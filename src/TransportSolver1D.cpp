#include "TransportSolver1D.h"

#include "Parallel.h"

#include "Eigen/SparseLU"

template <typename T>
TransportSolver1D<T>::TransportSolver1D(BrickMesh1D & mesh, const InputParameters & params, bool verbose, bool dsa, unsigned int num_threads)
  : _num_groups(params._num_e_groups),
    _mode(params._mode),
    _mesh(mesh),
    _angular_quad(2u * params._num_polar),
    _eq_system(),
    _verbose(verbose),
    _use_dsa(dsa),
    _sit(params._src_it_tol),
    _smi(params._num_src_it),
    _mgt(params._gs_tol),
    _mgi(params._num_mg_it),
    _k(1.0),
    _k_prev(0.0),
    _k_tol(params._k_tol),
    _t0(params._t0),
    _dt((params._t1 - params._t0) / static_cast<double>(params._num_steps)),
    _t_steps(params._num_steps),
    _ic(params._ic),
    _num_threads(num_threads),
    _diffusion_mat(_mesh._cells.size() + 1, _mesh._cells.size() + 1),
    _diffusion_src_vec(_mesh._cells.size() + 1),
    _diffusion_fluxes(_mesh._cells.size() + 1)
{
  _mesh._num_groups = _num_groups;
}

template <typename T>
bool
TransportSolver1D<T>::solveEigenvalue(const std::string & output_file_base)
{
  std::cout << "Solving..." << std::endl;
  initializeSolve();

  unsigned int mg_iteration = 0u;
  double current_shape_residual = 0.0;
  double current_k_diff = 0.0;
  do
  {
    std::cout << "Performing PI " << mg_iteration;
    if (mg_iteration != 0u)
      std::cout << ", shape residual = " << current_shape_residual
                << ", k residual = " << current_k_diff << std::endl
                << "k = " << _k << std::endl;
    else
      std::cout << std::endl;
    std::cout << "----------------------------------------------------------------------------------"
              << std::endl;

    for (unsigned int g = 0u; g < _num_groups; ++g)
    {
      updateMultigroupSourceEigen(g);
      auto res = sourceIteration(g);
      if (!res)
        return res;
    }

    updateEigenvalue();

    current_shape_residual = computeMGFluxResidual();
    current_k_diff = std::abs(_k - _k_prev);

    mg_iteration++;
  } while (mg_iteration < _mgi && (_mgt < current_shape_residual || _k_tol < current_k_diff));

  if (mg_iteration < _mgi)
  {
    std::cout << "PI converged after " << mg_iteration
              << " iterations with a flux shape residual of " << current_shape_residual
              << " and a k residual of " << current_k_diff << "." << std::endl
              << "Converged value of k = " << _k << std::endl;

    if (output_file_base != "")
      _mesh.dumpToTextFile(output_file_base);

    return true;
  }
  else
  {
    std::cout << "PI failed to converge after " << mg_iteration
              << " iterations with a flux shape residual of " << current_shape_residual
              << " and a k residual of " << current_k_diff << "." << std::endl;
    return false;
  }
}

template <typename T>
bool
TransportSolver1D<T>::solveTransient(const std::string & output_file_base)
{
  // Compute initial conditions.
  bool res = true;
  std::cout << "Setting up initial conditions..." << std::endl;
  switch (_ic)
  {
    case TransientIC::Zero:
      initZeroIC();
      break;
    case TransientIC::SteadyState:
      res = initSteadyIC();
      break;
    default:
      break;
  }

  // Failed to initialize with a steady-state calculation, abort.
  if (!res)
    return res;

  // Dump IC to file.
  _mesh.dumpToTextFile(output_file_base);

  double t = 0.0;
  for (unsigned int step = 0u; step < _t_steps; ++step)
  {
    std::cout << "Solving timestep " << step << "..." << std::endl;
    res = solveFixedSource("", t);

    // Failed to solve the steady-state problem at the current timestep, abort.
    if (!res)
      return res;

    // Step the delayed neutron precursors.
    stepDNPs();

    _mesh.dumpToTextFile(output_file_base + "_t" + std::to_string(step), true);
    _mesh.dumpDNPsToTextFile(output_file_base + "_t" + std::to_string(step));

    // Copy the curernt timestep's fluxes into the previous step's and zero the current iteration's fluxes.
    if (step != _t_steps - 1)
      updateStepFluxes();

    t += _dt;
  }

  return true;
}

template <typename T>
void
TransportSolver1D<T>::initZeroIC()
{
  #pragma omp parallel for
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
  {
    auto & cell = _mesh._cells[i];
    const auto & p = cell.getMatProps();

    cell._total_scalar_flux.resize(_num_groups, 0.0);
    cell._prev_mg_scalar_flux.resize(_num_groups, 0.0);
    cell._last_t_scalar_flux.resize(_num_groups, 0.0);
    cell._current_t_dnps.resize(p._num_d_groups, 0.0);
    cell._last_t_dnps.resize(p._num_d_groups, 0.0);
  }
}

template <typename T>
bool
TransportSolver1D<T>::initSteadyIC()
{
  initializeSolve();

  // Run a steady-state solve to compute initial conditions.
  auto res = solveFixedSource();

  // Failed to solve the steady-state problem at the current timestep, abort.
  if (!res)
    return res;

  #pragma omp parallel for
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
  {
    auto & cell = _mesh._cells[i];
    const auto & p = cell.getMatProps();

    // Init DNPs based on the steady-state equation.
    if (p._num_d_groups > 0)
    {
      cell._current_t_dnps.resize(p._num_d_groups, 0.0);
      cell._last_t_dnps.resize(p._num_d_groups, 0.0);

      for (unsigned int d = 0; d < p._num_d_groups; ++d)
      {
        for (unsigned int g = 0; g < _num_groups; ++g)
          cell._last_t_dnps[d] += p._g_prod[g] * cell._total_scalar_flux[g] * p._g_n_beta[g * p._num_d_groups + d];

        cell._last_t_dnps[d] /= p._n_lambda[d];
      }
    }

    // Copy the steady state solution into the previous timestep vector.
    std::copy(cell._total_scalar_flux.begin(),
              cell._total_scalar_flux.end(),
              std::back_inserter(cell._last_t_scalar_flux));

    for (unsigned int g = 0; g < _num_groups; ++g)
      cell._total_scalar_flux[g] = 0.0;

    cell._mg_source = 0.0;
    cell._current_iteration_source = 0.0;
    cell._current_scalar_flux = 0.0;
    cell._current_current = 0.0;
    cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)] = 0;
    cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)] = 0;
    for (unsigned int tid = 0; tid < _num_threads; ++tid)
    {
      cell.setSweptFlux(0.0, tid);
      cell.setAllInterfaceFluxes(0.0, tid);
      // For DSA.
      cell.setSweptCurrent(0.0, tid);
    }
  }
  #pragma omp parallel for
  for (unsigned int tid = 0; tid < _num_threads; ++tid)
  {
    for (unsigned int i = 0; i < _mesh._cells.size() + 1; ++i)
      _mesh._interface_scalar_fluxes[tid][i] = 0.0;
  }

  return res;
}

template <typename T>
void
TransportSolver1D<T>::updateStepFluxes()
{
  #pragma omp parallel for
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
  {
    auto & cell = _mesh._cells[i];
    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      cell._last_t_scalar_flux[g] = cell._total_scalar_flux[g];
      cell._total_scalar_flux[g] = 0.0;
    }

    cell._mg_source = 0.0;
    cell._current_iteration_source = 0.0;
    cell._current_scalar_flux = 0.0;
    cell._current_current = 0.0;
    cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)] = 0;
    cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)] = 0;
    for (unsigned int tid = 0; tid < _num_threads; ++tid)
    {
      cell.setSweptFlux(0.0, tid);
      cell.setAllInterfaceFluxes(0.0, tid);
      // For DSA.
      cell.setSweptCurrent(0.0, tid);
    }
  }
  #pragma omp parallel for
  for (unsigned int tid = 0; tid < _num_threads; ++tid)
  {
    for (unsigned int i = 0; i < _mesh._cells.size() + 1; ++i)
      _mesh._interface_scalar_fluxes[tid][i] = 0.0;
  }
}

template <typename T>
void
TransportSolver1D<T>::stepDNPs()
{
  #pragma omp parallel for
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
  {
    auto & cell = _mesh._cells[i];
    const auto & p = cell.getMatProps();

    for (unsigned int d = 0u; d < p._num_d_groups; ++d)
    {
      // Accumulate the DNP fission source for the current + previous timestep.
      double c_src = 0.0;
      for (unsigned int g = 0u; g < _num_groups; ++g)
        c_src += p._g_prod[g] * cell._total_scalar_flux[g] * p._g_n_beta[g * p._num_d_groups + d];

      // Compute the current step DNP concentrations. This expression is obtained by analytically
      // solving the precursor ODE while assuming the DNP source from fission is a constant between steps.
      const auto decay = std::exp(-p._n_lambda[d] * _dt);
      cell._current_t_dnps[d] = (c_src / p._n_lambda[d]) * (1.0 - decay) + cell._last_t_dnps[d] * decay;
      cell._last_t_dnps[d] = cell._current_t_dnps[d];
    }
  }
}

template <typename T>
bool
TransportSolver1D<T>::solveFixedSource(const std::string & output_file_base, const double & t)
{
  if (_mode != RunMode::Transient)
  {
    std::cout << "Solving..." << std::endl;
    initializeSolve();
  }

  unsigned int mg_iteration = 0u;
  double current_residual = 0.0;
  do
  {
    std::cout << "Performing MGI " << mg_iteration;
    if (mg_iteration != 0u)
      std::cout <<  ", residual = " << current_residual << std::endl;
    else
      std::cout << std::endl;
    std::cout << "----------------------------------------------------------------------------------"
              << std::endl;

    for (unsigned int g = 0u; g < _num_groups; ++g)
    {
      updateMultigroupSource(g, t);
      auto res = sourceIteration(g);
      if (!res)
        return res;
    }

    current_residual = computeMGFluxResidual();

    mg_iteration++;
  } while (mg_iteration < _mgi && _mgt < current_residual && _num_groups > 1u);

  if (mg_iteration < _mgi)
  {
    if (_num_groups > 1u)
    {
      std::cout << "MGI converged after " << mg_iteration
                << " iterations with a residual of " << current_residual << std::endl;
    }

    if (output_file_base != "")
      _mesh.dumpToTextFile(output_file_base);

    return true;
  }
  else
  {
    std::cout << "MGI failed to converge after " << mg_iteration
              << " iterations with a residual of " << current_residual << std::endl;
    return false;
  }
}

template <typename T>
double
TransportSolver1D<T>::computeMGFluxResidual()
{
  double diff_L2 = 0.0;
  double total_L2 = 0.0;

  #pragma omp parallel for reduction(+:diff_L2,total_L2)
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
  {
    const auto & cell = _mesh._cells[i];
    for (unsigned int g = 0u; g < _num_groups; ++g)
    {
      diff_L2  = diff_L2 + std::pow((cell._total_scalar_flux[g] - cell._prev_mg_scalar_flux[g]), 2.0) * cell._l_x;
      total_L2 = total_L2 + std::pow(cell._total_scalar_flux[g], 2.0) * cell._l_x;
    }
  }

  return total_L2 > 1e-8 ? std::sqrt(diff_L2) / std::sqrt(total_L2) : 0.0;
}

template <typename T>
void
TransportSolver1D<T>::updateMultigroupSource(unsigned int g, double t)
{
  #pragma omp parallel for
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
  {
    auto & cell = _mesh._cells[i];
    const auto & p = cell.getMatProps();

    cell._current_scalar_flux = 0.0;
    cell._current_current = 0.0;
    cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)] = 0;
    cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)] = 0;

    cell._mg_source = p._g_src.size() != 0u ? 0.5 * p._g_src[g] : 0.0;

    // Accumulate the transient step source.
    if (_mode == RunMode::Transient && cell.hasStepSource())
    {
      const auto & ss = cell.getSourceStep();
      switch (ss._type)
      {
        case StepType::Both:
          cell._mg_source += ss._insert_time <= t && t <= ss._remove_time ? 0.5 * ss._g_src[g] : 0.0;
          break;
        case StepType::Insert:
          cell._mg_source += ss._insert_time <= t ? 0.5 * ss._g_src[g] : 0.0;
          break;
        case StepType::Remove:
          cell._mg_source += t <= ss._remove_time ? 0.5 * ss._g_src[g] : 0.0;
          break;
      }
    }

    for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
    {
      // Accumulate the in-scattering contribution.
      if (g_prime != g)
        cell._mg_source += 0.5 * p._g_g_scatter_mat[g * _num_groups + g_prime] * cell._total_scalar_flux[g_prime];

      // Accumulate the prompt fission contribution.
      if (p._g_chi_p.size() > 0u && p._num_d_groups == 0u)
        cell._mg_source += 0.5 * p._g_chi_p[g] * p._g_prod[g_prime] * cell._total_scalar_flux[g_prime];
      else if (p._g_chi_p.size() > 0u && p._num_d_groups > 0u)
      {
        double g_beta = 0.0;
        for (unsigned int d = 0u; d < p._num_d_groups; ++d)
          g_beta += p._g_n_beta[g_prime * p._num_d_groups + d];

        cell._mg_source += 0.5 * (1.0 - g_beta) * p._g_chi_p[g] * p._g_prod[g_prime] * cell._total_scalar_flux[g_prime];
      }

      // Accumulate the contribution from delayed neutrons.
      if (p._num_d_groups > 0u && _mode == RunMode::Transient)
        for (unsigned int d = 0u; d < p._num_d_groups; ++d)
          cell._mg_source += 0.5 * p._n_g_chi_d[g * p._num_d_groups + d] * cell._current_t_dnps[d] * p._n_lambda[d];

      // Accumulate the transient source.
      if (_mode == RunMode::Transient)
        cell._mg_source += 0.5 * cell._last_t_scalar_flux[g] * p._g_inv_v[g] / _dt;
    }

    cell._prev_mg_scalar_flux[g] = cell._total_scalar_flux[g];
    cell._total_scalar_flux[g] = 0.0;
  }
}

template <typename T>
void
TransportSolver1D<T>::updateMultigroupSourceEigen(unsigned int g)
{
  #pragma omp parallel for
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
  {
    auto & cell = _mesh._cells[i];
    const auto & p = cell.getMatProps();

    cell._mg_source = 0.0;
    cell._current_scalar_flux = 0.0;
    cell._current_current = 0.0;
    cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)] = 0;
    cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)] = 0;

    for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
    {
      // Accumulate the in-scattering contribution.
      if (g_prime != g)
        cell._mg_source += 0.5 * p._g_g_scatter_mat[g * _num_groups + g_prime] * cell._total_scalar_flux[g_prime];

      // Accumulate the fission source scaled by k_{eff}.
      if (p._g_chi_p.size() > 0u)
        cell._mg_source += 0.5 * p._g_chi_p[g] * p._g_prod[g_prime] * cell._total_scalar_flux[g_prime] / _k;
    }

    cell._prev_mg_scalar_flux[g] = cell._total_scalar_flux[g];
    cell._total_scalar_flux[g] = 0.0;
  }
}

template <typename T>
void
TransportSolver1D<T>::updateEigenvalue()
{
  double num = 0.0;
  double den = 0.0;

  #pragma omp parallel for reduction(+:num,den)
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
  {
    const auto & cell = _mesh._cells[i];
    const auto & p = cell.getMatProps();
    if (p._g_prod.size() == 0u)
      continue;

    for (unsigned int g = 0u; g < _num_groups; ++g)
    {
      num = num + cell._l_x * p._g_prod[g] * cell._total_scalar_flux[g];
      den = den + cell._l_x * p._g_prod[g] * cell._prev_mg_scalar_flux[g];
    }
  }

  _k_prev = _k;
  _k = _k_prev * num / den;
}

template <typename T>
bool
TransportSolver1D<T>::sourceIteration(unsigned int g)
{
  // First iteration scattering source.
  updateScatteringSource(g);

  unsigned int source_iteration = 0u;
  double current_residual = 0.0;
  do
  {
    if (_verbose)
    {
      std::cout << "Performing SI " << source_iteration << " for G" << g;
      if (source_iteration != 0u)
        std::cout <<  ", residual = " << current_residual << std::endl;
      else
        std::cout << std::endl;
    }

    sweep(g);

    if (_use_dsa)
      syntheticAcceleration(g);

    current_residual = computeScatteringResidual(g);
    updateScatteringSource(g);

    source_iteration++;
  }
  while (source_iteration < _smi && _sit < current_residual);

  if (source_iteration < _smi)
  {
    std::cout << "SI converged after " << source_iteration
              << " iterations with a residual of " << current_residual << " for G"
              << g << std::endl;
    return true;
  }
  else
  {
    std::cout << "SI failed to converge after " << source_iteration
              << " iterations with a residual of " << current_residual << " for G"
              << g << std::endl;
    return false;
  }
}

template <typename T>
void
TransportSolver1D<T>::initializeSolve()
{
  std::cout << "Initializing the solver..." << std::endl;

  // Initializing the boundary condition data structure.
  for (unsigned int i = 0u; i < 2u; ++i)
    if (_mesh._bcs[i] != BoundaryCondition::Vacuum)
      _mesh._boundary_angular_fluxes[i].resize(_mesh._boundary_cells[i].size() * _angular_quad.order(), 0.0);

  #pragma omp parallel for
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
  {
    auto & cell = _mesh._cells[i];
    switch (_mode)
    {
      case RunMode::Eigen: cell._total_scalar_flux.resize(_num_groups, 1.0); break;
      default:             cell._total_scalar_flux.resize(_num_groups, 0.0); break;
    }
    cell._prev_mg_scalar_flux.resize(_num_groups, 0.0);

    cell._mg_source = 0.0;
    cell._current_iteration_source = 0.0;
    cell._current_scalar_flux = 0.0;
    cell._current_current = 0.0;
    cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)] = 0;
    cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)] = 0;
    for (unsigned int tid = 0; tid < _num_threads; ++tid)
    {
      cell.setSweptFlux(0.0, tid);
      cell.setAllInterfaceFluxes(0.0, tid);
      // For DSA.
      cell.setSweptCurrent(0.0, tid);
    }
  }
  #pragma omp parallel for
  for (unsigned int tid = 0; tid < _num_threads; ++tid)
  {
    for (unsigned int i = 0; i < _mesh._cells.size() + 1; ++i)
      _mesh._interface_scalar_fluxes[tid][i] = 0.0;
  }
}

/**
 * Based on "Diffusion Synthetic Acceleration Methods for the Diamond-Differenced Discrete-Ordinates Equations"
 * by R. E. Alcouffe.
 * https://doi.org/10.13182/NSE77-1
 */
template <typename T>
void
TransportSolver1D<T>::syntheticAcceleration(unsigned int g)
{
  _diffusion_mat_entries.clear();
  _diffusion_mat.setZero();

  const unsigned int g_prime = _num_groups * g + g;
  // Matrix interior is different between schemes. Boundary conditions are the same.
  switch (_mode)
  {
    // For eigenvalue calculations, we use the diffusion correction scheme unless it fails.
    // If it does we switch to removal correction.
    case RunMode::Eigen:
    {
      bool swap_to_removal = false;

      for (unsigned int i = 0; i < _mesh._cells.size() - 1; ++i)
      {
        // Current cell at i.
        const auto & c_i   = _mesh._cells[i];
        const auto & p_i   = c_i.getMatProps();
        // Neighboring cell at i + 1.
        const auto & c_i_1 = _mesh._cells[i + 1];
        const auto & p_i_1 = c_i_1.getMatProps();

        const unsigned int l_half = i;     // Left side of the tridiagonal
        const unsigned int m_half = i + 1; // Center of the tridiagonal.
        const unsigned int r_half = i + 2; // Right side of the tridiagonal.

        // Modified diffusion coefficients.
        double D_i = -1.0 * c_i._l_x * c_i._current_current;
        D_i /= (c_i._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)]
                - c_i._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)]);

        double D_i_1 = -1.0 * c_i_1._l_x * c_i_1._current_current;
        D_i_1 /= (c_i_1._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)]
                  - c_i_1._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)]);

        // If negative diffusion coefficients are encountered, exit early and use
        // the removal correction.
        if (D_i <= 0.0 || D_i_1 <= 0.0)
        {
          swap_to_removal = true;
          break;
        }

        // Eq. 23a
        double sigma_r_half = (p_i._g_total[g] - p_i._g_g_scatter_mat[g_prime]) * c_i._l_x * c_i._current_scalar_flux;
        sigma_r_half += (p_i_1._g_total[g] - p_i_1._g_g_scatter_mat[g_prime]) * c_i_1._l_x * c_i_1._current_scalar_flux;
        sigma_r_half *= 0.5 / c_i._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)];

        // Diffusion correction matrix.
        _diffusion_mat_entries.push_back(Eigen::Triplet<double>(m_half, l_half, -1.0 * D_i / c_i._l_x));
        _diffusion_mat_entries.push_back(Eigen::Triplet<double>(m_half, m_half, (D_i / c_i._l_x) + (D_i_1 / c_i_1._l_x) + sigma_r_half));
        _diffusion_mat_entries.push_back(Eigen::Triplet<double>(m_half, r_half, -1.0 * D_i_1 / c_i_1._l_x));

        // Eq. 23b. The source pre-multiplies by 0.5.
        const double qqh = (c_i._mg_source * c_i._l_x + c_i_1._mg_source * c_i_1._l_x);
        _diffusion_src_vec(m_half) = qqh;
      }

      if (swap_to_removal)
      {
        _diffusion_mat_entries.clear();
        _diffusion_mat.setZero();

        for (unsigned int i = 0; i < _mesh._cells.size() - 1; ++i)
        {
          // Current cell at i.
          const auto & c_i   = _mesh._cells[i];
          const auto & p_i   = c_i.getMatProps();
          // Neighboring cell at i + 1.
          const auto & c_i_1 = _mesh._cells[i + 1];
          const auto & p_i_1 = c_i_1.getMatProps();

          const unsigned int l_half = i;     // Left side of the tridiagonal
          const unsigned int m_half = i + 1; // Center of the tridiagonal.
          const unsigned int r_half = i + 2; // Right side of the tridiagonal.

          // Eqs. 21 and 23 (diffusion coefficients divided by h).
          const double Dh_i   = 1.0 / 3.0 / p_i._g_total[g] / c_i._l_x;
          const double Dh_i_1 = 1.0 / 3.0 / p_i_1._g_total[g] / c_i_1._l_x;
          // Eq. 23a
          double sigma_r_half = (p_i._g_total[g] - p_i._g_g_scatter_mat[g_prime]) * c_i._l_x * c_i._current_scalar_flux;
          sigma_r_half += (p_i_1._g_total[g] - p_i_1._g_g_scatter_mat[g_prime]) * c_i_1._l_x * c_i_1._current_scalar_flux;
          sigma_r_half *= 0.5 / c_i._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)];

          // Removal correction term.
          double removal_corr = c_i_1._current_current - c_i._current_current;
          removal_corr += Dh_i_1 * (c_i_1._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)]
                                    - c_i_1._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)]);
          removal_corr -= Dh_i * (c_i._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)]
                                  - c_i._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)]);
          removal_corr /= c_i._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)];

          // Diffusion correction matrix.
          _diffusion_mat_entries.push_back(Eigen::Triplet<double>(m_half, l_half, -1.0 * Dh_i));
          _diffusion_mat_entries.push_back(Eigen::Triplet<double>(m_half, m_half, Dh_i + Dh_i_1 + sigma_r_half + removal_corr));
          _diffusion_mat_entries.push_back(Eigen::Triplet<double>(m_half, r_half, -1.0 * Dh_i_1));

          // Eq. 23b. The source pre-multiplies by 0.5.
          const double qqh = (c_i._mg_source * c_i._l_x + c_i_1._mg_source * c_i_1._l_x);
          _diffusion_src_vec(m_half) = qqh;
        }
      }

      break;
    }
    // For fixed source calculations, we use the source correction scheme.
    case RunMode::FixedSrc:
    case RunMode::Transient:
    {
      for (unsigned int i = 0; i < _mesh._cells.size() - 1; ++i)
      {
        // Current cell at i.
        const auto & c_i   = _mesh._cells[i];
        const auto & p_i   = c_i.getMatProps();
        // Neighboring cell at i + 1.
        const auto & c_i_1 = _mesh._cells[i + 1];
        const auto & p_i_1 = c_i_1.getMatProps();

        const unsigned int l_half = i;     // Left side of the tridiagonal
        const unsigned int m_half = i + 1; // Center of the tridiagonal.
        const unsigned int r_half = i + 2; // Right side of the tridiagonal.

        // Eqs. 21 and 23 (diffusion coefficients divided by h).
        const double Dh_i   = 1.0 / 3.0 / p_i._g_total[g] / c_i._l_x;
        const double Dh_i_1 = 1.0 / 3.0 / p_i_1._g_total[g] / c_i_1._l_x;
        // Eq. 23a
        double sigma_r_half = (p_i._g_total[g] - p_i._g_g_scatter_mat[g_prime]) * c_i._l_x * c_i._current_scalar_flux;
        sigma_r_half += (p_i_1._g_total[g] - p_i_1._g_g_scatter_mat[g_prime]) * c_i_1._l_x * c_i_1._current_scalar_flux;
        sigma_r_half *= 0.5 / c_i._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)];

        // Eq. 23 proper.
        _diffusion_mat_entries.push_back(Eigen::Triplet<double>(m_half, l_half, -1.0 * Dh_i));
        _diffusion_mat_entries.push_back(Eigen::Triplet<double>(m_half, m_half, Dh_i + Dh_i_1 + sigma_r_half));
        _diffusion_mat_entries.push_back(Eigen::Triplet<double>(m_half, r_half, -1.0 * Dh_i_1));

        // Eq. 23b. The source pre-multiplies by 0.5.
        const double qqh = (c_i._mg_source * c_i._l_x + c_i_1._mg_source * c_i_1._l_x);

        // Eq. 22 (DSA source correction).
        double r = c_i_1._current_current - c_i._current_current;
        r += Dh_i_1 * (c_i_1._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)]
                       - c_i_1._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)]);
        r -= Dh_i * (c_i._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)]
                     - c_i._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)]);
        _diffusion_src_vec(m_half) = qqh + r;
      }
      break;
    }
  }

  /**
   * BCs are common to all schemes.
   * Left boundary condition:  Eq. 34a and Eq. 37a.
   * Right boundary condition: Eq. 34b and Eq. 37b.
   * "Unconditionally Stable Diffusion-Synthetic Acceleration Methods
   * for the Slab Geometry Discrete Ordinates Equations. Part I: Theory"
   * by E. W. Larsen. https://doi.org/10.13182/NSE82-1
   */
  {
    const auto & c_0 = _mesh._cells[0];
    // The source pre-multiplies by 0.5.
    _diffusion_src_vec(0) = c_0._mg_source;

    const auto & p_0 = c_0.getMatProps();
    const double half_d1_h   = 0.5 / 3.0 / p_0._g_total[g] / c_0._l_x;
    const double eight_sr1_h = 0.125 * (p_0._g_total[g] - p_0._g_g_scatter_mat[g_prime]) * c_0._l_x;

    double beta = 0.0;
    for (unsigned int n = 0u; n < _angular_quad.order() / 2u; ++n)
      beta += _angular_quad.direction(n) * _angular_quad.weight(n);

    _diffusion_mat_entries.push_back(Eigen::Triplet<double>(0, 0, beta + half_d1_h + eight_sr1_h));
    _diffusion_mat_entries.push_back(Eigen::Triplet<double>(0, 1, -1.0 * half_d1_h + eight_sr1_h));
  }
  {
    const auto & c_I = _mesh._cells[_mesh._cells.size() - 1];
    // The source pre-multiplies by 0.5.
    _diffusion_src_vec(_mesh._cells.size()) = -1.0 * c_I._mg_source;

    const auto & p_I = c_I.getMatProps();
    const double half_dI_h   = 0.5 / 3.0 / p_I._g_total[g] / c_I._l_x;
    const double eight_srI_h = 0.125 * (p_I._g_total[g] - p_I._g_g_scatter_mat[g_prime]) * c_I._l_x;

    double beta = 0.0;
    for (unsigned int n = _angular_quad.order() / 2u; n < _angular_quad.order(); ++n)
      beta += _angular_quad.direction(n) * _angular_quad.weight(n);

    _diffusion_mat_entries.push_back(Eigen::Triplet<double>(_mesh._cells.size(), _mesh._cells.size(), beta - half_dI_h - eight_srI_h));
    _diffusion_mat_entries.push_back(Eigen::Triplet<double>(_mesh._cells.size(), _mesh._cells.size() - 1, half_dI_h - eight_srI_h));
  }

  // Build the matrix.
  _diffusion_mat.setFromTriplets(_diffusion_mat_entries.begin(), _diffusion_mat_entries.end());
  _diffusion_mat.makeCompressed();

  // Solve it!
  _diffusion_solver.analyzePattern(_diffusion_mat);
  _diffusion_solver.factorize(_diffusion_mat);
  _diffusion_fluxes = _diffusion_solver.solve(_diffusion_src_vec);

  // Update transport guess with consistent diffusion fluxes.
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
    _mesh._cells[i]._current_scalar_flux = 0.5 * (_diffusion_fluxes(i) + _diffusion_fluxes(i + 1));
}

template <typename T>
void
TransportSolver1D<T>::updateScatteringSource(unsigned int g)
{
  #pragma omp parallel for
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
  {
    auto & cell = _mesh._cells[i];
    const auto & p = cell.getMatProps();

    cell._total_scalar_flux[g] = cell._current_scalar_flux;
    cell._current_iteration_source = cell._mg_source;
    cell._current_iteration_source += 0.5 * p._g_g_scatter_mat[g * _num_groups + g] * cell._current_scalar_flux;
    cell._current_scalar_flux = 0.0;
    cell._current_current = 0.0;
    cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)] = 0;
    cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)] = 0;
  }
}

template <typename T>
double
TransportSolver1D<T>::computeScatteringResidual(unsigned int g)
{
  double diff_L2 = 0.0;
  double total_L2 = 0.0;

  #pragma omp parallel for reduction(+:diff_L2,total_L2)
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
  {
    const auto & cell = _mesh._cells[i];
    diff_L2  = diff_L2 + std::pow(cell._current_scalar_flux - cell._total_scalar_flux[g], 2.0) * cell._l_x;
    total_L2 = total_L2 + std::pow(cell._current_scalar_flux, 2.0) * cell._l_x;
  }

  return total_L2 > 1e-8 ? std::sqrt(diff_L2) / std::sqrt(total_L2) : 0.0;
}

template <typename T>
void
TransportSolver1D<T>::sweep(unsigned int g)
{
  // Sweep +\mu.
  #pragma omp parallel for
  for (unsigned int n = 0u; n < _angular_quad.order() / 2u; ++n)
  {
    // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
    // if the upwind/downwind directions are aligned properly.
    sweepR(_angular_quad.direction(n), _angular_quad.weight(n), n, g, omp_get_thread_num());
  }

  // Sweep -\mu.
  #pragma omp parallel for
  for (unsigned int n = _angular_quad.order() / 2u; n < _angular_quad.order(); ++n)
  {
    // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
    // if the upwind/downwind directions are aligned properly.
    sweepL(_angular_quad.direction(n), _angular_quad.weight(n), n, g, omp_get_thread_num());
  }

  // Accumulate all of the swept flux moments across all threads.
  #pragma omp parallel for
  for (unsigned int i = 0; i < _mesh._cells.size(); ++i)
  {
    auto & cell = _mesh._cells[i];

    for (unsigned int tid = 0; tid < _num_threads; ++tid)
    {
      cell._current_scalar_flux += cell.getSweptFlux(tid);
      cell.setSweptFlux(0.0, tid);
      cell.setAllInterfaceFluxes(0.0, tid);
      // For DSA.
      cell._current_current += cell.getSweptCurrent(tid);
      cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Right)]
        += cell.getSweptInterfaceSF(CertesianFaceSide::Right, tid);
      cell._current_interface_sf[static_cast<unsigned int>(CertesianFaceSide::Left)]
        += cell.getSweptInterfaceSF(CertesianFaceSide::Left, tid);
      cell.setSweptCurrent(0.0, tid);
    }
  }
  #pragma omp parallel for
  for (unsigned int tid = 0; tid < _num_threads; ++tid)
  {
    for (unsigned int i = 0; i < _mesh._cells.size() + 1; ++i)
      _mesh._interface_scalar_fluxes[tid][i] = 0.0;
  }
}

template <typename T>
void
TransportSolver1D<T>::sweepR(const double & mu, const double & weight,
                             unsigned int ordinate_index, unsigned int g, unsigned int tid)
{
  for (unsigned int row = 0u; row < _mesh._tot_num_x; ++row)
  {
    auto & cell = _mesh._cells[row];
    _eq_system.solve(cell, weight, mu, ordinate_index, g,
                     CertesianFaceSide::Left, CertesianFaceSide::Right, tid);
  }
}

template <typename T>
void
TransportSolver1D<T>::sweepL(const double & mu, const double & weight,
                             unsigned int ordinate_index, unsigned int g, unsigned int tid)
{
  unsigned int row = _mesh._tot_num_x;
  while (row --> 0)
  {
    auto & cell =_mesh._cells[row];
    _eq_system.solve(cell, weight, mu, ordinate_index, g,
                     CertesianFaceSide::Right, CertesianFaceSide::Left, tid);
  }
}

template class TransportSolver1D<DiamondDifference1D>;
