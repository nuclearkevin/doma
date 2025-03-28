#include "TransportSolver1D.h"

#include "Parallel.h"

template <typename T>
TransportSolver1D<T>::TransportSolver1D(BrickMesh1D & mesh, const InputParameters & params, bool verbose, unsigned int num_threads)
  : _num_groups(params._num_e_groups),
    _mode(params._mode),
    _mesh(mesh),
    _angular_quad(2u * params._num_polar),
    _eq_system(),
    _verbose(verbose),
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
    _num_threads(num_threads)
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

    cell._current_iteration_source = 0.0;
    cell._current_scalar_flux = 0.0;
    for (unsigned int tid = 0; tid < _num_threads; ++tid)
    {
      cell.setSweptFlux(0.0, tid);
      cell.setAllInterfaceFluxes(0.0, tid);
    }
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

    cell._current_iteration_source = 0.0;
    cell._current_scalar_flux = 0.0;
    for (unsigned int tid = 0; tid < _num_threads; ++tid)
    {
      cell.setSweptFlux(0.0, tid);
      cell.setAllInterfaceFluxes(0.0, tid);
    }
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

    cell._current_iteration_source = p._g_src.size() != 0u ? 0.5 * p._g_src[g] : 0.0;

    // Accumulate the transient step source.
    if (_mode == RunMode::Transient && cell.hasStepSource())
    {
      const auto & ss = cell.getSourceStep();
      switch (ss._type)
      {
        case StepType::Both:
          cell._current_iteration_source += ss._insert_time <= t && t <= ss._remove_time ? 0.5 * ss._g_src[g] : 0.0;
          break;
        case StepType::Insert:
          cell._current_iteration_source += ss._insert_time <= t ? 0.5 * ss._g_src[g] : 0.0;
          break;
        case StepType::Remove:
          cell._current_iteration_source += t <= ss._remove_time ? 0.5 * ss._g_src[g] : 0.0;
          break;
      }
    }

    for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
    {
      // Accumulate the in-scattering contribution.
      if (g_prime != g)
        cell._current_iteration_source += 0.5 * p._g_g_scatter_mat[g * _num_groups + g_prime] * cell._total_scalar_flux[g_prime];

      // Accumulate the prompt fission contribution.
      if (p._g_chi_p.size() > 0u && p._num_d_groups == 0u)
        cell._current_iteration_source += 0.5 * p._g_chi_p[g] * p._g_prod[g_prime] * cell._total_scalar_flux[g_prime];
      else if (p._g_chi_p.size() > 0u && p._num_d_groups > 0u)
      {
        double g_beta = 0.0;
        for (unsigned int d = 0u; d < p._num_d_groups; ++d)
          g_beta += p._g_n_beta[g_prime * p._num_d_groups + d];

        cell._current_iteration_source += 0.5 * (1.0 - g_beta) * p._g_chi_p[g] * p._g_prod[g_prime] * cell._total_scalar_flux[g_prime];
      }

      // Accumulate the contribution from delayed neutrons.
      if (p._num_d_groups > 0u && _mode == RunMode::Transient)
        for (unsigned int d = 0u; d < p._num_d_groups; ++d)
          cell._current_iteration_source += 0.5 * p._n_g_chi_d[g * p._num_d_groups + d] * cell._current_t_dnps[d] * p._n_lambda[d];

      // Accumulate the transient source.
      if (_mode == RunMode::Transient)
        cell._current_iteration_source += 0.5 * cell._last_t_scalar_flux[g] * p._g_inv_v[g] / _dt;
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

    cell._current_scalar_flux = 0.0;

    for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
    {
      // Accumulate the in-scattering contribution.
      if (g_prime != g)
        cell._current_iteration_source += 0.5 * p._g_g_scatter_mat[g * _num_groups + g_prime] * cell._total_scalar_flux[g_prime];

      // Accumulate the fission source scaled by k_{eff}.
      if (p._g_chi_p.size() > 0u)
        cell._current_iteration_source += 0.5 * p._g_chi_p[g] * p._g_prod[g_prime] * cell._total_scalar_flux[g_prime] / _k;
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
      den = den + cell._l_x * p._g_prod[g] * cell._prev_mg_scalar_flux[g] / _k;
    }
  }

  _k_prev = _k;
  _k = num / den;
}

template <typename T>
bool
TransportSolver1D<T>::sourceIteration(unsigned int g)
{
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

    cell._current_iteration_source = 0.0;
    cell._current_scalar_flux = 0.0;
    for (unsigned int tid = 0; tid < _num_threads; ++tid)
    {
      cell.setSweptFlux(0.0, tid);
      cell.setAllInterfaceFluxes(0.0, tid);
    }
  }
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

    cell._total_scalar_flux[g] += cell._current_scalar_flux;
    cell._current_iteration_source = 0.5 * p._g_g_scatter_mat[g * _num_groups + g] * cell._current_scalar_flux;
    cell._current_scalar_flux = 0.0;
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
    diff_L2  = diff_L2 + std::pow(cell._current_scalar_flux, 2.0) * cell._l_x;
    total_L2 = total_L2 + std::pow(cell._total_scalar_flux[g] + cell._current_scalar_flux, 2.0) * cell._l_x;
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
    sweepR(std::abs(_angular_quad.direction(n)), _angular_quad.weight(n), n, g, omp_get_thread_num());
  }

  // Sweep -\mu.
  #pragma omp parallel for
  for (unsigned int n = _angular_quad.order() / 2u; n < _angular_quad.order(); ++n)
  {
    // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
    // if the upwind/downwind directions are aligned properly.
    sweepL(std::abs(_angular_quad.direction(n)), _angular_quad.weight(n), n, g, omp_get_thread_num());
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
    }
  }
}

template <typename T>
void
TransportSolver1D<T>::sweepR(const double & abs_mu, const double & weight,
                             unsigned int ordinate_index, unsigned int g, unsigned int tid)
{
  for (unsigned int row = 0u; row < _mesh._tot_num_x; ++row)
  {
    auto & cell = _mesh._cells[row];
    _eq_system.solve(cell, weight, abs_mu, ordinate_index, g,
                     CertesianFaceSide::Left, CertesianFaceSide::Right, tid);
  }
}

template <typename T>
void
TransportSolver1D<T>::sweepL(const double & abs_mu, const double & weight,
                             unsigned int ordinate_index, unsigned int g, unsigned int tid)
{
  unsigned int row = _mesh._tot_num_x;
  while (row --> 0)
  {
    auto & cell =_mesh._cells[row];
    _eq_system.solve(cell, weight, abs_mu, ordinate_index, g,
                     CertesianFaceSide::Right, CertesianFaceSide::Left, tid);
  }
}

template class TransportSolver1D<DiamondDifference1D>;
