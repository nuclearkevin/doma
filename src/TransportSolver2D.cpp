#include "TransportSolver2D.h"

template <typename T>
TransportSolver2D<T>::TransportSolver2D(BrickMesh2D & mesh, const InputParameters & params, bool verbose)
  : _num_groups(params._num_e_groups),
    _mode(params._mode),
    _mesh(mesh),
    _angular_quad(2u * params._num_azimuthal, 2u * params._num_polar, 2u),
    _eq_system(),
    _verbose(verbose),
    _sit(params._src_it_tol),
    _smi(params._num_src_it),
    _mgt(params._gs_tol),
    _mgi(params._num_mg_it),
    _t0(params._t0),
    _dt((params._t1 - params._t0) / static_cast<double>(params._num_steps)),
    _t_steps(params._num_steps),
    _ic(params._ic)
{
  _mesh._num_groups = _num_groups;
}

template <typename T>
bool
TransportSolver2D<T>::solveTransient(const std::string & output_file_base)
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
TransportSolver2D<T>::initZeroIC()
{
  for (auto & cell : _mesh._cells)
  {
    const auto & p = cell.getMatProps();

    cell._total_scalar_flux.resize(_num_groups, 0.0);
    cell._last_t_scalar_flux.resize(_num_groups, 0.0);
    cell._current_t_dnps.resize(p._num_d_groups, 0.0);
    cell._last_t_dnps.resize(p._num_d_groups, 0.0);

    cell._current_iteration_source = 0.0;
    cell._current_scalar_flux = 0.0;
    cell._interface_angular_fluxes.fill(0.0);
  }
}

template <typename T>
bool
TransportSolver2D<T>::initSteadyIC()
{
  initializeSolve();

  // Run a steady-state solve to compute initial conditions.
  auto res = solveFixedSource();

  // Failed to solve the steady-state problem at the current timestep, abort.
  if (!res)
    return res;

  for (auto & cell : _mesh._cells)
  {
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
    cell._interface_angular_fluxes.fill(0.0);
  }

  return res;
}

template <typename T>
void
TransportSolver2D<T>::updateStepFluxes()
{
  for (auto & cell : _mesh._cells)
  {
    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      cell._last_t_scalar_flux[g] = cell._total_scalar_flux[g];
      cell._total_scalar_flux[g] = 0.0;
    }
    cell._current_iteration_source = 0.0;
    cell._current_scalar_flux = 0.0;
    cell._interface_angular_fluxes.fill(0.0);
  }
}

template <typename T>
void
TransportSolver2D<T>::stepDNPs()
{
  for (auto & cell : _mesh._cells)
  {
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
TransportSolver2D<T>::solveFixedSource(const std::string & output_file_base, const double & t)
{
  if (_mode != RunMode::Transient)
  {
    std::cout << "Solving..." << std::endl;
    initializeSolve();
  }

  unsigned int mg_iteration = 0u;
  double current_residual = 0.0;
  double previous_norm = 1.0;
  double current_norm = 0.0;
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

    current_norm = computeMGFluxNorm();
    current_residual = std::abs(current_norm - previous_norm) / previous_norm;
    previous_norm = current_norm;

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
TransportSolver2D<T>::computeMGFluxNorm()
{
  double norm = 0.0;
  for (auto & cell : _mesh._cells)
    for (unsigned int g = 0u; g < _num_groups; ++g)
      norm += std::pow(cell._total_scalar_flux[g], 2.0) * cell._area;

  return std::sqrt(norm);
}

template <typename T>
void
TransportSolver2D<T>::updateMultigroupSource(unsigned int g, double t)
{
  for (auto & cell : _mesh._cells)
  {
    const auto & p = cell.getMatProps();

    cell._current_scalar_flux = 0.0;

    cell._current_iteration_source = p._g_src.size() != 0u ? 0.5 * p._g_src[g] / M_PI : 0.0;

    // Accumulate the transient step source.
    if (_mode == RunMode::Transient && cell.hasStepSource())
    {
      const auto & ss = cell.getSourceStep();
      switch (ss._type)
      {
        case StepType::Both:
          cell._current_iteration_source += ss._insert_time <= t && t <= ss._remove_time ? 0.5 * ss._g_src[g] / M_PI : 0.0;
          break;
        case StepType::Insert:
          cell._current_iteration_source += ss._insert_time <= t ? 0.5 * ss._g_src[g] / M_PI : 0.0;
          break;
        case StepType::Remove:
          cell._current_iteration_source += t <= ss._remove_time ? 0.5 * ss._g_src[g] / M_PI : 0.0;
          break;
      }
    }

    for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
    {
      // Accumulate the in-scattering contribution.
      if (g_prime != g)
        cell._current_iteration_source += 0.5 * p._g_g_scatter_mat[g * _num_groups + g_prime] * cell._total_scalar_flux[g_prime] / M_PI;

      // Accumulate the prompt fission contribution.
      if (p._g_chi_p.size() > 0u && p._num_d_groups == 0u)
        cell._current_iteration_source += 0.5 * p._g_chi_p[g] * p._g_prod[g_prime] * cell._total_scalar_flux[g_prime] / M_PI;
      else if (p._g_chi_p.size() > 0u && p._num_d_groups > 0u)
      {
        double g_beta = 0.0;
        for (unsigned int d = 0u; d < p._num_d_groups; ++d)
          g_beta += p._g_n_beta[g_prime * p._num_d_groups + d];

        cell._current_iteration_source += 0.5 * (1.0 - g_beta) * p._g_chi_p[g] * p._g_prod[g_prime] * cell._total_scalar_flux[g_prime] / M_PI;
      }

      // Accumulate the contribution from delayed neutrons.
      if (p._num_d_groups > 0u && _mode == RunMode::Transient)
        for (unsigned int d = 0u; d < p._num_d_groups; ++d)
          cell._current_iteration_source += 0.5 * p._n_g_chi_d[g * p._num_d_groups + d] * cell._current_t_dnps[d] * p._n_lambda[d] / M_PI;

      // Accumulate the transient source.
      if (_mode == RunMode::Transient)
        cell._current_iteration_source += 0.5 * cell._last_t_scalar_flux[g] * p._g_inv_v[g] / _dt / M_PI;
    }
    cell._total_scalar_flux[g] = 0.0;
  }
}

template <typename T>
bool
TransportSolver2D<T>::sourceIteration(unsigned int g)
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
TransportSolver2D<T>::initializeSolve()
{
  std::cout << "Initializing the solver..." << std::endl;

  // Initializing the boundary condition data structure.
  for (unsigned int i = 0u; i < 4u; ++i)
    if (_mesh._bcs[i] != BoundaryCondition::Vacuum)
      _mesh._boundary_angular_fluxes[i].resize(_mesh._boundary_cells[i].size() * _angular_quad.totalOrder(), 0.0);

  for (auto & cell : _mesh._cells)
  {
    const auto & p = cell.getMatProps();
    cell._total_scalar_flux.resize(_num_groups, 0.0);

    cell._current_iteration_source = 0.0;
    cell._current_scalar_flux = 0.0;
    cell._interface_angular_fluxes.fill(0.0);
  }
}

template <typename T>
void
TransportSolver2D<T>::updateScatteringSource(unsigned int g)
{
  for (auto & cell : _mesh._cells)
  {
    const auto & p = cell.getMatProps();

    cell._total_scalar_flux[g] += cell._current_scalar_flux;
    cell._current_iteration_source = 0.5 * p._g_g_scatter_mat[g * _num_groups + g] * cell._current_scalar_flux / M_PI;
    cell._current_scalar_flux = 0.0;

    cell._interface_angular_fluxes.fill(0.0);
  }
}

template <typename T>
double
TransportSolver2D<T>::computeScatteringResidual(unsigned int g)
{
  double diff_L2 = 0.0;
  double total_L2 = 0.0;
  for (auto & cell : _mesh._cells)
  {
    diff_L2 += std::pow(cell._current_scalar_flux, 2.0) * cell._area;
    total_L2 += std::pow(cell._total_scalar_flux[g] + cell._current_scalar_flux, 2.0) * cell._area;
  }

  return total_L2 > 1e-8 ? std::sqrt(diff_L2) / std::sqrt(total_L2) : 0.0;
}

template <typename T>
void
TransportSolver2D<T>::sweep(unsigned int g)
{
  double mu = 0.0;
  double eta = 0.0;
  double xi = 0.0;

  double abs_mu = 0.0;
  double abs_eta = 0.0;
  double abs_xi = 0.0;

  double weight = 0.0;

  // Sweep +\mu, +\eta, +\xi.
  {
    const auto current_oct = Octant::PPP;
    for (unsigned int n = 0u; n < _angular_quad.order(current_oct); ++n)
    {
      _angular_quad.direction(current_oct, n, mu, eta, xi);
      weight = _angular_quad.weight(current_oct, n);

      // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
      // if the upwind/downwind directions are aligned properly.
      abs_mu = std::abs(mu);
      abs_eta = std::abs(eta);
      abs_xi = std::abs(xi);

      sweepPPP(abs_mu, abs_eta, abs_xi, weight, n, g);
    }
  }

  // Sweep +\mu, -\eta, +\xi.
  {
    const auto current_oct = Octant::PMP;
    for (unsigned int n = 0u; n < _angular_quad.order(current_oct); ++n)
    {
      _angular_quad.direction(current_oct, n, mu, eta, xi);
      weight = _angular_quad.weight(current_oct, n);

      // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
      // if the upwind/downwind directions are aligned properly.
      abs_mu = std::abs(mu);
      abs_eta = std::abs(eta);
      abs_xi = std::abs(xi);

      sweepPMP(abs_mu, abs_eta, abs_xi, weight, n, g);
    }
  }

  // Sweep -\mu, +\eta, +\xi.
  {
    const auto current_oct = Octant::MPP;
    for (unsigned int n = 0u; n < _angular_quad.order(current_oct); ++n)
    {
      _angular_quad.direction(current_oct, n, mu, eta, xi);
      weight = _angular_quad.weight(current_oct, n);

      // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
      // if the upwind/downwind directions are aligned properly.
      abs_mu = std::abs(mu);
      abs_eta = std::abs(eta);
      abs_xi = std::abs(xi);

      sweepMPP(abs_mu, abs_eta, abs_xi, weight, n, g);
    }
  }

  // Sweep -\mu, -\eta, +\xi.
  {
    const auto current_oct = Octant::MMP;
    for (unsigned int n = 0u; n < _angular_quad.order(current_oct); ++n)
    {
      _angular_quad.direction(current_oct, n, mu, eta, xi);
      weight = _angular_quad.weight(current_oct, n);

      // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
      // if the upwind/downwind directions are aligned properly.
      abs_mu = std::abs(mu);
      abs_eta = std::abs(eta);
      abs_xi = std::abs(xi);

      sweepMMP(abs_mu, abs_eta, abs_xi, weight, n, g);
    }
  }
}

template <typename T>
void
TransportSolver2D<T>::sweepPPP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                               const double & weight, unsigned int ordinate_index, unsigned int g)
{
  // Y second.
  for (unsigned int column = 0u; column < _mesh._tot_num_y; ++column)
  {
    // X first.
    for (unsigned int row = 0u; row < _mesh._tot_num_x; ++row)
    {
      auto & cell = _mesh._cells[column * _mesh._tot_num_x + row];
      _eq_system.solve(cell, weight, abs_mu, abs_eta, abs_xi, ordinate_index, g,
                       CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                       CertesianFaceSide::Back, CertesianFaceSide::Front,
                       _mode, _dt); // z
    }
  }
}

template <typename T>
void
TransportSolver2D<T>::sweepPMP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                               const double & weight, unsigned int ordinate_index, unsigned int g)
{
  unsigned int column;
  // Y second.
  column = _mesh._tot_num_y;
  while (column --> 0)
  {
    // X first.
    for (unsigned int row = 0u; row < _mesh._tot_num_x; ++row)
    {
      auto & cell =_mesh._cells[column * _mesh._tot_num_x + row];
      _eq_system.solve(cell, weight, abs_mu, abs_eta, abs_xi, ordinate_index, g,
                       CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                       CertesianFaceSide::Front, CertesianFaceSide::Back,
                       _mode, _dt); // z
    }
  }
}

template <typename T>
void
TransportSolver2D<T>::sweepMPP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                               const double & weight, unsigned int ordinate_index, unsigned int g)
{
  unsigned int row;
  // Y second.
  for (unsigned int column = 0u; column < _mesh._tot_num_y; ++column)
  {
    // X first.
    row = _mesh._tot_num_x;
    while (row --> 0)
    {
      auto & cell =_mesh._cells[column * _mesh._tot_num_x + row];
      _eq_system.solve(cell, weight, abs_mu, abs_eta, abs_xi, ordinate_index, g,
                       CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                       CertesianFaceSide::Back, CertesianFaceSide::Front,
                       _mode, _dt); // z
    }
  }
}

template <typename T>
void
TransportSolver2D<T>::sweepMMP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                               const double & weight, unsigned int ordinate_index, unsigned int g)
{
  unsigned int column;
  unsigned int row;
  // Y second.
  column = _mesh._tot_num_y;
  while (column --> 0)
  {
    // X first.
    row = _mesh._tot_num_x;
    while (row --> 0)
    {
      auto & cell =_mesh._cells[column * _mesh._tot_num_x + row];
      _eq_system.solve(cell, weight, abs_mu, abs_eta, abs_xi, ordinate_index, g,
                       CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                       CertesianFaceSide::Front, CertesianFaceSide::Back,
                       _mode, _dt); // z
    }
  }
}

template class TransportSolver2D<TWDiamondDifference2D>;
template class TransportSolver2D<DiamondDifference2D>;
