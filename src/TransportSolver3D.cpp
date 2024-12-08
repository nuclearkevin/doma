#include "TransportSolver3D.h"

template <typename T>
TransportSolver3D<T>::TransportSolver3D(BrickMesh3D & mesh, unsigned int num_groups, unsigned int n_l, unsigned int n_c)
  : _num_groups(num_groups),
    _mesh(mesh),
    _angular_quad(2u * n_c, 2u * n_l, 3u),
    _eq_system()
{
  _mesh.initFluxes(_num_groups);
}

template <typename T>
bool
TransportSolver3D<T>::solveFixedSource(const double & sit, unsigned int smi, double mgt, unsigned int mgi)
{
  initializeSolve();

  std::cout << "Solving..." << std::endl;
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
      updateMultigroupSource(g);
      auto res = sourceIteration(sit, smi, g);
      if (!res)
        return res;
    }

    current_norm = computeMGFluxNorm();
    current_residual = std::abs(current_norm - previous_norm) / previous_norm;
    previous_norm = current_norm;

    mg_iteration++;
  } while (mg_iteration < mgi && mgt < current_residual && _num_groups > 1u);

  if (mg_iteration < mgi)
  {
    if (_num_groups > 1u)
    {
      std::cout << "MGI converged after " << mg_iteration
                << " iterations with a residual of " << current_residual << std::endl;
    }
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
TransportSolver3D<T>::computeMGFluxNorm()
{
  double norm = 0.0;
  for (auto & cell : _mesh._cells)
    for (unsigned int g = 0u; g < _num_groups; ++g)
      norm += std::pow(cell._total_scalar_flux[g] * cell._volume, 2.0);

  return std::sqrt(norm);
}

template <typename T>
void
TransportSolver3D<T>::updateMultigroupSource(unsigned int g)
{
  for (auto & cell : _mesh._cells)
  {
    const auto & p = cell.getMatProps();

    cell._current_iteration_source[g] = p._g_src.size() != 0u ? 0.25 * p._g_src[g] / M_PI : 0.0;
    cell._current_scalar_flux[g] = 0.0;

    for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
    {
      if (g_prime != g)
        cell._current_iteration_source[g] += 0.25 * p._g_g_scatter_mat[g * _num_groups + g_prime] * cell._total_scalar_flux[g_prime] / M_PI;

      if (p._g_chi_p.size() > 0u)
        cell._current_iteration_source[g] += 0.25 * p._g_chi_p[g] * p._g_prod[g_prime] * cell._total_scalar_flux[g_prime] / M_PI;
    }
    cell._total_scalar_flux[g] = 0.0;
  }
}

template <typename T>
bool
TransportSolver3D<T>::sourceIteration(const double & sit, unsigned int smi, unsigned int g)
{
  unsigned int source_iteration = 0u;
  double current_residual = 0.0;
  do
  {
    std::cout << "Performing SI " << source_iteration << " for G" << g;
    if (source_iteration != 0u)
      std::cout <<  ", residual = " << current_residual << std::endl;
    else
      std::cout << std::endl;

    sweep(g);
    current_residual = computeScatteringResidual(g);
    updateScatteringSource(g);

    source_iteration++;
  }
  while (source_iteration < smi &&  sit < current_residual);

  if (source_iteration < smi)
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
TransportSolver3D<T>::initializeSolve()
{
  std::cout << "Initializing the solver..." << std::endl;

  // Initializing the boundary condition data structure.
  for (unsigned int i = 0u; i < 8u; ++i)
    if (_mesh._bcs[i] != BoundaryCondition::Vacuum)
      _mesh._boundary_angular_fluxes[i].resize(_mesh._boundary_cells[i].size() * _angular_quad.totalOrder(), 0.0);

  for (auto & cell : _mesh._cells)
  {
    const auto & p = cell.getMatProps();
    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      cell._total_scalar_flux[g] = 0.0;
      cell._current_iteration_source[g] = 0.0;
      cell._current_scalar_flux[g] = 0.0;
    }
    cell._interface_angular_fluxes.fill(0.0);
  }
}

template <typename T>
void
TransportSolver3D<T>::updateScatteringSource(unsigned int g)
{
  for (auto & cell : _mesh._cells)
  {
    const auto & p = cell.getMatProps();

    cell._total_scalar_flux[g] += cell._current_scalar_flux[g];
    cell._current_iteration_source[g] = 0.25 * p._g_g_scatter_mat[g * _num_groups + g] * cell._current_scalar_flux[g] / M_PI;
    cell._current_scalar_flux[g] = 0.0;

    cell._interface_angular_fluxes.fill(0.0);
  }
}

template <typename T>
double
TransportSolver3D<T>::computeScatteringResidual(unsigned int g)
{
  double diff_L2 = 0.0;
  double total_L2 = 0.0;
  for (auto & cell : _mesh._cells)
  {
    diff_L2 += std::pow(cell._current_scalar_flux[g], 2.0) * cell._volume;
    total_L2 += std::pow(cell._total_scalar_flux[g] + cell._current_scalar_flux[g], 2.0) * cell._volume;
  }

  return total_L2 > 1e-8 ? std::sqrt(diff_L2) / std::sqrt(total_L2) : 0.0;
}

template <typename T>
void
TransportSolver3D<T>::sweep(unsigned int g)
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

  // Sweep +\mu, +\eta, -\xi.
  {
    const auto current_oct = Octant::PPM;
    for (unsigned int n = 0u; n < _angular_quad.order(current_oct); ++n)
    {
      _angular_quad.direction(current_oct, n, mu, eta, xi);
      weight = _angular_quad.weight(current_oct, n);

      // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
      // if the upwind/downwind directions are aligned properly.
      abs_mu = std::abs(mu);
      abs_eta = std::abs(eta);
      abs_xi = std::abs(xi);

      sweepPPM(abs_mu, abs_eta, abs_xi, weight, n, g);
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

  // Sweep +\mu, -\eta, -\xi.
  {
    const auto current_oct = Octant::PMM;
    for (unsigned int n = 0u; n < _angular_quad.order(current_oct); ++n)
    {
      _angular_quad.direction(current_oct, n, mu, eta, xi);
      weight = _angular_quad.weight(current_oct, n);

      // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
      // if the upwind/downwind directions are aligned properly.
      abs_mu = std::abs(mu);
      abs_eta = std::abs(eta);
      abs_xi = std::abs(xi);

      sweepPMM(abs_mu, abs_eta, abs_xi, weight, n, g);
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

  // Sweep -\mu, +\eta, -\xi.
  {
    const auto current_oct = Octant::MPM;
    for (unsigned int n = 0u; n < _angular_quad.order(current_oct); ++n)
    {
      _angular_quad.direction(current_oct, n, mu, eta, xi);
      weight = _angular_quad.weight(current_oct, n);

      // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
      // if the upwind/downwind directions are aligned properly.
      abs_mu = std::abs(mu);
      abs_eta = std::abs(eta);
      abs_xi = std::abs(xi);

      sweepMPM(abs_mu, abs_eta, abs_xi, weight, n, g);
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

  // Sweep -\mu, -\eta, -\xi.
  {
    const auto current_oct = Octant::MMM;
    for (unsigned int n = 0u; n < _angular_quad.order(current_oct); ++n)
    {
      _angular_quad.direction(current_oct, n, mu, eta, xi);
      weight = _angular_quad.weight(current_oct, n);

      // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
      // if the upwind/downwind directions are aligned properly.
      abs_mu = std::abs(mu);
      abs_eta = std::abs(eta);
      abs_xi = std::abs(xi);

      sweepMMM(abs_mu, abs_eta, abs_xi, weight, n, g);
    }
  }
}

template <typename T>
void
TransportSolver3D<T>::sweepPPP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                               const double & weight, unsigned int ordinate_index, unsigned int g)
{
  // Z third.
  for (unsigned int slice = 0u; slice < _mesh._tot_num_z; ++slice)
  {
    // Y second.
    for (unsigned int column = 0u; column < _mesh._tot_num_y; ++column)
    {
      // X first.
      for (unsigned int row = 0u; row < _mesh._tot_num_x; ++row)
      {
        auto & cell = _mesh._cells[slice * _mesh._tot_num_y * _mesh._tot_num_x + column * _mesh._tot_num_x + row];
        _eq_system.solve(cell, weight, abs_mu, abs_eta, abs_xi, ordinate_index, g,
                         CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                         CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                         CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
      }
    }
  }
}

template <typename T>
void
TransportSolver3D<T>::sweepPPM(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                               const double & weight, unsigned int ordinate_index, unsigned int g)
{
  // Z third.
  unsigned int slice = _mesh._tot_num_z;
  while (slice --> 0)
  {
    // Y second.
    for (unsigned int column = 0u; column < _mesh._tot_num_y; ++column)
    {
      // X first.
      for (unsigned int row = 0u; row < _mesh._tot_num_x; ++row)
      {
        auto & cell =_mesh._cells[slice * _mesh._tot_num_y * _mesh._tot_num_x + column * _mesh._tot_num_x + row];
        _eq_system.solve(cell, weight, abs_mu, abs_eta, abs_xi, ordinate_index, g,
                         CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                         CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                         CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
      }
    }
  }
}

template <typename T>
void
TransportSolver3D<T>::sweepPMP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                               const double & weight, unsigned int ordinate_index, unsigned int g)
{
  unsigned int column;
  // Z third.
  for (unsigned int slice = 0u; slice < _mesh._tot_num_z; ++slice)
  {
    // Y second.
    column = _mesh._tot_num_y;
    while (column --> 0)
    {
      // X first.
      for (unsigned int row = 0u; row < _mesh._tot_num_x; ++row)
      {
        auto & cell =_mesh._cells[slice * _mesh._tot_num_y * _mesh._tot_num_x + column * _mesh._tot_num_x + row];
        _eq_system.solve(cell, weight, abs_mu, abs_eta, abs_xi, ordinate_index, g,
                         CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                         CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                         CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
      }
    }
  }
}

template <typename T>
void
TransportSolver3D<T>::sweepPMM(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                               const double & weight, unsigned int ordinate_index, unsigned int g)
{
  unsigned int slice;
  unsigned int column;
  // Z third.
  slice = _mesh._tot_num_z;
  while (slice --> 0)
  {
    // Y second.
    column = _mesh._tot_num_y;
    while (column --> 0)
    {
      // X first.
      for (unsigned int row = 0u; row < _mesh._tot_num_x; ++row)
      {
        auto & cell =_mesh._cells[slice * _mesh._tot_num_y * _mesh._tot_num_x + column * _mesh._tot_num_x + row];
        _eq_system.solve(cell, weight, abs_mu, abs_eta, abs_xi, ordinate_index, g,
                         CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                         CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                         CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
      }
    }
  }
}

template <typename T>
void
TransportSolver3D<T>::sweepMPP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                               const double & weight, unsigned int ordinate_index, unsigned int g)
{
  unsigned int row;
  // Z third.
  for (unsigned int slice = 0u; slice < _mesh._tot_num_z; ++slice)
  {
    // Y second.
    for (unsigned int column = 0u; column < _mesh._tot_num_y; ++column)
    {
      // X first.
      row = _mesh._tot_num_x;
      while (row --> 0)
      {
        auto & cell =_mesh._cells[slice * _mesh._tot_num_y * _mesh._tot_num_x + column * _mesh._tot_num_x + row];
        _eq_system.solve(cell, weight, abs_mu, abs_eta, abs_xi, ordinate_index, g,
                         CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                         CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                         CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
      }
    }
  }
}

template <typename T>
void
TransportSolver3D<T>::sweepMPM(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                               const double & weight, unsigned int ordinate_index, unsigned int g)
{
  unsigned int row;
  // Z third.
  unsigned int slice = _mesh._tot_num_z;
  while (slice --> 0)
  {
    // Y second.
    for (unsigned int column = 0u; column < _mesh._tot_num_y; ++column)
    {
      // X first.
      row = _mesh._tot_num_x;
      while (row --> 0)
      {
        auto & cell =_mesh._cells[slice * _mesh._tot_num_y * _mesh._tot_num_x + column * _mesh._tot_num_x + row];
        _eq_system.solve(cell, weight, abs_mu, abs_eta, abs_xi, ordinate_index, g,
                         CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                         CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                         CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
      }
    }
  }
}

template <typename T>
void
TransportSolver3D<T>::sweepMMP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                               const double & weight, unsigned int ordinate_index, unsigned int g)
{
  unsigned int column;
  unsigned int row;
  // Z third.
  for (unsigned int slice = 0u; slice < _mesh._tot_num_z; ++slice)
  {
    // Y second.
    column = _mesh._tot_num_y;
    while (column --> 0)
    {
      // X first.
      row = _mesh._tot_num_x;
      while (row --> 0)
      {
        auto & cell =_mesh._cells[slice * _mesh._tot_num_y * _mesh._tot_num_x + column * _mesh._tot_num_x + row];
        _eq_system.solve(cell, weight, abs_mu, abs_eta, abs_xi, ordinate_index, g,
                         CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                         CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                         CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
      }
    }
  }
}

template <typename T>
void
TransportSolver3D<T>::sweepMMM(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                               const double & weight, unsigned int ordinate_index, unsigned int g)
{
  unsigned int column;
  unsigned int row;
  // Z third.
  unsigned int slice = _mesh._tot_num_z;
  while (slice --> 0)
  {
    // Y second.
    column = _mesh._tot_num_y;
    while (column --> 0)
    {
      // X first.
      row = _mesh._tot_num_x;
      while (row --> 0)
      {
        auto & cell =_mesh._cells[slice * _mesh._tot_num_y * _mesh._tot_num_x + column * _mesh._tot_num_x + row];
        _eq_system.solve(cell, weight, abs_mu, abs_eta, abs_xi, ordinate_index, g,
                         CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                         CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                         CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
      }
    }
  }
}

template class TransportSolver3D<TWDiamondDifference3D>;
template class TransportSolver3D<DiamondDifference3D>;
