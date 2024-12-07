#include "TransportSolver1D.h"

template <typename T>
TransportSolver1D<T>::TransportSolver1D(BrickMesh1D & mesh, unsigned int num_groups, unsigned int n_l)
  : _num_groups(num_groups),
    _mesh(mesh),
    _angular_quad(2u * n_l),
    _eq_system()
{
  _mesh.initFluxes(_num_groups);
}

template <typename T>
bool
TransportSolver1D<T>::solveFixedSource(const double & sit, unsigned int smi, double mgt, unsigned int mgi)
{
  initializeSolve();

  std::cout << "Solving..." << std::endl;
  unsigned int mg_iteration = 0u;
  double current_residual = 0.0;
  double previous_norm = 0.0;
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
  } while (mg_iteration < mgi && mgt < current_residual);

  if (mg_iteration < mgi)
  {
    std::cout << "MGI converged after " << mg_iteration
              << " iterations with a residual of " << current_residual << std::endl;
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
TransportSolver1D<T>::computeMGFluxNorm()
{
  double norm = 0.0;
  for (auto & cell : _mesh._cells)
    for (unsigned int g = 0u; g < _num_groups; ++g)
      norm += std::pow(cell._total_scalar_flux[g] * cell._l_x, 2.0);

  return std::sqrt(norm);
}

template <typename T>
void
TransportSolver1D<T>::updateMultigroupSource(unsigned int g)
{
  for (auto & cell : _mesh._cells)
  {
    const auto & p = cell.getMatProps();

    cell._current_iteration_source[g] = p._g_src.size() != 0u ? 0.5 * p._g_src[g] : 0.0;
    cell._current_scalar_flux[g] = 0.0;

    for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
    {
      if (g_prime != g)
        cell._current_iteration_source[g] += 0.5 * p._g_g_scatter_mat[g * _num_groups + g_prime] * cell._total_scalar_flux[g_prime];

      if (p._g_chi_p.size() > 0u)
        cell._current_iteration_source[g] += 0.5 * p._g_chi_p[g] * p._g_prod[g_prime] * cell._total_scalar_flux[g_prime];
    }
    cell._total_scalar_flux[g] = 0.0;
  }
}

template <typename T>
bool
TransportSolver1D<T>::sourceIteration(const double & sit, unsigned int smi, unsigned int g)
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
TransportSolver1D<T>::initializeSolve()
{
  std::cout << "Initializing the solver..." << std::endl;

  // Initializing the boundary condition data structure.
  for (unsigned int i = 0u; i < 2u; ++i)
    if (_mesh._bcs[i] != BoundaryCondition::Vacuum)
      _mesh._boundary_angular_fluxes[i].resize(_mesh._boundary_cells[i].size() * _angular_quad.order(), 0.0);

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
TransportSolver1D<T>::updateScatteringSource(unsigned int g)
{
  for (auto & cell : _mesh._cells)
  {
    const auto & p = cell.getMatProps();

    cell._total_scalar_flux[g] += cell._current_scalar_flux[g];
    cell._current_iteration_source[g] = 0.5 * p._g_g_scatter_mat[g * _num_groups + g] * cell._current_scalar_flux[g];
    cell._current_scalar_flux[g] = 0.0;

    cell._interface_angular_fluxes.fill(0.0);
  }
}

template <typename T>
double
TransportSolver1D<T>::computeScatteringResidual(unsigned int g)
{
  double diff_L2 = 0.0;
  double total_L2 = 0.0;
  for (auto & cell : _mesh._cells)
  {
    diff_L2 += std::pow(cell._current_scalar_flux[g], 2.0) * cell._l_x;
    total_L2 += std::pow(cell._total_scalar_flux[g] + cell._current_scalar_flux[g], 2.0) * cell._l_x;
  }

  return total_L2 > 1e-8 ? std::sqrt(diff_L2) / std::sqrt(total_L2) : 0.0;
}

template <typename T>
void
TransportSolver1D<T>::sweep(unsigned int g)
{
  double mu = 0.0;
  double abs_mu = 0.0;
  double weight = 0.0;

  // Sweep +\mu.
  {
    for (unsigned int n = 0u; n < _angular_quad.order() / 2u; ++n)
    {
      mu = _angular_quad.direction(n);
      weight = _angular_quad.weight(n);

      // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
      // if the upwind/downwind directions are aligned properly.
      abs_mu = std::abs(mu);

      sweepR(abs_mu, weight, n, g);
    }
  }

  // Sweep -\mu.
  {
    for (unsigned int n = _angular_quad.order() / 2u; n < _angular_quad.order(); ++n)
    {
      mu = _angular_quad.direction(n);
      weight = _angular_quad.weight(n);

      // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
      // if the upwind/downwind directions are aligned properly.
      abs_mu = std::abs(mu);

      sweepL(abs_mu, weight, n, g);
    }
  }
}

template <typename T>
void
TransportSolver1D<T>::sweepR(const double & abs_mu, const double & weight,
                             unsigned int ordinate_index, unsigned int g)
{
  for (unsigned int row = 0u; row < _mesh._tot_num_x; ++row)
  {
    auto & cell = _mesh._cells[row];
    _eq_system.solve(cell, weight, abs_mu, ordinate_index, g,
                     CertesianFaceSide::Left, CertesianFaceSide::Right);
  }
}

template <typename T>
void
TransportSolver1D<T>::sweepL(const double & abs_mu, const double & weight,
                             unsigned int ordinate_index, unsigned int g)
{
  unsigned int row = _mesh._tot_num_x;
  while (row --> 0)
  {
    auto & cell =_mesh._cells[row];
    _eq_system.solve(cell, weight, abs_mu, ordinate_index, g,
                     CertesianFaceSide::Right, CertesianFaceSide::Left);
  }
}

template class TransportSolver1D<DiamondDifference1D>;
