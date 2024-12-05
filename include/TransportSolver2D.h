#pragma once

#include <cmath>

#include "TransportBase.h"
#include "BrickMesh2D.h"
#include "GCAngularQuadrature.h"
#include "BrickCellEquation2D.h"

// The main class which solves the transport equation using either the upwinded
// diamond difference approximation or the upwinded step approximation.
// Implements a naive sweeper, no fancy wavefront sweeps with KBA (yet).
template <typename T>
class TransportSolver2D
{
public:
  TransportSolver2D(BrickMesh2D & mesh, unsigned int num_groups, unsigned int n_l, unsigned int n_c, unsigned int num_threads = 1u)
    : _num_threads(num_threads),
      _num_groups(num_groups),
      _mesh(mesh),
      _angular_quad(2u * n_c, 2u * n_l, 2u),
      _eq_system()
  {
    static_assert(std::is_base_of<BrickCellEquation2D, T>::value, "Class must derive from BrickCellEquation2D.");
  }

  bool solve(const double & sit, unsigned int smi, double mgt, unsigned int mmi)
  {
    initializeSolve();

    std::cout << "Solving..." << std::endl;
    // Currently the equivalent of solving a lower-triangular matrix.
    // TODO: Gauss-Seidel
    for (unsigned int g = 0u; g < _num_groups; ++g)
    {
      updateMultigroupSource(g);
      auto res = sourceIteration(sit, smi, g);
      if (!res)
        return res;
    }

    return true;
  }

private:
  void updateMultigroupSource(unsigned int g)
  {
    for (auto & cell : _mesh._cells)
    {
      const auto & p = cell.getMatProps();

      cell._total_scalar_flux[g] = 0.0;
      cell._current_iteration_source[g] = 0.5 * p._g_src[g] / M_PI;
      cell._current_scalar_flux[g] = 0.0;

      for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
      {
        if (g_prime != g)
          cell._current_iteration_source[g] += 0.5 * p._g_g_scatter_mat[g * _num_groups + g_prime] * cell._total_scalar_flux[g_prime] / M_PI;
      }
    }
  }

  bool sourceIteration(const double & sit, unsigned int smi, unsigned int g)
  {
    unsigned int source_iteration = 0u;
    double current_residual = 0.0;
    do
    {
      std::cout << "Performing source iteration " << source_iteration << " for group " << g;
      if (source_iteration != 0u)
        std::cout <<  ", current source iteration residual: " << current_residual << std::endl;
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
      std::cout << "Scattering source iteration converged after " << source_iteration
                << " iterations with a residual of " << current_residual << " for group "
                << g << std::endl;
      return true;
    }
    else
    {
      std::cout << "Scattering source iteration failed to convergence after " << source_iteration
                << " iterations with a residual of " << current_residual << " for group "
                << g << std::endl;
      return false;
    }
  }

  // Initialize the solver.
  void initializeSolve()
  {
    std::cout << "Initializing the solver..." << std::endl;

    // Initializing the boundary condition data structure.
    for (unsigned int i = 0u; i < 4u; ++i)
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

  // Update the angular flux boundary conditions.
  void updateBoundaryAngularFluxes(unsigned int ordinate_index, Octant current_oct,
                                   CertesianFaceSide x_uw, CertesianFaceSide x_dw,
                                   CertesianFaceSide y_uw, CertesianFaceSide y_dw)
  {
    switch (_mesh._bcs[static_cast<unsigned int>(x_dw)])
    {
      case BoundaryCondition::Reflective:
      {

        break;
      }
      default: break;
    }

    switch (_mesh._bcs[static_cast<unsigned int>(y_dw)])
    {
      case BoundaryCondition::Reflective:
      {

        break;
      }
      default: break;
    }
  }

  // Update the scattering source.
  void updateScatteringSource(unsigned int g)
  {
    for (auto & cell : _mesh._cells)
    {
      const auto & p = cell.getMatProps();

      cell._total_scalar_flux[g] += cell._current_scalar_flux[g];
      cell._current_iteration_source[g] = 0.5 * p._g_g_scatter_mat[g * _num_groups + g] * cell._current_scalar_flux[g] / M_PI;
      cell._current_scalar_flux[g] = 0.0;

      cell._interface_angular_fluxes.fill(0.0);
    }
  }

  // Compute the scattering residual to check for source iteration convergence.
  // We use a relative error metric in the L2 integral norm.
  double computeScatteringResidual(unsigned int g)
  {
    double diff_L2 = 0.0;
    double total_L2 = 0.0;
    for (auto & cell : _mesh._cells)
    {
      diff_L2 += std::pow(cell._current_scalar_flux[g], 2.0) * cell._area;
      total_L2 += std::pow(cell._total_scalar_flux[g] + cell._current_scalar_flux[g], 2.0) * cell._area;
    }

    return total_L2 > 1e-8 ? std::sqrt(diff_L2) / std::sqrt(total_L2) : 0.0;
  }

  // Invert the streaming and collision operator with a sweep.
  void sweep(unsigned int g)
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
        updateBoundaryAngularFluxes(n, current_oct,
                                    CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                                    CertesianFaceSide::Back, CertesianFaceSide::Front); // z
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
        updateBoundaryAngularFluxes(n, current_oct,
                                    CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                                    CertesianFaceSide::Front, CertesianFaceSide::Back); // z
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
        updateBoundaryAngularFluxes(n, current_oct,
                                    CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                                    CertesianFaceSide::Back, CertesianFaceSide::Front); // z
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
        updateBoundaryAngularFluxes(n, current_oct,
                                    CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                                    CertesianFaceSide::Front, CertesianFaceSide::Back); // z
      }
    }
  }

  // Individual sweeping functions for each octant.
  void sweepPPP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int ordinate_index, unsigned int g)
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
                         CertesianFaceSide::Back, CertesianFaceSide::Front); // z
      }
    }
  }

  void sweepPMP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int ordinate_index, unsigned int g)
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
                         CertesianFaceSide::Front, CertesianFaceSide::Back); // z
      }
    }
  }

  void sweepMPP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int ordinate_index, unsigned int g)
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
                         CertesianFaceSide::Back, CertesianFaceSide::Front); // z
      }
    }
  }

  void sweepMMP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int ordinate_index, unsigned int g)
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
                         CertesianFaceSide::Front, CertesianFaceSide::Back); // z
      }
    }
  }

  const unsigned int _num_threads;

  const unsigned int _num_groups;

  // The mesh to run the transport solver on.
  BrickMesh2D & _mesh;

  // The angular quadrature.
  const GCAngularQuadrature _angular_quad;

  // The templated equation system responsible for solving for cell-centered and interface fluxes.
  const T _eq_system;
}; // class TransportSolver
