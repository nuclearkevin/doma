#pragma once

#include <cmath>

#include "TransportBase.h"
#include "BrickMesh3D.h"
#include "GCAngularQuadrature.h"
#include "BrickCellEquation.h"

// The main class which solves the transport equation using either the upwinded
// diamond difference approximation or the upwinded step approximation.
// Implements a naive sweeper, no fancy wavefront sweeps with KBA (yet).
template <typename T>
class TransportSolver3D
{
public:
  TransportSolver3D(BrickMesh3D & mesh, unsigned int n_l, unsigned int n_c, unsigned int num_threads = 1u)
    : _num_threads(num_threads),
      _mesh(mesh),
      _angular_quad(2u * n_c, 2u * n_l),
      _eq_system()
  {
    static_assert(std::is_base_of<BrickCellEquation, T>::value, "Class must derive from BrickCellEquation.");
  }

  bool solve(const double & source_iteration_tolerance = 1e-5, unsigned int max_iterations = 1000u)
  {
    initializeSolve();

    std::cout << "Solving..." << std::endl;
    unsigned int source_iteration = 0u;
    double current_residual = 0.0;
    do
    {
      std::cout << "Performing source iteration " << source_iteration;
      if (source_iteration != 0u)
        std::cout <<  ", current source iteration residual: " << current_residual << std::endl;
      else
        std::cout << std::endl;

      sweep(source_iteration);
      current_residual = computeScatteringResidual();
      updateScatteringSource();

      source_iteration++;
    }
    while (source_iteration < max_iterations &&  source_iteration_tolerance < current_residual);

    if (source_iteration < max_iterations)
    {
      std::cout << "Scattering source iteration converged after " << source_iteration << " iterations with a residual of " << current_residual << std::endl;
      return true;
    }
    else
    {
      std::cout << "Scattering source iteration failed to convergence after " << source_iteration << " iterations with a residual of " << current_residual << std::endl;
      return false;
    }
  }

private:
  // Initialize the solver.
  void initializeSolve()
  {
    std::cout << "Initializing the solver..." << std::endl;
    for (auto & cell : _mesh._cells)
    {
      cell._total_scalar_flux = 0.0;
      cell._previous_scalar_flux = 0.0;
      cell._current_scalar_flux = 0.0;

      cell._interface_angular_fluxes.fill(0.0);
    }
  }

  // Invert the streaming and collision operator with a sweep.
  void sweep(unsigned int iteration)
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

        sweepPPP(abs_mu, abs_eta, abs_xi, weight, iteration);
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

        sweepPPM(abs_mu, abs_eta, abs_xi, weight, iteration);
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

        sweepPMP(abs_mu, abs_eta, abs_xi, weight, iteration);
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

        sweepPMM(abs_mu, abs_eta, abs_xi, weight, iteration);
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

        sweepMPP(abs_mu, abs_eta, abs_xi, weight, iteration);
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

        sweepMPM(abs_mu, abs_eta, abs_xi, weight, iteration);
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

        sweepMMP(abs_mu, abs_eta, abs_xi, weight, iteration);
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

        sweepMMM(abs_mu, abs_eta, abs_xi, weight, iteration);
      }
    }
  }

  // Update the scattering source.
  void updateScatteringSource()
  {
    for (auto & cell : _mesh._cells)
    {
      cell._total_scalar_flux += cell._current_scalar_flux;
      cell._previous_scalar_flux = cell._current_scalar_flux;
      cell._current_scalar_flux = 0.0;

      cell._interface_angular_fluxes.fill(0.0);
    }
  }

  // Compute the scattering residual to check for source iteration convergence.
  double computeScatteringResidual()
  {
    double residual = 0.0;
    for (auto & cell : _mesh._cells)
      residual += cell._current_scalar_flux * cell._volume;

    return std::abs(residual / _mesh._total_volume);
  }

  // Individual sweeping functions for each octant.
  void sweepPPP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
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
          _eq_system.solve(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                           CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                           CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                           CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
        }
      }
    }
  }

  void sweepPPM(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
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
          _eq_system.solve(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                           CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                           CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                           CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
        }
      }
    }
  }

  void sweepPMP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
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
          _eq_system.solve(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                           CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                           CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                           CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
        }
      }
    }
  }

  void sweepPMM(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
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
          _eq_system.solve(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                           CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                           CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                           CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
        }
      }
    }
  }

  void sweepMPP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
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
          _eq_system.solve(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                           CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                           CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                           CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
        }
      }
    }
  }

  void sweepMPM(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
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
          _eq_system.solve(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                           CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                           CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                           CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
        }
      }
    }
  }

  void sweepMMP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
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
          _eq_system.solve(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                           CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                           CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                           CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
        }
      }
    }
  }

  void sweepMMM(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
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
          _eq_system.solve(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                           CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                           CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                           CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
        }
      }
    }
  }

  const unsigned int _num_threads;

  // The mesh to run the transport solver on.
  BrickMesh3D & _mesh;

  // The 3D angular quadrature.
  const GCAngularQuadrature _angular_quad;

  // The templated equation system responsible for solving for cell-centered and interface fluxes.
  const T _eq_system;
}; // class TransportSolver
