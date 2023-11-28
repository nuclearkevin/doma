#include "TransportSolver.h"

#include <cmath>

TransportSolver::TransportSolver(BrickMesh3D & mesh, unsigned int n_l, unsigned int n_c, DiscretizationType disc_type)
  : _disc_type(disc_type),
    _mesh(mesh),
    _angular_quad(2u * n_c, 2u * n_l)
{ }

bool
TransportSolver::solve(const double & source_iteration_tolerance, unsigned int max_iterations)
{
  initializeSolve();

  std::cout << "Solving..." << std::endl;
  unsigned int iteration = 0u;
  double current_residual = 0.0;
  do
  {
    std::cout << "Performing source iteration " << iteration;
    if (iteration != 0u)
      std::cout <<  ", current source iteration residual: " << current_residual << std::endl;
    else
      std::cout << std::endl;

    sweep(iteration);
    current_residual = computeScatteringResidual();
    updateScatteringSource();

    iteration++;
  }
  while (iteration < max_iterations &&  source_iteration_tolerance < current_residual);

  if (iteration < max_iterations)
  {
    std::cout << "Scattering source iteration converged after " << iteration << " iterations with a residual of " << current_residual << std::endl;
    return true;
  }
  else
  {
    std::cout << "Scattering source iteration failed to convergence after " << iteration << " iterations with a residual of " << current_residual << std::endl;
    return false;
  }
}

Octant
TransportSolver::classifyDirection(const double & mu, const double & eta, const double & xi)
{
  if (mu > 0.0 && eta > 0.0 && xi > 0.0)
    return Octant::PPP;
  if (mu > 0.0 && eta > 0.0 && xi < 0.0)
    return Octant::PPM;
  if (mu > 0.0 && eta < 0.0 && xi > 0.0)
    return Octant::PMP;
  if (mu > 0.0 && eta < 0.0 && xi < 0.0)
    return Octant::PMM;
  if (mu < 0.0 && eta > 0.0 && xi > 0.0)
    return Octant::MPP;
  if (mu < 0.0 && eta > 0.0 && xi < 0.0)
    return Octant::MPM;
  if (mu < 0.0 && eta < 0.0 && xi > 0.0)
    return Octant::MMP;
  if (mu < 0.0 && eta < 0.0 && xi < 0.0)
    return Octant::MMM;

  std::cout << "Error: TransportSolver::classifyDirection(const double & mu, const double & eta, const double & xi)." << std::endl;
  exit(1);
  return Octant::PPP;
}

// Initialize the solver.
void
TransportSolver::initializeSolve()
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
void
TransportSolver::sweep(unsigned int iteration)
{
  double mu = 0.0;
  double eta = 0.0;
  double xi = 0.0;

  double abs_mu = 0.0;
  double abs_eta = 0.0;
  double abs_xi = 0.0;

  double weight = 0.0;

  for (unsigned int n = 0u; n < _angular_quad.totalOrder(); ++n)
  {
    _angular_quad.direction(n, mu, eta, xi);

    //std::cout << "Sweeping direction " << n << " (" << mu << ", " << eta << ", " << xi << ")." << std::endl;

    // We use the absolute value of each ordinate as the system of equations ends up being symmetrical
    // if the upwind/downwind directions are aligned properly.
    abs_mu = std::abs(mu);
    abs_eta = std::abs(eta);
    abs_xi = std::abs(xi);

    weight = _angular_quad.weight(n);

    switch (classifyDirection(mu, eta, xi))
    {
      case Octant::PPP: sweepPPP(abs_mu, abs_eta, abs_xi, weight, iteration); break; // Left->Right, Back->Front, Bottom->Top
      case Octant::PPM: sweepPPM(abs_mu, abs_eta, abs_xi, weight, iteration); break; // Left->Right, Back->Front, Top->Bottom
      case Octant::PMP: sweepPMP(abs_mu, abs_eta, abs_xi, weight, iteration); break; // Left->Right, Front->Back, Bottom->Top
      case Octant::PMM: sweepPMM(abs_mu, abs_eta, abs_xi, weight, iteration); break; // Left->Right, Front->Back, Top->Bottom
      case Octant::MPP: sweepMPP(abs_mu, abs_eta, abs_xi, weight, iteration); break; // Right->Left, Back->Front, Bottom->Top
      case Octant::MPM: sweepMPM(abs_mu, abs_eta, abs_xi, weight, iteration); break; // Right->Left, Back->Front, Top->Bottom
      case Octant::MMP: sweepMMP(abs_mu, abs_eta, abs_xi, weight, iteration); break; // Right->Left, Front->Back, Bottom->Top
      case Octant::MMM: sweepMMM(abs_mu, abs_eta, abs_xi, weight, iteration); break; // Right->Left, Front->Back, Top->Bottom
    }
  }
}

// Update the scattering source.
void
TransportSolver::updateScatteringSource()
{
  for (auto & cell : _mesh._cells)
  {
    cell._total_scalar_flux += cell._current_scalar_flux;
    cell._previous_scalar_flux = cell._current_scalar_flux;
    cell._current_scalar_flux = 0.0;

    cell._interface_angular_fluxes.fill(0.0);
  }
}

// Check for convergence.
double
TransportSolver::computeScatteringResidual()
{
  double residual = 0.0;
  for (auto & cell : _mesh._cells)
    residual += cell._current_scalar_flux * cell._volume;

  return std::abs(residual / _mesh._total_volume);
}

// Individual sweeping functions for each octant.
// Left->Right, Back->Front, Bottom->Top
void
TransportSolver::sweepPPP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
{
  switch (_disc_type)
  {
    case DiscretizationType::DiamondDifference:
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
            diamondDifference(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                              CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                              CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                              CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
          }
        }
      }
      break;

    case DiscretizationType::StepCharacteristics:
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
            stepCharacteristic(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                               CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                               CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                               CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
          }
        }
      }
      break;
  }
}

// Left->Right, Back->Front, Top->Bottom
void
TransportSolver::sweepPPM(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
{
  unsigned int slice;
  switch (_disc_type)
  {
    case DiscretizationType::DiamondDifference:
      // Z third.
      slice = _mesh._tot_num_z;
      while (slice --> 0)
      {
        // Y second.
        for (unsigned int column = 0u; column < _mesh._tot_num_y; ++column)
        {
          // X first.
          for (unsigned int row = 0u; row < _mesh._tot_num_x; ++row)
          {
            auto & cell =_mesh._cells[slice * _mesh._tot_num_y * _mesh._tot_num_x + column * _mesh._tot_num_x + row];
            diamondDifference(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                              CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                              CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                              CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
          }
        }
      }
      break;

    case DiscretizationType::StepCharacteristics:
      // Z third.
      slice = _mesh._tot_num_z;
      while (slice --> 0)
      {
        // Y second.
        for (unsigned int column = 0u; column < _mesh._tot_num_y; ++column)
        {
          // X first.
          for (unsigned int row = 0u; row < _mesh._tot_num_x; ++row)
          {
            auto & cell =_mesh._cells[slice * _mesh._tot_num_y * _mesh._tot_num_x + column * _mesh._tot_num_x + row];
            stepCharacteristic(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                               CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                               CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                               CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
          }
        }
      }
      break;
  }
}

// Left->Right, Front->Back, Bottom->Top
void
TransportSolver::sweepPMP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
{
  unsigned int column;
  switch (_disc_type)
  {
    case DiscretizationType::DiamondDifference:
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
            diamondDifference(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                              CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                              CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                              CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
          }
        }
      }
      break;

    case DiscretizationType::StepCharacteristics:
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
            stepCharacteristic(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                               CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                               CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                               CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
          }
        }
      }
      break;
  }
}

// Left->Right, Front->Back, Top->Bottom
void
TransportSolver::sweepPMM(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
{
  unsigned int slice;
  unsigned int column;
  switch (_disc_type)
  {
    case DiscretizationType::DiamondDifference:
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
            diamondDifference(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                              CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                              CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                              CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
          }
        }
      }
      break;

    case DiscretizationType::StepCharacteristics:
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
            stepCharacteristic(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                               CertesianFaceSide::Left, CertesianFaceSide::Right,  // x
                               CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                               CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
          }
        }
      }
      break;
  }
}

// Right->Left, Back->Front, Bottom->Top
void
TransportSolver::sweepMPP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
{
  unsigned int row;
  switch (_disc_type)
  {
    case DiscretizationType::DiamondDifference:
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
            diamondDifference(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                              CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                              CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                              CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
          }
        }
      }
      break;

    case DiscretizationType::StepCharacteristics:
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
            stepCharacteristic(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                               CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                               CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                               CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
          }
        }
      }
      break;
  }
}

// Right->Left, Back->Front, Top->Bottom
void
TransportSolver::sweepMPM(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
{
  unsigned int row;
  unsigned int slice;
  switch (_disc_type)
  {
    case DiscretizationType::DiamondDifference:
      // Z third.
      slice = _mesh._tot_num_z;
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
            diamondDifference(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                              CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                              CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                              CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
          }
        }
      }
      break;

    case DiscretizationType::StepCharacteristics:
      // Z third.
      slice = _mesh._tot_num_z;
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
            stepCharacteristic(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                               CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                               CertesianFaceSide::Back, CertesianFaceSide::Front,  // y
                               CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
          }
        }
      }
      break;
  }
}

// Right->Left, Front->Back, Bottom->Top
void
TransportSolver::sweepMMP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
{
  unsigned int column;
  unsigned int row;
  switch (_disc_type)
  {
    case DiscretizationType::DiamondDifference:
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
            diamondDifference(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                              CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                              CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                              CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
          }
        }
      }
      break;

    case DiscretizationType::StepCharacteristics:
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
            stepCharacteristic(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                               CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                               CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                               CertesianFaceSide::Bottom, CertesianFaceSide::Top); // z
          }
        }
      }
      break;
  }
}

// Right->Left, Front->Back, Top->Bottom
void
TransportSolver::sweepMMM(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration)
{
  unsigned int slice;
  unsigned int column;
  unsigned int row;
  switch (_disc_type)
  {
    case DiscretizationType::DiamondDifference:
      // Z third.
      slice = _mesh._tot_num_z;
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
            diamondDifference(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                              CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                              CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                              CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
          }
        }
      }
      break;

    case DiscretizationType::StepCharacteristics:
      // Z third.
      slice = _mesh._tot_num_z;
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
            stepCharacteristic(cell, weight, iteration == 0u, abs_mu, abs_eta, abs_xi,
                               CertesianFaceSide::Right, CertesianFaceSide::Left,  // x
                               CertesianFaceSide::Front, CertesianFaceSide::Back,  // y
                               CertesianFaceSide::Top, CertesianFaceSide::Bottom); // z
          }
        }
      }
      break;
  }
}

// Functions to compute the diamond difference or the step difference spatial discretization scheme.
void
TransportSolver::diamondDifference(CartesianCell3D & cell, const double & angular_weight, bool first_iteration,
                                   const double & abs_mu, const double & abs_eta, const double & abs_xi,
                                   CertesianFaceSide x_uw, CertesianFaceSide x_dw,
                                   CertesianFaceSide y_uw, CertesianFaceSide y_dw,
                                   CertesianFaceSide z_uw, CertesianFaceSide z_dw)
{
  // Apply the fixed source on the first iteration.
  double cell_source = first_iteration ? cell._fixed_src / (4.0 * M_PI) : 0.0;
  // Add the scattering source.
  cell_source += cell._sigma_s / (4.0 * M_PI) * cell._previous_scalar_flux;

  // Grab the upwind interfacing angular fluxes. The conditions handle the vacuum boundary conditions.
  const auto x_uw_af = cell.neighbor(x_uw) ? cell.neighbor(x_uw)->interfaceFlux(x_dw) : 0.0;
  const auto y_uw_af = cell.neighbor(y_uw) ? cell.neighbor(y_uw)->interfaceFlux(y_dw) : 0.0;
  const auto z_uw_af = cell.neighbor(z_uw) ? cell.neighbor(z_uw)->interfaceFlux(z_dw) : 0.0;

  // Diamond difference approximation.
  // Computing cell-centered angular fluxes.
  double center_af = cell_source + (2.0 * x_uw_af * abs_mu / cell._l_x) + (2.0 * y_uw_af * abs_eta / cell._l_y) + (2.0 * z_uw_af * abs_xi / cell._l_z);
  center_af /= (2.0 * abs_mu / cell._l_x) + (2.0 * abs_eta / cell._l_x) + (2.0 * abs_xi / cell._l_x) + cell._sigma_t;

  // Add this angular flux's contribution to the scalar flux.
  cell._current_scalar_flux += angular_weight * center_af;

  // Update the interface angular fluxes using the diamond difference closures. 2nd order accurate!
  cell.setInterfaceFlux(x_dw, 2.0 * center_af - x_uw_af);
  cell.setInterfaceFlux(y_dw, 2.0 * center_af - y_uw_af);
  cell.setInterfaceFlux(z_dw, 2.0 * center_af - z_uw_af);
}

void
TransportSolver::stepCharacteristic(CartesianCell3D & cell, const double & angular_weight, bool first_iteration,
                                    const double & abs_mu, const double & abs_eta, const double & abs_xi,
                                    CertesianFaceSide x_uw, CertesianFaceSide x_dw,
                                    CertesianFaceSide y_uw, CertesianFaceSide y_dw,
                                    CertesianFaceSide z_uw, CertesianFaceSide z_dw)
{
  // Apply the fixed source on the first iteration.
  double cell_source = first_iteration ? cell._fixed_src / (4.0 * M_PI) : 0.0;
  // Add the scattering source.
  cell_source += cell._sigma_s / (4.0 * M_PI) * cell._previous_scalar_flux;

  // Grab the upwind interfacing angular fluxes. The conditions handle the vacuum boundary conditions.
  const auto x_uw_af = cell.neighbor(x_uw) ? cell.neighbor(x_uw)->interfaceFlux(x_dw) : 0.0;
  const auto y_uw_af = cell.neighbor(y_uw) ? cell.neighbor(y_uw)->interfaceFlux(y_dw) : 0.0;
  const auto z_uw_af = cell.neighbor(z_uw) ? cell.neighbor(z_uw)->interfaceFlux(z_dw) : 0.0;

  // Step characteristics approximation.
  // Computing cell-centered angular fluxes.
  double center_af = cell_source + (x_uw_af * abs_mu / cell._l_x) + (y_uw_af * abs_eta / cell._l_y) + (z_uw_af * abs_xi / cell._l_z);
  center_af /= (abs_mu / cell._l_x) + (abs_eta / cell._l_x) + (abs_xi / cell._l_x) + cell._sigma_t;

  // Add this angular flux's contribution to the scalar flux.
  cell._current_scalar_flux += angular_weight * center_af;

  // Update the interface angular fluxes using the step closures. 1st order accurate!
  cell.setInterfaceFlux(x_dw, center_af);
  cell.setInterfaceFlux(y_dw, center_af);
  cell.setInterfaceFlux(z_dw, center_af);
}
