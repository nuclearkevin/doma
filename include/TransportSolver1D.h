#pragma once

#include <cmath>

#include "TransportBase.h"
#include "BrickMesh1D.h"
#include "GLAngularQuadrature.h"

#include "DiamondDifference1D.h"

template <typename T>
class TransportSolver1D;

typedef TransportSolver1D<DiamondDifference1D> DDTransportSolver1D;

// The main class which solves the transport equation using either the upwinded
// diamond difference approximation or the upwinded step approximation.
// Implements a naive sweeper, no fancy wavefront sweeps with KBA (yet).
template <typename T>
class TransportSolver1D
{
public:
  TransportSolver1D(BrickMesh1D & mesh, unsigned int num_groups, unsigned int n_l);

  // Solve the subcritical multiplication fixed source problem.
  bool solveFixedSource(const double & sit, unsigned int smi, double mgt, unsigned int mgi);

private:
  // Update the external multi-group sources (in-scattering, fission, and external sources)
  // between Gauss-Seidel iterations.
  void updateMultigroupSource(unsigned int g);

  // Solve the within-group equations for the scalar fluxes.
  bool sourceIteration(const double & sit, unsigned int smi, unsigned int g);

  // Initialize the solver.
  void initializeSolve();

  // Update the scattering source.
  void updateScatteringSource(unsigned int g);

  // Compute the scattering residual to check for source iteration convergence.
  // We use a relative error metric in the L2 integral norm.
  double computeScatteringResidual(unsigned int g);

  // Compute the L2 norm of the multi-group flux vectors. This is used to assess
  // the convergence of multi-group Gauss-Seidel iteration.
  double computeMGFluxNorm();

  // Invert the streaming and collision operator with a sweep.
  void sweep(unsigned int g);

  // Individual sweeping functions for each octant.
  void sweepR(const double & abs_mu, const double & weight, unsigned int ordinate_index,
              unsigned int g);
  void sweepL(const double & abs_mu, const double & weight, unsigned int ordinate_index,
              unsigned int g);

  // Number of neutron energy groups.
  const unsigned int _num_groups;

  // The mesh to run the transport solver on.
  BrickMesh1D & _mesh;

  // The angular quadrature.
  const GLAngularQuadrature _angular_quad;

  // The templated equation system responsible for solving for cell-centered and interface fluxes.
  const T _eq_system;
}; // class TransportSolver

extern template class TransportSolver1D<DiamondDifference1D>;
