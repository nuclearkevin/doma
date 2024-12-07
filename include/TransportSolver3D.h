#pragma once

#include <cmath>

#include "TransportBase.h"
#include "BrickMesh3D.h"
#include "GCAngularQuadrature.h"

#include "TWDiamondDifference3D.h"
#include "DiamondDifference3D.h"

template <typename T>
class TransportSolver3D;

typedef TransportSolver3D<TWDiamondDifference3D> TWDDTransportSolver3D;
typedef TransportSolver3D<DiamondDifference3D> DDTransportSolver3D;

// The main class which solves the transport equation using either the upwinded
// diamond difference approximation or the upwinded step approximation.
// Implements a naive sweeper, no fancy wavefront sweeps with KBA (yet).
template <typename T>
class TransportSolver3D
{
public:
  TransportSolver3D(BrickMesh3D & mesh, unsigned int n_l, unsigned int n_c);

  bool solveFixedSource(const double & sit, unsigned int smi);

private:
  // Initialize the solver.
  void initializeSolve();

  // Update the scattering source.
  void updateScatteringSource();

  // Compute the scattering residual to check for source iteration convergence.
  // We use a relative error metric in the L2 integral norm.
  double computeScatteringResidual();

  // Invert the streaming and collision operator with a sweep.
  void sweep();

  // Individual sweeping functions for each octant.
  void sweepPPP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index);
  void sweepPPM(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index);
  void sweepPMP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index);
  void sweepPMM(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index);
  void sweepMPP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index);
  void sweepMPM(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index);
  void sweepMMP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index);
  void sweepMMM(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index);

  // The mesh to run the transport solver on.
  BrickMesh3D & _mesh;

  // The angular quadrature.
  const GCAngularQuadrature _angular_quad;

  // The templated equation system responsible for solving for cell-centered and interface fluxes.
  const T _eq_system;
}; // class TransportSolver

extern template class TransportSolver3D<TWDiamondDifference3D>;
extern template class TransportSolver3D<DiamondDifference3D>;
