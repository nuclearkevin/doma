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
  TransportSolver3D(BrickMesh3D & mesh, unsigned int num_groups, unsigned int n_l, unsigned int n_c);

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
  void sweepPPP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index, unsigned int g);
  void sweepPPM(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index, unsigned int g);
  void sweepPMP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index, unsigned int g);
  void sweepPMM(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index, unsigned int g);
  void sweepMPP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index, unsigned int g);
  void sweepMPM(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index, unsigned int g);
  void sweepMMP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index, unsigned int g);
  void sweepMMM(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index, unsigned int g);

  // Number of neutron energy groups.
  const unsigned int _num_groups;

  // The mesh to run the transport solver on.
  BrickMesh3D & _mesh;

  // The angular quadrature.
  const GCAngularQuadrature _angular_quad;

  // The templated equation system responsible for solving for cell-centered and interface fluxes.
  const T _eq_system;
}; // class TransportSolver

extern template class TransportSolver3D<TWDiamondDifference3D>;
extern template class TransportSolver3D<DiamondDifference3D>;
