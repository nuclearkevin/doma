#pragma once

#include "TransportBase.h"
#include "BrickMesh3D.h"
#include "GCAngularQuadrature.h"

// The main class which solves the transport equation using either the upwinded
// diamond difference approximation or the upwinded step approximation.
// Implements a naive sweeper, no fancy wavefront sweeps with KBA (yet).
class TransportSolver
{
public:
  TransportSolver(BrickMesh3D & mesh, unsigned int n_l, unsigned int n_c, DiscretizationType disc_type, unsigned int num_threads = 1u);

  bool solve(const double & source_iteration_tolerance = 1e-5, unsigned int max_iterations = 1000u);

private:
  Octant classifyDirection(const double & mu, const double & eta, const double & xi);

  // Initialize the solver.
  void initializeSolve();
  // Invert the streaming and collision operator with a sweep.
  void sweep(unsigned int iteration);
  // Update the scattering source.
  void updateScatteringSource();
  // Compute the scattering residual to check for source iteration convergence.
  double computeScatteringResidual();

  // Individual sweeping functions for each octant.
  void sweepPPP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration);
  void sweepPPM(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration);
  void sweepPMP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration);
  void sweepPMM(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration);
  void sweepMPP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration);
  void sweepMPM(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration);
  void sweepMMP(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration);
  void sweepMMM(const double & abs_mu, const double & abs_eta, const double & abs_xi, const double & weight, unsigned int iteration);

  // Functions to compute the diamond difference or the step difference spatial discretization scheme.
  void diamondDifference(CartesianCell3D & cell, const double & angular_weight, bool first_iteration,
                         const double & abs_mu, const double & abs_eta, const double & abs_xi,
                         CertesianFaceSide x_uw, CertesianFaceSide x_dw,
                         CertesianFaceSide y_uw, CertesianFaceSide y_dw,
                         CertesianFaceSide z_uw, CertesianFaceSide z_dw);
  void stepCharacteristic(CartesianCell3D & cell, const double & angular_weight, bool first_iteration,
                          const double & abs_mu, const double & abs_eta, const double & abs_xi,
                          CertesianFaceSide x_uw, CertesianFaceSide x_dw,
                          CertesianFaceSide y_uw, CertesianFaceSide y_dw,
                          CertesianFaceSide z_uw, CertesianFaceSide z_dw);

  const DiscretizationType _disc_type;
  const unsigned int _num_threads;

  // The mesh to run the transport solver on.
  BrickMesh3D & _mesh;

  // The 3D angular quadrature.
  const GCAngularQuadrature _angular_quad;
}; // class TransportSolver
