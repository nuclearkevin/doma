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
  TransportSolver3D(BrickMesh3D & mesh, const InputParameters & params, bool verbose);

  // Solve the subcritical multiplication fixed source problem.
  bool solveFixedSource(const std::string & output_file_base = "", const double & t = 0.0);

  // TODO: solve eigenvalue problems.

  // Solve the subcritical transient fixed source problem with delayed neutron precursors.
  bool solveTransient(const std::string & output_file_base = "");

private:
  // Initialize the zero initial condition.
  void initZeroIC();

  // Initialize the steady-state initial condition.
  bool initSteadyIC();

  // Update the fluxes between timesteps.
  void updateStepFluxes();

  // Update the delayed neutron precursors.
  void stepDNPs();

  // Update the external multi-group sources (in-scattering, fission, transient, and external sources)
  // between Gauss-Seidel iterations.
  void updateMultigroupSource(unsigned int g, double t = 0.0);

  // Solve the within-group equations for the scalar fluxes.
  bool sourceIteration(unsigned int g);

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

  // The execution mode for the solver.
  const RunMode _mode;

  // The mesh to run the transport solver on.
  BrickMesh3D & _mesh;

  // The angular quadrature.
  const GCAngularQuadrature _angular_quad;

  // The templated equation system responsible for solving for cell-centered and interface fluxes.
  const T _eq_system;

  // Whether or not verbose output should be used.
  const bool _verbose;

  // The convergence tolerance for source iteration.
  const double _sit;
  // The maximum number of source iterations.
  const unsigned int _smi;

  // The convergence tolerance for multigroup iteration (Gauss-Seidel).
  const double _mgt;
  // The maxmimum number of multigroup iterations.
  const unsigned int _mgi;

  // The initial time (if running a transient solve).
  const double _t0;
  // The timestep size to use (if running a transient solve).
  const double _dt;
  // The number of timesteps to solve (if running a transient solve).
  const unsigned int _t_steps;

  // The initial condition type (if running a transient solve).
  const TransientIC _ic;
}; // class TransportSolver

extern template class TransportSolver3D<TWDiamondDifference3D>;
extern template class TransportSolver3D<DiamondDifference3D>;
