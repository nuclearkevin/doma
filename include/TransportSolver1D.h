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
  TransportSolver1D(BrickMesh1D & mesh, const InputParameters & params, bool verbose, unsigned int num_threads);

  // Solve the subcritical multiplication fixed source problem.
  bool solveFixedSource(const std::string & output_file_base = "", const double & t = 0.0);

  // Solve the fission eigenvalue problem.
  bool solveEigenvalue(const std::string & output_file_base = "");

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

  // Update the external multi-group source (in-scattering and eigenvalue-scaled fission) between power iterations.
  void updateMultigroupSourceEigen(unsigned int g);

  // Update k_{eff} at the end of a power iteration.
  void updateEigenvalue();

  // Solve the within-group equations for the scalar fluxes.
  bool sourceIteration(unsigned int g);

  // Initialize the solver.
  void initializeSolve();

  // Update the scattering source.
  void updateScatteringSource(unsigned int g);

  // Compute the scattering residual to check for source iteration convergence.
  // We use a relative error metric in the L2 integral norm.
  double computeScatteringResidual(unsigned int g);

  // Compute the residual of Gauss-Seidel iteration for the fixed-source flux shape.
  double computeMGFluxResidual();

  // Invert the streaming and collision operator with a sweep.
  void sweep(unsigned int g);

  // Individual sweeping functions for each octant.
  void sweepR(const double & abs_mu, const double & weight, unsigned int ordinate_index,
              unsigned int g, unsigned int tid);
  void sweepL(const double & abs_mu, const double & weight, unsigned int ordinate_index,
              unsigned int g, unsigned int tid);

  // Number of neutron energy groups.
  const unsigned int _num_groups;

  // The execution mode for the solver.
  const RunMode _mode;

  // The mesh to run the transport solver on.
  BrickMesh1D & _mesh;

  // The angular quadrature.
  const GLAngularQuadrature _angular_quad;

  // The templated equation system responsible for solving for cell-centered and interface fluxes.
  const T _eq_system;

  // Whether or not verbose output should be used.
  const bool _verbose;

  // The convergence tolerance for source iteration.
  const double _sit;
  // The maximum number of source iterations.
  const unsigned int _smi;

  // The convergence tolerance for multi-group iteration (Gauss-Seidel / power iteration).
  const double _mgt;
  // The maxmimum number of multi-group iterations.
  const unsigned int _mgi;

  // The criticality eigenvalue k_{eff} at the current power iteration.
  double _k;
  // The criticality eigenvalue k_{eff} at the previous power iteration.
  double _k_prev;
  // The convergence criteria for k_{eff}.
  const double _k_tol;

  // The initial time (if running a transient solve).
  const double _t0;
  // The timestep size to use (if running a transient solve).
  const double _dt;
  // The number of timesteps to solve (if running a transient solve).
  const unsigned int _t_steps;

  // The initial condition type (if running a transient solve).
  const TransientIC _ic;

  // The number of OpenMP threads to use.
  const unsigned int _num_threads;
}; // class TransportSolver

extern template class TransportSolver1D<DiamondDifference1D>;
