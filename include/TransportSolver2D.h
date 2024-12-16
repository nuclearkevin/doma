#pragma once

#include <cmath>

#include "TransportBase.h"
#include "BrickMesh2D.h"
#include "GCAngularQuadrature.h"

#include "TWDiamondDifference2D.h"
#include "DiamondDifference2D.h"

template <typename T>
class TransportSolver2D;

typedef TransportSolver2D<TWDiamondDifference2D> TWDDTransportSolver2D;
typedef TransportSolver2D<DiamondDifference2D> DDTransportSolver2D;

// The main class which solves the transport equation using either the upwinded
// diamond difference approximation or the upwinded step approximation.
// Implements a naive sweeper, no fancy wavefront sweeps with KBA (yet).
template <typename T>
class TransportSolver2D
{
public:
  TransportSolver2D(BrickMesh2D & mesh, const RunMode & mode, bool verbose, unsigned int num_groups, unsigned int n_l, unsigned int n_c);

  // Solve the subcritical multiplication fixed source problem.
  bool solveFixedSource(const double & sit, unsigned int smi, double mgt, unsigned int mgi, double t = 0.0);

  // Solve the subcritical transient fixed source problem with delayed neutron precursors.
  bool solveTransient(const double & t0, const double & dt, unsigned int t_steps, const TransientIC & ic,
                      const double & sit, unsigned int smi, double mgt, unsigned int mgi,
                      const std::string & output_file_base = "");

private:
  // Initialize the zero initial condition.
  void initZeroIC();

  // Initialize the steady-state initial condition.
  bool initSteadyIC(const double & sit, unsigned int smi, double mgt, unsigned int mgi);

  // Update the fluxes between timesteps.
  void updateStepFluxes();

  // Update the delayed neutron precursors.
  void stepDNPs();

  // Update the external multi-group sources (in-scattering, fission, and external sources)
  // between Gauss-Seidel iterations.
  void updateMultigroupSource(unsigned int g, double t = 0.0);

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
  void sweepPMP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index, unsigned int g);
  void sweepMPP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index, unsigned int g);
  void sweepMMP(const double & abs_mu, const double & abs_eta, const double & abs_xi,
                const double & weight, unsigned int ordinate_index, unsigned int g);

  // Number of neutron energy groups.
  const unsigned int _num_groups;

  // The execution mode for the solver.
  const RunMode _mode;

  // The dt to use (if running a transient solve).
  double _dt;

  // The mesh to run the transport solver on.
  BrickMesh2D & _mesh;

  // The angular quadrature.
  const GCAngularQuadrature _angular_quad;

  // The templated equation system responsible for solving for cell-centered and interface fluxes.
  const T _eq_system;

  // Whether or not verbose output should be used.
  const bool _verbose;
}; // class TransportSolver

extern template class TransportSolver2D<TWDiamondDifference2D>;
extern template class TransportSolver2D<DiamondDifference2D>;
