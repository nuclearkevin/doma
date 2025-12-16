#pragma once

#include <unordered_map>
#include <vector>
#include <string>

#include "TransportBase.h"

// A container to collect properties of step insertion / removal transients for external sources.
struct SourceStep
{
  // The type of step transient.
  StepType _type;
  // The time at which the source is inserted.
  double _insert_time;
  // The time at which the source is removed.
  double _remove_time;

  // Vector of energy-group-wise isotropic source intensities.
  std::vector<double> _g_src;
};

// A container to collect material properties.
struct MaterialProps
{
  // Vector of energy-group-wise total cross sections.
  std::vector<double> _g_total;
  // Vector of energy-group-wise fission production cross sections.
  std::vector<double> _g_prod;
  // Vector of energy-group-wise prompt fission spectra.
  std::vector<double> _g_chi_p;
  // Vector of energy-group-wise inverse velocities.
  std::vector<double> _g_inv_v;
  // Scattering matrix of energy-group to energy-group isotropic scattering cross sections flattened into a vector.
  std::vector<double> _g_g_scatter_mat;
  // Vector of energy-group-wise isotropic source intensities.
  std::vector<double> _g_src;

  // Number of delayed neutron precusor groups.
  unsigned int _num_d_groups;

  // The delayed-group to energy-group delayed fission spectra matrix flattened into a vector.
  std::vector<double> _n_g_chi_d;
  // The energy-group to delayed-group values of the delayed neutron fraction matrix flattened into a vector.
  std::vector<double> _g_n_beta;
  // Vector of delayed neutron precursor decay constants.
  std::vector<double> _n_lambda;
};

struct InputParameters
{
  // The execution mode of the simulation.
  RunMode _mode;

  // Number of energy groups.
  unsigned int _num_e_groups;

  // Initial time for transient simulations.
  double _t0;
  // Final time for transient simulations.
  double _t1;
  // Number of timesteps to take.
  unsigned int _num_steps;
  // Initial conditions to use.
  TransientIC _ic;

  // Source iteration tolerance.
  double _src_it_tol;
  // The multi-group iteration (Gauss-Seidel) tolerance.
  double _gs_tol;
  // The k_{eff} convergence criteria for eigenvalue calculations.
  double _k_tol;

  // Whether diffusion synthetic acceleration should be enabled or not.
  bool _dsa;

  // Maximum number of source iterations.
  unsigned int _num_src_it;
  // Maximum number of multi-group iterations.
  unsigned int _num_mg_it;

  // Angular quadrature parameters.
  unsigned int _num_polar;
  unsigned int _num_azimuthal;

  // The equation system to use.
  EquationType _eq_type;

  // Mesh settings.
  // Number of dimensions for the simulation.
  unsigned int _num_dims;

  // The uniform mesh refinement.
  unsigned int _refinement;

  // The number of intervals in x, y, and z.
  std::vector<unsigned int> _x_intervals;
  std::vector<unsigned int> _y_intervals;
  std::vector<unsigned int> _z_intervals;

  // The length of each x, y, and z interval.
  std::vector<double> _dx;
  std::vector<double> _dy;
  std::vector<double> _dz;

  // The blocks for each unique interval.
  std::vector<unsigned int> _blocks;

  // Boundary conditions, currently unused as the only supported BC is vacuum.
  std::unordered_map<CertesianFaceSide, BoundaryCondition> _bcs;

  // Materials.
  std::unordered_map<unsigned int, MaterialProps> _block_mat_info;

  // Source step insertion / removal transients.
  std::unordered_map<unsigned int, SourceStep> _block_step_src;
};

InputParameters parseInputParameters(const std::string & file_path);
