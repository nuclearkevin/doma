#pragma once

#include <vector>
#include <iostream>

#include "CartesianCell2D.h"

#include "InputParameters.h"

template <typename T>
class TransportSolver2D;

// A regular 2D cartesian mesh.
class BrickMesh2D
{
public:
  BrickMesh2D(const InputParameters & params, const std::array<BoundaryCondition, 4u> & bcs, unsigned int num_threads);

  // Make sure each block has material properties.
  void validateProps();

  // Returns true if the point exists on the mesh, false if it does not. The flux at that point
  // will be stored in 'returned_flux' if the point is on the mesh.
  bool fluxAtPoint(const double & x, const double & y, unsigned int g, double & returned_flux) const;

  // Dump the flux to a text file.
  void dumpToTextFile(const std::string & file_name, bool only_flux = false);

  // Dump the DNPs to a text file.
  void dumpDNPsToTextFile(const std::string & file_name);

  // A debug helper to print all of the mesh cells.
  void printAllBlocks();
  void printAllCoords();
  void printAllSideLengths();

  // Getters to fetch the properties of note.
  std::vector<CartesianCell2D *> & getBoundaryCells(CertesianFaceSide side)
  {
    return _boundary_cells[static_cast<unsigned int>(side)];
  }

private:
  friend class CartesianCell2D;

  template <typename T>
  friend class TransportSolver2D;

  // Number of brick cells per subdivision.
  const std::vector<unsigned int> _nx;
  const std::vector<unsigned int> _ny;

  unsigned int _tot_num_x;
  unsigned int _tot_num_y;

  // Length of each mesh subdivision.
  const std::vector<double> _dx;
  const std::vector<double> _dy;

  // Number of energy groups.
  unsigned int _num_groups;

  // The unique material blocks.
  const std::vector<unsigned int> _blocks;

  // The cells themselves.
  unsigned int _num_cells;
  std::vector<CartesianCell2D> _cells;

  // The total area of the mesh.
  double _total_area;

  // A list of boundary conditions. Organized in the following order: Front, Back, Right, Left.
  const std::array<BoundaryCondition, 4u> _bcs;

  // A list of boundary cells. Organized in the following order: Front, Back, Right, Left.
  std::array<std::vector<CartesianCell2D *>, 4u> _boundary_cells;

  // An array of boundary angular fluxes for reflective bcs.
  std::array<std::vector<double>, 4u> _boundary_angular_fluxes;

  // Material properties for this mesh.
  const std::unordered_map<unsigned int, MaterialProps> & _block_mat_info;

  // Source step transients for this mesh.
  const std::unordered_map<unsigned int, SourceStep> & _block_step_src;

  // The number of OpenMP threads to use.
  const unsigned int _num_threads;

  // Data structures that store the scalar flux as it's being swept.
  std::vector<std::vector<double>> _swept_scalar_flux;
  // The interface angular fluxes. Downwind fluxes are computed by the cell, upwind fluxes are pulled from neighboring cells.
  // Organized in the following order: Front, Back, Right, Left.
  std::vector<std::vector<std::array<double, 4>>> _interface_angular_fluxes;
}; // class BrickMesh2D
