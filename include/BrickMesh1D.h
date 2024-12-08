#pragma once

#include <vector>
#include <iostream>

#include "CartesianCell1D.h"

#include "InputParameters.h"

template <typename T>
class TransportSolver1D;

// A regular 2D cartesian mesh.
class BrickMesh1D
{
public:
  BrickMesh1D(const std::vector<unsigned int> & nx, const std::vector<double> & dx,
              const std::vector<unsigned int> & blocks, const std::array<BoundaryCondition, 2u> & bcs,
              const std::unordered_map<unsigned int, MaterialProps> & props);

  // Make sure each block has material properties.
  void validateProps();

  // A function to initialize the group-wise scalar fluxes in each cell.
  void initFluxes(unsigned int num_groups);

  // Returns true if the point exists on the mesh, false if it does not. The flux at that point
  // will be stored in 'returned_flux' if the point is on the mesh.
  bool fluxAtPoint(const double & x, unsigned int g, double & returned_flux) const;

  // Dump the flux to a text file.
  void dumpToTextFile(const std::string & file_name);

  // A debug helper to print all of the mesh cells.
  void printAllBlocks();
  void printAllCoords();
  void printAllSideLengths();

  // Getters to fetch the properties of note.
  std::vector<CartesianCell1D *> & getBoundaryCells(CertesianFaceSide side)
  {
    return _boundary_cells[static_cast<unsigned int>(side)];
  }

private:
  friend class CartesianCell1D;

  template <typename T>
  friend class TransportSolver1D;

  // Number of brick cells per subdivision.
  const std::vector<unsigned int> _nx;

  unsigned int _tot_num_x;

  // Length of each mesh subdivision.
  const std::vector<double> _dx;

  // Number of energy groups.
  unsigned int _num_groups;

  // The unique material blocks.
  const std::vector<unsigned int> _blocks;

  // The cells themselves.
  unsigned int _num_cells;
  std::vector<CartesianCell1D> _cells;

  // The total length of the mesh.
  double _total_length;

  // A list of boundary conditions. Organized in the following order: Front, Back, Right, Left.
  const std::array<BoundaryCondition, 2u> _bcs;

  // A list of boundary cells. Organized in the following order: Front, Back.
  std::array<std::vector<CartesianCell1D *>, 2u> _boundary_cells;

  // An array of boundary angular fluxes for reflective bcs.
  std::array<std::vector<double>, 2u> _boundary_angular_fluxes;

  // Material properties for this mesh.
  const std::unordered_map<unsigned int, MaterialProps> & _block_mat_info;
}; // class BrickMesh1D
