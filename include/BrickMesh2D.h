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
  BrickMesh2D(const std::vector<unsigned int> & nx, const std::vector<unsigned int> & ny,
              const std::vector<double> & dx, const std::vector<double> & dy,
              const std::vector<unsigned int> & blocks, const std::array<BoundaryCondition, 4u> & bcs);

  // Set material properties for each cell.
  void addPropsToBlock(unsigned int block, const MaterialProps & props);

  // A function to initialize the group-wise scalar fluxes in each cell.
  void initFluxes(unsigned int num_groups);

  // Returns true if the point exists on the mesh, false if it does not. The flux at that point
  // will be stored in 'returned_flux' if the point is on the mesh.
  bool fluxAtPoint(const double & x, const double & y, unsigned int g, double & returned_flux) const;

  // Dump the flux to a text file.
  void dumpToTextFile(const std::string & file_name, unsigned int g);

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
  std::unordered_map<unsigned int, MaterialProps> _block_mat_info;
}; // class BrickMesh2D
