#pragma once

#include <vector>
#include <iostream>

#include "CartesianCell3D.h"

template <typename T>
class TransportSolver3D;

// A regular 3D cartesian mesh.
class BrickMesh3D
{
public:
  BrickMesh3D(const std::vector<unsigned int> & nx, const std::vector<unsigned int> & ny, const std::vector<unsigned int> & nz,
              const std::vector<double> & dx, const std::vector<double> & dy, const std::vector<double> & dz,
              const std::vector<unsigned int> & blocks, const std::array<BoundaryCondition, 6u> & bcs);

  // Set material properties for each cell.
  void addPropsToBlock(unsigned int block, const double & sigma_total, const double & sigma_scattering, const double & fixed_source);

  // Returns true if the point exists on the mesh, false if it does not. The flux at that point
  // will be stored in 'returned_flux' if the point is on the mesh.
  bool fluxAtPoint(const double & x, const double & y, const double & z, double & returned_flux) const;

  // Dump the flux to a text file.
  void dumpToTextFile(const std::string & file_name);

  // A debug helper to print all of the mesh cells.
  void printAllBlocks();
  void printAllCoords();
  void printAllSideLengths();

    // Getters to fetch the properties of note.
  std::vector<CartesianCell3D *> & getBoundaryCells(CertesianFaceSide side)
  {
    return _boundary_cells[static_cast<unsigned int>(side)];
  }

private:
  friend class CartesianCell3D;

  template <typename T>
  friend class TransportSolver3D;

  // Number of brick cells per subdivision.
  const std::vector<unsigned int> _nx;
  const std::vector<unsigned int> _ny;
  const std::vector<unsigned int> _nz;

  unsigned int _tot_num_x;
  unsigned int _tot_num_y;
  unsigned int _tot_num_z;

  // Length of each mesh subdivision.
  const std::vector<double> _dx;
  const std::vector<double> _dy;
  const std::vector<double> _dz;

  // The unique material blocks.
  const std::vector<unsigned int> _blocks;

  // The cells themselves.
  unsigned int _num_cells;
  std::vector<CartesianCell3D> _cells;

  // The total volume of the mesh.
  double _total_volume;

  // A list of boundary conditions. Organized in the following order: Front, Back, Right, Left, Top, Bottom.
  const std::array<BoundaryCondition, 6u> _bcs;

  // A list of boundary cells. Organized in the following order: Front, Back, Right, Left, Top, Bottom.
  std::array<std::vector<CartesianCell3D *>, 6u> _boundary_cells;

  // An array of boundary angular fluxes for reflective bcs.
  std::array<std::vector<double>, 6u> _boundary_angular_fluxes;
}; // class BrickMesh3D
