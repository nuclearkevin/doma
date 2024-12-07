#pragma once

#include <array>

#include "TransportBase.h"
#include "InputParameters.h"

class BrickMesh1D;

template <typename T>
class TransportSolver1D;

// A 1D cartesian cell.
class CartesianCell1D
{
public:
  CartesianCell1D(const double & lx, const double & x_center, unsigned int cell_id,
                  unsigned int block_id, BrickMesh1D * parent_mesh);

  // Initialize the flux storage datastructures.
  void initFluxes(unsigned int num_groups);

  // Link the cell's neighbors to the cell.
  void addNeighbor(const CartesianCell1D * cell, CertesianFaceSide side);

  // Helper to fetch the interface flux.
  const double & interfaceFlux(CertesianFaceSide side) const { return _interface_angular_fluxes[static_cast<unsigned int>(side)]; }
  void setInterfaceFlux(CertesianFaceSide side, const double & val)
  {
    _interface_angular_fluxes[static_cast<unsigned int>(side)] = val;
  }

  // Helper to fetch BCs.
  double boundaryFlux(CertesianFaceSide side, unsigned int ordinate_index);

  // Helper to fetch cell neighbors.
  const CartesianCell1D * neighbor(CertesianFaceSide side) const { return _neighbors[static_cast<unsigned int>(side)]; }

  // Helper to fetch the material properties of this cell.
  const MaterialProps & getMatProps() const;

  // Check to see if a point is in the cell.
  bool pointInCell(const double & x) const;

  double getLength() const
  {
    return _l_x;
  }

  const unsigned int _cell_id;
  const unsigned int _block_id;

  // Geometric properties.
  const double _x_c;  // X component of the centroid.

  const double _l_x;  // X length of the rectangular prism.

  // Flux properties.
  std::vector<double> _total_scalar_flux;        // The sum of scalar fluxes from each scattering iteration.
  std::vector<double> _current_iteration_source; // The scattering source.
  std::vector<double> _current_scalar_flux;      // For accumulating the current iteration's scalar flux while the angular flux is being swept.

protected:
  friend class BrickMesh1D;

  template <typename T>
  friend class TransportSolver1D;

  // The mesh which owns this cell.
  BrickMesh1D * _parent_mesh;

  // The interface angular fluxes. Downwind fluxes are computed by the cell, upwind fluxes are pulled from neighboring cells.
  // Organized in the following order: Front, Back, Right, Left.
  std::array<double, 2u> _interface_angular_fluxes;

  // A list of neighboring cells. Organized in the following order: Front, Back.
  std::array<const CartesianCell1D *, 2u> _neighbors;
}; // class CartesianCell1D
