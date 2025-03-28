#pragma once

#include <array>

#include "Parallel.h"

#include "TransportBase.h"
#include "InputParameters.h"

class BrickMesh3D;

template <typename T>
class TransportSolver3D;

// A 3D cartesian cell.
class CartesianCell3D
{
public:
  CartesianCell3D(const double & lx, const double & ly, const double & lz,
                  const double & x_center, const double & y_center, const double & z_center,
                  unsigned int cell_id, unsigned int block_id, BrickMesh3D * parent_mesh);

  // Link the cell's neighbors to the cell.
  void addNeighbor(const CartesianCell3D * cell, CertesianFaceSide side);

  // Helper to fetch the interface flux.
  const double & interfaceFlux(CertesianFaceSide side, unsigned int tid) const;
  void setInterfaceFlux(CertesianFaceSide side, const double & val, unsigned int tid);
  void setAllInterfaceFluxes(const double & val, unsigned int tid);

  // Helper functions to get or modify the swept scalar flux.
  void accumulateSweptFlux(const double & val, unsigned int tid);
  void setSweptFlux(const double & val, unsigned int tid);
  double getSweptFlux(unsigned int tid) const;

  // Helper to fetch BCs.
  double boundaryFlux(CertesianFaceSide side, unsigned int ordinate_index);

  // Helper to fetch cell neighbors.
  const CartesianCell3D * neighbor(CertesianFaceSide side) const { return _neighbors[static_cast<unsigned int>(side)]; }

  // Helper to fetch the material properties of this cell.
  const MaterialProps & getMatProps() const;

  // Check to see if this cell has a step source.
  bool hasStepSource() const;

  // Helper to fetch the source step transients of this cell.
  const SourceStep & getSourceStep() const;

  // Check to see if a point is in the cell.
  bool pointInCell(const double & x, const double & y, const double & z) const;

  double getSideLength(CertesianAxis axis) const
  {
    switch (axis)
    {
      case CertesianAxis::X: return _l_x;
      case CertesianAxis::Y: return _l_y;
      case CertesianAxis::Z: return _l_z;
    }
    return 0.0;
  }

  double getFaceArea(CertesianFaceSide side) const
  {
    switch (side)
    {
      case CertesianFaceSide::Front:  return _l_x * _l_z;
      case CertesianFaceSide::Back:   return _l_x * _l_z;
      case CertesianFaceSide::Right:  return _l_y * _l_z;
      case CertesianFaceSide::Left:   return _l_y * _l_z;
      case CertesianFaceSide::Top:    return _l_x * _l_y;
      case CertesianFaceSide::Bottom: return _l_x * _l_y;
    }
    return 0.0;
  }

  double getVolume() const
  {
    return _volume;
  }

  const unsigned int _cell_id;
  const unsigned int _block_id;

  // Geometric properties.
  const double _x_c;    // X component of the centroid.
  const double _y_c;    // Y component of the centroid.
  const double _z_c;    // Z component of the centroid.

  const double _l_x;    // X length of the rectangular prism.
  const double _l_y;    // Y length of the rectangular prism.
  const double _l_z;    // Z length of the rectangular prism.

  const double _volume; // Volume of the rectangular prism.

  // Flux properties.
  std::vector<double> _total_scalar_flux;        // The accumulated scalar flux at the current multi-group iteration index.
  std::vector<double> _prev_mg_scalar_flux;      // The accumulated scalar flux at the previous multi-group iteration index.
  double              _current_iteration_source; // The current iteration scattering source.
  double              _current_scalar_flux;      // For accumulating the current iteration's scalar flux while the angular flux is being swept.

  // Previous timestep fluxes.
  std::vector<double> _last_t_scalar_flux;

  // Delayed neutron precursors.
  std::vector<double> _current_t_dnps;
  std::vector<double> _last_t_dnps;

protected:
  friend class BrickMesh3D;

  template <typename T>
  friend class TransportSolver3D;

  // The mesh which owns this cell.
  BrickMesh3D * _parent_mesh;

  // A list of neighboring cells. Organized in the following order: Front, Back, Right, Left, Top, Bottom.
  std::array<const CartesianCell3D *, 6u> _neighbors;
}; // class CartesianCell3D
