#include "CartesianCell3D.h"

#include <iostream>

#include "BrickMesh3D.h"

CartesianCell3D::CartesianCell3D(const double & lx, const double & ly, const double & lz,
                                 const double & x_center, const double & y_center, const double & z_center,
                                 unsigned int cell_id, unsigned int block_id, BrickMesh3D * parent_mesh)
  : _cell_id(cell_id),
    _block_id(block_id),
    _x_c(std::move(x_center)),
    _y_c(std::move(y_center)),
    _z_c(std::move(z_center)),
    _l_x(std::move(lx)),
    _l_y(std::move(ly)),
    _l_z(std::move(lz)),
    _volume(_l_x * _l_y * _l_z),
    _parent_mesh(parent_mesh),
    _interface_angular_fluxes({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}),
    _neighbors({nullptr, nullptr, nullptr, nullptr, nullptr, nullptr})
{ }

// Initialize the flux storage data structures.
void
CartesianCell3D::initFluxes(unsigned int num_groups)
{
  _total_scalar_flux.resize(num_groups, 0.0);
  _current_iteration_source.resize(num_groups, 0.0);
  _current_scalar_flux.resize(num_groups, 0.0);
}

void
CartesianCell3D::addNeighbor(const CartesianCell3D * cell, CertesianFaceSide side)
{
  _neighbors[static_cast<unsigned int>(side)] = cell;
}

double
CartesianCell3D::boundaryFlux(CertesianFaceSide side, unsigned int ordinate_index)
{
  if (_parent_mesh->_bcs[static_cast<unsigned int>(side)] == BoundaryCondition::Vacuum)
    return 0.0;

  return 0.0;
}

const MaterialProps &
CartesianCell3D::getMatProps() const
{
  return _parent_mesh->_block_mat_info.at(_block_id);
}

// Check to see if a point is in the cell.
bool
CartesianCell3D::pointInCell(const double & x, const double & y, const double & z) const
{
  const bool in_x = ((_x_c - (0.5 * _l_x)) <= x && (_x_c + (0.5 * _l_x)) >= x);
  const bool in_y = ((_y_c - (0.5 * _l_y)) <= y && (_y_c + (0.5 * _l_y)) >= y);
  const bool in_z = ((_z_c - (0.5 * _l_z)) <= z && (_z_c + (0.5 * _l_z)) >= z);

  return in_x && in_y && in_z;
}
