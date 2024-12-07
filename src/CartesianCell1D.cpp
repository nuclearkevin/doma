#include "CartesianCell1D.h"

#include "BrickMesh1D.h"

CartesianCell1D::CartesianCell1D(const double & lx, const double & x_center, unsigned int cell_id,
                                 unsigned int block_id, BrickMesh1D * parent_mesh)
  : _cell_id(cell_id),
    _block_id(block_id),
    _x_c(std::move(x_center)),
    _l_x(std::move(lx)),
    _total_scalar_flux(0.0),
    _current_iteration_source(0.0),
    _current_scalar_flux(0.0),
    _parent_mesh(parent_mesh),
    _interface_angular_fluxes({0.0, 0.0}),
    _neighbors({nullptr, nullptr})
{ }

void
CartesianCell1D::initFluxes(unsigned int num_groups)
{
  _total_scalar_flux.resize(num_groups, 0.0);
  _current_iteration_source.resize(num_groups, 0.0);
  _current_scalar_flux.resize(num_groups, 0.0);
}

void
CartesianCell1D::addNeighbor(const CartesianCell1D * cell, CertesianFaceSide side)
{
  _neighbors[static_cast<unsigned int>(side)] = cell;
}

double
CartesianCell1D::boundaryFlux(CertesianFaceSide side, unsigned int ordinate_index)
{
  if (_parent_mesh->_bcs[static_cast<unsigned int>(side)] == BoundaryCondition::Vacuum)
    return 0.0;

  return 0.0;
}

// Helper to fetch the material properties of this cell.
const MaterialProps &
CartesianCell1D::getMatProps() const
{
  return _parent_mesh->_block_mat_info[_block_id];
}

// Check to see if a point is in the cell.
bool
CartesianCell1D::pointInCell(const double & x) const
{
  return ((_x_c - (0.5 * _l_x)) <= x && (_x_c + (0.5 * _l_x)) >= x);
}
