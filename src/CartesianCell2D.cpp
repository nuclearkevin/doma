#include "CartesianCell2D.h"

#include "BrickMesh2D.h"

CartesianCell2D::CartesianCell2D(const double & lx, const double & ly,
                                 const double & x_center, const double & y_center,
                                 unsigned int cell_id, unsigned int block_id, BrickMesh2D * parent_mesh)
  : _cell_id(cell_id),
    _block_id(block_id),
    _x_c(std::move(x_center)),
    _y_c(std::move(y_center)),
    _l_x(std::move(lx)),
    _l_y(std::move(ly)),
    _area(_l_x * _l_y),
    _current_iteration_source(0.0),
    _current_scalar_flux(0.0),
    _parent_mesh(parent_mesh),
    _neighbors({nullptr, nullptr, nullptr, nullptr})
{ }

void
CartesianCell2D::addNeighbor(const CartesianCell2D * cell, CertesianFaceSide side)
{
  _neighbors[static_cast<unsigned int>(side)] = cell;
}

double
CartesianCell2D::boundaryFlux(CertesianFaceSide side, unsigned int ordinate_index)
{
  if (_parent_mesh->_bcs[static_cast<unsigned int>(side)] == BoundaryCondition::Vacuum)
    return 0.0;

  return 0.0;
}

// Helper to fetch the material properties of this cell.
const MaterialProps &
CartesianCell2D::getMatProps() const
{
  return _parent_mesh->_block_mat_info.at(_block_id);
}

bool
CartesianCell2D::hasStepSource() const
{
  return _parent_mesh->_block_step_src.count(_block_id) != 0u;
}

const SourceStep &
CartesianCell2D::getSourceStep() const
{
  return _parent_mesh->_block_step_src.at(_block_id);
}

// Check to see if a point is in the cell.
bool
CartesianCell2D::pointInCell(const double & x, const double & y) const
{
  const bool in_x = ((_x_c - (0.5 * _l_x)) <= x && (_x_c + (0.5 * _l_x)) >= x);
  const bool in_y = ((_y_c - (0.5 * _l_y)) <= y && (_y_c + (0.5 * _l_y)) >= y);

  return in_x && in_y;
}

const double &
CartesianCell2D::interfaceFlux(CertesianFaceSide side, unsigned int tid) const
{
  return _parent_mesh->_interface_angular_fluxes[tid][_cell_id][static_cast<unsigned int>(side)];
}

void
CartesianCell2D::setInterfaceFlux(CertesianFaceSide side, const double & val, unsigned int tid)
{
  _parent_mesh->_interface_angular_fluxes[tid][_cell_id][static_cast<unsigned int>(side)] = val;
}

void
CartesianCell2D::setAllInterfaceFluxes(const double & val, unsigned int tid)
{
  _parent_mesh->_interface_angular_fluxes[tid][_cell_id].fill(val);
}

void
CartesianCell2D::accumulateSweptFlux(const double & val, unsigned int tid)
{
  _parent_mesh->_swept_scalar_flux[tid][_cell_id] += val;
}

void
CartesianCell2D::setSweptFlux(const double & val, unsigned int tid)
{
  _parent_mesh->_swept_scalar_flux[tid][_cell_id] = val;
}

double
CartesianCell2D::getSweptFlux(unsigned int tid) const
{
  return _parent_mesh->_swept_scalar_flux[tid][_cell_id];
}
