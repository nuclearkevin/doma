#include "CartesianCell1D.h"

#include "BrickMesh1D.h"

CartesianCell1D::CartesianCell1D(const double & lx, const double & x_center, unsigned int cell_id,
                                 unsigned int block_id, BrickMesh1D * parent_mesh)
  : _cell_id(cell_id),
    _block_id(block_id),
    _x_c(std::move(x_center)),
    _l_x(std::move(lx)),
    _mg_source(0.0),
    _current_iteration_source(0.0),
    _current_scalar_flux(0.0),
    _current_current(0.0),
    _current_interface_sf({0.0, 0.0}),
    _parent_mesh(parent_mesh),
    _neighbors({nullptr, nullptr})
{ }

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
  return _parent_mesh->_block_mat_info.at(_block_id);
}

bool
CartesianCell1D::hasStepSource() const
{
  return _parent_mesh->_block_step_src.count(_block_id) != 0u;
}

const SourceStep &
CartesianCell1D::getSourceStep() const
{
  return _parent_mesh->_block_step_src.at(_block_id);
}

// Check to see if a point is in the cell.
bool
CartesianCell1D::pointInCell(const double & x) const
{
  return ((_x_c - (0.5 * _l_x)) <= x && (_x_c + (0.5 * _l_x)) >= x);
}

const double &
CartesianCell1D::interfaceFlux(CertesianFaceSide side, unsigned int tid) const
{
  return _parent_mesh->_interface_angular_fluxes[tid][_cell_id][static_cast<unsigned int>(side)];
}

void
CartesianCell1D::setInterfaceFlux(CertesianFaceSide side, const double & val, unsigned int tid)
{
  _parent_mesh->_interface_angular_fluxes[tid][_cell_id][static_cast<unsigned int>(side)] = val;
}

void
CartesianCell1D::setAllInterfaceFluxes(const double & val, unsigned int tid)
{
  _parent_mesh->_interface_angular_fluxes[tid][_cell_id].fill(val);
}

void
CartesianCell1D::accumulateSweptFlux(const double & val, unsigned int tid)
{
  _parent_mesh->_swept_scalar_flux[tid][_cell_id] += val;
}

void
CartesianCell1D::setSweptFlux(const double & val, unsigned int tid)
{
  _parent_mesh->_swept_scalar_flux[tid][_cell_id] = val;
}

double
CartesianCell1D::getSweptFlux(unsigned int tid) const
{
  return _parent_mesh->_swept_scalar_flux[tid][_cell_id];
}

void
CartesianCell1D::accumulateSweptCurrent(const double & val, unsigned int tid)
{
  _parent_mesh->_swept_net_current[tid][_cell_id] += val;
}

void
CartesianCell1D::setSweptCurrent(const double & val, unsigned int tid)
{
  _parent_mesh->_swept_net_current[tid][_cell_id] = val;
}

double
CartesianCell1D::getSweptCurrent(unsigned int tid)
{
  return _parent_mesh->_swept_net_current[tid][_cell_id];
}

void
CartesianCell1D::accumulateSweptInterfaceSF(const double & val, CertesianFaceSide side, unsigned int tid)
{
  const unsigned int side_offset = side == CertesianFaceSide::Right ? 1 : 0;
  _parent_mesh->_interface_scalar_fluxes[tid][_cell_id + side_offset] += val;
}

void
CartesianCell1D::setSweptInterfaceSF(const double & val, CertesianFaceSide side, unsigned int tid)
{
  const unsigned int side_offset = side == CertesianFaceSide::Right ? 1 : 0;
  _parent_mesh->_interface_scalar_fluxes[tid][_cell_id + side_offset] = val;
}

double
CartesianCell1D::getSweptInterfaceSF(CertesianFaceSide side, unsigned int tid)
{
  const unsigned int side_offset = side == CertesianFaceSide::Right ? 1 : 0;
  return _parent_mesh->_interface_scalar_fluxes[tid][_cell_id + side_offset];
}
