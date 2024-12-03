#include "CartesianCell2D.h"

#include "BrickMesh2D.h"

CartesianCell2D::CartesianCell2D(const double & lx, const double & ly,
                                 const double & x_center, const double & y_center,
                                 unsigned int cell_id, unsigned int block_id, BrickMesh2D * parent_mesh)
  : _cell_id(cell_id),
    _block_id(block_id),
    _sigma_t(0.0),
    _sigma_s(0.0),
    _fixed_src(0.0),
    _x_c(std::move(x_center)),
    _y_c(std::move(y_center)),
    _l_x(std::move(lx)),
    _l_y(std::move(ly)),
    _area(_l_x * _l_y),
    _total_scalar_flux(0.0),
    _current_iteration_source(0.0),
    _current_scalar_flux(0.0),
    _parent_mesh(parent_mesh),
    _interface_angular_fluxes({0.0, 0.0, 0.0, 0.0}),
    _neighbors({nullptr, nullptr, nullptr, nullptr})
{ }

void
CartesianCell2D::applyProperties(const double & sigma_total, const double & sigma_scattering, const double & fixed_source)
{
  _sigma_t = std::move(sigma_total);
  _sigma_s = std::move(sigma_scattering);
  _fixed_src = std::move(fixed_source);
}

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

// Check to see if a point is in the cell.
bool
CartesianCell2D::pointInCell(const double & x, const double & y) const
{
  const bool in_x = ((_x_c - (0.5 * _l_x)) <= x && (_x_c + (0.5 * _l_x)) >= x);
  const bool in_y = ((_y_c - (0.5 * _l_y)) <= y && (_y_c + (0.5 * _l_y)) >= y);

  return in_x && in_y;
}
