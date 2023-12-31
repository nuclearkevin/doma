#include "CartesianCell3D.h"

#include <iostream>

CartesianCell3D::CartesianCell3D(const double & lx, const double & ly, const double & lz,
                                 const double & x_center, const double & y_center, const double & z_center,
                                 unsigned int cell_id, unsigned int block_id)
  : _cell_id(cell_id),
    _block_id(block_id),
    _sigma_t(0.0),
    _sigma_s(0.0),
    _fixed_src(0.0),
    _x_c(std::move(x_center)),
    _y_c(std::move(y_center)),
    _z_c(std::move(z_center)),
    _l_x(std::move(lx)),
    _l_y(std::move(ly)),
    _l_z(std::move(lz)),
    _volume(_l_x * _l_y * _l_z),
    _total_scalar_flux(0.0),
    _previous_scalar_flux(0.0),
    _current_scalar_flux(0.0),
    _interface_angular_fluxes({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}),
    _neighbors({nullptr, nullptr, nullptr, nullptr, nullptr, nullptr})
{ }

void
CartesianCell3D::applyProperties(const double & sigma_total, const double & sigma_scattering, const double & fixed_source)
{
  _sigma_t = std::move(sigma_total);
  _sigma_s = std::move(sigma_scattering);
  _fixed_src = std::move(fixed_source);
}

void
CartesianCell3D::addNeighbor(const CartesianCell3D * cell, CertesianFaceSide side)
{
  _neighbors[static_cast<unsigned int>(side)] = cell;
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
