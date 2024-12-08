#pragma once

#include "TransportBase.h"
#include "CartesianCell3D.h"

#include <cmath>

// A class which implements the diamond difference cell balance equations.
class DiamondDifference3D
{
public:
  DiamondDifference3D() = default;

  void solve(CartesianCell3D & cell, const double & angular_weight,
             const double & abs_mu, const double & abs_eta,
             const double & abs_xi, unsigned int ordinate_index,
             unsigned int group, CertesianFaceSide x_uw,
             CertesianFaceSide x_dw, CertesianFaceSide y_uw,
             CertesianFaceSide y_dw, CertesianFaceSide z_uw,
             CertesianFaceSide z_dw) const
  {
    const auto & p = cell.getMatProps();

    // Grab the upwind interfacing angular fluxes. The conditions handle the vacuum boundary conditions.
    double x_uw_af = cell.neighbor(x_uw) ? cell.neighbor(x_uw)->interfaceFlux(x_dw) : cell.boundaryFlux(x_uw, ordinate_index);
    double y_uw_af = cell.neighbor(y_uw) ? cell.neighbor(y_uw)->interfaceFlux(y_dw) : cell.boundaryFlux(y_uw, ordinate_index);
    double z_uw_af = cell.neighbor(z_uw) ? cell.neighbor(z_uw)->interfaceFlux(z_dw) : cell.boundaryFlux(z_uw, ordinate_index);

    // Diamond difference approximation.
    // Computing cell-centered angular fluxes.
    double center_af = cell._current_iteration_source[group] + (2.0 * x_uw_af * abs_mu / cell._l_x) + (2.0 * y_uw_af * abs_eta / cell._l_y) + (2.0 * z_uw_af * abs_xi / cell._l_z);
    center_af /= (2.0 * abs_mu / cell._l_x) + (2.0 * abs_eta / cell._l_x) + (2.0 * abs_xi / cell._l_x) + p._g_total[group];

    double interface_af_x = 2.0 * center_af - x_uw_af;
    double interface_af_y = 2.0 * center_af - y_uw_af;
    double interface_af_z = 2.0 * center_af - z_uw_af;

    // Add this angular flux's contribution to the scalar flux.
    cell._current_scalar_flux[group] += angular_weight * center_af;

    // Update the interface angular fluxes using the diamond difference closures. 2nd order accurate!
    cell.setInterfaceFlux(x_dw, interface_af_x);
    cell.setInterfaceFlux(y_dw, interface_af_y);
    cell.setInterfaceFlux(z_dw, interface_af_z);
  }
}; // class DiamondDifference3D
