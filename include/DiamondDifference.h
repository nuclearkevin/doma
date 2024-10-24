#pragma once

#include <cmath>

#include "BrickCellEquation.h"

// A class which implements the diamond difference cell balance equations.
class DiamondDifference : public BrickCellEquation
{
public:
  DiamondDifference() = default;

  virtual void solve(CartesianCell3D & cell, const double & angular_weight, bool first_iteration,
                     const double & abs_mu, const double & abs_eta, const double & abs_xi,
                     CertesianFaceSide x_uw, CertesianFaceSide x_dw,
                     CertesianFaceSide y_uw, CertesianFaceSide y_dw,
                     CertesianFaceSide z_uw, CertesianFaceSide z_dw) const override final
  {
    // Apply the fixed source on the first iteration, otherwise apply the scattering source.
    double cell_source = first_iteration ? cell._fixed_src / (4.0 * M_PI) : cell._sigma_s / (4.0 * M_PI) * cell._previous_scalar_flux;

    // Grab the upwind interfacing angular fluxes. The conditions handle the vacuum boundary conditions.
    double x_uw_af = cell.neighbor(x_uw) ? cell.neighbor(x_uw)->interfaceFlux(x_dw) : cell.interfaceFlux(x_uw);
    double y_uw_af = cell.neighbor(y_uw) ? cell.neighbor(y_uw)->interfaceFlux(y_dw) : cell.interfaceFlux(y_uw);
    double z_uw_af = cell.neighbor(z_uw) ? cell.neighbor(z_uw)->interfaceFlux(z_dw) : cell.interfaceFlux(z_uw);

    // Diamond difference approximation.
    // Computing cell-centered angular fluxes.
    double center_af = cell_source + (2.0 * x_uw_af * abs_mu / cell._l_x) + (2.0 * y_uw_af * abs_eta / cell._l_y) + (2.0 * z_uw_af * abs_xi / cell._l_z);
    center_af /= (2.0 * abs_mu / cell._l_x) + (2.0 * abs_eta / cell._l_x) + (2.0 * abs_xi / cell._l_x) + cell._sigma_t;

    double interface_af_x = 2.0 * center_af - x_uw_af;
    double interface_af_y = 2.0 * center_af - y_uw_af;
    double interface_af_z = 2.0 * center_af - z_uw_af;

    // Add this angular flux's contribution to the scalar flux.
    cell._current_scalar_flux += angular_weight * center_af;

    // Update the interface angular fluxes using the diamond difference closures. 2nd order accurate!
    cell.setInterfaceFlux(x_dw, interface_af_x);
    cell.setInterfaceFlux(y_dw, interface_af_y);
    cell.setInterfaceFlux(z_dw, interface_af_z);
  }
}; // class DiamondDifference
