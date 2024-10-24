#pragma once

#include <cmath>
#include <iostream>

#include "BrickCellEquation.h"

// A class which implements the theta-weighted diamond difference cell balance equations.
// The values of theta are the ones used by TORT and DENOVO.
class TWDiamondDifference : public BrickCellEquation
{
public:
  TWDiamondDifference() = default;

  virtual void solve(CartesianCell3D & cell, const double & angular_weight, bool first_iteration,
                     const double & abs_mu, const double & abs_eta, const double & abs_xi,
                     CertesianFaceSide x_uw, CertesianFaceSide x_dw,
                     CertesianFaceSide y_uw, CertesianFaceSide y_dw,
                     CertesianFaceSide z_uw, CertesianFaceSide z_dw) const override final
  {
    constexpr double theta = 0.9;
    constexpr double theta_s = 0.9;

    // Apply the fixed source on the first iteration, otherwise apply the scattering source.
    double cell_source = first_iteration ? cell._fixed_src / (4.0 * M_PI) : cell._sigma_s / (4.0 * M_PI) * cell._previous_scalar_flux;

     // Grab the upwind interfacing angular fluxes. The branches handle the boundary conditions.
    double x_uw_af = cell.neighbor(x_uw) ? cell.neighbor(x_uw)->interfaceFlux(x_dw) : cell.interfaceFlux(x_uw);
    double y_uw_af = cell.neighbor(y_uw) ? cell.neighbor(y_uw)->interfaceFlux(y_dw) : cell.interfaceFlux(y_uw);
    double z_uw_af = cell.neighbor(z_uw) ? cell.neighbor(z_uw)->interfaceFlux(z_dw) : cell.interfaceFlux(z_uw);

    // Compute the weighting factors.
    const double A = cell._l_y * cell._l_z;
    const double B = cell._l_x * cell._l_z;
    const double C = cell._l_x * cell._l_y;
    const double V = cell._volume;

    double a = 0.5;
    if (x_uw_af > 1e-8)
    {
      a = cell_source * V * theta_s + (abs_eta * B * y_uw_af + abs_xi * C * z_uw_af) * theta + abs_mu * A * x_uw_af;
      a /= (cell._sigma_t * V + 2.0 * abs_eta * B + 2.0 * abs_xi * C) * x_uw_af;
      a = -1.0 * (a - 1.0);
    }

    double b = 0.5;
    if (y_uw_af > 1e-8)
    {
      b = cell_source * V * theta_s + (abs_mu * A * x_uw_af + abs_xi * C * z_uw_af) * theta + abs_eta * B * y_uw_af;
      b /= (cell._sigma_t * V + 2.0 * abs_mu * A + 2.0 * abs_xi * C) * y_uw_af;
      b = -1.0 * (b - 1.0);
    }

    double c = 0.5;
    if (z_uw_af > 1e-8)
    {
      c = cell_source * V * theta_s + (abs_mu * A * x_uw_af + abs_eta * B * y_uw_af) * theta + abs_xi * C * z_uw_af;
      c /= (cell._sigma_t * V + 2.0 * abs_mu * A + 2.0 * abs_eta * B) * z_uw_af;
      c = -1.0 * (c - 1.0);
    }

    // Clamp the weighting factors between 0.5 and 1.0.
    a = a < 0.5 ? 0.5 : a;
    a = a > 1.0 ? 1.0 : a;

    b = b < 0.5 ? 0.5 : b;
    b = b > 1.0 ? 1.0 : b;

    c = c < 0.5 ? 0.5 : c;
    c = c > 1.0 ? 1.0 : c;

    // Compute the center angular fluxes and the downwind interface angular fluxes using the weighted diamond difference approximation.
    double center_af = cell_source * V + (abs_mu * A / a * x_uw_af) + (abs_eta * B / b * y_uw_af) + (abs_xi * C / c * z_uw_af);
    center_af /= cell._sigma_t * V + (abs_mu * A / a) + (abs_eta * B / b) + (abs_xi * C / c);

    double interface_af_x = (center_af / a) - (((1.0 - a) / a) * x_uw_af);
    double interface_af_y = (center_af / b) - (((1.0 - b) / b) * y_uw_af);
    double interface_af_z = (center_af / c) - (((1.0 - c) / c) * z_uw_af);

    // Add this angular flux's contribution to the scalar flux.
    cell._current_scalar_flux += angular_weight * center_af;

    // Update the interface angular fluxes using the diamond difference closures. 2nd order accurate!
    cell.setInterfaceFlux(x_dw, interface_af_x);
    cell.setInterfaceFlux(y_dw, interface_af_y);
    cell.setInterfaceFlux(z_dw, interface_af_z);
  }
};
