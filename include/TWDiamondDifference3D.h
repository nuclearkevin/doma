#pragma once

#include "TransportBase.h"
#include "CartesianCell3D.h"

#include <cmath>

// A class which implements the theta-weighted diamond difference cell balance equations.
// The values of theta are the ones used by TORT and DENOVO.
class TWDiamondDifference3D
{
public:
  TWDiamondDifference3D() = default;

  void solve(CartesianCell3D & cell, const double & angular_weight,
             const double & abs_mu, const double & abs_eta,
             const double & abs_xi, unsigned int ordinate_index,
             unsigned int group, CertesianFaceSide x_uw,
             CertesianFaceSide x_dw, CertesianFaceSide y_uw,
             CertesianFaceSide y_dw, CertesianFaceSide z_uw,
             CertesianFaceSide z_dw) const
  {
    const auto & p = cell.getMatProps();

    constexpr double theta = 0.9;
    constexpr double theta_s = 0.9;

    // Grab the upwind interfacing angular fluxes. The branches handle the boundary conditions.
    double x_uw_af = cell.neighbor(x_uw) ? cell.neighbor(x_uw)->interfaceFlux(x_dw) : cell.boundaryFlux(x_uw, ordinate_index);
    double y_uw_af = cell.neighbor(y_uw) ? cell.neighbor(y_uw)->interfaceFlux(y_dw) : cell.boundaryFlux(y_uw, ordinate_index);
    double z_uw_af = cell.neighbor(z_uw) ? cell.neighbor(z_uw)->interfaceFlux(z_dw) : cell.boundaryFlux(z_uw, ordinate_index);

    // Compute the weighting factors.
    const double A = cell._l_y * cell._l_z;
    const double B = cell._l_x * cell._l_z;
    const double C = cell._l_x * cell._l_y;
    const double V = cell._volume;

    double a = 0.5;
    if (x_uw_af > 1e-8)
    {
      a = cell._current_iteration_source[group] * V * theta_s + (abs_eta * B * y_uw_af + abs_xi * C * z_uw_af) * theta + abs_mu * A * x_uw_af;
      a /= (p._g_total[group] * V + 2.0 * abs_eta * B + 2.0 * abs_xi * C) * x_uw_af;
      a = -1.0 * (a - 1.0);
    }

    double b = 0.5;
    if (y_uw_af > 1e-8)
    {
      b = cell._current_iteration_source[group] * V * theta_s + (abs_mu * A * x_uw_af + abs_xi * C * z_uw_af) * theta + abs_eta * B * y_uw_af;
      b /= (p._g_total[group] * V + 2.0 * abs_mu * A + 2.0 * abs_xi * C) * y_uw_af;
      b = -1.0 * (b - 1.0);
    }

    double c = 0.5;
    if (z_uw_af > 1e-8)
    {
      c = cell._current_iteration_source[group] * V * theta_s + (abs_mu * A * x_uw_af + abs_eta * B * y_uw_af) * theta + abs_xi * C * z_uw_af;
      c /= (p._g_total[group] * V + 2.0 * abs_mu * A + 2.0 * abs_eta * B) * z_uw_af;
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
    double center_af = cell._current_iteration_source[group] * V + (abs_mu * A / a * x_uw_af) + (abs_eta * B / b * y_uw_af) + (abs_xi * C / c * z_uw_af);
    center_af /= p._g_total[group] * V + (abs_mu * A / a) + (abs_eta * B / b) + (abs_xi * C / c);

    double interface_af_x = (center_af / a) - (((1.0 - a) / a) * x_uw_af);
    double interface_af_y = (center_af / b) - (((1.0 - b) / b) * y_uw_af);
    double interface_af_z = (center_af / c) - (((1.0 - c) / c) * z_uw_af);

    // Add this angular flux's contribution to the scalar flux.
    cell._current_scalar_flux[group] += angular_weight * center_af;

    // Update the interface angular fluxes using the diamond difference closures. 2nd order accurate!
    cell.setInterfaceFlux(x_dw, interface_af_x);
    cell.setInterfaceFlux(y_dw, interface_af_y);
    cell.setInterfaceFlux(z_dw, interface_af_z);
  }
};
