#pragma once

#include "TransportBase.h"
#include "CartesianCell1D.h"

#include <cmath>

// A class which implements the diamond difference cell balance equations.
class DiamondDifference1D
{
public:
  DiamondDifference1D() = default;

  void solve(CartesianCell1D & cell, const double & angular_weight,
             const double & abs_mu, unsigned int ordinate_index,
             unsigned int group, CertesianFaceSide x_uw,
             CertesianFaceSide x_dw, RunMode mode = RunMode::FixedSrc,
             const double & dt = 0.0) const
  {
    const auto & p = cell.getMatProps();
    auto tot = p._g_total[group];
    if (mode == RunMode::Transient)
      tot += p._g_inv_v[group] / dt;

    // Grab the upwind interfacing angular fluxes. The conditions handle the vacuum boundary conditions.
    double x_uw_af = cell.neighbor(x_uw) ? cell.neighbor(x_uw)->interfaceFlux(x_dw) : cell.boundaryFlux(x_uw, ordinate_index);

    // Diamond difference approximation.
    // Computing cell-centered angular fluxes.
    double center_af = cell._current_iteration_source + (2.0 * x_uw_af * abs_mu / cell._l_x);
    center_af /= (2.0 * abs_mu / cell._l_x) + tot;

    double interface_af_x = 2.0 * center_af - x_uw_af;

    // Add this angular flux's contribution to the scalar flux.
    cell._current_scalar_flux += angular_weight * center_af;

    // Update the interface angular fluxes using the diamond difference closures. 2nd order accurate!
    cell.setInterfaceFlux(x_dw, interface_af_x);
  }
}; // class DiamondDifference1D
