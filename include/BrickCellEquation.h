#pragma once

#include "TransportBase.h"
#include "CartesianCell3D.h"

// An abstract class for other equation systems to derive from.
class BrickCellEquation
{
public:
  BrickCellEquation() {}

  virtual void solve(CartesianCell3D & cell, const double & angular_weight,
                     const double & abs_mu, const double & abs_eta,
                     const double & abs_xi, unsigned int ordinate_index,
                     CertesianFaceSide x_uw, CertesianFaceSide x_dw,
                     CertesianFaceSide y_uw, CertesianFaceSide y_dw,
                     CertesianFaceSide z_uw, CertesianFaceSide z_dw) const = 0;
}; // class BrickCellEquation
