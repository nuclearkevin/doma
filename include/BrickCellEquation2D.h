#pragma once

#include "TransportBase.h"
#include "CartesianCell2D.h"

// An abstract class for other equation systems to derive from.
class BrickCellEquation2D
{
public:
  BrickCellEquation2D() {}

  virtual void solve(CartesianCell2D & cell, const double & angular_weight,
                     const double & abs_mu, const double & abs_eta,
                     const double & abs_xi, unsigned int ordinate_index,
                     unsigned int group, CertesianFaceSide x_uw,
                     CertesianFaceSide x_dw, CertesianFaceSide y_uw,
                     CertesianFaceSide y_dw) const = 0;
}; // class BrickCellEquation3D
