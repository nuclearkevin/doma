#pragma once

#include "TransportBase.h"
#include "AngularPolynomials.h"

#include <vector>
#include <array>

// Gauss-Legendre angular quadrature.
class GLAngularQuadrature
{
public:
  GLAngularQuadrature(unsigned int n_l);

  unsigned int order() const { return _n_l; }
  double direction(unsigned int n) const;
  const double & weight(unsigned int n) const ;

  LegendrePolynomial getPolarLegendre() const { return _polar_quadrature; }

private:
  // Number of Legendre quadrature points.
  const unsigned int _n_l;

  LegendrePolynomial _polar_quadrature;
}; // class GLAngularQuadrature
