#pragma once

#include "TransportBase.h"
#include "AngularPolynomials.h"

#include <vector>
#include <array>

//------------------------------------------------------------------------------
// Gauss-Legendre-Chebyshev product quadrature here.
//------------------------------------------------------------------------------
class GCAngularQuadrature
{
public:
  GCAngularQuadrature(unsigned int n_c, unsigned int n_l, unsigned int num_dims);

  unsigned int order(Octant oct) const;
  void direction(Octant oct, unsigned int n, double & mu, double & eta, double & xi) const;
  const double & weight(Octant oct, unsigned int n) const ;

  unsigned int totalOrder() const { return 2u * _n_c * _n_l; }

  unsigned int legendreOrder() const { return _n_l; }
  LegendrePolynomial getPolarLegendre() const { return _polar_quadrature; }

  unsigned int chebyshevOrder() const { return _n_c; }
  ChebyshevPolynomial getAzimuthalChebyshev() const { return _azimuthal_quadrature; }

private:
  Octant classifyDirection(const double & mu, const double & eta, const double & xi);

  // Number of Chebyshev quadrature points and number of Legendre quadrature
  // points, respectively.
  const unsigned int _n_c;
  const unsigned int _n_l;

  // Number of spatial dimensions of the problem.
  const unsigned int _n_d;

  LegendrePolynomial _polar_quadrature;
  ChebyshevPolynomial _azimuthal_quadrature;

  // Quadrature set directions indexed by octant of the unit sphere. Order is as follows:
  // PPP, PPM, PMP, PMM, MPP, MPM, MMP, MMM
  std::array<std::vector<double>, 8u> _quadrature_set_mu;
  std::array<std::vector<double>, 8u> _quadrature_set_eta;
  std::array<std::vector<double>, 8u> _quadrature_set_xi;
  std::array<std::vector<double>, 8u> _quadrature_set_weight;
}; // class GCAngularQuadrature
