#pragma once

#include <vector>

//------------------------------------------------------------------------------
// Legendre polynomials here.
//------------------------------------------------------------------------------
// Implementation inspired by:
// http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
class LegendrePolynomial
{
public:
  LegendrePolynomial(unsigned int degree);

  unsigned int degree() const { return _degree; }
  const double & root(unsigned int index) const { return _roots[index]; }
  const double & weight(unsigned int index) const { return _weights[index]; }

  const std::vector<double> & getRoots() const { return _roots; }
  const std::vector<double> & getWeights() const { return _weights; }

private:
  // Evaluate the value of a Legendre polynomial at x using the recurrance
  // relationship.
  double evaluate(double x);

  // Evaluate the value of a Legendre polynomial's derivative at x using the
  // recurrance relationship.
  double evaluateDerivative(double x);

  const unsigned int _degree;

  // Gauss-Legendre weights are easier to compute while finding zeros.
  std::vector<double> _roots;
  std::vector<double> _weights;
}; // class LegendrePolynomial

//------------------------------------------------------------------------------
// Chebyshev polynomials here.
//------------------------------------------------------------------------------
class ChebyshevPolynomial
{
public:
  ChebyshevPolynomial(unsigned int degree);

  unsigned int degree() const { return _degree; }
  const double & root(unsigned int index) const { return _roots[index]; }
  const double & angularRoot(unsigned int index) const { return _angular_roots[index]; }
  const double & weight(unsigned int index) const { return _weights[index]; }

  const std::vector<double> & getRoots() const { return _roots; }
  const std::vector<double> & getAngularRoots() const { return _angular_roots; }
  const std::vector<double> & getWeights() const { return _weights; }

private:
  const unsigned int _degree;

  // Storing weights to mirror LegendrePolynomial.
  // _roots stores the roots on -1 \leq y \leq 1
  // _angular_roots stores the roots on 0 \leq \omega \leq 2\pi
  std::vector<double> _roots;
  std::vector<double> _angular_roots;
  std::vector<double> _weights;
}; // class ChebyshevPolynomial

//------------------------------------------------------------------------------
// Gauss-Legendre-Chebyshev product quadrature here.
//------------------------------------------------------------------------------
class GCAngularQuadrature
{
public:
  GCAngularQuadrature(unsigned int n_c, unsigned int n_l);

  unsigned int totalOrder() const;
  void direction(unsigned int n, double & mu, double & eta, double & xi) const;
  const double & weight(unsigned int n) const ;
  const std::vector<double> & getWeights() const;

  const double & getPolarRoot(unsigned int n) const;
  const double & getAzimuthalAngularRoot(unsigned int n) const;

  unsigned int legendreOrder() const { return _n_l; }
  LegendrePolynomial getPolarLegendre() const { return _polar_quadrature; }

  unsigned int chebyshevOrder() const { return _n_c; }
  ChebyshevPolynomial getAzimuthalChebyshev() const { return _azimuthal_quadrature; }

private:
  // Generate a weight-ordinate pair for 3 dimensional problems.
  void generateWeightOrdiantePair(unsigned int i, unsigned int j);

  // Number of Chebyshev quadrature points and number of Legendre quadrature
  // points, respectively.
  const unsigned int _n_c;
  const unsigned int _n_l;

  LegendrePolynomial _polar_quadrature;
  ChebyshevPolynomial _azimuthal_quadrature;

  std::vector<double> _quadrature_set_mu;
  std::vector<double> _quadrature_set_eta;
  std::vector<double> _quadrature_set_xi;
  std::vector<double> _quadrature_set_weight;
}; // class GCAngularQuadrature
