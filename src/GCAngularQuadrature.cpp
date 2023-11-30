#include "GCAngularQuadrature.h"

#include <iostream>
#include <cmath>

//------------------------------------------------------------------------------
// Legendre polynomials here.
//------------------------------------------------------------------------------
LegendrePolynomial::LegendrePolynomial(unsigned int degree) : _degree(std::move(degree))
{
  // Precompute and store the roots of the polynomial using Newton's method.
  _roots.resize(degree, 0.0);
  _weights.resize(degree, 0.0);

  double dr = 0.0;
  double x = 0.0;
  double v = 0.0;
  double d = 0.0;
  for (unsigned int i = 1; i <= _degree; ++i)
  {
    dr = 1;

    x = std::cos(M_PI * (static_cast<double>(i) - 0.25) /
                 (static_cast<double>(_degree) + 0.5)); // Magical initial guess.
    v = evaluate(x);
    d = evaluateDerivative(x);

    do
    {
      dr = v / d;
      x -= dr; // x_{i+1} = x_{i} - f(x_{i}) / f'(x_{i})
      v = evaluate(x);
      d = evaluateDerivative(x);
    } while (std::abs(dr) > 2e-16);

    _roots[i - 1u] = x;
    _weights[i - 1u] = 2.0 / ((1.0 - x * x) * d * d);
  }
}

// Evaluate the value of a Legendre polynomial at x using the recurrance
// relationship.
double
LegendrePolynomial::evaluate(double x)
{
  if (_degree == 0u)
    return 1.0;
  if (_degree == 1u)
    return x;

  double v = 1.0;
  double v_sub_1 = x;
  double v_sub_2 = 1.0;

  for (unsigned int i = 2; i <= _degree; ++i)
  {
    v = ((2.0 * static_cast<double>(i) - 1.0) * x * v_sub_1 - (i - 1.0) * v_sub_2) /
        static_cast<double>(i);

    v_sub_2 = v_sub_1;
    v_sub_1 = v;
  }

  return v;
}

// Evaluate the value of a Legendre polynomial's derivative at x using the
// recurrance relationship.
double
LegendrePolynomial::evaluateDerivative(double x)
{
  if (_degree == 0u)
    return 0.0;
  if (_degree == 1u)
    return 1.0;

  double v = 1.0;
  double d = 0.0;
  double v_sub_1 = x;
  double v_sub_2 = 1.0;
  double f = 1.0 / (x * x - 1.0);

  for (unsigned int i = 2; i <= _degree; ++i)
  {
    v = ((2.0 * static_cast<double>(i) - 1.0) * x * v_sub_1 -
         (static_cast<double>(i) - 1.0) * v_sub_2) /
        static_cast<double>(i);
    d = static_cast<double>(i) * f * (x * v - v_sub_1);

    v_sub_2 = v_sub_1;
    v_sub_1 = v;
  }

  return d;
}

//------------------------------------------------------------------------------
// Chebyshev polynomials here.
//------------------------------------------------------------------------------
ChebyshevPolynomial::ChebyshevPolynomial(unsigned int degree) : _degree(std::move(degree))
{
  _roots.resize(_degree, 0.0);
  _angular_roots.resize(_degree, 0.0);
  _weights.resize(_degree, M_PI / static_cast<double>(_degree));

  for (unsigned int i = 1; i <= _degree; ++i)
  {
    _angular_roots[i - 1u] =
        (2 * static_cast<double>(i) - 1.0) / (2.0 * static_cast<double>(_degree)) * M_PI;
    _roots[i - 1u] = std::cos(_angular_roots[i - 1u]);
  }
}

//------------------------------------------------------------------------------
// Gauss-Legendre-Chebyshev product quadrature here.
//------------------------------------------------------------------------------
GCAngularQuadrature::GCAngularQuadrature(unsigned int n_c,
                                         unsigned int n_l)
  : _n_c(std::move(n_c)),
    _n_l(std::move(n_l)),
    _polar_quadrature(std::move(n_l)),
    _azimuthal_quadrature(std::move(n_c))
{
  // Generate the quadrature directions and weights in an unordered fashion.
  std::vector<double> unordered_mus;
  std::vector<double> unordered_etas;
  std::vector<double> unordered_xis;
  std::vector<double> unordered_weights;

  for (unsigned int i = 1u; i <= _n_l; ++i)
  {
    for (unsigned int j = 1u; j <= _n_c; ++j)
    {
      const auto weight = _polar_quadrature.weight(i - 1u) * _azimuthal_quadrature.weight(j - 1u);
      const auto & mu = _polar_quadrature.root(i - 1u);
      const auto & omega = _azimuthal_quadrature.angularRoot(j - 1u);

      // Positive octant.
      unordered_mus.emplace_back(mu);
      unordered_etas.emplace_back(std::sqrt(1.0 - (mu * mu)) * std::cos(omega));
      unordered_xis.emplace_back(std::sqrt(1.0 - (mu * mu)) * std::sin(omega));
      unordered_weights.emplace_back(weight);
      // Negative octant.
      unordered_mus.emplace_back(-1.0 * mu);
      unordered_etas.emplace_back(-1.0 * std::sqrt(1.0 - (mu * mu)) * std::cos(omega));
      unordered_xis.emplace_back(-1.0 * std::sqrt(1.0 - (mu * mu)) * std::sin(omega));
      unordered_weights.emplace_back(weight);
    }
  }

  // Sort the quadrature weights and directions into bins depending on their octant of residence.
  for (unsigned int n = 0u; n < unordered_weights.size(); ++n)
  {
    auto oct = classifyDirection(unordered_mus[n], unordered_etas[n], unordered_xis[n]);
    _quadrature_set_mu[static_cast<unsigned int>(oct)].emplace_back(unordered_mus[n]);
    _quadrature_set_eta[static_cast<unsigned int>(oct)].emplace_back(unordered_etas[n]);
    _quadrature_set_xi[static_cast<unsigned int>(oct)].emplace_back(unordered_xis[n]);
    _quadrature_set_weight[static_cast<unsigned int>(oct)].emplace_back(unordered_weights[n]);
  }
}

unsigned int
GCAngularQuadrature::order(Octant oct) const
{
  return _quadrature_set_weight[static_cast<unsigned int>(oct)].size();
}

void
GCAngularQuadrature::direction(Octant oct, unsigned int n, double & mu, double & eta, double & xi) const
{
  mu = _quadrature_set_mu[static_cast<unsigned int>(oct)][n];
  eta = _quadrature_set_eta[static_cast<unsigned int>(oct)][n];
  xi = _quadrature_set_xi[static_cast<unsigned int>(oct)][n];
}

const double &
GCAngularQuadrature::weight(Octant oct, unsigned int n) const
{
  return _quadrature_set_weight[static_cast<unsigned int>(oct)][n];
}

Octant
GCAngularQuadrature::classifyDirection(const double & mu, const double & eta, const double & xi)
{
  if (mu > 0.0 && eta > 0.0 && xi > 0.0)
    return Octant::PPP;
  if (mu > 0.0 && eta > 0.0 && xi < 0.0)
    return Octant::PPM;
  if (mu > 0.0 && eta < 0.0 && xi > 0.0)
    return Octant::PMP;
  if (mu > 0.0 && eta < 0.0 && xi < 0.0)
    return Octant::PMM;
  if (mu < 0.0 && eta > 0.0 && xi > 0.0)
    return Octant::MPP;
  if (mu < 0.0 && eta > 0.0 && xi < 0.0)
    return Octant::MPM;
  if (mu < 0.0 && eta < 0.0 && xi > 0.0)
    return Octant::MMP;
  if (mu < 0.0 && eta < 0.0 && xi < 0.0)
    return Octant::MMM;

  std::cout << "Error: TransportSolver::classifyDirection(const double & mu, const double & eta, const double & xi)." << std::endl;
  exit(1);
  return Octant::PPP;
}
