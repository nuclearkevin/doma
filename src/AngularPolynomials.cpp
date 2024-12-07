#include "AngularPolynomials.h"

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
