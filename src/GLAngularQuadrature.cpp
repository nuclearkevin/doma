#include "GLAngularQuadrature.h"

GLAngularQuadrature::GLAngularQuadrature(unsigned int n_l)
  : _n_l(std::move(n_l)),
    _polar_quadrature(std::move(n_l))
{

}

double
GLAngularQuadrature::direction(unsigned int n) const
{
  return _polar_quadrature.root(n);
}

const double &
GLAngularQuadrature::weight(unsigned int n) const
{
  return _polar_quadrature.weight(n);
}
