#include "GCAngularQuadrature.h"

#include <iostream>
#include <cmath>

//------------------------------------------------------------------------------
// Gauss-Legendre-Chebyshev product quadrature here.
//------------------------------------------------------------------------------
GCAngularQuadrature::GCAngularQuadrature(unsigned int n_c, unsigned int n_l, unsigned int num_dims)
  : _n_c(std::move(n_c)),
    _n_l(std::move(n_l)),
    _n_d(std::move(num_dims)),
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
      std::cout << mu << " " << i << std::endl;
      const auto & omega = _azimuthal_quadrature.angularRoot(j - 1u);

      // Positive octant.
      unordered_mus.emplace_back(mu);
      unordered_etas.emplace_back(std::sqrt(1.0 - (mu * mu)) * std::cos(omega));
      unordered_xis.emplace_back(std::sqrt(1.0 - (mu * mu)) * std::sin(omega));
      unordered_weights.emplace_back(weight);
      // Negative octant if in 3D.
      if (_n_d > 2)
      {
        unordered_mus.emplace_back(-1.0 * mu);
        unordered_etas.emplace_back(-1.0 * std::sqrt(1.0 - (mu * mu)) * std::cos(omega));
        unordered_xis.emplace_back(-1.0 * std::sqrt(1.0 - (mu * mu)) * std::sin(omega));
        unordered_weights.emplace_back(weight);
      }
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
