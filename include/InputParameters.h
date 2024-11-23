#pragma once

#include <array>
#include <vector>

#include "TransportBase.h"

struct InputParameters
{
  unsigned int _num_dims;

  std::vector<unsigned int> _x_intervals;
  std::vector<unsigned int> _y_intervals;
  std::vector<unsigned int> _z_intervals;

  std::vector<double> _dx;
  std::vector<double> _dy;
  std::vector<double> _dz;

  std::vector<unsigned int> _blocks;

  std::array<BoundaryCondition, 6u> _bcs;
};
