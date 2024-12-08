#include "BrickMesh1D.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>

BrickMesh1D::BrickMesh1D(const std::vector<unsigned int> & nx, const std::vector<double> & dx,
                         const std::vector<unsigned int> & blocks, const std::array<BoundaryCondition, 2u> & bcs,
                         const std::unordered_map<unsigned int, MaterialProps> & props)
  : _nx(nx),
    _tot_num_x(0u),
    _num_groups(1u),
    _dx(dx),
    _blocks(blocks),
    _num_cells(0u),
    _total_length(0.0),
    _bcs(bcs),
    _block_mat_info(props)
{
  if (_nx.size() != _blocks.size())
  {
    std::cerr << "Error: The blocks vector is not equal to the number of subdivisions!" << std::endl;
    std::exit(1);
  }

  for (auto nx : _nx)
    _tot_num_x += nx;

  std::vector<double> tot_dx;
  std::vector<unsigned int> tot_blocks;

  for (unsigned int i = 0u; i < _nx.size(); ++i)
  {
    const auto d_x = _dx[i] / static_cast<double>(_nx[i]);
    for (unsigned int j = 0u; j < _nx[i]; ++j)
      tot_dx.emplace_back(d_x);
  }

  // Unpack the material block array and expand it to include all subdivisions.
  for (unsigned int i = 0u; i < _nx.size(); ++i)
  {
    for (unsigned int ii = 0u; ii < _nx[i]; ++ii)
      tot_blocks.emplace_back(_blocks[i]);
  }

  // Build the mesh.
  double x_c = 0.5 * tot_dx[0];
  for (unsigned int i = 0u; i < tot_dx.size(); ++i)
  {
    _cells.emplace_back(tot_dx[i], x_c, _num_cells, tot_blocks[i], this);
    _total_length += _cells.back().getLength();
    _num_cells++;
    x_c += tot_dx[i];
  }

  // Link each cell to it's neighbors and handle boundary conditions.
  for (unsigned int i = 0u; i < tot_dx.size(); ++i)
  {
    auto & cell = _cells[i];

    // Neighbors and boundary cells in x.
    //----------------------------------------------------------------------------------------------------------
    // Left neighbor is a boundary condition when i = 0, otherwise we have a neighboring cell.
    if (i == 0u)
    {
      cell.addNeighbor(nullptr, CertesianFaceSide::Left);
      _boundary_cells[static_cast<unsigned int>(CertesianFaceSide::Left)].emplace_back(&cell);
    }
    else
      cell.addNeighbor(&_cells[i - 1u], CertesianFaceSide::Left);

    // Right neighbor is a boundary condition when i = NX - 1, otherwise we have a neighboring cell.
    if (i == tot_dx.size() - 1u)
    {
      cell.addNeighbor(nullptr, CertesianFaceSide::Right);
      _boundary_cells[static_cast<unsigned int>(CertesianFaceSide::Right)].emplace_back(&cell);
    }
    else
      cell.addNeighbor(&_cells[i + 1u], CertesianFaceSide::Right);
  }
}

void
BrickMesh1D::validateProps()
{
  for (const auto & cell : _cells)
  {
    if (_block_mat_info.count(cell._block_id) == 0)
    {
      std::cout << "Block " << cell._block_id << " (defined on the mesh) has no material properties!" << std::endl;
      std::exit(1);
    }
  }
}

// A function to initialize the group-wise scalar fluxes in each cell.
void
BrickMesh1D::initFluxes(unsigned int num_groups)
{
  _num_groups = num_groups;
  for (auto & cell : _cells)
    cell.initFluxes(num_groups);
}

// Returns true if the point exists on the mesh, false if it does not. The flux at that point
// will be stored in 'returned_flux' if the point is on the mesh.
bool
BrickMesh1D::fluxAtPoint(const double & x, unsigned int g, double & returned_flux) const
{
  for (auto & cell : _cells)
  {
    if (cell.pointInCell(x))
    {
      returned_flux = cell._total_scalar_flux[g];
      return true;
    }
  }
  return false;
}

// Dump the flux to a text file.
void
BrickMesh1D::dumpToTextFile(const std::string & file_name)
{
  std::ofstream dims(file_name + "_dims.txt", std::ofstream::out);
  std::ofstream blocks(file_name + "_blocks.txt", std::ofstream::out);
  std::ofstream x(file_name + "_meshx.txt", std::ofstream::out);
  x << std::setprecision(6);
  dims << "num_x: " << _tot_num_x << std::endl;
  dims << "num_g: " << _num_groups << std::endl;
  dims.close();

  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    std::ofstream flux(file_name + "_g" + std::to_string(g) + "_flux.txt", std::ofstream::out);

    flux << std::setprecision(6);
    for (const auto & cell : _cells)
    {
      flux << cell._total_scalar_flux[g] << std::endl;

      if (g == 0u)
      {
        x << cell._x_c << std::endl;
        blocks << cell._block_id << std::endl;
      }
    }
    flux.close();
  }
  blocks.close();
  x.close();
}

void
BrickMesh1D::printAllBlocks()
{
  for (unsigned int i = 0u; i < _tot_num_x; ++i)
  {
    std::cout << _cells[i]._block_id;
    std::cout << " ";
  }
  std::cout << std::flush;
}

void
BrickMesh1D::printAllCoords()
{
  std::cout << std::setprecision(6);
  for (unsigned int i = 0u; i < _tot_num_x; ++i)
  {
    const auto & cell = _cells[i];
    std::cout << "(" << cell._x_c << ")";
    std::cout << " ";
  }
  std::cout << std::flush;
}

void
BrickMesh1D::printAllSideLengths()
{
  std::cout << std::setprecision(6);
  for (unsigned int i = 0u; i < _tot_num_x; ++i)
  {
    const auto & cell = _cells[i];
    std::cout << "(" << cell.getLength() << ")";
    std::cout << " ";
  }
  std::cout << std::flush;
}
