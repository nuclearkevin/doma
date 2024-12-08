#include "BrickMesh3D.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>

BrickMesh3D::BrickMesh3D(const std::vector<unsigned int> & nx, const std::vector<unsigned int> & ny, const std::vector<unsigned int> & nz,
                         const std::vector<double> & dx, const std::vector<double> & dy, const std::vector<double> & dz,
                         const std::vector<unsigned int> & blocks, const std::array<BoundaryCondition, 6u> & bcs,
                         const std::unordered_map<unsigned int, MaterialProps> & props)
  : _nx(nx),
    _ny(ny),
    _nz(nz),
    _tot_num_x(0u),
    _tot_num_y(0u),
    _tot_num_z(0u),
    _dx(dx),
    _dy(dy),
    _dz(dz),
    _blocks(blocks),
    _num_cells(0u),
    _total_volume(0.0),
    _bcs(bcs),
    _block_mat_info(props)
{
  if (_nx.size() * _ny.size() * _nz.size() != _blocks.size())
  {
    std::cout << "Error: The blocks vector is not equal to the number of subdivisions!" << std::endl;
    exit(1);
  }

  for (auto nx : _nx)
    _tot_num_x += nx;
  for (auto n_y : _ny)
    _tot_num_y += n_y;
  for (auto n_z : _nz)
    _tot_num_z += n_z;

  std::vector<double> tot_dx;
  std::vector<double> tot_dy;
  std::vector<double> tot_dz;
  std::vector<unsigned int> tot_blocks;

  for (unsigned int i = 0u; i < _nx.size(); ++i)
  {
    const auto d_x = _dx[i] / static_cast<double>(_nx[i]);
    for (unsigned int j = 0u; j < _nx[i]; ++j)
      tot_dx.emplace_back(d_x);
  }
  for (unsigned int i = 0u; i < _ny.size(); ++i)
  {
    const auto d_y = _dy[i] / static_cast<double>(_ny[i]);
    for (unsigned int j = 0u; j < _ny[i]; ++j)
      tot_dy.emplace_back(d_y);
  }
  for (unsigned int i = 0u; i < _nz.size(); ++i)
  {
    const auto d_z = _dz[i] / static_cast<double>(_nz[i]);
    for (unsigned int j = 0u; j < _nz[i]; ++j)
      tot_dz.emplace_back(d_z);
  }

  // Unpack the material block array and expand it to include all subdivisions.
  for (unsigned int k = 0u; k < _nz.size(); ++k)
  {
    for (unsigned int kk = 0u; kk < _nz[k]; ++kk)
    {
      for (unsigned int j = 0u; j < _ny.size(); ++j)
      {
        for (unsigned int jj = 0u; jj < _ny[j]; ++jj)
        {
          for (unsigned int i = 0u; i < _nx.size(); ++i)
          {
            for (unsigned int ii = 0u; ii < _nx[i]; ++ii)
              tot_blocks.emplace_back(_blocks[k * _ny.size() * _nx.size() + j * _nx.size() + i]);
          }
        }
      }
    }
  }

  // Build the mesh.
  double x_c = 0.5 * tot_dx[0];
  double y_c = 0.5 * tot_dy[0];
  double z_c = 0.5 * tot_dz[0];
  for (unsigned int k = 0u; k < tot_dz.size(); ++k)
  {
    for (unsigned int j = 0u; j < tot_dy.size(); ++j)
    {
      for (unsigned int i = 0u; i < tot_dx.size(); ++i)
      {
        _cells.emplace_back(tot_dx[i], tot_dy[j], tot_dz[k], x_c, y_c, z_c, _num_cells, tot_blocks[k * tot_dy.size() * tot_dx.size() + j * tot_dx.size() + i], this);
        _total_volume += _cells.back().getVolume();
        _num_cells++;
        x_c += tot_dx[i];
      }
      y_c += tot_dy[j];
      x_c = 0.5 * tot_dx[0];
    }
    z_c += tot_dz[k];
    y_c = 0.5 * tot_dy[0];
  }

  // Link each cell to it's neighbors and handle boundary conditions.
  for (unsigned int k = 0u; k < tot_dz.size(); ++k)
  {
    for (unsigned int j = 0u; j < tot_dy.size(); ++j)
    {
      for (unsigned int i = 0u; i < tot_dx.size(); ++i)
      {
        auto & cell = _cells[k * _tot_num_y * _tot_num_x + j * _tot_num_x + i];

        // Neighbors and boundary cells in x.
        //----------------------------------------------------------------------------------------------------------
        // Left neighbor is a boundary condition when i = 0, otherwise we have a neighboring cell.
        if (i == 0u)
        {
          cell.addNeighbor(nullptr, CertesianFaceSide::Left);
          _boundary_cells[static_cast<unsigned int>(CertesianFaceSide::Left)].emplace_back(&cell);
        }
        else
          cell.addNeighbor(&_cells[k * _tot_num_y * _tot_num_x + j * _tot_num_x + i - 1u], CertesianFaceSide::Left);

        // Right neighbor is a boundary condition when i = NX - 1, otherwise we have a neighboring cell.
        if (i == tot_dx.size() - 1u)
        {
          cell.addNeighbor(nullptr, CertesianFaceSide::Right);
          _boundary_cells[static_cast<unsigned int>(CertesianFaceSide::Right)].emplace_back(&cell);
        }
        else
          cell.addNeighbor(&_cells[k * _tot_num_y * _tot_num_x + j * _tot_num_x + i + 1u], CertesianFaceSide::Right);

        // Neighbors and boundary cells in y.
        //----------------------------------------------------------------------------------------------------------
        // Back neighbor is a boundary condition when j = 0, otherwise we have a neighboring cell.
        if (j == 0u)
        {
          cell.addNeighbor(nullptr, CertesianFaceSide::Back);
          _boundary_cells[static_cast<unsigned int>(CertesianFaceSide::Back)].emplace_back(&cell);
        }
        else
          cell.addNeighbor(&_cells[k * _tot_num_y * _tot_num_x + (j - 1u) * _tot_num_x + i], CertesianFaceSide::Back);

        // Front neighbor is a boundary condition when j = NY - 1, otherwise we have a neighboring cell.
        if (j == tot_dy.size() - 1u)
        {
          cell.addNeighbor(nullptr, CertesianFaceSide::Front);
          _boundary_cells[static_cast<unsigned int>(CertesianFaceSide::Front)].emplace_back(&cell);
        }
        else
          cell.addNeighbor(&_cells[k * _tot_num_y * _tot_num_x + (j + 1u) * _tot_num_x + i], CertesianFaceSide::Front);

        // Neighbors and boundary cells in z.
        //----------------------------------------------------------------------------------------------------------
        // Bottom neighbor is a boundary condition when k = 0, otherwise we have a neighboring cell.
        if (k == 0u)
        {
          cell.addNeighbor(nullptr, CertesianFaceSide::Bottom);
          _boundary_cells[static_cast<unsigned int>(CertesianFaceSide::Bottom)].emplace_back(&cell);
        }
        else
          cell.addNeighbor(&_cells[(k - 1u) * _tot_num_y * _tot_num_x + j * _tot_num_x + i], CertesianFaceSide::Bottom);

        // Top neighbor is a boundary condition when k = NZ - 1, otherwise we have a neighboring cell.
        if (k == tot_dz.size() - 1u)
        {
          cell.addNeighbor(nullptr, CertesianFaceSide::Top);
          _boundary_cells[static_cast<unsigned int>(CertesianFaceSide::Top)].emplace_back(&cell);
        }
        else
          cell.addNeighbor(&_cells[(k + 1u) * _tot_num_y * _tot_num_x + j * _tot_num_x + i], CertesianFaceSide::Top);
      }
    }
  }
}

void
BrickMesh3D::validateProps()
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

void
BrickMesh3D::initFluxes(unsigned int num_groups)
{
  _num_groups = num_groups;
  for (auto & cell : _cells)
    cell.initFluxes(num_groups);
}

// Returns true if the point exists on the mesh, false if it does not. The flux at that point
// will be stored in 'returned_flux' if the point is on the mesh.
bool
BrickMesh3D::fluxAtPoint(const double & x, const double & y, const double & z, unsigned int g, double & returned_flux) const
{
  for (auto & cell : _cells)
  {
    if (cell.pointInCell(x, y, z))
    {
      returned_flux = cell._total_scalar_flux[g];
      return true;
    }
  }
  return false;
}

// Dump the flux to a text file.
void
BrickMesh3D::dumpToTextFile(const std::string & file_name)
{
  std::ofstream dims(file_name + "_dims.txt", std::ofstream::out);
  std::ofstream blocks(file_name + "_blocks.txt", std::ofstream::out);
  std::ofstream x(file_name + "_meshx.txt", std::ofstream::out);
  x << std::setprecision(6);
  std::ofstream y(file_name + "_meshy.txt", std::ofstream::out);
  y << std::setprecision(6);
  std::ofstream z(file_name + "_meshz.txt", std::ofstream::out);
  z << std::setprecision(6);
  dims << "num_x: " << _tot_num_x << std::endl;
  dims << "num_y: " << _tot_num_y << std::endl;
  dims << "num_z: " << _tot_num_z << std::endl;
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
        y << cell._y_c << std::endl;
        z << cell._z_c << std::endl;
        blocks << cell._block_id << std::endl;
      }
    }
    flux.close();
  }
  blocks.close();
  x.close();
  y.close();
  z.close();
}

void
BrickMesh3D::printAllBlocks()
{
  for (unsigned int k = 0u; k < _tot_num_z; ++k)
  {
    for (unsigned int j = 0u; j < _tot_num_y; ++j)
    {
      for (unsigned int i = 0u; i < _tot_num_x; ++i)
      {
        std::cout << _cells[k * _tot_num_y * _tot_num_x + j * _tot_num_x + i]._block_id;
        std::cout << " ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  std::cout << std::flush;
}

void
BrickMesh3D::printAllCoords()
{
  std::cout << std::setprecision(6);
  for (unsigned int k = 0u; k < _tot_num_z; ++k)
  {
    for (unsigned int j = 0u; j < _tot_num_y; ++j)
    {
      for (unsigned int i = 0u; i < _tot_num_x; ++i)
      {
        const auto & cell = _cells[k * _tot_num_y * _tot_num_x + j * _tot_num_x + i];
        std::cout << "(" << cell._x_c << ", " << cell._y_c << ", " << cell._z_c << ")";
        std::cout << " ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  std::cout << std::flush;
}

void
BrickMesh3D::printAllSideLengths()
{
  std::cout << std::setprecision(6);
  for (unsigned int k = 0u; k < _tot_num_z; ++k)
  {
    for (unsigned int j = 0u; j < _tot_num_y; ++j)
    {
      for (unsigned int i = 0u; i < _tot_num_x; ++i)
      {
        const auto & cell = _cells[k * _tot_num_y * _tot_num_x + j * _tot_num_x + i];
        std::cout << "(" << cell.getSideLength(CertesianAxis::X) << ", " << cell.getSideLength(CertesianAxis::Y) << ", " << cell.getSideLength(CertesianAxis::Z) << ")";
        std::cout << " ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  std::cout << std::flush;
}
