#include "InputParameters.h"

#include "TransportSolver1D.h"
#include "TransportSolver2D.h"
#include "TransportSolver3D.h"

#include "argparse/argparse.hpp"

#include <iostream>
#include <iomanip>
#include <filesystem>

int
main(int argc, char** argv)
{
  argparse::ArgumentParser args("Discrete Ordinates MiniApp");
  args.add_argument("-i", "--input_file")
      .required()
      .help("The path of a DOMA input file.");
  args.add_argument("--verbose")
      .help("Whether to use verbose output or not.")
      .default_value(false)
      .implicit_value(true);

  try
  {
    args.parse_args(argc, argv);
  }
  catch (const std::exception & err)
  {
    std::cerr << err.what() << std::endl;
    std::cerr << args;
    std::exit(1);
  }

  auto inp_path = std::filesystem::path(args.get<std::string>("--input_file"));
  const auto params = parseInputParameters(inp_path.string());

  if (params._num_dims == 1)
  {
    std::array<BoundaryCondition, 2u> boundary_conditions = { BoundaryCondition::Vacuum, BoundaryCondition::Vacuum };

    auto mesh = BrickMesh1D(params._x_intervals,
                            params._dx,
                            params._blocks,
                            boundary_conditions,
                            params._block_mat_info);
    mesh.validateProps();

    if (params._eq_type == EquationType::DD)
    {
      DDTransportSolver1D solver(mesh, params._num_e_groups, params._num_polar);
      if (solver.solveFixedSource(params._src_it_tol, params._num_src_it, params._gs_tol, params._num_mg_it))
        mesh.dumpToTextFile((inp_path.parent_path().string() / inp_path.stem()).string());
    }
    else if (params._eq_type == EquationType::TW_DD)
    {
      std::cerr << "At present 1D calculations do not support theta-weighted diamond differences!" << std::endl;
      std::exit(1);
    }
  }
  else if (params._num_dims == 2)
  {
    std::array<BoundaryCondition, 4u> boundary_conditions = { BoundaryCondition::Vacuum, BoundaryCondition::Vacuum,
                                                              BoundaryCondition::Vacuum, BoundaryCondition::Vacuum };

    auto mesh = BrickMesh2D(params, boundary_conditions);
    mesh.validateProps();

    if (params._eq_type == EquationType::DD)
    {
      DDTransportSolver2D solver(mesh, params, args.get<bool>("--verbose"));
      if (params._mode == RunMode::FixedSrc)
      {
        if (solver.solveFixedSource((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Fixed source solver finished executing." << std::endl;
          std::exit(0);
        }
      }
      else if (params._mode == RunMode::Transient)
      {
        if (solver.solveTransient((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Transient solver finished executing." << std::endl;
          std::exit(0);
        }
      }
      else
      {
        std::cerr << "Unsupported execution mode!" << std::endl;
        std::exit(1);
      }
    }
    else if (params._eq_type == EquationType::TW_DD)
    {
      TWDDTransportSolver2D solver(mesh, params, args.get<bool>("--verbose"));
      if (params._mode == RunMode::FixedSrc)
      {
        if (solver.solveFixedSource((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Fixed source solver finished executing." << std::endl;
          std::exit(0);
        }
      }
      else if (params._mode == RunMode::Transient)
      {
        if (solver.solveTransient((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Transient solver finished executing." << std::endl;
          std::exit(0);
        }
      }
      else
      {
        std::cerr << "Unsupported execution mode!" << std::endl;
        std::exit(1);
      }
    }
  }
  else if (params._num_dims == 3)
  {
    std::array<BoundaryCondition, 6u> boundary_conditions = { BoundaryCondition::Vacuum, BoundaryCondition::Vacuum,
                                                              BoundaryCondition::Vacuum, BoundaryCondition::Vacuum,
                                                              BoundaryCondition::Vacuum, BoundaryCondition::Vacuum };

    auto mesh = BrickMesh3D(params._x_intervals,
                            params._y_intervals,
                            params._z_intervals,
                            params._dx,
                            params._dy,
                            params._dz,
                            params._blocks,
                            boundary_conditions,
                            params._block_mat_info);
    mesh.validateProps();

    if (params._eq_type == EquationType::DD)
    {
      DDTransportSolver3D solver(mesh, params._num_e_groups, params._num_polar, params._num_azimuthal);
      if (solver.solveFixedSource(params._src_it_tol, params._num_src_it, params._gs_tol, params._num_mg_it))
        mesh.dumpToTextFile((inp_path.parent_path().string() / inp_path.stem()).string());
    }
    else if (params._eq_type == EquationType::TW_DD)
    {
      TWDDTransportSolver3D solver(mesh, params._num_e_groups, params._num_polar, params._num_azimuthal);
      if (solver.solveFixedSource(params._src_it_tol, params._num_src_it, params._gs_tol, params._num_mg_it))
        mesh.dumpToTextFile((inp_path.parent_path().string() / inp_path.stem()).string());
    }
  }
  else
  {
    std::cerr << "Invalid number of dimensions: " << params._num_dims << std::endl;
    std::exit(1);
  }
}
