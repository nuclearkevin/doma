#include "InputParameters.h"

#include "TransportSolver1D.h"
#include "TransportSolver2D.h"
#include "TransportSolver3D.h"

#include "argparse/argparse.hpp"

#include <iostream>
#include <iomanip>
#include <filesystem>
#include <chrono>

#include "Parallel.h"

int
main(int argc, char** argv)
{
  argparse::ArgumentParser args("Discrete Ordinates Mini-App");
  args.add_argument("-i", "--input_file")
      .required()
      .help("The path of a DOMA input file.");
  args.add_argument("--verbose")
      .help("Whether to use verbose output or not.")
      .default_value(false)
      .implicit_value(true);
  args.add_argument("-n", "--n-threads")
      .help("The number of OpenMP threads to use. Default is 1.")
      .default_value(1)
      .scan<'d', int>();

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
  const int threads = args.get<int>("--n-threads");

  omp_set_num_threads(threads);

  auto start = std::chrono::high_resolution_clock::now();

  if (params._num_dims == 1)
  {
    std::array<BoundaryCondition, 2u> boundary_conditions = { BoundaryCondition::Vacuum, BoundaryCondition::Vacuum };

    auto mesh = BrickMesh1D(params, boundary_conditions, threads);
    mesh.validateProps();

    if (params._eq_type == EquationType::DD)
    {
      DDTransportSolver1D solver(mesh, params, args.get<bool>("--verbose"), threads);
      if (params._mode == RunMode::FixedSrc)
      {
        if (solver.solveFixedSource((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Fixed source solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
          std::exit(0);
        }
      }
      else if (params._mode == RunMode::Transient)
      {
        if (solver.solveTransient((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Transient solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
          std::exit(0);
        }
      }
      else if (params._mode == RunMode::Eigen)
      {
        if (solver.solveEigenvalue((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Eigenvalue solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
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
      std::cerr << "Theta-Weighted Diamond Differences are currently not supported in slab geometries!" << std::endl;
      std::exit(1);
    }
  }
  else if (params._num_dims == 2)
  {
    std::array<BoundaryCondition, 4u> boundary_conditions = { BoundaryCondition::Vacuum, BoundaryCondition::Vacuum,
                                                              BoundaryCondition::Vacuum, BoundaryCondition::Vacuum };

    auto mesh = BrickMesh2D(params, boundary_conditions, threads);
    mesh.validateProps();

    if (params._eq_type == EquationType::DD)
    {
      DDTransportSolver2D solver(mesh, params, args.get<bool>("--verbose"), threads);
      if (params._mode == RunMode::FixedSrc)
      {
        if (solver.solveFixedSource((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Fixed source solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
          std::exit(0);
        }
      }
      else if (params._mode == RunMode::Transient)
      {
        if (solver.solveTransient((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Transient solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
          std::exit(0);
        }
      }
      else if (params._mode == RunMode::Eigen)
      {
        if (solver.solveEigenvalue((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Eigenvalue solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
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
      TWDDTransportSolver2D solver(mesh, params, args.get<bool>("--verbose"), threads);
      if (params._mode == RunMode::FixedSrc)
      {
        if (solver.solveFixedSource((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Fixed source solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
          std::exit(0);
        }
      }
      else if (params._mode == RunMode::Transient)
      {
        if (solver.solveTransient((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Transient solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
          std::exit(0);
        }
      }
      else if (params._mode == RunMode::Eigen)
      {
        if (solver.solveEigenvalue((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Eigenvalue solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
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

    auto mesh = BrickMesh3D(params, boundary_conditions, threads);
    mesh.validateProps();

    if (params._eq_type == EquationType::DD)
    {
      DDTransportSolver3D solver(mesh, params, args.get<bool>("--verbose"), threads);
      if (params._mode == RunMode::FixedSrc)
      {
        if (solver.solveFixedSource((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Fixed source solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
          std::exit(0);
        }
      }
      else if (params._mode == RunMode::Transient)
      {
        if (solver.solveTransient((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Transient solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
          std::exit(0);
        }
      }
      else if (params._mode == RunMode::Eigen)
      {
        if (solver.solveEigenvalue((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Eigenvalue solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
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
      TWDDTransportSolver3D solver(mesh, params, args.get<bool>("--verbose"), threads);
      if (params._mode == RunMode::FixedSrc)
      {
        if (solver.solveFixedSource((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Fixed source solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
          std::exit(0);
        }
      }
      else if (params._mode == RunMode::Transient)
      {
        if (solver.solveTransient((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Transient solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
          std::exit(0);
        }
      }
      else if (params._mode == RunMode::Eigen)
      {
        if (solver.solveEigenvalue((inp_path.parent_path().string() / inp_path.stem()).string()))
        {
          std::cout << "Eigenvalue solver finished executing." << std::endl;
          auto stop = std::chrono::high_resolution_clock::now();
          std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
                    << " (s)." << std::endl;
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
  else
  {
    std::cerr << "Invalid number of dimensions: " << params._num_dims << std::endl;
    std::exit(1);
  }
}
