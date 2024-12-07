#include "TransportSolver2D.h"

#include "DiamondDifference3D.h"
#include "TWDiamondDifference3D.h"
#include "TransportSolver3D.h"

#include "InputParameters.h"

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

  try
  {
    args.parse_args(argc, argv);
  }
  catch (const std::exception& err)
  {
    std::cerr << err.what() << std::endl;
    std::cerr << args;
    std::exit(1);
  }

  auto inp_path = std::filesystem::path(args.get<std::string>("--input_file"));
  const auto params = parseInputParameters(inp_path.string());

  if (params._num_dims == 1)
  {
    std::cerr << "1D transport solves are currently not supported." << std::endl;
    std::exit(1);
  }
  else if (params._num_dims == 2)
  {
    std::array<BoundaryCondition, 4u> boundary_conditions = { BoundaryCondition::Vacuum, BoundaryCondition::Vacuum,
                                                              BoundaryCondition::Vacuum, BoundaryCondition::Vacuum };

    auto mesh = BrickMesh2D(params._x_intervals,
                            params._y_intervals,
                            params._dx,
                            params._dy,
                            params._blocks,
                            boundary_conditions);

    for (auto & [block, props] : params._block_mat_info)
      mesh.addPropsToBlock(block, props);

    if (params._eq_type == EquationType::DD)
    {
      TWDDTransportSolver2D solver(mesh, params._num_e_groups, params._num_polar, params._num_azimuthal);
      if (solver.solveFixedSource(params._src_it_tol, params._num_src_it, params._gs_tol, params._num_mg_it))
        mesh.dumpToTextFile((inp_path.parent_path().string() / inp_path.stem()).string());
    }
    else if (params._eq_type == EquationType::TW_DD)
    {
      DDTransportSolver2D solver(mesh, params._num_e_groups, params._num_polar, params._num_azimuthal);
      if (solver.solveFixedSource(params._src_it_tol, params._num_src_it, params._gs_tol, params._num_mg_it))
        mesh.dumpToTextFile((inp_path.parent_path().string() / inp_path.stem()).string());
    }
  }
  else if (params._num_dims == 3)
  {
    if (params._num_e_groups > 1)
    {
      std::cerr << "At present the 3D transport solver only supports monoenergetic calculations." << std::endl;
      std::exit(1);
    }

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
                            boundary_conditions);

    for (auto & [block, props] : params._block_mat_info)
      mesh.addPropsToBlock(block, props._g_total[0], props._g_g_scatter_mat[0], props._g_src[0]);

    if (params._eq_type == EquationType::DD)
    {
      TransportSolver3D<DiamondDifference3D> solver(mesh, params._num_polar, params._num_azimuthal);
      if (solver.solve(params._src_it_tol, params._num_src_it))
        mesh.dumpToTextFile((inp_path.parent_path().string() / inp_path.stem()).string());
    }
    else if (params._eq_type == EquationType::TW_DD)
    {
      TransportSolver3D<TWDiamondDifference3D> solver(mesh, params._num_polar, params._num_azimuthal);
      if (solver.solve(params._src_it_tol, params._num_src_it))
        mesh.dumpToTextFile((inp_path.parent_path().string() / inp_path.stem()).string());
    }
  }
  else
  {
    std::cerr << "Invalid number of dimensions: " << params._num_dims << std::endl;
    std::exit(1);
  }
}
