#include "TransportSolver3D.h"
#include "DiamondDifference.h"
#include "TWDiamondDifference.h"

#include "InputParameters.h"

#include "argparse/argparse.hpp"

#include <iostream>
#include <iomanip>

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

  auto inp_file = args.get<std::string>("--input_file");
  auto params = parseInputParameters(inp_file);

  if (params._num_dims == 1)
  {
    std::cerr << "1D transport solves are currently not supported." << std::endl;
    std::exit(1);
  }
  else if (params._num_dims == 2)
  {
    std::cerr << "2D transport solves are currently not supported." << std::endl;
    std::exit(1);
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
                            boundary_conditions);

    for (auto & [block, props] : params._block_mat_info)
      mesh.addPropsToBlock(block, props._g_total[0], props._g_g_scatter_mat[0], props._g_src[0]);

    if (params._eq_type == EquationType::DD)
    {
      TransportSolver3D<DiamondDifference> solver(mesh, params._num_polar, params._num_azimuthal);
      if (solver.solve(params._src_it_tol, params._num_src_it))
        mesh.dumpToTextFile(inp_file.substr(0, inp_file.find_first_of(".")));
    }
    else if (params._eq_type == EquationType::TW_DD)
    {
      TransportSolver3D<TWDiamondDifference> solver(mesh, params._num_polar, params._num_azimuthal);
      if (solver.solve(params._src_it_tol, params._num_src_it))
        mesh.dumpToTextFile(inp_file.substr(0, inp_file.find_first_of(".")));
    }
  }
  else
  {
    std::cerr << "Invalid number of dimensions: " << params._num_dims << std::endl;
    std::exit(1);
  }
}
