#include "TransportSolver3D.h"
#include "DiamondDifference.h"
#include "TWDiamondDifference.h"

#include <iostream>
#include <iomanip>

int
main(int argc, char** argv)
{
  const unsigned int refinement = 4u;
  const unsigned int nl = 10u;
  const unsigned int nc = 10u;
  //----------------------------------------------------------------------------------------------------------
  // The first Kobayashi benchmark.
  {
    std::cout << "Building the First Kobayashi Benchmark.\n";
    std::cout << "---------------------------------------\n";
    std::vector<unsigned int> ix{5 * refinement, 4 * refinement, 2 * refinement, 4 * refinement, 5 * refinement};
    std::vector<unsigned int> iy{5 * refinement, 4 * refinement, 2 * refinement, 4 * refinement, 5 * refinement};
    std::vector<unsigned int> iz{5 * refinement, 4 * refinement, 2 * refinement, 4 * refinement, 5 * refinement};

    std::vector<double> dx{50.0, 40.0, 20.0, 40.0, 50.0};
    std::vector<double> dy{50.0, 40.0, 20.0, 40.0, 50.0};
    std::vector<double> dz{50.0, 40.0, 20.0, 40.0, 50.0};

    std::vector<unsigned int> blocks{2, 2, 2, 2, 2,
                                     2, 2, 2, 2, 2,
                                     2, 2, 2, 2, 2,
                                     2, 2, 2, 2, 2,
                                     2, 2, 2, 2, 2,

                                     2, 2, 2, 2, 2,
                                     2, 1, 1, 1, 2,
                                     2, 1, 1, 1, 2,
                                     2, 1, 1, 1, 2,
                                     2, 2, 2, 2, 2,

                                     2, 2, 2, 2, 2,
                                     2, 1, 1, 1, 2,
                                     2, 1, 0, 1, 2,
                                     2, 1, 1, 1, 2,
                                     2, 2, 2, 2, 2,

                                     2, 2, 2, 2, 2,
                                     2, 1, 1, 1, 2,
                                     2, 1, 1, 1, 2,
                                     2, 1, 1, 1, 2,
                                     2, 2, 2, 2, 2,

                                     2, 2, 2, 2, 2,
                                     2, 2, 2, 2, 2,
                                     2, 2, 2, 2, 2,
                                     2, 2, 2, 2, 2,
                                     2, 2, 2, 2, 2};

    std::array<BoundaryCondition, 6u> boundary_conditions = { BoundaryCondition::Vacuum, BoundaryCondition::Vacuum,
                                                              BoundaryCondition::Vacuum, BoundaryCondition::Vacuum,
                                                              BoundaryCondition::Vacuum, BoundaryCondition::Vacuum };

    auto mesh = BrickMesh3D(ix, iy, iz, dx, dy, dz, blocks, boundary_conditions);
    mesh.addPropsToBlock(0, 1e-1, 5e-2, 1.0);
    mesh.addPropsToBlock(1, 1e-4, 5e-5, 0.0);
    mesh.addPropsToBlock(2, 1e-1, 5e-2, 0.0);
    // mesh.addPropsToBlock(0, 1e-1, 0.0, 1.0);
    // mesh.addPropsToBlock(1, 1e-4, 0.0, 0.0);
    // mesh.addPropsToBlock(2, 1e-1, 0.0, 0.0);

    {
      std::cout << "Solving with Theta-Weighted Diamond Differences." << std::endl;
      std::cout << "---------------------------------------\n";
      TransportSolver3D<TWDiamondDifference> solver(mesh, nl, nc);
      if (solver.solve(1e-8, 1000u))
      {
        mesh.dumpToTextFile("Kobayashi_1_diamond");
        std::cout << "Solved with Theta-Weighted Diamond Differences." << std::endl;
      }
      std::cout << "---------------------------------------\n";
    }
  }
  //----------------------------------------------------------------------------------------------------------
}
