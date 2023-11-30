#include "TransportSolver.h"

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

    auto mesh = BrickMesh3D(ix, iy, iz, dx, dy, dz, blocks);
    mesh.addPropsToBlock(0, 1e-1, 5e-2, 1.0);
    mesh.addPropsToBlock(1, 1e-4, 5e-5, 0.0);
    mesh.addPropsToBlock(2, 1e-1, 5e-2, 0.0);

    {
      std::cout << "Solving with Diamond Differences." << std::endl;
      std::cout << "---------------------------------------\n";
      TransportSolver solver(mesh, nl, nc, DiscretizationType::DiamondDifference);
      if (solver.solve(1e-8, 1000u))
      {
        mesh.dumpToTextFile("Kobayashi_1_diamond");

        std::cout << "Results for Diamond Differences." << std::endl;
        std::cout << "---------------------------------------\n";
        double one_a_fluxes[10u] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        mesh.fluxAtPoint(105.0, 105.0, 105.0, one_a_fluxes[0]);
        mesh.fluxAtPoint(105.0, 115.0, 105.0, one_a_fluxes[1]);
        mesh.fluxAtPoint(105.0, 125.0, 105.0, one_a_fluxes[2]);
        mesh.fluxAtPoint(105.0, 135.0, 105.0, one_a_fluxes[3]);
        mesh.fluxAtPoint(105.0, 145.0, 105.0, one_a_fluxes[4]);
        mesh.fluxAtPoint(105.0, 155.0, 105.0, one_a_fluxes[5]);
        mesh.fluxAtPoint(105.0, 165.0, 105.0, one_a_fluxes[6]);
        mesh.fluxAtPoint(105.0, 175.0, 105.0, one_a_fluxes[7]);
        mesh.fluxAtPoint(105.0, 185.0, 105.0, one_a_fluxes[8]);
        mesh.fluxAtPoint(105.0, 195.0, 105.0, one_a_fluxes[9]);
        std::cout << "1-A fluxes:\n" << std::setprecision(6);
        for (unsigned int i = 0u; i < 10; ++i)
          std::cout << one_a_fluxes[i] << "\n";

        double one_b_fluxes[10u] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        mesh.fluxAtPoint(105.0, 105.0, 105.0, one_b_fluxes[0]);
        mesh.fluxAtPoint(115.0, 115.0, 115.0, one_b_fluxes[1]);
        mesh.fluxAtPoint(125.0, 125.0, 125.0, one_b_fluxes[2]);
        mesh.fluxAtPoint(135.0, 135.0, 135.0, one_b_fluxes[3]);
        mesh.fluxAtPoint(145.0, 145.0, 145.0, one_b_fluxes[4]);
        mesh.fluxAtPoint(155.0, 155.0, 155.0, one_b_fluxes[5]);
        mesh.fluxAtPoint(165.0, 165.0, 165.0, one_b_fluxes[6]);
        mesh.fluxAtPoint(175.0, 175.0, 175.0, one_b_fluxes[7]);
        mesh.fluxAtPoint(185.0, 185.0, 185.0, one_b_fluxes[8]);
        mesh.fluxAtPoint(195.0, 195.0, 195.0, one_b_fluxes[9]);
        std::cout << "1-B fluxes:\n" << std::setprecision(6);
        for (unsigned int i = 0u; i < 10; ++i)
          std::cout << one_b_fluxes[i] << "\n";

        double one_c_fluxes[10u] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        mesh.fluxAtPoint(105.0, 155.0, 105.0, one_c_fluxes[0]);
        mesh.fluxAtPoint(115.0, 155.0, 105.0, one_c_fluxes[1]);
        mesh.fluxAtPoint(125.0, 155.0, 105.0, one_c_fluxes[2]);
        mesh.fluxAtPoint(135.0, 155.0, 105.0, one_c_fluxes[3]);
        mesh.fluxAtPoint(145.0, 155.0, 105.0, one_c_fluxes[4]);
        mesh.fluxAtPoint(155.0, 155.0, 105.0, one_c_fluxes[5]);
        mesh.fluxAtPoint(165.0, 155.0, 105.0, one_c_fluxes[6]);
        mesh.fluxAtPoint(175.0, 155.0, 105.0, one_c_fluxes[7]);
        mesh.fluxAtPoint(185.0, 155.0, 105.0, one_c_fluxes[8]);
        mesh.fluxAtPoint(195.0, 155.0, 105.0, one_c_fluxes[9]);
        std::cout << "1-C fluxes:\n" << std::setprecision(6);
        for (unsigned int i = 0u; i < 10; ++i)
          std::cout << one_c_fluxes[i] << "\n";
        std::cout << std::flush;
      }
      std::cout << "---------------------------------------\n";
    }
  }
  //----------------------------------------------------------------------------------------------------------
}
