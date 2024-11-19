// #include "output_writer/write_paraview_output.h"
// #include "storage/array2d.h"
// #include "storage/FieldVariable.h"
// #include "discretization/discretization.h"
// #include "discretization/centralDifferences.h"
// #include "pressureSolver/pressureSolver.h"
// #include "pressureSolver/gaussSeidel.h"

#include "computation/computation.h"
#include <iostream>
#include <cstdlib>
#include <memory>
#include <chrono>

int main(int argc, char *argv[])
{

  // std::unique_ptr<PressureSolver> solver = std::make_unique<GaussSeidel>(disc, epsilon, iteration);
  // solver->solver();

  // for (int j = disc->pJEnd(); j >= 0; --j) {
  //     for (int i = 0; i < disc->pIEnd()+1; ++i) {
  //       std::cout << disc->p(i, j) << " ";
  //     }
  //     std::cout << std::endl;
  // }

  //   // double interpolate = myArray.interpolateAt(0.15, 0.15);
  //   // std::cout << interpolate << std::endl; 



  // // // write 5 output files
  // // for (int i = 0; i < 5; i++)
  // // {
  // //   writeParaviewOutput(i);
  // // }

  Computation comp;
  comp.initialize(argc, argv);

  // first measures
  auto start = std::chrono::high_resolution_clock::now();

  comp.runSimulation();

  auto stop = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  // int n = 5;

  // for(int i = 0; i<n-1; i++) {
  //   start = std::chrono::high_resolution_clock::now();

  //   comp.runSimulation();

  //   stop = std::chrono::high_resolution_clock::now();

  //   duration += std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  // }

  std::cout << "Time needed: " << duration.count() << " ms" << std::endl;



  return EXIT_SUCCESS;
}
