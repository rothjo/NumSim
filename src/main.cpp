#include "output_writer/write_paraview_output.h"
#include "storage/array2d.h"
#include "storage/FieldVariable.h"
#include "discretization/discretization.h"
#include "discretization/centralDifferences.h"
#include "pressureSolver/pressureSolver.h"
#include "pressureSolver/gaussSeidel.h"

#include <iostream>
#include <cstdlib>
#include <memory>

int main(int argc, char *argv[])
{

  // Initialize the array with size 3x3
    std::array<int, 2> size = {12, 12};
    std::array<double, 2> origin = {-0.05, -0.05};
    std::array<double, 2> meshWidth = {0.1, 0.1};
    std::array<int, 2> nCells = {6, 6};
    double epsilon = 1e-5;
    int iteration = 1e4;
    // std::cout << size[0] << std::endl;
    FieldVariable myArray(size, origin, meshWidth);
    std::shared_ptr<Discretization> disc = std::make_shared<CentralDifferences>(nCells, meshWidth);
    GaussSeidel solver(disc, epsilon, iteration);

    // Test pressure boundaries

    for (int j = disc->pJEnd(); j >= 0; --j) {
        for (int i = 0; i < disc->pIEnd()+1; ++i) {
          if (i == 0 || i == disc->pIEnd() ) {
            disc->p(i,j) = 0;
          } else if (j == 0 || j == disc->pJEnd()){
            disc->p(i,j) = 0;
          } else{
            disc->p(i,j) = 1;
          }
          std::cout << disc->p(i, j) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    
    for (int j = disc->pJBegin(); j < disc->pJEnd(); j++) {
        disc->p(0, j) = disc->p(1, j);
        disc->p(disc->pIEnd(), j) = disc->p(disc->pIEnd() - 1, j);
    }

    // Set pressure boundary values for the upper and lower side of the grid
    for (int i = disc->pIBegin(); i < disc->pIEnd(); i++) {
        disc->p(i, 0) = disc->p(i, 1);
        disc->p(i, disc->pJEnd()) = disc->p(i, disc->pJEnd() - 1);
    }

    for (int j = disc->pJEnd(); j >= 0; --j) {
        for (int i = 0; i < disc->pIEnd()+1; ++i) {
          std::cout << disc->p(i, j) << " ";
          disc->rhs(i,j) = -0.001*disc->p(i,j);
        }
        std::cout << std::endl;
    }

    for (int j = disc->pJEnd(); j >= 0; --j) {
        for (int i = 0; i < disc->pIEnd()+1; ++i) {
          std::cout << disc->p(i, j) << " ";
          disc->p(i,j) = 0;
        }
        std::cout << std::endl;
    }

    for (int j = disc->pJEnd(); j >= 0; --j) {
        for (int i = 0; i < disc->pIEnd()+1; ++i) {
          std::cout << disc->rhs(i, j) << " ";
        }
        std::cout << std::endl;
    }

    // TEST SOLVER

    solver.solve();

    for (int j = disc->pJEnd(); j >= 0; --j) {
        for (int i = 0; i < disc->pIEnd()+1; ++i) {
          std::cout << disc->p(i, j) << " ";
        }
        std::cout << std::endl;
    }

    // double interpolate = myArray.interpolateAt(0.15, 0.15);
    // std::cout << interpolate << std::endl; 



  // // write 5 output files
  // for (int i = 0; i < 5; i++)
  // {
  //   writeParaviewOutput(i);
  // }

  return EXIT_SUCCESS;
}
