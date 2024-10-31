#include "output_writer/write_paraview_output.h"
#include "storage/array2d.h"
#include "storage/FieldVariable.h"
#include "discretization/discretization.h"
#include "discretization/centralDifferences.h"
#include "pressureSolver/pressureSolver.h"

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
    double epsilon = 0.1;
    int iteration = 10;
    // std::cout << size[0] << std::endl;
    FieldVariable myArray(size, origin, meshWidth);
    std::shared_ptr<Discretization> disc = std::make_shared<CentralDifferences>(nCells, meshWidth);
    //PressureSolver solver(disc, epsilon, iteration);



    for (int j = nCells[1]-1; j >= 0; --j) {
        for (int i = 0; i < nCells[0]; ++i) {
          if (i == 0 || i == nCells[0] - 1) {
            disc->p(i,j) = 0;
          } else if (j == 0 || j == nCells[1] - 1){
            disc->p(i,j) = 0;
          } else{
            disc->p(i,j) = 1;
          }
          std::cout << disc->p(i, j) << " ";
        }
        std::cout << std::endl;
    }

    for (int j = nCells[1]-1; j >= 0; --j) {
        for (int i = 0; i < nCells[0]; ++i) {
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
