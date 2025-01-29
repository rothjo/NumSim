// #include "output_writer/write_paraview_output.h"
#include "storage/array2d.h"
#include "storage/FieldVariable.h"
#include "discretization/discretization.h"
#include "discretization/centralDifferences.h"
#include "pressureSolver/pressureSolver.h"
#include "pressureSolver/gaussSeidel.h"
#include "pressureSolver/parallelPressureSolver.h"
#include "pressureSolver/parallelGaussSeidel.h"
#include "computation/parallelComputation.h"

#include "computation/computation.h"
#include <iostream>
#include <cstdlib>
#include <memory>
#include <mpi.h>
#include "partitioning/partitioning.h"
#include <array>
#include <unistd.h>
#include <chrono>

int main(int argc, char *argv[])
{
  // sleep(15);
  // Computation comp;
  // comp.initialize(argc, argv);
  // comp.runSimulation();

  


  
  MPI_Init(&argc, &argv);
  // MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  ParallelComputation parallelcomp;
  parallelcomp.initialize(argc, argv);

  double start_time = MPI_Wtime();

  parallelcomp.runSimulation();

  double end_time = MPI_Wtime();
  double elapsed_time = end_time - start_time;
  

  // Rank 0 prints the elapsed time
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
      std::cout << "Total elapsed time: " << elapsed_time << " seconds" << std::endl;
  }

  MPI_Finalize();

  return EXIT_SUCCESS;
}
