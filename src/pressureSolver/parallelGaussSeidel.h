#pragma once

#include "gaussSeidel.h"
#include "parallelPressureSolver.h"
#include "partitioning/partitioning.h"


/**
 * Implementation of the parallel Gauss-Seidel pressure solver.
 * Implements the solve method of the ParallelPressureSolver interface.
 */
class ParallelGaussSeidel : public ParallelPressureSolver {
public:
    // constructor that calls ParallelPressureSolver constructor
    ParallelGaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     * Applying suitable communication for parallelization
     */
    void solve() override;
};
  
