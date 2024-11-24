#pragma once

#include "gaussSeidel.h"
#include "parallelPressureSolver.h"
#include "partitioning/partitioning.h"


/**
 * Implementation of the Gauss-Seidel pressure solver.
 * Implements the solve method of the PressureSolver interface.
 */
class ParallelGaussSeidel : public ParallelPressureSolver {
public:
    // constructor that calls ParallelPressureSolver constructor
    ParallelGaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    void solve() override;
};
  
