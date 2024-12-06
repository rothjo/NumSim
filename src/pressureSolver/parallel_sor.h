#pragma once

#include "sor.h"
#include "parallelPressureSolver.h"
#include "partitioning/partitioning.h"


/**
 * Implementation of the parallel SOR pressure solver.
 * Implements the solve method of the ParallelPressureSolver interface.
 */
class ParallelSOR : public ParallelPressureSolver {
public:
    // constructor that calls ParallelPressureSolver constructor
    ParallelSOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega, std::shared_ptr<Partitioning> partitioning);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     * Applying suitable communication for parallelization
     */
    void solve() override;

protected:
    double omega_;
};
  
