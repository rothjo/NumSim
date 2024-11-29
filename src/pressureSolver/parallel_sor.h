#pragma once

#include "sor.h"
#include "parallelPressureSolver.h"
#include "partitioning/partitioning.h"


/**
 * Implementation of the parallel sor pressure solver.
 * Implements the solve method of the PressureSolver interface.
 */
class ParallelSOR : public ParallelPressureSolver {
public:
    // constructor that calls ParallelPressureSolver constructor
    ParallelSOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega, std::shared_ptr<Partitioning> partitioning);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    void solve() override;

protected:
    double omega_;
};
  
