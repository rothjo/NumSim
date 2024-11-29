#pragma once

#include "parallelPressureSolver.h"

/**
 * Implementation of the Conjugate Gradient (CG) pressure solver.
 * Implements the solve method of the PressureSolver interface.
 */
class ParallelCG : public ParallelPressureSolver {
public:
    ParallelCG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    void solve() override;
};
