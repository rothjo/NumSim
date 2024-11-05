#pragma once

#include "pressureSolver.h"


/**
 * Implementation of the Gauss-Seidel pressure solver.
 * Implements the solve method of the PressureSolver interface.
 */
class GaussSeidel : public PressureSolver {
public:
    GaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    void solve() override;
};
  
