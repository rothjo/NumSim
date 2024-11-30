#pragma once

#include "pressureSolver.h"


/**
 * Implementation of the CG pressure solver.
 * Implements the solve method of the PressureSolver interface.
 */
class CG : public PressureSolver {
public:
    CG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    void solve() override;
};
  
