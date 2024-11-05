#pragma once

#include "pressureSolver.h"



class GaussSeidel : public PressureSolver {
public:
    GaussSeidel(std::shared_ptr<Discretization>, double epsilon, int maximumNumberOfIterations);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    void solve() override;
};
  
