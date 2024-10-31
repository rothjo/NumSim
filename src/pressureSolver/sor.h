#pragma once

#include "pressureSolver.h"

/**
 * Implementation of the sor pressure solver. 
 */
class SOR : public PressureSolver {
public:
    SOR(std::shared_ptr<Discretization>, double epsilon, int maximumNumberOfIterations, double omega);
    
    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    void solve() override;

private:
    double omega_;
    void computeResidualNorm() override;
};
  
