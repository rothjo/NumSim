#pragma once

#include "gaussSeidel.h"
#include "parallelPressureSolver.h"


/**
 * Implementation of the Gauss-Seidel pressure solver.
 * Implements the solve method of the PressureSolver interface.
 */
class ParallelGaussSeidel : public ParallelPressureSolver {
public:
    // constructor that calls ParallelPressureSolver constructor
    ParallelGaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    void solve() override;

protected:
    /**
     * communicates values between subdomains after each half iteration
     * TODO: optimize by only sending half the pressure data (initially send all to simplify)
     */
    void pressureCommunication();
};
  
