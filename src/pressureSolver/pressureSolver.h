#pragma once

#include <memory>
#include <discretization/discretization.h>


/**
 * Interface for the pressure solver. It computes the pressure field variable such that the continuity equation is fulfilled.
 */
class PressureSolver {
public:
    PressureSolver(std::shared_ptr<Discretization>, double epsilon, int maximumNumberOfIterations);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    virtual void solve() = 0;

protected:

    /**
     * Set boundary values for the pressure, needs to be called after every iteration
     */
    void setBoundaryValues();

    /**
     * Compute residual norm for possible exit condition
     */
    virtual void computeResidualNorm();

    /**
     * Get the euclidean residual norm
     */
    double residualNorm();

    std::shared_ptr<Discretization> discretization_;
    double epsilon_;
    int maximumNumberOfIterations_;
    double residualNorm2_;
    int numberOfIterations_;






};
