#pragma once

#include <memory>
#include "../discretization/discretization.h"
#include "../partitioning/partitioning.h"
#include "pressureSolver.h"


/**
 * Interface for the pressure solver. It computes the pressure field variable such that the continuity equation is fulfilled.
 */
class ParallelPressureSolver : public PressureSolver {
public:
    ParallelPressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations,
                           std::shared_ptr<Partitioning> partitioning);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    virtual void solve() = 0;

protected:

    /**
     * Set boundary values at the bottom boundary
     */
    void setBoundaryValuesBottom();

    /**
     * Set boundary values at the top boundary
     */
    void setBoundaryValuesTop();

    /**
     * Set boundary values at the left boundary
     */
    void setBoundaryValuesLeft();

    /**
     * Set boundary values at the right boundary
     */
    void setBoundaryValuesRight();
    
    /**
     * compute residualNorm
     */
    void computeResidualNorm() override;

protected:
    std::shared_ptr<Partitioning> partitioning_;
};
