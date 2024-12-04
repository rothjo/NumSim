#pragma once

#include <memory>
#include "discretization/discretization.h"
#include "partitioning/partitioning.h"
#include "pressureSolver.h"
#include "mpi.h"


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
     * Set boundary values for the pressure, needs to be called after every iteration
     * Also, send necessary information to neighbouring ranks
     */
    void communicateAndBoundaries();
    
    /**
     * compute residualNorm
     */
    void computeResidualNorm() override;

protected:
    std::shared_ptr<Partitioning> partitioning_;
};
