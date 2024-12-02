#pragma once

#include "pressureSolver.h"

/**
 * Implementation of the Conjugate Gradient (CG) pressure solver.
 * Implements the solve method of the PressureSolver interface.
 */
class CG : public PressureSolver {
public:
    CG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    void solve() override;

private:
    void communicateAndBoundariesD();

    double LaplaceP(int i, int j) const;

    double LaplaceD(int i, int j) const;

    // alpha, beta, res_old, res_new, d, Ad,
    FieldVariable r_;
    FieldVariable d_;
    FieldVariable Ad_;
    double res_old2_;
    double res_new2_;
    double dx2_;
    double dy2_; 
    
};
