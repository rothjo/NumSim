#pragma once

#include "pressureSolver.h"

/**
 * Implementation of the Conjugate Gradient (CG) pressure solver.
 * Implements the solve method of the PressureSolver interface.
 */
class CG : public PressureSolver {
public:
    CG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    void solve() override;

private:
    /**
     * Set boundary values for the search direction d, to be called after avery iteration
     */
    void setBoundariesD();

    /**
     * get the discrete laplace operator applied to the pressure field variable at position (i, j)
     */
    double LaplaceP(int i, int j) const;

    /**
     * get the discrete laplace operator applied to the search direction field variable at position (i, j)
     */
    double LaplaceD(int i, int j) const;

    FieldVariable r_;
    FieldVariable d_;
    FieldVariable Ad_;

    double res_old2_;
    double res_new2_;
    double dx2_;
    double dy2_; 
};
