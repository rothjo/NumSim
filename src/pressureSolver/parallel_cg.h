#pragma once

#include "parallelPressureSolver.h"

/**
 * Implementation of the parallel Conjugate Gradient (CG) pressure solver with Jacobi precondition.
 * Implements the solve method of the ParallelPressureSolver interface.
 */
class ParallelCG : public ParallelPressureSolver {
public:
    ParallelCG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     * See Wikipedia CG with preconditioning
     */
    void solve() override;

private:
    /**
     * Communicate boundaries and ghost layers for the search direction d
     * Can be derived out of Neumann restriction for p at time (n+1) and using CG iteration rule
     * p^(n)_0,j + alpha^(n)*d^(n)_0,j = p^(n)_1,j + alpha^(n)*d^(n)_1,j 
     * Since we have Neumann for p^(n) and alpha is const, we use the same Neumann BC for d
     */
    void communicateAndBoundariesD();

    /**
     * Compute Laplace operator for P
     */
    double LaplaceP(int i, int j) const;

    /**
     * Compute Laplace operator for P
     */
    double LaplaceD(int i, int j) const;

    // alpha, beta, res_old, res_new, d, Ad, dx^2, dy^2
    FieldVariable r_;
    FieldVariable d_;
    FieldVariable Ad_;
    double res_old2_;
    double res_new2_;
    double dx2_;
    double dy2_; 
    
};
