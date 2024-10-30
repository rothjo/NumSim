#pragma once
#include <memory>
#include <discretization/discretization.h>


class PressureSolver {
public:
    PressureSolver(std::shared_ptr<Discretization>, double epsilon, int maximumNumberOfIterations);

    /**
     * Solve poisson problem for the pressure
     */
    virtual void solve() = 0;

protected:

    /**
     * Set boundary values for the pressure
     */
    void setBoundaryValues();

    std::shared_ptr<Discretization> discretization_;
    double epsilon_;
    int maximumNumberOfIterations_;






}
