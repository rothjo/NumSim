#include "pressureSolver/pressureSolver.h"
#include <cmath>
#include "discretization/discretization.h"
#include <iostream>


PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations) :
    discretization_(discretization),
    epsilon_(epsilon),
    maximumNumberOfIterations_(maximumNumberOfIterations) {}
    
void PressureSolver::setBoundaryValues() {
    // Set pressure boundary values for the upper and lower side of the grid
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        discretization_->p(i, discretization_->pJBegin() - 1) = discretization_->p(i, discretization_->pJBegin());
        discretization_->p(i, discretization_->pJEnd()) = discretization_->p(i, discretization_->pJEnd() - 1);
    }

    // Set pressure boundary values for the left and right side of the grid
    for (int j = discretization_->pJBegin() - 1; j < discretization_->pJEnd() + 1; j++) {
        discretization_->p(discretization_->pIBegin() - 1, j) = discretization_->p(discretization_->pIBegin(), j);
        discretization_->p(discretization_->pIEnd(), j) = discretization_->p(discretization_->pIEnd() - 1, j);
    }
}

void PressureSolver::computeResidualNorm() {
    // Define needed parameters
    double residualNorm2 = 0.0;
    const double dx_2 = discretization_->dx() * discretization_->dx();
    const double dy_2 = discretization_->dy() * discretization_->dy();
    const int N = (discretization_->pIEnd() - discretization_->pIBegin()) * (discretization_->pJEnd() - discretization_->pJBegin()); // Number of points used to compute norm

    // Compute l^2 norm of residual
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            const double p_xx = (discretization_->p(i + 1, j) - 2.0 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / dx_2;
            const double p_yy = (discretization_->p(i, j + 1) - 2.0 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / dy_2;
            residualNorm2 += (discretization_->rhs(i, j) - p_xx - p_yy) * (discretization_->rhs(i, j) - p_xx - p_yy);
        }
    }
    residualNorm2_ = residualNorm2 / N;
}   

double PressureSolver::residualNorm() {
    return residualNorm2_;
}

double PressureSolver::laplaceP(int i, int j) const {
    const double dx_2 = discretization_->dx() * discretization_->dx();
    const double dy_2 = discretization_->dy() * discretization_->dy();
    return ((discretization_->p(i + 1, j) - 2.0 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / dx_2) + ((discretization_->p(i, j + 1) - 2.0 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / dy_2);
}
int PressureSolver::numberOfIterations() {
    return numberOfIterations_;
}