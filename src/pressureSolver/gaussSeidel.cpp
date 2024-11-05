#include "gaussSeidel.h"
#include <iostream>	


GaussSeidel::GaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations) {}

void GaussSeidel::solve() {
    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();
    const double k = (dx2 * dy2) / (2.0 * (dx2 + dy2));
    const double eps2 = epsilon_ * epsilon_;

    int iteration = 0;
    // applyBoundaryValues(); already set at t = 0
    computeResidualNorm();

    while (residualNorm2_ > eps2 && iteration < maximumNumberOfIterations_) {
        ++iteration;
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
                double px = (discretization_->p(i + 1, j) + discretization_->p(i - 1, j)) / dx2;
                double py = (discretization_->p(i, j + 1) + discretization_->p(i, j - 1)) / dy2;
                discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
            }
        }
        setBoundaryValues();
        computeResidualNorm();

        
    }
    this->numberOfIterations_ = iteration;

}
    
