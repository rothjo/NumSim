#include "sor.h"
#include <iostream>

SOR::SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), omega_(omega) {}

void SOR::solve() {
    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();
    const double k = (dx2 * dy2) / (2.0 * (dx2 + dy2));
    const double eps2 = epsilon_ * epsilon_;

    int iteration = 0;

    computeResidualNorm();

    while (residualNorm2_ > eps2 && iteration < maximumNumberOfIterations_) {
        ++iteration;
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
                const double px = (discretization_->p(i + 1, j) + discretization_->p(i - 1, j)) / dx2;
                const double py = (discretization_->p(i, j + 1) + discretization_->p(i, j - 1)) / dy2;

                const double correction = k * (px + py - discretization_->rhs(i, j)) - discretization_->p(i, j);
                discretization_->p(i,j) += omega_ * correction;
            }
        }
        setBoundaryValues();
        computeResidualNorm();
    }
    this->numberOfIterations_ += iteration;
}
    
