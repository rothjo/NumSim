#include "parallelGaussSeidel.h"
#include <iostream>	


ParallelGaussSeidel::ParallelGaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning) :
    ParallelPressureSolver(discretization, epsilon, maximumNumberOfIterations, partitioning) {}

void ParallelGaussSeidel::solve() {
    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();
    const double k = (dx2 * dy2) / (2.0 * (dx2 + dy2));
    const double eps2 = epsilon_ * epsilon_;

    int iteration = 0;
    // applyBoundaryValues(); already set at t = 0
    computeResidualNorm();

    while (residualNorm2_ > eps2 && iteration < maximumNumberOfIterations_) {
        ++iteration;
        
        if (partitioning_->nodeOffsetSum() % 2 == 0) {

            // Go through all cells beginning in the bottom left corner
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
                for (int i = discretization_-> pIBegin() + (j - discretization_->pJBegin()) % 2; i < discretization_->pIEnd(); i += 2) {
                    double px = (discretization_->p(i + 1, j) + discretization_->p(i - 1, j)) / dx2;
                    double py = (discretization_->p(i, j + 1) + discretization_->p(i, j - 1)) / dy2;
                    discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
                }
            }
            communicateAndBoundaries();

            // Go through all cells beginning one right to the bottom left corner
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
                for (int i = discretization_-> pIBegin() + (j - discretization_->pJBegin() + 1) % 2; i < discretization_->pIEnd(); i += 2) {
                    double px = (discretization_->p(i + 1, j) + discretization_->p(i - 1, j)) / dx2;
                    double py = (discretization_->p(i, j + 1) + discretization_->p(i, j - 1)) / dy2;
                    discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
                }
            }
        } else {

            // Go through all cells beginning one right to the bottom left corner
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
                for (int i = discretization_-> pIBegin() + (j - discretization_->pJBegin() + 1) % 2; i < discretization_->pIEnd(); i += 2) {
                    double px = (discretization_->p(i + 1, j) + discretization_->p(i - 1, j)) / dx2;
                    double py = (discretization_->p(i, j + 1) + discretization_->p(i, j - 1)) / dy2;
                    discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
                }
            }
            communicateAndBoundaries();

            // Go through all cells beginning in the bottom left corner
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
                for (int i = discretization_-> pIBegin() + (j - discretization_->pJBegin()) % 2; i < discretization_->pIEnd(); i += 2) {
                    double px = (discretization_->p(i + 1, j) + discretization_->p(i - 1, j)) / dx2;
                    double py = (discretization_->p(i, j + 1) + discretization_->p(i, j - 1)) / dy2;
                    discretization_->p(i, j) = k * (px + py - discretization_->rhs(i, j));
                }
            } 
        }  

        communicateAndBoundaries();
        computeResidualNorm(); 
    }
    this->numberOfIterations_ = iteration;
}
    