#include "pressureSolver/parallel_sor.h"

#include <cmath>

/**
 * Implementation of the red-black solver, a parallelisized version of the SOR solver.
 * @param discretization pointer to the implementation of the discretization
 * @param epsilon error tolerance below which we consider the solver to be converged
 * @param maximumNumberOfIterations when this number is reached, the solver stops without converging
 * @param partitioning information about subdomain
 */

ParallelSOR::ParallelSOR(std::shared_ptr<Discretization> discretization,
                         double epsilon,
                         int maximumNumberOfIterations,
                         double omega,
                         std::shared_ptr<Partitioning> partitioning) :
        ParallelPressureSolver(discretization, epsilon, maximumNumberOfIterations, partitioning),
        omega_(omega) {

}

/**
 * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
 */
void ParallelSOR::solve() {
    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();
    const double k = (dx2 * dy2) / (2.0 * (dx2 + dy2));
    const double eps2 = epsilon_ * epsilon_;

    int iteration = 0;
    // computeResidualNorm();
    residualNorm2_ = 1.0;
    while (residualNorm2_ > eps2 && iteration < maximumNumberOfIterations_) {
        ++iteration;
        
        // Initalize chessboard pattern, assure we start partition with the right pattern
        if (partitioning_->nodeOffsetSum() % 2 == 0) {

            // Go through all cells beginning in the bottom left corner
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
                for (int i = discretization_-> pIBegin() + (j - discretization_->pJBegin()) % 2; i < discretization_->pIEnd(); i += 2) {
                    const double px = (discretization_->p(i + 1, j) + discretization_->p(i - 1, j)) / dx2;
                    const double py = (discretization_->p(i, j + 1) + discretization_->p(i, j - 1)) / dy2;
                    const double correction = k * (px + py - discretization_->rhs(i, j)) - discretization_->p(i, j);
                    discretization_->p(i,j) += omega_ * correction;
                }
            }
            communicateAndBoundaries();

            // Go through all cells beginning one right to the bottom left corner
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
                for (int i = discretization_-> pIBegin() + (j - discretization_->pJBegin() + 1) % 2; i < discretization_->pIEnd(); i += 2) {
                    const double px = (discretization_->p(i + 1, j) + discretization_->p(i - 1, j)) / dx2;
                    const double py = (discretization_->p(i, j + 1) + discretization_->p(i, j - 1)) / dy2;
                    const double correction = k * (px + py - discretization_->rhs(i, j)) - discretization_->p(i, j);
                    discretization_->p(i,j) += omega_ * correction;
                }
            }
        } else {

            // Go through all cells beginning one right to the bottom left corner
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
                for (int i = discretization_-> pIBegin() + (j - discretization_->pJBegin() + 1) % 2; i < discretization_->pIEnd(); i += 2) {
                    const double px = (discretization_->p(i + 1, j) + discretization_->p(i - 1, j)) / dx2;
                    const double py = (discretization_->p(i, j + 1) + discretization_->p(i, j - 1)) / dy2;
                    const double correction = k * (px + py - discretization_->rhs(i, j)) - discretization_->p(i, j);
                    discretization_->p(i,j) += omega_ * correction;
                }
            }
            communicateAndBoundaries();

            // Go through all cells beginning in the bottom left corner
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
                for (int i = discretization_-> pIBegin() + (j - discretization_->pJBegin()) % 2; i < discretization_->pIEnd(); i += 2) {
                    const double px = (discretization_->p(i + 1, j) + discretization_->p(i - 1, j)) / dx2;
                    const double py = (discretization_->p(i, j + 1) + discretization_->p(i, j - 1)) / dy2;
                    const double correction = k * (px + py - discretization_->rhs(i, j)) - discretization_->p(i, j);
                    discretization_->p(i,j) += omega_ * correction;
                }
            } 
        }  

        communicateAndBoundaries();
        // Compute res only every 20 steps to decrease communication time
        if (iteration % 20 == 0) {
            computeResidualNorm();
        }
    }
    this->numberOfIterations_ = iteration;
}