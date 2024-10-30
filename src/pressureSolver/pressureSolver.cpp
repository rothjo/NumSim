#include "pressureSolver/pressureSolver.h"
#include <cmath>


/**
 * Implementation of the pressure solver. 
 */
PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations) :
    
    discretization_(discretization),
    epsilon_(epsilon),
    maximumNumberOfIterations_(maximumNumberOfIterations) {


}

void PressureSolver::setBoundaryValues() {

}

/**
 * Compute euclidean residual norm
 */
void PressureSolver::computeResidualNorm() {
    // Define needed parameters
    double residualNorm2 = 0.0;
    const double dx_2 = pow(discretization_->dx(), 2);
    const double dy_2= pow(discretization_->dy(), 2);
    const int N = discretization_->nCells()[0] * discretization_->nCells()[1]; // Number of points used to compute norm

    // Compute l^2 norm by taking into account rhs and u_xx, u_yy

    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {

            double p_xx = (discretization_->p(i + 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / dx_2;
            double p_yy = (discretization_->p(i, j + 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / dy_2;
            residualNorm2 += pow(discretization_->rhs(i, j) - p_xx - p_yy, 2);
        }
    }

    residualNorm2_ = residualNorm2 / N;
}


double PressureSolver::residualNorm() {
    return residualNorm2_;
}
