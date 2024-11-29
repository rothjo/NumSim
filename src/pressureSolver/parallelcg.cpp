#include "parallelcg.h"
#include <numeric>
#include <iostream>

ParallelCG::ParallelCG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning)
    : ParallelPressureSolver(discretization, epsilon, maximumNumberOfIterations, partitioning) {}


void ParallelCG::solve() {
    const double eps2 = epsilon_ * epsilon_;

    const int gridWidth = discretization_->pIEnd() - discretization_->pIBegin();
    const int gridHeight = discretization_->pJEnd() - discretization_->pJBegin();

    // Initialize vectors for r and p directly in the grid
    double rsold = 0.0;

    // Compute initial residual r = b - A * p and initialize p = r
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); ++i) {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); ++j) {
            const double p_xx = (discretization_->p(i + 1, j) - 2.0 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / (discretization_->dx() * discretization_->dx());
            const double p_yy = (discretization_->p(i, j + 1) - 2.0 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / (discretization_->dy() * discretization_->dy());
            const double residual = discretization_->rhs(i, j) - (p_xx + p_yy);
            discretization_->p(i, j) = residual; // p = r initially
            rsold += residual * residual;
        }
    }
    rsold = partitioning_->globalSum(rsold);

    if (rsold < eps2) {
        return;
    }

    for (int k = 0; k < maximumNumberOfIterations_; ++k) {
        double ApSum = 0.0;

        // Compute Ap = A * p
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); ++i) {
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); ++j) {
                const double p_xx = (discretization_->p(i + 1, j) - 2.0 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / (discretization_->dx() * discretization_->dx());
                const double p_yy = (discretization_->p(i, j + 1) - 2.0 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / (discretization_->dy() * discretization_->dy());
                const double Ap = p_xx + p_yy;
                ApSum += Ap * discretization_->p(i, j);
            }
        }
        ApSum = partitioning_->globalSum(ApSum);

        // Compute alpha = rsold / (p^T * Ap)
        const double alpha = rsold / ApSum;

        double rsnew = 0.0;
        // Update pressure field and residual
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); ++i) {
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); ++j) {
                discretization_->p(i, j) += alpha * discretization_->p(i, j);
                const double residual = discretization_->rhs(i, j) - alpha * ApSum;
                rsnew += residual * residual;
            }
        }
        rsnew = partitioning_->globalSum(rsnew);

        if (rsnew < eps2) {
            break;
        }

        // Update p = r + (rsnew / rsold) * p
        const double beta = rsnew / rsold;
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); ++i) {
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); ++j) {
                discretization_->p(i, j) = beta * discretization_->p(i, j);
            }
        }
        rsold = rsnew;

        // Synchronize boundary values
        communicateAndBoundaries();
    }
}