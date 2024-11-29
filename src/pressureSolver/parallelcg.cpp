#include "parallelCG.h"
#include <cmath>
#include <vector>
#include <iostream>

ParallelCG::ParallelCG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning) :
    ParallelPressureSolver(discretization, epsilon, maximumNumberOfIterations, partitioning) {}

void ParallelCG::solve() {
    const double eps2 = epsilon_ * epsilon_;
    const int size = discretization_->nCellsGlobal()[0] * discretization_->nCellsGlobal()[1];

    // Initialize vectors
    std::vector<double> r(size), p(size), Ap(size);

    // Compute initial residual r = b - A * p
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); ++i) {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); ++j) {
            const double p_xx = (discretization_->p(i + 1, j) - 2.0 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / (discretization_->dx() * discretization_->dx());
            const double p_yy = (discretization_->p(i, j + 1) - 2.0 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / (discretization_->dy() * discretization_->dy());
            r[i * discretization_->nCellsGlobal()[1] + j] = discretization_->rhs(i, j) - (p_xx + p_yy);
        }
    }
    communicateAndBoundaries();

    // Copy r into p
    std::copy(r.begin(), r.end(), p.begin());

    double rsold = partitioning_->globalSum(std::inner_product(r.begin(), r.end(), r.begin(), 0.0));
    if (rsold < eps2) {
        return;
    }

    for (int k = 0; k < maximumNumberOfIterations_; ++k) {
        // Compute Ap
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); ++i) {
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); ++j) {
                const double p_xx = (p[(i + 1) * discretization_->nCellsGlobal()[1] + j] - 2.0 * p[i * discretization_->nCellsGlobal()[1] + j] + p[(i - 1) * discretization_->nCellsGlobal()[1] + j]) / (discretization_->dx() * discretization_->dx());
                const double p_yy = (p[i * discretization_->nCellsGlobal()[1] + (j + 1)] - 2.0 * p[i * discretization_->nCellsGlobal()[1] + j] + p[i * discretization_->nCellsGlobal()[1] + (j - 1)]) / (discretization_->dy() * discretization_->dy());
                Ap[i * discretization_->nCellsGlobal()[1] + j] = p_xx + p_yy;
            }
        }
        communicateAndBoundaries();

        // Alpha = rsold / (p^T * Ap)
        double alpha = rsold / partitioning_->globalSum(std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0));

        // Update pressure field and residual
        for (int i = 0; i < size; ++i) {
            discretization_->p(i) += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rsnew = partitioning_->globalSum(std::inner_product(r.begin(), r.end(), r.begin(), 0.0));
        if (rsnew < eps2) {
            break;
        }

        // Update p = r + (rsnew / rsold) * p
        double beta = rsnew / rsold;
        for (int i = 0; i < size; ++i) {
            p[i] = r[i] + beta * p[i];
        }
        rsold = rsnew;
    }
}
