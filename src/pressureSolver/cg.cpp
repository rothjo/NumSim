#include "cg.h"
#include <iostream>

CG::CG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations)
    : PressureSolver(discretization, epsilon, maximumNumberOfIterations),
    dx2_(discretization->dx() * discretization->dx()),
    dy2_(discretization->dy() * discretization->dy()),
    r_(FieldVariable({discretization->nCells()[0] + 3, discretization_->nCells()[1] + 3}, {-1.5 * discretization_->dx(), -1.5 * discretization_->dy()}, discretization_->meshWidth())),
    d_(FieldVariable({discretization->nCells()[0] + 3, discretization_->nCells()[1] + 3}, {-1.5 * discretization_->dx(), -1.5 * discretization_->dy()}, discretization_->meshWidth())),
    Ad_(FieldVariable({discretization->nCells()[0] + 3, discretization_->nCells()[1] + 3}, {-1.5 * discretization_->dx(), -1.5 * discretization_->dy()}, discretization_->meshWidth()))
    {}


void CG::solve() {
    // Initialize constants and loop limits
    const double eps2 = epsilon_ * epsilon_;
    const int N = (*discretization_).nCells()[0] * (*discretization_).nCells()[1];
    const double N_eps2 = N * eps2;

    int pIBegin = discretization_->pIBegin();
    int pIEnd = discretization_->pIEnd();
    int pJBegin = discretization_->pJBegin();
    int pJEnd = discretization_->pJEnd();

    // Initialize first time step
    res_old2_ = 0.0;
    const double M_inv = 1.0 / (2.0 / dx2_ + 2.0 / dy2_); // Diagonal element approximation, Jacobi
    // Compute initial residual r and apply preconditioner, i.e. solve Mz_0 = r_0
    for (int i = pIBegin; i < pIEnd; i++) {
        for (int j = pJBegin; j < pJEnd; j++) {
            double res_ij = discretization_->rhs(i, j) - laplaceP(i, j);
            r_(i, j) = res_ij;

            // Apply preconditioner (Jacobi: divide by diagonal element)
            d_(i, j) = res_ij * M_inv;

            res_old2_ += r_(i, j) * d_(i, j); // Compute preconditioned residual
        }
    }
    setBoundariesD();

    if (res_old2_ < N_eps2) {
        return;
    }

    numberOfIterations_ = maximumNumberOfIterations_;
    double dAd_ = 0.0;

    // Time loop
    for (int k = 0; k < maximumNumberOfIterations_; ++k) {
        dAd_ = 0.0;

        // Compute Ad and partial dAd_ locally
        for (int i = pIBegin; i < pIEnd; i++) {
            for (int j = pJBegin; j < pJEnd; j++) {
                const double laplace_d = LaplaceD(i, j);
                Ad_(i, j) = laplace_d;
                dAd_ += d_(i, j) * laplace_d;
            }
        }
        // Compute alpha
        double alpha = res_old2_ / dAd_;

        res_new2_ = 0.0;

        // Update p, r and compute local res_new2_
        for (int i = pIBegin; i < pIEnd; i++) {
            for (int j = pJBegin; j < pJEnd; j++) {
                (*discretization_).p(i, j) += alpha * d_(i, j);
                r_(i, j) -= alpha * Ad_(i, j);
                // Apply preconditioner to update residual
                double z_ij = r_(i, j) * M_inv;
                res_new2_ += r_(i, j) * z_ij;
                d_(i, j) = z_ij + (res_new2_ / res_old2_) * d_(i, j);
            }
        }


        // Check if new residuum is lower than tolerance
        if (res_new2_ < N_eps2) {
            numberOfIterations_ = k;
            std::cout << "Converged after " << k << " iterations." << std::endl;
            break;
        }

        res_old2_ = res_new2_;
        setBoundariesD();
    }

    setBoundaryValues();
}

void CG::setBoundariesD() {
    // Set d boundary values for the upper and lower side of the grid
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        d_(i, discretization_->pJBegin() - 1) = d_(i, discretization_->pJBegin());
        d_(i, discretization_->pJEnd()) = d_(i, discretization_->pJEnd() - 1);
    }

    // Set d boundary values for the left and right side of the grid
    for (int j = discretization_->pJBegin() - 1; j < discretization_->pJEnd() + 1; j++) {
        d_(discretization_->pIBegin() - 1, j) = d_(discretization_->pIBegin(), j);
        d_(discretization_->pIEnd(), j) = d_(discretization_->pIEnd() - 1, j);
    }

}



double CG::LaplaceD(int i, int j) const {
    return ((d_(i + 1, j) - 2.0 * d_(i, j) + d_(i - 1, j)) / dx2_) + ((d_(i, j + 1) - 2.0 * d_(i, j) + d_(i, j - 1)) / dy2_);
}