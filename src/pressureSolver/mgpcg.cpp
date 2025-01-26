#include "mgpcg.h"
#include <cmath>
#include <iostream>

// Constructor
MGPCG::MGPCG(std::shared_ptr<Discretization> baseDiscretization, double epsilon,
                    int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning)
    : CG(baseDiscretization, epsilon, maximumNumberOfIterations),
      partitioning_(partitioning),
      baseRHS_(FieldVariable({baseDiscretization->nCells()[0] + 3, baseDiscretization->nCells()[1] + 3}, {-1.5 * baseDiscretization->meshWidth()[0], -1.5 * baseDiscretization->meshWidth()[1]}, baseDiscretization->meshWidth())),
      baseP_(FieldVariable({baseDiscretization->nCells()[0] + 3, baseDiscretization->nCells()[1] + 3}, {-1.5 * baseDiscretization->meshWidth()[0], -1.5 * baseDiscretization->meshWidth()[1]}, baseDiscretization->meshWidth())),
      z_(FieldVariable({baseDiscretization->nCells()[0] + 3, baseDiscretization->nCells()[1] + 3}, {-1.5 * baseDiscretization->dx(), -1.5 * baseDiscretization->dy()}, baseDiscretization->meshWidth())),
      levels_(std::log2(baseDiscretization->nCells()[0])) {
    assert((baseDiscretization->nCells()[0] == baseDiscretization->nCells()[1] && baseDiscretization->nCells()[0] > 0 &&  (baseDiscretization->nCells()[0] & (baseDiscretization->nCells()[0] - 1)) == 0) &&
           "Number of cells in each direction must be a power of 2.");
    grids_.push_back(baseDiscretization);

    // Build hierarchy of grids
    for (int l = 1; l < levels_; ++l) {
        grids_.push_back(coarsenGrid(l-1));
    }

    for (size_t l = 0; l < grids_.size() - 1; ++l) { // Exclude the last grid
        // const auto& grid = grids_[l];
        smoothers_.emplace_back(std::make_shared<GaussSeidel>(grids_[l], epsilon, 2));
    }
    smoothers_.emplace_back(std::make_shared<GaussSeidel>(grids_[grids_.size() - 1], epsilon, 4));

}

// Solve method
void MGPCG::solve() {
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
    // std::cout << discretization_->rhs(6, 5) << std::endl;
    // const double M_inv = 1.0 / (2.0 / dx2_ + 2.0 / dy2_); // Diagonal element approximation, Jacobi
    // Compute initial residual r and apply preconditioner, i.e. solve Mz_0 = r_0
    for (int i = pIBegin; i < pIEnd; i++) {
        for (int j = pJBegin; j < pJEnd; j++) {
            double res_ij = discretization_->rhs(i, j) - laplaceP(i, j);
            r_(i, j) = res_ij;
            discretization_->rhs(i, j) = r_(i, j);
            res_old2_ += r_(i, j) * r_(i, j); 

        }
    }
    // std::cout << discretization_->rhs(6, 5) << std::endl;

    // Perform a multigrid V-cycle
    // std::cout << discretization_->rhs(5,5) << std::endl;
    vCycle();

    // for (int i = pIBegin; i < pIEnd; i++) {
    //     for (int j = pJBegin; j < pJEnd; j++) {
    //         z_(i,j) = M_inv * r_(i,j);
    //     }
    // }

    // std::cout << "sigde3" << std::endl;

    double alpha_top = 0.0;

    for (int i = pIBegin; i < pIEnd; i++) {
        for (int j = pJBegin; j < pJEnd; j++) {
            d_(i, j) = z_(i, j);
            // std::cout << "is this zero?" <<d_(i,j) << std::endl;

            alpha_top += r_(i, j) * z_(i, j); // Compute preconditioned residual

            // discretization_->rhs(i, j) = baseRHS_(i, j);
            // discretization_->p(i, j) = baseP_(i, j);
        }
    }
    // std::cout << "sigde4" << std::endl;
    setBoundariesD();

    // std::cout << res_old2_ << std::endl;
    // std::cout << res_old2_ << std::endl;
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
        double alpha = alpha_top / dAd_;

        res_new2_ = 0.0;

        // std::cout << "sigde5" << std::endl;

         // Update p, r and compute local res_new2_
        for (int i = pIBegin; i < pIEnd; i++) {
            for (int j = pJBegin; j < pJEnd; j++) {
                (*discretization_).p(i, j) += alpha * d_(i, j);
                r_(i, j) -= alpha * Ad_(i, j);
                (*discretization_).rhs(i, j) = r_(i, j);
                res_new2_ += r_(i, j) * r_(i, j);
            }
        }

        // Precondition multigrid cycle
        vCycle();
        // for (int i = pIBegin; i < pIEnd; i++) {
        //     for (int j = pJBegin; j < pJEnd; j++) {
        //         z_(i,j) = M_inv * r_(i,j);
        //     }
        // }

        // std::cout << "sigde6" << std::endl;

        double beta_top = 0.0;

        // Update p, r and compute local res_new2_
        for (int i = pIBegin; i < pIEnd; i++) {
            for (int j = pJBegin; j < pJEnd; j++) {
                // z_(i, j) = discretization_->p(i, j); // Leave out for jacobi
                // discretization_->p(i, j) = baseP_(i, j); // Leave out for jacobi
                // std::cout << baseP_(i, j) << alpha_top <<std::endl;
                beta_top += r_(i, j) * z_(i, j);
            }
        }

        for (int i = pIBegin; i < pIEnd; i++) {
            for (int j = pJBegin; j < pJEnd; j++) {
                d_(i, j) = z_(i, j) + (beta_top / alpha_top) * d_(i, j);
            }
        }

        alpha_top = beta_top;



        // std::cout << "sigde7" << std::endl;


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

// Compute residual (void version)
void MGPCG::computeResidual(int level) {
    for (int i = grids_[level]->pIBegin(); i < grids_[level]->pIEnd(); ++i) {
        for (int j = grids_[level]->pJBegin(); j < grids_[level]->pJEnd(); ++j) {
            const double dx_2 = grids_[level]->dx() * grids_[level]->dx();
            const double dy_2 = grids_[level]->dy() * grids_[level]->dy();
            double laplaceP = ((grids_[level]->p(i + 1, j) - 2.0 * grids_[level]->p(i, j) + grids_[level]->p(i - 1, j)) / dx_2) + ((grids_[level]->p(i, j + 1) - 2.0 * grids_[level]->p(i, j) + grids_[level]->p(i, j - 1)) / dy_2);
            grids_[level]->rhs(i, j) = grids_[level]->rhs(i, j) - laplaceP;
        }
    }
}


void MGPCG::vCycle() {
    // Perform a MGPCG V-cycle
    // computeResidual(0);

    for (int i = grids_[0]->pIBegin(); i < grids_[0]->pIEnd(); ++i) {
        for (int j = grids_[0]->pJBegin(); j < grids_[0]->pJEnd(); ++j) {
            baseRHS_(i, j) = discretization_->rhs(i, j);
            baseP_(i, j) = discretization_->p(i, j);
            discretization_->p(i, j) = 0.0;
        }
    }

    // std::cout << discretization_->rhs(6, 5) << std::endl;

    // Go down the V-Cycle: Restrict residuals to coarser grids
    for (int l = 0; l < levels_ - 1; ++l) {
        smoothers_[l]->solve();  // Pre-smoothing
        computeResidual(l);
        restrictToCoarserGrid(l);
    }

    // Solve on the coarsest grid
    smoothers_[levels_ - 1]->solve();

    // Go up the V-Cycle: Prolong errors_ to finer grids
    for (int l = levels_ - 2; l >= 0; --l) {
        prolongAndCorrectToFinerGrid(l+1);
        smoothers_[l]->solve();  // Post-smoothing

    }

    for (int i = grids_[0]->pIBegin(); i < grids_[0]->pIEnd(); ++i) {
        for (int j = grids_[0]->pJBegin(); j < grids_[0]->pJEnd(); ++j) {
            // std::cout << discretization_->p(i, j) << std::endl;
            z_(i, j) = discretization_->p(i, j);
            discretization_->rhs(i, j) = baseRHS_(i, j);
            discretization_->p(i, j) = baseP_(i, j);
        }
    }

    
    // std::cout << discretization_->rhs(5, 5) << std::endl;

    // for (int i = discretization_->pIBegin(); discretization_->pIEnd(); ++i) {
    //     for (int j = discretization_->pJBegin(); discretization_->pJEnd(); ++j) {
    //         std::cout << "i = " << i << ", j = " << j << std::endl;
    //         discretization_->rhs(i, j) = baseRHS_(i, j);
    //     }
    // }

    // full_solver->solve();

    // std::cout << "MGPCG reached maximum iterations without convergence.\n";
}

// Coarsen the grid
std::shared_ptr<Discretization> MGPCG::coarsenGrid(int fineLevel) {
    int nx = grids_[fineLevel]->nCells()[0] / 2;
    double dx = grids_[fineLevel]->dx() * 2;

    //! TODO: change 0.5 parameter to settings_.omega
    return std::make_shared<DonorCell>(std::array<int, 2>{nx, nx}, std::array<double, 2>{dx, dx}, partitioning_, 0.5);
}
// Restrict residual to coarser grid (void version)
void MGPCG::restrictToCoarserGrid(int fineLevel) {
    int coarseLevel = fineLevel + 1; 
    for (int i = grids_[coarseLevel]->pIBegin(); i < grids_[coarseLevel]->pIEnd(); ++i) {
        for (int j = grids_[coarseLevel]->pJBegin(); j < grids_[coarseLevel]->pJEnd(); ++j) {
            int i_fine = 2 * i - grids_[coarseLevel]->pIBegin();
            int j_fine = 2 * j - grids_[coarseLevel]->pJBegin();
            grids_[coarseLevel]->rhs(i, j) = 0.25 * (grids_[fineLevel]->rhs(i_fine, j_fine) +
                                                  grids_[fineLevel]->rhs(i_fine + 1, j_fine) +
                                                  grids_[fineLevel]->rhs(i_fine, j_fine + 1) +
                                                  grids_[fineLevel]->rhs(i_fine + 1, j_fine + 1));
        }
    }
}
//! TODO: use bilinear interpolation instead of constant interpolation
// Prolongate correction to finer grid (void version)
void MGPCG::prolongAndCorrectToFinerGrid(int coarseLevel) {
    int fineLevel = coarseLevel - 1;
    for (int i = grids_[coarseLevel]->pIBegin(); i < grids_[coarseLevel]->pIEnd(); ++i) {
        for (int j = grids_[coarseLevel]->pJBegin(); j < grids_[coarseLevel]->pJEnd(); ++j) {
            int i_fine = 2 * i - grids_[coarseLevel]->pIBegin();
            int j_fine = 2 * j - grids_[coarseLevel]->pJBegin();
            grids_[fineLevel]->p(i_fine, j_fine) += grids_[coarseLevel]->p(i, j);
            grids_[fineLevel]->p(i_fine + 1, j_fine) += grids_[coarseLevel]->p(i, j);
            grids_[fineLevel]->p(i_fine, j_fine + 1) += grids_[coarseLevel]->p(i, j);
            grids_[fineLevel]->p(i_fine + 1, j_fine + 1) += grids_[coarseLevel]->p(i, j);
        }
    }
}

