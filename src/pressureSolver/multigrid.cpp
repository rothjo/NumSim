#include "multigrid.h"
#include <cmath>
#include <iostream>

// Constructor
Multigrid::Multigrid(std::shared_ptr<Discretization> baseDiscretization, double epsilon,
                    int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning)
    : PressureSolver(baseDiscretization, epsilon, maximumNumberOfIterations),
      epsilon_(epsilon), maxIterations_(maximumNumberOfIterations),
      partitioning_(partitioning),
      baseRHS_(FieldVariable({baseDiscretization->nCells()[0] + 3, baseDiscretization->nCells()[1] + 3}, {-1.5 * baseDiscretization->meshWidth()[0], -1.5 * baseDiscretization->meshWidth()[1]}, baseDiscretization->meshWidth())),
      full_solver(std::make_shared<CG>(baseDiscretization, epsilon, maximumNumberOfIterations)),
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
    smoothers_.emplace_back(std::make_shared<GaussSeidel>(grids_[grids_.size() - 1], epsilon, 2));

    // for (int i = discretization_->pIBegin(); discretization_->pIEnd(); ++i) {
    //     for (int j = discretization_->pJBegin(); discretization_->pJEnd(); ++j) {

    //         baseRHS_(i, j) = discretization_->rhs(i, j);
    //     }
    // }
}

// Solve method
void Multigrid::solve() {
    // Perform a multigrid V-cycle
    // computeResidual(0);

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

    // for (int i = discretization_->pIBegin(); discretization_->pIEnd(); ++i) {
    //     for (int j = discretization_->pJBegin(); discretization_->pJEnd(); ++j) {
    //         std::cout << "i = " << i << ", j = " << j << std::endl;
    //         discretization_->rhs(i, j) = baseRHS_(i, j);
    //     }
    // }

    // full_solver->solve();

    // std::cout << "Multigrid reached maximum iterations without convergence.\n";
}

// Coarsen the grid
std::shared_ptr<Discretization> Multigrid::coarsenGrid(int fineLevel) {
    int nx = grids_[fineLevel]->nCells()[0] / 2;
    double dx = grids_[fineLevel]->dx() * 2;

    //! TODO: change 0.5 parameter to settings_.omega
    return std::make_shared<DonorCell>(std::array<int, 2>{nx, nx}, std::array<double, 2>{dx, dx}, partitioning_, 0.5);
}

// Compute residual (void version)
void Multigrid::computeResidual(int level) {
    for (int i = grids_[level]->pIBegin(); i < grids_[level]->pIEnd(); ++i) {
        for (int j = grids_[level]->pJBegin(); j < grids_[level]->pJEnd(); ++j) {
            const double dx_2 = grids_[level]->dx() * grids_[level]->dx();
            const double dy_2 = grids_[level]->dy() * grids_[level]->dy();
            double laplaceP = ((grids_[level]->p(i + 1, j) - 2.0 * grids_[level]->p(i, j) + grids_[level]->p(i - 1, j)) / dx_2) + ((grids_[level]->p(i, j + 1) - 2.0 * grids_[level]->p(i, j) + grids_[level]->p(i, j - 1)) / dy_2);
            grids_[level]->rhs(i, j) = grids_[level]->rhs(i, j) - laplaceP;
        }
    }
}

// Restrict residual to coarser grid (void version)
void Multigrid::restrictToCoarserGrid(int fineLevel) {
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
void Multigrid::prolongAndCorrectToFinerGrid(int coarseLevel) {
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


// Apply correction
// void Multigrid::applyCorrection(int level) {
//     for (int i = grids_[level]->pIBegin(); i < grids_[level]->pIEnd(); ++i) {
//         for (int j = grids_[level]->pJBegin(); j < grids_[level]->pJEnd(); ++j) {
//             grids_[level]->p(i, j) += errors_[level](i, j);
        
//     }
// }
// }
