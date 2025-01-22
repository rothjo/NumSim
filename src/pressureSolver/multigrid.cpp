#include "multigrid.h"
#include <cmath>
#include <iostream>

// Constructor
Multigrid::Multigrid(std::shared_ptr<Discretization> baseDiscretization,
                    double epsilon, int maximumNumberOfIterations)
    : PressureSolver(baseDiscretization, epsilon, maximumNumberOfIterations),
      epsilon_(epsilon), maxIterations_(maximumNumberOfIterations),
      full_solver(std::make_shared<CG>(baseDiscretization, epsilon, maximumNumberOfIterations)),
      levels_(std::log2(baseDiscretization->nCells()[0])) {
    assert((baseDiscretization->nCells()[0] == baseDiscretization->nCells()[1] && baseDiscretization->nCells()[0] > 0 &&  (baseDiscretization->nCells()[0] & (baseDiscretization->nCells()[0] - 1)) == 0) &&
           "Number of cells in each direction must be a power of 2.");
    grids_.push_back(baseDiscretization);

    // Build hierarchy of grids
    for (int l = 1; l < levels_; ++l) {
        grids_.push_back(coarsenGrid(grids_[l - 1]));
    }

    // Create a smoother for each grid level
    for (const auto& grid : grids_) {
        smoothers_.emplace_back(std::make_shared<CG>(grid, epsilon, 2));
        errors_.emplace_back(FieldVariable({grid->nCells()[0] + 3, grid->nCells()[1] + 3}, {-1.5 * grid->dx(), -1.5 * grid->dy()}, grid->meshWidth()));
    }
}

// Solve method
void Multigrid::solve() {
    // Perform a multigrid V-cycle
    computeResidual(grids_[0], errors_[0]);

    // Downward sweep: Restrict residuals to coarser grids
    for (int l = 0; l < levels_ - 1; ++l) {
        smoothers_[l]->solve();  // Pre-smoothing
        computeResidual(grids_[l], errors_[l]);
        restrictToCoarserGrid(errors_[l], errors_[l + 1]);
    }

    // Solve on the coarsest grid
    smoothers_[levels_ - 1]->solve();

    // Upward sweep: Prolong errors_ to finer grids
    for (int l = levels_ - 2; l >= 0; --l) {
        prolongToFinerGrid(errors_[l + 1], errors_[l]);
        applyCorrection(grids_[l], errors_[l]);
        smoothers_[l]->solve();  // Post-smoothing
    }
    full_solver->solve();

    std::cout << "Multigrid reached maximum iterations without convergence.\n";
}

// Coarsen the grid
std::shared_ptr<Discretization> Multigrid::coarsenGrid(std::shared_ptr<Discretization> fineGrid) {
    int nx = fineGrid->nCells()[0] / 2;
    double dx = fineGrid->dx() * 2;

    return std::make_shared<DonorCell>(std::array<int, 2>{nx, nx}, std::array<double, 2>{dx, dx}, partitioning_, 0.5);
}

// Compute residual (void version)
void Multigrid::computeResidual(std::shared_ptr<Discretization> grid, FieldVariable &residual) {
    for (int i = grid->pIBegin(); i < grid->pIEnd(); ++i) {
        for (int j = grid->pJBegin(); j < grid->pJEnd(); ++j) {
            residual(i, j) = grid->rhs(i, j) - laplaceP(i, j);
        }
    }
}

// Restrict residual to coarser grid (void version)
void Multigrid::restrictToCoarserGrid(const FieldVariable &fineResidual, FieldVariable &coarseResidual) {
    for (int i = 1; i < fineResidual.size()[0] - 1; i += 2) {
        for (int j = 1; j < fineResidual.size()[1] - 1; j += 2) {
            coarseResidual(i / 2, j / 2) = 0.25 * (fineResidual(i, j) +
                                                  fineResidual(i + 1, j) +
                                                  fineResidual(i, j + 1) +
                                                  fineResidual(i + 1, j + 1));
        }
    }
}

// Prolongate correction to finer grid (void version)
void Multigrid::prolongToFinerGrid(const FieldVariable &coarseCorrection, FieldVariable &fineCorrection) {
    for (int i = 0; i < coarseCorrection.size()[0]; ++i) {
        for (int j = 0; j < coarseCorrection.size()[1]; ++j) {
            fineCorrection(2 * i, 2 * j) = coarseCorrection(i, j);
            fineCorrection(2 * i + 1, 2 * j) = coarseCorrection(i, j);
            fineCorrection(2 * i, 2 * j + 1) = coarseCorrection(i, j);
            fineCorrection(2 * i + 1, 2 * j + 1) = coarseCorrection(i, j);
        }
    }
}

// Apply correction
void Multigrid::applyCorrection(std::shared_ptr<Discretization> grid, const FieldVariable &correction) {
    for (int i = grid->pIBegin(); i < grid->pIEnd(); ++i) {
        for (int j = grid->pJBegin(); j < grid->pJEnd(); ++j) {
            grid->p(i, j) += correction(i, j);
        }
    }
}
