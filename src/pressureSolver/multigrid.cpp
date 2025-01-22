#include "multigrid.h"
#include <cmath>
#include <iostream>

// Constructor
Multigrid::Multigrid(std::shared_ptr<Discretization> baseDiscretization,
                     int levels, double epsilon, int maximumNumberOfIterations)
    : PressureSolver(baseDiscretization, epsilon, maximumNumberOfIterations),
      levels_(levels), epsilon_(epsilon), maxIterations_(maximumNumberOfIterations) {
    grids_.push_back(baseDiscretization);

    // Build hierarchy of grids
    for (int l = 1; l < levels; ++l) {
        grids_.push_back(coarsenGrid(grids_[l - 1]));
    }

    // Create a smoother for each grid level
    for (const auto &grid : grids_) {
        smoothers_.emplace_back(std::make_shared<CG>(grid, epsilon, maximumNumberOfIterations));
    }
}

// Solve method
void Multigrid::solve() {
    for (int k = 0; k < maxIterations_; ++k) {
        // Perform a multigrid V-cycle
        FieldVariable residual = computeResidual(grids_[0]);
        std::vector<FieldVariable> corrections(levels_);

        // Downward sweep: Restrict residuals to coarser grids
        for (int l = 0; l < levels_ - 1; ++l) {
            smoothers_[l]->solve();  // Pre-smoothing
            corrections[l] = computeResidual(grids_[l]);
            corrections[l + 1] = restrictToCoarserGrid(corrections[l]);
        }

        // Solve on the coarsest grid
        smoothers_[levels_ - 1]->solve();

        // Upward sweep: Prolong corrections to finer grids
        for (int l = levels_ - 2; l >= 0; --l) {
            FieldVariable correctionFiner = prolongToFinerGrid(corrections[l + 1]);
            applyCorrection(grids_[l], correctionFiner);
            smoothers_[l]->solve();  // Post-smoothing
        }

        // Check convergence
        if (checkConvergence(grids_[0])) {
            std::cout << "Multigrid converged after " << k + 1 << " iterations.\n";
            return;
        }
    }

    std::cout << "Multigrid reached maximum iterations without convergence.\n";
}

// Coarsen the grid
std::shared_ptr<Discretization> Multigrid::coarsenGrid(std::shared_ptr<Discretization> fineGrid) {
    // Create a coarser grid by reducing the number of cells and adjusting the mesh width
    int nx = fineGrid->nx() / 2;
    int ny = fineGrid->ny() / 2;
    double dx = fineGrid->dx() * 2;
    double dy = fineGrid->dy() * 2;

    return std::make_shared<Discretization>(nx, ny, dx, dy);
}

// Compute residual
FieldVariable Multigrid::computeResidual(std::shared_ptr<Discretization> grid) {
    FieldVariable residual(grid->size(), grid->origin(), grid->meshWidth());

    for (int i = grid->pIBegin(); i < grid->pIEnd(); ++i) {
        for (int j = grid->pJBegin(); j < grid->pJEnd(); ++j) {
            residual(i, j) = grid->rhs(i, j) - (grid->LaplaceP(i, j));
        }
    }

    return residual;
}

// Restrict residual to coarser grid
FieldVariable Multigrid::restrictToCoarserGrid(const FieldVariable &fineResidual) {
    FieldVariable coarseResidual(fineResidual.size() / 2, fineResidual.origin(), fineResidual.meshWidth() * 2);

    for (int i = 1; i < fineResidual.size()[0] - 1; i += 2) {
        for (int j = 1; j < fineResidual.size()[1] - 1; j += 2) {
            coarseResidual(i / 2, j / 2) = 0.25 * (fineResidual(i, j) +
                                                  fineResidual(i + 1, j) +
                                                  fineResidual(i, j + 1) +
                                                  fineResidual(i + 1, j + 1));
        }
    }

    return coarseResidual;
}

// Prolongate correction to finer grid
FieldVariable Multigrid::prolongToFinerGrid(const FieldVariable &coarseCorrection) {
    FieldVariable fineCorrection(coarseCorrection.size() * 2, coarseCorrection.origin(), coarseCorrection.meshWidth() / 2);

    for (int i = 0; i < coarseCorrection.size()[0]; ++i) {
        for (int j = 0; j < coarseCorrection.size()[1]; ++j) {
            fineCorrection(2 * i, 2 * j) = coarseCorrection(i, j);
            fineCorrection(2 * i + 1, 2 * j) = coarseCorrection(i, j);
            fineCorrection(2 * i, 2 * j + 1) = coarseCorrection(i, j);
            fineCorrection(2 * i + 1, 2 * j + 1) = coarseCorrection(i, j);
        }
    }

    return fineCorrection;
}

// Apply correction
void Multigrid::applyCorrection(std::shared_ptr<Discretization> grid, const FieldVariable &correction) {
    for (int i = grid->pIBegin(); i < grid->pIEnd(); ++i) {
        for (int j = grid->pJBegin(); j < grid->pJEnd(); ++j) {
            grid->p(i, j) += correction(i, j);
        }
    }
}

// Check convergence
bool Multigrid::checkConvergence(std::shared_ptr<Discretization> grid) {
    double residualNorm = 0.0;

    for (int i = grid->pIBegin(); i < grid->pIEnd(); ++i) {
        for (int j = grid->pJBegin(); j < grid->pJEnd(); ++j) {
            double res = grid->rhs(i, j) - grid->LaplaceP(i, j);
            residualNorm += res * res;
        }
    }

    double globalResidualNorm;
    MPI_Allreduce(&residualNorm, &globalResidualNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return std::sqrt(globalResidualNorm) < epsilon_;
}
