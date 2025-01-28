#include "multigrid.h"
#include <cmath>
#include <iostream>

// Constructor
Multigrid::Multigrid(std::shared_ptr<Discretization> baseDiscretization, double epsilon,
                    int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning)
    : PressureSolver(baseDiscretization, epsilon, maximumNumberOfIterations),
      partitioning_(partitioning),
      levels_(std::log2(baseDiscretization->nCells()[0])) {
    assert((baseDiscretization->nCells()[0] == baseDiscretization->nCells()[1] && baseDiscretization->nCells()[0] > 0 &&  (baseDiscretization->nCells()[0] & (baseDiscretization->nCells()[0] - 1)) == 0) &&
           "Number of cells in each direction must be a power of 2.");
}

// Solve method
void Multigrid::solve() {

    for (int i = 0; i < 10; ++i) {
        vCycle(discretization_);
    }

}
//! TODO: implement with pointers
void Multigrid::vCycle(std::shared_ptr<Discretization> discretization) {

    if (discretization->nCells()[0] == 2) {
        GaussSeidel coarsesmoother = GaussSeidel(discretization, epsilon_, maximumNumberOfIterations_);
        coarsesmoother.solve();
        return;
    }
    // Pre-smoothing
    // Init gaussseidel, solve on discretization, eps, maximumNumberIter
    GaussSeidel smoother = GaussSeidel(discretization, epsilon_, 2);
    smoother.solve();

    // Compute residual with p, rhs
    FieldVariable residual = discretization->p();
    computeResidual(discretization, residual);

    // Restrict residual to coarser grid -> New Fieldvariable rhs
    std::shared_ptr<Discretization> coarseDiscretization = std::make_shared<DonorCell>(std::array<int, 2>{discretization->nCells()[0] / 2 + 3, discretization->nCells()[1] / 2 + 3}, std::array<double, 2>{discretization->meshWidth()[0] * 2, discretization->meshWidth()[1] * 2}, partitioning_, 0.5);
    FieldVariable coarseResidual = FieldVariable(coarseDiscretization->nCells(), {-1.5 * coarseDiscretization->meshWidth()[0], -1.5 * coarseDiscretization->meshWidth()[1]}, coarseDiscretization->meshWidth());
    restrictToCoarserGrid(residual, coarseResidual, coarseDiscretization);

    for (int i = coarseDiscretization->pIBegin(); i < coarseDiscretization->pIEnd(); ++i) {
        for (int j = coarseDiscretization->pJBegin(); j < coarseDiscretization->pJEnd(); ++j) {
            coarseDiscretization->rhs(i, j) = coarseResidual(i, j);
        }
    }

    // if smallest grid size, do coarse level solve
    // else vCycle(error, rhs)
    vCycle(coarseDiscretization);

    FieldVariable correction = discretization->p();

    // Prologate error to finer grid -> New Fieldvariable error p = p + prolongate(error)
    prolongation(coarseDiscretization, correction);

    // add error
    for (int i = discretization->pIBegin(); i < discretization->pIEnd(); ++i) {
        for (int j = discretization->pJBegin(); j < discretization->pJEnd(); ++j) {
            discretization->p(i, j) += correction(i, j);
        }
    }

    

    // Post-smoothing p, rhs
    smoother.solve();

}

// Compute residual (void version)
void Multigrid::computeResidual(std::shared_ptr<Discretization> discretization, FieldVariable& residual) {
    for (int i = discretization->pIBegin(); i < discretization->pIEnd(); ++i) {
        for (int j = discretization->pJBegin(); j < discretization->pJEnd(); ++j) {
            const double dx_2 = discretization->dx() * discretization->dx();
            const double dy_2 = discretization->dy() * discretization->dy();
            double laplaceP = ((discretization->p(i + 1, j) - 2.0 * discretization->p(i, j) + discretization->p(i - 1, j)) / dx_2) + ((discretization->p(i, j + 1) - 2.0 * discretization->p(i, j) + discretization->p(i, j - 1)) / dy_2);
            residual(i, j) = discretization->rhs(i, j) - laplaceP;
        }
    }
}

void Multigrid::copyFieldVariable(FieldVariable& source, FieldVariable& target, std::shared_ptr<Discretization> discretization) {
    for (int i = discretization->pIBegin(); i < discretization->pIEnd(); ++i) {
        for (int j = discretization->pJBegin(); j < discretization->pJEnd(); ++j) {
            target(i, j) = source(i, j);
        }
    }
}

// Restrict residual to coarser grid (void version)
void Multigrid::restrictToCoarserGrid(FieldVariable& fineResidual, FieldVariable& coarseResidual, std::shared_ptr<Discretization> coarseDiscretization) {
    for (int i = coarseDiscretization->pIBegin(); i < coarseDiscretization->pIEnd(); ++i) {
        for (int j = coarseDiscretization->pJBegin(); j < coarseDiscretization->pJEnd(); ++j) {
            int i_fine = 2 * i - coarseDiscretization->pIBegin();
            int j_fine = 2 * j - coarseDiscretization->pJBegin();
            coarseResidual(i, j) = 0.25 * (fineResidual(i_fine, j_fine) +
                                                  fineResidual(i_fine + 1, j_fine) +
                                                  fineResidual(i_fine, j_fine + 1) +
                                                  fineResidual(i_fine + 1, j_fine + 1));
        }
    }
}
//! TODO: use bilinear interpolation instead of constant interpolation
// Prolongate correction to finer grid (void version)
void Multigrid::prolongation(std::shared_ptr<Discretization> coarseDiscretization, FieldVariable& correction) {
    for (int i = coarseDiscretization->pIBegin(); i < coarseDiscretization->pIEnd(); ++i) {
        for (int j = coarseDiscretization->pJBegin(); j < coarseDiscretization->pJEnd(); ++j) {
            int i_fine = 2 * i - coarseDiscretization->pIBegin();
            int j_fine = 2 * j - coarseDiscretization->pJBegin();
            correction(i_fine, j_fine) = coarseDiscretization->p(i, j);
            correction(i_fine + 1, j_fine) = coarseDiscretization->p(i, j);
            correction(i_fine, j_fine + 1) = coarseDiscretization->p(i, j);
            correction(i_fine + 1, j_fine + 1) = coarseDiscretization->p(i, j);
        }
    }
}
