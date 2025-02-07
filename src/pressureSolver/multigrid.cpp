#include "multigrid.h"
#include <cmath>
#include <iostream>

/**
 * @brief Constructor for the Multigrid solver.
 * @param baseDiscretization The finest-level discretization.
 * @param epsilon Convergence tolerance.
 * @param maximumNumberOfIterations Maximum iterations for the multigrid cycle.
 * @param cycle The type of cycle to use (V or W).
 * @param lowestLevel The lowest level of the multigrid hierarchy.
 * @param partitioning The partitioning object.
 * @param smoothingIterations Number of smoothing iterations.
 * @param coarseGridIterations Number of iterations on the coarsest grid.
 */
Multigrid::Multigrid(std::shared_ptr<Discretization> baseDiscretization, double epsilon,
                    int maximumNumberOfIterations, std::string cycle, int lowestLevel, std::shared_ptr<Partitioning> partitioning,
                    int smoothingIterations, int coarseGridIterations)
    : PressureSolver(baseDiscretization, epsilon, maximumNumberOfIterations),
        cycle_(cycle),
        lowestLevel_(lowestLevel),
        partitioning_(partitioning),
        smoothingIterations_(smoothingIterations),
        coarseGridIterations_(coarseGridIterations){
    // assert that the number of cells in each direction is a power of 2
    assert((baseDiscretization->nCells()[0] == baseDiscretization->nCells()[1] && baseDiscretization->nCells()[0] > 0 &&  (baseDiscretization->nCells()[0] & (baseDiscretization->nCells()[0] - 1)) == 0) &&
           "Number of cells in each direction must be a power of 2.");
    // assert that lowestLevel is smaller than the maximum level
    assert((lowestLevel_ <= std::log2(baseDiscretization->nCells()[0])) && "lowestLevel must be smaller than maxLevel!");
}

/**
 * @brief Solve the pressure equation using the multigrid method.
 */
void Multigrid::solve() {
    int maxCycles = 100;
    int iteration = 0;
    const double eps2 = epsilon_ * epsilon_;
    computeResidualNorm();

    if (cycle_ == "V") {
        while (residualNorm2_ > eps2 && iteration < maxCycles) {
        ++iteration;
        vCycle(discretization_);
        computeResidualNorm();
        }
    } else if (cycle_ == "W") {
        while (residualNorm2_ > eps2 && iteration < maxCycles) {
        ++iteration;
        wCycle(discretization_);
        computeResidualNorm();
        residualNormVector_.push_back(residualNorm2_);
        }
    }
    numberOfIterations_ = iteration;
}
/**
 * @brief Perform a V-cycle on a given discretization.
 * @param discretization The grid to perform the V-cycle on.
 */
void Multigrid::vCycle(std::shared_ptr<Discretization> discretization) {
    // if smallest grid size, do coarse level solve
    if (discretization->nCells()[0] == std::pow(2, lowestLevel_)) {
        GaussSeidel coarsesmoother = GaussSeidel(discretization, epsilon_, coarseGridIterations_);
        coarsesmoother.solve();
        solverIterations_.push_back(coarsesmoother.numberOfIterations());
        return;
    }
    // Smoothing
    GaussSeidel smoother = GaussSeidel(discretization, epsilon_, smoothingIterations_);
    smoother.solve();
    // Saving number of GaussSeidel Iterations used for smoothing
    solverIterations_.push_back(smoother.numberOfIterations());

    // Compute residual with p, rhs
    FieldVariable residual = discretization->p();
    computeResidual(discretization, residual);

    // Restrict residual to coarser grid -> New Fieldvariable rhs
    std::shared_ptr<Discretization> coarseDiscretization = std::make_shared<DonorCell>(std::array<int, 2>{discretization->nCells()[0] / 2, discretization->nCells()[1] / 2}, std::array<double, 2>{discretization->meshWidth()[0] * 2, discretization->meshWidth()[1] * 2}, partitioning_, 0.5);
    FieldVariable coarseResidual = coarseDiscretization->p();
    restrictToCoarserGrid(residual, coarseResidual, coarseDiscretization);

    for (int i = coarseDiscretization->pIBegin(); i < coarseDiscretization->pIEnd(); ++i) {
        for (int j = coarseDiscretization->pJBegin(); j < coarseDiscretization->pJEnd(); ++j) {
            coarseDiscretization->rhs(i, j) = coarseResidual(i, j);
        }
    }

    // enter recursion loop
    vCycle(coarseDiscretization);

    FieldVariable correction = discretization->p();

    // Prologate error to finer grid -> New Fieldvariable error p = p + prolongate(error)
    prolongation(coarseDiscretization, correction);

    // add error up
    for (int i = discretization->pIBegin(); i < discretization->pIEnd(); ++i) {
        for (int j = discretization->pJBegin(); j < discretization->pJEnd(); ++j) {
            discretization->p(i, j) += correction(i, j);
        }
    }

    // Post-smoothing p, rhs
    smoother.solve();
    solverIterations_.push_back(smoother.numberOfIterations());

}
/**
 * @brief Perform a W-cycle on a given discretization.
 * @param discretization The grid to perform the W-cycle on.
 */
void Multigrid::wCycle(std::shared_ptr<Discretization> discretization) {
    // if smallest grid size, do coarse level solve
    if (discretization->nCells()[0] == std::pow(2, lowestLevel_)) {
        GaussSeidel coarsesmoother = GaussSeidel(discretization, epsilon_, coarseGridIterations_);
        coarsesmoother.solve();
        solverIterations_.push_back(coarsesmoother.numberOfIterations());
        return;
    }
    // Pre-smoothing
    GaussSeidel smoother = GaussSeidel(discretization, epsilon_, smoothingIterations_);
    smoother.solve();
    solverIterations_.push_back(smoother.numberOfIterations());

    // Compute residual with p, rhs
    FieldVariable residual = discretization->p();
    computeResidual(discretization, residual);

    // Restrict residual to coarser grid -> New Fieldvariable rhs
    std::shared_ptr<Discretization> coarseDiscretization = std::make_shared<DonorCell>(std::array<int, 2>{discretization->nCells()[0] / 2, discretization->nCells()[1] / 2}, std::array<double, 2>{discretization->meshWidth()[0] * 2, discretization->meshWidth()[1] * 2}, partitioning_, 0.5);
    FieldVariable coarseResidual = coarseDiscretization->p();
    restrictToCoarserGrid(residual, coarseResidual, coarseDiscretization);

    for (int i = coarseDiscretization->pIBegin(); i < coarseDiscretization->pIEnd(); ++i) {
        for (int j = coarseDiscretization->pJBegin(); j < coarseDiscretization->pJEnd(); ++j) {
            coarseDiscretization->rhs(i, j) = coarseResidual(i, j);
        }
    }

    // enter recursion loop
    wCycle(coarseDiscretization);

    FieldVariable correction = discretization->p();

    // Prologate error to finer grid -> New Fieldvariable error p = p + prolongate(error)
    prolongation(coarseDiscretization, correction);

    // add error
    for (int i = discretization->pIBegin(); i < discretization->pIEnd(); ++i) {
        for (int j = discretization->pJBegin(); j < discretization->pJEnd(); ++j) {
            discretization->p(i, j) += correction(i, j);
        }
    }

    smoother.solve();

    // Compute residual with p, rhs
    residual = discretization->p();
    computeResidual(discretization, residual);

    // Restrict residual to coarser grid -> New Fieldvariable rhs
    coarseResidual = coarseDiscretization->p();
    restrictToCoarserGrid(residual, coarseResidual, coarseDiscretization);

    for (int i = coarseDiscretization->pIBegin(); i < coarseDiscretization->pIEnd(); ++i) {
        for (int j = coarseDiscretization->pJBegin(); j < coarseDiscretization->pJEnd(); ++j) {
            coarseDiscretization->rhs(i, j) = coarseResidual(i, j);
        }
    }
    // enter recursion loop
    wCycle(coarseDiscretization);

    correction = discretization->p();

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
    solverIterations_.push_back(smoother.numberOfIterations());

}

/**
 * @brief Compute the residual for a given grid level.
 * @param grid The grid to compute the residual on.
 * @param residual The residual field variable to populate.
 */
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

/**
 * @brief Restrict a residual from a fine grid to a coarser grid.
 * @param fineResidual The residual on the finer grid.
 * @param coarseResidual The residual on the coarser grid to populate.
 */
void Multigrid::restrictToCoarserGrid(FieldVariable& fineResidual, FieldVariable& coarseResidual, std::shared_ptr<Discretization> coarseDiscretization) {
    for (int i = coarseDiscretization->pIBegin(); i < coarseDiscretization->pIEnd(); ++i) {
        for (int j = coarseDiscretization->pJBegin(); j < coarseDiscretization->pJEnd(); ++j) {
            int i_fine = 2 * i - coarseDiscretization->pIBegin();
            int j_fine = 2 * j - coarseDiscretization->pJBegin();
            // std::cout << "i:" <<i << "j:"<< j << std::endl;
            // std::cout << coarseResidual.size()[0] << "loopsize" << coarseDiscretization->pJEnd() - coarseDiscretization->pJBegin() << std::endl;
            coarseResidual(i, j) = 0.25 * (fineResidual(i_fine, j_fine) +
                                                  fineResidual(i_fine + 1, j_fine) +
                                                  fineResidual(i_fine, j_fine + 1) +
                                                  fineResidual(i_fine + 1, j_fine + 1));
        }
    }
}
/**
 * @brief Prolongate a correction from a coarser grid to a finer grid.
 * @param coarseDiscretization The coarser grid.
 * @param correction The correction to prolongate.
 */
void Multigrid::prolongation(std::shared_ptr<Discretization> coarseDiscretization, FieldVariable& correction) {
    double x_origin = -1.5 * coarseDiscretization->meshWidth()[0];
    double y_origin = -1.5 * coarseDiscretization->meshWidth()[0];
    double half_dx = 0.5*coarseDiscretization->meshWidth()[0];
    for (int i = coarseDiscretization->pIBegin(); i < coarseDiscretization->pIEnd(); ++i) {
        for (int j = coarseDiscretization->pJBegin(); j < coarseDiscretization->pJEnd(); ++j) {
            int i_fine = 2 * i - coarseDiscretization->pIBegin();
            int j_fine = 2 * j - coarseDiscretization->pJBegin();
            double x = x_origin + i*coarseDiscretization->meshWidth()[0];
            double y = y_origin + j*coarseDiscretization->meshWidth()[1]; // Transformation to point in coarse grid
            correction(i_fine, j_fine) = coarseDiscretization->p().interpolateAt(x - half_dx , y - half_dx);
            correction(i_fine + 1, j_fine) = coarseDiscretization->p().interpolateAt(x + half_dx, y - half_dx);
            correction(i_fine, j_fine + 1) = coarseDiscretization->p().interpolateAt(x - half_dx, y + half_dx);
            correction(i_fine + 1, j_fine + 1) = coarseDiscretization->p().interpolateAt(x + half_dx, y + half_dx);
        }
    }
}

