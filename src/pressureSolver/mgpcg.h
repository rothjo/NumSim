#pragma once

#include "cg.h"
#include "gaussSeidel.h"
#include "sor.h"
#include "discretization/donorCell.h"
#include "partitioning/partitioning.h"
#include <memory>
#include <vector>

/**
 * Implementation of the Conjugate Gradient (CG) pressure solver.
 * Implements the solve method of the PressureSolver interface.
 */
class MGPCG : public CG {
public:
    MGPCG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning);

    /**
     * Solve poisson problem for the pressure, using the rhs and p field variables in staggeredGrid
     */
    void solve() override;

private:

    int levels_; ///< Number of levels in the multigrid hierarchy.

    std::vector<std::shared_ptr<Discretization>> grids_; ///< Hierarchy of grids (finest to coarsest).
    std::vector<std::shared_ptr<PressureSolver>> smoothers_; ///< Smoothers for each level.
    // std::vector<FieldVariable> errors_; ///< Error grid for each level.

    /**
     * @brief Coarsen the given grid to create the next level in the hierarchy.
     * @param fineGrid The finer grid to coarsen.
     * @return Coarsened grid.
     */
    std::shared_ptr<Discretization> coarsenGrid(int fineLevel);

    /**
     * @brief Compute the residual for a given grid level.
     * @param grid The grid to compute the residual on.
     * @param residual The residual field variable to populate.
     */
    void vCycle();	

    /**
     * @brief Compute the residual for a given grid level.
     * @param grid The grid to compute the residual on.
     * @param residual The residual field variable to populate.
     */
    void computeResidual(int level);

    /**
     * @brief Restrict a residual from a fine grid to a coarser grid.
     * @param fineResidual The residual on the finer grid.
     * @param coarseResidual The residual on the coarser grid to populate.
     */
    void restrictToCoarserGrid(int fineLevel);

    /**
     * @brief Prolongate a correction from a coarser grid to a finer grid.
     * @param coarseCorrection The correction on the coarser grid.
     * @param fineCorrection The correction on the finer grid to populate.
     */
    void prolongAndCorrectToFinerGrid(int coarseLevel);

    /**
     * @brief Apply a correction to the current grid level.
     * @param grid The grid to apply the correction to.
     * @param correction The correction field variable.
     */
    // void applyCorrection(int level);

    FieldVariable baseP_;
    FieldVariable baseRHS_;
    std::shared_ptr<Partitioning> partitioning_;
    FieldVariable z_;
};
