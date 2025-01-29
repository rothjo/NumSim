#include "pressureSolver.h"
#include "cg.h"
#include "gaussSeidel.h"
#include "sor.h"
#include "discretization/donorCell.h"
#include "partitioning/partitioning.h"
#include <memory>
#include <vector>

/**
 * @brief Multigrid solver class, derived from ParallelPressureSolver.
 * Implements the geometric multigrid method with V-cycles.
 */
class Multigrid : public PressureSolver {
public:
    /**
     * @brief Constructor for the Multigrid solver.
     * @param baseDiscretization The finest-level discretization.
     * @param levels Number of levels in the multigrid hierarchy.
     * @param epsilon Convergence tolerance.
     * @param maximumNumberOfIterations Maximum iterations for the multigrid cycle.
     */
    Multigrid(std::shared_ptr<Discretization> baseDiscretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning);

    /**
     * @brief Solve the pressure equation using the multigrid method.
     */
    void solve() override;

private:
    int levels_; ///< Number of levels in the multigrid hierarchy.

    std::vector<std::shared_ptr<Discretization>> grids_; ///< Hierarchy of grids (finest to coarsest).
    std::vector<std::shared_ptr<PressureSolver>> smoothers_; ///< Smoothers for each level.
    // std::vector<FieldVariable> errors_; ///< Error grid for each level.


    void vCycle(std::shared_ptr<Discretization> discretization);

    void wCycle(std::shared_ptr<Discretization> discretization);

    /**
     * @brief Compute the residual for a given grid level.
     * @param grid The grid to compute the residual on.
     * @param residual The residual field variable to populate.
     */
    void computeResidual(std::shared_ptr<Discretization> discretization, FieldVariable& residual);

    /**
     * @brief Copy the values of a source field variable to a target field variable.
     * @param source The source field variable.
     * @param target The target field variable.
     */
    void copyFieldVariable(FieldVariable& source, FieldVariable& target, std::shared_ptr<Discretization> discretization);

    /**
     * @brief Restrict a residual from a fine grid to a coarser grid.
     * @param fineResidual The residual on the finer grid.
     * @param coarseResidual The residual on the coarser grid to populate.
     */
    void restrictToCoarserGrid(FieldVariable& fineResidual, FieldVariable& coarseResidual, std::shared_ptr<Discretization> coarseDiscretization);

    void prolongation(std::shared_ptr<Discretization> coarseDiscretization, FieldVariable& correction);  

    /**
     * @brief Apply a correction to the current grid level.
     * @param grid The grid to apply the correction to.
     * @param correction The correction field variable.
     */
    // void applyCorrection(int level);

    std::shared_ptr<Partitioning> partitioning_;
};
