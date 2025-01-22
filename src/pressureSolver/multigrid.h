#include "pressureSolver.h"
#include "cg.h"
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
    Multigrid(std::shared_ptr<Discretization> baseDiscretization, double epsilon, int maximumNumberOfIterations);

    /**
     * @brief Solve the pressure equation using the multigrid method.
     */
    void solve() override;

private:
    int levels_; ///< Number of levels in the multigrid hierarchy.
    double epsilon_; ///< Convergence tolerance.
    int maxIterations_; ///< Maximum number of iterations for V-cycles.

    std::vector<std::shared_ptr<Discretization>> grids_; ///< Hierarchy of grids (finest to coarsest).
    std::vector<std::shared_ptr<CG>> smoothers_; ///< Smoothers for each level.
    std::vector<FieldVariable> errors_; ///< Error grid for each level.

    /**
     * @brief Coarsen the given grid to create the next level in the hierarchy.
     * @param fineGrid The finer grid to coarsen.
     * @return Coarsened grid.
     */
    std::shared_ptr<Discretization> coarsenGrid(std::shared_ptr<Discretization> fineGrid);

    /**
     * @brief Compute the residual for a given grid level.
     * @param grid The grid to compute the residual on.
     * @param residual The residual field variable to populate.
     */
    void computeResidual(std::shared_ptr<Discretization> grid, FieldVariable &residual);

    /**
     * @brief Restrict a residual from a fine grid to a coarser grid.
     * @param fineResidual The residual on the finer grid.
     * @param coarseResidual The residual on the coarser grid to populate.
     */
    void restrictToCoarserGrid(const FieldVariable &fineResidual, FieldVariable &coarseResidual);

    /**
     * @brief Prolongate a correction from a coarser grid to a finer grid.
     * @param coarseCorrection The correction on the coarser grid.
     * @param fineCorrection The correction on the finer grid to populate.
     */
    void prolongToFinerGrid(const FieldVariable &coarseCorrection, FieldVariable &fineCorrection);

    /**
     * @brief Apply a correction to the current grid level.
     * @param grid The grid to apply the correction to.
     * @param correction The correction field variable.
     */
    void applyCorrection(std::shared_ptr<Discretization> grid, const FieldVariable &correction);

    std::shared_ptr<CG> full_solver;
    std::shared_ptr<Partitioning> partitioning_;
};
