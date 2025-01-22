
#include "pressureSolver.h"
#include "cg.h"
#include "discretization/discretization.h"
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
    Multigrid(std::shared_ptr<Discretization> baseDiscretization,
              int levels, double epsilon, int maximumNumberOfIterations);

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

    /**
     * @brief Coarsen the given grid to create the next level in the hierarchy.
     * @param fineGrid The finer grid to coarsen.
     * @return Coarsened grid.
     */
    std::shared_ptr<Discretization> coarsenGrid(std::shared_ptr<Discretization> fineGrid);

    /**
     * @brief Compute the residual for a given grid level.
     * @param grid The grid to compute the residual on.
     * @return Residual field variable.
     */
    FieldVariable computeResidual(std::shared_ptr<Discretization> grid);

    /**
     * @brief Restrict a residual from a fine grid to a coarser grid.
     * @param fineResidual The residual on the finer grid.
     * @return Restricted residual on the coarser grid.
     */
    FieldVariable restrictToCoarserGrid(const FieldVariable &fineResidual);

    /**
     * @brief Prolongate a correction from a coarser grid to a finer grid.
     * @param coarseCorrection The correction on the coarser grid.
     * @return Prolongated correction on the finer grid.
     */
    FieldVariable prolongToFinerGrid(const FieldVariable &coarseCorrection);

    /**
     * @brief Apply a correction to the current grid level.
     * @param grid The grid to apply the correction to.
     * @param correction The correction field variable.
     */
    void applyCorrection(std::shared_ptr<Discretization> grid, const FieldVariable &correction);

     /**
     * @brief Check if the solution has converged on the given grid level.
     * @param grid The grid to check convergence on.
     * @return True if converged, false otherwise.
     */
    bool checkConvergence(std::shared_ptr<Discretization> grid);
};
