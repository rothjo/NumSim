#include "pressureSolver.h"
#include "cg.h"
#include "gaussSeidel.h"
#include "discretization/donorCell.h"
#include "partitioning/partitioning.h"
#include <memory>
#include <vector>
#include <settings.h>

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
    Multigrid(std::shared_ptr<Discretization> baseDiscretization, double epsilon, int maximumNumberOfIterations, std::string cycle, int lowestLevel, std::shared_ptr<Partitioning> partitioning, int smoothingIterations, int coarseGridIterations);

    /**
     * @brief Solve the pressure equation using the multigrid method.
     */
    void solve() override;


private:
    /**
     * @brief Perform a V-cycle on a given discretization.
     * @param discretization The grid to perform the V-cycle on.
     */
    void vCycle(std::shared_ptr<Discretization> discretization);

    /**
     * @brief Perform a W-cycle on a given discretization.
     * @param discretization The grid to perform the W-cycle on.
     */
    void wCycle(std::shared_ptr<Discretization> discretization);

    /**
     * @brief Compute the residual for a given grid level.
     * @param grid The grid to compute the residual on.
     * @param residual The residual field variable to populate.
     */
    void computeResidual(std::shared_ptr<Discretization> discretization, FieldVariable& residual);

    /**
     * @brief Restrict a residual from a fine grid to a coarser grid.
     * @param fineResidual The residual on the finer grid.
     * @param coarseResidual The residual on the coarser grid to populate.
     */
    void restrictToCoarserGrid(FieldVariable& fineResidual, FieldVariable& coarseResidual, std::shared_ptr<Discretization> coarseDiscretization);

    /**
     * @brief Prolongate a correction from a coarser grid to a finer grid.
     * @param coarseDiscretization The coarser grid.
     * @param correction The correction to prolongate.
     */
    void prolongation(std::shared_ptr<Discretization> coarseDiscretization, FieldVariable& correction);  

    std::shared_ptr<Partitioning> partitioning_; ///< The partitioning object.
    std::string cycle_; ///< The type of cycle to use (V or W).
    int lowestLevel_; ///< The lowest level in the multigrid hierarchy.
    int smoothingIterations_; ///< Number of smoothing iterations.
    int coarseGridIterations_; ///< Number of iterations on the coarsest grid.

};
