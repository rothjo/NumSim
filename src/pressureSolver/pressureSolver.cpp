#include "pressureSolver/pressureSolver.h"


/**
 * Implementation of the pressure solver. 
 */
PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations) :
    
    discretization_(discretization),
    epsilon_(epsilon),
    maximumNumberOfIterations_(maximumNumberOfIterations) {


}

void PressureSolver::setBoundaryValues() {

}

void PressureSolver::computeResidualNorm() {


}

