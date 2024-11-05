#include "computation.h"

/**
 * Initialize the computation object, parse the settings from the file that is given as the only command line argument
 * 
 */
void Computation::initialize(int argc, char* argv[]) {
    // Parse the settings from the file that is given as the only command line argument
}


/**
 * run the whole simulation until t_end
 */
void Computation::runSimulation() {

    // Loop over all time steps until t_end is reached
    
    computeTimeStepWidth();

    for (int time = 0; time * dt_ < settings_.endTime; ++time) {

        applyBoundaryValues();

        computePreliminaryVelocities();

        computeRightHandSide();

        computePressure();

        computeVelocities();

        outputWriterParaview_->writeFile(time * dt_); // Output
    }


}

void Computation::computeTimeStepWidth() {
    // Compute the time step width dt from maximum velocities
    // Idea: CFL-condition satisfied with factor 0.5
    double dx2 = discretization_->dx() * discretization_->dx();
    double dy2 = discretization_->dy() * discretization_->dy();

    // Compute CFL condition from diffusion operator
    double dt_diffusion = (settings_.re / 2.0) * (dx2 * dy2)/(dx2 + dy2);
    // Compute CFL condition from convection operator, use |u_max|  = 1, |v_max| = 1
    double dt_convection = 0.5 * std::min(discretization_->dx(), discretization_->dy());

    dt_ = std::min(dt_diffusion, dt_convection);
} 

    

void Computation::applyBoundaryValues() {
    // Set boundary values of u and v to correct values

    // loop over top boundary: 
    // set u,v to dirichletBcTop
    // set u at uJend() to dirichtletBcTop[0]: u at ujend() = 2 * dirichletBcTop[0] - u at ujend() -1
    // set v at vJend() to dirichtletBcTop[1]

    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        discretization_->u(i, discretization_->uJEnd()) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 1);
    }

    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        discretization_->v(i, discretization_->vJEnd()) = settings_.dirichletBcTop[1];
    }
    // loop over bottom boundary
    // set u,v to dirichletBcBottom
    // set u at ujbegin() - 1: u at ujbegin() - 1 = 2 * dirichletBcBottom[0] - u at ujbegin() 
    // set v at vjbegin() - 1 to dirichletBcBottom[1]
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        discretization_->u(i, discretization_->uJBegin() - 1) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin());
    }

    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        discretization_->v(i, discretization_->vJBegin() - 1) = settings_.dirichletBcBottom[1];
    }
    // loop over left boundary
    // set u,v to dirichletBcLeft
    // set u at uIbegin() - 1 to dirichletBcLeft[0]
    // set v at vIbegin() - 1 to 2 * dirichletBcLeft[1] - v at vIbegin()
    for (int j = (discretization_->uJBegin() - 1); j < (discretization_->uJEnd() + 1); j++) {
        discretization_->u(discretization_->uIBegin() - 1, j) = settings_.dirichletBcLeft[0];
    }

    for (int j = (discretization_->vJBegin() - 1); j < (discretization_->vJEnd() + 1); j++) {
        discretization_->v(discretization_->vIBegin() - 1, j) = 2 * settings_.dirichletBcLeft[1] -  discretization_->v(discretization_->vIBegin(), j);
    }

    // loop over right boundary
    // set u,v to dirichletBcRight
    // set u at uIend() to dirichletBcRight[0]
    // set v at vIend() to 2 * dirichletBcRight[1] - v at vIend() - 1
    for (int j = (discretization_->uJBegin() - 1); j < (discretization_->uJEnd() + 1); j++) {
        discretization_->u(discretization_->uIEnd(), j) = settings_.dirichletBcRight[0];
    }

    for (int j = (discretization_->vJBegin() - 1); j < (discretization_->vJEnd() + 1); j++) {
        discretization_->v(discretization_->vIEnd(), j) = 2 * settings_.dirichletBcRight[1] -  discretization_->v(discretization_->vIEnd() - 1, j);
    }

}

void Computation::computePreliminaryVelocities() {
    // Compute the preliminary velocities, F and G
}

void Computation::computeRightHandSide() {
    // Compute the right hand side of the Poisson equation for the pressure
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {

            double F_x = (discretization_->f(i, j) - discretization_->f(i - 1, j)) / discretization_->dx();
            double G_y = (discretization_->g(i, j) - discretization_->g(i, j - 1)) / discretization_->dy();

            discretization_->rhs(i, j) = (F_x + G_y) / dt_;
        }
    }
}

void Computation::computePressure() {
    // Solve the Poisson equation for the pressure
    pressureSolver_->solve();
}

void Computation::computeVelocities() {
    // Compute the new velocities, u,v from the preliminary velocities, F,G and the pressure, p
    
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            discretization_->u(i, j) = discretization_->f(i, j) - dt_ * discretization_->computeDpDx(i, j);
        }
    }

    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            discretization_->v(i, j) = discretization_->g(i, j) - dt_ * discretization_->computeDpDy(i, j);
        }
    }
}