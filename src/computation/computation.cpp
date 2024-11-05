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

    // Loop over all time steps untill t_end is reached


    // Set boundary values of u and v to correct values
    applyBoundaryValues();


    // Compute the time step width dt from maximum velocities
    computeTimeStepWidth();



    // Compute the preliminary velocities, F and G 
    // SET BOUNDARIES FOR F,G HERE???
    computePreliminaryVelocities();

    // Compute the right hand side of the Poisson equation for the pressure
    computeRightHandSide();

    // Solve the Poisson equation for the pressure
    computePressure();

    // Compute the new velocities, u,v from the preliminary velocities, F,G and the pressure, p
    computeVelocities();


    //output_uvp
}

void Computation::computeTimeStepWidth() {
    // Compute the time step width dt from maximum velocities
}

void Computation::applyBoundaryValues() {
    // Set boundary values of u and v to correct values
}

void Computation::computePreliminaryVelocities() {
    // Compute the preliminary velocities, F and G
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            // Compute F
            double f_diffusion_term = (discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j))/settings_.re;
            double f_convection_term = (discretization_->computeDu2Dx(i,j) + discretization_->computeDuvDy(i,j));
            double discretization_->f(i,j) = discretization_->u(i,j) + dt*(f_diffusion_term - f_convection_term);
        }
    }
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            // Compute G
            double g_diffusion_term = (discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j))/settings_.re;
            double g_convection_term = (discretization_->computeDuvDx(i,j) + discretization_->computeDv2Dy(i,j));
            double discretization_->g(i,j) = discretization_->v(i,j) + dt*(g_diffusion_term - g_convection_term);
        }
    }     
}

void Computation::computeRightHandSide() {
    // Compute the right hand side of the Poisson equation for the pressure
}

void Computation::computePressure() {
    // Solve the Poisson equation for the pressure
}

void Computation::computeVelocities() {
    // Compute the new velocities, u,v from the preliminary velocities, F,G and the pressure, p
}