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