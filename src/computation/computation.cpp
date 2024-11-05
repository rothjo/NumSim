#include "computation.h"

/**
 * Initialize the computation object, parse the settings from the file that is given as the only command line argument
 * 
 */
void Computation::initialize(int argc, char* argv[]) {
    // load settigns from file
    settings_ = Settings();
    settings_.loadFromFile("parameters.txt");

    // init meshWidth
    meshWidth_ = {settings_.physicalSize[0] / settings_.nCells[0], settings_.physicalSize[1] / settings_.nCells[1]};

    // init discretization
    if (settings_.useDonorCell) {
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    }

    // init output writers
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
    outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);

    // init pressure solvers
    if (settings_.pressureSolver == "SOR") {
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    } else if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    } else {
        std::cerr << "Unknown pressure solver: " << settings_.pressureSolver << std::endl;
        std::exit(1);
    }
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