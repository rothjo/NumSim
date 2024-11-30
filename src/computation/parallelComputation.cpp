#include "parallelComputation.h"
#include <cmath>
#include <unistd.h>


void ParallelComputation::initialize(int argc, char* argv[]) {
    // load settigns from file
    std::string filename = argv[1];
    settings_.loadFromFile(filename);

    meshWidth_ = {settings_.physicalSize[0] / settings_.nCells[0], settings_.physicalSize[1] / settings_.nCells[1]};

    partitioning_ = std::make_shared<Partitioning>();
    partitioning_->initialize(settings_.nCells);
    std::array<int,2> localCells = partitioning_->nCellsLocal();

    // init discretization
    if (settings_.useDonorCell) {
        discretization_ = std::make_shared<DonorCell>(localCells, meshWidth_, partitioning_, settings_.alpha);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(localCells, meshWidth_, partitioning_);
    }

    // init output writers
    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, *partitioning_);
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discretization_, *partitioning_);

    // init pressure solvers
    if (settings_.pressureSolver == "SOR") {
        pressureSolver_ = std::make_unique<ParallelSOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega, partitioning_);
    } else if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<ParallelGaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_);
    } else {
        std::cerr << "Unknown pressure solver: " << settings_.pressureSolver << std::endl;
        std::exit(1);
    }
}

void ParallelComputation::runSimulation() {
    applyInitalBoundaryValues();

    double time = 0.0;
    double time_epsilon = 1e-8;
    int output = 1;

    
    // Loop over all time steps until t_end is reached
    while (time < (settings_.endTime - time_epsilon)) {

        applyBoundaryValues();


        computeTimeStepWidth();
  

        // Check to choose the last time step differently to hit t_end exactly
        if (time + dt_ > settings_.endTime - time_epsilon) {
            dt_ = settings_.endTime - time;
        }
        time += dt_;

        // Output
        if (time >= output) {
            outputWriterParaview_->writeFile(time); // Output
            outputWriterText_->writeFile(time); // Output
            output++;
        }

    	
        computePreliminaryVelocities();
    

        computeRightHandSide();
        
      

        computePressure();
        

        computeVelocities();
        
    }
}

//! TODO: Maybe change that u_max and v_max is communicated instead of dt_min
void ParallelComputation::computeTimeStepWidth() {
    const double dx = discretization_->dx();
    const double dy = discretization_->dy();

    const double dx2 = dx * dx;
    const double dy2 = dy * dy;

    const double dt_diffusion = (settings_.re / 2.0) * (dx2 * dy2)/(dx2 + dy2);

    double u_max = partitioning_->globalMax(discretization_->u().computeMaxAbs());
    double v_max = partitioning_->globalMax(discretization_->v().computeMaxAbs());
    
    // compute local timesteps
    const double dt_convection_x = dx / u_max;
    const double dt_convection_y = dy / v_max;
    const double dt_convection = std::min(dt_convection_x, dt_convection_y);


    const double dt = settings_.tau * std::min(dt_diffusion, dt_convection);

    if (dt > settings_.maximumDt) {
        std::cout << "Warning: Time step width is larger than maximum time step width. Using maximum time step width instead." << std::endl;
        dt_ = settings_.maximumDt;
    } else {
        dt_ = dt;
    }
}

void ParallelComputation::applyInitalBoundaryValues() {
    // u and f boundary values

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = (discretization_->uJBegin() - 1); j < (discretization_->uJEnd() + 1); j++) {
            // left boundary
            discretization_->u(discretization_->uIBegin() - 1, j) = settings_.dirichletBcLeft[0];
            discretization_->f(discretization_->uIBegin() - 1, j) = settings_.dirichletBcLeft[0];
        }
    } 

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = (discretization_->uJBegin() - 1); j < (discretization_->uJEnd() + 1); j++) {
            // right boundary
            discretization_->u(discretization_->uIEnd(), j) = settings_.dirichletBcRight[0];
            discretization_->f(discretization_->uIEnd(), j) = settings_.dirichletBcRight[0];
        }
    }

    // v and g boundary values

    if (partitioning_->ownPartitionContainsTopBoundary()){
        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
            // top boundary
            discretization_->v(i, discretization_->vJEnd()) = settings_.dirichletBcTop[1];
            discretization_->g(i, discretization_->vJEnd()) = settings_.dirichletBcTop[1];
        }
    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        // bottom boundary
        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
            discretization_->v(i, discretization_->vJBegin() - 1) = settings_.dirichletBcBottom[1];
            discretization_->g(i, discretization_->vJBegin() - 1) = settings_.dirichletBcBottom[1];
        }
    }
}
//! TODO: Concate u and v for communication
void ParallelComputation::applyBoundaryValues() {

    // define u buffers
    std::vector<double> uTopBuffer(discretization_->uIEnd() - discretization_->uIBegin(), 0.0);
    std::vector<double> uBottomBuffer(discretization_->uIEnd() - discretization_->uIBegin(), 0.0);
    std::vector<double> uRightBuffer(discretization_->uJEnd() - discretization_->uJBegin(), 0.0);
    std::vector<double> uLeftBuffer(discretization_->uJEnd() - discretization_->uJBegin(), 0.0);

    // define v buffers
    std::vector<double> vTopBuffer(discretization_->vIEnd() - discretization_->vIBegin(), 0.0);
    std::vector<double> vBottomBuffer(discretization_->vIEnd() - discretization_->vIBegin(), 0.0);
    std::vector<double> vRightBuffer(discretization_->vJEnd() - discretization_->vJBegin(), 0.0);
    std::vector<double> vLeftBuffer(discretization_->vJEnd() - discretization_->vJBegin(), 0.0);


    // define requests
    MPI_Request requestuTop, requestuBottom, requestuLeft, requestuRight;
    MPI_Request requestvTop, requestvBottom, requestvLeft, requestvRight;

    
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
            discretization_->u(i, discretization_->uJEnd()) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 1);
        }
    } else {
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
            uTopBuffer[i - discretization_->uIBegin()] = discretization_->u(i, discretization_->uJEnd() - 1);
        }
        partitioning_->send(uTopBuffer, partitioning_->topNeighbourRankNo(), requestuTop);

        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
            vTopBuffer[i - discretization_->vIBegin()] = discretization_->v(i, discretization_->vJEnd() - 2);
        }
        partitioning_->send(vTopBuffer, partitioning_->topNeighbourRankNo(), requestvTop);

        partitioning_->receive(uTopBuffer, partitioning_->topNeighbourRankNo(), requestuTop);
        partitioning_->receive(vTopBuffer, partitioning_->topNeighbourRankNo(), requestvTop);
    }

    // Communication to bottom neighbour
    if (partitioning_->ownPartitionContainsBottomBoundary()){
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
            discretization_->u(i, discretization_->uJBegin() - 1) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin());
        }
    } else {
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
            uBottomBuffer[i - discretization_->uIBegin()] = discretization_->u(i, discretization_->uJBegin());
        }
        partitioning_->send(uBottomBuffer, partitioning_->bottomNeighbourRankNo(), requestuBottom);

        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
            vBottomBuffer[i - discretization_->vIBegin()] = discretization_->v(i, discretization_->vJBegin() + 1);
        }
        partitioning_->send(vBottomBuffer, partitioning_->bottomNeighbourRankNo(), requestvBottom);

        partitioning_->receive(uBottomBuffer, partitioning_->bottomNeighbourRankNo(), requestuBottom);
        partitioning_->receive(vBottomBuffer, partitioning_->bottomNeighbourRankNo(), requestvBottom);
    }

    // Communication to left neighbour
    if (partitioning_->ownPartitionContainsLeftBoundary()){
        for (int j = (discretization_->vJBegin() - 1); j < (discretization_->vJEnd() + 1); j++) {
            discretization_->v(discretization_->vIBegin() - 1, j) = 2 * settings_.dirichletBcLeft[1] -  discretization_->v(discretization_->vIBegin(), j);
        }   
    } else {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            uLeftBuffer[j - discretization_->uJBegin()] = discretization_->u(discretization_->uIBegin() + 1, j);
        }
        partitioning_->send(uLeftBuffer, partitioning_->leftNeighbourRankNo(), requestuLeft);

        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            vLeftBuffer[j - discretization_->vJBegin()] = discretization_->v(discretization_->vIBegin(), j);
        }
        partitioning_->send(vLeftBuffer, partitioning_->leftNeighbourRankNo(), requestvLeft);

        partitioning_->receive(uLeftBuffer, partitioning_->leftNeighbourRankNo(), requestuLeft);
        partitioning_->receive(vLeftBuffer, partitioning_->leftNeighbourRankNo(), requestvLeft);
    }

    // Communication to right neighbour
    if (partitioning_->ownPartitionContainsRightBoundary()){
        for (int j = (discretization_->vJBegin() - 1); j < (discretization_->vJEnd() + 1); j++) {
            discretization_->v(discretization_->vIEnd(), j) = 2 * settings_.dirichletBcRight[1] -  discretization_->v(discretization_->vIEnd() - 1, j);
        }
    } else {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            uRightBuffer[j - discretization_->uJBegin()] = discretization_->u(discretization_->uIEnd() - 2, j);
        }
        partitioning_->send(uRightBuffer, partitioning_->rightNeighbourRankNo(), requestuRight);

        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            vRightBuffer[j - discretization_->vJBegin()] = discretization_->v(discretization_->vIEnd() - 1, j);
        }
        partitioning_->send(vRightBuffer, partitioning_->rightNeighbourRankNo(), requestvRight);

        partitioning_->receive(uRightBuffer, partitioning_->rightNeighbourRankNo(), requestuRight);
        partitioning_->receive(vRightBuffer, partitioning_->rightNeighbourRankNo(), requestvRight);
    }

    // Set the received values
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestvTop, MPI_STATUS_IGNORE);
        MPI_Wait(&requestuTop, MPI_STATUS_IGNORE);
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
            discretization_->u(i, discretization_->uJEnd()) = uTopBuffer[i - discretization_->uIBegin()];
        }
        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
            discretization_->v(i, discretization_->vJEnd()) = vTopBuffer[i - discretization_->vIBegin()];
        }
    }

    if(!partitioning_->ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&requestvBottom, MPI_STATUS_IGNORE);
        MPI_Wait(&requestuBottom, MPI_STATUS_IGNORE);
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
            discretization_->u(i, discretization_->uJBegin() -1) = uBottomBuffer[i - discretization_->uIBegin()];
        }
        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
            discretization_->v(i, discretization_->vJBegin() -1) = vBottomBuffer[i - discretization_->vIBegin()];
        }
    }

    if(!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestvLeft, MPI_STATUS_IGNORE);
        MPI_Wait(&requestuLeft, MPI_STATUS_IGNORE);
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            discretization_->u(discretization_->uIBegin()-1, j) = uLeftBuffer[j - discretization_->uJBegin()];
        }
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            discretization_->v(discretization_->vIBegin()-1, j) = vLeftBuffer[j - discretization_->vJBegin()];
        }
    }


    if(!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestvRight, MPI_STATUS_IGNORE);
        MPI_Wait(&requestuRight, MPI_STATUS_IGNORE);
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            discretization_->u(discretization_->uIEnd(), j) = uRightBuffer[j - discretization_->uJBegin()];
        }
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            discretization_->v(discretization_->vIEnd(), j) = vRightBuffer[j - discretization_->vJBegin()];
        }
    }
}