#include "parallelComputation.h"
#include <cmath>
#include <unistd.h>


void ParallelComputation::initialize(int argc, char* argv[]) {
    // load settigns from file
    std::string filename = argv[1];
    settings_.loadFromFile(filename);

    meshWidth_ = {settings_.physicalSize[0] / settings_.nCells[0], settings_.physicalSize[1] / settings_.nCells[1]};

    partitioning_ = std::make_shared<Partitioning>();
    (*partitioning_).initialize(settings_.nCells);
    std::array<int,2> localCells = (*partitioning_).nCellsLocal();

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
    // if (settings_.pressureSolver == "SOR") {
    //     pressureSolver_ = std::make_unique<ParallelSOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega, partitioning_);
    // } else if (settings_.pressureSolver == "GaussSeidel") {
    //     pressureSolver_ = std::make_unique<ParallelGaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_);
    // } else if (settings_.pressureSolver == "CG") {
    //     pressureSolver_ = std::make_unique<ParallelCG>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_);
    // } else {
    //     std::cerr << "Unknown pressure solver: " << settings_.pressureSolver << std::endl;
    //     std::exit(1);
    // }

    // hard coded CG
    pressureSolver_ = std::make_unique<ParallelCG>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_);
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

    	
        computePreliminaryVelocities();
    

        computeRightHandSide();
        

        computePressure();
        

        computeVelocities();

        // Output
        if (time >= output) {
            (*outputWriterParaview_).writeFile(time); // Output
            // outputWriterText_->writeFile(time); // Output
            output++;
        }
        
    }
}

void ParallelComputation::computeTimeStepWidth() {
    const double dx = (*discretization_).dx();
    const double dy = (*discretization_).dy();

    const double dx2 = dx * dx;
    const double dy2 = dy * dy;

    const double dt_diffusion = (settings_.re / 2.0) * (dx2 * dy2)/(dx2 + dy2);

    const double dt_convection_x = dx / (*discretization_).u().computeMaxAbs();
    const double dt_convection_y = dy / (*discretization_).v().computeMaxAbs();
    const double dt_convection_local = std::min(dt_convection_x, dt_convection_y);

    MPI_Request timerequest;

    double dt_convection = 0.0;

    MPI_Iallreduce(&dt_convection_local, &dt_convection, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, &timerequest);
    MPI_Wait(&timerequest, MPI_STATUS_IGNORE);

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
    // left boundary
    if ((*partitioning_).ownPartitionContainsLeftBoundary()) {
        for (int j = ((*discretization_).uJBegin() - 1); j < ((*discretization_).uJEnd() + 1); j++) {
            (*discretization_).u((*discretization_).uIBegin() - 1, j) = settings_.dirichletBcLeft[0];
            (*discretization_).f((*discretization_).uIBegin() - 1, j) = settings_.dirichletBcLeft[0];
        }
    } 

    // right boundary
    if ((*partitioning_).ownPartitionContainsRightBoundary()) {
        for (int j = ((*discretization_).uJBegin() - 1); j < ((*discretization_).uJEnd() + 1); j++) {
            (*discretization_).u((*discretization_).uIEnd(), j) = settings_.dirichletBcRight[0];
            (*discretization_).f((*discretization_).uIEnd(), j) = settings_.dirichletBcRight[0];
        }
    }

    // v and g boundary values
    // top boundary
    if ((*partitioning_).ownPartitionContainsTopBoundary()){
        for (int i = (*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++) {
            (*discretization_).v(i, (*discretization_).vJEnd()) = settings_.dirichletBcTop[1];
            (*discretization_).g(i, (*discretization_).vJEnd()) = settings_.dirichletBcTop[1];
        }
    }

    // bottom boundary
    if ((*partitioning_).ownPartitionContainsBottomBoundary()) {
        for (int i = (*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++) {
            (*discretization_).v(i, (*discretization_).vJBegin() - 1) = settings_.dirichletBcBottom[1];
            (*discretization_).g(i, (*discretization_).vJBegin() - 1) = settings_.dirichletBcBottom[1];
        }
    }
}

void ParallelComputation::applyBoundaryValues() {
    // define loop limits
    int uIBegin = (*discretization_).uIBegin();
    int uIEnd = (*discretization_).uIEnd();
    int uJBegin = (*discretization_).uJBegin();
    int uJEnd = (*discretization_).uJEnd();

    int vIBegin = (*discretization_).vIBegin();
    int vIEnd = (*discretization_).vIEnd();
    int vJBegin = (*discretization_).vJBegin();
    int vJEnd = (*discretization_).vJEnd();

    // define buffers
    std::vector<double> topBuffer(uIEnd - uIBegin + vIEnd - vIBegin, 0.0);
    std::vector<double> bottomBuffer(uIEnd - uIBegin + vIEnd - vIBegin, 0.0);
    std::vector<double> leftBuffer(uJEnd - uJBegin + vJEnd - vJBegin, 0.0);
    std::vector<double> rightBuffer(uJEnd - uJBegin + vJEnd - vJBegin, 0.0);

    // define requests
    MPI_Request requestTop, requestBottom, requestLeft, requestRight;

    // top boundary condition
    if ((*partitioning_).ownPartitionContainsTopBoundary()) {
        for (int i = uIBegin; i < uIEnd; i++) {
            (*discretization_).u(i, uJEnd) = 2 * settings_.dirichletBcTop[0] - (*discretization_).u(i, uJEnd - 1);
        }
    } else { // communication to top neighbor
        for (int i = uIBegin; i < uIEnd; i++) {
            topBuffer[i - uIBegin] = (*discretization_).u(i, uJEnd - 1);
        }

        for (int i = vIBegin; i < vIEnd; i++) {
            topBuffer[i - vIBegin + uIEnd - uIBegin] = (*discretization_).v(i, vJEnd - 2);
        }

        MPI_Isend(topBuffer.data(), topBuffer.size(), MPI_DOUBLE, (*partitioning_).topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestTop);
        MPI_Irecv(topBuffer.data(), topBuffer.size(), MPI_DOUBLE, (*partitioning_).topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestTop);
    }

    // Bottom boundary condition
    if ((*partitioning_).ownPartitionContainsBottomBoundary()){
        for (int i = uIBegin; i < uIEnd; i++) {
            (*discretization_).u(i, uJBegin - 1) = 2 * settings_.dirichletBcBottom[0] - (*discretization_).u(i, uJBegin);
        }
    } else { // communication to bottom neighbor
        for (int i = uIBegin; i < uIEnd; i++) {
            bottomBuffer[i - uIBegin] = (*discretization_).u(i, uJBegin);
        }

        for (int i = vIBegin; i < vIEnd; i++) {
            bottomBuffer[i - vIBegin + uIEnd - uIBegin] = (*discretization_).v(i, vJBegin + 1);
        }

        MPI_Isend(bottomBuffer.data(), bottomBuffer.size(), MPI_DOUBLE, (*partitioning_).bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestBottom);
        MPI_Irecv(bottomBuffer.data(), bottomBuffer.size(), MPI_DOUBLE, (*partitioning_).bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestBottom);
    }

    // Left boundary condition
    if ((*partitioning_).ownPartitionContainsLeftBoundary()){
        for (int j = (vJBegin - 1); j < (vJEnd + 1); j++) {
            (*discretization_).v(vIBegin - 1, j) = 2 * settings_.dirichletBcLeft[1] -  (*discretization_).v(vIBegin, j);
        }   
    } else { // communication to left neighbor
        for (int j = uJBegin; j < uJEnd; j++) {
            leftBuffer[j - uJBegin] = (*discretization_).u(uIBegin + 1, j);
        }

        for (int j = vJBegin; j < vJEnd; j++) {
            leftBuffer[j - vJBegin + uJEnd - uJBegin] = (*discretization_).v(vIBegin, j);
        }

        MPI_Isend(leftBuffer.data(), leftBuffer.size(), MPI_DOUBLE, (*partitioning_).leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestLeft);
        MPI_Irecv(leftBuffer.data(), leftBuffer.size(), MPI_DOUBLE, (*partitioning_).leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestLeft);
    }

    // Right boundary condition
    if ((*partitioning_).ownPartitionContainsRightBoundary()){
        for (int j = (vJBegin - 1); j < (vJEnd + 1); j++) {
            (*discretization_).v(vIEnd, j) = 2 * settings_.dirichletBcRight[1] -  (*discretization_).v(vIEnd - 1, j);
        }
    } else { // communication to right neighbor
        for (int j = uJBegin; j < uJEnd; j++) {
            rightBuffer[j - uJBegin] = (*discretization_).u(uIEnd - 2, j);
        }

        for (int j = vJBegin; j < vJEnd; j++) {
            rightBuffer[j - vJBegin + uJEnd - uJBegin] = (*discretization_).v(vIEnd - 1, j);
        }

        MPI_Isend(rightBuffer.data(), rightBuffer.size(), MPI_DOUBLE, (*partitioning_).rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestRight);
        MPI_Irecv(rightBuffer.data(), rightBuffer.size(), MPI_DOUBLE, (*partitioning_).rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestRight);
    }

    // Set the received values
    if (!(*partitioning_).ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestTop, MPI_STATUS_IGNORE);

        for (int i = uIBegin; i < uIEnd; i++) {
            (*discretization_).u(i, uJEnd) = topBuffer[i - uIBegin];
        }

        for (int i = vIBegin; i < vIEnd; i++) {
            (*discretization_).v(i, vJEnd) = topBuffer[i - vIBegin + uIEnd - uIBegin];
        }
    }

    if(!(*partitioning_).ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&requestBottom, MPI_STATUS_IGNORE);

        for (int i = uIBegin; i < uIEnd; i++) {
            (*discretization_).u(i, uJBegin -1) = bottomBuffer[i - uIBegin];
        }

        for (int i = vIBegin; i < vIEnd; i++) {
            (*discretization_).v(i, vJBegin -1) = bottomBuffer[i - vIBegin + uIEnd - uIBegin];
        }
    }

    if(!(*partitioning_).ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestLeft, MPI_STATUS_IGNORE);

        for (int j = uJBegin; j < uJEnd; j++) {
            (*discretization_).u(uIBegin-1, j) = leftBuffer[j - uJBegin];
        }

        for (int j = vJBegin; j < vJEnd; j++) {
            (*discretization_).v(vIBegin-1, j) = leftBuffer[j - vJBegin + uJEnd - uJBegin];
        }
    }


    if(!(*partitioning_).ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestRight, MPI_STATUS_IGNORE);

        for (int j = uJBegin; j < uJEnd; j++) {
            (*discretization_).u(uIEnd, j) = rightBuffer[j - uJBegin];
        }

        for (int j = vJBegin; j < vJEnd; j++) {
            (*discretization_).v(vIEnd, j) = rightBuffer[j - vJBegin + uJEnd - uJBegin];
        }
    }
}