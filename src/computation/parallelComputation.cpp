#include "parallelComputation.h"
#include <cmath>

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
        discretization_ = std::make_shared<DonorCell>(localCells, meshWidth_, settings_.alpha);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(localCells, meshWidth_);
    }

    // init output writers
    // outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
    // outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discretization_, partitioning_);

    // init pressure solvers
    if (settings_.pressureSolver == "SOR") {
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
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
    int t_iter = 0;
    double time_epsilon = 1e-8;
    
    // Loop over all time steps until t_end is reached
    while (time < (settings_.endTime - time_epsilon)) {

        applyBoundaryValues();
        std::cout << "1\n" << std::endl;
        
        computeTimeStepWidth();
        std::cout << "12\n" << std::endl;

        // Check to choose the last time step differently to hit t_end exactly
        if (time + dt_ > settings_.endTime - time_epsilon) {
            dt_ = settings_.endTime - time;
        }
        time += dt_;


        computePreliminaryVelocities();
        std::cout << "123\n" << std::endl;

        communicatePreliminaryVelocities();
        std::cout << "1234\n" << std::endl;

        computeRightHandSide();
        std::cout << "12345\n" << std::endl;

        computePressure();
        std::cout << "123456\n" << std::endl;

        computeVelocities();
        std::cout << "1234567\n" << std::endl;

        // outputWriterParaview_->writeFile(time); // Output
        // outputWriterText_->writeFile(time); // Output
    }
}

//! TODO: Maybe change that u_max and v_max is communicated instead of dt_min
void ParallelComputation::computeTimeStepWidth() {
    const double dx = discretization_->dx();
    const double dy = discretization_->dy();

    const double dx2 = dx * dx;
    const double dy2 = dy * dy;

    const double dt_diffusion = (settings_.re / 2.0) * (dx2 * dy2)/(dx2 + dy2);
    
    // compute local timesteps
    const double dt_convection_x = dx / discretization_->u().computeMaxAbs();
    const double dt_convection_y = dy / discretization_->v().computeMaxAbs();
    const double dt_convection_local = std::min(dt_convection_x, dt_convection_y);

    std::cout << "dt_convection_local: " << dt_convection_local << std::endl;


    // compute global timestep
    double dt_convection = partitioning_->globalMin(dt_convection_local);


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
    std::vector<double> senduTopBuffer(discretization_->uIEnd() - discretization_->uIBegin(), 0);
    std::vector<double> senduBottomBuffer(discretization_->uIEnd() - discretization_->uIBegin(), 0);
    std::vector<double> senduRightBuffer(discretization_->uJEnd() - discretization_->uJBegin(), 0);
    std::vector<double> senduLeftBuffer(discretization_->uJEnd() - discretization_->uJBegin(), 0);

    std::vector<double> receiveuTopBuffer(discretization_->uIEnd() - discretization_->uIBegin(), 0);
    std::vector<double> receiveuBottomBuffer(discretization_->uIEnd() - discretization_->uIBegin(), 0);
    std::vector<double> receiveuLeftBuffer(discretization_->uJEnd() - discretization_->uJBegin(), 0);
    std::vector<double> receiveuRightBuffer(discretization_->uJEnd() - discretization_->uJBegin(), 0);

    // define v buffers
    std::vector<double> sendvTopBuffer(discretization_->vIEnd() - discretization_->vIBegin(), 0);
    std::vector<double> sendvBottomBuffer(discretization_->vIEnd() - discretization_->vIBegin(), 0);
    std::vector<double> sendvRightBuffer(discretization_->vJEnd() - discretization_->vJBegin(), 0);
    std::vector<double> sendvLeftBuffer(discretization_->vJEnd() - discretization_->vJBegin(), 0);

    std::vector<double> receivevTopBuffer(discretization_->vIEnd() - discretization_->vIBegin(), 0);
    std::vector<double> receivevBottomBuffer(discretization_->vIEnd() - discretization_->vIBegin(), 0);
    std::vector<double> receivevLeftBuffer(discretization_->vJEnd() - discretization_->vJBegin(), 0);
    std::vector<double> receivevRightBuffer(discretization_->vJEnd() - discretization_->vJBegin(), 0);

    // define requests
    MPI_Request requestuSendTop, requestuSendBottom, requestuSendLeft, requestuSendRight;
    MPI_Request requestuReceiveTop, requestuReceiveBottom, requestuReceiveLeft, requestuReceiveRight;

    MPI_Request requestvSendTop, requestvSendBottom, requestvSendLeft, requestvSendRight;
    MPI_Request requestvReceiveTop, requestvReceiveBottom, requestvReceiveLeft, requestvReceiveRight;

    std::vector<MPI_Request> requests = {requestuSendTop, requestuSendBottom, requestuSendLeft, requestuSendRight, requestuReceiveTop, requestuReceiveBottom, requestuReceiveLeft, requestuReceiveRight, requestvSendTop, requestvSendBottom, requestvSendLeft, requestvSendRight, requestvReceiveTop, requestvReceiveBottom, requestvReceiveLeft, requestvReceiveRight};

    // Communication to top neighbour
    if (partitioning_->ownPartitionContainsTopBoundary()){
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
            discretization_->u(i, discretization_->uJEnd()) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 1);
        }
    } else {
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
            senduTopBuffer[i - discretization_->uIBegin()] = discretization_->u(i, discretization_->uJEnd() - 1);
        }
        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
            sendvTopBuffer[i - discretization_->vIBegin()] = discretization_->v(i, discretization_->vJEnd() - 1);
        }
        partitioning_->communicate(senduTopBuffer, receiveuTopBuffer, partitioning_->topNeighbourRankNo(), requestuSendTop, requestuReceiveTop);
        partitioning_->communicate(sendvTopBuffer, receivevTopBuffer, partitioning_->topNeighbourRankNo(), requestvSendTop, requestvReceiveTop);  
    }

    // Communication to bottom neighbour
    if (partitioning_->ownPartitionContainsBottomBoundary()){
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
            discretization_->u(i, discretization_->uJBegin() - 1) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin());
        }
    } else {
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
            senduBottomBuffer[i - discretization_->uIBegin()] = discretization_->u(i, discretization_->uJBegin());
        }
        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
            sendvBottomBuffer[i - discretization_->vIBegin()] = discretization_->v(i, discretization_->vJBegin());
        }
        partitioning_->communicate(senduBottomBuffer, receiveuBottomBuffer, partitioning_->bottomNeighbourRankNo(), requestuSendBottom, requestuReceiveBottom);
        partitioning_->communicate(sendvBottomBuffer, receivevBottomBuffer, partitioning_->bottomNeighbourRankNo(), requestvSendBottom, requestvReceiveBottom);  
    }

    // Communication to left neighbour
    if (partitioning_->ownPartitionContainsLeftBoundary()){
        for (int j = (discretization_->vJBegin() - 1); j < (discretization_->vJEnd() + 1); j++) {
            discretization_->v(discretization_->vIBegin() - 1, j) = 2 * settings_.dirichletBcLeft[1] -  discretization_->v(discretization_->vIBegin(), j);
        }   
    } else {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            senduLeftBuffer[j - discretization_->uJBegin()] = discretization_->u(discretization_->uIBegin(), j);
        }
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            sendvLeftBuffer[j - discretization_->vJBegin()] = discretization_->v(discretization_->vIBegin(), j);
        }
        partitioning_->communicate(senduLeftBuffer, receiveuLeftBuffer, partitioning_->leftNeighbourRankNo(), requestuSendLeft, requestuReceiveLeft);
        partitioning_->communicate(sendvLeftBuffer, receivevLeftBuffer, partitioning_->leftNeighbourRankNo(), requestvSendLeft, requestvReceiveLeft);   
    }

    // Communication to right neighbour
    if (partitioning_->ownPartitionContainsRightBoundary()){
        for (int j = (discretization_->vJBegin() - 1); j < (discretization_->vJEnd() + 1); j++) {
            discretization_->v(discretization_->vIEnd(), j) = 2 * settings_.dirichletBcRight[1] -  discretization_->v(discretization_->vIEnd() - 1, j);
        }
    } else {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            senduRightBuffer[j - discretization_->uJBegin()] = discretization_->u(discretization_->uIEnd(), j);
        }
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            sendvRightBuffer[j - discretization_->vJBegin()] = discretization_->v(discretization_->vIEnd(), j);
        }
        partitioning_->communicate(senduRightBuffer, receiveuRightBuffer, partitioning_->rightNeighbourRankNo(), requestuSendRight, requestuReceiveRight);
        partitioning_->communicate(sendvRightBuffer, receivevRightBuffer, partitioning_->rightNeighbourRankNo(), requestvSendRight, requestvReceiveRight);   
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    std::cout << receiveuRightBuffer.size() << std::endl;

    // Set the received values
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
            discretization_->u(i, discretization_->uJBegin() - 1) = receiveuBottomBuffer[i - discretization_->uIBegin()];
            discretization_->u(i, discretization_->uJEnd()) = receiveuTopBuffer[i - discretization_->uIBegin()];
        }
    }

    if(!partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            discretization_->u(discretization_->pIBegin() - 1, j) = receiveuLeftBuffer[j - discretization_->uJBegin()];
            discretization_->u(discretization_->uIEnd(), j) = receiveuRightBuffer[j - discretization_->uJBegin()];
        }
    }

    if(!partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
            discretization_->v(i, discretization_->vJBegin() - 1) = receivevBottomBuffer[i - discretization_->vIBegin()];
            discretization_->v(i, discretization_->vJEnd()) = receivevTopBuffer[i - discretization_->vIBegin()];
        }
    }


    if(!partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            discretization_->v(discretization_->vIBegin() - 1, j) = receivevLeftBuffer[j - discretization_->vJBegin()];
            discretization_->v(discretization_->vIEnd(), j) = receivevRightBuffer[j - discretization_->vJBegin()];
        }
    }


}

void ParallelComputation::communicatePreliminaryVelocities() {
    // define buffers
    std::vector<double> sendTopBuffer(discretization_->vIEnd() - discretization_->vIBegin(), 0);
    std::vector<double> sendRightBuffer(discretization_->uJEnd() - discretization_->uJBegin(), 0);

    std::vector<double> receiveBottomBuffer(discretization_->vIEnd() - discretization_->vIBegin(), 0);
    std::vector<double> receiveLeftBuffer(discretization_->uJEnd() - discretization_->uJBegin(), 0);

    MPI_Request requestSendTop, requestSendRight;
    MPI_Request requestReceiveBottom, requestReceiveLeft;
    std::vector<MPI_Request> requests = {requestSendTop, requestSendRight, requestReceiveBottom, requestReceiveLeft};

    // Send values to top neighbour
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
            sendTopBuffer[i - discretization_->vIBegin()] = discretization_->g(i, discretization_->vJEnd());
        }
        partitioning_->send(sendTopBuffer, partitioning_->topNeighbourRankNo(), requestSendTop);
    }
    // send values to right neighbour
    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            sendRightBuffer[j - discretization_->uJBegin()] = discretization_->f(discretization_->uIEnd(), j);
        }
        partitioning_->send(sendRightBuffer, partitioning_->rightNeighbourRankNo(), requestSendRight);
    }


    // receive value from bottom neighbour
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        partitioning_->receive(receiveBottomBuffer, partitioning_->bottomNeighbourRankNo(), requestReceiveBottom);
        MPI_Wait(&requestReceiveBottom, MPI_STATUS_IGNORE);
        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
            discretization_->g(i, discretization_->vJBegin() - 1) = receiveBottomBuffer[i - discretization_->vIBegin()];
        }
    }

    // receive value from left neighbour
    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        partitioning_->receive(receiveLeftBuffer, partitioning_->leftNeighbourRankNo(), requestReceiveLeft);
        MPI_Wait(&requestReceiveLeft, MPI_STATUS_IGNORE);
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            discretization_->f(discretization_->uIBegin() - 1, j) = receiveLeftBuffer[j - discretization_->uJBegin()];
        }
    } 
}
