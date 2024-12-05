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

//! TODO: Maybe change that u_max and v_max is communicated instead of dt_min
void ParallelComputation::computeTimeStepWidth() {
    const double dx = (*discretization_).dx();
    const double dy = (*discretization_).dy();

    const double dx2 = dx * dx;
    const double dy2 = dy * dy;

    const double dt_diffusion = (settings_.re / 2.0) * (dx2 * dy2)/(dx2 + dy2);

    // double u_max = (*partitioning_).globalMax((*discretization_).u().computeMaxAbs());
    // double v_max = (*partitioning_).globalMax((*discretization_).v().computeMaxAbs());

    
    
    // compute local timesteps
    // const double dt_convection_x = dx / u_max;
    // const double dt_convection_y = dy / v_max;
    // const double dt_convection = std::min(dt_convection_x, dt_convection_y);

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

    if ((*partitioning_).ownPartitionContainsLeftBoundary()) {
        for (int j = ((*discretization_).uJBegin() - 1); j < ((*discretization_).uJEnd() + 1); j++) {
            // left boundary
            (*discretization_).u((*discretization_).uIBegin() - 1, j) = settings_.dirichletBcLeft[0];
            (*discretization_).f((*discretization_).uIBegin() - 1, j) = settings_.dirichletBcLeft[0];
        }
    } 

    if ((*partitioning_).ownPartitionContainsRightBoundary()) {
        for (int j = ((*discretization_).uJBegin() - 1); j < ((*discretization_).uJEnd() + 1); j++) {
            // right boundary
            (*discretization_).u((*discretization_).uIEnd(), j) = settings_.dirichletBcRight[0];
            (*discretization_).f((*discretization_).uIEnd(), j) = settings_.dirichletBcRight[0];
        }
    }

    // v and g boundary values

    if ((*partitioning_).ownPartitionContainsTopBoundary()){
        for (int i = (*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++) {
            // top boundary
            (*discretization_).v(i, (*discretization_).vJEnd()) = settings_.dirichletBcTop[1];
            (*discretization_).g(i, (*discretization_).vJEnd()) = settings_.dirichletBcTop[1];
        }
    }

    if ((*partitioning_).ownPartitionContainsBottomBoundary()) {
        // bottom boundary
        for (int i = (*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++) {
            (*discretization_).v(i, (*discretization_).vJBegin() - 1) = settings_.dirichletBcBottom[1];
            (*discretization_).g(i, (*discretization_).vJBegin() - 1) = settings_.dirichletBcBottom[1];
        }
    }
}
//! TODO: Concate u and v for communication
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

    // // define u buffers
    // std::vector<double> uTopBuffer(uIEnd - uIBegin, 0.0);
    // std::vector<double> uBottomBuffer(uIEnd - uIBegin, 0.0);
    // std::vector<double> uRightBuffer(uJEnd - uJBegin, 0.0);
    // std::vector<double> uLeftBuffer(uJEnd - uJBegin, 0.0);

    // // define v buffers
    // std::vector<double> vTopBuffer(vIEnd - vIBegin, 0.0);
    // std::vector<double> vBottomBuffer(vIEnd - vIBegin, 0.0);
    // std::vector<double> vRightBuffer(vJEnd - vJBegin, 0.0);
    // std::vector<double> vLeftBuffer(vJEnd - vJBegin, 0.0);

    // define buffers
    std::vector<double> sendTopBuffer(uIEnd - uIBegin + vIEnd - vIBegin, 0.0);
    std::vector<double> sendBottomBuffer(uIEnd - uIBegin + vIEnd - vIBegin, 0.0);
    std::vector<double> sendLeftBuffer(uJEnd - uJBegin + vJEnd - vJBegin, 0.0);
    std::vector<double> sendRightBuffer(uJEnd - uJBegin + vJEnd - vJBegin, 0.0);

    // std::vector<double> recvTopBuffer(uIEnd - uIBegin + vIEnd - vIBegin, 0.0);
    // std::vector<double> recvBottomBuffer(uIEnd - uIBegin + vIEnd - vIBegin, 0.0);
    // std::vector<double> recvLeftBuffer(uJEnd - uJBegin + vJEnd - vJBegin, 0.0);
    // std::vector<double> recvRightBuffer(uJEnd - uJBegin + vJEnd - vJBegin, 0.0);


    // // define requests
    // MPI_Request requestuTop, requestuBottom, requestuLeft, requestuRight;
    // MPI_Request requestvTop, requestvBottom, requestvLeft, requestvRight;

    // define requests
    MPI_Request requestTop, requestBottom, requestLeft, requestRight;

    if ((*partitioning_).ownPartitionContainsTopBoundary()) {
        for (int i = uIBegin; i < uIEnd; i++) {
            (*discretization_).u(i, uJEnd) = 2 * settings_.dirichletBcTop[0] - (*discretization_).u(i, uJEnd - 1);
        }
    } else {
        // for (int i = uIBegin; i < uIEnd; i++) {
        //     uTopBuffer[i - uIBegin] = (*discretization_).u(i, uJEnd - 1);
        // }
        // (*partitioning_).send(uTopBuffer, (*partitioning_).topNeighbourRankNo(), requestuTop);

        // for (int i = vIBegin; i < vIEnd; i++) {
        //     vTopBuffer[i - vIBegin] = (*discretization_).v(i, vJEnd - 2);
        // }
        // (*partitioning_).send(vTopBuffer, (*partitioning_).topNeighbourRankNo(), requestvTop);

        // (*partitioning_).receive(uTopBuffer, (*partitioning_).topNeighbourRankNo(), requestuTop);
        // (*partitioning_).receive(vTopBuffer, (*partitioning_).topNeighbourRankNo(), requestvTop);


        for (int i = uIBegin; i < uIEnd; i++) {
            sendTopBuffer[i - uIBegin] = (*discretization_).u(i, uJEnd - 1);
        }

        for (int i = vIBegin; i < vIEnd; i++) {
            sendTopBuffer[i - vIBegin + uIEnd - uIBegin] = (*discretization_).v(i, vJEnd - 2);
        }

        // (*partitioning_).send(topBuffer, (*partitioning_).topNeighbourRankNo(), requestTop);
        // (*partitioning_).receive(topBuffer, (*partitioning_).topNeighbourRankNo(), requestTop);
        MPI_Isend(sendTopBuffer.data(), sendTopBuffer.size(), MPI_DOUBLE, (*partitioning_).topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestTop);
        MPI_Irecv(sendTopBuffer.data(), sendTopBuffer.size(), MPI_DOUBLE, (*partitioning_).topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestTop);
    }

    // Communication to bottom neighbour
    if ((*partitioning_).ownPartitionContainsBottomBoundary()){
        for (int i = uIBegin; i < uIEnd; i++) {
            (*discretization_).u(i, uJBegin - 1) = 2 * settings_.dirichletBcBottom[0] - (*discretization_).u(i, uJBegin);
        }
    } else {
        // for (int i = uIBegin; i < uIEnd; i++) {
        //     uBottomBuffer[i - uIBegin] = (*discretization_).u(i, uJBegin);
        // }
        // (*partitioning_).send(uBottomBuffer, (*partitioning_).bottomNeighbourRankNo(), requestuBottom);

        // for (int i = vIBegin; i < vIEnd; i++) {
        //     vBottomBuffer[i - vIBegin] = (*discretization_).v(i, vJBegin + 1);
        // }
        // (*partitioning_).send(vBottomBuffer, (*partitioning_).bottomNeighbourRankNo(), requestvBottom);

        // (*partitioning_).receive(uBottomBuffer, (*partitioning_).bottomNeighbourRankNo(), requestuBottom);
        // (*partitioning_).receive(vBottomBuffer, (*partitioning_).bottomNeighbourRankNo(), requestvBottom);
        for (int i = uIBegin; i < uIEnd; i++) {
            sendBottomBuffer[i - uIBegin] = (*discretization_).u(i, uJBegin);
        }

        for (int i = vIBegin; i < vIEnd; i++) {
            sendBottomBuffer[i - vIBegin + uIEnd - uIBegin] = (*discretization_).v(i, vJBegin + 1);
        }

        // (*partitioning_).send(bottomBuffer, (*partitioning_).bottomNeighbourRankNo(), requestBottom);
        // (*partitioning_).receive(bottomBuffer, (*partitioning_).bottomNeighbourRankNo(), requestBottom);
        MPI_Isend(sendBottomBuffer.data(), sendBottomBuffer.size(), MPI_DOUBLE, (*partitioning_).bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestBottom);
        MPI_Irecv(sendBottomBuffer.data(), sendBottomBuffer.size(), MPI_DOUBLE, (*partitioning_).bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestBottom);
    }

    // Communication to left neighbour
    if ((*partitioning_).ownPartitionContainsLeftBoundary()){
        for (int j = (vJBegin - 1); j < (vJEnd + 1); j++) {
            (*discretization_).v(vIBegin - 1, j) = 2 * settings_.dirichletBcLeft[1] -  (*discretization_).v(vIBegin, j);
        }   
    } else {
        // for (int j = uJBegin; j < uJEnd; j++) {
        //     uLeftBuffer[j - uJBegin] = (*discretization_).u(uIBegin + 1, j);
        // }
        // (*partitioning_).send(uLeftBuffer, (*partitioning_).leftNeighbourRankNo(), requestuLeft);

        // for (int j = vJBegin; j < vJEnd; j++) {
        //     vLeftBuffer[j - vJBegin] = (*discretization_).v(vIBegin, j);
        // }
        // (*partitioning_).send(vLeftBuffer, (*partitioning_).leftNeighbourRankNo(), requestvLeft);

        // (*partitioning_).receive(uLeftBuffer, (*partitioning_).leftNeighbourRankNo(), requestuLeft);
        // (*partitioning_).receive(vLeftBuffer, (*partitioning_).leftNeighbourRankNo(), requestvLeft);

        for (int j = uJBegin; j < uJEnd; j++) {
            sendLeftBuffer[j - uJBegin] = (*discretization_).u(uIBegin + 1, j);
        }

        for (int j = vJBegin; j < vJEnd; j++) {
            sendLeftBuffer[j - vJBegin + uJEnd - uJBegin] = (*discretization_).v(vIBegin, j);
        }

        // (*partitioning_).send(leftBuffer, (*partitioning_).leftNeighbourRankNo(), requestLeft);
        // (*partitioning_).receive(leftBuffer, (*partitioning_).leftNeighbourRankNo(), requestLeft);
        MPI_Isend(sendLeftBuffer.data(), sendLeftBuffer.size(), MPI_DOUBLE, (*partitioning_).leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestLeft);
        MPI_Irecv(sendLeftBuffer.data(), sendLeftBuffer.size(), MPI_DOUBLE, (*partitioning_).leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestLeft);
    }

    // Communication to right neighbour
    if ((*partitioning_).ownPartitionContainsRightBoundary()){
        for (int j = (vJBegin - 1); j < (vJEnd + 1); j++) {
            (*discretization_).v(vIEnd, j) = 2 * settings_.dirichletBcRight[1] -  (*discretization_).v(vIEnd - 1, j);
        }
    } else {
        // for (int j = uJBegin; j < uJEnd; j++) {
        //     uRightBuffer[j - uJBegin] = (*discretization_).u(uIEnd - 2, j);
        // }
        // (*partitioning_).send(uRightBuffer, (*partitioning_).rightNeighbourRankNo(), requestuRight);

        // for (int j = vJBegin; j < vJEnd; j++) {
        //     vRightBuffer[j - vJBegin] = (*discretization_).v(vIEnd - 1, j);
        // }
        // (*partitioning_).send(vRightBuffer, (*partitioning_).rightNeighbourRankNo(), requestvRight);

        // (*partitioning_).receive(uRightBuffer, (*partitioning_).rightNeighbourRankNo(), requestuRight);
        // (*partitioning_).receive(vRightBuffer, (*partitioning_).rightNeighbourRankNo(), requestvRight);
        for (int j = uJBegin; j < uJEnd; j++) {
            sendRightBuffer[j - uJBegin] = (*discretization_).u(uIEnd - 2, j);
        }

        for (int j = vJBegin; j < vJEnd; j++) {
            sendRightBuffer[j - vJBegin + uJEnd - uJBegin] = (*discretization_).v(vIEnd - 1, j);
        }

        // (*partitioning_).send(rightBuffer, (*partitioning_).rightNeighbourRankNo(), requestRight);
        // (*partitioning_).receive(rightBuffer, (*partitioning_).rightNeighbourRankNo(), requestRight);
        MPI_Isend(sendRightBuffer.data(), sendRightBuffer.size(), MPI_DOUBLE, (*partitioning_).rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestRight);
        MPI_Irecv(sendRightBuffer.data(), sendRightBuffer.size(), MPI_DOUBLE, (*partitioning_).rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestRight);
    }

    // Set the received values
    if (!(*partitioning_).ownPartitionContainsTopBoundary()) {
        // MPI_Wait(&requestvTop, MPI_STATUS_IGNORE);
        // MPI_Wait(&requestuTop, MPI_STATUS_IGNORE);
        MPI_Wait(&requestTop, MPI_STATUS_IGNORE);
        for (int i = uIBegin; i < uIEnd; i++) {
            (*discretization_).u(i, uJEnd) = sendTopBuffer[i - uIBegin];
        }
        for (int i = vIBegin; i < vIEnd; i++) {
            (*discretization_).v(i, vJEnd) = sendTopBuffer[i - vIBegin + uIEnd - uIBegin];
        }
    }

    if(!(*partitioning_).ownPartitionContainsBottomBoundary()) {
        // MPI_Wait(&requestvBottom, MPI_STATUS_IGNORE);
        // MPI_Wait(&requestuBottom, MPI_STATUS_IGNORE);
        MPI_Wait(&requestBottom, MPI_STATUS_IGNORE);
        for (int i = uIBegin; i < uIEnd; i++) {
            (*discretization_).u(i, uJBegin -1) = sendBottomBuffer[i - uIBegin];
        }
        for (int i = vIBegin; i < vIEnd; i++) {
            (*discretization_).v(i, vJBegin -1) = sendBottomBuffer[i - vIBegin + uIEnd - uIBegin];
        }
    }

    if(!(*partitioning_).ownPartitionContainsLeftBoundary()) {
        // MPI_Wait(&requestvLeft, MPI_STATUS_IGNORE);
        // MPI_Wait(&requestuLeft, MPI_STATUS_IGNORE);
        MPI_Wait(&requestLeft, MPI_STATUS_IGNORE);
        for (int j = uJBegin; j < uJEnd; j++) {
            (*discretization_).u(uIBegin-1, j) = sendLeftBuffer[j - uJBegin];
        }
        for (int j = vJBegin; j < vJEnd; j++) {
            (*discretization_).v(vIBegin-1, j) = sendLeftBuffer[j - vJBegin + uJEnd - uJBegin];
        }
    }


    if(!(*partitioning_).ownPartitionContainsRightBoundary()) {
        // MPI_Wait(&requestvRight, MPI_STATUS_IGNORE);
        // MPI_Wait(&requestuRight, MPI_STATUS_IGNORE);
        MPI_Wait(&requestRight, MPI_STATUS_IGNORE);
        for (int j = uJBegin; j < uJEnd; j++) {
            (*discretization_).u(uIEnd, j) = sendRightBuffer[j - uJBegin];
        }
        for (int j = vJBegin; j < vJEnd; j++) {
            (*discretization_).v(vIEnd, j) = sendRightBuffer[j - vJBegin + uJEnd - uJBegin];
        }
    }
}