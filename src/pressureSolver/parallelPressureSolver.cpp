#include "pressureSolver/parallelPressureSolver.h"
#include "pressureSolver/pressureSolver.h"
#include <cmath>
#include "discretization/discretization.h"
#include "partitioning/partitioning.h"
#include <iostream>


ParallelPressureSolver::ParallelPressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), 
    partitioning_(partitioning) {}

void ParallelPressureSolver::communicateAndBoundaries() {

    // Initialize buffers and loop limits

    int pIBegin = discretization_->pIBegin();
    int pIEnd = discretization_->pIEnd();
    int pJBegin = discretization_->pJBegin();
    int pJEnd = discretization_->pJEnd();

    std::vector<double> sendTopBuffer(pIEnd - pIBegin, 0.0);
    std::vector<double> sendBottomBuffer(pIEnd - pIBegin, 0.0);
    std::vector<double> sendLeftBuffer(pJEnd - pJBegin, 0.0);
    std::vector<double> sendRightBuffer(pJEnd - pJBegin, 0.0);

    std::vector<double> recvTopBuffer(pIEnd - pIBegin, 0.0);
    std::vector<double> recvBottomBuffer(pIEnd - pIBegin, 0.0);
    std::vector<double> recvLeftBuffer(pJEnd - pJBegin, 0.0);
    std::vector<double> recvRightBuffer(pJEnd - pJBegin, 0.0);


    

    MPI_Request requestTop, requestBottom, requestLeft, requestRight;

    
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        // Set boundary values
        for (int i = pIBegin; i < pIEnd; i++) {
            discretization_->p(i, pJEnd) = discretization_->p(i, pJEnd - 1); 
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int i = pIBegin; i < pIEnd; i++) {
            sendTopBuffer[i - pIBegin] = discretization_->p(i, pJEnd - 1);
        }
        // partitioning_->send(topBuffer, partitioning_->topNeighbourRankNo(), requestTop);
        // partitioning_->receive(topBuffer, partitioning_->topNeighbourRankNo(), requestTop);
        MPI_Isend(sendTopBuffer.data(), sendTopBuffer.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestTop);
        MPI_Irecv(recvTopBuffer.data(), recvTopBuffer.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestTop);
        // partitioning_->communicate(sendTopBuffer, receiveTopBuffer, partitioning_->topNeighbourRankNo(), requestSendTop, requestReceiveTop);

    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = pIBegin; i < pIEnd; i++) {
            discretization_->p(i, pJBegin - 1) = discretization_->p(i, pJBegin);
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int i = pIBegin; i < pIEnd; i++) {
            sendBottomBuffer[i - pIBegin] = discretization_->p(i, pJBegin);
        }
        // partitioning_->communicate(sendBottomBuffer, receiveBottomBuffer, partitioning_->bottomNeighbourRankNo(), requestSendBottom, requestReceiveBottom); 
        // partitioning_->send(bottomBuffer, partitioning_->bottomNeighbourRankNo(), requestBottom);
        // partitioning_->receive(bottomBuffer, partitioning_->bottomNeighbourRankNo(), requestBottom);  
        MPI_Isend(sendBottomBuffer.data(), sendBottomBuffer.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestBottom);
        MPI_Irecv(recvBottomBuffer.data(), recvBottomBuffer.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestBottom);
    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = pJBegin - 1; j < pJEnd + 1; j++) {
            discretization_->p(pIBegin - 1, j) = discretization_->p(pIBegin, j);
        } 
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int j = pJBegin; j < pJEnd; j++) {
            sendLeftBuffer[j - pJBegin] = discretization_->p(pIBegin, j);
        }
        // partitioning_->communicate(sendLeftBuffer, receiveLeftBuffer, partitioning_->leftNeighbourRankNo(), requestSendLeft, requestReceiveLeft);
        // partitioning_->send(leftBuffer, partitioning_->leftNeighbourRankNo(), requestLeft);
        // partitioning_->receive(leftBuffer, partitioning_->leftNeighbourRankNo(), requestLeft); 
        MPI_Isend(sendLeftBuffer.data(), sendLeftBuffer.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestLeft);
        MPI_Irecv(recvLeftBuffer.data(), recvLeftBuffer.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestLeft);
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = pJBegin - 1; j < pJEnd + 1; j++) {
            discretization_->p(pIEnd, j) = discretization_->p(pIEnd - 1, j);
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int j = pJBegin; j < pJEnd; j++) {
            sendRightBuffer[j - pJBegin] = discretization_->p(pIEnd - 1, j);
        }
        // partitioning_->communicate(sendRightBuffer, receiveRightBuffer, partitioning_->rightNeighbourRankNo(), requestSendRight, requestReceiveRight);
        // partitioning_->send(rightBuffer, partitioning_->rightNeighbourRankNo(), requestRight);
        // partitioning_->receive(rightBuffer, partitioning_->rightNeighbourRankNo(), requestRight);
        MPI_Isend(sendRightBuffer.data(), sendRightBuffer.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestRight);
        MPI_Irecv(recvRightBuffer.data(), recvRightBuffer.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestRight);
    }

    // Set the buffers to the suiting columns/rows
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestTop, MPI_STATUS_IGNORE);
        for (int i = pIBegin; i < pIEnd; i++) {
            discretization_->p(i, pJEnd) = recvTopBuffer[i - pIBegin];
        }   
    }

    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&requestBottom, MPI_STATUS_IGNORE);
        for (int i = pIBegin; i < pIEnd; i++) {
            discretization_->p(i, pJBegin - 1) = recvBottomBuffer[i - pIBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestLeft, MPI_STATUS_IGNORE);
        for (int j = pJBegin; j < pJEnd; j++) {
            discretization_->p(pIBegin - 1, j) = recvLeftBuffer[j - pJBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestRight, MPI_STATUS_IGNORE);
        for (int j = pJBegin; j < pJEnd; j++) {
            discretization_->p(pIEnd, j) = recvRightBuffer[j - pJBegin];
        }
    }

}


// TODO: Change /N factor for the subdomains in total  
void ParallelPressureSolver::computeResidualNorm() {
    // Define needed parameters
    double residualNorm2 = 0.0;
    const double dx_2 = discretization_->dx() * discretization_->dx();
    const double dy_2 = discretization_->dy() * discretization_->dy();
    const int N = partitioning_->nCellsGlobal()[0] * partitioning_->nCellsGlobal()[1]; // Number of points used to compute norm

    // Compute l^2 norm of residual
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            const double p_xx = (discretization_->p(i + 1, j) - 2.0 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / dx_2;
            const double p_yy = (discretization_->p(i, j + 1) - 2.0 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / dy_2;
            residualNorm2 += (discretization_->rhs(i, j) - p_xx - p_yy) * (discretization_->rhs(i, j) - p_xx - p_yy);
        }
    }
    residualNorm2_ = partitioning_->globalSum(residualNorm2) / N;
}   