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

    //! TODO: Put into else statements afterwards, optimize sizes of buffers
    std::vector<double> topBuffer(discretization_->pIEnd() - discretization_->pIBegin(), 0.0);
    std::vector<double> bottomBuffer(discretization_->pIEnd() - discretization_->pIBegin(), 0.0);
    std::vector<double> leftBuffer(discretization_->pJEnd() - discretization_->pJBegin(), 0.0);
    std::vector<double> rightBuffer(discretization_->pJEnd() - discretization_->pJBegin(), 0.0);

    MPI_Request requestTop, requestBottom, requestLeft, requestRight;
    
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        // Set boundary values
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            discretization_->p(i, discretization_->pJEnd()) = discretization_->p(i, discretization_->pJEnd() - 1); 
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            topBuffer[i - discretization_->pIBegin()] = discretization_->p(i, discretization_->pJEnd() - 1);
        }
        partitioning_->send(topBuffer, partitioning_->topNeighbourRankNo(), requestTop);
        partitioning_->receive(topBuffer, partitioning_->topNeighbourRankNo(), requestTop);
        // partitioning_->communicate(sendTopBuffer, receiveTopBuffer, partitioning_->topNeighbourRankNo(), requestSendTop, requestReceiveTop);

    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            discretization_->p(i, discretization_->pJBegin() - 1) = discretization_->p(i, discretization_->pJBegin());
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            bottomBuffer[i - discretization_->pIBegin()] = discretization_->p(i, discretization_->pJBegin());
        }
        // partitioning_->communicate(sendBottomBuffer, receiveBottomBuffer, partitioning_->bottomNeighbourRankNo(), requestSendBottom, requestReceiveBottom); 
        partitioning_->send(bottomBuffer, partitioning_->bottomNeighbourRankNo(), requestBottom);
        partitioning_->receive(bottomBuffer, partitioning_->bottomNeighbourRankNo(), requestBottom);  
    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = discretization_->pJBegin() - 1; j < discretization_->pJEnd() + 1; j++) {
            discretization_->p(discretization_->pIBegin() - 1, j) = discretization_->p(discretization_->pIBegin(), j);
        } 
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            leftBuffer[j - discretization_->pJBegin()] = discretization_->p(discretization_->pIBegin(), j);
        }
        // partitioning_->communicate(sendLeftBuffer, receiveLeftBuffer, partitioning_->leftNeighbourRankNo(), requestSendLeft, requestReceiveLeft);
        partitioning_->send(leftBuffer, partitioning_->leftNeighbourRankNo(), requestLeft);
        partitioning_->receive(leftBuffer, partitioning_->leftNeighbourRankNo(), requestLeft); 
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = discretization_->pJBegin() - 1; j < discretization_->pJEnd() + 1; j++) {
            discretization_->p(discretization_->pIEnd(), j) = discretization_->p(discretization_->pIEnd() - 1, j);
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            rightBuffer[j - discretization_->pJBegin()] = discretization_->p(discretization_->pIEnd() - 1, j);
        }
        // partitioning_->communicate(sendRightBuffer, receiveRightBuffer, partitioning_->rightNeighbourRankNo(), requestSendRight, requestReceiveRight);
        partitioning_->send(rightBuffer, partitioning_->rightNeighbourRankNo(), requestRight);
        partitioning_->receive(rightBuffer, partitioning_->rightNeighbourRankNo(), requestRight);
    }

    // Set the buffers to the suiting columns/rows
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestTop, MPI_STATUS_IGNORE);
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            discretization_->p(i, discretization_->pJEnd()) = topBuffer[i - discretization_->pIBegin()];
        }   
    }

    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&requestBottom, MPI_STATUS_IGNORE);
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            discretization_->p(i, discretization_->pJBegin() - 1) = bottomBuffer[i - discretization_->pIBegin()];
        }
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestLeft, MPI_STATUS_IGNORE);
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            discretization_->p(discretization_->pIBegin() - 1, j) = leftBuffer[j - discretization_->pJBegin()];
        }
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestRight, MPI_STATUS_IGNORE);
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            discretization_->p(discretization_->pIEnd(), j) = rightBuffer[j - discretization_->pJBegin()];
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