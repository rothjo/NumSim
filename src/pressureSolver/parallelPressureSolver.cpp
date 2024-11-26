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
    std::vector<double> sendTopBuffer(discretization_->pIEnd() - discretization_->pIBegin());
    std::vector<double> sendBottomBuffer(discretization_->pIEnd() - discretization_->pIBegin());
    std::vector<double> sendLeftBuffer(discretization_->pJEnd() - discretization_->pJBegin());
    std::vector<double> sendRightBuffer(discretization_->pJEnd() - discretization_->pJBegin());

    std::vector<double> receiveTopBuffer(discretization_->pIEnd() - discretization_->pIBegin(), 0);
    std::vector<double> receiveBottomBuffer(discretization_->pIEnd() - discretization_->pIBegin(), 0);
    std::vector<double> receiveLeftBuffer(discretization_->pJEnd() - discretization_->pJBegin(), 0);
    std::vector<double> receiveRightBuffer(discretization_->pJEnd() - discretization_->pJBegin(), 0);

    MPI_Request requestSendTop, requestSendBottom, requestSendLeft, requestSendRight;
    MPI_Request requestReceiveTop, requestReceiveBottom, requestReceiveLeft, requestReceiveRight;
    std::vector<MPI_Request> requests = {requestSendTop, requestSendBottom, requestSendLeft, requestSendRight, requestReceiveTop, requestReceiveBottom, requestReceiveLeft, requestReceiveRight};


    if (partitioning_->ownPartitionContainsTopBoundary()) {
        // Set boundary values
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            discretization_->p(i, discretization_->pJEnd()) = discretization_->p(i, discretization_->pJEnd() - 1);
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            sendTopBuffer[i - discretization_->pIBegin()] = discretization_->p(i, discretization_->pJEnd() - 1);
        }
        partitioning_->communicate(sendTopBuffer, receiveTopBuffer, partitioning_->topNeighbourRankNo(), requestSendTop, requestReceiveTop);   

    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            discretization_->p(i, discretization_->pJBegin() - 1) = discretization_->p(i, discretization_->pJBegin());
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            sendBottomBuffer[i - discretization_->pIBegin()] = discretization_->p(i, discretization_->pJBegin());
        }
        partitioning_->communicate(sendBottomBuffer, receiveBottomBuffer, partitioning_->bottomNeighbourRankNo(), requestSendBottom, requestReceiveBottom);   

    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = discretization_->pJBegin() - 1; j < discretization_->pJEnd() + 1; j++) {
            discretization_->p(discretization_->pIBegin() - 1, j) = discretization_->p(discretization_->pIBegin(), j);
        } 
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            sendLeftBuffer[j - discretization_->pJBegin()] = discretization_->p(discretization_->pIBegin(), j);
        }
        partitioning_->communicate(sendLeftBuffer, receiveLeftBuffer, partitioning_->leftNeighbourRankNo(), requestSendLeft, requestReceiveLeft);   
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = discretization_->pJBegin() - 1; j < discretization_->pJEnd() + 1; j++) {
            discretization_->p(discretization_->pIEnd(), j) = discretization_->p(discretization_->pIEnd() - 1, j);
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            sendRightBuffer[j - discretization_->pJBegin()] = discretization_->p(discretization_->pIEnd() - 1, j);
        }
        partitioning_->communicate(sendRightBuffer, receiveRightBuffer, partitioning_->rightNeighbourRankNo(), requestSendRight, requestReceiveRight);   
    }

    // Set the buffers to the suiting columns/rows
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestReceiveTop, MPI_STATUS_IGNORE);
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            discretization_->p(i, discretization_->pJEnd()) = receiveTopBuffer[i - discretization_->pIBegin()];
        }   
    }

    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestReceiveBottom, MPI_STATUS_IGNORE);
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            discretization_->p(i, discretization_->pJBegin() - 1) = receiveBottomBuffer[i - discretization_->pIBegin()];
        }
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestReceiveLeft, MPI_STATUS_IGNORE);
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            discretization_->p(discretization_->pIBegin() - 1, j) = receiveLeftBuffer[j - discretization_->pJBegin()];
        }
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestReceiveRight, MPI_STATUS_IGNORE);
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            discretization_->p(discretization_->pIEnd(), j) = receiveRightBuffer[j - discretization_->pJBegin()];
        }
    }

    //! TODO: Move to the place where needed
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);

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
    std::cout  << " Residual: " << residualNorm2 << std::endl;
    std::cout << "N: " << N << std::endl;
    residualNorm2_ = partitioning_->globalSum(residualNorm2) / N;
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Residual norm: " << residualNorm2_ << std::endl;
}   