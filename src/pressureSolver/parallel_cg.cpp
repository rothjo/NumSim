#include "parallel_cg.h"
#include <numeric>
#include <iostream>

ParallelCG::ParallelCG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning)
    : ParallelPressureSolver(discretization, epsilon, maximumNumberOfIterations, partitioning),
    dx2_(discretization->dx() * discretization->dx()),
    dy2_(discretization->dy() * discretization->dy()),
    r_(FieldVariable({partitioning->nCellsLocal()[0] + 3, partitioning->nCellsLocal()[1] + 3}, {-1.5 * discretization_->dx(), -1.5 * discretization_->dy()}, discretization_->meshWidth())),
    d_(FieldVariable({partitioning->nCellsLocal()[0] + 3, partitioning->nCellsLocal()[1] + 3}, {-1.5 * discretization_->dx(), -1.5 * discretization_->dy()}, discretization_->meshWidth())),
    Ad_(FieldVariable({partitioning->nCellsLocal()[0] + 3, partitioning->nCellsLocal()[1] + 3}, {-1.5 * discretization_->dx(), -1.5 * discretization_->dy()}, discretization_->meshWidth()))
    {}


void ParallelCG::solve() {
    // init res, d
    // while loop
    // compute alpha: res_old / (d^T * Ad)
    // using Ad = d_xx + d_yy
    // communicate alpha (only alpha or nominator and demon separately?)
    // update p += alpha * d
    // update res_new = res_old - alpha * Ad
    // calculate res_new^T * res_new by summing up all res_new^2
    // global Sum res_new
    // if res_new^T * res_new < eps2, break
    // compute beta = res_new / res_old
    // d  = res_new + beta * d
    // res_old = res_new
    // communicate d boundaries

    // Initialize constants
    const double eps2 = epsilon_ * epsilon_;
    const int N = partitioning_->nCellsGlobal()[0] * partitioning_->nCellsGlobal()[1]; // Number of points used to compute norm
    const double N_eps2 = N * eps2;

    int pIBegin = discretization_->pIBegin();
    int pIEnd = discretization_->pIEnd();
    int pJBegin = discretization_->pJBegin();
    int pJEnd = discretization_->pJEnd();

    // Set initial data
    res_old2_ = 0.0;
    for(int i = pIBegin; i < pIEnd; i++) {
        for(int j = pJBegin; j < pJEnd; j++) {
            double res_ij = discretization_->rhs(i, j) - LaplaceP(i, j);
            r_(i,j) = res_ij;
            res_old2_ += res_ij * res_ij;
            d_(i, j) = res_ij;
        }
    }

    

    res_old2_ = partitioning_->globalSum(res_old2_);
    


    if (res_old2_ < N_eps2) {
        return;
    }
    // res_new2_ = res_old2_;
    numberOfIterations_ = maximumNumberOfIterations_;
    double dAd_ = 0.0;
    communicateAndBoundariesD();

    for (int k = 0; k < maximumNumberOfIterations_; ++k) {
        dAd_ = 0.0;
        // double laplace_d = 0.0;
        // Compute Ad = A * d
        for(int i = pIBegin; i < pIEnd; i++) {
            for(int j = pJBegin; j < pJEnd; j++) {
                const double laplace_d = LaplaceD(i,j);
                Ad_(i,j) = laplace_d;
                dAd_ += d_(i, j) * laplace_d;
            }
        }
        // std::cout << "before dAd" << dAd_ << std::endl;
        dAd_ = partitioning_->globalSum(dAd_);
        // std::cout << "after dAd" << dAd_ << std::endl;

        const double alpha = res_old2_ / dAd_;

        res_new2_ = 0.0;

        // Update pressure field and residual
        for(int i = pIBegin; i < pIEnd; i++) {
            for(int j = pJBegin; j < pJEnd; j++) {
                discretization_->p(i, j) += alpha * d_(i, j);
                r_(i, j) -= alpha * Ad_(i, j);

                res_new2_ += r_(i, j) * r_(i, j);
            }
        }
        res_new2_ = partitioning_->globalSum(res_new2_);
    	
        if (res_new2_ < N_eps2) {
            numberOfIterations_ = k;
            break;
        }

        double beta = res_new2_ / res_old2_;

        for(int i = pIBegin; i < pIEnd; i++) {
            for(int j = pJBegin; j < pJEnd; j++) {
                d_(i, j) = r_(i,j) + beta * d_(i, j);
            }
        }
        res_old2_ = res_new2_;

        // Synchronize boundary values
        communicateAndBoundariesD();
    }

    communicateAndBoundaries();

}

void ParallelCG::communicateAndBoundariesD() {
    // Initialize buffers and loop limits
    std::vector<double> topBuffer(discretization_->pIEnd() - discretization_->pIBegin(), 0.0);
    std::vector<double> bottomBuffer(discretization_->pIEnd() - discretization_->pIBegin(), 0.0);
    std::vector<double> leftBuffer(discretization_->pJEnd() - discretization_->pJBegin(), 0.0);
    std::vector<double> rightBuffer(discretization_->pJEnd() - discretization_->pJBegin(), 0.0);

    MPI_Request requestTop, requestBottom, requestLeft, requestRight;

    int pIBegin = discretization_->pIBegin();
    int pIEnd = discretization_->pIEnd();
    int pJBegin = discretization_->pJBegin();
    int pJEnd = discretization_->pJEnd();
    
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        // Set boundary values
        for (int i = pIBegin; i < pIEnd; i++) {
            d_(i, pJEnd) = d_(i, pJEnd - 1); 
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int i = pIBegin; i < pIEnd; i++) {
            topBuffer[i - pIBegin] = d_(i, pJEnd - 1);
        }
        // partitioning_->send(topBuffer, partitioning_->topNeighbourRankNo(), requestTop);
        // partitioning_->receive(topBuffer, partitioning_->topNeighbourRankNo(), requestTop);
        MPI_Isend(topBuffer.data(), topBuffer.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestTop);
        MPI_Irecv(topBuffer.data(), topBuffer.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestTop);
        // partitioning_->communicate(sendTopBuffer, receiveTopBuffer, partitioning_->topNeighbourRankNo(), requestSendTop, requestReceiveTop);

    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = pIBegin; i < pIEnd; i++) {
            d_(i, pJBegin - 1) = d_(i, pJBegin);
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int i = pIBegin; i < pIEnd; i++) {
            bottomBuffer[i - pIBegin] = d_(i, pJBegin);
        }
        // partitioning_->communicate(sendBottomBuffer, receiveBottomBuffer, partitioning_->bottomNeighbourRankNo(), requestSendBottom, requestReceiveBottom); 
        // partitioning_->send(bottomBuffer, partitioning_->bottomNeighbourRankNo(), requestBottom);
        // partitioning_->receive(bottomBuffer, partitioning_->bottomNeighbourRankNo(), requestBottom);  
        MPI_Isend(bottomBuffer.data(), bottomBuffer.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestBottom);
        MPI_Irecv(bottomBuffer.data(), bottomBuffer.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestBottom);
    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = pJBegin - 1; j < pJEnd + 1; j++) {
            d_(pIBegin - 1, j) = d_(pIBegin, j);
        } 
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int j = pJBegin; j < pJEnd; j++) {
            leftBuffer[j - pJBegin] = d_(pIBegin, j);
        }
        // partitioning_->communicate(sendLeftBuffer, receiveLeftBuffer, partitioning_->leftNeighbourRankNo(), requestSendLeft, requestReceiveLeft);
        // partitioning_->send(leftBuffer, partitioning_->leftNeighbourRankNo(), requestLeft);
        // partitioning_->receive(leftBuffer, partitioning_->leftNeighbourRankNo(), requestLeft); 
        MPI_Isend(leftBuffer.data(), leftBuffer.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestLeft);
        MPI_Irecv(leftBuffer.data(), leftBuffer.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestLeft);
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = pJBegin - 1; j < pJEnd + 1; j++) {
            d_(pIEnd, j) = d_(pIEnd - 1, j);
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int j = pJBegin; j < pJEnd; j++) {
            rightBuffer[j - pJBegin] = d_(pIEnd - 1, j);
        }
        // partitioning_->communicate(sendRightBuffer, receiveRightBuffer, partitioning_->rightNeighbourRankNo(), requestSendRight, requestReceiveRight);
        // partitioning_->send(rightBuffer, partitioning_->rightNeighbourRankNo(), requestRight);
        // partitioning_->receive(rightBuffer, partitioning_->rightNeighbourRankNo(), requestRight);
        MPI_Isend(rightBuffer.data(), rightBuffer.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestRight);
        MPI_Irecv(rightBuffer.data(), rightBuffer.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestRight);
    }

    // Set the buffers to the suiting columns/rows
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestTop, MPI_STATUS_IGNORE);
        for (int i = pIBegin; i < pIEnd; i++) {
            d_(i, pJEnd) = topBuffer[i - pIBegin];
        }   
    }

    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&requestBottom, MPI_STATUS_IGNORE);
        for (int i = pIBegin; i < pIEnd; i++) {
            d_(i, pJBegin - 1) = bottomBuffer[i - pIBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestLeft, MPI_STATUS_IGNORE);
        for (int j = pJBegin; j < pJEnd; j++) {
            d_(pIBegin - 1, j) = leftBuffer[j - pJBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestRight, MPI_STATUS_IGNORE);
        for (int j = pJBegin; j < pJEnd; j++) {
            d_(pIEnd, j) = rightBuffer[j - pJBegin];
        }
    }
}

double ParallelCG::LaplaceP(int i, int j) const {
    return ((discretization_->p(i + 1, j) - 2.0 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / dx2_) + ((discretization_->p(i, j + 1) - 2.0 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / dy2_);
}

double ParallelCG::LaplaceD(int i, int j) const {
    return ((d_(i + 1, j) - 2.0 * d_(i, j) + d_(i - 1, j)) / dx2_) + ((d_(i, j + 1) - 2.0 * d_(i, j) + d_(i, j - 1)) / dy2_);
}