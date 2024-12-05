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
    const double eps2 = epsilon_ * epsilon_;
    const int N = partitioning_->nCellsGlobal()[0] * partitioning_->nCellsGlobal()[1];
    const double N_eps2 = N * eps2;

    int pIBegin = discretization_->pIBegin();
    int pIEnd = discretization_->pIEnd();
    int pJBegin = discretization_->pJBegin();
    int pJEnd = discretization_->pJEnd();

    res_old2_ = 0.0;
    for (int i = pIBegin; i < pIEnd; i++) {
        for (int j = pJBegin; j < pJEnd; j++) {
            double res_ij = discretization_->rhs(i, j) - LaplaceP(i, j);
            r_(i, j) = res_ij;
            res_old2_ += res_ij * res_ij;
            d_(i, j) = res_ij;
        }
    }

    // Non-blocking reduction for res_old2_
    MPI_Request req_res_old2;
    MPI_Iallreduce(MPI_IN_PLACE, &res_old2_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req_res_old2);

    communicateAndBoundariesD();

    MPI_Wait(&req_res_old2, MPI_STATUS_IGNORE); // Wait for res_old2_ reduction
    if (res_old2_ < N_eps2) {
        return;
    }

    numberOfIterations_ = maximumNumberOfIterations_;
    double dAd_ = 0.0;

    for (int k = 0; k < maximumNumberOfIterations_; ++k) {
        dAd_ = 0.0;

        // Compute Ad and partial dAd_ locally
        for (int i = pIBegin; i < pIEnd; i++) {
            for (int j = pJBegin; j < pJEnd; j++) {
                const double laplace_d = LaplaceD(i, j);
                Ad_(i, j) = laplace_d;
                dAd_ += d_(i, j) * laplace_d;
            }
        }

        // Non-blocking reduction for dAd_
        MPI_Request req_dAd;
        MPI_Iallreduce(MPI_IN_PLACE, &dAd_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req_dAd);

        // Perform computations while dAd_ reduction progresses
        double alpha;
        MPI_Wait(&req_dAd, MPI_STATUS_IGNORE); // Wait for dAd_ reduction
        alpha = res_old2_ / dAd_;

        res_new2_ = 0.0;

        // Update p, r and compute local res_new2_
        for (int i = pIBegin; i < pIEnd; i++) {
            for (int j = pJBegin; j < pJEnd; j++) {
                discretization_->p(i, j) += alpha * d_(i, j);
                r_(i, j) -= alpha * Ad_(i, j);
                res_new2_ += r_(i, j) * r_(i, j);
            }
        }

        // Non-blocking reduction for res_new2_
        MPI_Request req_res_new2;
        MPI_Iallreduce(MPI_IN_PLACE, &res_new2_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req_res_new2);

        MPI_Wait(&req_res_new2, MPI_STATUS_IGNORE); // Wait for res_new2_ reduction
        if (res_new2_ < N_eps2) {
            numberOfIterations_ = k;
            break;
        }

        double beta = res_new2_ / res_old2_;

        // Update search direction d
        for (int i = pIBegin; i < pIEnd; i++) {
            for (int j = pJBegin; j < pJEnd; j++) {
                d_(i, j) = r_(i, j) + beta * d_(i, j);
            }
        }

        res_old2_ = res_new2_;
        communicateAndBoundariesD(); // Synchronize boundary values
    }

    communicateAndBoundaries();
}

void ParallelCG::communicateAndBoundariesD() {
    //! TODO: Put into else statements afterwards, optimize sizes of buffers
    std::vector<double> topBuffer(discretization_->pIEnd() - discretization_->pIBegin(), 0.0);
    std::vector<double> bottomBuffer(discretization_->pIEnd() - discretization_->pIBegin(), 0.0);
    std::vector<double> leftBuffer(discretization_->pJEnd() - discretization_->pJBegin(), 0.0);
    std::vector<double> rightBuffer(discretization_->pJEnd() - discretization_->pJBegin(), 0.0);

    MPI_Request requestTop, requestBottom, requestLeft, requestRight;
    
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        // Set boundary values
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            d_(i, discretization_->pJEnd()) = d_(i, discretization_->pJEnd() - 1); 
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            topBuffer[i - discretization_->pIBegin()] = d_(i, discretization_->pJEnd() - 1);
        }
        partitioning_->send(topBuffer, partitioning_->topNeighbourRankNo(), requestTop);
        partitioning_->receive(topBuffer, partitioning_->topNeighbourRankNo(), requestTop);
        // partitioning_->communicate(sendTopBuffer, receiveTopBuffer, partitioning_->topNeighbourRankNo(), requestSendTop, requestReceiveTop);

    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            d_(i, discretization_->pJBegin() - 1) = d_(i, discretization_->pJBegin());
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            bottomBuffer[i - discretization_->pIBegin()] = d_(i, discretization_->pJBegin());
        }
        // partitioning_->communicate(sendBottomBuffer, receiveBottomBuffer, partitioning_->bottomNeighbourRankNo(), requestSendBottom, requestReceiveBottom); 
        partitioning_->send(bottomBuffer, partitioning_->bottomNeighbourRankNo(), requestBottom);
        partitioning_->receive(bottomBuffer, partitioning_->bottomNeighbourRankNo(), requestBottom);  
    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = discretization_->pJBegin() - 1; j < discretization_->pJEnd() + 1; j++) {
            d_(discretization_->pIBegin() - 1, j) = d_(discretization_->pIBegin(), j);
        } 
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            leftBuffer[j - discretization_->pJBegin()] = d_(discretization_->pIBegin(), j);
        }
        // partitioning_->communicate(sendLeftBuffer, receiveLeftBuffer, partitioning_->leftNeighbourRankNo(), requestSendLeft, requestReceiveLeft);
        partitioning_->send(leftBuffer, partitioning_->leftNeighbourRankNo(), requestLeft);
        partitioning_->receive(leftBuffer, partitioning_->leftNeighbourRankNo(), requestLeft); 
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = discretization_->pJBegin() - 1; j < discretization_->pJEnd() + 1; j++) {
            d_(discretization_->pIEnd(), j) = d_(discretization_->pIEnd() - 1, j);
        }
    } else {
        // Send top boundary values to top neighbour
        // Receive top boundary values from top neighbour
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            rightBuffer[j - discretization_->pJBegin()] = d_(discretization_->pIEnd() - 1, j);
        }
        // partitioning_->communicate(sendRightBuffer, receiveRightBuffer, partitioning_->rightNeighbourRankNo(), requestSendRight, requestReceiveRight);
        partitioning_->send(rightBuffer, partitioning_->rightNeighbourRankNo(), requestRight);
        partitioning_->receive(rightBuffer, partitioning_->rightNeighbourRankNo(), requestRight);
    }

    // Set the buffers to the suiting columns/rows
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestTop, MPI_STATUS_IGNORE);
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            d_(i, discretization_->pJEnd()) = topBuffer[i - discretization_->pIBegin()];
        }   
    }

    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&requestBottom, MPI_STATUS_IGNORE);
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            d_(i, discretization_->pJBegin() - 1) = bottomBuffer[i - discretization_->pIBegin()];
        }
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestLeft, MPI_STATUS_IGNORE);
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            d_(discretization_->pIBegin() - 1, j) = leftBuffer[j - discretization_->pJBegin()];
        }
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestRight, MPI_STATUS_IGNORE);
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            d_(discretization_->pIEnd(), j) = rightBuffer[j - discretization_->pJBegin()];
        }
    }

}

double ParallelCG::LaplaceP(int i, int j) const {
    return ((discretization_->p(i + 1, j) - 2.0 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / dx2_) + ((discretization_->p(i, j + 1) - 2.0 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / dy2_);
}

double ParallelCG::LaplaceD(int i, int j) const {
    return ((d_(i + 1, j) - 2.0 * d_(i, j) + d_(i - 1, j)) / dx2_) + ((d_(i, j + 1) - 2.0 * d_(i, j) + d_(i, j - 1)) / dy2_);
}