#include "parallel_cg.h"
#include <numeric>
#include <iostream>

ParallelCG::ParallelCG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning)
    : ParallelPressureSolver(discretization, epsilon, maximumNumberOfIterations, partitioning),
    dx2_((*discretization_).dx() * (*discretization_).dx()),
    dy2_(discretization->dy() * discretization->dy()),
    r_(FieldVariable({(*partitioning_).nCellsLocal()[0] + 3, (*partitioning_).nCellsLocal()[1] + 3}, {-1.5 * discretization_->dx(), -1.5 * discretization_->dy()}, discretization_->meshWidth())),
    d_(FieldVariable({(*partitioning_).nCellsLocal()[0] + 3, (*partitioning_).nCellsLocal()[1] + 3}, {-1.5 * discretization_->dx(), -1.5 * discretization_->dy()}, discretization_->meshWidth())),
    Ad_(FieldVariable({(*partitioning_).nCellsLocal()[0] + 3, (*partitioning_).nCellsLocal()[1] + 3}, {-1.5 * discretization_->dx(), -1.5 * discretization_->dy()}, discretization_->meshWidth()))
    {}


void ParallelCG::solve() {
    // Initialize constants and loop limits
    const double eps2 = epsilon_ * epsilon_;
    const int N = (*partitioning_).nCellsGlobal()[0] * (*partitioning_).nCellsGlobal()[1];
    const double N_eps2 = N * eps2;

    int pIBegin = discretization_->pIBegin();
    int pIEnd = discretization_->pIEnd();
    int pJBegin = discretization_->pJBegin();
    int pJEnd = discretization_->pJEnd();

    // Initialize first time step
    res_old2_ = 0.0;
    const double M_inv = 1.0 / (2.0 / dx2_ + 2.0 / dy2_); // Diagonal element approximation, Jacobi
    // Compute initial residual r and apply preconditioner, i.e. solve Mz_0 = r_0
    for (int i = pIBegin; i < pIEnd; i++) {
        for (int j = pJBegin; j < pJEnd; j++) {
            double res_ij = discretization_->rhs(i, j) - LaplaceP(i, j);
            r_(i, j) = res_ij;

            // Apply preconditioner (Jacobi: divide by diagonal element)
            d_(i, j) = res_ij * M_inv;

            res_old2_ += r_(i, j) * d_(i, j); // Compute preconditioned residual
        }
    }

    // Compute residual
    MPI_Request req_res_old2;
    MPI_Iallreduce(MPI_IN_PLACE, &res_old2_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req_res_old2);

    communicateAndBoundariesD();

    MPI_Wait(&req_res_old2, MPI_STATUS_IGNORE);

    if (res_old2_ < N_eps2) {
        return;
    }

    numberOfIterations_ = maximumNumberOfIterations_;
    double dAd_ = 0.0;

    // Time loop
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
        double alpha;
        MPI_Wait(&req_dAd, MPI_STATUS_IGNORE);

        // Compute alpha
        alpha = res_old2_ / dAd_;

        res_new2_ = 0.0;

        // Update p, r and compute local res_new2_
        for (int i = pIBegin; i < pIEnd; i++) {
            for (int j = pJBegin; j < pJEnd; j++) {
                (*discretization_).p(i, j) += alpha * d_(i, j);
                r_(i, j) -= alpha * Ad_(i, j);
                // Apply preconditioner to update residual
                double z_ij = r_(i, j) * M_inv;
                res_new2_ += r_(i, j) * z_ij;
                d_(i, j) = z_ij + (res_new2_ / res_old2_) * d_(i, j);
            }
        }

        MPI_Request req_res_new2;
        MPI_Iallreduce(MPI_IN_PLACE, &res_new2_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req_res_new2);
        MPI_Wait(&req_res_new2, MPI_STATUS_IGNORE); 

        // Check if new residuum is lower than tolerance
        if (res_new2_ < N_eps2) {
            numberOfIterations_ = k;
            break;
        }

        res_old2_ = res_new2_;
        communicateAndBoundariesD();
    }

    communicateAndBoundaries();
}

void ParallelCG::communicateAndBoundariesD() {

    // Initialize buffers and loop limits
    int pIBegin = (*discretization_).pIBegin();
    int pIEnd = (*discretization_).pIEnd();
    int pJBegin = (*discretization_).pJBegin();
    int pJEnd = (*discretization_).pJEnd();

    std::vector<double> topBuffer(pIEnd - pIBegin, 0.0);
    std::vector<double> bottomBuffer(pIEnd - pIBegin, 0.0);
    std::vector<double> leftBuffer(pJEnd - pJBegin, 0.0);
    std::vector<double> rightBuffer(pJEnd - pJBegin, 0.0);

    // Requests
    MPI_Request requestTop, requestBottom, requestLeft, requestRight;

    // Go through all boundaries and send/recv all needed values
    if ((*partitioning_).ownPartitionContainsTopBoundary()) {
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
        MPI_Isend(topBuffer.data(), topBuffer.size(), MPI_DOUBLE, (*partitioning_).topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestTop);
        MPI_Irecv(topBuffer.data(), topBuffer.size(), MPI_DOUBLE, (*partitioning_).topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestTop);
    }

    if ((*partitioning_).ownPartitionContainsBottomBoundary()) {
        for (int i = pIBegin; i < pIEnd; i++) {
            d_(i, pJBegin - 1) = d_(i, pJBegin);
        }
    } else {
        for (int i = pIBegin; i < pIEnd; i++) {
            bottomBuffer[i - pIBegin] = d_(i, pJBegin);
        } 
        MPI_Isend(bottomBuffer.data(), bottomBuffer.size(), MPI_DOUBLE, (*partitioning_).bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestBottom);
        MPI_Irecv(bottomBuffer.data(), bottomBuffer.size(), MPI_DOUBLE, (*partitioning_).bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestBottom);
    }

    if ((*partitioning_).ownPartitionContainsLeftBoundary()) {
        for (int j = pJBegin - 1; j < pJEnd + 1; j++) {
            d_(pIBegin - 1, j) = d_(pIBegin, j);
        } 
    } else {
        for (int j = pJBegin; j < pJEnd; j++) {
            leftBuffer[j - pJBegin] = d_(pIBegin, j);
        }
        MPI_Isend(leftBuffer.data(), leftBuffer.size(), MPI_DOUBLE, (*partitioning_).leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestLeft);
        MPI_Irecv(leftBuffer.data(), leftBuffer.size(), MPI_DOUBLE, (*partitioning_).leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestLeft);
    }

    if ((*partitioning_).ownPartitionContainsRightBoundary()) {
        for (int j = pJBegin - 1; j < pJEnd + 1; j++) {
            d_(pIEnd, j) = d_(pIEnd - 1, j);
        }
    } else {
        for (int j = pJBegin; j < pJEnd; j++) {
            rightBuffer[j - pJBegin] = d_(pIEnd - 1, j);
        }
        MPI_Isend(rightBuffer.data(), rightBuffer.size(), MPI_DOUBLE, (*partitioning_).rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestRight);
        MPI_Irecv(rightBuffer.data(), rightBuffer.size(), MPI_DOUBLE, (*partitioning_).rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestRight);
    }

    // Set the buffers to the suiting columns/rows, go through all boundaries
    if (!(*partitioning_).ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestTop, MPI_STATUS_IGNORE);
        for (int i = pIBegin; i < pIEnd; i++) {
            d_(i, pJEnd) = topBuffer[i - pIBegin];
        }   
    }

    if (!(*partitioning_).ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&requestBottom, MPI_STATUS_IGNORE);
        for (int i = pIBegin; i < pIEnd; i++) {
            d_(i, pJBegin - 1) = bottomBuffer[i - pIBegin];
        }
    }

    if (!(*partitioning_).ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestLeft, MPI_STATUS_IGNORE);
        for (int j = pJBegin; j < pJEnd; j++) {
            d_(pIBegin - 1, j) = leftBuffer[j - pJBegin];
        }
    }

    if (!(*partitioning_).ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestRight, MPI_STATUS_IGNORE);
        for (int j = pJBegin; j < pJEnd; j++) {
            d_(pIEnd, j) = rightBuffer[j - pJBegin];
        }
    }

}

double ParallelCG::LaplaceP(int i, int j) const {
    return (((*discretization_).p(i + 1, j) - 2.0 * (*discretization_).p(i, j) + (*discretization_).p(i - 1, j)) / dx2_) + (((*discretization_).p(i, j + 1) - 2.0 * (*discretization_).p(i, j) + (*discretization_).p(i, j - 1)) / dy2_);
}

double ParallelCG::LaplaceD(int i, int j) const {
    return ((d_(i + 1, j) - 2.0 * d_(i, j) + d_(i - 1, j)) / dx2_) + ((d_(i, j + 1) - 2.0 * d_(i, j) + d_(i, j - 1)) / dy2_);
}