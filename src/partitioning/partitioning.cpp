#include "partitioning.h"

//! TODO: implement correct decomposition method
void Partitioning::initialize(std::array<int,2> nCellsGlobal) {
    nSubdomains_ = {0, 0};
    nCellsGlobal_ = nCellsGlobal; 

    // int ranks = 2;

    // nCellsLocal_[0] = nCellsGlobal[0] / ranks;
    // nCellsLocal_[1] = nCellsGlobal[1];

    MPI_Comm_size(MPI_COMM_WORLD, &nRanks_);
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo_);

    MPI_Dims_create(nRanks_, 2, nSubdomains_.data());

    nCellsLocal_[0] = nCellsGlobal[0] / nSubdomains_[0];
    nCellsLocal_[1] = nCellsGlobal[1] / nSubdomains_[1];

    const int currentColumn = ownRankNo_ % nSubdomains_[0];
    const int currentRow = ownRankNo_ / nSubdomains_[0];
    
    nodeOffset_[0] = currentColumn * nCellsLocal_[0];
    nodeOffset_[1] = currentRow * nCellsLocal_[1];

    const int xRemainder = nCellsGlobal[0] % nSubdomains_[0];
    const int yRemainder = nCellsGlobal[1] % nSubdomains_[1];

    // split n remaining cells among the first n subdomains
    if (currentColumn < xRemainder) {
        nCellsLocal_[0] += 1;
        
        // increase offset with the number of cells in the previous subdomains
        nodeOffset_[0] += currentColumn;
    } else {
        // increase offset with the total number of remaining cells
        nodeOffset_[0] += xRemainder;
    }

    // analogous for y direction
    if (currentRow < yRemainder) {
        nCellsLocal_[1] += 1;
        nodeOffset_[1] += currentRow;
    } else {
        nodeOffset_[1] += yRemainder;
    }

    nodeOffsetSum_ = nodeOffset_[0] + nodeOffset_[1];

    // std::cout << "Rank: " << ownRankNo_ << std::endl;
    // std::cout << " nCellsLocal: " << nCellsLocal_[0] << " " << nCellsLocal_[1] << std::endl;
    // std::cout << "nodeOffset: " << nodeOffset_[0] << " " << nodeOffset_[1] << '\n' << std::endl;
}

std::array<int,2> Partitioning::nCellsLocal() const {
    return nCellsLocal_;
}

std::array<int,2> Partitioning::nCellsGlobal() const {
    return nCellsGlobal_;
}

int Partitioning::ownRankNo() const {
    return ownRankNo_;
}

int Partitioning::nRanks() const {
    return nRanks_;
}

bool Partitioning::ownPartitionContainsBottomBoundary() const {
    return (ownRankNo_ < nSubdomains_[0]);
}

bool Partitioning::ownPartitionContainsTopBoundary() const {
    return (ownRankNo_ >= nRanks_ - nSubdomains_[0]);
}

bool Partitioning::ownPartitionContainsLeftBoundary() const {
    return (ownRankNo_ % nSubdomains_[0] == 0);
}

bool Partitioning::ownPartitionContainsRightBoundary() const {
    return (ownRankNo_ % nSubdomains_[0] == nSubdomains_[0] - 1);
}

int Partitioning::leftNeighbourRankNo() const {
    // assert(ownRankNo_ % nSubdomains_[0] != 0);
    assert(!ownPartitionContainsLeftBoundary());
    return ownRankNo_ - 1;
}

int Partitioning::rightNeighbourRankNo() const {
    // assert(ownRankNo_ % nSubdomains_[0] != nSubdomains_[0] - 1);
    assert(!ownPartitionContainsRightBoundary());
    return ownRankNo_ + 1;
}

int Partitioning::topNeighbourRankNo() const {
    // assert(ownRankNo_ < nRanks_ - nSubdomains_[0]);
    assert(!ownPartitionContainsTopBoundary());
    return ownRankNo_ + nSubdomains_[0];
}

int Partitioning::bottomNeighbourRankNo() const {
    // assert(ownRankNo_ >= nSubdomains_[0]);
    assert(!ownPartitionContainsBottomBoundary());
    return ownRankNo_ - nSubdomains_[0];
}

std::array<int,2> Partitioning::nodeOffset() const {
    return nodeOffset_;
}

int Partitioning::nodeOffsetSum() const {
    return nodeOffsetSum_;
}

double Partitioning::globalSum(double localValue) const {
    double globalValue = 0.0;
    MPI_Allreduce(&localValue, &globalValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return globalValue;
}


double Partitioning::globalMin(double localValue) const {
    double globalValue = 0.0;
    MPI_Allreduce(&localValue, &globalValue, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return globalValue;
}


void Partitioning::send(std::vector<double> valuesToSend, int neighbourRankNo, MPI_Request &requestSend) {
    const int nValuesToSend = valuesToSend.size();
    MPI_Isend(valuesToSend.data(), nValuesToSend, MPI_DOUBLE, neighbourRankNo, 0, MPI_COMM_WORLD, &requestSend);
}

void Partitioning::receive(std::vector<double> &valuesToReceive, int neighbourRankNo, MPI_Request &requestReceive) {
    const int nValuesToReceive = valuesToReceive.size();
    MPI_Irecv(valuesToReceive.data(), nValuesToReceive, MPI_DOUBLE, neighbourRankNo, 0, MPI_COMM_WORLD, &requestReceive);
}

// MPI_wait is still required after this function
void Partitioning::communicate(std::vector<double> valuesToSend, std::vector<double> &valuesToReceive, int neighbourRankNo, MPI_Request &requestSend, MPI_Request &requestReceive) {
    this->send(valuesToSend, neighbourRankNo, requestSend);
    this->receive(valuesToReceive, neighbourRankNo, requestReceive);
}