#pragma once

#include <array>
#include <mpi.h>
#include <cassert>
#include <iostream>
#include <vector>

class Partitioning
{
public:

  //! compute partitioning, set internal variables
  void initialize(std::array<int,2> nCellsGlobal);

  //! get the local number of cells in the own subdomain
  std::array<int,2> nCellsLocal() const;

  //! get the global number of cells in the whole computational domain
  //! used in OutputWriterParaviewParallel
  std::array<int,2> nCellsGlobal() const;

  //! get the own MPI rank no
  //! used in OutputWriterParaviewParallel and OutputWriterTextParallel
  int ownRankNo() const;

  //! number of MPI ranks
  int nRanks() const;

  //! if the own partition has part of the bottom boundary of the whole domain
  bool ownPartitionContainsBottomBoundary() const;

  //! if the own partition has part of the top boundary of the whole domain
  //! used in OutputWriterParaviewParallel
  bool ownPartitionContainsTopBoundary() const;

  //! if the own partition has part of the left boundary of the whole domain
  bool ownPartitionContainsLeftBoundary() const;

  //! if the own partition has part of the right boundary of the whole domain
  //! used in OutputWriterParaviewParallel
  bool ownPartitionContainsRightBoundary() const;

  //! get the rank no of the left neighbouring rank
  int leftNeighbourRankNo() const;

  //! get the rank no of the right neighbouring rank
  int rightNeighbourRankNo() const;

  //! get the rank no of the top neighbouring rank
  int topNeighbourRankNo() const;

  //! get the rank no of the bottom neighbouring rank
  int bottomNeighbourRankNo() const;

  //! get the offset values for counting local nodes in x and y direction. 
  //! (i_local,j_local) + nodeOffset = (i_global,j_global)
  //! used in OutputWriterParaviewParallel
  std::array<int,2> nodeOffset() const;

  //! get the sum of the offset in x and y direction
  int nodeOffsetSum() const;

  //! get the global sum over all ranks of a local value
  double globalSum(double localValue) const;

  //! 
  void communicate(std::vector<double> valuesToSend, std::vector<double> &valuesToReceive, int neighbourRankNo, MPI_Request &requestSend, MPI_Request &requestReceive);


  //! TODO: sendToTop, Bottom, Left, Right

private:
  std::array<int,2> nCellsGlobal_;
  std::array<int,2> nCellsLocal_;
  int ownRankNo_;
  int nRanks_;
  int nodeOffsetSum_;
  std::array<int,2> nSubdomains_;
  std::array<int,2> nodeOffset_;
};
