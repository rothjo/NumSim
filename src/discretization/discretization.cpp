#include "discretization.h"

Discretization::Discretization(std::array<int, 2> nCells, std::array<double,2> meshWidth, std::shared_ptr<Partitioning> partitioning_) : 
    StaggeredGrid(nCells, meshWidth, partitioning_) {}

double Discretization::computeD2uDx2(int i, int j) const {
    return (u(i+1,j) - 2.0 * u(i,j) + u(i-1,j)) / (dx() * dx());
}

double Discretization::computeD2uDy2(int i, int j) const {
    return (u(i,j+1) - 2.0 * u(i,j) + u(i,j-1)) / (dy() * dy());
}

double Discretization::computeD2vDx2(int i, int j) const {
    return (v(i+1,j) - 2.0 * v(i,j) + v(i-1,j)) / (dx() * dx());
}

double Discretization::computeD2vDy2(int i, int j) const {
    return (v(i,j+1) - 2.0 * v(i,j) + v(i,j-1)) / (dy() * dy());
}

double Discretization::computeDpDx(int i, int j) const {
    return (p(i+1,j) - p(i,j)) / dx();
}

double Discretization::computeDpDy(int i, int j) const {
    return (p(i,j+1) - p(i,j)) / dy();
}