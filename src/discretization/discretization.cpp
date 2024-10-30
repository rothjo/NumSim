#include "discretization.h"

Discretization::Discretization(std::array<int, 2> nCells, std::array<double,2> meshWidth) : 
    StaggeredGrid(nCells, meshWidth) {};

/**
 * Compute the derivative of u^2 in x-direction
 * @param i index in x-direction
 * @param j index in y-direction
 * @return approximated derivative of u^2 in x-direction
 */
double Discretization::computeD2uDx2(int i, int j) const {
    return (u(i+1,j) - 2*u(i,j) + u(i-1,j)) / (dx() * dx());
}

/**
 * Compute the derivative of u^2 in y-direction
 * @param i index in x-direction
 * @param j index in y-direction
 * @return approximated derivative of u^2 in y-direction
 */
double Discretization::computeD2uDy2(int i, int j) const {
    return (u(i,j+1) - 2*u(i,j) + u(i,j-1)) / (dy() * dy());
}

/**
 * Compute the derivative of v^2 in x-direction
 * @param i index in x-direction
 * @param j index in y-direction
 * @return approximated derivative of v^2 in x-direction
 */
double Discretization::computeD2vDx2(int i, int j) const {
    return (v(i+1,j) - 2*v(i,j) + v(i-1,j)) / (dx() * dx());
}

/**
 * Compute the derivative of v^2 in y-direction
 * @param i index in x-direction
 * @param j index in y-direction
 * @return approximated derivative of v^2 in y-direction
 */
double Discretization::computeD2vDy2(int i, int j) const {
    return (v(i,j+1) - 2*v(i,j) + v(i,j-1)) / (dy() * dy());
}

/**
 * Compute the derivative of u*v in x-direction
 * @param i index in x-direction
 * @param j index in y-direction
 * @return approximated derivative of u*v in x-direction
 */
double Discretization::computeDpDx(int i, int j) const {
    return (p(i+1,j) - p(i,j)) / dx();
}

/**
 * Compute the derivative of u*v in y-direction
 * @param i index in x-direction
 * @param j index in y-direction
 * @return approximated derivative of u*v in y-direction
 */
double Discretization::computeDpDy(int i, int j) const {
    return (p(i,j+1) - p(i,j)) / dy();
}



