#include "centralDifferences.h"
#include <iostream>
#include <cmath>

CentralDifferences::CentralDifferences(std::array<int, 2> nCells, std::array<double,2> meshWidth) : 
    Discretization(nCells, meshWidth) {}

double CentralDifferences::computeDu2Dx(int i, int j) const {
    const double u_left = (u(i,j) + u(i-1,j)) / 2.0;
    const double u_right = (u(i,j) + u(i+1,j)) / 2.0;
    return (u_right*u_right - u_left*u_left) / dx();
}

double CentralDifferences::computeDv2Dy(int i, int j) const {
    const double v_bottom = (v(i,j) + v(i,j-1)) / 2.0;
    const double v_top = (v(i,j) + v(i,j+1)) / 2.0;
    return (v_top*v_top - v_bottom*v_bottom) / dy();
}

double CentralDifferences::computeDuvDx(int i, int j) const {
    // CONVENTION: DIRECTION_POSITION
    const double v_left = (v(i,j) + v(i-1,j)) / 2.0;
    const double v_right = (v(i,j) + v(i+1,j)) / 2.0;

    const double u_top = (u(i,j) + u(i,j+1)) / 2.0;
    const double u_top_left = (u(i-1,j) + u(i-1,j+1)) / 2.0;

    return (u_top*v_right - u_top_left*v_left) / dx();
}

double CentralDifferences::computeDuvDy(int i, int j) const {
    const double u_bottom = (u(i,j) + u(i,j-1)) / 2.0;
    const double u_top = (u(i,j) + u(i,j+1)) / 2.0;

    const double v_right = (v(i,j) + v(i+1,j)) / 2.0;
    const double v_right_down = (v(i,j-1) + v(i+1,j-1)) / 2.0;

    return (u_top*v_right - u_bottom*v_right_down) / dy();
}
