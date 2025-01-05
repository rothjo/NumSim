#include "centralDifferences.h"
#include <iostream>
#include <cmath>

CentralDifferences::CentralDifferences(std::array<int, 2> nCells, std::array<double, 2> meshWidth, std::shared_ptr<Partitioning> partitioning_) :
    Discretization(nCells, meshWidth, partitioning_) {}

double CentralDifferences::computeDu2Dx(int i, int j) const {
    const double u_left = (u(i, j) + u(i - 1, j)) / 2.0;
    const double u_right = (u(i, j) + u(i + 1, j)) / 2.0;
    return (u_right * u_right - u_left * u_left) / dx();
}

double CentralDifferences::computeDv2Dy(int i, int j) const {
    const double v_bottom = (v(i, j) + v(i, j - 1)) / 2.0;
    const double v_top = (v(i, j) + v(i, j + 1)) / 2.0;
    return (v_top * v_top - v_bottom * v_bottom) / dy();
}

double CentralDifferences::computeDuvDx(int i, int j) const {
    // CONVENTION: DIRECTION_POSITION
    const double v_left = (v(i, j) + v(i - 1, j)) / 2.0;
    const double v_right = (v(i, j) + v(i + 1, j)) / 2.0;

    const double u_top = (u(i, j) + u(i, j + 1)) / 2.0;
    const double u_top_left = (u(i - 1, j) + u(i - 1, j + 1)) / 2.0;

    return (u_top * v_right - u_top_left * v_left) / dx();
}

double CentralDifferences::computeDuvDy(int i, int j) const {
    const double u_bottom = (u(i, j) + u(i, j - 1)) / 2.0;
    const double u_top = (u(i, j) + u(i, j + 1)) / 2.0;

    const double v_right = (v(i, j) + v(i + 1, j)) / 2.0;
    const double v_right_down = (v(i, j - 1) + v(i + 1, j - 1)) / 2.0;

    return (u_top * v_right - u_bottom * v_right_down) / dy();
}

double CentralDifferences::computeDuTDx(int i, int j) const {
    const double t_right_sum = (t(i,j) + t(i+1,j)) / 2.0;
    const double t_left_sum = (t(i-1,j) + t(i,j)) / 2.0;
    const double cd_term = (u(i,j)*t_right_sum - u(i-1,j)*t_left_sum) / dx();
    return cd_term;
}

double CentralDifferences::computeDvTDy(int i, int j) const {
    const double t_top_sum = (t(i,j) + t(i,j+1)) / 2.0;
    const double t_bottom_sum = (t(i,j-1) + t(i,j)) / 2.0;
    const double cd_term = (v(i,j)*t_top_sum - v(i, j-1)*t_bottom_sum) / dy();
    return cd_term;
}
