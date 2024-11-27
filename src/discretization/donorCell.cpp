#include "donorCell.h"
#include <cmath>
#include <iostream>

DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double,2> meshWidth, std::shared_ptr<Partitioning> partitioning_, double alpha) : 
    Discretization(nCells, meshWidth, partitioning_), alpha_(alpha) {}

double DonorCell::computeDu2Dx(int i, int j) const {
    const double u_left_sum = (u(i,j) + u(i-1,j)) / 2.0;
    const double u_left_diff = (u(i-1,j) - u(i,j)) / 2.0;
    const double u_right_sum = (u(i,j) + u(i+1,j)) / 2.0;
    const double u_right_diff = (u(i,j) - u(i+1,j)) / 2.0;

    const double cd_term = ((u_right_sum*u_right_sum) - (u_left_sum*u_left_sum)) / dx();
    const double donor_term = (std::fabs(u_right_sum) * u_right_diff - std::fabs(u_left_sum) * u_left_diff) / dx();
    
    return cd_term + alpha_ * donor_term;
}

double DonorCell::computeDv2Dy(int i, int j) const {
    const double v_bottom_sum = (v(i,j) + v(i,j-1)) / 2.0;
    const double v_bottom_diff = (v(i,j-1) - v(i,j)) / 2.0;
    const double v_top_sum = (v(i,j) + v(i,j+1)) / 2.0;
    const double v_top_diff = (v(i,j) - v(i,j+1)) / 2.0;

    const double cd_term = ((v_top_sum*v_top_sum) - (v_bottom_sum*v_bottom_sum)) / dy();
    const double donor_term = (std::fabs(v_top_sum) * v_top_diff - std::fabs(v_bottom_sum) * v_bottom_diff) / dy();

    return cd_term + alpha_ * donor_term;
}

double DonorCell::computeDuvDx(int i, int j) const {
    // CONVENTION: DIRECTION_POSITION
    const double u_top_sum = (u(i,j) + u(i,j+1)) / 2.0;
    const double v_right_sum = (v(i,j) + v(i+1,j)) / 2.0;
    const double u_top_left_sum = (u(i-1,j) + u(i-1,j+1)) / 2.0;
    const double v_left_sum = (v(i-1,j) + v(i,j)) / 2.0;
    const double v_right_diff = (v(i,j) - v(i+1,j)) / 2.0;
    const double v_left_diff = (v(i-1,j) - v(i,j)) / 2.0;

    const double cd_term = (u_top_sum*v_right_sum - u_top_left_sum*v_left_sum) / dx();
    const double donor_term = (std::fabs(u_top_sum) * v_right_diff - std::fabs(u_top_left_sum) * v_left_diff) / dx();

    return cd_term + alpha_ * donor_term;
}

double DonorCell::computeDuvDy(int i, int j) const {
    const double v_right_sum = (v(i,j) + v(i+1,j)) / 2.0;
    const double u_top_sum = (u(i,j) + u(i,j+1)) / 2.0;
    const double v_right_bottom_sum = (v(i,j-1) + v(i+1,j-1)) / 2.0;
    const double u_bottom_sum = (u(i,j-1) + u(i,j)) / 2.0;
    const double u_top_diff = (u(i,j) - u(i,j+1)) / 2.0;
    const double u_bottom_diff = (u(i,j-1) - u(i,j)) / 2.0;

    const double cd_term = (v_right_sum*u_top_sum - v_right_bottom_sum*u_bottom_sum) / dy();
    const double donor_term = (std::fabs(v_right_sum) * u_top_diff - std::fabs(v_right_bottom_sum) * u_bottom_diff) / dy();

    return cd_term + alpha_ * donor_term;
}