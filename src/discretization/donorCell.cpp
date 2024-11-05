#include "donorCell.h"
#include <cmath>

DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double,2> meshWidth, double alpha) : 
    Discretization(nCells, meshWidth), alpha_(alpha) {};

double DonorCell::computeDu2Dx(int i, int j) const {
    const double u_left_sum = (u(i,j) + u(i-1,j)) / 2.0;
    const double u_left_diff = (u(i,j) - u(i-1,j)) / 2.0;
    const double u_right_sum = (u(i,j) + u(i+1,j)) / 2.0;
    const double u_right_diff = (u(i,j) - u(i+1,j)) / 2.0;

    const double diff_of_squares = ((u_right_sum*u_right_sum) - (u_left_sum*u_left_sum)) / dx();
    const double alpha_term = alpha_ * (std::abs(u_right_sum) * u_right_diff - std::abs(u_left_sum) * u_left_diff) / dx();
    
    return diff_of_squares + alpha_term;
}

double DonorCell::computeDv2Dy(int i, int j) const {
    const double v_bottom_sum = (v(i,j) + v(i,j-1)) / 2.0;
    const double v_bottom_diff = (v(i,j) - v(i,j-1)) / 2.0;
    const double v_top_sum = (v(i,j) + v(i,j+1)) / 2.0;
    const double v_top_diff = (v(i,j) - v(i,j+1)) / 2.0;

    const double cd_term = ((v_top_sum*v_top_sum) - (v_bottom_sum*v_bottom_sum)) / dy();
    const double alpha_term = alpha_ * (std::abs(v_top_sum) * v_top_diff - std::abs(v_bottom_sum) * v_bottom_diff) / dy();

    return cd_term + alpha_term;
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
    const double alpha_term = alpha_ * (std::abs(u_top_sum) * v_right_diff - std::abs(u_top_left_sum) * v_left_diff) / dx();

    return cd_term + alpha_term;
}

double DonorCell::computeDuvDy(int i, int j) const {
    const double v_right_sum = (v(i,j) + v(i+1,j)) / 2.0;
    const double u_top_sum = (u(i,j) + u(i,j+1)) / 2.0;
    const double v_right_bottom_sum = (v(i,j-1) + v(i+1,j-1)) / 2.0;
    const double u_bottom_sum = (u(i,j-1) + u(i,j)) / 2.0;
    const double u_top_diff = (u(i,j) - u(i,j+1)) / 2.0;
    const double u_bottom_diff = (u(i,j-1) - u(i,j)) / 2.0;

    const double cd_term = (v_right_sum*u_top_sum - v_right_bottom_sum*u_bottom_sum) / dy();
    const double alpha_term = alpha_ * (std::abs(v_right_sum) * u_top_diff - std::abs(v_right_bottom_sum) * u_bottom_diff) / dy();

    return cd_term + alpha_term;
}