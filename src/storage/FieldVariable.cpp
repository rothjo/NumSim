#include "FieldVariable.h"
#include <cassert>
#include <cmath>
#include <iostream>

FieldVariable::FieldVariable(std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth) :
    Array2D(size), origin_(origin), meshWidth_(meshWidth) {}


/**
 * Get the bilinear interpolated value of the discrete field variable
 * @param x x-coordinate of the point
 * @param y y-coordinate of the point
 * @return interpolated value of the discrete field variable at position (x,y)
 */
double FieldVariable::interpolateAt(double x, double y) const {
    
    // Check if x and y are located in the domain
    assert(0.0 <= x && x < size_[0]*meshWidth_[0]);
    assert(0.0 <= y && y < size_[1]*meshWidth_[1]);

    // Find the indicies i and j of the corresponding cells, we are looking for the bottom left point of the cell
    int i = (x - origin_[0]) / meshWidth_[0];
    int j = (y - origin_[1]) / meshWidth_[1];

    // Find the values of the fieldvariable for the interpolation at the neighbouring points
    double valueLeftBottom = (*this)(i,j);
    double valueRightBottom = (*this)(i + 1,j);
    double valueLeftTop = (*this)(i, j + 1);
    double valueRightTop = (*this)(i + 1, j + 1);

    // Find the x and y coordinates of the neighbouring points

    double xLeft = origin_[0] + i * meshWidth_[0];
    double xRight = xLeft + meshWidth_[0];
    double yBottom = origin_[1] + j * meshWidth_[1];
    double yTop = yBottom + meshWidth_[1];


    // Use bilinear interpolation to find the interpolated value

    // First, interpolate in x-direction for the bottom horizontal neighboring line 
    double interpolateBottom = (xRight - x) / (xRight - xLeft) * valueLeftBottom + (x - xLeft) / (xRight - xLeft) * valueRightBottom;
    
    // Interpolate in x-direction for the top horizontal neighboring line 
    double interpolateTop = (xRight - x) / (xRight - xLeft) * valueLeftTop + (x - xLeft) / (xRight - xLeft) * valueRightTop;

    // Interpolate in y-direction
    double interpolation = (yTop - y) / (yTop - yBottom) * interpolateBottom + (y - yBottom) / (yTop - yBottom) * interpolateTop;

    return interpolation;

}

double FieldVariable::computeMaxAbs() const {
    double max_value = 0.0;

    // Go through inner points
    for (int i = 1; i < (size_[0] - 1); i++) {
        for (int j = 1; j < (size_[1] - 1); j++) {
            max_value = std::max(max_value, std::fabs((*this)(i,j)));
        }     
    }

    return max_value;
}
    









