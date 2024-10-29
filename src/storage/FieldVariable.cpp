#include "FieldVariable.h"
#include <cassert>

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

    // Find the indicies i and j of the corresponding cells, we are looking for the buttom left point of the cell
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
    double yButtom = origin_[1] + j * meshWidth_[1];
    double yTop = yButtom + meshWidth_[1];


    // Use bilinear interpolation to find the interpolated value

    // First, interpolate in x-direction for the bottom horizontal neighboring line 
    double interpolateButtom = (xRight - x) / (xRight - xLeft) * valueLeftBottom + (x - xLeft) / (xRight - xLeft) * valueRightBottom;
    
    // Interpolate in x-direction for the top horizontal neighboring line 
    double interpolateTop = (xRight - x) / (xRight - xLeft) * valueLeftTop + (x - xLeft) / (xRight - xLeft) * valueRightTop;

    // Interpolate in y-direction
    double interpolation = (yTop - y) / (yTop - yButtom) * interpolateButtom + (y - yButtom) / (yTop - yButtom) * interpolateTop;

    return interpolation;

}

    









