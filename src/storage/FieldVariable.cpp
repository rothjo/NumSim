#include "FieldVariable.h"
#include <cassert>
#include <cmath>
#include <iostream>

FieldVariable::FieldVariable(std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth) :
    Array2D(size), origin_(origin), meshWidth_(meshWidth) {}

FieldVariable& FieldVariable::operator=(const FieldVariable& other) {
    if (this != &other) {
        Array2D::operator=(other); // Copy base class part
        origin_[0] = other.origin_[0];
        origin_[1] = other.origin_[1];
        meshWidth_[0] = other.meshWidth_[0];
        meshWidth_[1] = other.meshWidth_[1];
    }
    return *this;
}

// works only for constant dx, dy; rectangular grid
double FieldVariable::interpolateAt(double x, double y) const {
    // Check if x and y are located in the domain
    assert(0.0 <= x && x < size_[0]*meshWidth_[0]);
    assert(0.0 <= y && y < size_[1]*meshWidth_[1]);

    // Find the indicies i and j of the corresponding cells, we are looking for the bottom left point of the cell
    const int i = static_cast<int>((x - origin_[0]) / meshWidth_[0]);
    const int j = static_cast<int>((y - origin_[1]) / meshWidth_[1]);

    const double xLeft = origin_[0] + i * meshWidth_[0];
    const double yBottom = origin_[1] + j * meshWidth_[1];

    const double valueLeftBottom = (*this)(i,j);
    const double valueRightBottom = (*this)(i + 1,j);
    const double valueLeftTop = (*this)(i, j + 1);
    const double valueRightTop = (*this)(i + 1, j + 1);

    const double interpolateBottom = valueLeftBottom + (valueRightBottom - valueLeftBottom) * (x - xLeft) / meshWidth_[0];
    const double interpolateTop = valueLeftTop + (valueRightTop - valueLeftTop) * (x - xLeft) / meshWidth_[0];
    
    return interpolateBottom + (interpolateTop - interpolateBottom) * (y - yBottom) / meshWidth_[1];
}

double FieldVariable::computeMaxAbs() const {
    double max_value = 0.0;

    // Go through all points
    for (int i = 0; i < size_[0]; i++) {
        for (int j = 0; j < size_[1]; j++) {
            double abs_val = fabs((*this)(i,j));
            if (abs_val > max_value) {
                max_value = abs_val;
            }
        }     
    }

    return max_value;
}