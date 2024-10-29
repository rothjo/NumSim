#include "FieldVariable.h"

FieldVariable::FieldVariable(std::array<int,2> size, std::array<int,2> origin, std::array<int,2> meshWidth) :
    Array2D(size), origin_(origin), meshWidth_(meshWidth) {}

double FieldVariable::interpolateAt(double x, double y) const {
    // eval f(x) with given origin, meshWidth and array2d
}