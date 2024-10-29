#pragma once

#include "array2d.h"

/** A field variable is the discretization of a scalar function f(x) with x in the computational domain.
 *  More specifically, a scalar value is stored at discrete nodes/points. The nodes are arranged in an equidistant mesh
 *  with specified mesh width.
 */

class FieldVariable : public Array2D {
public:
    //! constructor
    FieldVariable(std::array<int,2> size, std::array<int,2> origin, std::array<int,2> meshWidth);

    /**
     * Interpolate FieldVariable at given position x and y
     */
    double interpolateAt(double x, double y) const;

private:
    const std::array<int,2> origin_; //< origin of computational domain
    const std::array<int,2> meshWidth_; //< meshWidth in x and y direction
};
  
