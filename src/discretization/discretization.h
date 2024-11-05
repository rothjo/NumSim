#pragma once

#include <array>
#include "staggeredGrid.h"

/**
 * Calculate derivatives
 */

class Discretization : public StaggeredGrid {
public:
    Discretization(std::array<int, 2> nCells, std::array<double,2> meshWidth);

    /**
     * Compute the derivative of u^2 in x-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return derivative of u^2 in x-direction
     */

    virtual double computeDu2Dx (int i, int j) const = 0;

    /**
     * Compute the derivative of u*v in x-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return derivative of u*v in x-direction
     */

    virtual double computeDuvDx (int i, int j) const = 0;

    /**
     * Compute the derivative of p in x-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return derivative of p in x-direction
     */

    virtual double computeDpDx (int i, int j) const;

    /**
     * Compute the derivative of v^2 in y-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return derivative of v^2 in y-direction
     */

    virtual double computeDv2Dy (int i, int j) const = 0;

    /**
     * Compute the derivative of u*v in y-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return derivative of u*v in y-direction
     */

    virtual double computeDuvDy (int i, int j) const = 0;

    /**
     * Compute the derivative of p in y-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return derivative of p in y-direction
     */

    virtual double computeDpDy (int i, int j) const;

    /**
     * Compute the second derivative of u in x-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return second derivative of u in x-direction
     */

    virtual double computeD2uDx2 (int i, int j) const;

    /**
     * Compute the second derivative of u in y-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return second derivative of u in y-direction
     */

    virtual double computeD2uDy2 (int i, int j) const;

    /**
     * Compute the second derivative of v in x-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return second derivative of v in x-direction
     */

    virtual double computeD2vDx2 (int i, int j) const;

    /**
     * Compute the second derivative of v in y-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return second derivative of v in y-direction
     */

    virtual double computeD2vDy2 (int i, int j) const;

};