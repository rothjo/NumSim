#pragma once

#include "discretization.h"

/**
 * Calculate derivatives needed for donor cell
 */
class DonorCell : public Discretization {
public:
    DonorCell(std::array<int, 2> nCells, std::array<double,2> meshWidth, std::shared_ptr<Partitioning> partitioning_, double alpha);

    /**
     * Compute the derivative of u^2 in x-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return derivative of u^2 in x-direction
     */
    virtual double computeDu2Dx (int i, int j) const;

    /**
     * Compute the derivative of u*v in x-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return derivative of u*v in x-direction
     */
    virtual double computeDuvDx (int i, int j) const;

    /**
     * Compute the derivative of v^2 in y-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return derivative of v^2 in y-direction
     */
    virtual double computeDv2Dy (int i, int j) const;

    /**
     * Compute the derivative of u*v in y-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return derivative of u*v in y-direction
     */
    virtual double computeDuvDy (int i, int j) const;

    /**
     * Compute the derivative of u*T in x-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return derivative of u*v in y-direction
     */
    virtual double computeDuTDx (int i, int j) const;

    /**
     * Compute the derivative of v*T in y-direction
     * 
     * @param i index in x-direction
     * @param j index in y-direction
     * @return derivative of u*v in y-direction
     */
    virtual double computeDvTDy (int i, int j) const;

private:
    double alpha_;
    double gamma_;
};