#pragma once

#include <vector>
#include <array>
#include "../storage/FieldVariable.h"

class StaggeredGrid
{
public:
    //! constructor
    StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth);

    //! get the mesh width in x and y direction
    const std::array<double,2> meshWidth() const;

    //! get the number of cells in x and y direction
    const std::array<int,2> nCells() const;

    //! get a reference to the field variable u
    const FieldVariable &u() const;

    //! get a reference to the field variable v
    const FieldVariable &v() const;

    //! get a reference to the field variable p
    const FieldVariable &p() const;

    //! get value of u at index i,j
    double u(int i, int j) const;

    //! get reference to u at index i,j
    double &u(int i, int j);

    //! get value of v at index i,j
    double v(int i, int j) const;

    //! get reference to v at index i,j
    double &v(int i, int j);
 
    //! get value of p at index i,j
    double p(int i, int j) const;

    //! get reference to p at index i,j
    double &p(int i, int j);
 
    //! get reference to the rhs at index i,j
    double &rhs(int i, int j);

    //! get reference to F at index i,j
    double &f(int i, int j);

    //! get reference to G at index i,j
    double &g(int i, int j);

    //! get the mesh width in x direction
    double dx() const;

    //! get the mesh width in y direction
    double dy() const;

    //! get the first valid index for u in x direction
    int uIBegin() const;

    //! get the last valid index for u in x direction
    int uIEnd() const;

    //! get the first valid index for u in y direction
    int uJBegin() const;

    //! get the last valid index for u in y direction
    int uJEnd() const;

    //! get the first valid index for v in x direction
    int vIBegin() const;

    //! get the last valid index for v in x direction
    int vIEnd() const;

    //! get the first valid index for v in y direction
    int vJBegin() const;

    //! get the last valid index for v in y direction
    int vJEnd() const;

    //! get the first valid index for p in x direction
    int pIBegin() const;

    //! get the last valid index for p in x direction
    int pIEnd() const;

    //! get the first valid index for p in y direction
    int pJBegin() const;

    //! get the last valid index for p in y direction
    int pJEnd() const;


protected:
    FieldVariable u_;
    FieldVariable v_;
    FieldVariable p_;
    FieldVariable rhs_;
    FieldVariable f_;
    FieldVariable g_;

    const std::array< int, 2 > nCells_;
    const std::array< double, 2 > meshWidth_;
};