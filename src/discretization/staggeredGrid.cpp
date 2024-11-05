#include "staggeredGrid.h"

// nCells has to be N + 2, where N is the number of cells in x or y direction
StaggeredGrid::StaggeredGrid(std::array<int, 2> nCells, std::array<double, 2> meshWidth)
    : nCells_(nCells), meshWidth_(meshWidth),
      u_(FieldVariable({nCells[0] + 2, nCells[1] + 2}, {0.0, -0.5 * meshWidth[1]}, meshWidth)),
      v_(FieldVariable({nCells[0] + 2, nCells[1] + 2}, {-0.5 * meshWidth[0], 0.0}, meshWidth)),
      p_(FieldVariable({nCells[0] + 2, nCells[1] + 2}, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth)),
      rhs_(FieldVariable({nCells[0] + 2, nCells[1] + 2}, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth)),
      f_(FieldVariable({nCells[0] + 2, nCells[1] + 2}, {0.0, -0.5 * meshWidth[1]}, meshWidth)),
      g_(FieldVariable({nCells[0] + 2, nCells[1] + 2}, {-0.5 * meshWidth[0], 0.0}, meshWidth)) {
}

// StaggeredGrid::StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth) :
//     nCells_(nCells), meshWidth_(meshWidth) {
//         std::array<int, 2> size = {nCells[0] + 2, nCells[1] + 2};

//         u_ = FieldVariable(size, {0.0, -0.5 * meshWidth[1]}, meshWidth);
//         v_ = FieldVariable(size, {-0.5 * meshWidth[0], 0.0}, meshWidth);
//         p_ = FieldVariable(size, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth);
//         rhs_ = FieldVariable(size, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth);
//         f_ = FieldVariable(size, {0.0, -0.5 * meshWidth[1]}, meshWidth);
//         g_ = FieldVariable(size, {-0.5 * meshWidth[0], 0.0}, meshWidth);
//     }

const std::array<double,2> StaggeredGrid::meshWidth() const {
    return meshWidth_;
}

const std::array<int,2> StaggeredGrid::nCells() const {
    return nCells_;
}

const FieldVariable& StaggeredGrid::u() const {
    return u_;
}

const FieldVariable& StaggeredGrid::v() const {
    return v_;
}

const FieldVariable& StaggeredGrid::p() const {
    return p_;
}

double StaggeredGrid::u(int i, int j) const {
    return u_(i, j);
}

double& StaggeredGrid::u(int i, int j) {
    return u_(i, j);
}

double StaggeredGrid::v(int i, int j) const {
    return v_(i, j);
}

double& StaggeredGrid::v(int i, int j) {
    return v_(i, j);
}

double StaggeredGrid::p(int i, int j) const {
    return p_(i, j);
}

double& StaggeredGrid::p(int i, int j) {
    return p_(i, j);
}

double& StaggeredGrid::rhs(int i, int j) {
    return rhs_(i, j);
}

double& StaggeredGrid::f(int i, int j) {
    return f_(i, j);
}

double& StaggeredGrid::g(int i, int j) {
    return g_(i, j);
}

double StaggeredGrid::dx() const {
    return meshWidth_[0];
}

double StaggeredGrid::dy() const {
    return meshWidth_[1];
}

int StaggeredGrid::uIBegin() const {
    return 1;
}

int StaggeredGrid::uIEnd() const {
    return nCells_[0];
}

int StaggeredGrid::uJBegin() const {
    return 1;
}

int StaggeredGrid::uJEnd() const {
    return nCells_[1] + 1;
}

int StaggeredGrid::vIBegin() const {
    return 1;
}

int StaggeredGrid::vIEnd() const {
    return nCells_[0] + 1;
}

int StaggeredGrid::vJBegin() const {
    return 1;
}

int StaggeredGrid::vJEnd() const {
    return nCells_[1];
}

int StaggeredGrid::pIBegin() const {
    return 1;
}

int StaggeredGrid::pIEnd() const {
    return nCells_[0] + 1;
}

int StaggeredGrid::pJBegin() const {
    return 1;
}

int StaggeredGrid::pJEnd() const {
    return nCells_[1] + 1;
}