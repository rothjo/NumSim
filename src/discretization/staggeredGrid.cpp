#include "staggeredGrid.h"

// nCells has to be N + 2, where N is the number of cells in x or y direction
// StaggeredGrid::StaggeredGrid(std::array<int, 2> nCells, std::array<double, 2> meshWidth)
//     : nCells_(nCells), meshWidth_(meshWidth),
//       u_(FieldVariable({nCells[0] + 2, nCells[1] + 2}, {0.0, -0.5 * meshWidth[1]}, meshWidth)),
//       v_(FieldVariable({nCells[0] + 2, nCells[1] + 2}, {-0.5 * meshWidth[0], 0.0}, meshWidth)),
//       p_(FieldVariable({nCells[0] + 2, nCells[1] + 2}, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth)),
//       rhs_(FieldVariable({nCells[0] + 2, nCells[1] + 2}, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth)),
//       f_(FieldVariable({nCells[0] + 2, nCells[1] + 2}, {0.0, -0.5 * meshWidth[1]}, meshWidth)),
//       g_(FieldVariable({nCells[0] + 2, nCells[1] + 2}, {-0.5 * meshWidth[0], 0.0}, meshWidth)) {
// }

StaggeredGrid::StaggeredGrid(std::array<int, 2> nCells, std::array<double, 2> meshWidth, std::shared_ptr<Partitioning> partitioning)
    : nCells_(nCells), meshWidth_(meshWidth), partitioning_(partitioning),
      u_(FieldVariable({nCells[0] + 3, nCells[1] + 3}, {0.0, -0.5 * meshWidth[1]}, meshWidth)),
      v_(FieldVariable({nCells[0] + 3, nCells[1] + 3}, {-0.5 * meshWidth[0], 0.0}, meshWidth)),
      p_(FieldVariable({nCells[0] + 3, nCells[1] + 3}, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth)),
      rhs_(FieldVariable({nCells[0] + 3, nCells[1] + 3}, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth)),
      f_(FieldVariable({nCells[0] + 3, nCells[1] + 3}, {0.0, -0.5 * meshWidth[1]}, meshWidth)),
      g_(FieldVariable({nCells[0] + 3, nCells[1] + 3}, {-0.5 * meshWidth[0], 0.0}, meshWidth)) {
        containsLeftBoundary = partitioning_->ownPartitionContainsLeftBoundary();
        containsRightBoundary = partitioning_->ownPartitionContainsRightBoundary();
        containsTopBoundary = partitioning_->ownPartitionContainsTopBoundary();
        containsBottomBoundary = partitioning_->ownPartitionContainsBottomBoundary();
}

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
    if(containsLeftBoundary) {
        return 2;
    }
    return 1;
}

int StaggeredGrid::uIEnd() const {
    if(containsRightBoundary) {
        return nCells_[0] + 1;
    }
    return nCells_[0] + 2;
}

int StaggeredGrid::uJBegin() const {
    // return 1;
    return 2;
}

int StaggeredGrid::uJEnd() const {
    // return nCells_[1] + 1;
    return nCells_[1] + 2;
}

int StaggeredGrid::vIBegin() const {
    // return 1;
    return 2;
}

int StaggeredGrid::vIEnd() const {
    // return nCells_[0] + 1;
    return nCells_[0] + 2;
}

int StaggeredGrid::vJBegin() const {
    if (containsBottomBoundary) {
        return 2;
    }    
    return 1;
}

int StaggeredGrid::vJEnd() const {
    if (containsTopBoundary) {
        return nCells_[1] + 1;
    }
    return nCells_[1] + 2;
}

int StaggeredGrid::pIBegin() const {
    // return 1;
    return 2;
}

int StaggeredGrid::pIEnd() const {
    // return nCells_[0] + 1;
    return nCells_[0] + 2;
}

int StaggeredGrid::pJBegin() const {
    // return 1;
    return 2;
}

int StaggeredGrid::pJEnd() const {
    // return nCells_[1] + 1;
    return nCells_[1] + 2;
}