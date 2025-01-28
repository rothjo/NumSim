#include "array2d.h"

#include <stdexcept>
#include <cassert>

Array2D::Array2D(std::array<int,2> size) :
  size_(size)
{
  // allocate data, initialize to 0
  data_.resize(size_[0]*size_[1], 0.0);
}

//! get the size
std::array<int,2> Array2D::size() const
{
  return size_;
}

double &Array2D::operator()(int i, int j)
{
  const int index = j*size_[0] + i;

  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());

  return data_[index];
}

double Array2D::operator()(int i, int j) const
{
  const int index = j*size_[0] + i;

  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());

  return data_[index];
}

void Array2D::setToZero() {
  std::fill(data_.begin(), data_.end(), 0.0);
}

void* Array2D::data() {
  return data_.data();
}

Array2D& Array2D::operator=(const Array2D& other) {
    // Check if the sizes match
    if (size_ != other.size_) {
        throw std::runtime_error("Cannot assign Array2D objects of different sizes.");
    }

    // Copy the data
    data_ = other.data_;

    // Return a reference to this
    return *this;
}