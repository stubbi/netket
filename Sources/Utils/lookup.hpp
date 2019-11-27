// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef NQS_LOOKUP_HPP
#define NQS_LOOKUP_HPP

#include <cassert>
#include <vector>

#include <Eigen/Core>

namespace nqs {

// Generic look-up table
template <class T>
class Lookup {
  using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

  std::vector<std::vector<VectorType>> vv_;
  std::vector<VectorType> v_;
  std::vector<MatrixType> m_;

 public:
  // Lookup<T> &operator=(const Lookup<T> lt) {}

  int AddVV(int a) {
    vv_.push_back(std::vector<VectorType>(a));
    return vv_.size() - 1;
  }

  int AddVector(int a) {
    v_.push_back(VectorType(a));
    return v_.size() - 1;
  }

  int AddMatrix(int a, int b) {
    m_.push_back(MatrixType(a, b));
    return m_.size() - 1;
  }

  int VVSize() { return vv_.size(); }

  int VectorSize() { return v_.size(); }

  int MatrixSize() { return m_.size(); }

  std::vector<VectorType> &VV(std::size_t i) {
    assert(i < vv_.size() && i >= 0);
    return vv_[i];
  }

  const std::vector<VectorType> &VV(std::size_t i) const {
    assert(i < vv_.size() && i >= 0);
    return vv_[i];
  }

  VectorType &V(std::size_t i) {
    assert(i < v_.size() && i >= 0);
    return v_[i];
  }

  const VectorType &V(std::size_t i) const {
    assert(i < v_.size() && i >= 0);
    return v_[i];
  }

  MatrixType &M(std::size_t i) {
    assert(i < m_.size() && i >= 0);
    return m_[i];
  }

  const MatrixType &M(std::size_t i) const {
    assert(i < m_.size() && i >= 0);
    return m_[i];
  }
};

}  // namespace nqs

#endif  // NQS_LOOKUP_HPP
