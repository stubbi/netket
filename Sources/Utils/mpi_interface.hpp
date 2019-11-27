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

#ifndef NQS_MPI_INTERFACE_HPP
#define NQS_MPI_INTERFACE_HPP

#include <cassert>
#include <complex>
#include <vector>

#include <mpi.h>
#include <Eigen/Dense>

#include "common_types.hpp"

namespace nqs {

inline void SendToAll(double &val, int sendnode = 0,
                      const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Bcast(&val, 1, MPI_DOUBLE, sendnode, comm);
}
inline void SendToAll(int &val, int sendnode = 0,
                      const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Bcast(&val, 1, MPI_INT, sendnode, comm);
}
inline void SendToAll(Complex &val, int sendnode = 0,
                      const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Bcast(&val, 1, MPI_DOUBLE_COMPLEX, sendnode, comm);
}

inline void SendToAll(std::vector<int> &value, int root = 0,
                      const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Bcast(&value[0], value.size(), MPI_INT, root, comm);
}
inline void SendToAll(std::vector<unsigned int> &value, int root = 0,
                      const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Bcast(&value[0], value.size(), MPI_UNSIGNED, root, comm);
}
inline void SendToAll(std::vector<unsigned long> &value, int root = 0,
                      const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Bcast(&value[0], value.size(), MPI_UNSIGNED_LONG, root, comm);
}
inline void SendToAll(std::vector<double> &value, int root = 0,
                      const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Bcast(&value[0], value.size(), MPI_DOUBLE, root, comm);
}
inline void SendToAll(std::vector<Complex> &value, int root = 0,
                      const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Bcast(&value[0], value.size(), MPI_DOUBLE_COMPLEX, root, comm);
}
inline void SendToAll(Eigen::VectorXi &value, int root = 0,
                      const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Bcast(value.data(), value.size(), MPI_INT, root, comm);
}
inline void SendToAll(Eigen::VectorXd &value, int root = 0,
                      const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Bcast(value.data(), value.size(), MPI_DOUBLE, root, comm);
}
inline void SendToAll(Eigen::VectorXcd &value, int root = 0,
                      const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Bcast(value.data(), value.size(), MPI_DOUBLE_COMPLEX, root, comm);
}

// Accumulates the sum of val collected from all nodes and the sum is
// distributed back to all processors
inline void SumOnNodes(double &val, double &sum,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(&val, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
}

inline void SumOnNodes(int &val, int &sum,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(&val, &sum, 1, MPI_INT, MPI_SUM, comm);
}
inline void SumOnNodes(int &value, const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_SUM, comm);
}

inline void SumOnNodes(Complex &val, Complex &sum,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(&val, &sum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
}

inline void SumOnNodes(std::vector<double> &val, std::vector<double> &sum,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(&val[0], &sum[0], val.size(), MPI_DOUBLE, MPI_SUM, comm);
}

inline void SumOnNodes(std::vector<double> &value,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(MPI_IN_PLACE, &value[0], value.size(), MPI_DOUBLE, MPI_SUM,
                comm);
}

#if 0
inline void SumOnNodes(std::valarray<double> &val, std::valarray<double> &sum,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(&val[0], &sum[0], val.size(), MPI_DOUBLE, MPI_SUM, comm);
}

inline void SumOnNodes(std::valarray<double> &value,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(MPI_IN_PLACE, &value[0], value.size(), MPI_DOUBLE, MPI_SUM,
                comm);
}
#endif

inline void SumOnNodes(Eigen::VectorXd &val, Eigen::VectorXd &sum,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  assert(sum.size() >= val.size());
  MPI_Allreduce(val.data(), sum.data(), val.size(), MPI_DOUBLE, MPI_SUM, comm);
}

inline void SumOnNodes(Eigen::VectorXd &value,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(MPI_IN_PLACE, value.data(), value.size(), MPI_DOUBLE, MPI_SUM,
                comm);
}

inline void SumOnNodes(Eigen::MatrixXd &value,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(MPI_IN_PLACE, value.data(), value.size(), MPI_DOUBLE, MPI_SUM,
                comm);
}

inline void SumOnNodes(Eigen::MatrixXcd &value,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(MPI_IN_PLACE, value.data(), value.size(), MPI_DOUBLE_COMPLEX,
                MPI_SUM, comm);
}

inline void SumOnNodes(double &value, const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_SUM, comm);
}

inline void SumOnNodes(double *value, int size,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(MPI_IN_PLACE, value, size, MPI_DOUBLE, MPI_SUM, comm);
}

inline void SumOnNodes(Complex *value, int size,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(MPI_IN_PLACE, value, size, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
}

inline void SumOnNodes(std::vector<Complex> &val, std::vector<Complex> &sum,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(&val[0], &sum[0], val.size(), MPI_DOUBLE_COMPLEX, MPI_SUM,
                comm);
}

inline void SumOnNodes(std::vector<Complex> &value,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(MPI_IN_PLACE, &value[0], value.size(), MPI_DOUBLE_COMPLEX,
                MPI_SUM, comm);
}

inline void SumOnNodes(Complex &value, const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
}

inline void SumOnNodes(Eigen::VectorXcd &value,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  MPI_Allreduce(MPI_IN_PLACE, value.data(), value.size(), MPI_DOUBLE_COMPLEX,
                MPI_SUM, comm);
}

inline void SumOnNodes(Eigen::VectorXcd &val, Eigen::VectorXcd &sum,
                       const MPI_Comm comm = MPI_COMM_WORLD) {
  assert(sum.size() >= val.size());
  MPI_Allreduce(val.data(), sum.data(), val.size(), MPI_DOUBLE_COMPLEX, MPI_SUM,
                comm);
}

template <typename T>
inline void MeanOnNodes(T &value, const MPI_Comm comm = MPI_COMM_WORLD) {
  int size;
  MPI_Comm_size(comm, &size);
  SumOnNodes(value);
  value /= size;
}

struct MPIHelpers {
  static int MPIRank() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
  }

  static int MPISize() {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
  }
};

}  // namespace nqs

#endif
