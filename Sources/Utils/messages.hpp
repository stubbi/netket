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

#ifndef NQS_MESSAGES_HPP
#define NQS_MESSAGES_HPP

#include <mpi.h>
#include <iostream>
#include <streambuf>
#include <string>

namespace nqs {

class NullBuffer : public std::streambuf {
 public:
  int overflow(int c) { return c; }
};

class NullStream : public std::ostream {
 public:
  NullStream() : std::ostream(&m_sb) {}

 private:
  NullBuffer m_sb;
};

inline std::ostream& InfoMessage(const std::string& comment = "# ") {
  // null stream
  static NullStream nullstream;

  // get MPI rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    return std::cout << comment;
  } else {
    return nullstream;
  }
}

inline std::ostream& WarningMessage(
    const std::string& comment = "# WARNING: ") {
  // null stream
  static NullStream nullstream;

  // get MPI rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    return std::cerr << comment;
  } else {
    return nullstream;
  }
}

inline std::ostream& ErrorMessage(const std::string& comment = "# ERROR: ") {
  // null stream
  static NullStream nullstream;

  // get MPI rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    return std::cerr << comment;
  } else {
    return nullstream;
  }
}

inline std::ostream& DebugMessage(const std::string& comment = "# DEBUG: ") {
  // null stream
  static NullStream nullstream;

#ifdef NDEBUG
  (void)comment;
  return nullstream;
#else
  // get MPI rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    return std::cout << comment;
  } else {
    return nullstream;
  }
#endif
}

}  // namespace nqs

#endif  // NQS_MESSAGES_HPP
