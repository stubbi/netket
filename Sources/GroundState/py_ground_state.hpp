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

#ifndef NQS_PYGROUND_STATE_HPP
#define NQS_PYGROUND_STATE_HPP

#include <pybind11/pybind11.h>
#include "py_exact.hpp"
#include "py_variational_montecarlo.hpp"

namespace py = pybind11;

namespace nqs {

void AddGroundStateModule(py::module &m) {
  AddVariationalMonteCarloModule(m);
  AddExactModule(m);
}

}  // namespace nqs

#endif
