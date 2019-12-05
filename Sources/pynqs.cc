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

#ifndef NQS_PYNQS_CC
#define NQS_PYNQS_CC

#include <nqs.hpp>
#include "NQS/py_nqs.hpp"
#include "Dynamics/py_dynamics.hpp"
#include "Graph/py_graph.hpp"
#include "GroundState/py_ground_state.hpp"
#include "Hilbert/py_hilbert.hpp"
#include "Machine/py_machine.hpp"
#include "Operator/py_operator.hpp"
#include "Optimizer/py_optimizer.hpp"
#include "Output/py_output.hpp"
#include "Sampler/py_sampler.hpp"
#include "Stats/py_stats.hpp"
#include "Supervised/py_supervised.hpp"
#include "Unsupervised/py_unsupervised.hpp"
#include "Utils/mpi_interface.hpp"  // for MPIInitializer
#include "Utils/py_utils.hpp"
#include "Utils/pybind_helpers.hpp"

namespace nqs {

PYBIND11_MODULE(_C_nqs, m) {
  AddDynamicsModule(m);
  AddGraphModule(m);
  AddGroundStateModule(m);
  AddHilbertModule(m);
  AddMachineModule(m);
  AddNQSModule(m);
  AddOperatorModule(m);
  AddOptimizerModule(m);
  AddOutputModule(m);
  AddSamplerModule(m);
  AddStatsModule(m);
  AddUtilsModule(m);
  AddSupervisedModule(m);
  AddUnsupervisedModule(m);
}  // PYBIND11_MODULE

}  // namespace nqs


#endif
