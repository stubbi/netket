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

#ifndef NQS_PYSUPERVISED_HPP
#define NQS_PYSUPERVISED_HPP

#include <mpi.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <complex>
#include <vector>
#include "supervised.hpp"

namespace py = pybind11;

namespace nqs {

void AddSupervisedModule(py::module &m) {
  auto subm = m.def_submodule("supervised");

  py::class_<Supervised>(
      subm, "Supervised",
      R"EOF(Supervised learning scheme to learn data, i.e. the given state, by stochastic gradient descent with log overlap loss or MSE loss.)EOF")
      .def(py::init([](AbstractMachine &ma,
                       MetropolisLocal &sa,
                       AbstractOptimizer &opt,
                       int batch_size,
                       std::vector<Eigen::VectorXd> trainingSamples,
                       std::vector<Eigen::VectorXcd> trainingTargets,
                       std::vector<Eigen::VectorXd> testSamples,
                       std::vector<Eigen::VectorXcd> testTargets,
                       bool sr,
                       double diag_shift,
                       bool use_iterative, bool use_cholesky) {
             return Supervised{ma,
                               sa,
                               opt,
                               batch_size,
                               std::move(trainingSamples),
                               std::move(trainingTargets),
                               std::move(testSamples),
                               std::move(testTargets),
                               sr,
                               diag_shift,
                               use_iterative,
                               use_cholesky};
           }),
           py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
           py::arg("machine"), py::arg("sampler"), py::arg("opt"),
           py::arg("batch_size"), py::arg("trainingSamples"),
           py::arg("trainingTargets"),
           py::arg("testSamples"), py::arg("testTargets"), py::arg("sr"),
           py::arg("diag_shift") = 0.01, py::arg("use_iterative") = false,
           py::arg("use_cholesky") = true,
           R"EOF(
           Construct a Supervised object given a machine, an optimizer, batch size and
           data, including samples and targets.

           Args:
               machine: The machine representing the wave function.
               optimizer: The optimizer object that determines how the SGD optimization.
               sampler: The sampler.
               batch_size: The batch size used in SGD.
               samples: The input data, i.e. many-body basis.
               targets: The output label, i.e. amplitude of the corresponding basis.
               targetNormalisation: max abs of target psi
               method: The chosen method to learn the parameters of the
                   wave-function. Possible choices are `Gd` (Regular Gradient descent),
                   and `Sr` (Stochastic reconfiguration a.k.a. natural gradient).
               diag_shift: The regularization parameter in stochastic
                   reconfiguration. The default is 0.01.
               use_iterative: Whether to use the iterative solver in the Sr
                   method (this is extremely useful when the number of
                   parameters to optimize is very large). The default is false.
               use_cholesky: Whether to use cholesky decomposition. The default
                   is true.

           )EOF")
      .def("run", &Supervised::Run, py::arg("n_iter"),
           py::arg("earlyStopping"),
           py::arg("output_prefix") = "output",
           py::arg("save_params_every") = 50, R"EOF(
           Run supervised learning.

           Args:
               n_iter: The number of iterations for running.
               earlyStopping: If true, training will be stopped as soon as testing error increases during training iterations
               output_prefix: The output file name, without extension.
               save_params_every: Frequency to dump wavefunction parameters. The default is 50.

           )EOF");
}

}  // namespace nqs

#endif
