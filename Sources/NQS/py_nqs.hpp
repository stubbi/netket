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

#ifndef NQS_PYNQS_HPP
#define NQS_PYNQS_HPP

#include <mpi.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <complex>
#include <vector>
#include "nqs.hpp"

namespace py = pybind11;

namespace nqs {

void AddNQSModule(py::module &m) {
  auto subm = m.def_submodule("nqs");

  py::class_<NQS>(
      subm, "NQS",
      R"EOF(Neural Quantum state for classical simulation of quantum circuits using RBMs.)EOF")
    .def(py::init([](int nqubits, int initialHidden, int sampleSteps, int randomRestarts, bool earlyStopping, std::string &optimizer) {
                    return NQS{nqubits, initialHidden, sampleSteps, randomRestarts, earlyStopping, optimizer};
           }),
           py::arg("nqubits"),
           py::arg("initialHidden"),
           py::arg("sampleSteps"),
           py::arg("randomRestarts"),
           py::arg("earlyStopping"),
           py::arg("optimizer"),
           R"EOF(
           Construct a NQS object with the given number of qubits.

           Args:
               nqubits: Number of Qubits of the circuit to be simulated.
               initialHidden: Initial number of hidden units of the underlying Boltzmann machine.
               sampleSteps: The number of sample steps in the MCMC for the Hadamard. 
               randomRestarts: If true, small random permutations will be applied to the RBM's parameters. The training will be performed X times and the parameters with lowest testing error will be set. 
               earlyStopping: If true, the training will stop as soon as the testing error increases from one iteration to the next one. 
               optimizer: The optimizer. One of XXX. 
           )EOF")
      .def("applyHadamard", &NQS::applyHadamard, py::arg("qubit"),
                            py::arg("numSamples"), py::arg("numIterations"),
            R"EOF(
           Apply Hadamard gate to qubit.

           Args:
               qubit: The index of the qubit the gate will be applied to
           )EOF")
      .def("applySqrtX", &NQS::applySqrtX, py::arg("qubit"),
                            py::arg("numSamples"), py::arg("numIterations"),
            R"EOF(
           Apply sqrtX gate to qubit.

           Args:
               qubit: The index of the qubit the gate will be applied to
           )EOF")
      .def("applySqrtY", &NQS::applySqrtY, py::arg("qubit"),
                            py::arg("numSamples"), py::arg("numIterations"),
            R"EOF(
           Apply sqrtY gate to qubit.

           Args:
               qubit: The index of the qubit the gate will be applied to
           )EOF")
      .def("applyPauliX", &NQS::applyPauliX, py::arg("qubit"),
            R"EOF(
           Apply Pauli X gate to qubit.

           Args:
               qubit: The index of the qubit the gate will be applied to
           )EOF")
      .def("applyPauliY", &NQS::applyPauliY, py::arg("qubit"),
            R"EOF(
           Apply Pauli Y gate to qubit.

           Args:
               qubit: The index of the qubit the gate will be applied to
           )EOF")
      .def("applyPauliZ", &NQS::applyPauliZ, py::arg("qubit"),
            R"EOF(
           Apply Pauli Z gate to qubit.

           Args:
               qubit: The index of the qubit the gate will be applied to
           )EOF")
      .def("applySingleZRotation", &NQS::applySingleZRotation,
           py::arg("qubit"), py::arg("theta"), R"EOF(
           Apply Z rotation to qubit.

           Args:
               qubit: The index of the qubit the gate will be applied to
               theta: angle of the rotation

           )EOF")
      .def("applyControlledZRotation", &NQS::applyControlledZRotation, py::arg("controlQubit"),
            py::arg("qubit"), py::arg("theta"), py::arg("numSamples"), py::arg("numIterations"),
            R"EOF(
           Apply Controlled Z rotation to qubit.

           Args:
               qubit: The index of the qubit the gate will be applied to
               controlQubit: The index of the qubit depending on which value the rotation will be applied
               theta: angle of the rotation
               numSamples: number of training samples
               numIterations: number of training iterations
           )EOF")
      .def("applyT", &NQS::applyT, py::arg("qubit"),
            R"EOF(
           Apply T gate as defined in Nielsen and Chuang to qubit.

           Args:
               qubit: The index of the qubit the gate will be applied to
           )EOF")
      .def("applyTDagger", &NQS::applyTDagger, py::arg("qubit"),
            R"EOF(
           Apply dagger of T gate as defined in Nielsen and Chuang to qubit.

           Args:
               qubit: The index of the qubit the gate will be applied to
           )EOF")
      .def("sample", &NQS::sample,
            R"EOF(
           Sample from the nqs.
           )EOF")
      .def("save", &NQS::savePsiParams, py::arg("filename"),
           R"EOF(
                 Member function to save the machine parameters.

                 Args:
                     filename: name of file to save parameters to.
           )EOF")
      .def("load", &NQS::loadPsiParams, py::arg("filename"),
           R"EOF(
                 Member function to load machine parameters from a json file.

                 Args:
                     filename: name of file to load parameters from.
           )EOF")
      .def("psi", &NQS::psi, py::arg("v"),
            R"EOF(
           Get the psi value of a given configuration of qubits.

           Args:
               v: configuration of qubits to calculate wave function for.
           )EOF");
}

}  // namespace nqs

#endif
