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

#ifndef NETKET_PYNQS_HPP
#define NETKET_PYNQS_HPP

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

namespace netket {

void AddNQSModule(py::module &m) {
  auto subm = m.def_submodule("nqs");

  py::class_<NQS>(
      subm, "NQS",
      R"EOF(Neural Quantum state for classical simulation of quantum circuits using RBMs.)EOF")
      .def(py::init([](int nqubits) {
             return NQS{nqubits};
           }),
           py::arg("nqubits"),
           R"EOF(
           Construct a NQS object with the given number of qubits.

           Args:
               nqubits: Number of Qubits of the circuit to be simulated.
           )EOF")
      .def("applyHadamard", &NQS::applyHadamard, py::arg("qubit"),
                            py::arg("numSamples"), py::arg("numIterations"),
            R"EOF(
           Apply Hadamard gate to qubit.

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
            py::arg("qubit"), py::arg("theta"),
            R"EOF(
           Apply Controlled Z rotation to qubit.

           Args:
               qubit: The index of the qubit the gate will be applied to
               controlQubit: The index of the qubit depending on which value the rotation will be applied
               theta: angle of the rotation
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
      .def("applyToffoli", &NQS::applyToffoli,
            py::arg("qubit1"), py::arg("qubit2"), py::arg("qubit3"), py::arg("numSamples"), py::arg("numIterations"),
            R"EOF(
           Apply dagger of T gate as defined in Nielsen and Chuang to qubit.

           Args:
               qubit1: The index of the first qubit the Toffoli gate will be applied to
               qubit2: The index of the second qubit the Toffoli gate will be applied to
               qubit3: The index of the third qubit the Toffoli gate will be applied to
               numSamples: Number of samples to draw for Hadamard gates being applied along the way
               numIterations: Number of training iterations for Hadamard gates being applied along the way
           )EOF")
      .def("sample", &NQS::sample,
            R"EOF(
           Sample from the nqs.
           )EOF")
      .def("getPsiParams", &NQS::getPsiParams,
            R"EOF(
           Get parameters of the underlying Boltzmann machine.
           )EOF")
      .def("psi", &NQS::psi, py::arg("v"),
            R"EOF(
           Get the psi value of a given configuration of qubits.

           Args:
               v: configuration of qubits to calculate wave function for.
           )EOF");
}

}  // namespace netket

#endif
