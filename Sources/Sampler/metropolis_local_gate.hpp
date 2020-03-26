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

#ifndef NQS_METROPOLISLOCALHADAMARD_HPP
#define NQS_METROPOLISLOCALHADAMARD_HPP

#include <mpi.h>
#include <Eigen/Dense>
#include <iostream>
#include <limits>
#include "Utils/parallel_utils.hpp"
#include "Utils/random_utils.hpp"
#include "abstract_sampler.hpp"

namespace nqs {

// Metropolis sampling generating local moves in hilbert space on RBM after Hadamard gate has been applied
class MetropolisLocalGate : public AbstractSampler {
  // number of visible units
  const int nv_;

  // states of visible units
  Eigen::VectorXd v_;

  Eigen::VectorXd accept_;
  Eigen::VectorXd moves_;

  //matrix representation of the gate
  using MatrixType = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;
  MatrixType gatematrix_;

  int mynode_;
  int totalnodes_;

  // Look-up tables
  typename AbstractMachine::LookupType lt_;

  int nstates_;
  std::vector<double> localstates_;

  int sweep_size_;

 public:
  explicit MetropolisLocalGate(AbstractMachine& psi, MatrixType gate_matrix, int sweep_size)
      : AbstractSampler(psi), nv_(GetHilbert().Size()), gatematrix_(gate_matrix), sweep_size_(sweep_size) {
    Init();
  }

  void Init() {
    v_.resize(nv_);

    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode_);

    if (!GetHilbert().IsDiscrete()) {
      throw InvalidInputError(
          "Local Metropolis Gate sampler works only for discrete "
          "Hilbert spaces");
    }

    accept_.resize(1);
    moves_.resize(1);

    nstates_ = GetHilbert().LocalSize();
    localstates_ = GetHilbert().LocalStates();

    Reset(true);
  }

  void Reset(bool initrandom) override {
    if (initrandom) {
      GetHilbert().RandomVals(v_, this->GetRandomEngine());
    }

    GetMachine().InitLookup(v_, lt_);

    accept_ = Eigen::VectorXd::Zero(1);
    moves_ = Eigen::VectorXd::Zero(1);
  }

  void Sweep() override {}

  std::complex<double> PsiAfterGate(Eigen::VectorXd v, int qubit) {
    double valueOfQubit = v(qubit);
    //set qubit to 0
    v(qubit) = 0.0;
    std::complex<double> psi0 = std::exp(GetMachine().LogVal(v));
    //set qubit to 1
    v(qubit) = 1.0;
    std::complex<double> psi1 = std::exp(GetMachine().LogVal(v));
    v(qubit) = valueOfQubit;

    std::complex<double> psi;

    if(valueOfQubit == 0.0) {
      psi = gatematrix_(0,0) * psi0 + gatematrix_(1,0) * psi1;
    } else {
      psi = gatematrix_(0,1) * psi0 + gatematrix_(1,1) *  psi1;
    }
    return psi;
  }


  void Sweep(int qubit) {
    std::vector<int> tochange(1);
    std::vector<double> newconf(1);

    std::uniform_real_distribution<double> distu;
    std::uniform_int_distribution<int> distrs(0, nv_ - 1);
    std::uniform_int_distribution<int> diststate(0, nstates_ - 1);

    for (int i = 0; i < sweep_size_; i++) {
      // picking a random site to be changed
      int si = distrs(this->GetRandomEngine());
      assert(si < nv_);
      tochange[0] = si;

      // picking a random state
      int newstate = diststate(this->GetRandomEngine());
      newconf[0] = localstates_[newstate];

      // make sure that the new state is not equal to the current one
      while (std::abs(newconf[0] - v_(si)) <
             std::numeric_limits<double>::epsilon()) {
        newstate = diststate(this->GetRandomEngine());
        newconf[0] = localstates_[newstate];
      }

      double psiBefore = std::norm(PsiAfterGate(v_, qubit));
      psiBefore = psiBefore > 0 ? psiBefore : 0.00000001;

      double valueOfQubitToChange = v_(tochange[0]);
      v_(tochange[0]) = newconf[0];
      double psiAfter = std::norm(PsiAfterGate(v_, qubit));
      psiAfter = psiAfter > 0 ? psiAfter : 0.00000001;

      //reset
      v_(tochange[0]) = valueOfQubitToChange;

      double ratio;
      if(psiBefore != 0.0) {
        ratio = psiAfter/ psiBefore;
      } else {
        ratio = 1.0;
      }

#ifndef NDEBUG
      const auto psival1 = GetMachine().LogVal(v_);
      if (std::abs(
              std::exp(GetMachine().LogVal(v_) - GetMachine().LogVal(v_, lt_)) -
              1.) > 1.0e-8) {
        std::cerr << GetMachine().LogVal(v_) << "  and LogVal with Lt is "
                  << GetMachine().LogVal(v_, lt_) << std::endl;
        std::abort();
      }
#endif

      // Metropolis acceptance test
      if (ratio > distu(this->GetRandomEngine())) {
        accept_[0] += 1;
        GetMachine().UpdateLookup(v_, tochange, newconf, lt_);
        GetHilbert().UpdateConf(v_, tochange, newconf);

#ifndef NDEBUG
        const auto psival2 = GetMachine().LogVal(v_);
        if (std::abs(std::exp(psival2 - psival1 - lvd) - 1.) > 1.0e-8) {
          std::cerr << psival2 - psival1 << " and logvaldiff is " << lvd
                    << std::endl;
          std::cerr << psival2 << " and LogVal with Lt is "
                    << GetMachine().LogVal(v_, lt_) << std::endl;
          std::abort();
        }
#endif
      }
      moves_[0] += 1;
    }
  }

  const Eigen::VectorXd& Visible() const noexcept override { return v_; }

  void SetVisible(const Eigen::VectorXd& v) override { v_ = v; }

  AbstractMachine::VectorType DerLogVisible() override {
    return GetMachine().DerLog(v_, lt_);
  }

  Eigen::VectorXd Acceptance() const override {
    Eigen::VectorXd acc = accept_;
    for (int i = 0; i < 1; i++) {
      acc(i) /= moves_(i);
    }
    return acc;
  }
};

}  // namespace nqs

#endif
