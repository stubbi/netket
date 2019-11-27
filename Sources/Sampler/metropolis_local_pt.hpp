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

#ifndef NQS_METROPOLISFLIPT_HPP
#define NQS_METROPOLISFLIPT_HPP

#include <mpi.h>
#include <Eigen/Dense>
#include <iostream>
#include "Utils/parallel_utils.hpp"
#include "Utils/random_utils.hpp"
#include "abstract_sampler.hpp"

namespace nqs {

// Metropolis sampling generating local changes
// Parallel tempering is also used
class MetropolisLocalPt : public AbstractSampler {
  // number of visible units
  const int nv_;

  // states of visible units
  // for each sampled temperature
  std::vector<Eigen::VectorXd> v_;

  Eigen::VectorXd accept_;
  Eigen::VectorXd moves_;

  int mynode_;
  int totalnodes_;

  // clusters to do updates
  std::vector<std::vector<int>> clusters_;

  // Look-up tables
  std::vector<typename AbstractMachine::LookupType> lt_;

  int nrep_;

  std::vector<double> beta_;

  int nstates_;
  std::vector<double> localstates_;

  int sweep_size_;

 public:
  // Constructor with one replica by default
  explicit MetropolisLocalPt(AbstractMachine& psi, int nreplicas = 1)
      : AbstractSampler(psi), nv_(GetHilbert().Size()), nrep_(nreplicas) {
    Init();
  }

  void Init() {
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode_);

    nstates_ = GetHilbert().LocalSize();
    localstates_ = GetHilbert().LocalStates();

    SetNreplicas(nrep_);

    // Always use odd sweep size to avoid possible ergodicity problems
    if (nv_ % 2 == 0) {
      sweep_size_ = nv_ + 1;
    } else {
      sweep_size_ = nv_;
    }

    InfoMessage() << "Metropolis sampler with parallel tempering is ready "
                  << std::endl;
    InfoMessage() << "Nreplicas is equal to " << nrep_ << std::endl;
  }

  void SetNreplicas(int nrep) {
    nrep_ = nrep;
    v_.resize(nrep_);
    for (int i = 0; i < nrep_; i++) {
      v_[i].resize(nv_);
    }

    for (int i = 0; i < nrep_; i++) {
      beta_.push_back(1. - double(i) / double(nrep_));
    }

    lt_.resize(nrep_);

    accept_.resize(2 * nrep_);
    moves_.resize(2 * nrep_);

    Reset(true);
  }

  void Reset(bool initrandom = false) override {
    if (initrandom) {
      for (int i = 0; i < nrep_; i++) {
        GetHilbert().RandomVals(v_[i], this->GetRandomEngine());
      }
    }

    for (int i = 0; i < nrep_; i++) {
      GetMachine().InitLookup(v_[i], lt_[i]);
    }

    accept_ = Eigen::VectorXd::Zero(2 * nrep_);
    moves_ = Eigen::VectorXd::Zero(2 * nrep_);
  }

  // Exchange sweep at given temperature
  void LocalSweep(int rep) {
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
      while (std::abs(newconf[0] - v_[rep](si)) <
             std::numeric_limits<double>::epsilon()) {
        newstate = diststate(this->GetRandomEngine());
        newconf[0] = localstates_[newstate];
      }

      const auto lvd =
          GetMachine().LogValDiff(v_[rep], tochange, newconf, lt_[rep]);
      double ratio = this->GetMachineFunc()(std::exp(beta_[rep] * lvd));

#ifndef NDEBUG
      const auto psival1 = GetMachine().LogVal(v_[rep]);
      if (std::abs(std::exp(GetMachine().LogVal(v_[rep]) -
                            GetMachine().LogVal(v_[rep], lt_[rep])) -
                   1.) > 1.0e-8) {
        std::cerr << GetMachine().LogVal(v_[rep]) << "  and LogVal with Lt is "
                  << GetMachine().LogVal(v_[rep], lt_[rep]) << std::endl;
        std::abort();
      }
#endif
      // Metropolis acceptance test
      if (ratio > distu(this->GetRandomEngine())) {
        accept_(rep) += 1;

        GetMachine().UpdateLookup(v_[rep], tochange, newconf, lt_[rep]);
        GetHilbert().UpdateConf(v_[rep], tochange, newconf);

#ifndef NDEBUG
        const auto psival2 = GetMachine().LogVal(v_[rep]);
        if (std::abs(std::exp(psival2 - psival1 - lvd) - 1.) > 1.0e-8) {
          std::cerr << psival2 - psival1 << " and logvaldiff is " << lvd
                    << std::endl;
          std::cerr << psival2 << " and LogVal with Lt is "
                    << GetMachine().LogVal(v_[rep], lt_[rep]) << std::endl;
          std::abort();
        }
#endif
      }
      moves_(rep) += 1;
    }
  }

  void Sweep() override {
    // First we do local sweeps
    for (int i = 0; i < nrep_; i++) {
      LocalSweep(i);
    }

    // Tempearture exchanges
    std::uniform_real_distribution<double> distribution(0, 1);

    for (int r = 1; r < nrep_; r += 2) {
      if (ExchangeProb(r, r - 1) > distribution(this->GetRandomEngine())) {
        Exchange(r, r - 1);
        accept_(nrep_ + r) += 1.;
        accept_(nrep_ + r - 1) += 1;
      }
      moves_(nrep_ + r) += 1.;
      moves_(nrep_ + r - 1) += 1;
    }

    for (int r = 2; r < nrep_; r += 2) {
      if (ExchangeProb(r, r - 1) > distribution(this->GetRandomEngine())) {
        Exchange(r, r - 1);
        accept_(nrep_ + r) += 1.;
        accept_(nrep_ + r - 1) += 1;
      }
      moves_(nrep_ + r) += 1.;
      moves_(nrep_ + r - 1) += 1;
    }
  }

  // computes the probability to exchange two replicas
  double ExchangeProb(int r1, int r2) {
    const double lf1 = 2 * std::real(GetMachine().LogVal(v_[r1], lt_[r1]));
    const double lf2 = 2 * std::real(GetMachine().LogVal(v_[r2], lt_[r2]));

    return std::exp((beta_[r1] - beta_[r2]) * (lf2 - lf1));
  }

  void Exchange(int r1, int r2) {
    std::swap(v_[r1], v_[r2]);
    std::swap(lt_[r1], lt_[r2]);
  }

  const Eigen::VectorXd& Visible() const noexcept override { return v_[0]; }

  void SetVisible(const Eigen::VectorXd& v) override { v_[0] = v; }

  AbstractMachine::VectorType DerLogVisible() override {
    return GetMachine().DerLog(v_[0], lt_[0]);
  }

  Eigen::VectorXd Acceptance() const override {
    Eigen::VectorXd acc = accept_;
    for (int i = 0; i < acc.size(); i++) {
      acc(i) /= moves_(i);
    }
    return acc;
  }
};

}  // namespace nqs

#endif
