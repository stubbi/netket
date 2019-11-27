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

#include "custom_hilbert.hpp"

namespace nqs {

CustomHilbert::CustomHilbert(const AbstractGraph &graph,
                             const std::vector<double> &localstates)
    : graph_(graph), local_(localstates) {
  size_ = graph.Size();
  nstates_ = local_.size();
}

bool CustomHilbert::IsDiscrete() const { return true; }

int CustomHilbert::LocalSize() const { return nstates_; }

int CustomHilbert::Size() const { return size_; }

std::vector<double> CustomHilbert::LocalStates() const { return local_; }

void CustomHilbert::RandomVals(Eigen::Ref<Eigen::VectorXd> state,
                               nqs::default_random_engine &rgen) const {
  std::uniform_int_distribution<int> distribution(0, nstates_ - 1);

  assert(state.size() == size_);

  // unconstrained random
  for (int i = 0; i < state.size(); i++) {
    state(i) = local_[distribution(rgen)];
  }
}

void CustomHilbert::UpdateConf(Eigen::Ref<Eigen::VectorXd> v,
                               nonstd::span<const int> tochange,
                               nonstd::span<const double> newconf) const {
  assert(v.size() == size_);

  int i = 0;
  for (auto sf : tochange) {
    v(sf) = newconf[i];
    i++;
  }
}

const AbstractGraph &CustomHilbert::GetGraph() const noexcept { return graph_; }

}  // namespace nqs
