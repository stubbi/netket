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

#include "rbm_nqs.hpp"

#include "Utils/json_utils.hpp"
#include "Utils/messages.hpp"

namespace nqs {

RbmNQS::RbmNQS(std::shared_ptr<const AbstractHilbert> hilbert, int nhidden,
                 int alpha, bool usea, bool useb)
    : AbstractMachine(hilbert), nv_(hilbert->Size()), usea_(usea), useb_(useb) {
  nh_ = std::max(nhidden, alpha * nv_);
  Init();
}

int RbmNQS::Nvisible() const { return nv_; }

int RbmNQS::Npar() const { return npar_; }

void RbmNQS::addHidden() {
  b_.conservativeResize(b_.size() + 1);
  b_(b_.size() - 1) = 0;

  W_.conservativeResize(W_.rows(), W_.cols() + 1);
  W_.col(W_.cols()-1).setZero();
  
  nh_ += 1;

  thetas_.resize(nh_);
  lnthetas_.resize(nh_);
  thetasnew_.resize(nh_);
  lnthetasnew_.resize(nh_);

  npar_ += W_.col(W_.cols() - 1).size() + 1;
}

void RbmNQS::Init() {
  W_.resize(nv_, nh_);
  a_.resize(nv_);
  b_.resize(nh_);

  thetas_.resize(nh_);
  lnthetas_.resize(nh_);
  thetasnew_.resize(nh_);
  lnthetasnew_.resize(nh_);

  npar_ = nv_ * nh_;

  if (usea_) {
    npar_ += nv_;
  } else {
    a_.setZero();
  }

  if (useb_) {
    npar_ += nh_;
  } else {
    b_.setZero();
  }

  InfoMessage() << "RBM Initizialized with nvisible = " << nv_
                << " and nhidden = " << nh_ << std::endl;
  InfoMessage() << "Using visible bias = " << usea_ << std::endl;
  InfoMessage() << "Using hidden bias  = " << useb_ << std::endl;
}

void RbmNQS::InitRandomPars(int seed, double sigma) {
  VectorType par(npar_);

  nqs::RandomGaussian(par, seed, sigma);

  SetParameters(par);
}

void RbmNQS::InitLookup(VisibleConstType v, LookupType &lt) {
  if (lt.VectorSize() == 0) {
    lt.AddVector(b_.size());
  }
  if (lt.V(0).size() != b_.size()) {
    lt.V(0).resize(b_.size());
  }

  lt.V(0) = (W_.transpose() * v + b_);
}

void RbmNQS::UpdateLookup(VisibleConstType v, const std::vector<int> &tochange,
                           const std::vector<double> &newconf, LookupType &lt) {
  if (tochange.size() != 0) {
    for (std::size_t s = 0; s < tochange.size(); s++) {
      const int sf = tochange[s];
      lt.V(0) += W_.row(sf) * (newconf[s] - v(sf));
    }
  }
}

RbmNQS::VectorType RbmNQS::DerLog(VisibleConstType v) {
  LookupType ltnew;
  InitLookup(v, ltnew);
  return DerLog(v, ltnew);
}

RbmNQS::VectorType RbmNQS::DerLog(VisibleConstType v, const LookupType &lt) {
  VectorType der(npar_);

  if (usea_) {
    der.head(nv_) = v;
  }

  RbmNQS::sigmoid(lt.V(0), lnthetas_);

  if (useb_) {
    der.segment(usea_ * nv_, nh_) = lnthetas_;
  }

  MatrixType wder = (v * lnthetas_.transpose());
  der.tail(nv_ * nh_) = Eigen::Map<VectorType>(wder.data(), nv_ * nh_);

  return der;
}

RbmNQS::VectorType RbmNQS::GetParameters() {
  VectorType pars(npar_);

  if (usea_) {
    pars.head(nv_) = a_;
  }

  if (useb_) {
    pars.segment(usea_ * nv_, nh_) = b_;
  }

  pars.tail(nv_ * nh_) = Eigen::Map<VectorType>(W_.data(), nv_ * nh_);

  return pars;
}

void RbmNQS::SetParameters(VectorConstRefType pars) {
  if (usea_) {
    a_ = pars.head(nv_);
  }

  if (useb_) {
    b_ = pars.segment(usea_ * nv_, nh_);
  }

  VectorType Wpars = pars.tail(nv_ * nh_);

  W_ = Eigen::Map<MatrixType>(Wpars.data(), nv_, nh_);
}

// Value of the logarithm of the wave-function
Complex RbmNQS::LogVal(VisibleConstType v) {
  RbmNQS::softplus(W_.transpose() * v + b_, lnthetas_);

  return (v.dot(a_) + lnthetas_.sum());
}

// Value of the logarithm of the wave-function
// using pre-computed look-up tables for efficiency
Complex RbmNQS::LogVal(VisibleConstType v, const LookupType &lt) {
  RbmNQS::softplus(lt.V(0), lnthetas_);

  return (v.dot(a_) + lnthetas_.sum());
}

// Difference between logarithms of values, when one or more visible variables
// are being flipped
RbmNQS::VectorType RbmNQS::LogValDiff(
    VisibleConstType v, const std::vector<std::vector<int>> &tochange,
    const std::vector<std::vector<double>> &newconf) {
  const std::size_t nconn = tochange.size();
  VectorType logvaldiffs = VectorType::Zero(nconn);

  thetas_ = (W_.transpose() * v + b_);
  RbmNQS::softplus(thetas_, lnthetas_);

  Complex logtsum = lnthetas_.sum();

  for (std::size_t k = 0; k < nconn; k++) {
    if (tochange[k].size() != 0) {
      thetasnew_ = thetas_;

      for (std::size_t s = 0; s < tochange[k].size(); s++) {
        const int sf = tochange[k][s];

        logvaldiffs(k) += a_(sf) * (newconf[k][s] - v(sf));

        thetasnew_ += W_.row(sf) * (newconf[k][s] - v(sf));
      }

      RbmNQS::softplus(thetasnew_, lnthetasnew_);
      logvaldiffs(k) += lnthetasnew_.sum() - logtsum;
    }
  }
  return logvaldiffs;
}

// Difference between logarithms of values, when one or more visible variables
// are being flipped Version using pre-computed look-up tables for efficiency
// on a small number of spin flips
Complex RbmNQS::LogValDiff(VisibleConstType v,
                            const std::vector<int> &tochange,
                            const std::vector<double> &newconf,
                            const LookupType &lt) {
  Complex logvaldiff = 0.;

  if (tochange.size() != 0) {
    RbmNQS::softplus(lt.V(0), lnthetas_);

    thetasnew_ = lt.V(0);

    for (std::size_t s = 0; s < tochange.size(); s++) {
      const int sf = tochange[s];

      logvaldiff += a_(sf) * (newconf[s] - v(sf));

      thetasnew_ += W_.row(sf) * (newconf[s] - v(sf));
    }

    RbmNQS::softplus(thetasnew_, lnthetasnew_);
    logvaldiff += (lnthetasnew_.sum() - lnthetas_.sum());
  }
  return logvaldiff;
}

void RbmNQS::Save(const std::string &filename) const {
  json state;
  state["Name"] = "RbmNQS";
  state["Nvisible"] = nv_;
  state["Nhidden"] = nh_;
  state["UseVisibleBias"] = usea_;
  state["UseHiddenBias"] = useb_;
  state["a"] = a_;
  state["b"] = b_;
  state["W"] = W_;
  WriteJsonToFile(state, filename);
}

void RbmNQS::Load(const std::string &filename) {
  auto const pars = ReadJsonFromFile(filename);
  std::string name = FieldVal<std::string>(pars, "Name");
  if (name != "RbmNQS") {
    throw InvalidInputError(
        "Error while constructing RbmNQS from input parameters");
  }

  if (FieldExists(pars, "Nvisible")) {
    nv_ = FieldVal<int>(pars, "Nvisible");
  }
  if (nv_ != GetHilbert().Size()) {
    throw InvalidInputError(
        "Number of visible units is incompatible with given "
        "Hilbert space");
  }

  if (FieldExists(pars, "Nhidden")) {
    nh_ = FieldVal<int>(pars, "Nhidden");
  } else {
    nh_ = nv_ * double(FieldVal<double>(pars, "Alpha"));
  }

  usea_ = FieldOrDefaultVal(pars, "UseVisibleBias", true);
  useb_ = FieldOrDefaultVal(pars, "UseHiddenBias", true);

  Init();

  // Loading parameters, if defined in the input
  if (FieldExists(pars, "a")) {
    a_ = FieldVal<VectorType>(pars, "a");
  } else {
    a_.setZero();
  }

  if (FieldExists(pars, "b")) {
    b_ = FieldVal<VectorType>(pars, "b");
  } else {
    b_.setZero();
  }
  if (FieldExists(pars, "W")) {
    W_ = FieldVal<MatrixType>(pars, "W");
  }
}

bool RbmNQS::IsHolomorphic() const noexcept { return true; }

}  // namespace nqs
