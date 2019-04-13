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

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "Utils/all_utils.hpp"
#include "Utils/lookup.hpp"

#ifndef NETKET_RBM_SPIN_REAL_HPP
#define NETKET_RBM_SPIN_REAL_HPP

namespace netket {

/** Restricted Boltzmann machine class with spin 1/2 hidden units.
and real-valued weights
 *
 */
class RbmSpinReal : public AbstractMachine {
  const AbstractHilbert &hilbert_;

  // number of visible units
  int nv_;

  // number of hidden units
  int nh_;

  // number of parameters
  int npar_;

  // weights
  RealMatrixType W_;

  // visible units bias
  RealVectorType a_;

  // hidden units bias
  RealVectorType b_;

  RealVectorType thetas_;
  RealVectorType lnthetas_;
  RealVectorType thetasnew_;
  RealVectorType lnthetasnew_;

  bool usea_;
  bool useb_;

 public:
  explicit RbmSpinReal(const AbstractHilbert &hilbert, int nhidden = 0,
                       int alpha = 0, bool usea = true, bool useb = true)
      : hilbert_(hilbert), nv_(hilbert.Size()), usea_(usea), useb_(useb) {
    nh_ = std::max(nhidden, alpha * nv_);

    Init();
  }

  void Init() {
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

    InfoMessage() << "Real-valued RBM Initizialized with nvisible = " << nv_
                  << " and nhidden = " << nh_ << std::endl;
    InfoMessage() << "Using visible bias = " << usea_ << std::endl;
    InfoMessage() << "Using hidden bias  = " << useb_ << std::endl;
  }

  int Nvisible() const override { return nv_; }

  int Nhidden() const { return nh_; }

  int Npar() const override { return npar_; }

  void InitRandomPars(int seed, double sigma) override {
    RealVectorType par(npar_);

    netket::RandomGaussian(par, seed, sigma);

    SetParameters(VectorType(par));
  }

  void InitLookup(VisibleConstType v, LookupType &lt) override {
    if (lt.VectorSize() == 0) {
      lt.AddVector(b_.size());
    }
    if (lt.V(0).size() != b_.size()) {
      lt.V(0).resize(b_.size());
    }

    lt.V(0) = (W_.transpose() * v + b_);
  }

  void UpdateLookup(VisibleConstType v, const std::vector<int> &tochange,
                    const std::vector<double> &newconf,
                    LookupType &lt) override {
    if (tochange.size() != 0) {
      for (std::size_t s = 0; s < tochange.size(); s++) {
        const int sf = tochange[s];
        lt.V(0) += W_.row(sf) * (newconf[s] - v(sf));
      }
    }
  }

  VectorType DerLog(VisibleConstType v) override {
    VectorType der(npar_);

    int k = 0;

    if (usea_) {
      for (; k < nv_; k++) {
        der(k) = v(k);
      }
    }

    RbmSpin::tanh(W_.transpose() * v + b_, lnthetas_);

    if (useb_) {
      for (int p = 0; p < nh_; p++) {
        der(k) = lnthetas_(p);
        k++;
      }
    }

    for (int i = 0; i < nv_; i++) {
      for (int j = 0; j < nh_; j++) {
        der(k) = lnthetas_(j) * v(i);
        k++;
      }
    }
    return der;
  }

  VectorType GetParameters() override {
    VectorType pars(npar_);

    int k = 0;

    if (usea_) {
      for (; k < nv_; k++) {
        pars(k) = a_(k);
      }
    }

    if (useb_) {
      for (int p = 0; p < nh_; p++) {
        pars(k) = b_(p);
        k++;
      }
    }

    for (int i = 0; i < nv_; i++) {
      for (int j = 0; j < nh_; j++) {
        pars(k) = W_(i, j);
        k++;
      }
    }

    return pars;
  }

  void SetParameters(VectorConstRefType pars) override {
    int k = 0;

    if (usea_) {
      for (; k < nv_; k++) {
        a_(k) = std::real(pars(k));
      }
    }

    if (useb_) {
      for (int p = 0; p < nh_; p++) {
        b_(p) = std::real(pars(k));
        k++;
      }
    }

    for (int i = 0; i < nv_; i++) {
      for (int j = 0; j < nh_; j++) {
        W_(i, j) = std::real(pars(k));
        k++;
      }
    }
  }

  // Value of the logarithm of the wave-function
  Complex LogVal(VisibleConstType v) override {
    RbmSpin::lncosh(W_.transpose() * v + b_, lnthetas_);

    return (v.dot(a_) + lnthetas_.sum());
  }

  // Value of the logarithm of the wave-function
  // using pre-computed look-up tables for efficiency
  Complex LogVal(VisibleConstType v, const LookupType &lt) override {
    RbmSpin::lncosh(lt.V(0).real(), lnthetas_);

    return (v.dot(a_) + lnthetas_.sum());
  }

  // Difference between logarithms of values, when one or more visible variables
  // are being flipped
  VectorType LogValDiff(
      VisibleConstType v, const std::vector<std::vector<int>> &tochange,
      const std::vector<std::vector<double>> &newconf) override {
    const std::size_t nconn = tochange.size();
    VectorType logvaldiffs = VectorType::Zero(nconn);

    thetas_ = (W_.transpose() * v + b_);
    RbmSpin::lncosh(thetas_, lnthetas_);

    Complex logtsum = lnthetas_.sum();

    for (std::size_t k = 0; k < nconn; k++) {
      if (tochange[k].size() != 0) {
        thetasnew_ = thetas_;

        for (std::size_t s = 0; s < tochange[k].size(); s++) {
          const int sf = tochange[k][s];

          logvaldiffs(k) += a_(sf) * (newconf[k][s] - v(sf));

          thetasnew_ += W_.row(sf) * (newconf[k][s] - v(sf));
        }

        RbmSpin::lncosh(thetasnew_, lnthetasnew_);
        logvaldiffs(k) += lnthetasnew_.sum() - logtsum;
      }
    }
    return logvaldiffs;
  }

  // Difference between logarithms of values, when one or more visible variables
  // are being flipped Version using pre-computed look-up tables for efficiency
  // on a small number of spin flips
  Complex LogValDiff(VisibleConstType v, const std::vector<int> &tochange,
                     const std::vector<double> &newconf,
                     const LookupType &lt) override {
    Complex logvaldiff = 0.;

    if (tochange.size() != 0) {
      RbmSpin::lncosh(lt.V(0).real(), lnthetas_);

      thetasnew_ = lt.V(0).real();

      for (std::size_t s = 0; s < tochange.size(); s++) {
        const int sf = tochange[s];

        logvaldiff += a_(sf) * (newconf[s] - v(sf));

        thetasnew_ += W_.row(sf) * (newconf[s] - v(sf));
      }

      RbmSpin::lncosh(thetasnew_, lnthetasnew_);
      logvaldiff += (lnthetasnew_.sum() - lnthetas_.sum());
    }
    return logvaldiff;
  }

  inline static double lncosh(double x) {
    const double xp = std::abs(x);
    if (xp <= 12.) {
      return std::log(std::cosh(xp));
    } else {
      const static double log2v = std::log(2.);
      return xp - log2v;
    }
  }

  const AbstractHilbert &GetHilbert() const noexcept override {
    return hilbert_;
  }

  void to_json(json &j) const override {
    j["Name"] = "RbmSpinReal";
    j["Nvisible"] = nv_;
    j["Nhidden"] = nh_;
    j["UseVisibleBias"] = usea_;
    j["UseHiddenBias"] = useb_;
    j["a"] = a_;
    j["b"] = b_;
    j["W"] = W_;
  }

  void from_json(const json &pars) override {
    std::string name = FieldVal<std::string>(pars, "Name");
    if (name != "RbmSpinReal") {
      throw InvalidInputError(
          "Error while constructing RbmSpinReal from input parameters");
    }

    if (FieldExists(pars, "Nvisible")) {
      nv_ = FieldVal<int>(pars, "Nvisible");
    }
    if (nv_ != hilbert_.Size()) {
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
      a_ = FieldVal<RealVectorType>(pars, "a");
    } else {
      a_.setZero();
    }

    if (FieldExists(pars, "b")) {
      b_ = FieldVal<RealVectorType>(pars, "b");
    } else {
      b_.setZero();
    }
    if (FieldExists(pars, "W")) {
      W_ = FieldVal<RealMatrixType>(pars, "W");
    }
  }

  bool IsHolomorphic() override { return false; }
};

}  // namespace netket

#endif
