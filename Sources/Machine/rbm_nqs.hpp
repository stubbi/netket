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

#ifndef NETKET_RBM_NQS_HPP
#define NETKET_RBM_NQS_HPP

#include <cmath>

#include "Machine/abstract_machine.hpp"

namespace netket {

/** Restricted Boltzmann machine class with spin 1/2 hidden units.
 *
 */
class RbmNQS : public AbstractMachine {
  // number of visible units
  int nv_;

  // number of hidden units
  int nh_;

  // number of parameters
  int npar_;

  // weights
  MatrixType W_;

  // visible units bias
  VectorType a_;

  // hidden units bias
  VectorType b_;

  VectorType thetas_;
  VectorType lnthetas_;
  VectorType thetasnew_;
  VectorType lnthetasnew_;

  bool usea_;
  bool useb_;

 public:
  RbmNQS(std::shared_ptr<const AbstractHilbert> hilbert, int nhidden = 0,
          int alpha = 0, bool usea = true, bool useb = true);

  int Nvisible() const override;
  int Npar() const override;
  int Nhidden() const noexcept { return nh_; }
  void addHidden();

  void InitRandomPars(int seed, double sigma) override;
  void InitLookup(VisibleConstType v, LookupType &lt) override;
  void UpdateLookup(VisibleConstType v, const std::vector<int> &tochange,
                    const std::vector<double> &newconf,
                    LookupType &lt) override;
  VectorType DerLog(VisibleConstType v) override;
  VectorType DerLog(VisibleConstType v, const LookupType &lt) override;
  Complex LogVal(VisibleConstType v) override;
  Complex LogVal(VisibleConstType v, const LookupType &lt) override;

  VectorType GetParameters() override;
  void SetParameters(VectorConstRefType pars) override;

  // Difference between logarithms of values, when one or more visible variables
  // are being flipped
  VectorType LogValDiff(
      VisibleConstType v, const std::vector<std::vector<int>> &tochange,
      const std::vector<std::vector<double>> &newconf) override;

  // Difference between logarithms of values, when one or more visible variables
  // are being flipped Version using pre-computed look-up tables for efficiency
  // on a small number of spin flips
  Complex LogValDiff(VisibleConstType v, const std::vector<int> &tochange,
                     const std::vector<double> &newconf,
                     const LookupType &lt) override;

  void Save(const std::string &filename) const override;
  void Load(const std::string &filename) override;

  bool IsHolomorphic() const noexcept override;


  static void sigmoid(VectorConstRefType x, VectorType &y) {
    assert(y.size() >= x.size());
    y = 1.0 / (1.0 + Eigen::exp(-1.0 * x.array()));
  }

  static void sigmoid(RealVectorConstRefType x, RealVectorType &y) {
    assert(y.size() >= x.size());
    y = 1.0 / (1.0 + Eigen::exp(-1.0 * x.array()));
  }

  static double softplus(double x) {
    const double xp = std::abs(x);
    if (xp <= 12.) {
      return std::log(1.0 + std::exp(xp));
    } else {
      const static double log2v = std::log(2.);
      return xp - log2v;
    }
  }

  // softplus(x) for std::complex argument
  // the modulus is computed by means of the previously defined function
  // for real argument
  static Complex softplus(Complex x) {
    return std::log(1.0 + std::exp(x));
  }

  static void softplus(VectorConstRefType x, VectorType &y) {
    assert(y.size() >= x.size());
    for (int i = 0; i < x.size(); i++) {
      y(i) = softplus(x(i));
    }
  }

  static void softplus(RealVectorConstRefType x, RealVectorType &y) {
    assert(y.size() >= x.size());
    for (int i = 0; i < x.size(); i++) {
      y(i) = softplus(x(i));
    }
  }


 private:
  inline void Init();
};

}  // namespace netket

#endif
