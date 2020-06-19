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

#ifndef NQS_SUPERVISED_CC
#define NQS_SUPERVISED_CC

#include <bitset>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "Machine/machine.hpp"
#include "Optimizer/optimizer.hpp"
#include "Output/json_output_writer.hpp"
#include "Stats/stats.hpp"
#include "Utils/parallel_utils.hpp"
#include "Utils/random_utils.hpp"

namespace nqs {

class Supervised {
  using MatrixT = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;

  AbstractMachine &psi_;
  AbstractOptimizer &opt_;
  MetropolisLocal& sa_;

  SR sr_;
  bool dosr_;

  // Total batchsize
  int batchsize_;
  // Batchsize per node
  int batchsize_node_;

  // Total number of computational nodes to run on
  int totalnodes_;
  int mynode_;

  // Number of parameters of the machine
  int npar_;
  // Stores the gradient of loss w.r.t. the machine parameters
  Eigen::VectorXcd grad_;
  Eigen::VectorXcd grad_part_1_;
  Eigen::VectorXcd grad_part_2_;
  Complex grad_part_3_;
  Complex grad_num_1_;
  Complex grad_num_2_;
  Complex grad_num_3_;

  // Training samples and targets
  std::vector<Eigen::VectorXd> trainingSamples_;
  std::vector<Eigen::VectorXcd> trainingTargets_;
  // Test samples and targets
  std::vector<Eigen::VectorXd> testSamples_;
  std::vector<Eigen::VectorXcd> testTargets_;

  // All loss function is real
  double training_loss_log_overlap_;
  double training_loss_mse_;
  double training_loss_mse_log_;
  double test_loss_log_overlap_;
  double test_loss_mse_;
  double test_loss_mse_log_;

  std::uniform_int_distribution<int> distribution_uni_;
  std::discrete_distribution<> distribution_phi_;

  MatrixT Ok_;

 protected:
  // Random number generator with correct seeding for parallel processes
  default_random_engine &GetRandomEngine() { return engine_.Get(); }

 private:
  DistributedRandomEngine engine_;

 public:
  Supervised(AbstractMachine &psi,
             MetropolisLocal &sa,
             AbstractOptimizer &opt,
             int batchsize,
             std::vector<Eigen::VectorXd> trainingSamples,
             std::vector<Eigen::VectorXcd> trainingTargets,
             std::vector<Eigen::VectorXd> testSamples,
             std::vector<Eigen::VectorXcd> testTargets,
             bool sr,
             double diag_shift = 0.01,
             bool use_iterative = false, bool use_cholesky = true)
      : psi_(psi),
        sa_(sa),
        opt_(opt),
        trainingSamples_(trainingSamples),
        trainingTargets_(trainingTargets),
        testSamples_(testSamples),
        testTargets_(testTargets),
        dosr_(sr) {
    npar_ = psi_.Npar();


    grad_.resize(npar_);
    grad_part_1_.resize(npar_);
    grad_part_2_.resize(npar_);

    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode_);

    batchsize_ = batchsize;
    batchsize_node_ = int(std::ceil(double(batchsize_) / double(totalnodes_)));

    // Initialize a uniform distribution to draw training samples from
    distribution_uni_ =
        std::uniform_int_distribution<int>(0, trainingSamples_.size() - 1);

    opt_.Init(npar_, psi_.IsHolomorphic());

    if (dosr_) {
      setSrParameters(diag_shift, use_iterative, use_cholesky);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }

  /// Computes the gradient estimate of the derivative of negative log
  /// of wavefunction overlap, with index i sampled from unifrom[1, N]
  void DerLogOverlap_uni(std::vector<Eigen::VectorXd> &batchSamples,
                         std::vector<Eigen::VectorXcd> &batchTargets) {
    // ... and zero them out
    grad_.setZero(psi_.Npar());
    grad_part_1_.setZero(psi_.Npar());
    grad_part_2_.setZero(psi_.Npar());
    grad_num_1_ = 0;
    grad_num_2_ = 0;

    double max_log_psi = 0;
    /// [TODO] avoid going through psi twice.
    for (int i = 0; i < batchsize_node_; i++) {
      Complex value(psi_.LogVal(batchSamples[i]));
      if (max_log_psi < value.real()) {
        max_log_psi = value.real();
      }
    }

    Ok_.resize(batchsize_node_, psi_.Npar());

    // For each sample in the batch
    for (int i = 0; i < batchsize_node_; i++) {
      // Extract log(config)
      Eigen::VectorXd sample(batchSamples[i]);
      // And the corresponding target
      Eigen::VectorXcd target(batchTargets[i]);
      Complex t(target[0].real(), target[0].imag());
      // Undo log
      // already in exp form t = exp(t);

      Complex value(psi_.LogVal(sample));
      // Undo Log
      value = value;// - max_log_psi;
      value = exp(value);

      // Compute derivative of log
      auto der = psi_.DerLog(sample);
      Ok_.row(i) = der;

      der = der.conjugate();

      grad_part_1_ = grad_part_1_ + der * pow(abs(value), 2);
      grad_num_1_ = grad_num_1_ + pow(abs(value), 2);
      grad_part_2_ = grad_part_2_ + der * std::conj(value) * t;
      grad_num_2_ = grad_num_2_ + std::conj(value) * t;
    }

    SumOnNodes(grad_part_1_);
    SumOnNodes(grad_num_1_);
    SumOnNodes(grad_part_2_);
    SumOnNodes(grad_num_2_);
    /// No need to divide by totalnodes_
    grad_ = grad_part_1_ / grad_num_1_ - grad_part_2_ / grad_num_2_;
  }


  /// Runs the supervised learning on the training samples and targets
  void Run(int n_iter, bool early_stopping,
           const std::string &output_prefix = "output",
           int save_params_every = 500) {
    assert(n_iter > 0);
    assert(save_params_every > 0);

    // for early stopping
    double last_test_log_overlap = std::numeric_limits<double>::infinity();

    /// Writer to the output
    /// This optional will contain a value iff the MPI rank is 0.
    nonstd::optional<JsonOutputWriter> writer;
    if (mynode_ == 0) {
      /// Initializes the output log and wavefunction files
      writer.emplace(output_prefix + ".log", output_prefix + ".wf",
                     save_params_every);
    }

    opt_.Reset();

    std::vector<Eigen::VectorXd> batchSamples(batchsize_node_);
    std::vector<Eigen::VectorXcd> batchTargets(batchsize_node_);

    for (int i = 0; i < n_iter; i++) {

      int index;

      // Randomly select a batch of training data
      for (int k = 0; k < batchsize_node_; k++) {
        // Draw from the distribution using the nqs random number generator
        index = distribution_uni_(this->GetRandomEngine());
        batchSamples[k] = trainingSamples_[index];
        batchTargets[k] = trainingTargets_[index];
      }

      DerLogOverlap_uni(batchSamples, batchTargets);
      UpdateParameters();
      ComputeLosses();
      // writer.has_value() iff the MPI rank is 0, so the output is only
      // written once
      
      if (writer.has_value()) {
        json out_data;
        out_data["training_log_overlap"] = training_loss_log_overlap_;
        out_data["training_mse"] = training_loss_mse_;
        out_data["training_mse_log"] = training_loss_mse_log_;
        out_data["test_log_overlap"] = test_loss_log_overlap_;
        out_data["test_mse"] = test_loss_mse_;
        out_data["test_mse_log"] = test_loss_mse_log_;

        writer->WriteLog(i, out_data);
        writer->WriteState(i, psi_);
      }

      // early stopping
      if (early_stopping && last_test_log_overlap < test_loss_log_overlap_) {
        break;
      }

      last_test_log_overlap = test_loss_log_overlap_;
      
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  /// Updates the machine parameters with the current gradient
  void UpdateParameters() {
    auto pars = psi_.GetParameters();

    Eigen::VectorXcd deltap(npar_);

    if (dosr_) {
      sr_.ComputeUpdate(Ok_, grad_, deltap);
    } else {
      deltap = grad_;
    }

    opt_.Update(deltap, pars);
    SendToAll(pars);

    psi_.SetParameters(pars);
  }

  void ComputeLosses() {
    ComputeMse(trainingSamples_, trainingTargets_, training_loss_mse_, training_loss_mse_log_);
    ComputeMse(testSamples_, testTargets_, test_loss_mse_, test_loss_mse_log_);
    ComputeLogOverlap(trainingSamples_, trainingTargets_, training_loss_log_overlap_);
    ComputeLogOverlap(testSamples_, testTargets_, test_loss_log_overlap_);
  }

  /// Compute the MSE of psi and the MSE of log(psi)
  /// for monitoring the convergence.
  void ComputeMse(std::vector<Eigen::VectorXd> samples, std::vector<Eigen::VectorXcd> targets, double& loss_mse, double& loss_mse_log) {
    const int numSamples = samples.size();

    double mse_log = 0.0;
    double mse = 0.0;
    for (int i = 0; i < numSamples; i++) {
      Eigen::VectorXd sample = samples[i];
      Complex value(psi_.LogVal(sample));

      Eigen::VectorXcd target = targets[i];
      Complex t(target[0].real(), target[0].imag());

      mse_log += 0.5 * std::norm(value - t);
      mse += 0.5 * std::norm(exp(value) - exp(t));
    }

    loss_mse = mse / numSamples;
    loss_mse_log = mse_log / numSamples;
  }

  double GetTrainingMse() const { return training_loss_mse_; }
  double GetTestMse() const { return test_loss_mse_; }

  double GetTrainingMseLog() const { return training_loss_mse_log_; }
  double GetTestMseLog() const { return test_loss_mse_log_; }

  void ComputeLogOverlap(std::vector<Eigen::VectorXd> samples, std::vector<Eigen::VectorXcd> targets, double& loss_log_overlap) {
    const int numSamples = samples.size();

    // Allocate vectors for storing the derivatives ...
    Complex num1(0.0, 0.0);
    Complex num2(0.0, 0.0);
    Complex num3(0.0, 0.0);
    Complex num4(0.0, 0.0);

    Eigen::VectorXcd logpsi(numSamples);
    double max_log_psi = -std::numeric_limits<double>::infinity();

    for (int i = 0; i < numSamples; i++) {
      logpsi(i) = psi_.LogVal(samples[i]);
      if (std::real(logpsi(i)) > max_log_psi) {
        max_log_psi = std::real(logpsi(i));
      }
    }

    for (int i = 0; i < numSamples; i++) {
      // Extract log(config)
      Eigen::VectorXd sample(samples[i]);
      // And the corresponding target
      Eigen::VectorXcd target(targets[i]);

      // Cast value and target to Complex and undo logs
      Complex value(psi_.LogVal(sample));
      value = exp(value - max_log_psi);
      Complex t(target[0].real(), target[0].imag());
      // setting t to exp already in NQS.hpp t = exp(t);

      num1 += std::conj(value) * t;
      num2 += value * std::conj(t);
      num3 += std::norm(value);
      num4 += std::norm(t);
    }

    Complex complex_log_overlap_ =
        -(log(num1) + log(num2) - log(num3) - log(num4));
    assert(std::abs(complex_log_overlap_.imag()) < 1e-8);
    loss_log_overlap = complex_log_overlap_.real();
  }

  double GetTrainingLogOverlap() const { return training_loss_log_overlap_; }
  double GetTestLogOverlap() const { return test_loss_log_overlap_; }

  void setSrParameters(double diag_shift = 0.01, bool use_iterative = false,
                       bool use_cholesky = true) {
    sr_.setParameters(diag_shift, use_iterative, use_cholesky,
                      psi_.IsHolomorphic());
  }
};

}  // namespace nqs

#endif
