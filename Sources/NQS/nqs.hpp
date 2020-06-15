#include <complex>
#include "Machine/rbm_nqs.hpp"
#include "Graph/hypercube.hpp"
#include "Hilbert/spins.hpp"
#include "Sampler/metropolis_local_gate.hpp"
#include "Supervised/supervised.hpp"
#include "Operator/local_operator.hpp"
#include "Optimizer/abstract_optimizer.hpp"
#include "Utils/parallel_utils.hpp"
#include "Utils/random_utils.hpp"
#include <vector>
#include <iostream>
#include <chrono>
#include <math.h>
#include <string>


namespace nqs {

class NQS {

    using VectorType = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;
    using MatrixType = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;

    int nqubits_;
    Hypercube& g_;
    Spin& hi_;
    RbmNQS& psi_;
    MetropolisLocal& sa_;
    std::string& optimizer_;
    std::shared_ptr<AbstractOptimizer> opt_;

    int totalnodes_;
    int samplesteps_;
    int gateNo_;
    int randomRestarts_;
    bool earlyStopping_;
    bool sr;

    public:

  NQS(int nqubits, int initialHidden, int sampleSteps, int randomRestarts, bool earlyStopping, std::string &optimizer)
            : nqubits_(nqubits), g_(*new Hypercube(nqubits,1,false)), samplesteps_(sampleSteps),
            hi_(*new Spin(g_, 0.5)), psi_(*new RbmNQS(std::make_shared<Spin>(hi_), initialHidden, 0, true, true)),
              sa_(*new MetropolisLocal(psi_)), optimizer_(optimizer), gateNo_(0), randomRestarts_(randomRestarts), earlyStopping_(earlyStopping), sr(false) {
                VectorType a = getPsi_a();
                VectorType b = getPsi_b();
                MatrixType W = getPsi_W();
                
                a.setZero();
                b.setZero();
                W.setZero();

                setPsiParams(a,b,W);

                // InfoMessage() << "optimizer_ " << optimizer_ << std::endl;

                if (optimizer_ == "AdaDelta") {
                  opt_ = std::make_shared<AdaDelta>();
                  InfoMessage() << "initialized AdaDelta()" << std::endl;
                } else if (optimizer_ == "AdaGrad") {
                  opt_ = std::make_shared<AdaGrad>();
                  InfoMessage() << "initialized AdaGrad()" << std::endl;
                } else if (optimizer_ == "AdaMax") {
                  opt_ = std::make_shared<AdaMax>();
                  InfoMessage() << "initialized AdaMax()" << std::endl;
                } else if (optimizer_ == "AMSGrad") {
                  opt_ = std::make_shared<AMSGrad>();
                  InfoMessage() << "initialized AMSGrad()" << std::endl;
                } else if (optimizer_ == "Momentum") {
                  opt_ = std::make_shared<Momentum>();
                  InfoMessage() << "initialized Momentum()" << std::endl;
                } else if (optimizer_ == "RMSProp") {
                  opt_ = std::make_shared<RMSProp>();
                  InfoMessage() << "initialized RMSProb()" << std::endl;
                } else if (optimizer_ == "Sgd") {
                  opt_ = std::make_shared<Sgd>();
                  InfoMessage() << "initialized Sgd()" << std::endl;
                } else {
                  sr = true;
                  // avoid deferencing nullptr later on
                  opt_ = std::make_shared<AdaDelta>();
                  InfoMessage() << "initialized Sgd() [for stochastic reconfiguration]" << std::endl;
                }

                MPI_Comm_size(MPI_COMM_WORLD, &totalnodes_);
        }


        void learnGate(int qubit1, int qubit2, int numSamples, int numIterations, MatrixType gateMatrix) {
            gateNo_++;
            MetropolisLocalGate sampler = MetropolisLocalGate(psi_, gateMatrix, samplesteps_);

            int numSamplesNode = int(std::ceil(double(numSamples) / double(totalnodes_)));
            int batchSize = numSamplesNode < 10 ? numSamplesNode : 10; //https://arxiv.org/pdf/1606.02228.pdf 3.7 => prefer smaller batch sizes

            std::vector<Eigen::VectorXd> trainingSamples;
            std::vector<Eigen::VectorXcd> trainingTargets;

            // independently sampled test set. Same size as training set
            std::vector<Eigen::VectorXd> testSamples;
            std::vector<Eigen::VectorXcd> testTargets;

            generateSamples(qubit1, qubit2, numSamplesNode, sampler, trainingSamples, trainingTargets);
            generateSamples(qubit1, qubit2, numSamplesNode, sampler, testSamples, testTargets);

            InfoMessage() << "trainingSamples.size()" << trainingSamples.size() << std::endl;
            InfoMessage() << "trainingTargets.size()" << trainingTargets.size() << std::endl;
            InfoMessage() << "testSamples.size()" << testSamples.size() << std::endl;
            InfoMessage() << "testTargets.size()" << testTargets.size() << std::endl;

            InfoMessage() << "randomRstearts_ " << randomRestarts_ << std::endl;

            if (randomRestarts_ <= 0) {
              Supervised spvsd = Supervised(psi_, sa_, *opt_, batchSize, trainingSamples, trainingTargets, testSamples, testTargets, sr);
              InfoMessage() << "initialized spvsd()" << std::endl;
              spvsd.Run(numIterations, earlyStopping_, std::to_string(gateNo_));
              InfoMessage() << "ran spvsd" << std::endl;
            } else {
              savePsiParams("supervised_gate_" + std::to_string(gateNo_) + "_random_restarts_" + std::to_string(0) + ".json");

              double minLogOverlap = std::numeric_limits<double>::infinity();
              int rMinLogOverlap = 0;

              for (int r = 0; r < randomRestarts_; r++) {
                loadPsiParams("supervised_gate_" + std::to_string(gateNo_) + "_random_restarts_" + std::to_string(0) + ".json");

                RbmNQS::VectorType pars = psi_.GetParameters();
                nqs::RandomGaussianPermutation(pars, r, 0.01);
                psi_.SetParameters(pars);

                Supervised spvsd = Supervised(psi_, sa_, *opt_, batchSize, trainingSamples, trainingTargets, testSamples, testTargets, sr);
                InfoMessage() << "initialized spvsd()" << std::endl;
                spvsd.Run(numIterations, earlyStopping_, std::to_string(gateNo_));
                InfoMessage() << "ran spvsd" << std::endl;

                if (spvsd.GetTestLogOverlap() < minLogOverlap) {
                  minLogOverlap = spvsd.GetTestLogOverlap();
                  rMinLogOverlap = r;
                }

                savePsiParams("supervised_gate_" + std::to_string(gateNo_) + "_random_restarts_" + std::to_string(r) + ".json");
              }

              loadPsiParams("supervised_gate_" + std::to_string(gateNo_) + "_random_restarts_" + std::to_string(rMinLogOverlap) + ".json");
            }

        }


        void applyHadamard(int qubit, int numSamples, int numIterations) {
            MatrixType H(2,2);
            H << 1,1,1,-1;
            H = 1.0/sqrt(2.0) * H;
            learnGate(qubit, -1, numSamples, numIterations, H);
        }

        void applySqrtX(int qubit, int numSamples, int numIterations) {
            MatrixType sqrtX(2,2);
            sqrtX << std::complex<double>(1, 1), std::complex<double>(1, -1), std::complex<double>(1, -1), std::complex<double>(1, 1);
            learnGate(qubit, -1, numSamples, numIterations, 0.5*sqrtX);
        }
        
        void applySqrtY(int qubit, int numSamples, int numIterations) {
            MatrixType sqrtY(2,2);
            sqrtY << std::complex<double>(1, 1), std::complex<double>(-1, -1), std::complex<double>(1, 1), std::complex<double>(1, 1);
            learnGate(qubit, -1, numSamples, numIterations, 0.5*sqrtY);
        }

        void applyPauliX(int qubit){
            VectorType a = getPsi_a();
            VectorType b = getPsi_b();
            MatrixType W = getPsi_W();
            
            a(qubit) = -a(qubit);
            for(int k = 0; k < b.size(); k++) {
                b(k) += W(qubit, k);
            }
            W.row(qubit) = -W.row(qubit);

            setPsiParams(a,b,W);
        }

        void applyPauliY(int qubit){
            VectorType a = getPsi_a();
            VectorType b = getPsi_b();
            MatrixType W = getPsi_W();
            
            a(qubit) = -a(qubit) + (0, M_PI/2.0);
            for(int k = 0; k < b.size(); k++) {
                b(k) += W(qubit, k);
            }
            W.row(qubit) = -W.row(qubit);

            setPsiParams(a,b,W);
        }

        void applyPauliZ(int qubit){
            VectorType a = getPsi_a();
            VectorType b = getPsi_b();
            MatrixType W = getPsi_W();
            
            a(qubit) = a(qubit) + std::complex<double>(0, M_PI);

            setPsiParams(a,b,W);
        }

        void applySingleZRotation(int qubit, double theta) {
            VectorType a = getPsi_a();
            VectorType b = getPsi_b();
            MatrixType W = getPsi_W();
            
            a(qubit) = a(qubit) + std::complex<double>(0, theta);

            setPsiParams(a,b,W);
        }

        void applyControlledZRotation(int controlQubit, int qubit, double theta, int numSamples, int numIterations) {
            MatrixType cZ(4,4);
            cZ << 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,std::exp(std::complex<double>(0, theta));
            learnGate(controlQubit, qubit, numSamples, numIterations, cZ);
        }

        /**
         * T gate as in Nielsen and Chuang. 
         * 
         * [ 1 , 0 ; 0 , exp(i*pi/4.0) ]
        **/
        void applyT(int qubit) {
            applySingleZRotation(qubit, M_PI/4.0);
        }

        void applyTDagger(int qubit) {
            applySingleZRotation(qubit, -M_PI/4.0);
        }

        const Eigen::VectorXd& sample() {
            sa_.Reset(true);
            sa_.Sweep();
            return sa_.Visible();
        }

        void savePsiParams(const std::string &filename) {
            psi_.Save(filename);
        }

        void loadPsiParams(const std::string &filename) {
            psi_.Load(filename);
        }

        Complex psi(AbstractMachine::VisibleConstType v) {
            return std::exp(psi_.LogVal(v));
        }

    private:

        void generateSamples(int qubit1, int qubit2, int numSamples, MetropolisLocalGate& sampler, std::vector<Eigen::VectorXd>& samples, std::vector<Eigen::VectorXcd>& targets) {

            for(int i = 0; i < numSamples; i++) {
                sampler.Reset(true);
                sampler.Sweep(qubit1, qubit2);

                samples.push_back(sampler.Visible());

                Eigen::VectorXcd target(1);
                target(0) = sampler.PsiAfterGate(sampler.Visible(), qubit1, qubit2);
                targets.push_back(target);

            }
        }

        void convertToBinary(unsigned int n, std::vector<int> &v) {
            int remainder, i = 1, step = 0;

            while (n!=0) {
                v[step++] =n%2;
                n /= 2;
            }
            while(step < v.size()) {
                v[step++] = 0;
            }
        }
     
        VectorType getPsi_a() {
            RbmNQS::VectorType pars = psi_.GetParameters();
            //always "use_a" & "use_b"
            return pars.head(psi_.Nvisible());
        }

        VectorType getPsi_b() {
            RbmNQS::VectorType pars = psi_.GetParameters();
            //always "use_a" & "use_b"
            return pars.segment(psi_.Nvisible(), psi_.Nhidden());
        }

        MatrixType getPsi_W() {
            RbmNQS::VectorType pars = psi_.GetParameters();
            VectorType Wpars = pars.tail(psi_.Nvisible() * psi_.Nhidden());
            return Eigen::Map<MatrixType>(Wpars.data(), psi_.Nvisible(), psi_.Nhidden());
        }

        void setPsiParams(RbmNQS::VectorType a,
                            RbmNQS::VectorType b,
                            RbmNQS::MatrixType W) {
            VectorType pars(psi_.Npar());
            pars.head(psi_.Nvisible()) = a;
            pars.segment(psi_.Nvisible(), psi_.Nhidden()) = b;
            pars.tail(psi_.Nvisible() * psi_.Nhidden()) = Eigen::Map<VectorType>(W.data(), psi_.Nvisible() * psi_.Nhidden());
            psi_.SetParameters(pars);
        }
        
    };

}  // namespace nqs
