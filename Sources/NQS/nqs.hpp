#include <complex>
#include "Machine/rbm_nqs.hpp"
#include "Graph/hypercube.hpp"
#include "Hilbert/spins.hpp"
#include "Optimizer/ada_max.hpp"
#include "Sampler/metropolis_local_gate.hpp"
#include "Supervised/supervised.hpp"
#include "Operator/local_operator.hpp"
#include "Utils/parallel_utils.hpp"
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
    AdaMax& op_;

    // Total number of computational nodes to run on
    int totalnodes_;
    int samplesteps_;
    int gateNo_;

    public:

        NQS(int nqubits, int initialHidden, int sampleSteps)
            : nqubits_(nqubits), g_(*new Hypercube(nqubits,1,false)), samplesteps_(sampleSteps),
            hi_(*new Spin(g_, 0.5)), psi_(*new RbmNQS(std::make_shared<Spin>(hi_), nqubits*(nqubits-1)/2, 0, true, true)),
            sa_(*new MetropolisLocal(psi_)),
            op_(*new AdaMax()), gateNo_(0) {
                VectorType a = getPsi_a();
                VectorType b = getPsi_b();
                MatrixType W = getPsi_W();
                
                a.setZero();
                b.setZero();
                W.setZero();

                setPsiParams(a,b,W);

                for(int j = 0; j < initialHidden; j++) {
                    psi_.addHidden();
                }

                MPI_Comm_size(MPI_COMM_WORLD, &totalnodes_);
            }

        void learnGate(int qubit1, int qubit2, int numSamples, int numIterations, MatrixType gateMatrix) {
            gateNo_++;
            MetropolisLocalGate sampler = MetropolisLocalGate(psi_, gateMatrix, samplesteps_);

            int numSamples_ = int(std::ceil(double(numSamples) / double(totalnodes_)));
            std::vector<Eigen::VectorXd> trainingSamples;
            std::vector<Eigen::VectorXcd> trainingTargets;
            double maxTrainingTarget(0);

            int count00 = 0;
            int count01 = 0;
            int count10 = 0;
            int count11 = 0;

            for(int i = 0; i < numSamples_; i++) {
                sampler.Reset(true);
                sampler.Sweep(qubit1, qubit2);

                trainingSamples.push_back(sampler.Visible());

                Eigen::VectorXcd target(1);
                target(0) = sampler.PsiAfterGate(sampler.Visible(), qubit1, qubit2);
                trainingTargets.push_back(target);

                double targetAbs = abs(target(0));
                if (maxTrainingTarget < targetAbs) {
                    maxTrainingTarget = targetAbs;
                }

                int valueQubit2 = qubit2;

                if(qubit2 == -1) {
                    qubit2 = 0;
                }

                if(sampler.Visible()(qubit1) == 0 && sampler.Visible()(qubit2) == 0) {
                    count00++;
                }

                if(sampler.Visible()(qubit1) == 0 && sampler.Visible()(qubit2) == 1) {
                    count01++;
                }

                if(sampler.Visible()(qubit1) == 1 && sampler.Visible()(qubit2) == 0) {
                    count10++;
                }

                if(sampler.Visible()(qubit1) == 1 && sampler.Visible()(qubit2) == 1) {
                    count11++;
                }

                qubit2 = valueQubit2;
            }


            // in these cases, the gradient factors out and collapses
            if(count00 == numSamples_ || count01 == numSamples_ || count10 == numSamples_ || count11 == numSamples_) {
                // we have to add more samples
                for(int i = 0; i < numSamples_; i++) {
                    sampler.Reset(true);

                    auto sample = sampler.Visible();
                    trainingSamples.push_back(sample);

                    Eigen::VectorXcd target(1);
                    target(0) = sampler.PsiAfterGate(sample, qubit1, qubit2);
                    trainingTargets.push_back(target);

                    double targetAbs = abs(target(0));
                    if (maxTrainingTarget < targetAbs) {
                        maxTrainingTarget = targetAbs;
                    }
                }
            }

            Supervised spvsd = Supervised(psi_, op_, sa_, int(std::ceil(double(trainingSamples.size())/5.0)), trainingSamples, trainingTargets, maxTrainingTarget, "SR");
            spvsd.Run(numIterations, "Overlap_uni", std::to_string(gateNo_));
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

        void applyToffoli(int qubit1, int qubit2, int qubit3, int numSamples = 100, int numIterations = 1000) {
            applyControlledZRotation(qubit2, qubit3, M_PI, numSamples, numIterations);
            applyHadamard(qubit3, numSamples, numIterations);

            applyTDagger(qubit3);

            applyHadamard(qubit3, numSamples, numIterations);
            applyControlledZRotation(qubit1, qubit3, M_PI, numSamples, numIterations);
            applyHadamard(qubit3, numSamples, numIterations);

            applyT(qubit3);

            applyHadamard(qubit3, numSamples, numIterations);
            applyControlledZRotation(qubit2, qubit3, M_PI, numSamples, numIterations);
            applyHadamard(qubit3, numSamples, numIterations);

            applyTDagger(qubit3);

            applyHadamard(qubit3, numSamples, numIterations);
            applyControlledZRotation(qubit1, qubit3, M_PI, numSamples, numIterations);
            applyHadamard(qubit3, numSamples, numIterations);

            applyT(qubit2);
            applyT(qubit3);

            applyHadamard(qubit2, numSamples, numIterations);
            applyControlledZRotation(qubit1, qubit2, M_PI, numSamples, numIterations);
            applyHadamard(qubit2, numSamples, numIterations);

            applyHadamard(qubit3, numSamples, numIterations);

            applyT(qubit1);
            applyTDagger(qubit2);

            applyHadamard(qubit2, numSamples, numIterations);
            applyControlledZRotation(qubit1, qubit2, M_PI, numSamples, numIterations);
            applyHadamard(qubit2, numSamples, numIterations);
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
