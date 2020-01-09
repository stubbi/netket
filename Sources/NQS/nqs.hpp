#include <complex>
#include "Machine/rbm_nqs.hpp"
#include "Graph/hypercube.hpp"
#include "Hilbert/spins.hpp"
#include "Optimizer/ada_max.hpp"
#include "Sampler/metropolis_local_hadamard.hpp"
#include "Supervised/supervised.hpp"
#include "Operator/local_operator.hpp"
#include "Utils/parallel_utils.hpp"
#include <vector>
#include <iostream>
#include <chrono>


namespace nqs {

class NQS {

    using VectorType = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;
    using MatrixType = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;

    int nqubits_;
    Hypercube& g_;
    Spin& hi_;
    RbmNQS& psi_;
    MetropolisLocal& sa_;
    MetropolisLocalHadamard& saHadamard_;
    AdaMax& op_;

    // Total number of computational nodes to run on
    int totalnodes_;

    public:

        NQS(int nqubits)
            : nqubits_(nqubits), g_(*new Hypercube(nqubits,1,false)),
            hi_(*new Spin(g_, 0.5)), psi_(*new RbmNQS(std::make_shared<Spin>(hi_), 0, 0, true, true)),
            sa_(*new MetropolisLocal(psi_)),
            saHadamard_(*new MetropolisLocalHadamard(psi_)),
            op_(*new AdaMax()) {
                VectorType a = getPsi_a();
                VectorType b = getPsi_b();
                MatrixType W = getPsi_W();
                
                for(int k = 0; k < a.size(); k++) {
                    a(k) = 0.0;
                }

                setPsiParams(a,b,W);

                MPI_Comm_size(MPI_COMM_WORLD, &totalnodes_);
            }

        void applyHadamard(int qubit, int numSamples = 100, int numIterations = 1000) {
            int numSamples_ = int(std::ceil(double(numSamples) / double(totalnodes_)));
            std::vector<Eigen::VectorXd> trainingSamples(numSamples_);
            std::vector<Eigen::VectorXcd> trainingTargets(numSamples_);
            
            int countOne = 0;

            auto start = std::chrono::high_resolution_clock::now();
            for(int i = 0; i < numSamples_; i++) {
                saHadamard_.Reset(true);
                saHadamard_.Sweep(qubit);

                trainingSamples[i] = saHadamard_.Visible();

                Eigen::VectorXcd target(1);
                target(0) = std::log(saHadamard_.PsiValueAfterHadamard(saHadamard_.Visible(), qubit));
                trainingTargets[i] = target;

                //InfoMessage() << saHadamard_.Visible() << " " << target << std::endl;

                if(saHadamard_.Visible()(qubit) == 1) {
                    countOne++;
                }
            }

            // in these cases, the gradient factors out and collapses
            if(countOne == 0 || countOne == numSamples_) {
                // we have to add more samples, say 1%
                for(int i = 0; i < numSamples_/100.0; i++) {
                    saHadamard_.Reset(true);
                    saHadamard_.Sweep(qubit);

                    auto sample = saHadamard_.Visible();
                    sample(qubit) = 1.0 - sample(qubit);
                    trainingSamples[numSamples_ + i] = sample;

                    Eigen::VectorXcd target(1);
                    target(0) = std::log(saHadamard_.PsiValueAfterHadamard(sample, qubit));
                    trainingTargets[numSamples_ + i] = target;

                    //InfoMessage() << sample << " " << target << " " << saHadamard_.PsiValueAfterHadamard(sample, qubit) << std::endl;
                }
            }
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            InfoMessage() << "Elapsed time sampling: " << elapsed.count() << std::endl;

            Supervised spvsd = *new Supervised(psi_, op_, sa_, int(std::ceil(double(trainingSamples.size())/10.0)), trainingSamples, trainingTargets);
            start = std::chrono::high_resolution_clock::now();
            spvsd.Run(numIterations, "Overlap_phi");
            finish = std::chrono::high_resolution_clock::now();
            elapsed = finish - start;
            InfoMessage() << "Elapsed time supervised: " << elapsed.count() << std::endl;
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

        void applyControlledZRotation(int controlQubit, int qubit, double theta) {
            std::complex<double> A_theta = std::acosh(std::exp(std::complex<double>(0, -theta/2.0)));

            psi_.addHidden();

            VectorType a = getPsi_a();
            VectorType b = getPsi_b();
            MatrixType W = getPsi_W();

            W(controlQubit, W.cols()-1) = -2.0 * A_theta;
            W(qubit, W.cols()-1) = 2.0 * A_theta;
        
            a(controlQubit) += std::complex<double>(0, theta/2.0) + A_theta;
            a(qubit) += std::complex<double>(0, theta/2.0) - A_theta;

            setPsiParams(a,b,W);
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
            applyControlledZRotation(qubit2, qubit3, M_PI);
            applyHadamard(qubit3, numSamples, numIterations);

            applyTDagger(qubit3);

            applyHadamard(qubit3, numSamples, numIterations);
            applyControlledZRotation(qubit1, qubit3, M_PI);
            applyHadamard(qubit3, numSamples, numIterations);

            applyT(qubit3);

            applyHadamard(qubit3, numSamples, numIterations);
            applyControlledZRotation(qubit2, qubit3, M_PI);
            applyHadamard(qubit3, numSamples, numIterations);

            applyTDagger(qubit3);

            applyHadamard(qubit3, numSamples, numIterations);
            applyControlledZRotation(qubit1, qubit3, M_PI);
            applyHadamard(qubit3, numSamples, numIterations);

            applyT(qubit2);
            applyT(qubit3);

            applyHadamard(qubit2, numSamples, numIterations);
            applyControlledZRotation(qubit1, qubit2, M_PI);
            applyHadamard(qubit2, numSamples, numIterations);

            applyHadamard(qubit3, numSamples, numIterations);

            applyT(qubit1);
            applyTDagger(qubit2);

            applyHadamard(qubit2, numSamples, numIterations);
            applyControlledZRotation(qubit1, qubit2, M_PI);
            applyHadamard(qubit2, numSamples, numIterations);
        }

        const Eigen::VectorXd& sample() {
            sa_.Reset(true);
            sa_.Sweep();
            return sa_.Visible();
        }

        VectorType getPsiParams() {
            return psi_.GetParameters();
        }

        Complex psi(AbstractMachine::VisibleConstType v) {
            return std::exp(psi_.LogVal(v));
        }

        /*
        *
        * This function is for debug purpose only. It prints the truth table for
        * the current circuit and has an exponential runtime!
        */
        void truthTable() {
            std::map<std::vector<int>, float> occurences;
            InfoMessage() << "Truth table:" << std::endl;
            for(int i = 0; i < 1000; i++) {
               sa_.Reset(true);
               sa_.Sweep();
               std::vector<int> v(sa_.Visible().data(), sa_.Visible().data() + sa_.Visible().size());
               try {
                   occurences.at(v) = occurences.at(v) + 1.0/1000.0;
               } catch (...) {
                   occurences.insert(std::pair<std::vector<int>,float>(v,1.0/1000.0));
               }
            }

            std::vector<int> v(sa_.Visible().data(), sa_.Visible().data() + sa_.Visible().size());
            for(int i = 0; i < std::pow(2.0, sa_.Visible().size()); i++) {
                convertToBinary(i, v);
                std::ostringstream vts;
                std::copy(v.begin(), v.end(), std::ostream_iterator<int>(vts, "")); 
                try {
                    float o = occurences.at(v);
                    InfoMessage() << vts.str() << " " << o << std::endl;
               } catch (...) {
                    InfoMessage() << vts.str() << " 0.0" << std::endl;
               }
            }
            InfoMessage() << std::endl;
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