#include <complex>
#include "Machine/rbm_nqs.hpp"
#include "Graph/hypercube.hpp"
#include "Hilbert/spins.hpp"
#include "Optimizer/ada_max.hpp"
#include "Sampler/metropolis_local_hadamard.hpp"
#include "Supervised/supervised.hpp"
#include "Operator/local_operator.hpp"
#include "Unsupervised/quantum_state_reconstruction.hpp"
#include <vector>


namespace netket {

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
            }

        void applyHadamard(int qubit, int numSamples = 100, int numIterations = 1000) {
            std::vector<Eigen::VectorXd> trainingSamples;
            std::vector<Eigen::VectorXcd> trainingTargets;

            int countOne = 0;

            for(int i = 0; i < numSamples; i++) {
                saHadamard_.Reset(true);
                saHadamard_.Sweep(qubit);

                trainingSamples.push_back(saHadamard_.Visible());

                Eigen::VectorXcd target(1);
                target(0) = std::log(saHadamard_.PsiValueAfterHadamard(saHadamard_.Visible(), qubit));

                trainingTargets.push_back(target);

                if(saHadamard_.Visible()(qubit) == 1) {
                    countOne++;
                }

                InfoMessage() << saHadamard_.Visible() << " " << target << std::endl;
            }

            /*
            // in these cases, the gradient factors out and collapses
            if(countOne == 0 || countOne == numSamples) {
                // we have to add more samples, say 1%
                for(int i = 0; i < numSamples; i++) {
                    saHadamard_.Reset(true);
                    saHadamard_.Sweep(qubit);

                    auto sample = saHadamard_.Visible();
                    sample(qubit) = 1.0 - sample(qubit);

                    trainingSamples.push_back(sample);

                    Eigen::VectorXcd target(1);
                    target(0) = std::log(saHadamard_.PsiValueAfterHadamard(sample, qubit));

                    trainingTargets.push_back(target);

                    InfoMessage() << sample << " " << target << std::endl;
                }

            }
            */

            InfoMessage() << "ones " << countOne << std::endl;
            InfoMessage() << "#samples " << trainingSamples.size() << std::endl;

            LocalOperator* o = new LocalOperator(std::make_shared<Spin>(hi_), 1.0);
            std::vector<AbstractOperator *> rotations;
            rotations.push_back(o);

            std::vector<int> trainingBases;
            trainingBases.push_back(0);
    
            QuantumStateReconstruction qsr = *new QuantumStateReconstruction(sa_, op_,
                                                trainingSamples.size()/10.0, trainingSamples.size(),
                                                rotations, trainingSamples, trainingBases);

            qsr.Run("out", numIterations);

            //Supervised spvsd = *new Supervised(psi_, op_, trainingSamples.size(), trainingSamples, trainingTargets);
            //spvsd.Run(numIterations, "Overlap_uni");
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

        void applySingleZRotation(int qubit, double theta){
            VectorType a = getPsi_a();
            VectorType b = getPsi_b();
            MatrixType W = getPsi_W();
            
            a(qubit) = a(qubit) + std::complex<double>(0, theta);

            setPsiParams(a,b,W);
        }

        void applyControlledZRotation(int controlQubit, int qubit, double theta){
            std::complex<double> A_theta = std::acosh(std::exp(std::complex<double>(0, -theta/2.0)));

            InfoMessage() << "A_theta " << A_theta << std::endl;

            psi_.addHidden();

            VectorType a = getPsi_a();
            VectorType b = getPsi_b();
            MatrixType W = getPsi_W();

            W(qubit, W.cols()-1) = -2.0 * A_theta;
            W(controlQubit, W.cols()-1) = 2.0 * A_theta;
        
            a(qubit) += std::complex<double>(0, theta/2.0) + A_theta;
            a(controlQubit) += std::complex<double>(0, theta/2.0) - A_theta;

            setPsiParams(a,b,W);
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


    private:
     
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

}  // namespace netket