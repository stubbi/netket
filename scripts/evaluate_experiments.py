import subprocess, os, errno

experiment_name = 'new-plotting-test'
noctua_user = 'hpc-prf-nqs'
email = 'stubbi@mail.upb.de'

# parameters to be tested
number_of_qubits = [5, 7, 10]
number_of_cycles = [5, 7, 10]
number_of_circuits = 10 #number of random circuits with same number of qubits and cycles

number_of_nodes = [1]
number_of_tasks_per_node = [1]
number_of_omp_threads = [1]

number_of_training_samples = [100,300,500,700,1000]
number_of_training_iterations = [10000,30000,50000,70000,100000]

number_of_initial_hidden_units = [0]
number_of_sample_steps = [4,6,8,11,13]
number_of_runs = 10 #number of runs for a specific circuit


epxperiment_folder = "{pc2pfs}/{noctua_user}/{experiment_name}".format(noctua_user=noctua_user,
                                    pc2pfs=os.environ["PC2PFS"],
                                    experiment_name=experiment_name)


batch_script ="""#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -J {experiment_name}-evaluation
#SBATCH -A {noctua_user}
#SBATCH -p short
#SBATCH -t 00:10:00
#SBATCH --mail-type fail
#SBATCH --mail-user {email}

module reset
module load vis/matplotlib
python $HOME/nqs/scripts/evaluation.py {epxperiment_folder} {number_of_qubits} {number_of_cycles} {number_of_circuits} {listOMPNodes} {listOMPTasks} {listOMPThreads} {listSamples} {listIterations} {listInitialHidden} {listSampleSteps} {numRuns} > evaluation_out 2> evaluation_err""".format(
                        epxperiment_folder=epxperiment_folder,
                        experiment_name=experiment_name,
                        noctua_user=noctua_user,
                        email=email,
                        number_of_qubits=','.join(map(str, number_of_qubits)),
                        number_of_cycles=','.join(map(str, number_of_cycles)),
                        number_of_circuits=number_of_circuits,
                        listOMPNodes=','.join(map(str, number_of_nodes)),
                        listOMPTasks=','.join(map(str, number_of_tasks_per_node)),
                        listOMPThreads=','.join(map(str, number_of_omp_threads)),
                        listSamples=','.join(map(str, number_of_training_samples)),
                        listIterations=','.join(map(str, number_of_training_iterations)),
                        listInitialHidden=','.join(map(str, number_of_initial_hidden_units)),
                        listSampleSteps=','.join(map(str, number_of_sample_steps)),
                        numRuns=number_of_runs
                        )

f = open("{epxperiment_folder}/evaluation.slurm".format(epxperiment_folder=epxperiment_folder),'w')
print >>f, batch_script

bashCommand = "sbatch -D {epxperiment_folder} {epxperiment_folder}/evaluation.slurm".format(epxperiment_folder=epxperiment_folder)
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)