import subprocess, os, errno

experiment_name = 'bell-test-qasm'
noctua_user = 'hpc-prf-nqs'
email = 'stubbi@mail.upb.de'

# parameters which should be evaluated
number_of_nodes = [1]
number_of_tasks_per_node = [1]
number_of_omp_threads = [1]

number_of_training_samples = [100]
number_of_training_iterations = [100000]

number_of_initial_hidden_units = [0]
number_of_sample_steps = [0]
number_of_runs = 10


epxperiment_folder = "{pc2pfs}/{noctua_user}/{experiment_name}".format(noctua_user=noctua_user,
                                    pc2pfs=os.environ["PC2PFS"],
                                    experiment_name=experiment_name)


batch_script ="""#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -J {experiment_name}-evaluation
#SBATCH -A {noctua_user}
#SBATCH -p short
#SBATCH -t 00:30:00
#SBATCH --mail-type fail
#SBATCH --mail-user {email}

python $HOME/nqs/scripts/evaluation.py {epxperiment_folder} 2 1 1 1 100 100000 0 0 10 > evaluation_out 2> evaluation_err""".format(
                        epxperiment_folder=epxperiment_folder,
                        experiment_name=experiment_name,
                        noctua_user=noctua_user,
                        email=email)

f = open("{epxperiment_folder}/evaluation.slurm".format(epxperiment_folder=epxperiment_folder),'w')
print >>f, batch_script

bashCommand = "sbatch -D {epxperiment_folder} {epxperiment_folder}/job.slurm".format(epxperiment_folder=epxperiment_folder)
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)