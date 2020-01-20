import subprocess, os, errno

circuit = 'bell.qc'
experiment_name = 'bell-test-qasm'
noctua_partition = 'short'
max_wall_time = '00:30:00'
email = 'stubbi@mail.upb.de'
noctua_user = 'hpc-prf-nqs'
singularity_image_location = "{pc2pfs}/{noctua_user}/nqs.sif".format(
                        noctua_user=noctua_user,
                        pc2pfs=os.environ["PC2PFS"])

# parameters to be tested
number_of_nodes = [1]
number_of_tasks_per_node = [1]
number_of_omp_threads = [1]

number_of_training_samples = [100]
number_of_training_iterations = [100000]

number_of_initial_hidden_units = [0]
number_of_sample_steps = [0]
number_of_runs = 10


directory = "{pc2pfs}/{noctua_user}/{experiment_name}".format(noctua_user=noctua_user,
                                    pc2pfs=os.environ["PC2PFS"],
                                    experiment_name=experiment_name)
try: os.makedirs(directory)
except OSError, err:
    # Reraise the error unless it's about an already existing directory 
    if err.errno != errno.EEXIST or not os.path.isdir(directory): 
        raise

batch_script ="""#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -J {experiment_name}-exact
#SBATCH -A {noctua_user}
#SBATCH -p {noctua_partition}
#SBATCH -t {max_wall_time}
#SBATCH --mail-type fail
#SBATCH --mail-user {email}

module reset
module load singularity
module load mpi/OpenMPI/3.1.4-GCC-8.3.0
export OMP_NUM_THREADS=1

cp $HOME/nqs/scripts/{circuit} {directory}/in.qc
mpirun -mca pml cm -mca mtl psm2 --report-bindings singularity exec {singularity_image_location} python2.7 $HOME/nqs/scripts/qasm_reader.py 0 0 0 0 exact > out 2> err""".format(
                        experiment_name=experiment_name,
                        noctua_user=noctua_user,
                        noctua_partition=noctua_partition,
                        max_wall_time=max_wall_time,
                        email=email,
                        circuit=circuit,
                        singularity_image_location=singularity_image_location,
                        directory=directory
                    )

f = open("job.slurm",'w')
print >>f, batch_script

bashCommand = "sbatch -D {pc2pfs}/{noctua_user}/{experiment_name} job.slurm".format(noctua_user=noctua_user,
                                    pc2pfs=os.environ["PC2PFS"],
                                    experiment_name=experiment_name)
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

for nodes in number_of_nodes:
    for tasks in number_of_tasks_per_node:
        for threads in number_of_omp_threads:
            for samples in number_of_training_samples:
                for iterations in number_of_training_iterations:
                    for initial_hidden in number_of_initial_hidden_units:
                        for sample_steps in number_of_sample_steps:
                            for run in range(number_of_runs):

                                directory = "{pc2pfs}/{noctua_user}/{experiment_name}/{nodes}nodes/{tasks}tasks/{threads}threads/{samples}samples/{iterations}iterations/{initial_hidden}initialHidden/{sample_steps}sampleSteps/run{run}".format(
                                    noctua_user=noctua_user,
                                    pc2pfs=os.environ["PC2PFS"],
                                    experiment_name=experiment_name,
                                    nodes=nodes,
                                    tasks=tasks,
                                    threads=threads,
                                    samples=samples,
                                    iterations=iterations,
                                    initial_hidden=initial_hidden,
                                    sample_steps=sample_steps,
                                    run=run
                                )

                                try: os.makedirs(directory)
                                except OSError, err:
                                    # Reraise the error unless it's about an already existing directory 
                                    if err.errno != errno.EEXIST or not os.path.isdir(directory): 
                                        raise

                                batch_script ="""#!/bin/bash
#SBATCH -N {nodes}
#SBATCH --ntasks-per-node={tasks}
#SBATCH -J {experiment_name}-{nodes}nodes-{tasks}tasks-{threads}threads-{samples}samples-{iterations}iterations-run{run}
#SBATCH -A {noctua_user}
#SBATCH -p {noctua_partition}
#SBATCH -t {max_wall_time}
#SBATCH --mail-type fail
#SBATCH --mail-user {email}

module reset
module load singularity
module load mpi/OpenMPI/3.1.4-GCC-8.3.0
export OMP_NUM_THREADS={threads}

cp $HOME/nqs/scripts/{circuit} {directory}/in.qc
mpirun -mca pml cm -mca mtl psm2 --report-bindings singularity exec {singularity_image_location} python2.7 $HOME/nqs/scripts/qasm_reader.py {samples} {iterations} {initial_hidden} {sample_steps} nqs > out 2> err""".format(
                        nodes=nodes,
                        experiment_name=experiment_name,
                        tasks=tasks,
                        threads=threads,
                        samples=samples,
                        initial_hidden = initial_hidden,
                        sample_steps = sample_steps,
                        iterations=iterations,
                        noctua_user=noctua_user,
                        noctua_partition=noctua_partition,
                        max_wall_time=max_wall_time,
                        email=email,
                        circuit=circuit,
                        singularity_image_location=singularity_image_location,
                        run=run,
                        pc2pfs=os.environ["PC2PFS"],
                        directory=directory
                    )

                                f = open("job.slurm",'w')
                                print >>f, batch_script

                                bashCommand = "sbatch -D {directory} job.slurm".format(directory=directory)
                                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)                

                                print "started job {experiment_name} for {nodes}nodes {tasks}tasks {threads}threads {samples}samples {iterations}iterations {initial_hidden}initialHidden {sample_steps}sampleSteps run {run}".format(
                                    experiment_name=experiment_name,
                                    nodes=nodes,
                                    tasks=tasks,
                                    threads=threads,
                                    samples=samples,
                                    iterations=iterations,
                                    initial_hidden=initial_hidden,
                                    sample_steps=sample_steps,
                                    run=run
                                )