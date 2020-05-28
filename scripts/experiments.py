import subprocess, os, errno, time
from experiments_settings import number_of_qubits, number_of_cycles, number_of_circuits, number_of_nodes, number_of_tasks_per_node, number_of_omp_threads, number_of_training_samples, number_of_training_iterations, number_of_initial_hidden_units, number_of_sample_steps, number_of_runs, experiment_name, circuit_generator_script

noctua_partition = 'batch'
max_wall_time = '12:00:00'
noctua_user = 'hpc-prf-nqs'
singularity_image_location = "{pc2pfs}/{noctua_user}/nqs.sif".format(
                        noctua_user=noctua_user,
                        pc2pfs=os.environ["PC2PFS"])



jobDirs = []

def wait_for_job_queue():
    out = subprocess.Popen(['squeue'], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT)
    stdout,stderr = out.communicate()
    running_jobs = len(stdout.split('\n'))
    while(running_jobs > 99):
        time.sleep(10)
        out = subprocess.Popen(['squeue'], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT)
        stdout,stderr = out.communicate()
        running_jobs = len(stdout.split('\n'))

for qubits in number_of_qubits:
    for cycles in number_of_cycles:
        for circuit in range(number_of_circuits):

            circuitDirectory = "{pc2pfs}/{noctua_user}/{experiment_name}/{qubits}qubits/{cycles}cycles/circuit{circuit}".format(noctua_user=noctua_user,
                                                pc2pfs=os.environ["PC2PFS"],
                                                experiment_name=experiment_name,
                                                qubits=qubits,
                                                cycles=cycles,
                                                circuit=circuit)
            try: os.makedirs(circuitDirectory)
            except OSError, err:
                # Reraise the error unless it's about an already existing directory 
                if err.errno != errno.EEXIST or not os.path.isdir(circuitDirectory): 
                    raise

            wait_for_job_queue()
            bashCommand = "python {home}/nqs/scripts/{circuit_generator_script} {qubits} {cycles} {circuitDirectory}/in.qc".format(home=os.environ["HOME"], circuit_generator_script=circuit_generator_script, qubits=qubits, cycles=cycles, circuitDirectory=circuitDirectory)
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

            batch_script ="""#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -J {experiment_name}-exact
#SBATCH -A {noctua_user}
#SBATCH -p {noctua_partition}
#SBATCH -t {max_wall_time}

module reset
module load singularity
module load mpi/OpenMPI/3.1.4-GCC-8.3.0
export OMP_NUM_THREADS=1

mpirun -mca pml cm -mca mtl psm2 --report-bindings singularity exec {singularity_image_location} python2.7 $HOME/nqs/scripts/qasm_reader.py 0 0 0 0 exact > out 2> err""".format(
                        experiment_name=experiment_name,
                        noctua_user=noctua_user,
                        noctua_partition=noctua_partition,
                        max_wall_time=max_wall_time,
                        singularity_image_location=singularity_image_location
                    )

            f = open("{circuitDirectory}/job.slurm".format(circuitDirectory=circuitDirectory),'w')
            print >>f, batch_script

            for nodes in number_of_nodes:
                for tasks in number_of_tasks_per_node:
                    for threads in number_of_omp_threads:
                        for samples in number_of_training_samples:
                            for iterations in number_of_training_iterations:
                                for initial_hidden in number_of_initial_hidden_units:
                                    step_size = int(len(number_of_sample_steps)/len(number_of_qubits))
                                    index_sample_steps = number_of_qubits.index(qubits) * step_size
                                    for sample_steps in number_of_sample_steps[index_sample_steps:index_sample_steps + step_size]:
                                        for run in range(number_of_runs):


                                            directory = "{circuitDir}/{nodes}nodes/{tasks}tasks/{threads}threads/{samples}samples/{iterations}iterations/{initial_hidden}initialHidden/{sample_steps}sampleSteps/run{run}".format(
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
                                                run=run,
                                                circuitDir=circuitDirectory
                                            )

                                            try: os.makedirs(directory)
                                            except OSError, err:
                                                # Reraise the error unless it's about an already existing directory 
                                                if err.errno != errno.EEXIST or not os.path.isdir(directory): 
                                                    raise

                                            batch_script ="""#!/bin/bash
#SBATCH -N {nodes}
#SBATCH --ntasks-per-node={tasks}
#SBATCH -J {experiment_name}-{qubits}qubits-{cycles}cycles-circuit{circuit}-{nodes}nodes-{tasks}tasks-{threads}threads-{samples}samples-{iterations}iterations-run{run}
#SBATCH -A {noctua_user}
#SBATCH -p {noctua_partition}
#SBATCH -t {max_wall_time}

module reset
module load singularity
module load mpi/OpenMPI/3.1.4-GCC-8.3.0
export OMP_NUM_THREADS={threads}

cp {circuitDirectory}/in.qc {directory}/in.qc
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
                        singularity_image_location=singularity_image_location,
                        run=run,
                        pc2pfs=os.environ["PC2PFS"],
                        directory=directory,
                        circuitDirectory=circuitDirectory,
                        qubits=qubits,
                        cycles=cycles,
                        circuit=circuit
                    )

                                            f = open("{directory}/job.slurm".format(directory=directory),'w')
                                            print >>f, batch_script
                                            jobDirs.append(directory)

            jobDirs.append(circuitDirectory)

total_number_of_jobs = len(jobDirs)
job_number = 0
for directory in jobDirs:
    job_number += 1
    print('submitting job {} of {}'.format(job_number, total_number_of_jobs))
    wait_for_job_queue()
    bashCommand = "sbatch -D {directory} {directory}/job.slurm".format(directory=directory)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)  
    print "started job in {}".format(directory)           
