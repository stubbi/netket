import subprocess, os, errno

script = 'bell.py'
experiment_name = 'bell-test-scaling'
noctua_partition = 'short'
max_wall_time = '00:30:00'
email = 'stubbi@mail.upb.de'
noctua_user = 'hpc-prf-nqs'
singularity_image_location = "{pc2pfs}/{noctua_user}/nqs.sif".format(
                        noctua_user=noctua_user,
                        pc2pfs=os.environ["PC2PFS"])

# parameters to be tested
number_of_nodes = range(1,3)
number_of_tasks_per_node = range(1,10)
number_of_omp_threads = range(1,10)

number_of_training_samples = [1000]#range(1000, 100000, 1000)
number_of_training_iterations = [1000]#range(1000, 100000, 1000)

for nodes in number_of_nodes:
    for tasks in number_of_tasks_per_node:
        for threads in number_of_omp_threads:
            for samples in number_of_training_samples:
                for iterations in number_of_training_iterations:

                    batch_script ="""#!/bin/bash
#SBATCH -N {nodes}
#SBATCH --ntasks-per-node={tasks}
#SBATCH -J {experiment_name}-{nodes}nodes-{tasks}tasks-{threads}threads-{samples}samples-{iterations}iterations
#SBATCH -A {noctua_user}
#SBATCH -p {noctua_partition}
#SBATCH -t {max_wall_time}
#SBATCH --mail-type fail
#SBATCH --mail-user {email}

module reset
module load singularity
module load mpi/OpenMPI/3.1.4-GCC-8.3.0
export OMP_NUM_THREADS={threads}

mpirun -mca pml cm -mca mtl psm2 --report-bindings singularity exec {singularity_image_location} python2.7 $HOME/nqs/scripts/{script} {samples} {iterations}> out 2> err""".format(
                        nodes=nodes,
                        experiment_name=experiment_name,
                        tasks=tasks,
                        threads=threads,
                        samples=samples,
                        iterations=iterations,
                        noctua_user=noctua_user,
                        noctua_partition=noctua_partition,
                        max_wall_time=max_wall_time,
                        email=email,
                        script=script,
                        singularity_image_location=singularity_image_location
                    )

                    f = open('job.slurm','w')
                    print >>f, batch_script

                    directory = "{pc2pfs}/{noctua_user}/{experiment_name}/{nodes}nodes/{tasks}tasks/{threads}threads/{samples}samples/{iterations}iterations".format(
                        noctua_user=noctua_user,
                        pc2pfs=os.environ["PC2PFS"],
                        experiment_name=experiment_name,
                        nodes=nodes,
                        tasks=tasks,
                        threads=threads,
                        samples=samples,
                        iterations=iterations
                    )

                    try: os.makedirs(directory)
                    except OSError, err:
                        # Reraise the error unless it's about an already existing directory 
                        if err.errno != errno.EEXIST or not os.path.isdir(directory): 
                            raise

                    bashCommand = "sbatch -D {directory} job.slurm".format(directory=directory)
                    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

                    print "started job {experiment_name} for {nodes}nodes {tasks}tasks {threads}threads {samples}samples {iterations}iterations job.slurm".format(
                        experiment_name=experiment_name,
                        nodes=nodes,
                        tasks=tasks,
                        threads=threads,
                        samples=samples,
                        iterations=iterations
                    )