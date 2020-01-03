import subprocess

script = 'bell.py'
experiment_name = 'bell-test-scaling'
noctua_partition = 'batch'
max_wall_time = '00:10:00'
email = 'stubbi@mail.upb.de'

# parameters to be tested
number_of_nodes = range(1,40)
number_of_tasks_per_node = range(1,10)
number_of_omp_threads = range(1,40)

number_of_training_samples = [1000]#range(1000, 100000, 1000)
number_of_training_iterations = [1000]#range(1000, 100000, 1000)

for nodes in number_of_nodes:
    for tasks in number_of_tasks_per_node:
        for threads in number_of_omp_threads:
            for samples in number_of_training_samples:
                for iterations in number_of_training_iterations:

                    batch_script = f'''
                    #!/bin/bash
                    #SBATCH -N {number_of_nodes}
                    #SBATCH -J {experiment_name}-{nodes}nodes-{tasks}tasks-{threads}threads-{samples}samples-{iterations}iterations
                    #SBATCH -A hpc-prf-nqs
                    #SBATCH -p {noctua_partition}
                    #SBATCH -t {max_wall_time}
                    #SBATCH --mail-type all
                    #SBATCH --mail-user {email}

                    module reset
                    module load singularity
                    module load mpi/OpenMPI/3.1.4-GCC-8.3.0
                    export OMP_NUM_THREADS={number_of_omp_threads}

                    singularity pull --name nqs.sif shub://stubbi/nqs
                    mpirun singularity exec nqs.sif python2.7 $HOME/nqs/scripts/{script} > out 2> err'''

                    f = open('job.slurm','w')
                    print(batch_script, file=f)

                    bashCommand = "sbatch -D $PC2PFS/hpc-prf-nqs/{experiment_name}/{nodes}nodes/{tasks}tasks/{threads}threads/{samples}samples/{iterations}iterations job.slurm"
                    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)