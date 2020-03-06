experiment_name = 'rcs'
# parameters to be tested
number_of_qubits = [4,6,8,10]
number_of_cycles = [4,6,8,10]
number_of_circuits = 10 #number of random circuits with same number of qubits and cycles

number_of_nodes = [1]
number_of_tasks_per_node = [1]
number_of_omp_threads = [1]

number_of_training_samples = [100 + i * 200 for i in range(5)] 
number_of_training_iterations = [10000 + i * 20000 for i in range(5)]

number_of_initial_hidden_units = [0]
number_of_sample_steps = [q if q%2 != 0 else q+1 for q in number_of_qubits] #[9,11,13]#[3,5,7] #size must be multiple of qubits (n*size), each n entries will be used for corresponding qubits
number_of_runs = 3 #number of runs for a specific circuit
