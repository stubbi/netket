# parameters to be tested
number_of_qubits = range(5,16,2)
number_of_cycles = range(10,21,5)
number_of_circuits = 5 #number of random circuits with same number of qubits and cycles

number_of_nodes = range(1,3)
number_of_tasks_per_node = range(1,3)
number_of_omp_threads = range(1,3)

number_of_training_samples = [100 + i * 200 for i in range(5)] 
number_of_training_iterations = [10000 + i * 20000 for i in range(5)]

number_of_initial_hidden_units = [0]
number_of_sample_steps = [q if q%2 != 0 else q+1 for q in number_of_qubits] #[9,11,13]#[3,5,7] #size must be multiple of qubits (n*size), each n entries will be used for corresponding qubits
number_of_runs = 5 #number of runs for a specific circuit
