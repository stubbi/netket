experiment_name = 'rcs-16-12'
circuit_generator_script = 'random_circuit.py'
# parameters to be tested
number_of_qubits = [16]
number_of_cycles = [12]
number_of_circuits = 3 #number of random circuits with same number of qubits and cycles

number_of_nodes = [1]
number_of_tasks_per_node = [1]
number_of_omp_threads = [1]

number_of_training_samples = [100 + i * 400 for i in range(3)] 
number_of_training_iterations = [10000 + i * 40000 for i in range(3)]

number_of_initial_hidden_units = [0]
number_of_sample_steps = [0,7,15]#[q if q%2 != 0 else q+1 for q in number_of_qubits]
number_of_runs = 3 #number of runs for a specific circuit
