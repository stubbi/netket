experiment_name = 'rcsIterations'
circuit_generator_script = 'random_circuit.py'
# parameters to be tested
number_of_qubits = [6]
number_of_cycles = [8,10]
number_of_circuits = 5 #number of random circuits with same number of qubits and cycles

number_of_nodes = [1]
number_of_tasks_per_node = [1]
number_of_omp_threads = [1]

number_of_training_samples = [100 + i * 400 for i in range(3)] 
number_of_training_iterations = [1000 + i * 4000 for i in range(3)]

number_of_initial_hidden_units = [0]
number_of_sample_steps = [q if q%2 != 0 else q+1 for q in number_of_qubits]
number_of_runs = 1 #number of runs for a specific circuit
