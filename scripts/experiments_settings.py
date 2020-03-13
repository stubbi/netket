experiment_name = 'debug_learning'
circuit_generator_script = 'debug_circuit.py'
# parameters to be tested
number_of_qubits = [1]
number_of_cycles = [1,2]
number_of_circuits = 1 #number of random circuits with same number of qubits and cycles

number_of_nodes = [1]
number_of_tasks_per_node = [1]
number_of_omp_threads = [1]

number_of_training_samples = [100 + i * 200 for i in range(3)] 
number_of_training_iterations = [10000 + i * 20000 for i in range(3)]

number_of_initial_hidden_units = [1]
number_of_sample_steps = [q if q%2 != 0 else q+1 for q in number_of_qubits]
number_of_runs = 3 #number of runs for a specific circuit
