experiment_name = 'StochasticReconfiguration-5-10-5-False'
circuit_generator_script = 'random_circuit.py'
# parameters to be tested
number_of_qubits = [2]
number_of_cycles = [2]
number_of_circuits = 1 #number of random circuits with same number of qubits and cycles

number_of_nodes = [4]
number_of_tasks_per_node = [4]
number_of_omp_threads = [4]

number_of_training_samples = [500] 
number_of_training_iterations = [500]

number_of_initial_hidden_units = [int(q*(q-1)/2.0) for q in number_of_qubits]
number_of_sample_steps = [0]#q if q%2 != 0 else q+1 for q in number_of_qubits]
number_of_runs = 3 #number of runs for a specific circuit

randomRestarts = 5
earlyStopping = False
# AdaDelta, AdaGrad, AdaMax, AMSGrad, Momentum, RMSProp, Sgd, StochasticReconfiguration
optimizer = 'StochasticReconfiguration'
learnCZ = True
