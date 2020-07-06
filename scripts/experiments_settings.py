import math


def samples(n):
    return int(n * (math.log(n) + 0.577216) + 1.0/2.0)


experiment_name = 'final-4-SR-restarts-learned'
circuit_generator_script = 'random_circuit.py'
# parameters to be tested
number_of_qubits = [4]
n = 4
number_of_cycles = [5, 10, 15, 20]
number_of_circuits = 20 #number of random circuits with same number of qubits and cycles

number_of_nodes = [1]
number_of_tasks_per_node = [1]
number_of_omp_threads = [1]


number_of_training_samples = [samples(n), samples(n**2), samples(n**3), samples(2**n), samples(0.95* 2**n), samples(0.9 * 2**n), samples(0.85 * 2**n)] 
number_of_training_iterations = [1000, 10000, 100000]

number_of_initial_hidden_units = [int(q*(q-1)/2.0) for q in number_of_qubits]
number_of_sample_steps = [0]#q if q%2 != 0 else q+1 for q in number_of_qubits]
number_of_runs = 5 #number of runs for a specific circuit

randomRestarts = 5
earlyStopping = False
# AdaDelta, AdaGrad, AdaMax, AMSGrad, Momentum, RMSProp, Sgd, StochasticReconfiguration
optimizer = 'StochasticReconfiguration'
learnCZ = True


