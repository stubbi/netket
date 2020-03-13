#program for Random Circuit creation.
# follows algorithm from https://arxiv.org/pdf/1608.00263.pdf
#expected file name is 'in.qc'
import random
import sys
import math

numQubits = int(sys.argv[1])
numCycles = int(sys.argv[2])
qasmFile = sys.argv[3]

singleGates = ['sqrt_X']

root = math.floor(math.sqrt(numQubits))
single_qubit_cycles = [['n/a' for _ in range(numQubits)]]
def generateCycle(m):
    cycle = 'cycle {}\n'.format(m)

    for q in range(numQubits):
        gate = random.choice(singleGates)
        cycle = cycle + "{} q[{}]\n".format(gate, q)

    cycle = cycle + '\n\n'

    return cycle
    


# qubits will be initialized in |+> state
circuit = """
version 1.0
qubits {}
H q[0]
""".format(numQubits)

for m in range(numCycles):
    circuit = circuit + generateCycle(m)

f = open(qasmFile,'w')
print >>f, circuit
f.close()