#program for Random Circuit creation.
#expected file name is 'in.qc'
import random
import sys

numQubits = sys.argv[1]
numCycles = sys.argv[2]
qasmFile = sys.argv[3]

singleGates = ['X','Y','Z','H','T','Tdag']


def generateCycle(m):
    cycle = ''
    for q in range(numQubits):
        gate = random.choice(singleGates)
        cycle = cycle + "{} q[{}]\n".format(gate, q)
    
    #simulate Google's approach to connect each qubits to its for
    #neighbours one by one. Or at least kind of
    controls = [q for q in range(numQubits) if q%4 == m%4]
    for c in controls:
        cycle = cycle + "CZ q[{}], q[{}]\n".format(c, (c+1)%numQubits)
    



circuit = """
version 1.0
qubits {}
""".format(numQubits)

for m in range(numCycles):
    circuit = circuit + generateCycle(m)

f = open(qasmFile,'w')
print >>f, program
f.close()