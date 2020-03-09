#program for Random Circuit creation.
# follows algorithm from https://arxiv.org/pdf/1608.00263.pdf
#expected file name is 'in.qc'
import random
import sys
import math

numQubits = int(sys.argv[1])
numCycles = int(sys.argv[2])
qasmFile = sys.argv[3]

singleGates = ['sqrt_Y','sqrt_X','T']

root = math.floor(math.sqrt(numQubits))
single_qubit_cycles = [['n/a' for _ in range(numQubits)]]
def generateCycle(m):

    #for conversion to 2d lattice of qubits
    def row(qubit):
        return qubit % root

    def column(qubit):
        return math.floor(qubit/root)
    
    def pairIndex(index):
        if(index % 2 != 0):
            index -=1
        return index/2.0 

    def even(i):
        return i % 2 == 0

    czs = []
    #figure 6 of the paper involve all qubits for faster conversion to porter-thomas distribution
    if(m % 4 == 0):
        czs.extend([(i,j) for i in range(numQubits) for j in range(numQubits) if row(i) == row(j) and column(i) + 1 == column(j) and ((even(pairIndex(column(i))) and even(pairIndex(column(j))) and not even(row(i))) or (not even(pairIndex(column(i))) and not even(pairIndex(column(j))) and even(row(i))))])
        czs.extend([(i,j) for i in range(numQubits) for j in range(numQubits) if row(i) == row(j) and column(i) + 1 == column(j) and ((not even(pairIndex(column(i))) and not even(pairIndex(column(j))) and not even(row(i))) or (even(pairIndex(column(i))) and even(pairIndex(column(j))) and even(row(i))))])
    elif(m % 4 == 1):
        czs.extend([(i,j) for i in range(numQubits) for j in range(numQubits) if column(i) == column(j) and row(i) + 1 == row(j) and ((not even(pairIndex(row(i))) and even(pairIndex(row(j))) and even(column(i))) or (even(pairIndex(row(i))) and not even(pairIndex(row(j))) and not even(column(i))))])
        czs.extend([(i,j) for i in range(numQubits) for j in range(numQubits) if column(i) == column(j) and row(i) + 1 == row(j) and ((not even(pairIndex(row(i))) and even(pairIndex(row(j))) and not even(column(i))) or (even(pairIndex(row(i))) and not even(pairIndex(row(j))) and even(column(i))))])
    elif(m % 4 == 2):
        czs.extend([(i,j) for i in range(numQubits) for j in range(numQubits) if row(i) == row(j) and column(i) + 1 == column(j) and ((not even(pairIndex(column(i))) and even(pairIndex(column(j))) and not even(row(i))) or (even(pairIndex(column(i))) and not even(pairIndex(column(j))) and even(row(i))))])
        czs.extend([(i,j) for i in range(numQubits) for j in range(numQubits) if row(i) == row(j) and column(i) + 1 == column(j) and ((not even(pairIndex(column(i))) and even(pairIndex(column(j))) and even(row(i))) or (even(pairIndex(column(i))) and not even(pairIndex(column(j))) and not even(row(i))))])
    elif(m % 4 == 3):
        czs.extend([(i,j) for i in range(numQubits) for j in range(numQubits) if column(i) == column(j) and row(i) + 1 == row(j) and ((even(pairIndex(row(i))) and even(pairIndex(row(j))) and even(column(i))) or (not even(pairIndex(row(i))) and not even(pairIndex(row(j))) and not even(column(i))))])
        czs.extend([(i,j) for i in range(numQubits) for j in range(numQubits) if column(i) == column(j) and row(i) + 1 == row(j) and ((even(pairIndex(row(i))) and even(pairIndex(row(j))) and not even(column(i))) or (not even(pairIndex(row(i))) and not even(pairIndex(row(j))) and even(column(i))))])


    cycle = 'cycle {}\n'.format(m)

    for c in czs:
        cycle = cycle + "CZ q[{}], q[{}]\n".format(c[0], c[1])
    
    cycle = cycle + '\n'

    single_qubit_gates = []
    for q in range(numQubits):
        gate = random.choice([g for g in singleGates if g != single_qubit_cycles[m][q]])
        single_qubit_gates.append(gate)
        cycle = cycle + "{} q[{}]\n".format(gate, q)

    single_qubit_cycles.append(single_qubit_gates)
    cycle = cycle + '\n\n'

    return cycle
    


# qubits will be initialized in |+> state
circuit = """
version 1.0
qubits {}
""".format(numQubits)

for m in range(numCycles):
    circuit = circuit + generateCycle(m)

f = open(qasmFile,'w')
print >>f, circuit
f.close()