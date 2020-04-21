#program for Random Circuit creation.
# follows algorithm from https://arxiv.org/pdf/1608.00263.pdf
#expected file name is 'in.qc'
import random
import sys
import math

qasmFile = sys.argv[3]

# qubits will be initialized in |+> state
circuit = """
version 1.0
qubits 2

CZ q[0], q[1]
H q[1]
"""

f = open(qasmFile,'w')
print >>f, circuit
f.close()