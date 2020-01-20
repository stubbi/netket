#program for Toffolli gate.
#expected file name is 'in.qc'

program = """version 1.0
qubits 3
Toffoli q[0], q[1], q[2]
"""

f = open('in.qc','w')
print >>f, program
f.close()