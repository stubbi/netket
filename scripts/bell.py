#program for Bell state creation.
#expected file name is 'in.qc'

program = """version 1.0
qubits 2

# Create the Bell state with CZ gate
H q[0:1]
CZ q[0], q[1]
H q[1]"""

f = open('in.qc','w')
print >>f, program
f.close()