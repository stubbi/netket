version 1.0
qubits 2

# Create the Bell state with CZ gate
H q[0:1]
CZ q[0], q[1]
H q[1]