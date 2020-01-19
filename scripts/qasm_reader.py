import cmath
import numpy as np
import nqs as nq
import collections
import json
import sys

samples = int(sys.argv[1])
epochs = int(sys.argv[2])
initialHidden = int(sys.argv[3])
sampleSteps = int(sys.argv[4])
shots = 1000#int(sys.argv[5])

class QASMReader:
    def __init__(self, filename, numSamples, numIterations, numInitialHidden, numSampleSteps):
        self.filename = filename
        self.numSamples = numSamples
        self.numIterations = numIterations
        self.numInitialHidden = numInitialHidden
        self.numSampleSteps = numSampleSteps

    def buildCircuit(self):
        f = open(self.filename)
        for line in f:
            self.circuitFromLine(line)
        f.close()

    def circuitFromLine(self, line):
        def readQubits(self, line):
            line = line.strip(' q[]')
            qubits = []
            for s in line.split(','):
                if(len(s.split(':')) > 1):
                    for i in range(int(s.split(':')[0]), int(s.split(':')[1]) + 1):
                        qubits.append(i)
                else:
                    qubits.append(int(s))
            return qubits

        line = line.split('#')[0].strip()
        if(line.startswith('qubits')):
            qubits = int(line[7:])
            self.nqs = nq.nqs.NQS(qubits, self.numInitialHidden, self.numSampleSteps)
            for q in qubits:
                self.nqs.applyHadamard(q, self.numSamples, self.numIterations)
        #elif(line.startswith('prep_z')):
        #    for q in self.readQubits(line.strip('prep_z')):
        #        #self.nqs = nq.nqs.NQS(int(line[7:]), self.numInitialHidden, self.numSampleSteps)
        elif(line.startswith('X')):
            for q in self.readQubits(line.strip('X')):
                self.nqs.applyPauliX(q)
        elif(line.startswith('Y')):
            for q in self.readQubits(line.strip('Y')):
                self.nqs.applyPauliZ(q)
        elif(line.startswith('Z')):
            for q in self.readQubits(line.strip('Z')):
                self.nqs.applyPauliZ(q)
        elif(line.startswith('H')):
            for q in self.readQubits(line.strip('H')):
                self.nqs.applyHadamard(q, self.numSamples, self.numIterations)
        elif(line.startswith('Rz')):
            for q in self.readQubits(line.split(',')[0].strip('Rz')):
                self.nqs.applySingleZRotation(q, float(line.split(',')[1]))
        elif(line.startswith('T')):
            for q in self.readQubits(line.strip('T')):
                self.nqs.applyT(q)
        elif(line.startswith('Tdag')):
            for q in self.readQubits(line.strip('Tdag')):
                self.nqs.applyTDagger(q)
        elif(line.startswith('CZ')):
            q0 = self.readQubits(line.strip('CZ').split(',')[0])[0]
            q1 = self.readQubits(line.strip('CZ').split(',')[1])[0]
            self.nqs.applyControlledZRotation(q0, q1, cmath.pi)
        else:
            # no CNOT, named and repeated subcircuits among others...
            pass


    def display(self, suffix):
        #TODO use suffix
        def toDecimal(self, sample):
            return int(bin(int(sample)), 2)
        raw_data = (toDecimal(self.nqs.sample()) for _ in range(shots))
        histogram = collections.Counter(raw_data)

        with open('raw_data.json', 'w') as f:
            json.dump(raw_data, f)

        with open('histogram.json', 'w') as f:
            json.dump(histogram, f)



program = """version 1.0
# number of qubits for this backend is limited to 5 Qubits
qubits 2

# Create the Bell state with CZ gate
H q[0:1]
CZ q[0], q[1]
H q[1]"""

f = open('in.qc','w')
print >>f, program

qasm = QASMReader("in.qc", 100, 10000, 0, 0)
qasm.buildCircuit()
qasm.display("")
        

