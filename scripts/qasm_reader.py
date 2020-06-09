import cmath
import nqs as nq
import collections
import json
import sys
import pickle
from qupy import Qubits
from qupy.operator import H, X, Y, Z, T, Tdag, rz, sqrt_X, sqrt_Y
import time

samples = int(sys.argv[1])
epochs = int(sys.argv[2])
initialHidden = int(sys.argv[3])
sampleSteps = int(sys.argv[4])
randomRestarts = int(sys.argv[5])
earlyStopping = int(sys.argv[6])
optimizer = int(sys.argv[7])
method = str(sys.argv[8])
shots = 1000


class QASMReader:
    def __init__(self, method, numSamples, numIterations, numInitialHidden, numSampleSteps, numRandomRestarts, earlyStopping, optimizer):
        assert(method == 'nqs' or method == 'exact')
        self.numSamples = numSamples
        self.numIterations = numIterations
        self.numInitialHidden = numInitialHidden
        self.numSampleSteps = numSampleSteps
        self.numRandomRestarts = numRandomRestarts
        self.earlyStopping = earlyStopping
        self.optimizer = optimizer
        self.nqs = None
        self.exact = None
        self.method = method
        self.start = None
        self.end = None
        self.gateNo = 0

    def is_nqs(self):
        return self.method == 'nqs'

    def buildCircuit(self, filename):
        f = open(filename)
        oldGateNo = self.gateNo #to avoid storing data without changes
        for line in f:
            self.circuitFromLine(line)
            if(oldGateNo != self.gateNo):
                if(self.is_nqs()):
                    self.nqs.save('parameters_gate_{}.json'.format(self.gateNo))
                else:
                    with open('exact_gate_{}.json'.format(self.gateNo), 'wb') as f:
                        pickle.dump(self.exact.get_state(), f)
            oldGateNo = self.gateNo
        f.close()

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

    def circuitFromLine(self, line):
        self.gateNo += 1
        line = line.split('#')[0].strip()
        if(line.startswith('qubits')):
            qubits = int(line[7:])
            if(self.is_nqs()):
                self.nqs = nq.nqs.NQS(qubits, self.numInitialHidden, self.numSampleSteps, self.numRandomRestarts, self.earlyStopping, self.optimizer)
            else:
                self.exact = Qubits(qubits)
                for q in range(qubits):
                    self.exact.gate(H, target=q)

        elif(line.startswith('X')):
            for q in self.readQubits(line.strip('X')):
                if(self.is_nqs()):
                    self.nqs.applyPauliX(q)
                else:
                    self.exact.gate(X, target=q)

        elif(line.startswith('Y')):
            for q in self.readQubits(line.strip('Y')):
                if(self.is_nqs()):
                    self.nqs.applyPauliY(q)
                else:
                    self.exact.gate(Y, target=q)

        elif(line.startswith('Z')):
            for q in self.readQubits(line.strip('Z')):
                if(self.is_nqs()):
                    self.nqs.applyPauliZ(q)
                else:
                    self.exact.gate(Z, target=q)

        elif(line.startswith('H')):
            for q in self.readQubits(line.strip('H')):
                if(self.is_nqs()):
                    self.nqs.applyHadamard(q, self.numSamples, self.numIterations)
                else:
                    self.exact.gate(H, target=q)

        elif(line.startswith('Rz')):
            for q in self.readQubits(line.split(',')[0].strip('Rz')):
                if(self.is_nqs()):
                    self.nqs.applySingleZRotation(q, float(line.split(',')[1]))
                else:
                    self.exact.gate(rz(float(line.split(',')[1])), target=q)

        elif(line.startswith('Toffoli')):
            q0 = self.readQubits(line.strip('Toffoli').split(',')[0])[0]
            q1 = self.readQubits(line.strip('Toffoli').split(',')[1])[0]
            q2 = self.readQubits(line.strip('Toffoli').split(',')[2])[0]
            if(self.is_nqs()):
                self.nqs.applyToffoli(q0, q1, q2, self.numSamples, self.numIterations)
            else:
                #https://arxiv.org/pdf/0803.2316.pdf
                self.exact.gate(H, target = q2)
                self.exact.gate(X, target = q2, control=q1)
                self.exact.gate(Tdag, target = q2)
                self.exact.gate(X, target = q2, control=q0)
                self.exact.gate(T, target = q2)
                self.exact.gate(X, target = q2, control=q1)
                self.exact.gate(Tdag, target = q2)
                self.exact.gate(X, target = q2, control=q0)
                self.exact.gate(T, target = q1)
                self.exact.gate(T, target = q2)
                self.exact.gate(H, target = q2)
                self.exact.gate(X, target = q1, control=q0)
                self.exact.gate(T, target = q0)
                self.exact.gate(Tdag, target = q1)
                self.exact.gate(X, target = q1, control=q0)

        elif(line.startswith('Tdag')):
            for q in self.readQubits(line.strip('Tdag')):
                if(self.is_nqs()):
                    self.nqs.applyTDagger(q)
                else:
                    self.exact.gate(T, target=q)

        elif(line.startswith('T')):
            for q in self.readQubits(line.strip('T')):
                if(self.is_nqs()):
                    self.nqs.applyT(q)
                else:
                    self.exact.gate(T, target=q)

        elif(line.startswith('sqrt_X')):
            for q in self.readQubits(line.strip('sqrt_X')):
                if(self.is_nqs()):
                    self.nqs.applySqrtX(q, self.numSamples, self.numIterations)
                else:
                    self.exact.gate(sqrt_X, target=q)

        elif(line.startswith('sqrt_Y')):
            for q in self.readQubits(line.strip('sqrt_Y')):
                if(self.is_nqs()):
                    self.nqs.applySqrtY(q, self.numSamples, self.numIterations)
                else:
                    self.exact.gate(sqrt_Y, target=q)
                    
        elif(line.startswith('CZ')):
            q0 = self.readQubits(line.strip('CZ').split(',')[0])[0]
            q1 = self.readQubits(line.strip('CZ').split(',')[1])[0]
            if(self.is_nqs()):
                self.nqs.applyControlledZRotation(q0, q1, cmath.pi, self.numSamples, self.numIterations)
            else:
                self.exact.gate(Z, target=q0, control=q1)

        else:
            # no CNOT, named and repeated subcircuits among others...
            self.gateNo -= 1

    def toDecimal(self, sample):
        return int("".join(str(int(x)) for x in sample), 2)

    def display(self):
        if(self.is_nqs()):
            startSampling = time.time()
            raw_data = [self.toDecimal(self.nqs.sample()) for _ in range(shots)]
            histogram = collections.Counter(raw_data)
            endSampling = time.time()

            with open('raw_data.json', 'w') as f:
                json.dump(raw_data, f)

            with open('histogram.json', 'w') as f:
                json.dump(histogram, f)

            with open('duration_sampling.time', 'w') as f:
                f.write(str(endSampling-startSampling))

            self.nqs.save('parameters.json')
        else:
            with open('exact.json', 'wb') as f:
                pickle.dump(self.exact.get_state(), f)

        with open('duration.time', 'w') as f:
            f.write(str(self.end-self.start))

qasm = QASMReader(method, samples, epochs, initialHidden, sampleSteps, randomRestarts, earlyStopping, optimizer)
qasm.start = time.time()
qasm.buildCircuit("in.qc")
qasm.end = time.time()
qasm.display()
