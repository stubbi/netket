import nqs as nq
import cmath
import numpy as np
import time
import sys

samples = int(sys.argv[1])
epochs = int(sys.argv[2])
initialHidden = int(sys.argv[3])
sampleSteps = int(sys.argv[4])

for x in range(2):
    for y in range(2):

        nqs = nq.nqs.NQS(2, initialHidden, sampleSteps)

        nqs.applyHadamard(0,samples,epochs)
        nqs.applyHadamard(1,samples,epochs)


        if(x == 1):
            nqs.applyPauliX(0)
        if(y == 1):
            nqs.applyPauliX(1)

        nqs.truthTable()
        nqs.applyHadamard(1,samples,epochs)
        nqs.applyControlledZRotation(0, 1, cmath.pi)
        nqs.applyHadamard(1,samples,epochs)
        nqs.truthTable()