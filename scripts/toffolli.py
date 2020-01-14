#create Bell state
import nqs as nq
import cmath
import numpy as np
import time

samples = int(sys.argv[1])
epochs = int(sys.argv[2])
initialHidden = int(sys.argv[3])
sampleSteps = int(sys.argv[4])

for x in range(2):
    for y in range(2):
        for z in range(2):
            nqs = nq.nqs.NQS(3, initialHidden, sampleSteps)

            # init to 000
            nqs.applyHadamard(0,samples,epochs)
            nqs.applyHadamard(1,samples,epochs)
            nqs.applyHadamard(2,samples,epochs)

            if(x == 1):
                nqs.applyPauliX(2)
            if(y == 1):
                nqs.applyPauliX(1)
            if(z == 1):
                nqs.applyPauliX(0)

            nqs.truthTable()
            start = time.time()
            nqs.applyToffoli(0,1,2,samples,epochs)
            end = time.time()
            print(end - start)
            nqs.truthTable()