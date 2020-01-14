#create Bell state
import nqs as nq
import cmath
import numpy as np
import time
import sys

samples = int(sys.argv[1])
epochs = int(sys.argv[2])
initialHidden = int(sys.argv[3])
sampleSteps = int(sys.argv[4])

nqs = nq.nqs.NQS(2, initialHidden, sampleSteps)
nqs.applyControlledZRotation(0, 1, cmath.pi)
start = time.time()
nqs.applyHadamard(1, samples, epochs)
end = time.time()
nqs.truthTable()
print(end - start)