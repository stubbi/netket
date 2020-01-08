#create Bell state
import nqs as nq
import cmath
import numpy as np
import time


nqs = nq.nqs.NQS(2)
nqs.applyControlledZRotation(0, 1, cmath.pi)
nqs.applyHadamard(1, 100, 100000)
start = time.time()
nqs.applyHadamard(1, 100, 100000)
end = time.time()
nqs.truthTable()
print(end - start)