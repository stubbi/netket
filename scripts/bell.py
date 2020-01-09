#create Bell state
import nqs as nq
import cmath
import numpy as np
import time
import sys


nqs = nq.nqs.NQS(2)
nqs.applyControlledZRotation(0, 1, cmath.pi)
start = time.time()
nqs.applyHadamard(1, int(sys.argv[1]), int(sys.argv[2]))
end = time.time()
nqs.truthTable()
print(end - start)