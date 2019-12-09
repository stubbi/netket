#create Bell state
import nqs as nq
import cmath
import numpy as np

nqs = nq.nqs.NQS(2)
nqs.applyControlledZRotation(0, 1, cmath.pi)
nqs.applyHadamard(1, 100, 100000)

nqs.truthTable()