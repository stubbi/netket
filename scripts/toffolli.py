#create Bell state
import nqs as nq
import cmath
import numpy as np


for x in range(2):
    for y in range(2):
        for z in range(2):
            nqs = nq.nqs.NQS(3)

            # init to 000
            nqs.applyHadamard(0,100,100000)
            nqs.applyHadamard(1,100,100000)
            nqs.applyHadamard(2,100,100000)

            if(x == 1):
                nqs.applyPauliX(0)
            if(y == 1):
                nqs.applyPauliX(1)
            if(z == 1):
                nqs.applyPauliX(2)

            nqs.truthTable()
            nqs.applyToffoli(0,1,2,100,100000)
            nqs.truthTable()