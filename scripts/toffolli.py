#create Bell state
import nqs as nq
import cmath
import numpy as np


for x in range(2):
    for y in range(2):
        for z in range(2):
            for samples in [10,100,1000,10000]:
                for epochs in [100,1000,10000,100000,1000000,10000000]:
                    nqs = nq.nqs.NQS(3)

                    # init to 000
                    nqs.applyHadamard(0,samples,epochs)
                    nqs.applyHadamard(1,samples,epochs)
                    nqs.applyHadamard(2,samples,epochs)

                    if(x == 1):
                        nqs.applyPauliX(0)
                    if(y == 1):
                        nqs.applyPauliX(1)
                    if(z == 1):
                        nqs.applyPauliX(2)

                    nqs.truthTable()
                    nqs.applyToffoli(0,1,2,samples,epochs)
                    nqs.truthTable()