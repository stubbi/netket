#create Bell state
import nqs
import cmath
import numpy as np

nqs = nqs.nqs.NQS(2)
nqs.applyControlledZRotation(0, 1, cmath.pi)
nqs.applyHadamard(1, 100, 100000)

a = 0.0
b = 0.0
c = 0.0
d = 0.0

for i in range(10000):
    sample = nqs.sample()
    if(sample[0] == 0 and sample[1] == 0):
        a = a + 1.0
    if(sample[0] == 0 and sample[1] == 1):
        b = b + 1.0
    if(sample[0] == 1 and sample[1] == 0):
        c = c + 1.0
    if(sample[0] == 1 and sample[1] == 1):
        d = d + 1.0
        
print(a/10000.0)
print(b/10000.0)
print(c/10000.0)
print(d/10000.0)