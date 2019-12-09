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

            a = 0.0
            b = 0.0
            c = 0.0
            d = 0.0
            e = 0.0
            f = 0.0
            g = 0.0
            h = 0.0


            for i in range(10000):
                sample = nqs.sample()
                if(sample[0] == 0 and sample[1] == 0 and sample[2] == 0):
                    a = a + 1.0
                if(sample[0] == 0 and sample[1] == 0 and sample[2] == 1):
                    b = b + 1.0
                if(sample[0] == 0 and sample[1] == 1 and sample[2] == 0):
                    c = c + 1.0
                if(sample[0] == 0 and sample[1] == 1 and sample[2] == 1):
                    d = d + 1.0
                if(sample[0] == 1 and sample[1] == 0 and sample[2] == 0):
                    e = e + 1.0
                if(sample[0] == 1 and sample[1] == 0 and sample[2] == 1):
                    f = f + 1.0
                if(sample[0] == 1 and sample[1] == 1 and sample[2] == 0):
                    g = g + 1.0
                if(sample[0] == 1 and sample[1] == 1 and sample[2] == 1):
                    h = h + 1.0
                    
            print("init:")        
            print(a/10000.0)
            print(b/10000.0)
            print(c/10000.0)
            print(d/10000.0)
            print(e/10000.0)
            print(f/10000.0)
            print(g/10000.0)
            print(h/10000.0)
            print()



            nqs.applyToffoli(0,1,2,100,100000)

            a = 0.0
            b = 0.0
            c = 0.0
            d = 0.0
            e = 0.0
            f = 0.0
            g = 0.0
            h = 0.0


            for i in range(10000):
                sample = nqs.sample()
                if(sample[0] == 0 and sample[1] == 0 and sample[2] == 0):
                    a = a + 1.0
                if(sample[0] == 0 and sample[1] == 0 and sample[2] == 1):
                    b = b + 1.0
                if(sample[0] == 0 and sample[1] == 1 and sample[2] == 0):
                    c = c + 1.0
                if(sample[0] == 0 and sample[1] == 1 and sample[2] == 1):
                    d = d + 1.0
                if(sample[0] == 1 and sample[1] == 0 and sample[2] == 0):
                    e = e + 1.0
                if(sample[0] == 1 and sample[1] == 0 and sample[2] == 1):
                    f = f + 1.0
                if(sample[0] == 1 and sample[1] == 1 and sample[2] == 0):
                    g = g + 1.0
                if(sample[0] == 1 and sample[1] == 1 and sample[2] == 1):
                    h = h + 1.0
                    
            print("out")
            print(a/10000.0)
            print(b/10000.0)
            print(c/10000.0)
            print(d/10000.0)
            print(e/10000.0)
            print(f/10000.0)
            print(g/10000.0)
            print(h/10000.0)
            print()