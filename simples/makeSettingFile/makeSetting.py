import numpy as np

f = open("params.txt","w")

f.write("#NAME MIN MAX TARGET(sim) LOG\n")
for i in range(0, 72):
    if i < 36:
        f.write("weight%d\t0\t1\t0.1\t0\n" % (i+1))
    else:
        f.write("delay%d\t0.1\t10\t5\t0\n" % (i-35))
f.close()
