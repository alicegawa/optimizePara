import numpy as np

f = open("params.txt","w")

f.write("#NAME MIN MAX TARGET(sim) LOG\n")
for i in range(0, 50):
    if i < 25:
        f.write("weight%d\t0\t1\t0.1\t0\n" % (i+1))
    else:
        f.write("delay%d\t0.1\t10\t1\t0\n" % (i-24))
f.close()
