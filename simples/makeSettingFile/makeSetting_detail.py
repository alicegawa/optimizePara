import numpy as np
import sys

argvs = sys.argv
argc = len(argvs)
if argc < 2:
    print 'python makeSetting.py <NCELL>'
    quit()

f = open("params.txt","w")

NCELL = int(argvs[1])
print 'NCELL = ', NCELL

NCON = NCELL * NCELL
f.write("#NAME MIN MAX TARGET(sim) LOG\n")

for i in range(0, NCON * 2):
    if i < NCON:
        f.write("weight%d\t0\t1\t0.1\t0\n" % (i+1))
    else:
        f.write("pos%d\t0\t10\t5\t0\n" % (i+1))

# for i in range(0, 72):
#     if i < 36:
#         f.write("weight%d\t0\t1\t0.1\t0\n" % (i+1))
#     else:
#         f.write("delay%d\t0.1\t10\t5\t0\n" % (i-35))
f.close()

f = open("rev_potential.dat", "w")
for i in  range(0, NCELL):
    for j in range(0, NCELL):
        if j < (NCELL - 1):
            f.write("%d\t" % ( (j%2) * -40))
        else:
            f.write("0")
    f.write("\n")
f.close()
