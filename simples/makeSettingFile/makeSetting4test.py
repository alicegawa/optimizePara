import numpy as np
import sys
import math as m

argvs = sys.argv
argc = len(argvs)
if argc < 2:
    print 'python makeSetting.py <ndim>'
    quit()

f = open("params_test.txt","w")

ndim = int(argvs[1])
print 'ndim = ', ndim

f.write("#NAME MIN MAX TARGET(sim) LOG\n")

for i in range(0, ndim):
    f.write("param%d\t-1\t1\t0\t0\n" % (i+1))
f.close()

judge_sq = m.sqrt(ndim) - int(m.sqrt(ndim))
if(judge_sq > 0):
    print "judge_sq is not a squared num\n"
else:
    f = open("conMat_test.txt", "w")
    print(m.sqrt(ndim))
    for i in range(0, int(m.sqrt(ndim))):
        for j in range(0, int(m.sqrt(ndim))):
            f.write("1\t")
        f.write("\n")
f.close()
