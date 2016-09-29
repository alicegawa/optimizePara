import numpy as np
import os
import sys

argvs = sys.argv
argc = len(argvs)
#num_lines = sum(1 for line in open('./swc/0008.swc'))
if argc > 2:
    files = os.listdir(argvs[1])
    filename_tmp = argvs[1]
else:
    files = os.listdir("./swc/")
    filename_tmp = "./swc/"
#print files
#print files[0]
#filename_tmp = "./swc/"
filename_out = "./sectionData.txt"
filename_out_params = "./params.txt"
f = open(filename_out, 'w')
f_params = open(filename_out_params, 'w')
f_params.write("#NAME MIN MAX TARGET(sim) LOG\n")
loop = 0
for file in files:
    filename = filename_tmp + file
    print filename
    data = np.loadtxt(filename, comments='#', delimiter=' ')
    print data.shape[0]
    sectionNum = data.shape[0] - 1
    f.write("%s: %d\n"%(filename, sectionNum))
    f_params.write("pos%d\t0\t%d\t%lf\t0\n"%(loop, sectionNum, sectionNum / 2))
    loop = loop + 1
f.close()
f_params.close()
