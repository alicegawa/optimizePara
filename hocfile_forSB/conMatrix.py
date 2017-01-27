#!/usr/bin/python

import csv
import sys

#--- command line ---
argvs = sys.argv
argc = len(argvs)
if (argc != 3):
    print 'Usage: # python %s datafile matrixfile' % argvs[0]
    quit()
datafile = argvs[1]
matrixfile = argvs[2]
       
print "data:",datafile, " , matrix:",matrixfile
 
outputReader = csv.reader(open(datafile, 'rb'))
outputReader.next()
outputReader.next()
outputReader.next()

matrixWriter = csv.writer(open(matrixfile, 'wb'))

for row_o in outputReader:
    if row_o[0] != '':
        print "output:", row_o[0]
        matrixWriter.writerow(['output', row_o[0]])

        inputReader = csv.reader(open(datafile, 'rb'))
        inputReader.next()
        inputReader.next()
        inputReader.next()

        for row_i in inputReader:
            if row_i[0] != '':
                R_oLAL = float(row_o[9])  * float(row_i[19])
                R_iLAL = float(row_o[10]) * float(row_i[20])
                R_oVPC = float(row_o[11]) * float(row_i[21])
                R_iVPC = float(row_o[12]) * float(row_i[22])
                R_aiVPC = float(row_o[13]) * float(row_i[23])
                L_oLAL = float(row_o[14]) * float(row_i[24])
                L_iLAL = float(row_o[15]) * float(row_i[25])
                L_oVPC = float(row_o[16]) * float(row_i[26])
                L_iVPC = float(row_o[17]) * float(row_i[27])
                L_aiVPC = float(row_o[18]) * float(row_i[28])
                print "input :", row_i[0],\
                    R_oLAL, R_iLAL, R_oVPC, R_iVPC, R_aiVPC,\
                    L_oLAL, L_iLAL, L_oVPC, L_iVPC, L_aiVPC
                matrixWriter.writerow(['input', row_i[0],
                                       R_oLAL, R_iLAL, R_oVPC, R_iVPC, R_aiVPC,
                                       L_oLAL, L_iLAL, L_oVPC, L_iVPC, L_aiVPC])
                


            
