#***********************
# getSynList.py
# 2012/01/22
# Yusuke Mori
#***********************

#***********************
# get synapse list
#***********************

#--- import module ---
from mpi4py import *
from neuron import *
from random import *
import csv
import math
import re
import os
import shutil


#--- object ---
connectlist = []
candlist = []
randsynlist = []

#--- import ---
connectfile = 'list/connectlist.csv'
datafile = 'data.csv'


#-----------------------
# Make Random Synapse List
# 
#-----------------------

#----------------------------------
# data
#
# row = (#0, #1, #2, #3, #4)
# #0. pre  synapse
# #1. post synapse
# #2. distance (square)
# #3. volume of pre compartment
# #4. volume of post compartment
def readDataCsv():
    fh = open(datafile, 'rb')
    reader = csv.reader(fh)
    for row in reader:
        row.append(1) #--- row[5] means synapse strength
        candlist.append(row)
    fh.close()

#----------------------------------
# get Probability
#
#  prob_func(dist) : probability distribution function 
#  probability(dist_threshold) :  
#
def prob_func(distance, volume1, volume2):
    d = math.sqrt(float(distance)) # d = real distance
    v1 = float(volume1)
    v2 = float(volume2)
    sigma = 10.0
    #--- Gaussian : f(x) = ( 1 / sqrt(2*pi*sigma^2) ) exp (-(x-u)^2 / 2 * sigma^2)
    Gaussian_0 = math.exp(-math.pow(0.0-d,2)/(2*math.pow(sigma, 2))) / math.sqrt(2*math.pi*math.pow(sigma,2)) 
    return Gaussian_0 * (v1*v2)

def probability(dist_threshold):
    global candlist
    sum = 0.0

    for row in candlist:
        #--- if distance is too large, exclude the row ---
        if (float(row[2]) > dist_threshold * dist_threshold):
            candlist.remove(row)
        else:
            #--- sum for normalization ---
            sum += prob_func(row[2], row[3], row[4])

    #--- calc probability ---
    for row in candlist:
        tmp = prob_func(row[2], row[3], row[4]) / sum
        row.append(tmp) #--- row[6] means probability

#----------------------------------
# export Random Synapse List
#
def randSynList(N, ratio, synType, filename = "synList.csv"):
    #--- global access (for reset) ---
    global candlist
    global randsynlist

    synfile = filename

    #--- file open ---
    fw = open(synfile, 'wb')
    writer = csv.writer(fw)

    #--- calc ---
    num = 0
    duplicate = -1

    print "--- candidate list ---"
    for row in candlist:
        print row

    while num < N:
        tmp = random()
        for row in candlist:
            # print row
            tmp = tmp - row[6]
            if tmp < 0 :
                #--- selected row ---
                num += 1
                #--- check duplicate ---
                if (randsynlist == []):
                    randsynlist.append(row[:])
                else:
                    duplicate = -1                    
                    for exist_row in randsynlist:
                        if row[0] == exist_row[0] or row[1] == exist_row[1]:
                            exist_row[5] = int(exist_row[5]) + int(row[5])
                            duplicate = 1
                            break
                    #--- no dupulication ---
                    if duplicate == -1:
                        randsynlist.append(row[:])
                break
        #--- no candidate left
        if (candlist == []):
            break

    #--- export data ---
    print "--- synapse list ---"
    for row in randsynlist:
        synParam(row, ratio, synType)
        writer.writerow(row)
        print row
    #--- reset lists (global access) ---
    candlist = []
    randsynlist = []
    #--- file close ---
    fw.flush()
    fw.close()

#----------------------------------
# Chemical and Electrical synapse Ratio
#
def synParam(elem, ratio, synType):
    #---  ---
    tmp = random()
    if(tmp < ratio):
        elem.append('Chemi')
        chemiType(elem, synType)
    else:
        elem.append('Gap')
        gapType(elem, synType)
    return

#----------------------------------
# Synapse type (Chemi)
#
def chemiType(elem, synType):
    if(elem[7] != 'Chemi'):
        return
    elem.append(synType)

#----------------------------------
# Synapse type (Gap)
#
def gapType(elem, synType):
    if(elem[7] != 'Gap'):
        return
    elem.append('nil')

#----------------------------------
# randSynExe
#
# row = (#0, #1, #2, #3, #4, #5, #6, #7, #8)
# #0. pre  synapse
# #1. post synapse
# #2. distance (square)
# #3. volume of pre compartment
# #4. volume of post compartment
# #5. synapse strength (how many synapses are there)
# #6. probability
# #7. synapse type (Chemical or Electrical)
# #8. synapse type 
#

def randSynExe(N = 20, dist_threshold = 1000, ratio = 0.8, synType = "GABA", filename = "synList.csv"):
    readDataCsv()
    probability(dist_threshold)
    randSynList(N, ratio, synType, filename)



#----------------------------------
# csv2list
#
def csv2list(synlistname, outputfile = 'output.txt'):
    i = 0

    #--- open file ---
    fh = open(synlistname, 'rb')
    fw = open(outputfile, 'wb')
    reader = csv.reader(fh)

    #--- regular expression ---
    pattern = re.compile(r'CellSwc\[(\d+)\]\.Dend\[(\d+)\]')

    #--- first row / second and later row ---
    for row in reader:
        if (i == 0):
            fprintrow(row, pattern, fw, 1)
        else:
            fprintrow(row, pattern, fw, 0)
        i += 1

    fh.close()
    fw.close()

    
def fprintrow(row, pattern, fw, first = 0):
    tmplist = []

    #--- pre ---
    tmpstr = row[0]
    result = pattern.search(tmpstr)
    tmpsearch = result.groups()
    tmplist.append(tmpsearch[1])
    pre = tmpsearch[0]

    #--- post ---
    tmpstr = row[1]
    result = pattern.search(tmpstr)
    tmpsearch = result.groups()
    tmplist.append(tmpsearch[1])
    post = tmpsearch[0]    

    #--- synapse type ---
    tmplist.append(row[7])
    tmplist.append(row[8])
    #--- synapse strength ---
    tmplist.append(row[5])    

    if(first):
        fw.write('CellSwc[' + str(pre) + '] ')
        fw.write('CellSwc[' + str(post) + ']\n')
        fw.write('----------\n')

    for elem in tmplist:
        fw.write(elem)
        fw.write(' ')
    fw.write('\n')
           
    
def mklists():
    #--- open connect list---
    confile = open('list/connectlist.csv', 'rb')
    exehocfile = open('importSynList.hoc', 'wb')
    conReader = csv.reader(confile)
    dirlist = os.listdir("./synlist/")

    i = 0

    exehocfile.write('//this file is created by csv2list.py \n')
    exehocfile.write('//------------------------------------\n')


    for elem in dirlist:
        if('.csv' in elem):
            print elem
            i += 1
            inputname = ('synlist/'+ elem)
            outputname = ('synlist/output_' + str(i) + '.txt')
            csv2list(inputname, outputname)
            exehocfile.write('readsl("' + outputname + '")')
            exehocfile.write('\n')

#-----------------------
# Main
# 
#-----------------------

def main(process = 1, myid = 0):

    #--- rand ---
    seed()

    #--- make cells ---
    h('objref pc')
    h('pc = new ParallelContext()')
    h.load_file("ParaCellSwc.hoc")
    h.load_file("getdistance.hoc")
    h.p_mkcells("list/neuronlist.txt", 1) # make cells in serial mode


    if(os.path.exists("synlist") == False):
        os.mkdir("synlist")

    i = 0
    fcon = open(connectfile, 'rb')
    reader = csv.reader(fcon)
    for row in reader:
        i += 1

        if( i % process == myid):
        #--- get synapse & rename file ---
            csvfilename = 'synlist/synList' + str(i) + '.csv'

            pre = row[0]
            post = row[1]

            getcandidate_exe = 'getcandidate(' + pre + ', ' + post + ')'
            h(getcandidate_exe)

            randSynExe(N = int(row[2]), dist_threshold = 1000, ratio = int(row[3])*0.01, synType = row[4], filename = csvfilename)

    fcon.flush()
    fcon.close()
    mklists()


#--- execution ---
if (len(sys.argv) > 1):
    main(int(sys.argv[1]),int(sys.argv[2]))
else:
    main()


