#!/usr/bin/python


def findname():
    cellname = 'LAL-bilateral0'
    
    flist = open("list/neuronlist.txt")
    lines = flist.readlines()
    flist.close()
    
    count = -2
    
    for line in lines:
        line = line.rstrip()
        if line != '' and list(line)[0] != '#':
            count += 1
            
        if line.find(cellname) != -1:
            # print line
            return 'CellSwc[%d]' % count
    
name = findname()
print name
