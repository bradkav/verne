import numpy as np

def ReadParam(filename, param_name):
    f = open(filename)
    for line in f:
        vals = np.array(map(str, line.split()))
        if (len(vals) > 0):
            if (vals[0] == param_name):
                return vals[1]


#f=open("Xenon.txt")
#for line in f:
#    vals = np.array(map(str, line.split()))
#    print vals
#for line in f:
#    print line[0]