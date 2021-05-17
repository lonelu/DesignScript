import numpy as np

def readCA(file, coorType = 1):
    f = open(file, 'r')
    lines = []
    if coorType == 1:
        for line in f:
            if 'CA' in line:
                lines.append([float(n) for n in line.split()[6:9]])
    elif coorType == 0:
        for line in f:
            lines.append(([float(n) for n in line.split()]))
    return np.array(lines)
                

