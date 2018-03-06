
# coding: utf-8


import matplotlib.pyplot as plt
import numpy as np
import sys


file = sys.argv[1]
name = file.split('.')[0]
pngname = name + '.png'

breaks = {}
with open(file,'r') as bf:
    for line in bf:
        dsbcount = int(line.strip().split()[3])
        if dsbcount not in breaks.keys():
            breaks[dsbcount] = 1
        else:
            breaks[dsbcount] += 1


lists = sorted(breaks.items())
x,y = zip(*lists)

plt.plot(x,y,'r')
plt.yscale('log')
plt.xscale('log')
plt.title(name)
plt.ylabel('number of locations (log)')
plt.xlabel('number of breaks (log)')
plt.savefig(pngname)


