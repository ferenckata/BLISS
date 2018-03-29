# --------------- Dependencies ---------------
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# --------------- Inputs ---------------
indir = sys.argv[1] # where the original outputs (ending with "UMI.bed" can be found)
runid = sys.argv[2] # the ID of the run

# --------------- Output ---------------
pdfname = indir + runid + '.pdf'

# --------------- RUN ---------------
# currently with 5 colors for 5 datasets, colors are from the viridis scale
vir = ['#440154FF','#404788FF','#287D8EFF','#95D840FF','#FDE725FF']

cl = -1
for filename in os.listdir(indir):
    if filename.endswith('UMI.bed'):
        cl += 1
        print(filename)
        breaks = {}
        with open(indir + filename,'r') as bf:
            for line in bf:
                dsbcount = int(line.strip().split()[3])
                if dsbcount not in breaks.keys():
                    breaks[dsbcount] = 1
                else:
                    breaks[dsbcount] += 1
        lists = sorted(breaks.items())
        x,y = zip(*lists)
        plt.scatter(x,y,c=vir[cl],marker='.', label=filename.split('_')[1])
plt.legend(bbox_to_anchor=(0.75, 0.9), loc=2, borderaxespad=0.)
plt.yscale('log')
plt.xscale('log')
plt.xlim(0.7,)
plt.ylim(0.5,)
plt.title(runid)
plt.ylabel('number of locations (log)')
plt.xlabel('number of breaks (log)')
plt.savefig(pdfname)
