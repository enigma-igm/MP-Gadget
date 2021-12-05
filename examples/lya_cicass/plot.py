import numpy as np
import matplotlib.pyplot as plt
import bigfile

f = bigfile.File('output/PART_008')

# particle IDs: 0 = baryons, 1 = dark matter
dm = f['1/Position'][:]
bar = f['0/Position'][:]

fig,ax = plt.subplots(1,1,figsize=(6.5,6))
mask = (dm[:,2] > 490.0) & (dm[:,2] < 510.0)
plt.scatter(dm[mask,0],dm[mask,1],c='b',s=20,alpha=0.3,lw=0)
mask = (bar[:,2] > 490.0) & (bar[:,2] < 510.0)
plt.scatter(bar[mask,0],bar[mask,1],c='r',s=20,alpha=0.3,lw=0)
plt.xlim(0,1000)
plt.ylim(0,1000)
plt.tight_layout()
plt.show()

