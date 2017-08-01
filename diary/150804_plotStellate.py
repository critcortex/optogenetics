import LFPy
import Neuron
import os
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
import pylab
import numpy as np

"""

For some reason, this wouldn't work when pasted in a iPython notebook ... 

"""

#os.chdir('/home/sjarvis1/workspace/co_optogenetics')
cellStell = Neuron.SHStellate()

rotation = {'x' : 1.7, 'y' : 0, 'z' : 0}


cellStell.set_rotation(**rotation)
zips = []
for x, z in cellStell.get_idx_polygons():
    zips.append(zip(x, z))

polycol = PolyCollection(zips,
                          linewidths=2,
                         edgecolors='black',
                         facecolors='black')

fig = pylab.figure(figsize=(15,15))
ax = fig.add_subplot(111)

ax.add_collection(polycol)
ax.axis(ax.axis('equal'))

pylab.savefig("Fig3A_plotStellate.png",dpi=300)
pylab.savefig("Fig3A_plotStellate.svg",dpi=300)

