
from neuron import h
import Neuron
import numpy as np
import pylab

h.load_file('stdlib.hoc', 'String') 
h.load_file('stdrun.hoc')
h.load_file("import3d.hoc")

"""
pp = {'num_split': 2, 
      'dend0': {'mechanisms': ['pas_dend', 'hh'], 'Ra': 100, 'cm': 2.0}, 
       'dend1': {'mechanisms': ['pas_dend', 'hh'], 'Ra': 100, 'cm': 2.0}, 
       'soma': {'mechanisms': ['pas_soma', 'hh'], 'cm': 1.0, }, 
       'mechanisms': {'pas_dend': ('pas', {'e': -70}), 'pas_soma': ('pas', {'e': -70})}, 
       'num_base': 2, 
       'hh': ('hh',), 
       'defaultlevel': {'mechanisms': ['pas_dend', 'hh'], 'Ra': 100, 'cm': 2.0}, 
       'num_levels': 3}
"""

pp = {'num_split': 2, 
      'dend0': {'mechanisms': ['pas_dend', 'hh'],'diam':1.,'L':600}, 
       'dend1': {'mechanisms': ['pas_dend', 'hh'],'diam':1.,'L':600}, 
       'soma': {'mechanisms': ['pas_soma', 'hh'], 'diam':30,'L':30, 'Ra': 200, 'cm':1}, 
       'mechanisms': {'pas_dend': ('pas', {'e': -45,'g':0.0005}), 'pas_soma': ('pas', {'e': -65,'g':0.0005})}, 
       'num_base': 2, 
       'hh': ('hh',), 
       'defaultlevel': {'mechanisms': ['pas_dend', 'hh'],'diam':1.,'L':600}, 
       'num_levels': 3}


fn = Neuron.FractalNeuron(pp)


h.refrac_SpikeOut = 5
h.vrefrac_SpikeOut = -60
h.thresh_SpikeOut = -50

spkout = h.SpikeOut(0.5,sec=fn.soma)

stim = h.IClamp()
stim.loc(0.5, sec=fn.dend0)
stim.dur = 400
setattr(stim, 'del',200)
stim.amp = 19000.19


time = h.Vector()
vmsoma = h.Vector()

time.record(h._ref_t)
vmsoma.record (fn.soma(0.5)._ref_v)


h.init()
h.tstop = 700
h.run()

vms = np.array(vmsoma)
tt = np.array(time)
pylab.plot(tt,vms)
pylab.savefig('test_iclamp%g.png'%stim.amp)
