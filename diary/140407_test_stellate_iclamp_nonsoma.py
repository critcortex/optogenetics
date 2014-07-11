
from neuron import h
import Neuron
import numpy as np
import pylab

h.load_file('stdlib.hoc', 'String') 
h.load_file('stdrun.hoc')
h.load_file("import3d.hoc")

cell = Neuron.Stellate()
#cell = Neuron.L5PC()

h.refrac_SpikeOut = 5
h.vrefrac_SpikeOut = -60
h.thresh_SpikeOut = -50

spkout = h.SpikeOut(0.5,sec=cell.soma)

#clamp_site = cell.select_section_posn_bydistance('apic')[0]
# for dendrites @ hill
clamp_site = cell.select_section_posn_bydistance('myelin')[2]
# for soma
#clamp_site = (0, cell.soma,0.5)
stim = h.IClamp()
stim.loc(clamp_site[2], sec=clamp_site[1])
stim.dur = 400
setattr(stim, 'del',200)
stim.amp = 10.


time = h.Vector()
vmsoma = h.Vector()

time.record(h._ref_t)
vmsoma.record (cell.soma(0.5)._ref_v)


h.init()
h.tstop = 700
h.run()

vms = np.array(vmsoma)
tt = np.array(time)
pylab.plot(tt,vms)
pylab.savefig('test_iclamp%g.png'%stim.amp)


