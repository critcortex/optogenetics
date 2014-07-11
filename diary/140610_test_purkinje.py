from neuron import h
import Neuron


###############################################
# Version 1: worked fine
h.load_file('stdlib.hoc') 
h.load_file('stdrun.hoc')
h.load_file("import3d.hoc")
h.xopen("purkinje.hoc")


cell = h.purkinje_miyasho
cell.soma

