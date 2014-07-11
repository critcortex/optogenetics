import run_experiments as re
import run_stimulation as rs

"""

pp = {}
pp['cell'] = ['Neuron', 'FractalNeuron']

es = re.ExpSetup()
dp = es.get_default_params()
dp.update(pp)

NE = rs.NeuronExperiment()
NE.params.update(dp)
NE.setup()
NE.printstatus('Setting up')
NE.setup()
NE.set_experiment_type()
NE.printstatus('Setting stimulus')
NE.set_stimulus()
NE.new_set_stimulus()
NE.printstatus('Setting optogenetics')
#NE.add_optogenetics()
#NE.add_optogenetics_old(NE.params['opdict'])
NE.printstatus('Setting recording')
NE.setup_record()
NE.printstatus('Simulating')
#NE.simulate_exp()
NE._run_simulation(NE.cell)
NE.printstatus('Saving data')
NE.save_data()
NE.printstatus('Running plots')

"""

"""
from  Neuron2013lib import BipolarNeuron
bn = BipolarNeuron()
NE = rs.NeuronExperiment()
NE._run_simulation(bn)

"""

"""
from Neuron import BipolarNeuron
bn = BipolarNeuron()
NE = rs.NeuronExperiment()
NE._run_simulation(bn)
"""

from Neuron import SimpleBipolar, FractalNeuron
#bn = SimpleBipolar()
bn = FractalNeuron()
print bn.xmid
gg = bn.select_segments_bydistance('dend0')
#print len(gg[0])
dd = bn.get_section_byname('dend0_0','dend0')
print dd, dd[0][1].name()
NE = rs.NeuronExperiment()
NE._run_simulation(bn)
