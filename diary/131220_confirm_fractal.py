import Neuron
from run_experiments import ExpSetup

expbase = '131220_confirmFractal'


num_base = range(1,5)
num_child = range(1,5)
num_levels = range(1,5)
"""
num_base = [1]
num_child = [1]
num_levels = [3]

"""
es = ExpSetup()

for nb in num_base:
    for nc in num_child:
        for nl in num_levels:
            pp = {}
            pp['cell'] = ['Neuron', 'FractalNeuron']
            pp['cell_params'] = {'num_base': nb,
                  'num_split': nc,
                  'num_levels':nl}
            
            pp.update({'expname':expbase,
                       'description':'_nb%g_ns%g_nl%g'%(nb,nc,nl)})
            
            #print 'going to run on cluster'
            es.run_single_experiment(expbase, 'local', pp)
            
            