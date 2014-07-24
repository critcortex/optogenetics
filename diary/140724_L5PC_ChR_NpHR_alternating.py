import sys
import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import run_stimulation
import pylab
import file_io as fio
import run_analysis
import NeuroTools.stgen as ntst 



expbase = '140724_L5PC_ChR_NpHR_alternating'


stg = ntst.StGen()
es = ExpSetup()
 
light_on = 100
light_dur = 5
tstop = 800

# Optogenetics parameters
factors = [0,1.]
irrs = [ coeff*power for power in [0.001,0.01,0.1,1.] for coeff in [1,2,5] ]

# In vitro parameters
iclamp_amps = np.arange(0.0,5.51,0.5)

# In vivo parameters
freqs = range(20,150,10)
freqs = range(15,150,10) + [150]
#freqs = range(15,151,5)
Js = [2.]
nsite_range = [80]
dist = [100,-1]




L5PC_areas = ['soma', 'axon', 'apic','dend']
#L5PC_areas = ['soma','dend']
SHStellate_areas = ['soma','dendrite','axon']

TOTAL_NUM_AREAS = {'L5PC': 4,
                   'SHStellate':2 }
areas = {'L5PC': L5PC_areas,
         'SHStellate':SHStellate_areas}

celltype = 'L5PC'


#######################
factors = [0.125,0.25,0.5,1.,2.,4.,8.]
factors = [0.125,0.25,0.5,0.75,1.,1.5,2.] # for vitro
factors = [0.125,0.25,0.5,0.75,1.] # for vivo

irrs = [0.01,0.02]

factor = 0.7
irr = 0.02
n_pulses = 20
ipi = 65 # interpulse interval
pulse_offset = (light_dur+ipi)/2


def get_spiketimes(rate,num_inputs):
    spikes = []
    for i in range(num_inputs):
        # TODO: am guessing there's a more efficient way of implementing this
        ss = stg.poisson_generator(rate,t_stop=tstop,array=True)
        spikes.append(ss)
    #print 'Spikes = ',spikes
    return spikes

def alternating_pulses():
 
    pp = {}
    # neuron model and params
    pp['cell'] = ['Neuron',celltype]
    pp['cell_params'] = {}
    
    # opsin
    chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':ipi,  'n_pulses':n_pulses}
    hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on+pulse_offset,'interpulse_interval':ipi,  'n_pulses':n_pulses}
    
    pp['opsindict'] = {}
    pp['opsindict']['ChR'] =  {}
    for area in areas[celltype]:
        pp['opsindict']['ChR'][area] = chrdict    
    pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
    
    pp['opsindict']['NpHR'] =  {}
    for area in areas[celltype]:
        pp['opsindict']['NpHR'][area] = hrdict
    pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
    
    # general settings 
    pp['experiment_type'] = 'opsinonly'
    pp['savedata'] = True # False #True
    
    pp['tstart'] = 0
    pp['tstop'] = tstop
    
    pp['mark_loc'] = {}
    pp['mark_loc']['names'] = ['mysoma']
    pp['mark_loc']['sections'] = ['soma']
    pp['mark_loc']['ids'] = [(0,0.5)]
                 
    
    pp['record_loc'] = {}
    pp['record_loc']['v'] = ['mysoma'] #+labels
    pp['record_loc']['ina'] = ['mysoma']
    pp['record_loc']['ik'] = ['mysoma']
    
    #print labels[0]
    #pp['spiketrains'] = [{'tstims': get_spiketimes(freq,nsites), 'locations': labels, 'weights':np.ones(nsites)*J,  'el': 0.1}]
    
    
    vplots_soma = [['mysoma','v','k']]

    #iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
    pp['plot'] = {1:vplots_soma} 
        
    
    pp['num_threads'] = 1
                               
    pp.update({'expname':expbase,
               'description':'irr%.3f_factor%.2f_ipi%g'%(irr,factor,ipi)})

    es.run_single_experiment(expbase, 'local', pp)
    


 
if __name__ == '__main__':
    alternating_pulses()
                                            
                           