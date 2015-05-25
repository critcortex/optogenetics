
import Neuron
from run_experiments import ExpSetup
import run_analysis
import NeuroTools.stgen as ntst 
import numpy as np

expbase = '150518_attentionGating'

"""
Purpose:

Model for having input of three different types to Pyramidal cell (L5PC) neuron:
- feedforward (input)
- recurrent 
- feedback (attention)
and for:
- keeping input rate constant
- varying recurrent and feedback rates
- with and without optogenetics

Intent: that the change in gain modulation will be feed back into changes 
in recurrent population -- which will in turn modulate this cell

To consider:
-- how many inputs for basal vs apic - correct for 

--- TODO: move to include inhibitory inputs (10% of each)

"""

stg = ntst.StGen()
es = ExpSetup()


tstop = 2200
light_dur = 1000
light_on = 200


L5PC_areas = ['soma', 'apic','dend'] # NB: we don't include 'axon'

TOTAL_NUM_AREAS = {'L5PC': 3 }
areas = {'L5PC': L5PC_areas}

celltype = 'L5PC' 


gradients = np.arange(0.0,0.00101,0.0001)

num_sites = 500 # total number of connections
num_ff = int(0.1*num_sites) # which will all be in the apical dendrites
num_fb = int(0.1*num_sites) # which will all be in the basal dendrites
num_rr = int(0.8*num_sites) # which will be split between the apical and dendrites
# how we split our recurrent connections between apical and basal
rr_apic = int(0.5*num_rr)
rr_basal = num_rr - rr_apic
# so then the total number of synaptic sites on the apical and basal dendrites are:  
num_apic = num_ff + rr_apic
num_basal = num_fb + rr_basal

J = 0.5


def get_spiketimes(rate,num_inputs):
    spikes = []
    for i in range(num_inputs):
        ss = stg.poisson_generator(rate,t_stop=tstop,array=True)
        spikes.append(ss)
    #print 'Spikes = ',spikes
    return spikes


def _run_simulation(irr,factor,freq_ff,freq_rr,freq_fb,adist,bdist):
    """
    @params
        irr       irradiance
        factor    NpHR factor
    """
    
    pp = {}
    pp['cell'] = ['Neuron',celltype]
    pp['cell_params'] = {}
    pp['tstart'] = 0
    pp['tstop'] = tstop
    
    # opsin
    chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':0,  'n_pulses':1}
    hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':0,  'n_pulses':1}
    
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
    pp['experiment_type'] = 'opsinonly' #TODO change
    pp['savedata'] = True 
    
    
    #
    pp['stim_spiketrains'] = True                            
    labels = ['stim%g'%i for i in range(num_sites)]
    tstims = get_spiketimes(freq_ff,num_ff) + get_spiketimes(freq_rr,num_rr) + get_spiketimes(freq_fb,num_fb)
    pp['spiketrains'] = [{'tstims': tstims,  'locations': labels, 'weights':np.ones(num_sites)*J,  'el': 0.02}]
                                
    
    # mark locations for PSP throughout dendritic trees, along with soma
    pp['mark_loc'] = {}
    pp['mark_loc']['names'] = ['mysoma']+labels
    pp['mark_loc']['sections'] = ['soma']+['apic']*num_apic+['dend']*num_basal
    pp['mark_loc']['ids'] = [(0,0.5)] + [('select_section_posn_bydistance',{'sectionarea':'apic','mindist':adist[0],'maxdist':adist[1]}) for i in range(num_apic)] + [('select_section_posn_bydistance',{'sectionarea':'dend','mindist':bdist[0],'maxdist':bdist[1]}) for i in range(num_basal)] 
                                
                                
    # set up recordings
    # - record voltage from soma
    pp['record_loc'] = {}
    pp['record_loc']['v'] = ['mysoma']
    
    pp['plot'] = { 'voltages': [ ['mysoma', 'v','k-'] ]} 
    
    pp.update({'expname':expbase,
                   'description':'_irr%.3f_factor%.2f_sites%g_ff%g_rr%g_fb%g_J%.1f'%(irr,factor,num_sites,freq_ff,freq_rr,freq_fb,J)})
    es.run_single_experiment(expbase, 'cluster', pp)




ff_rates = [10,20,50,60,80]
scan_rates = range(20,81,5) # for fb and rr
irr = 0.001
factor = 0.6

for ffr in ff_rates:
    for rrr in scan_rates:
        for fbr in scan_rates:
            _run_simulation(irr,factor,ffr,rrr,fbr,[100,-1],[100,-1])   







