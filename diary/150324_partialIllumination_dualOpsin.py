
import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np

expbase = '150324_partialIllumination_dualOpsin'


stg = ntst.StGen()
es = ExpSetup()


tstop = 1500
light_dur = 1100
light_on = 200


L5PC_areas = ['soma', 'axon', 'apic','dend']

TOTAL_NUM_AREAS = {'L5PC': 4 }
areas = {'L5PC': L5PC_areas}

celltype = 'L5PC' 



def setup_test_chronicPulse(damp_chr,damp_hr,irr,exp_desc):
 
    pp = {}
    # neuron model and params
    pp['cell'] = ['Neuron',celltype]
    pp['cell_params'] = {}
    # opsin
    chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':0,  'n_pulses':1}
    chrdict.update(damp_chr)
    hrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':0,  'n_pulses':1}
    hrdict.update(damp_hr)
    
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
    pp['savedata'] = True 
    
    pp['mark_loc'] = {}
    pp['mark_loc']['names'] = ['mysoma']
    pp['mark_loc']['sections'] = ['soma']
    pp['mark_loc']['ids'] = [(0,0.5)]

    pp['record_loc'] = {}
    pp['record_loc']['v'] = ['mysoma']
    vplots_soma = [['mysoma','v','k']]
    pp['plot'] = {1:vplots_soma}
    
    
    
    pp['tstart'] = 0
    pp['tstop'] = tstop
    
    pp['num_threads'] = 1
                               
    pp.update({'expname':expbase,
               'description':exp_desc})

    es.run_single_experiment(expbase, 'missing', pp)
    

def run_irr_gradient():
    irrs = [0.02,0.05]
    gradients = np.arange(0.0,0.00101,0.0001)
    count=0
    for irr in irrs:
        for ch_gra in gradients:
            for nphr_gra in gradients:
                damp_chr = {'irr_gradient':ch_gra, 'irr_surface':1.0, 'projection': 'y'}
                damp_hr = {'irr_gradient':nphr_gra, 'irr_surface':1.0, 'projection': 'y'}
                setup_test_chronicPulse(damp_chr, damp_hr, irr, '_irr%.2f_ChRgrad%.4f_NpHRgrad%.4f'%(irr,ch_gra,nphr_gra))
                count+=1
    print("%g jobs submitted"%count)
            
        
        
run_irr_gradient()