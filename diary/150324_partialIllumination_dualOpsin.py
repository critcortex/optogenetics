
import Neuron
from run_experiments import ExpSetup
import run_analysis
import NeuroTools.stgen as ntst 
import numpy as np

expbase = '150324_partialIllumination_dualOpsin'

"""
history:

- ran with no factor specified in exp. description hence no factor implies equal irradiance i.e. factor=1. 
- reran with factor = 0.7. Factor specified in exp. description

Next:
- can't see any difference in different gradients. Suspect this is due to the fact that we're too far from Vrest i.e. too much illumination 
Aim: extend so that background input i.e. in vivo like input (frozen input)
- Kinda of surprising that there was no difference as the gradient changed 

25/03/15 12.20pm: just finished investigating to check that gradient was working. It wasn't. The bug was in opsin.py and we weren't pulling projection out
correctly, therefore it was breaking in the try block --> no gradient calculated.
Have fixed this, and will rerun yesterday's run to check that there is a difference 

"""


stg = ntst.StGen()
es = ExpSetup()


tstop = 1500
light_dur = 1100
light_on = 200


L5PC_areas = ['soma', 'axon', 'apic','dend']

TOTAL_NUM_AREAS = {'L5PC': 4 }
areas = {'L5PC': L5PC_areas}

celltype = 'L5PC' 

gradients = np.arange(0.0,0.00101,0.0001)

def setup_test_chronicPulse(damp_chr,damp_hr,irr,factor,exp_desc):
 
    pp = {}
    # neuron model and params
    pp['cell'] = ['Neuron',celltype]
    pp['cell_params'] = {}
    # opsin
    chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':0,  'n_pulses':1}
    chrdict.update(damp_chr)
    hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':0,  'n_pulses':1}
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
    irrs = [0.01,0.02,0.05]
    factors = [0.1,0.2,0.5]
    #gradients = np.arange(0.0,0.00101,0.0002)#1)
    
    count=0
    for irr in irrs:
        for fac in factors:
            for ch_gra in gradients:
                for nphr_gra in gradients:
                    damp_chr = {'irr_gradient':ch_gra, 'irr_surface':1.0, 'projection': 'y'}
                    damp_hr = {'irr_gradient':nphr_gra, 'irr_surface':1.0, 'projection': 'y'}
                    setup_test_chronicPulse(damp_chr, damp_hr, irr, fac,'_irr%.2f_factor%g_ChRgrad%.4f_NpHRgrad%.4f'%(irr,fac,ch_gra,nphr_gra))
                    count+=1
                    
    print("%g jobs submitted"%count)
    
    
def analyse():
    
    af = run_analysis.AnalyticFrame()
    af.update_params({'label_format':'_irr%.2f_factor%g_ChRgrad%.4f_NpHRgrad%.4f_NpHR_whole_ChR_whole'})
    af.update_params({'tstart':light_on,'tstop':light_on+light_dur,
                                  'tstart_bg': 0,'tstop_bg':light_on,
                                  'tstart_post':light_on+light_dur,'tstop_post':tstop})
    irrs = [0.01,0.02,0.05]
    factors = [0.1,0.2,0.5]
    #gradients = np.arange(0.0,0.00101,0.0002)
    for irr in irrs:
        for fac in factors:
            for ch_gra in gradients:
                #for nphr_gra in gradients:
                exp_comp_list = [['_irr%.2f'%irr+'_factor%g'%fac+'_ChRgrad%.4f'%ch_gra+'_NpHRgrad%.4f_NpHR_whole_ChR_whole','..']]
                expss = [ec[0] for ec in exp_comp_list]
                explabels = [ec[1] for ec in exp_comp_list]
                print explabels
    
                af.populate_expset(expbase,expss,explabels, [gradients])

                af.submenu_extractSpikes()
                
                af.submenu_runFI()
                for exp in af.experimentset:
                    exp.calculate_responses('FI')
                    exp.calculate_responses('FI_bg')
                    exp.calculate_responses('FI_post')
#     
                af.submenu_save()
    af.submenu_print()


def plot_gradient_matrices():
    """
    Plot how FR changes for ChR2 gradient vs NpHR gradient for partial illumination
    """
    
    irrs = [0.01,0.02,0.05]
    factors = [0.1,0.2,0.5]
    #gradients = np.arange(0.0,0.00101,0.0002)
    for irr in irrs:
        for fac in factors:
            af = run_analysis.AnalyticFrame()
            af.update_params({'label_format':'_irr%.2f_factor%g_ChRgrad%.4f_NpHRgrad%.4f_NpHR_whole_ChR_whole'})
            af.update_params({'tstart':light_on,'tstop':light_on+light_dur,
                                  'tstart_bg': 0,'tstop_bg':light_on,
                                  'tstart_post':light_on+light_dur,'tstop_post':tstop})
            
            exp_comp_list = [['_irr%.2f'%irr+'_factor%g'%fac+'_ChRgrad%.4f'+'_NpHRgrad%.4f_NpHR_whole_ChR_whole','..']]
            expss = [ec[0] for ec in exp_comp_list]
            explabels = [ec[1] for ec in exp_comp_list]
            af.populate_expset(expbase,expss,explabels, [gradients,gradients])
            af.submenu_load()
            af.submenu_print()
            labels = ['%.1f'%(g*1000) for g in gradients]
            af.submenu_plot(11, 'FI_gain_irr%.2f_factor%g_varyGradients'%(irr,fac),axestitles=['NpHR gradient (a.u.)','ChR2 gradient (a.u.)'],axeslabels=[labels,labels],plottitle='Graded illumination (irr=%.2f, E:I opsin factor=%.1f)'%(irr,fac))
            #af.submenu_plot(0, expbase+'FI_gain_irr%.2f_factor%g_varyGradients'%(irr,fac))

                

        
#run_irr_gradient()
#analyse()
plot_gradient_matrices()