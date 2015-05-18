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

def get_optdescript(irr,factor):
    opt_des = []
    if irr > 0 :
        opt_des.append('whole')
    else:
        opt_des.append('none')
        
    if irr > 0 and factor > 0:
        opt_des.append('whole')
    else:
        opt_des.append('none')
    return opt_des

# Purpose - redraw in vivo / in vitro IV curves but with better fit

expbase = '140722_cell_irr_ChRonly'
save_appendix = '_run141031'
               
factors = [0]
irrs = [ coeff*power for power in [0.001,0.01,0.1,1.] for coeff in [1,2,5] ]

light_on = 1050
light_dur = 1000
tstop = light_on + light_dur + 50 # as buffer
Js = [2.]
nsite_range = [80]
dist = [100,-1]
celltype = 'SHStellate'

freqs = range(20,150,10)

for f in factors:
    af = run_analysis.AnalyticFrame()
    
    af.update_params({'tstart':light_on,'tstop':light_on+light_dur,
                              'tstart_bg': 50,'tstop_bg':light_on,
                              'tstart_post':light_on+light_dur,'tstop_post':tstop})
    #optlog = get_optdescript(irr,1)
    exp_comp_list = [['irr%.3f_'%irr+'factor%.2f_'%f+'freq%g'+'_J%g_nsites%g_celltype%s'%(Js[0],nsite_range[0],celltype)+'_NpHR_%s_ChR_%s'%(get_optdescript(irr,f)[1],get_optdescript(irr,f)[0]),'=%.3f'%irr] for irr in irrs]
    
    print exp_comp_list
  
    expss = [ec[0] for ec in exp_comp_list]
    explabels = [ec[1] for ec in exp_comp_list]
    print explabels
    
    af.populate_expset(expbase,expss,explabels, [freqs])
        
    
    af.submenu_extractSpikes()
    
    
    af.submenu_runFI()
    for exp in af.experimentset:
        exp.calculate_responses('FI')
        exp.calculate_responses('FI_bg')
        exp.calculate_responses('FI_post')
    
    
    af.submenu_print()
    #af.submenu_save()
    
    #af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
    af.submenu_plot(5, expbase+'FI_gain_varyIrr_invivo'+save_appendix)
    af.submenu_plot(0, expbase+'FI_gain_varyIrr_invivo'+save_appendix)
    af.submenu_plot(10, expbase+'FI_gain_varyIrr_invivo'+save_appendix)

