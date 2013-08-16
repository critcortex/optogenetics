import run_experiments
import run_analysis
import numpy as np
import sys

def poly(x,a,b,c,d):
    print type(x)
    return a*x**3 + b*x**2 + c*x + d

expbase = '130709_increasing_uneq_irradiance_factor_somacurr'

current_amps = np.arange(0.1,2.,.1) # 2.,.1
current_amps = np.arange(-2.,0.01,0.1)
current_amps = np.arange(-2.,2.01,0.1)
#current_amps = [2.]

irradiances = np.arange(0,5.1,0.5)
irradiances = [1.0,5.0]
irradiances = [1.]
factors = [1.,2.,4.,0.5,0.25]
orig_factors = [2.,4.,8.,16.,32.]
factors = orig_factors + [1./x for x in orig_factors] + [1]

runon = 'cluster'

#explist = [['apical','none'],['whole','none']] ,['none','apical'],['none','whole'],['apical','apical'],['whole','whole'],['none','none']]
explist = [['apical','apical'],['whole','whole']] 
explist = [['apical','apical'],['whole','whole']] 


light_on = 700
curr_on = 0
curr_dur = 2500
light_dur = 1000
tstop = 2500


es = run_experiments.ExpSetup()
areas = es.get_areas_section()
params = {}

def run_experiments():
    for exparea in explist:
        params['ChR_areas'] = {exparea[0]:areas[exparea[0]]}
        params['NpHR_areas'] = {exparea[1]:areas[exparea[1]]}
        
        params['tstop'] = tstop

        for irr in irradiances:
            for factor in factors:
                for amp in current_amps:
                    params['experiment_type'] = 'BAC'
                    params['iclamp_amp']= amp
                    params['iclamp_start']= curr_on
                    params['iclamp_duration'] = curr_dur                    
                    opsin_expression_ch = [5e-4,irr,light_on,light_dur,light_dur,600,1] 
                    opsin_expression_np = [5e-4,irr*factor,light_on,light_dur,light_dur,600,1] 
                    params.update({'NpHR_times':opsin_expression_np,'ChR_times':opsin_expression_ch,'description':'_irr%.1f_factor%.2f_Isoma%.1f'%(irr,factor,amp)})
                    es.run_single_experiment(expbase, runon, params,checkdatexists=True)


def analyse():

    af = run_analysis.AnalyticFrame()
    af.update_params({'tstart':light_on,'tstop':light_on+light_dur})
    
    #irradiances = [0.5]
    #factors = [0.25] #[4.,0.25]
    irradiances = [1]
    
    exp_comp_list = [['_irr%.1f_factor%.2f_Isoma%.1f_'+'NpHR_%s_ChR_%s'%(exp[1],exp[0]),'ChR %s, NpHR %s'%(exp[0],exp[1])] for exp in explist ]
    print exp_comp_list
    
                 
    expss = [ec[0] for ec in exp_comp_list]
    explabels = [ec[1] for ec in exp_comp_list]
    af.populate_expset(expbase,expss,explabels,[irradiances,factors,current_amps])
    #af.submenu_load()
    #return af
    
    #af.run_analysis_menu()
    
    
    af.submenu_runFI()
    af.submenu_print()
    af.submenu_save()    
    
    
    #af.submenu_plot(3, expbase+'_FI_compare')
    #af.submenu_Iphoto()
    #af.submenu_print()
    #af.submenu_save()
    
    
def plot_gain():
    """
    for irr in [1,5]:
        for factor in factors:
            af = run_analysis.AnalyticFrame()
            af.update_params({'tstart':light_on,'tstop':light_on+light_dur})
            
            exp_comp_list = [['_irr%.1f_factor%.2f'%(irr,factor)+'_Isoma%.1f_'+'NpHR_%s_ChR_%s'%(exp[1],exp[0]),'ChR %s, NpHR %s'%(exp[0],exp[1])] for exp in explist ]
            print exp_comp_list
            
            expss = [ec[0] for ec in exp_comp_list]
            explabels = [ec[1] for ec in exp_comp_list]
            af.populate_expset(expbase,expss,explabels,[current_amps])
            
            af.submenu_load()
            af.submenu_print()
            af.submenu_plot(0, expbase+'FI_gain_irr%.1f_factor%.2f'%(irr,factor))
     """
     
    explist = [['whole','whole']] 
    explist = [['apical','apical'],['whole','whole'] ] 
    factors.sort()

    for irr in [1]: 
        for exp in explist:
            af = run_analysis.AnalyticFrame()
            af.update_params({'tstart':light_on,'tstop':light_on+light_dur})
            
            exp_comp_list = [['_irr%.1f'%(irr)+'_factor%.2f'%(factor)+'_Isoma%.1f_'+'NpHR_%s_ChR_%s'%(exp[1],exp[0]),'%g'%factor] for factor in factors]
            print exp_comp_list
            
            expss = [ec[0] for ec in exp_comp_list]
            explabels = [ec[1] for ec in exp_comp_list]
            af.populate_expset(expbase,expss,explabels,[current_amps])
            
            af.submenu_load()
            af.submenu_print()
            #af.submenu_plot(0, expbase+'FI_gain_irr%.1f_varyFactor_exp%s%s_'%(irr,exp[0],exp[1]))
            af.submenu_plot(5, expbase+'FI_gain_irr%.1f_varyFactor_exp%s%s_'%(irr,exp[0],exp[1]))






        
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        pass
    elif sys.argv[1] == 'run':
        run_experiments()        
    elif sys.argv[1] == 'analyse':
        analyse()  
    else:
        plot_gain()
              
    