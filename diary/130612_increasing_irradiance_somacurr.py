import run_experiments
import run_analysis
import numpy as np
import sys

expbase = '130612_increasing_irradiance_somacurr'

current_amps = np.arange(0.,2.,.1) # 2.,.1
irradiances = np.arange(0,10.1,0.5)

runon = 'cluster'

explist = [['apical','none'],['whole','none'],['none','apical'],['none','whole'],['apical','apical'],['whole','whole'],['none','none']]

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
            for amp in current_amps:
                params['experiment_type'] = 'BAC'
                params['iclamp_amp']= amp
                params['iclamp_start']= curr_on
                params['iclamp_duration'] = curr_dur
                opsin_expression = [5e-4,irr,light_on,light_dur,light_dur,600,1] 
                params.update({'NpHR_times':opsin_expression,'ChR_times':opsin_expression,'description':'_irr%.1f_Isoma%.1f'%(irr,amp)})
                es.run_single_experiment(expbase, runon, params)


def analyse():

    af = run_analysis.AnalyticFrame()
    af.update_params({'tstart':light_on,'tstop':light_on+light_dur})
    
    exp_comp_list = [['_irr%.1f_Isoma%.1f_'+'NpHR_%s_ChR_%s'%(exp[1],exp[0]),'ChR %s, NpHR %s'%(exp[0],exp[1])] for exp in explist ]
    print exp_comp_list
    
    
    
    '''
    exp_comp_list = [ ['_Isoma%.1f_NpHR_whole_ChR_none','ChR2 none,NpHR whole'],
                      ['_Isoma%.1f_NpHR_apical_ChR_none','ChR2 none,NpHR apical'],
                      ['_Isoma%.1f_NpHR_none_ChR_apical','ChR2 apical,NpHR none'],
                      ['_Isoma%.1f_NpHR_none_ChR_whole','ChR2 whole,NpHR none'],
                      #['_Isoma%.1f_NpHR_apical_ChR_basal','ChR2 basal,NpHR apical'],
                      #['_Isoma%.1f_NpHR_basal_ChR_apical','ChR2 apical,NpHR basal'],
                      ['_Isoma%.1f_NpHR_whole_ChR_whole','ChR2 whole,NpHR whole'],
                      ['_Isoma%.1f_NpHR_apical_ChR_apical','ChR2 apical,NpHR apical'], 
                      ['_Isoma%.1f_NpHR_none_ChR_none','ChR2 none,NpHR none'] ]
    '''
                 
    expss = [ec[0] for ec in exp_comp_list]
    explabels = [ec[1] for ec in exp_comp_list]
    af.populate_expset(expbase,expss,explabels,[irradiances,current_amps])
    af.submenu_load()
    #return af
    af.submenu_plot(4, 'test_vmax30')
    
    
    
    #af.submenu_plot(0, expbase+'FI_gain')
    #return 
    #af.run_analysis_menu()
    
    
    #af.submenu_runFI()
    #af.submenu_plot(3, expbase+'_FI_compare')
    #af.submenu_Iphoto()
    #af.submenu_print()
    #af.submenu_save()
    
    
def plot_gain():
    
    for irr in [1,5,10]:
    
        af = run_analysis.AnalyticFrame()
        af.update_params({'tstart':light_on,'tstop':light_on+light_dur})
        
        exp_comp_list = [['_irr%.1f'%irr+'_Isoma%.1f_'+'NpHR_%s_ChR_%s'%(exp[1],exp[0]),'ChR %s, NpHR %s'%(exp[0],exp[1])] for exp in explist ]
        print exp_comp_list
        
        expss = [ec[0] for ec in exp_comp_list]
        explabels = [ec[1] for ec in exp_comp_list]
        af.populate_expset(expbase,expss,explabels,[current_amps])
        af.submenu_load()
        
        af.submenu_plot(5, expbase+'FI_gain_irr%.1f_varyFactor_exp%s%s_'%(irr,exp[0],exp[1]))
        #af.submenu_plot(0, expbase+'FI_gain_irr%.1f'%irr)
    
        
        
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        pass
    elif sys.argv[1] == 'run':
        run_experiments()        
    elif sys.argv[1] == 'analyse':
        analyse()     
    elif sys.argv[1] == 'plotgain':
        plot_gain()
    