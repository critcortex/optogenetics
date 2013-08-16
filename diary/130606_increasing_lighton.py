import run_experiments
import run_analysis
import numpy as np
import sys

expbase = '130606_incr_lighton'
expbase = '130624_incr_lighton'

DEFAULT_FREQ = [0.5,1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,20,30,35,40,45,50]
current_amps = np.arange(0.,2.01,.1)

runon = 'cluster'
n_pulses = 100
transient=200
pulsewidth=10
offset_phase=0
tstop_long = 2000
tstop_short=500
BAC_start = 100 # should be before transient for opsins

explist = [['apical','none'],['whole','none'],['none','apical'],['none','whole'],['apical','apical'],['whole','whole'],['none','none']]

light_on = 100
curr_on = 200
curr_dur = 500
light_off = 1000

# for 130624_incr_lighton
light_on = 100
curr_on = 200
curr_dur = 1000
light_off = 1300

es = run_experiments.ExpSetup()
areas = es.get_areas_section()
params = {}

def run_experiments():
    for exparea in explist:
        params['ChR_areas'] = {exparea[0]:areas[exparea[0]]}
        params['NpHR_areas'] = {exparea[1]:areas[exparea[1]]}
        
        params['tstop'] = light_off+100

        for amp in current_amps:
            params['experiment_type'] = 'BAC'
            params['iclamp_amp']= amp
            params['iclamp_start']= curr_on
            params['iclamp_duration'] = curr_dur
            opsin_expression = [5e-4,10,light_on,light_off,light_off,600,1] 
            params.update({'NpHR_times':opsin_expression,'ChR_times':opsin_expression,'description':'_Isoma%.1f'%amp})
            es.run_single_experiment(expbase, runon, params)


def analyse():

    af = run_analysis.AnalyticFrame()
    af.update_params({'tstop':curr_dur,'tstart':curr_on})
    
    
    exp_comp_list = [ ['_Isoma%.1f_NpHR_whole_ChR_none','ChR2 none,NpHR whole'],
                      ['_Isoma%.1f_NpHR_apical_ChR_none','ChR2 none,NpHR apical'],
                      ['_Isoma%.1f_NpHR_none_ChR_apical','ChR2 apical,NpHR none'],
                      ['_Isoma%.1f_NpHR_none_ChR_whole','ChR2 whole,NpHR none'],
                      #['_Isoma%.1f_NpHR_apical_ChR_basal','ChR2 basal,NpHR apical'],
                      #['_Isoma%.1f_NpHR_basal_ChR_apical','ChR2 apical,NpHR basal'],
                      ['_Isoma%.1f_NpHR_whole_ChR_whole','ChR2 whole,NpHR whole'],
                      ['_Isoma%.1f_NpHR_apical_ChR_apical','ChR2 apical,NpHR apical'], 
                      ['_Isoma%.1f_NpHR_none_ChR_none','ChR2 none,NpHR none'] ]
                 
    expss = [ec[0] for ec in exp_comp_list]
    explabels = [ec[1] for ec in exp_comp_list]
    af.populate_expset(expbase,expss,explabels,current_amps)
    af.submenu_runFI()
    af.submenu_Iphoto()
    af.submenu_print()
    af.submenu_save()
    af.submenu_plot(0, expbase+'_FI_test2')
    af.submenu_plot(1, expbase+'_voltageclamp_test2')
    af.submenu_plot(2, expbase+'_FI_current_test2')
    
    
    
        
if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run_experiments()        
    elif sys.argv[1] == 'analyse':
        analyse()        
