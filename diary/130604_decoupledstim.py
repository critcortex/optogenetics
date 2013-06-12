import run_experiments
import run_analysis
import numpy as np
import sys

DEFAULT_FREQ = [0.5,1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,20,30,35,40,45,50]
epsp_amps = np.arange(1,10.5,1)

runon = 'cluster'
n_pulses = 100
transient=200
pulsewidth=10
offset_phase=0
tstop_long = 2000
tstop_short=500
BAC_start = 100 # should be before transient for opsins

explist = [['apical','none'],['whole','none'],['none','apical'],['none','whole'],['apical','basal'],['basal','apical']] 



es = run_experiments.ExpSetup()
areas = es.get_areas_section()
params = {}

def run_experiments():
    for exparea in explist:
        params['ChR_areas'] = {exparea[0]:areas[exparea[0]]}
        params['NpHR_areas'] = {exparea[1]:areas[exparea[1]]}
        
        params['tstop'] = tstop_long
        '''
        params['experiment_type'] = 'opsinonly'
        es.run_frequency_range('130606_decoupledapical_t2k',freqs=DEFAULT_FREQ,pulsewidth=pulsewidth,transient=transient,n_pulses=n_pulses,params=params,runon=runon,offset_phase=offset_phase)
        '''
        '''
        params['experiment_type'] = 'BAC'
        params['iclamp_start']= BAC_start
        params['iclamp_duration'] = tstop_long - BAC_start
        es.run_frequency_range('130606_decoupledapical_BAC_t2k',freqs=DEFAULT_FREQ,pulsewidth=pulsewidth,transient=transient,n_pulses=n_pulses,params=params,runon=runon,offset_phase=offset_phase)
        '''
        '''
        params['experiment_type'] = 'BAC'
        params['iclamp_start']= BAC_start
        params['iclamp_duration'] = tstop_long - BAC_start
        es.run_frequency_range('130606_coupledapical_BAC_t2k',freqs=DEFAULT_FREQ,pulsewidth=pulsewidth,transient=transient,n_pulses=n_pulses,params=params,runon=runon,offset_phase=0.5)
    '''
        '''
        for epsp in epsp_amps:
            params['tstop'] = tstop_short
            params['experiment_type'] = 'BAC'
            params['EPSP_transient'] = 200
            params['EPSP_amp'] = epsp
            params['iclamp_amp']= 0
            params['iclamp_duration'] = 0
            opsin_expression = [5e-4,10,100,tstop_short,tstop_short,600,1] 
            params.update({'NpHR_times':opsin_expression,'ChR_times':opsin_expression,'description':'_Jex%.1f'%epsp})
            es.run_single_experiment('130606_decoupled_singlestim', runon, params)
        '''
        
def get_baseline():
    DEFAULT_FREQ = [1]
    explist = [['none','none']] 
    #use these values when calculating deviation

        
def analyse1():
    expbase = '130606_decoupledapical_BAC_t2k'
    af = run_analysis.AnalyticFrame()
    af.update_params({'tstop':2000,'tstart':180,'baseline':'_freq1_pw10_NpHR_none_ChR_none'})
    
    '''
    exp_comp_list = [ ['_freq%s_pw10_NpHR_whole_ChR_none','ChR2 none,NpHR whole'],
                      ['_freq%s_pw10_NpHR_apical_ChR_none','ChR2 none,NpHR apical'],
                      ['_freq%s_pw10_NpHR_none_ChR_apical','ChR2 apical,NpHR none'],
                      ['_freq%s_pw10_NpHR_none_ChR_whole','ChR2 whole,NpHR none'],
                      ['_freq%s_pw10_NpHR_apical_ChR_basal','ChR2 basal,NpHR apical'],
                      ['_freq%s_pw10_NpHR_basal_ChR_apical','ChR2 apical,NpHR basal'] ]
                 
    expss = [ec[0] for ec in exp_comp_list]
    explabels = [ec[1] for ec in exp_comp_list]
    af.populate_expset(expbase,expss,explabels,DEFAULT_FREQ)
    af.submenu_runFI()
    af.submenu_save()
    af.submenu_plot(0, expbase+'_FI')
    af.submenu_plot(1, expbase+'_voltage')
    '''
    
    for exparea in explist:
        af = run_analysis.AnalyticFrame()
        af.update_params({'tstop':2000,'tstart':180,'baseline':'_freq1_pw10_NpHR_none_ChR_none'})
        expbases = ['130606_decoupledapical_BAC_t2k','130606_coupledapical_BAC_t2k']
        subselect = ['_freq%s_pw10'+'_NpHR_%s_ChR_%s'%(exparea[1],exparea[0]) for i in expbases]
        labels = ['ChR %s, NpHR %s, %s'%(exparea[0],exparea[1],i) for i in ['inphase','antiphase']]
        print expbases, subselect, labels
        
        af.populate_expset(expbases, subselect,labels, DEFAULT_FREQ)
        
        af.submenu_runFI()
        af.submenu_save()
        af.submenu_plot(0, '130606_comparecouple_FI'+'_NpHR_%s_ChR_%s'%(exparea[1],exparea[0]))
        af.submenu_plot(1, '130606_comparecouple_voltage'+'_NpHR_%s_ChR_%s'%(exparea[1],exparea[0]))

 
    
def analyse2():
    expbase = '130606_decoupled_singlestim'
    af = run_analysis.AnalyticFrame()
    af.update_params({'tstop':tstop_long,'tstart':180})
    exp_comp_list = [ ['_Jex%.1f_NpHR_whole_ChR_none','ChR2 none,NpHR whole'],
                      ['_Jex%.1f_NpHR_apical_ChR_none','ChR2 none,NpHR apical'],
                      ['_Jex%.1f_NpHR_none_ChR_apical','ChR2 apical,NpHR none'],
                      ['_Jex%.1f_NpHR_none_ChR_whole','ChR2 whole,NpHR none'],
                      ['_Jex%.1f_NpHR_apical_ChR_basal','ChR2 basal,NpHR apical'],
                      ['_Jex%.1f_NpHR_basal_ChR_apical','ChR2 apical,NpHR basal'] ]
    expss = [ec[0] for ec in exp_comp_list]
    explabels = [ec[1] for ec in exp_comp_list]
    af.populate_expset(expbase,expss,explabels,epsp_amps)
    af.submenu_runFI()
    af.submenu_save()
    af.submenu_plot(0, expbase+'_FI')
    af.submenu_plot(1, expbase+'_voltage')
    
    
    
    
if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run_experiments()
    elif sys.argv[1] == 'analyse':
        analyse1()
    elif sys.argv[1] == 'analyse_ss':
        analyse2()
    else:
        print 'add run or analyse as an arg'
        
    