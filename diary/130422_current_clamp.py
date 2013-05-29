import sys
import numpy as np
import run_experiments
import run_analysis


ie_range = [0.0] + [coeff*ten for coeff in [1,2,5] for ten in [1e-2,1e-1,1]]
max_amp = 2.0
ie_range = np.concatenate((np.arange(-1,0.,0.1) , np.arange(0,max_amp+0.01,0.1)))
#ie_range  = np.arange(0.4,0.61,0.02)
ie_range = np.concatenate((np.arange(-1,0.4,0.1),np.arange(0.4,0.61,0.02),np.arange(0.7,1.51,0.1)))

expbase = '130422_currentclamp'
tstart = 400
tstop  = 800
tsim = 1000.

num_dec_place = 2 ################# 1

def run_experiment():

    es = run_experiments.ExpSetup()
    areas = es.get_areas_section()
    
    params = {'NpHR_areas': {'none':areas['none']}, #'whole':areas['whole']}, #,, 
              'ChR_areas' : {'none':areas['none']},
              'experiment_type': 'BAC',
              'iclamp_start'      : tstart,
              'iclamp_duration'   : tstop-tstart,
              'NpHR_times': [5e-4,10,0,tsim,350,600,1],
              'tstop'     : tsim}
    
    for ie in ie_range:
        params['iclamp_amp'] = ie
        params['description'] = '_ie%.2f'%(ie)
        es.run_single_experiment(expbase,runon='cluster',params=params)



def analyse():

    af = run_analysis.AnalyticFrame()
    af.update_params({'tstop':tstop,'tstart':tstart})
    
    
    exp_comp_list = [ ['_ie%s_NpHR_whole_ChR_none','ChR none,NpHR whole'],
                      ['_ie%s_NpHR_none_ChR_none','ChR none,NpHR none']]
                 
    expss = [ec[0] for ec in exp_comp_list]
    explabels = [ec[1] for ec in exp_comp_list]
    af.populate_expset(expbase,expss,explabels,['%.2f'%d for d in ie_range])
    af.run_analysis_menu()
    
    
    
    
if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run_experiment()
    elif sys.argv[1] == 'analyse':
        analyse()
    else:
        print 'add run or analyse as an arg'