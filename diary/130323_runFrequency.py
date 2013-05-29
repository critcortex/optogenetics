import run_experiments
DEFAULT_FREQ = [0.5,1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,20,30,35,40,45,50]

runon = 'cluster'

es = run_experiments.ExpSetup()
areas = es.get_areas_section()
params = {}
params['tstop'] = 2000
params['NpHR_areas'] = {'none':areas['none'],'soma':areas['soma'],'apical':areas['apical'],'whole':areas['whole']} # 
params['ChR_areas'] = {'whole':areas['whole']}
es.run_frequency_range('130419_FIcurve_t2k',freqs=DEFAULT_FREQ,pulsewidth=10.,transient=200.,n_pulses=100,params=params,runon=runon)
params['experiment_type'] = 'BAC'
es.run_frequency_range('130419_FIcurveBAC_t2k',freqs=DEFAULT_FREQ,pulsewidth=10.,transient=200.,n_pulses=100,params=params,runon=runon)
