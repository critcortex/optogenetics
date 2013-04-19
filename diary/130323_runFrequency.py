import run_experiments
DEFAULT_FREQ = [0.5,1,1.5,2,2.5,4,5,6,7,8,9,10,15,20,30,35,40,45,50]


es = run_experiments.ExpSetup()
es.run_frequency_range('130418_FIcurve',freqs=DEFAULT_FREQ,pulsewidth=10.,transient=200.,n_pulses=100)
es.run_frequency_range('130418_FIcurveBAC',freqs=DEFAULT_FREQ,pulsewidth=10.,transient=200.,n_pulses=100,params={'experiment_type':'BAC'})


#es.run_test_exp('130417_test',params={'num_threads':2})