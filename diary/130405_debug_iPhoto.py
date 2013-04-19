import run_stimulation as run_stim
import run_experiments as run_exp

opsin_expression= [5e-4,10,200,50,400,600,10]
opsin2_expression= [5e-4,10,400,50,400,600,10]
exp_type = 'BAC'
tstart = 0
tstop = 2000


params = run_exp.get_default_params()
params['experiment_type'] = exp_type
params['expname'] = 'debug'
params['ChR_areas'] = {'soma' : ['soma']}
params['opdict'] = {'ChR': [['soma', opsin_expression]], 'NpHR': [['soma', opsin2_expression]]}
print params


run_stim.setup(tstart,tstop)
run_stim.set_experiment_type(exp_type)
run_stim.set_stimulus()
run_stim.add_optogenetics(params['opdict'])
run_stim.setup_record(params['opdict'])
run_stim.setup_plot()
run_stim.simulate_exp()
run_stim.save_data(params['expname'],True,params['opdict'])
run_stim.run_plots(params)