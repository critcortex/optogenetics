import Neuron 
import run_stimulation
from run_experiments import ExpSetup
expbase = '131104_test_newopsincode'


runon = 'local'

explist = [['whole','whole']] 
#explist = [['whole','none']] 
explist = [['none','whole']] 

light_on = 700
curr_on = 0
curr_dur = 2500
light_dur = 1000
tstop = 2500

tstop = 500

es = ExpSetup()
#areas = es.get_areas_section()
params = {}

for exparea in explist:
    for pd in range(1,2):
        print 'hereeeeeeeeeeeeeeeeeeeeeeee'
#            params['ChR_areas'] = {exparea[0]:areas[exparea[0]]}
#            params['NpHR_areas'] = {exparea[1]:areas[exparea[1]]}
        params['num_threads'] = 1
        
        params['tstart'] = 0
        params['tstop'] = tstop
        params['stim_iclamp'] = False # True
        params['iclamp'] = [{'amp':.5,'tstim':500., 'duration':500.,'location':'proximal'}]
        params['stim_epsp'] =  False
        params['epsp'] =  [{'tstim':200, 'EPSPamp':10.1, 'location':'mysoma','risetau':0.5, 'decaytau':5., 'BACdt':0.}]
        params['stim_spiketrains'] = False
        params['spiketrains'] = [{'tstims': [[200, 600],[400, 600,500]],  'locations': ['mysoma','myapic'], 'weights':[4.5,3.2],  'el': 0.02}]
        #params['spiketrains'] = [{'tstims': [[200, 600]],  'locations': [['mysoma',0.5]], 'weights':[50.],  'el': 0.02}]

        params['experiment_type'] = 'opsinonly'
        params['savedata'] = True
        
        params['mark_loc'] = {}
        params['mark_loc']['names'] = ['mysoma','myapic','proximal','distal']
        params['mark_loc']['sections'] = ['soma','apical','apical','apical']
        #params['mark_loc']['distances'] = [0.5]
        params['mark_loc']['ids'] = [(0,0.5),(0,0.5),(0,0.0753),(0,0.972326)]
        
        params['record_loc'] = {}
        params['record_loc']['v'] = ['mysoma']
        """
        params['record_loc']['v'] = ['mysoma', 'myapic','proximal','distal']
        params['record_loc']['ik'] = ['mysoma','myapic']
        params['record_loc']['ina'] = ['mysoma','myapic']
        """
        params['record_loc']['iphoto'] = ['mysoma']
        
        params['plot'] = { 'voltages': [ ['mysoma', 'v','k-'] ]} #, ['myapic','v','b-'] ] } #['proximal','v','r-'],['distal','v','b-'] ] }
        chrdict =  {'exp':5e-4, 'irradiance' :1.*pd, 'pulsewidth': 50,'lightdelay':200,'interpulse_interval':50,  'n_pulses':2}
        hrdict =  {'exp':5e-4, 'irradiance' :1.*pd, 'pulsewidth': 50,'lightdelay':250,'interpulse_interval':50,  'n_pulses':2}
        params['opsindict'] = { 'ChR' : {'soma': chrdict,
                                         'dend':chrdict,
                                         'axon': chrdict,
                                         'apic':chrdict} ,
                                'NpHR' : {'soma': hrdict,
                                         'dend' : hrdict,
                                         'axon':hrdict,
                                         'apic':hrdict}}
        
        
        params.update({#'NpHR_times':opsin_expression_np,'ChR_times':opsin_expression_ch,
                       'expname':expbase,
                       'description':'_irr%.1f_factor%.2f'%(1.,pd)})
        params['cell'] = ['Neuron','SimpleBipolar']
        params['cell'] = ['Neuron','L5PC']
        params['cell_description'] = 'test'
        pp = es.get_default_params()
        NE = run_stimulation.NeuronExperiment()
        pp.update(params)
        #es.run_single_experiment(expbase, runon, params)
        NE.main(pp)