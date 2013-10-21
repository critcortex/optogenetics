import Neuron 
import run_stimulation_old
import run_experiments
import sys


expbase = '131021_testSimpleNeuron'


runon = 'local'

explist = [['whole','whole']] 
explist = [['soma','soma']] 


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
        params['num_threads'] = 1
        
        params['tstart'] = 0
        params['tstop'] = tstop
        params['stim_iclamp'] = False # True
        params['iclamp'] = [{'amp':.5,'tstim':500., 'duration':500.,'location':'proximal'}]
        params['stim_epsp'] =  False
        params['epsp'] =  [{'tstim':200, 'EPSPamp':10.1, 'location':'mysoma','risetau':0.5, 'decaytau':5., 'BACdt':0.}]
        params['stim_spiketrains'] = True
        params['spiketrains'] = [{'tstims': [[1200, 1600],[1400, 1600,1900]],  'locations': ['mysoma','myapic'], 'weights':[0.5,0.2],  'el': 0.02}]
        #params['spiketrains'] = [{'tstims': [[200, 600]],  'locations': [['mysoma',0.5]], 'weights':[50.],  'el': 0.02}]

        params['experiment_type'] = 'opsinonly'
        params['savedata'] = True
        
        params['mark_loc'] = {}
        params['mark_loc']['names'] = ['mysoma','myapic','proximal','distal']
        params['mark_loc']['sections'] = ['soma','apical','apical','apical']
        #params['mark_loc']['distances'] = [0.5]
        params['mark_loc']['ids'] = [(0,0.5),(0,0.5),(0,0.0753),(0,0.972326)]
        
        params['record_loc'] = {}
        params['record_loc']['v'] = ['mysoma', 'proximal','distal']
        params['record_loc']['i_k'] = ['mysoma','myapic']
        params['record_loc']['i_na'] = ['mysoma','myapic']
        
        params['plot'] = { 'voltages': [ ['mysoma', 'v','k-'],['proximal','v','r-'],['distal','v','b-'] ] }
        params['opdict'] = {}
        
        irr = 1.
        factor = 1.
        opsin_expression_ch = [5e-4,irr,light_on,light_dur,light_dur,600,1] 
        opsin_expression_np = [5e-4,irr*factor,light_on,light_dur,light_dur,600,1] 
        params.update({'NpHR_times':opsin_expression_np,'ChR_times':opsin_expression_ch,'description':'_irr%.1f_factor%.2f'%(irr,factor)})
        params['cell'] = ['Neuron','SimpleBipolar']
        params['cell'] = ['Neuron','L5PC']
        params['cell_description'] = 'SimpleBipolar'
        pp = es.get_default_params()
        pp.update(params)
        es.run_single_experiment(expbase, runon, params)
        #run_stimulation_old.main(pp)



        
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        pass
    elif sys.argv[1] == 'run':
        run_experiments()        