import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import run_stimulation
import pylab
import file_io as fio

expbase = '140212_testActive'


pp = {}
# neuron model and params
pp['cell'] = ['Neuron', 'FractalNeuron']
pp['cell_params'] = {'num_base': 2,
      'num_split': 3,
      'num_levels':2,
      'defaultlevel':{'mechanisms':['pas_dend'],#,'ih','ca_lvast','ca_hva','sk3','ske2','na_tat','ca_dyn','im'],
                      'cm':2.,
                      'Ra':100},
      'dend0':{'mechanisms':['pas_dend'],#'ih','ca_lvast','ca_hva','sk3','ske2','na_tat','ca_dyn','im'],
               'cm':2.},
      'soma':{'mechanisms':['pas_soma','hh'],#'ih','ca_lvast','ca_hva','sk3','ske2','ktst','kpst','na_et2','na_tat','ca_dyn'],
              'cm':1.},
      'mechanisms': {'pas_soma':('pas',{'e':-70}),#,'g_pas':1./20}),
                     'pas_dend':('pas',{'e':-70}),#,,'g_pas':1./10}),
                      'ih':('Ih',{}),
                     'ca_lvast':('Ca_LVAst',{}),
                     'ca_hva':('Ca_HVA',{}),
                     'sk3':('SKv3_1',{}),
                     'ske2':('SK_E2',{}),
                     'ktst':('K_Tst',{}),
                     'kpst':('K_Pst',{}),
                     'na_et2':('Nap_Et2',{}),
                     'na_tat':('NaTa_t',{}),
                     'ca_dyn':('CaDynamics_E2',{}),
                     'im': ('Im',{})}
                     }


pp['opsindict'] = {}
irr = 50.
factor = 0.1
light_on = 250
light_dur = 50 #0
chrdict =  {'exp':100e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
hrdict =  {'exp':100e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}

"""
pp['opsindict']['NpHR'] =  {'soma': hrdict,
                            'dend0': hrdict}
                                """
pp['NpHR_areas'] = {'none'      : [None]}
pp['ChR_areas'] = {'none'      : [None]}
pp['experiment_type'] = 'opsinonly'
pp['savedata'] = True 
                    
pp['tstart'] = 0
pp['tstop'] = 1000
pp['mark_loc'] = {}
pp['mark_loc']['names'] = ['mysoma']
pp['mark_loc']['sections'] = ['soma']
pp['mark_loc']['ids'] = [(0,0.5)]
pp['record_loc'] = {}
pp['record_loc']['v'] = ['mysoma']
pp['record_loc']['ina'] = ['mysoma']
#pp['record_loc']['ica'] = ['mysoma']
pp['record_loc']['iphoto'] = ['mysoma']
#pp['record_loc']['iChR'] = ['mysoma']
pp['record_loc']['ik'] = ['mysoma']
pp['plot'] = {1:[['mysoma','v','k-']]}
            #,
            #  'curr': [['mysoma','ina','k-'],['mysoma','ica','b-'],['mysoma','ik','r-'],['mysoma','i_ChR','g-']]}                    
def full():
    es = ExpSetup()
    irrs = np.arange(5,50.1,5)
    expressions = np.arange(100e-4,500e-4,50e-4)
    for irr in irrs:
        for express in expressions:
            pp['expname'] = expbase
            pp['description'] = '_irr%g_exp%g'%(irr,express)
            factor = 0.1
            light_on = 250
            light_dur = 500
            chrdict =  {'exp':express, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
            hrdict =  {'exp':express, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
            pp['opsindict']['ChR'] =  {'soma': chrdict,
                                       'dend0': chrdict}
    
            es.run_single_experiment(expbase, 'cluster', pp)

def patial():
    es = ExpSetup()
    irrs = np.arange(5,50.1,5)
    expressions = np.arange(1000e-4,10010e-4,1000e-4)
    expressions = np.arange(100e-4,1001e-4,100e-4)
    for irr in irrs:
        for express in expressions:
            pp['expname'] = expbase
            pp['description'] = 'soma_irr%g_exp%g_pulsed'%(irr,express)
            factor = 0.1
            light_on = 250
            light_dur = 500
            n_pulses = 1#00
            chrdict =  {'exp':express, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':light_dur,  'n_pulses':n_pulses}
            hrdict =  {'exp':express, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':light_dur,  'n_pulses':n_pulses}
            pp['opsindict']['ChR'] =  {'soma': chrdict}#,
                                      # 'dend0': chrdict}
    
            es.run_single_experiment(expbase, 'cluster', pp)


def patial_slowpulsed():
    es = ExpSetup()
    irrs = np.arange(5,50.1,5)
    irrs = np.arange(1,5.,1)
    expressions = np.arange(1000e-4,10010e-4,1000e-4)
    expressions = np.arange(100e-4,1001e-4,100e-4)
    for irr in irrs:
        for express in expressions:
            pp['expname'] = expbase
            pp['description'] = 'soma_irr%g_exp%g_slowpulsed'%(irr,express)
            factor = 0.1
            light_on = 250
            light_dur = 5
            n_pulses = 40
            interpulse = 20
            chrdict =  {'exp':express, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':interpulse,  'n_pulses':n_pulses}
            hrdict =  {'exp':express, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':interpulse,  'n_pulses':n_pulses}
            pp['opsindict']['ChR'] =  {'soma': chrdict}#,
                                      # 'dend0': chrdict}
    
            es.run_single_experiment(expbase, 'cluster', pp)
            
def large_full_tree():
    es = ExpSetup()
    nb = 2
    ns = 3
    nl = 2
    pp['cell_params'] = {'num_base': 2,
      'num_split': 3,
      'num_levels':3}
    irrs = np.arange(1,10.1,0.5)
    expressions = np.arange(1000e-4,10010e-4,1000e-4)
    expressions = np.arange(100e-4,1001e-4,100e-4)
    expressions = np.arange(5e-4,101e-4,5e-4)
    for irr in irrs:
        for express in expressions:
            pp['expname'] = expbase
            pp['description'] = 'largetree_sub_irr%g_exp%g_nb%g_ns%g_nl%g'%(irr,express,nb,ns,nl)
            factor = 0.1
            light_on = 250
            light_dur = 500
            n_pulses = 1#00
            chrdict =  {'exp':express, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':light_dur,  'n_pulses':n_pulses}
            hrdict =  {'exp':express, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':light_dur,  'n_pulses':n_pulses}
            pp['opsindict']['ChR'] =  {'soma': chrdict,
                                       #'dend1': chrdict,
                                       'dend0': chrdict}
    
            #es.run_single_experiment(expbase, 'cluster', pp)
            NE = run_stimulation.NeuronExperiment()
            dp = es.get_default_params()
            dp.update(pp)
            NE.main(dp)
            return

large_full_tree()