import sys
import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import run_stimulation
import pylab
import file_io as fio

expbase = '140305_testiclamp'
tstop = 1000

es = ExpSetup()

def test_fractal():
    pp = {}
    pp['cell'] = ['Neuron', 'FractalNeuron']
    pp['cell_params'] = {'num_base': 2,
                  'num_split': 4,
                  'num_levels':3,
                  'defaultlevel':{'mechanisms':['pas_dend','hh'],#,'ih','ca_lvast','ca_hva','sk3','ske2','na_tat','ca_dyn','im'],
                                  'cm':2.,
                                  'Ra':100},
                  'soma':{'mechanisms':['pas_soma','hh'],#'ih','ca_lvast','ca_hva','sk3','ske2','ktst','kpst','na_et2','na_tat','ca_dyn'],
                          'cm':1.},
                  'dend0':{'mechanisms':['pas_dend','hh'],
                           'cm':2.,
                           'Ra':100},
                  'dend1':{'mechanisms':['pas_dend','hh'],
                           'cm':2.,
                           'Ra':100},               
                  'mechanisms': {'pas_soma':('pas',{'e':-70}),#,'g_pas':1./20}),
                                 'pas_dend':('pas',{'e':-70})},#,,'g_pas':1./10}),}
                                 'hh':('hh',),
            }
    
    clamp_neuron(pp,'fractal')
    
def test_bipolar():
    pp = {}
    pp['cell'] = ['Neuron', 'SimpleBipolar']
    pp['cell_params'] = {
                  'defaultlevel':{'mechanisms':['pas_dend','hh'],#,'ih','ca_lvast','ca_hva','sk3','ske2','na_tat','ca_dyn','im'],
                                  'cm':2.,
                                  'Ra':100},
                  'soma':{'mechanisms':['pas_soma','hh'],#'ih','ca_lvast','ca_hva','sk3','ske2','ktst','kpst','na_et2','na_tat','ca_dyn'],
                          'cm':1.},
                  'dend0':{'mechanisms':['pas_dend','hh'],
                           'cm':2.,
                           'Ra':100},
                  'dend1':{'mechanisms':['pas_dend','hh'],
                           'cm':2.,
                           'Ra':100},               
                  'mechanisms': {'pas_soma':('pas',{'e':-70}),#,'g_pas':1./20}),
                                 'pas_dend':('pas',{'e':-70})},#,,'g_pas':1./10}),}
                                 'hh':('hh',),
            }
    
    clamp_neuron(pp,'bipolar')
    
    
def clamp_neuron(pp,neurontype):
    Ias = np.arange(0.1,5.5,0.1)
    
    pp['experiment_type'] = 'opsinonly'
    pp['savedata'] = True # False #True
    
    pp['tstart'] = 0
    pp['tstop'] = tstop
    
    
    pp['mark_loc'] = {}
    pp['mark_loc']['names'] = ['mysoma','mydend0']
    pp['mark_loc']['sections'] = ['soma','dend0']
    pp['mark_loc']['ids'] = [(0,0.5),(0,0.5)] 
    
    pp['record_loc'] = {}
    pp['record_loc']['v'] = ['mysoma','mydend0']
    pp['record_loc']['ina'] = ['mysoma','mydend0']
    pp['record_loc']['ik'] = ['mysoma','mydend0']
    vplots_soma = [['mysoma','v','k']]
    pp['plot'] = {1:vplots_soma}
    pp['NpHR_areas'] = {'none'      : [None]}
    pp['ChR_areas'] = {'none'      : [None]} 
    pp['num_threads'] = 1
                                                       
    """
    pp['stim_iclamp'] = True
    for Ia in Ias:
        print Ia
        pp['iclamp'] = [{'amp':Ia,'tstim':100., 'duration':tstop,'location':'mydend0'}]
    
        pp.update({'expname':expbase,
               'description':'neurontype_%s_iclamp_%.1f'%(neurontype,Ia)})
        es.run_single_experiment(expbase, 'cluster', pp)
    """
    pp['stim_spiketrains'] = True
    for Ia in Ias:
        print Ia
        pp['spiketrains'] = [{'tstims': [np.arange(100.,700.,5)],  'locations': ['mydend0'], 'weights':[Ia],  'el': 0.1}]
    
        pp.update({'expname':expbase,
               'description':'neurontype_%s_spikestim_%.2f'%(neurontype,Ia)})
        es.run_single_experiment(expbase, 'cluster', pp)
        
def test_iclamp_train(pp,neurontype):
    
    
    pp = {}
    pp['cell'] = ['Neuron', 'FractalNeuron']
    pp['cell_params'] = {'num_base': 2,
                  'num_split': 4,
                  'num_levels':3,
                  'defaultlevel':{'mechanisms':['pas_dend','hh'],#,'ih','ca_lvast','ca_hva','sk3','ske2','na_tat','ca_dyn','im'],
                                  'cm':2.,
                                  'Ra':100},
                  'soma':{'mechanisms':['pas_soma','hh'],#'ih','ca_lvast','ca_hva','sk3','ske2','ktst','kpst','na_et2','na_tat','ca_dyn'],
                          'cm':1.},
                  'dend0':{'mechanisms':['pas_dend','hh'],
                           'cm':2.,
                           'Ra':100},
                  'dend1':{'mechanisms':['pas_dend','hh'],
                           'cm':2.,
                           'Ra':100},               
                  'mechanisms': {'pas_soma':('pas',{'e':-70}),#,'g_pas':1./20}),
                                 'pas_dend':('pas',{'e':-70})},#,,'g_pas':1./10}),}
                                 'hh':('hh',),
            }
    
    
    Ias = np.arange(0.1,5.5,0.1)
    
    pp['experiment_type'] = 'opsinonly'
    pp['savedata'] = True # False #True
    
    pp['tstart'] = 0
    pp['tstop'] = tstop
    
    
    pp['mark_loc'] = {}
    pp['mark_loc']['names'] = ['mysoma','mydend0']
    pp['mark_loc']['sections'] = ['soma','dend0']
    pp['mark_loc']['ids'] = [(0,0.5),(0,0.5)] 
    
    pp['record_loc'] = {}
    pp['record_loc']['v'] = ['mysoma','mydend0']
    pp['record_loc']['ina'] = ['mysoma','mydend0']
    pp['record_loc']['ik'] = ['mysoma','mydend0']
    vplots_soma = [['mysoma','v','k']]
    pp['plot'] = {1:vplots_soma}
    pp['NpHR_areas'] = {'none'      : [None]}
    pp['ChR_areas'] = {'none'      : [None]} 
    pp['num_threads'] = 1
                                                       
    """
    pp['stim_iclamp'] = True
    for Ia in Ias:
        print Ia
        pp['iclamp'] = [{'amp':Ia,'tstim':100., 'duration':tstop,'location':'mydend0'}]
    
        pp.update({'expname':expbase,
               'description':'neurontype_%s_iclamp_%.1f'%(neurontype,Ia)})
        es.run_single_experiment(expbase, 'cluster', pp)
    """
    pp['stim_iclamp_train'] = True
    for Ia in Ias:
        print Ia
        pp['iclamp_train'] = [{'tstims': np.arange(200.,700.,100),  'location': 'mydend0', 'weights':[Ia],'amp':Ia, 'dur':90,}]
    
        pp.update({'expname':expbase,
               'description':'neurontype_%s_spikestim_%.2f'%(neurontype,Ia)})
        #es.run_single_experiment(expbase, 'cluster', pp)        
        es.run_single_experiment(expbase, 'local', pp)
        return        
        
                            
#test_fractal()
#test_bipolar()     
test_iclamp_train({},'fractalclamp')     
     
     
     
     
     
     