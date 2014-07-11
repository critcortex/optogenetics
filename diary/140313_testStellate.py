import sys
import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import run_stimulation
import pylab
import file_io as fio

expbase = '140313_tagStellate'
tstop = 1000
num_proc = 4
es = ExpSetup()
stg = ntst.StGen()


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
    pp['cell'] = ['Neuron', 'Stellate']
    pp['cell_params'] = {}
    
    
    Ias = np.arange(50.,100.,1)
    
    pp['experiment_type'] = 'opsinonly'
    pp['savedata'] = True # False #True
    
    pp['tstart'] = 0
    pp['tstop'] = tstop
    
    
    pp['mark_loc'] = {}
    pp['mark_loc']['names'] = ['mysoma','mydend0']
    pp['mark_loc']['sections'] = ['soma','soma']
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
    count = 0
    pp['stim_iclamp_train'] = True
    for Ia in Ias:
        print Ia
        pp['iclamp_train'] = [{'tstims': np.arange(200.,700.,100),  'location': 'mydend0', 'weights':[Ia],'amp':Ia, 'dur':90,}]
    
        pp.update({'expname':expbase,
               'description':'neurontype_%s_spikestim_%.2f'%(neurontype,Ia)})
        es.run_single_experiment(expbase, 'cluster', pp)        
        #es.run_single_experiment(expbase, 'local', pp)
        count += 1
        return
        if count >= num_proc:
            return

def get_spiketimes(rate,num_inputs):
    spikes = []
    for i in range(num_inputs):
        ss = stg.poisson_generator(rate,t_stop=tstop,array=True)
        spikes.append(ss)
    #print 'Spikes = ',spikes
    return spikes

def test_stellate_train(pp,neurontype):
    
    count = 0
    Js = [30000]
    freq = 80
    pp = {}
    pp['cell'] = ['Neuron', 'Stellate']
    pp['cell_params'] = {}
    
    
    
    pp['experiment_type'] = 'opsinonly'
    pp['savedata'] = True # False #True
    
    pp['tstart'] = 0
    pp['tstop'] = tstop
    
    
    pp['mark_loc'] = {}
    pp['mark_loc']['names'] = ['mysoma','mydend0']
    pp['mark_loc']['sections'] = ['soma','soma']
    
    nsites = 1
    dist = [50,100]
    
    labels = ['stim%g'%i for i in range(nsites)]
    pp['mark_loc'] = {}
    pp['mark_loc']['names'] = ['mysoma']+labels
    pp['mark_loc']['sections'] = ['soma']+['a1']
    pp['mark_loc']['ids'] = [(0,0.5)] + [('select_section_posn_bydistance',{'sectionarea':'a1'}) for i in range(nsites)] 
    
    pp['record_loc'] = {}
    pp['record_loc']['v'] = ['mysoma']

    vplots_soma = [['mysoma','v','k']]
    pp['plot'] = {1:vplots_soma}
    pp['NpHR_areas'] = {'none'      : [None]}
    pp['ChR_areas'] = {'none'      : [None]} 
    pp['num_threads'] = 1
    
    """
    pp['stim_spiketrains'] = True 
    for J in Js: 
        pp['spiketrains'] = [{'tstims': get_spiketimes(freq,nsites),  'locations': labels, 'weights':np.ones(nsites)*J,  'el': 0.02}]
        pp.update({'expname':expbase,
               'description':'neurontype_%s_spiketrain%g_freq%g'%(neurontype,J,freq)})
        #es.run_single_experiment(expbase, 'cluster', pp)        
        es.run_single_experiment(expbase, 'local', pp)
        return
        count += 1
    print count,'jobs submitted'
    """
        
    pp['stim_iclamp'] = True 
    for J in Js: 
        pp['iclamp'] = [{'amp':J,'tstim':100., 'duration':100,'location':labels[0]}]
        #pp['iclamp'] = [{'amp':J,'tstim':100., 'duration':100,'location':'mysoma'}]
        pp.update({'expname':expbase,
               'description':'neurontype_%s_spiketrain%g_freq%g'%(neurontype,J,freq)})
        #es.run_single_experiment(expbase, 'cluster', pp)        
        es.run_single_experiment(expbase, 'local', pp)
        return
        count += 1
    print count,'jobs submitted'
    

def create_run_stellate():
    from neuron import h
    
    cell = Neuron.Stellate()
    clamp_site = cell.select_section_posn_bydistance('a2')[0]
    iclamp = h.IClamp()
    iclamp.loc(clamp_site[2], sec=clamp_site[1])
    
    return cell

    
        
#test_fractal()
#test_bipolar()     
test_stellate_train({},'test')     
     
     
     
     
     
     