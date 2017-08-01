'''
Created on 19 Oct 2015

@author: sjarvis

Run L5PC experiment but with template of neuron that has Ca/K/Na/Ih 
ion channels removed and replaced with plain HH channels

'''

#import sys
#import Neuron
from run_experiments import ExpSetup
#import NeuroTools.stgen as ntst 
import numpy as np
#import run_stimulation
#import pylab
#import file_io as fio
#import run_analysis
import NeuroTools.stgen as ntst 
#import math

stg = ntst.StGen()


expbase = '151019_removeL5PCChannels'
es = ExpSetup()
cellsections = ['apic','dend']
irrs = [1] # np.arange(0,5.1,0.5)
factors = [0.25] #[0.25,0.5,1.,2.]
Ias = [30.] #np.arange(50,1000.1,50.)

light_on = 150
light_dur = 200
tstop = 100
nsites = 1

dist_vitro = [[101,102]]
inj_locations = "apic"

nsite_range = [300]
dist_vivo = [20,-1] # [100,-1]
freqs = [705] #range(150,301,50)
J = 1550.


def get_spiketimes(rate,num_inputs):
    spikes = []
    for i in range(num_inputs):
        # TODO: am guessing there's a more efficient way of implementing this
        ss = stg.poisson_generator(rate,t_stop=tstop,array=True)
        spikes.append(ss)
    #print 'Spikes = ',spikes
    return spikes

def run_L5PC_vitro():
    
    for irr in irrs:
        for factor in factors:
            for Ia in Ias:
                
                pp = setup_L5PC_vitro()
                # set opsin
                chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                
                pp['opsindict'] = {}
                if irr > 0 :
                    pp['opsindict']['ChR'] = {'soma': chrdict}
                    for subdomain in cellsections:
                        pp['opsindict']['ChR'][subdomain] = chrdict
                    pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
                else:
                    pp['ChR_areas'] = {'none'      : [None]}
                    
                if irr > 0 and factor > 0:
                    pp['opsindict']['NpHR'] =  {'soma': hrdict}
                    for subdomain in cellsections:
                        pp['opsindict']['NpHR'][subdomain] = hrdict
                    pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
                else:
                    pp['NpHR_areas'] = {'none'      : [None]}
    
    
                # Setup current injection
                pp['stim_iclamp'] = False #True
                #print labels[0]
                #pp['spiketrains'] = [{'tstims': get_spiketimes(freq,nsites), 'locations': labels, 'weights':np.ones(nsites)*J,  'el': 0.1}]
                stimloc = 'stimloc' # for dendritic location
                stimloc = 'mysoma' # for somatic location
                pp['iclamp'] = [{'tstim':50,  'location': stimloc, 'amp':Ia, 'duration':tstop-100}]
                
                # Set name and submit
                pp.update({'expname':expbase,
                           'description':'irr%.3f_factor%.2f_I%.2f_stimloc_%s'%(irr,factor,Ia,'stimloc')})

                es.run_single_experiment(expbase, 'local', pp)
                #return
    
def run_L5PC_vivo():
    
    for irr in irrs:
        for factor in factors:
            for freq in freqs:
                for nsites_vivo in nsite_range:
                
                    pp = setup_L5PC_vitro()
                    # set opsin
                    chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                    hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                    
                    pp['opsindict'] = {}
                    if irr > 0 :
                        pp['opsindict']['ChR'] = {'soma': chrdict}
                        for subdomain in cellsections:
                            pp['opsindict']['ChR'][subdomain] = chrdict
                        pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
                    else:
                        pp['ChR_areas'] = {'none'      : [None]}
                        
                    if irr > 0 and factor > 0:
                        pp['opsindict']['NpHR'] =  {'soma': hrdict}
                        for subdomain in cellsections:
                            pp['opsindict']['NpHR'][subdomain] = hrdict
                        pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
                    else:
                        pp['NpHR_areas'] = {'none'      : [None]}
        
                    
                    labels = ['stim%g'%i for i in range(nsites_vivo*2)]
                            
                    
                    pp['mark_loc'] = {}
                    pp['mark_loc']['names'] = ['mysoma']+labels
                    pp['mark_loc']['sections'] = ['soma']+['apic']*nsites_vivo+['dend']*nsites_vivo
                    pp['mark_loc']['ids'] = [(0,0.5)] \
                        + [('select_section_posn_bydistance',{'sectionarea':'apic','mindist':dist_vivo[0],'maxdist':dist_vivo[1]}) for i in range(nsites_vivo)] \
                        + [('select_section_posn_bydistance',{'sectionarea':'dend','mindist':dist_vivo[0],'maxdist':dist_vivo[1]}) for i in range(nsites_vivo)] 
                    
                            
                    
                    
                    
                    pp['stim_spiketrains'] = True
                                    #locations = ['mysoma','myapic']
                    
                    
                    
                    pp['spiketrains'] = [{'tstims': get_spiketimes(freq,nsites_vivo),  'locations': labels, 'weights':np.ones(nsites_vivo)*J,  'el': 0.02}]
                    
                    
                    
                    
                    
                    # Set name and submit
                    pp.update({'expname':expbase,
                               'description':'irr%.3f_factor%.2f_freq%g_nsites%g_J%g'%(irr,factor,freq,nsites_vivo,J)})
        
                    es.run_single_experiment(expbase, 'local', pp)
                    #return

    
    

def setup_L5PC_vitro():
  
    pp = {}
    # neuron model and params
    pp['cell'] = ['Neuron', "L5PC"]
    pp['cell_params'] = {'model':'L5PC_biophys_plainHH.hoc'}
    
    
    # general settings 
    pp['experiment_type'] = 'opsinonly'
    pp['savedata'] = True # False #True
    
    pp['tstart'] = 0
    pp['tstop'] = tstop
    
    
    labels = ['stim%g'%i for i in range(nsites)]
    
    
    pp['mark_loc'] = {}
    pp['mark_loc']['names'] = ['mysoma']+['stimloc']
    pp['mark_loc']['sections'] = ['soma']+[inj_locations]
    pp['mark_loc']['ids'] = [(0,0.5)] + [('select_section_posn_bydistance',{'sectionarea':inj_locations,'mindist':dists[0],'maxdist':dists[1]}) for (i,dists) in enumerate(dist_vitro)] 
    #TODO: update
    
    pp['record_loc'] = {}
    pp['record_loc']['v'] = ['mysoma'] #+labels
    
    vplots_soma = [['mysoma','v','b']]

    #iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
    pp['plot'] = {1:vplots_soma} 
        
    pp['num_threads'] = 1
    return pp
    
run_L5PC_vivo()    