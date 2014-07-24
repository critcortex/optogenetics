import sys
import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import run_stimulation
import pylab
import math
import src.file_io as fio
import pickle

class Helper():
    
    def _save_spiketime(self,ss,filename):
        # can't save txt as the lengths will be unequal
        ## np.savetxt(filename+'.txt', ss)
        output = open(filename+'.pkl','wb')
        pickle.dump(ss,output)
        output.close()
        print "Saved spiketimes to %s.pkl"%filename
        
    
    def _load_spiketimes(self,filename):
        print "Trying to load spiketimes from %s.pkl"%filename, 
        #ss = np.loadtxt(filename)
        pkl_file = open('%s.pkl'%filename, 'rb')
        ss = pickle.load(pkl_file)
        pkl_file.close()
        print "... done!"
        return ss
    
    def generate_spike_times(self,freq,trainType,nOutput=1000,tmax=None,**gen_kw):
        """
        @param freq      : frequency of output 
        @param trainType : type of process to be generated (string). 
                Permissible values: 'poisson', 'inhomog_poisson', or 'gamma'  
        
        """
        self.stg = ntst.StGen()
        if tmax is None:
            tmax = self.tstop
            
        if trainType =='poisson':
            print freq, tmax
            #print self.stg.poisson_generator(freq,t_stop=tmax,array=True)
            
            ss = [self.stg.poisson_generator(freq,t_stop=tmax,array=True) for i in range(nOutput)]
        elif trainType =='inhomog_poisson':
            inhomg_rate = self.GAMMA
            bg_r = freq
            high_r = 2*freq
            T = 2*math.pi
            w = 1000./inhomg_rate/T
            tt = np.arange(0,tmax)
            ss = [self.stg.inh_poisson_generator((high_r-bg_r)*np.sin(1./w*tt)+bg_r,tt,1000,array=True) for i in range(nOutput)]
        elif trainType =='gamma':
            ss = [self.stg.inh_gamma_generator(t_stop=tmax, array=True,**gen_kw) for i in range(nOutput)]
        else:
            print "Unknown trainType:", trainType
            return -1
        return ss
    
    def _load_input_sites(self,filename):
        print "Trying to load input sites from %s.pkl"%filename,
        #ss = np.loadtxt(filename)
        pkl_file = open('%s.pkl'%filename, 'rb')
        idx = pickle.load(pkl_file)
        pkl_file.close()
        print "... done!"
        return idx
    
########################################################## End of Helper
    

hh = Helper()
# Params ========================================================
expbase = '140718_replaylong'
tstop = 5000
light_on = 3000
light_dur = 600
celltype = 'L5PC'
# Optogene related
irrs = [0.03]
irrs = [0.002] # ,0.02]
irrs = [0.003,0.004]
irrs = np.arange(0.005,0.011,0.001)
factor = 0.1
# Current injection related
Ias = range(1,10)
# Spike input related
J = 1.
syn_factor =  -0.1
nsites_e = 200
nsites_i = 120
nsites = nsites_e+nsites_i
freq = 50
shuffle_id = 5
spikefile = 'spikes_t10k_f%g_nsites%g'%(freq,nsites) # filename for input spikes

# Used parameters ---------------
# for CNS Talk 
irrs = [0.007,0.008]
Ias = [5]
# -------------------------------

# ===============================================================

def generate_spikes():
    # Step 1: generate spike times
    ss = hh.generate_spike_times(freq,'poisson',nsites,tstop)
    hh._save_spiketime(ss,spikefile)



def setup(pp, replay_type):
    
    
    # Cell type
    pp['cell'] = ['Neuron', 'L5PC']
    pp['cell_params'] = {}
    
    # Stim spiketrains
    pp['stim_spiketrains'] = True
    
    locnames = ['idx_%g'%i for i in range(nsites_e+nsites_i)]
    #syn_factor = -1.*nsites_e/nsites_i
    weights = np.concatenate((np.ones(nsites_e)*J ,  np.ones(nsites_i)*syn_factor*J)) 
    spikes = hh._load_spiketimes(spikefile)
    spikelist = spikes[:nsites]
    pp['spiketrains'] = [{'tstims': spikelist,  'locations':  locnames, 'weights':weights,  'el': 0.02}]
    
    idx_a = hh._load_input_sites('sites_%s_apic_mindist%g_shuffle%g'%(celltype,50,shuffle_id))
    idx_d = hh._load_input_sites('sites_%s_dend_mindist%g_shuffle%g'%(celltype,50,shuffle_id))
    nsites_e_a = 296
    nsites_e_d = 88
    nsites_ex = nsites_e_a + nsites_e_d
    otherbaseslabels = ['idx_%g'%i for i in range(nsites_ex)]
    othersections = ['apic']*nsites_e_a + ['dend']*nsites_e_d
    idxzip_e = [otherbaseslabels,np.concatenate((idx_a[:nsites_e_a],idx_d[:nsites_e_d]))]
    
    idx_a = hh._load_input_sites('sites_%s_apic_mindist%g_shuffle%g'%(celltype,150,shuffle_id))
    idx_d = hh._load_input_sites('sites_%s_dend_mindist%g_shuffle%g'%(celltype,150,shuffle_id))
    nsites_i_a = 74
    nsites_i_d = 22
    nsites_ix = nsites_i_a + nsites_i_d
    otherbaseslabels_i = ['idx_%g'%i for i in range(nsites_ex,nsites_ex+nsites_ix)]
    othersections_i = ['apic']*nsites_i_a + ['dend']*nsites_i_d
    idxzip_i = [otherbaseslabels_i,np.concatenate((idx_a[:nsites_i_a],idx_d[:nsites_i_d]))]
        
    #print olocs
    pp['mark_loc'] = {}
    pp['mark_loc']['names'] = ['mysoma']+['tmpA','tmpB']#+otherbaseslabels+otherbaseslabels_i
    pp['mark_loc']['sections'] = ['soma']+['tmpA','tmpB']#+othersections+othersections_i
    #pp['mark_loc']['ids'] = [(0,0.5)] + [('get_section_byname',{'sectioname':loc,'subdomain':'dend1'}) for loc in otherbaseslabels]+[(0,loc[1]) for loc in olocs ]
    pp['mark_loc']['ids'] = [(0,0.5)] + [{'idx':idxzip_e, 'sections':othersections}] + [{'idx':idxzip_i, 'sections':othersections_i}] 
    
    
    pp['record_loc'] = {}
    pp['record_loc']['v'] = ['mysoma']#+otherbaseslabels

    vplots_soma = [['mysoma','v','k']]
    pp['plot'] = {1:vplots_soma}
    
    pp['num_threads'] = 1
    
    pp['experiment_type'] = 'opsinonly'
    pp['savedata'] = True # False #True
        
    pp['tstart'] = 0
    pp['tstop'] = tstop
    
    pp.update({'expname':expbase,
               'description':'_replay_%s_freq%g'%(replay_type,freq)})
    
    es = ExpSetup()
    #es.run_single_experiment(expbase, 'cluster', pp)
    es.run_single_experiment(expbase, 'local', pp)

def simulate_normal():
    """
    We don't need to do anything apart from let the simulation run 'as is'
    """
    pp = {}
    pp['opsindict'] = {}
    pp['opsindict']['ChR'] = {}
    pp['opsindict']['NpHR'] = {}
    pp['ChR_areas'] = {'none'      : [None]}
    pp['NpHR_areas'] = {'none'      : [None]}
    setup(pp,'normal')

def simulation_injectedcurrent():
    pp = {}
    pp['opsindict'] = {}
    pp['opsindict']['ChR'] = {}
    pp['opsindict']['NpHR'] = {}
    pp['ChR_areas'] = {'none'      : [None]}
    pp['NpHR_areas'] = {'none'      : [None]}
    pp['stim_iclamp'] = True
    stimloc = 'mysoma'
    for Ia in Ias:
        pp['iclamp'] = [{'tstim': light_on,  'location': stimloc, 'amp':Ia, 'duration':light_dur}]
        setup(pp,'current_Ia%g'%Ia)

def simulation_injectedcurrent_biphase():
    # WARNING: not verified to work correctly atm 22/07/14
    pp = {}
    pp['opsindict'] = {}
    pp['opsindict']['ChR'] = {}
    pp['opsindict']['NpHR'] = {}
    pp['ChR_areas'] = {'none'      : [None]}
    pp['NpHR_areas'] = {'none'      : [None]}
    pp['stim_iclamp'] = True
    stimloc = 'mysoma'
    lightdur2 = light_dur/2
        
    for Ia in Ias:
        pp['iclamp'] = [{'tstim': light_on,  'location': stimloc, 'amp':Ia, 'duration':lightdur2}]#,
                       # {'tstim': light_on+lightdur2,  'location': stimloc, 'amp':-Ia, 'duration':lightdur2}]
        setup(pp,'current_biphasictt_Ia%g'%Ia)



def simulate_gainmod():
    """
    Add in the optogenetics section and let run
    
    """
    pp = {}
    pp['opsindict'] = {}
    pp['opsindict']['ChR'] = {}
    pp['opsindict']['NpHR'] = {}
        
    areas = ['soma', 'apic','dend']
    
    for irr in irrs:
        chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
        hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
        for area in areas:
            pp['opsindict']['ChR'][area] =  chrdict
            pp['opsindict']['NpHR'][area] = hrdict
        pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
        pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
        setup(pp,'optogen_irr%.3f_factor%.2f'%(irr,factor))






#generate_spikes()
#simulate_normal()
#simulation_injectedcurrent()
#simulation_injectedcurrent_biphase()
#simulate_gainmod()



