import sys
import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import run_stimulation
import file_io as fio
import run_analysis
import math
import pickle
import random

conv_inho_rate = { 20: {'bg_r':10,'high_r':50},
                   50: {'bg_r':50,'high_r':100},
                   80: {'bg_r':80,'high_r':160}}


cell_stim_params = {'SHSteallate':  {'nsites':125,
                                     'sectionlist': ['dendrite']},
                    'L5PC':         {'nsites':125,
                                     'sectionlist': ['apic','dend']}}

class Experiment_120614:
    GAMMA = 40
    
    def __init__(self,whole,stimulus,celltype,spiketrain_type):
        self.es = ExpSetup()
        
        self.stg = ntst.StGen()
        self.whole = (whole == 'whole')
        print 'whole = ', self.whole
        if not self.whole:
            self.explist = ['partialSame']
        else:
            self.explist = ['whole']
            
        self.stimulus = stimulus
        self.celltype = celltype
        self.spiketrain_type = spiketrain_type
        
        self.expbase = '140612_replay_input'
    
    
        
        self.stim_location_fracts = [['peripheral',[0.5,1.0]],
                 ['distributed',[0,1.0]],
                 ['proximal',[0,.5]]]
    
        self.tstop = 4000
    
        self.freq = 50.
        self.dist_frac_sections = 50.

           
        self.light_on = 2000
        self.light_dur = 1000
        
        #self.factors = [0,0.125,0.25,1.,4.,8.]
        self.factors = [0.,0.125,0.25,0.325,0.5,0.75,1.,1.5,2.,3.,4.,6.,8.]
        
        #self.clamp_amps = np.arange(0.2,5.,0.2)
        #self.clamp_amps = [1.5]
        
        
        
        #self.factors = [0.5,2,1.5,1.,0.75,1.25]
        #self.factors = [0,0.5,2,1.5,1.,0.75,1.25]
        #self.factors = [10.]
        #self.irrs = [0.05]#,0.02]
        #self.irrs = [0.0005,0.005]#,0.02]
        self.irrs = [0.015]#,0.02]
        #self.irrs = np.arange(0.005,0.0101,0.001)
        #self.freqs = range(50,81,10)#101,10)
        self.freqs = [20,50,80]
        self.freqs = [20,50]
        self.freqs = [50]
        self.freq_factor = 2
        
        
        #self.Js = np.arange(0.1,1.1,0.1) # for soma simulation
        self.Js = [0.5,0.7]
        
        #self.syn_factors = [-190./50, -2.]
        self.syn_factors = [-1.2]
        
        self.nsites = 1000
        
        #self.shuffles = range(10)
        self.shuffles = range(5)
    
    
    def run_experiments(self):
        count = 0
        for freq in self.freqs:
            spikes = self._load_spiketimes('spikes_poisson_freq%g'%freq)
            spikes_i = self._load_spiketimes('spikes_poisson_freq%g'%(freq*self.freq_factor))
            for factor in self.factors:
                for irr in self.irrs: 
                    for (i,J) in enumerate(self.Js):
                        for syn_factor in self.syn_factors:
                            for shuf in self.shuffles:
                                if freq == 0 and i > 0:
                                    continue
                                
                                pp = {}
                                pp['stim_spiketrains'] = True
                                nsites_e = 190
                                nsites_i = 50
                                locnames = ['idx_%g'%i for i in range(nsites_e+nsites_i)]
                                #syn_factor = -1.*nsites_e/nsites_i
                                weights = np.concatenate((np.ones(nsites_e)*J ,  np.ones(nsites_i)*syn_factor*J)) 
                                print len(spikes[:nsites_e]), spikes[0]
                                spikelist = spikes[:nsites_e] + spikes_i[:nsites_i]
                                
                                pp['spiketrains'] = [{'tstims': spikelist,  'locations':  locnames, 'weights':weights,  'el': 0.02}]
                                pp.update({'expname':self.expbase,
                                           'description':'irr%.3f_factor%.2f_spikes%g_J%.1f_cell%s_inh%.2f_freqfactor%g_shuffle%g'%(irr,factor,freq,J,self.celltype,syn_factor,self.freq_factor,shuf)})
                                    
                                self._run_single_exp(pp, irr, factor, self.whole,shuffle_id=shuf)
                                count += 1
                                return
                        
        print '%g jobs submitted'%count
        
    def _run_single_exp(self,pp,irr,factor,whole,other_tag_locations=None,runon=True,shuffle_id=0):
        
        # neuron model and params
        pp['cell'] = ['Neuron', self.celltype]
        pp['cell_params'] = {}
        # opsin
        
        #print locs
        chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
        hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
        
        
        pp['opsindict'] = {}
        areas = ['soma', 'dendrite']
        pp['opsindict']['ChR'] = {}
        pp['opsindict']['NpHR'] = {}
        if irr > 0 and self.whole:
            print 'Case 1'
            for area in areas:
                pp['opsindict']['ChR'][area] =  chrdict
                pp['opsindict']['NpHR'][area] = hrdict
            
            pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
            pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
        elif irr > 0 and not self.whole:
            
            print 'Case 2'
            areas = ['dendrite']
            for area in areas:
                pp['opsindict']['ChR'][area] =  chrdict
                pp['opsindict']['NpHR'][area] = hrdict
            
            pp['NpHR_areas'] = {'partialSame'      : pp['opsindict']['NpHR'].keys()}
            pp['ChR_areas'] = {'partialSame'      : pp['opsindict']['ChR'].keys()}
        elif irr==0 and self.whole:
            print 'Case 3'
            pp['NpHR_areas'] = {'none'      : [None]}
            pp['ChR_areas'] = {'none'      : [None]}
        else: # irr = 0 and partial only --> so why bother simulating
            print 'Case 4'
            return
        """
        if factor ==0:
            pp['NpHR_areas'] = {'none'      : [None]}
            pp['opsindict']['NpHR'] = {}
        """
        
        # general settings 
        pp['experiment_type'] = 'opsinonly'
        pp['savedata'] = True # False #True
        
        pp['tstart'] = 0
        pp['tstop'] = self.tstop
        
        
        
        idx = self._load_input_sites('sites_%s_mindist%g_shuffle%g'%(self.celltype,50,shuffle_id))
        nsites_e = 190
        otherbaseslabels = ['idx_%g'%i for i in range(nsites_e)]
        othersections = ['dendrite']*nsites_e
        idxzip_e = [otherbaseslabels,idx]
        
        idx = self._load_input_sites('sites_%s_mindist%g_shuffle%g'%(self.celltype,150,shuffle_id))
        nsites_i = 50
        otherbaseslabels_i = ['idx_%g'%i for i in range(nsites_e,nsites_e+nsites_i)]
        othersections_i = ['dendrite']*nsites_i
        idxzip_i = [otherbaseslabels_i,idx]
        
        
        #print olocs
        pp['mark_loc'] = {}
        pp['mark_loc']['names'] = ['mysoma']+['tmpA','tmpB']#+otherbaseslabels+otherbaseslabels_i
        pp['mark_loc']['sections'] = ['soma']+['tmpA','tmpB']#+othersections+othersections_i
        #pp['mark_loc']['ids'] = [(0,0.5)] + [('get_section_byname',{'sectioname':loc,'subdomain':'dend1'}) for loc in otherbaseslabels]+[(0,loc[1]) for loc in olocs ]
        pp['mark_loc']['ids'] = [(0,0.5)] + [{'idx':idxzip_e, 'sections':othersections}] + [{'idx':idxzip_i, 'sections':othersections_i}] 
        
        
        pp['record_loc'] = {}
        pp['record_loc']['v'] = ['mysoma']#+otherbaseslabels
        pp['record_loc']['ina'] = ['mysoma']#+otherbaseslabels
        pp['record_loc']['ik'] = ['mysoma']#+otherbaseslabels
        
        vplots_soma = [['mysoma','v','k']]
        """
        iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
        vplots = [['mysoma','v','k']]
        iplots_k = [['mysoma','ik','b']]
        iplots_na = [['mysoma','ina','g']]
        
        iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
        """
        pp['plot'] = {#1:vplots,
                      2:vplots_soma}#,
        """
                      #3:iplots_soma,
                      #4:iplots_k,
                      #5:iplots_na}
        """
        pp['num_threads'] = 1
        if runon:
            self.es.run_single_experiment(self.expbase, 'cluster', pp)
            #self.es.run_single_experiment(self.expbase, 'missing', pp)
            
        else:
            self.es.run_single_experiment(self.expbase, 'local', pp)
    
    
    
    
    
    def generate_spike_times(self,freq,trainType,nOutput=1000,tmax=None,**gen_kw):
        """
        @param freq      : frequency of output 
        @param trainType : type of process to be generated (string). 
                Permissible values: 'poisson', 'inhomog_poisson', or 'gamma'  
        
        """
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
    
    def generate_spike_files(self):
        ttypes = ['poisson','inhomog_poisson']
        
        for ttype in ttypes:
            for freq in self.freqs:
                ss = self.generate_spike_times(freq, ttype, self.nsites, self.tstop)
                #print ss.shape()
                self._save_spiketime(ss, 'spikes_%s_freq%g'%(ttype,freq))
                # and for corresponding inhibitory spikes
                ssi = self.generate_spike_times(freq*self.freq_factor, ttype, self.nsites, self.tstop)
                #print ss.shape()
                self._save_spiketime(ssi, 'spikes_%s_freq%g'%(ttype,freq*self.freq_factor))

    def _save_input_sites(self,filename,idx):
        output = open(filename+'.pkl','wb')
        pickle.dump(idx,output)
        output.close()
        print "Saved input sites to %s.pkl"%filename
        
        
    def _load_input_sites(self,filename):
        print "Trying to load input sites from %s.pkl"%filename,
        #ss = np.loadtxt(filename)
        pkl_file = open('%s.pkl'%filename, 'rb')
        idx = pickle.load(pkl_file)
        pkl_file.close()
        print "... done!"
        return idx

    def get_input_sites(self,mindist=150,shuffles=10):
        
        __import__('Neuron') 
        neuronclass = sys.modules['Neuron']
        CellClass = getattr(neuronclass,self.celltype)
        self.cell = CellClass()
        # get relevant section for each type of cell we expect to encounter
        if self.celltype == 'SHStellate':
            section ='dendrite' 
        elif self.celltype == 'L5PC':
            section = 'apical'
        idx = self.cell.select_idx_posn_bydistance(section,mindist=mindist)[0]
        for i in range(shuffles):
            random.shuffle(idx)
            self._save_input_sites('sites_%s_mindist%g_shuffle%g'%(self.celltype,mindist,i),idx) 
        


####################################################################################
if __name__ == '__main__':
    print sys.argv
    if len(sys.argv) <= 1:
        print 'Need to specify an input mode'
    elif sys.argv[1] == 'inputs':
        Experiment_120614('a','b','c','d').generate_spike_files()
    elif sys.argv[1] == 'sites':
        celltype = sys.argv[2]
        Experiment_120614('a','b',celltype,'d').get_input_sites()
    elif sys.argv[1] == 'run':
        print 'run'
        wholecell,inputtype,celltype,spiketrain_type = sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]
        Experiment_120614(wholecell,inputtype,celltype,spiketrain_type).run_experiments()
        