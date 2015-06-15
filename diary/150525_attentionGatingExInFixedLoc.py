
import Neuron
from run_experiments import ExpSetup
import run_analysis
import NeuroTools.stgen as ntst 
import numpy as np
import pickle
import math
import os.path

expbase = '150525_attentionGatingExInFixedLoc'

"""
Purpose:

Model for having input of three different types to Pyramidal cell (L5PC) neuron:
- feedforward (input)
- recurrent 
- feedback (attention)
and for:
- keeping input rate constant
- varying recurrent and feedback rates
- with and without optogenetics

Intent: that the change in gain modulation will be feed back into changes 
in recurrent population -- which will in turn modulate this cell

Modifications:
- updated from 150519_attentionGatingExIn--> so that the same sites are chosen each time
- combined with 140718_replay_long, where sites and spikes times for each site were exported to pkl, and reloaded for each trial
- pkl for spike times and sites are now located in folder simdata/
- sites created in diary/140616_replay_input_apic.py and method get_input_sites()

To consider:
-- how many inputs for basal vs apic - correct for 
-- synchrony between input signals?

"""

stg = ntst.StGen()
es = ExpSetup()


tstop = 2500# 2200
light_dur = 1000
light_on = 1200


L5PC_areas = ['soma', 'apic','dend'] # NB: we don't include 'axon'

TOTAL_NUM_AREAS = {'L5PC': 3 }
areas = {'L5PC': L5PC_areas}

celltype = 'L5PC' 

spikefile = 'simdata/spikes_t10k_f%g_nsites%g' # filename for input spikes


gradients = np.arange(0.0,0.00101,0.0001)

class Helper():
    # First created in 140718_replay_long
    
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
    
    def checkSpikeFile(self,rate,nOutput,filebase,tmax=2500):
        # check if filename already exists
        if os.path.isfile(filebase%(rate,nOutput)+'.pkl'):
            print("Spike file already exists - "+filebase%(rate,nOutput))
            return
        # if not, creates it
        print("Spike file doesn't exist for rate %g, so creating it"%rate)
        ss = self.generate_spike_times(rate,'poisson',nOutput,tmax)
        self._save_spiketime(ss,filebase%(rate,nOutput))
        
    
    def _load_input_sites(self,filename):
        print "Trying to load input sites from %s.pkl"%filename,
        #ss = np.loadtxt(filename)
        pkl_file = open('%s.pkl'%filename, 'rb')
        idx = pickle.load(pkl_file)
        pkl_file.close()
        print "... done!"
        return idx
    



class AttentionGating():

    def get_populations(self,num_sites_ex,ratioEI,ratioRRSplit=0.5,percentFF=0.1,percentFB=0.1):
        """
        
        @params
            num_sites_ex
            ratioEI
            ratioRRSplit        What percentage of the recurrent connections should go to the apical dendrites
                                (1.0 --> all apic, 0.0 --> all basal; default = 0.5)
            percentFF           percentage of incoming connections that are feedforward* 
            percentFB           percentage of incoming connections that are feedback*
                                * recurrent connections are the remainder
        """
        self.num_sites_ex = num_sites_ex
        self.ratioEI = ratioEI
        self.num_sites_in = int(self.num_sites_ex*self.ratioEI)
        self.num_sites_total = self.num_sites_ex + self.num_sites_in
    
        ############ FOR Excitatory
        self.num_ff = int(0.1*num_sites_ex) # which will all be in the apical dendrites
        self.num_fb = int(0.1*num_sites_ex) # which will all be in the basal dendrites
        self.num_rr = self.num_sites_ex - self.num_ff - self.num_fb # which should be == int(0.8*num_sites_ex) 
                                                                    # and which will be split between the apical and dendrites
        # how we split our recurrent connections between apical and basal
        self.rr_apic = int(ratioRRSplit*self.num_rr)
        self.rr_basal = self.num_rr - self.rr_apic
        # so then the total number of synaptic sites on the apical and basal dendrites are:  
        self.num_apic = self.num_ff + self.rr_apic
        self.num_basal = self.num_fb + self.rr_basal
        ############ FOR Inhibitory
        self.i_num_ff = int(0.1*self.num_sites_in) # which will all be in the apical dendrites
        self.i_num_fb = int(0.1*self.num_sites_in) # which will all be in the basal dendrites
        self.i_num_rr = self.num_sites_in - self.i_num_ff - self.i_num_fb #int(0.8*self.num_sites_in) # which will be split between the apical and dendrites
        # how we split our recurrent connections between apical and basal
        self.i_rr_apic = int(ratioRRSplit*self.i_num_rr)
        self.i_rr_basal = self.i_num_rr - self.i_rr_apic
        # so then the total number of synaptic sites on the apical and basal dendrites are:  
        self.i_num_apic = self.i_num_ff + self.i_rr_apic
        self.i_num_basal = self.i_num_fb + self.i_rr_basal
        
        self.index_ff = 0
        self.index_rr = self.index_ff+self.num_ff
        self.index_fb = self.index_rr+self.num_rr
        self.index_ff_i = self.index_fb  +self.num_fb
        self.index_rr_i = self.index_ff_i+self.i_num_ff
        self.index_fb_i = self.index_rr_i+self.i_num_rr
        
    
    
    
    
    def get_spiketimes(self,rate,num_inputs):
        spikes = []
        for i in range(num_inputs):
            ss = stg.poisson_generator(rate,t_stop=tstop,array=True)
            spikes.append(ss)
        #print 'Spikes = ',spikes
        return spikes
    
    
    def _run_simulation(self,irr,irr_factor,freq_ff,freq_rr,freq_fb,adist,bdist,Jex,Jin,hh,shuffle_id):
        """
        @params
            irr       irradiance
            factor    NpHR factor
        """
        
        pp = {}
        pp['cell'] = ['Neuron',celltype]
        pp['cell_params'] = {}
        pp['tstart'] = 0
        pp['tstop'] = tstop
        
        # opsin
        chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':0,  'n_pulses':1}
        hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':0,  'n_pulses':1}
        
        pp['opsindict'] = {}
        pp['opsindict']['ChR'] =  {}
        for area in areas[celltype]:
            pp['opsindict']['ChR'][area] = chrdict    
        pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
        
        pp['opsindict']['NpHR'] =  {}
        for area in areas[celltype]:
            pp['opsindict']['NpHR'][area] = hrdict
        pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
        
        # general settings 
        pp['experiment_type'] = 'opsinonly' #TODO change
        pp['savedata'] = True 
        pp['logdata'] = False
        
        
        
        locnames = ['idx_e_%g'%i for i in range(self.num_sites_ex)] + ['idx_i_%g'%i for i in range(self.num_sites_in)]
        
        weights = np.concatenate((np.ones(self.num_sites_ex)*Jex ,  np.ones(self.num_sites_in)*Jin)) 
        spikefile_ff = hh._load_spiketimes(spikefile%(freq_ff,MAX_SITES))
        spikes_ff = spikefile_ff[self.index_ff:self.index_ff+self.num_ff]
        spikes_ff_i = spikefile_ff[self.index_ff_i:self.index_ff_i+self.i_num_ff]
        # spikes for recurrent population
        spikefile_rr = hh._load_spiketimes(spikefile%(freq_rr,MAX_SITES))
        spikes_rr = spikefile_rr[self.index_rr:self.index_rr+self.num_rr]
        spikes_rr_i = spikefile_rr[self.index_rr_i:self.index_rr_i+self.i_num_rr]
        # spikes for fb population
        spikefile_fb = hh._load_spiketimes(spikefile%(freq_fb,MAX_SITES))
        spikes_fb = spikefile_fb[self.index_fb:self.index_fb+self.num_fb]
        spikes_fb_i = spikefile_fb[self.index_fb_i:self.index_fb_i+self.i_num_fb]
        print len(spikefile_fb)
        
        
        
        
        pp['stim_spiketrains'] = True   
        spikelist = spikes_ff+spikes_rr+spikes_fb +  spikes_ff_i+spikes_rr_i+spikes_fb_i 
        """
        print "Length spiketrains"
        print len(spikelist)
        print len(spikes_ff)
        print len(spikes_rr)
        print len(spikes_fb)
        print len(spikes_ff_i)
        print len(spikes_rr_i)
        print len(spikes_fb_i)        
        print len(locnames)       
        """
        pp['spiketrains'] = [{'tstims': spikelist,  'locations':  locnames, 'weights':weights,  'el': 0.02}]
        
        idx_a = list(hh._load_input_sites('simdata/sites_%s_apic_mindist%g_shuffle%g'%(celltype,50,shuffle_id)))
        idx_d = list(hh._load_input_sites('simdata/sites_%s_dend_mindist%g_shuffle%g'%(celltype,50,shuffle_id)))
        otherbaseslabels = ['idx_e_%g'%i for i in range(self.num_sites_ex)]
        othersections = ['apic']*self.num_apic + ['dend']*self.num_basal
        # check that there are enough sites for our desired number of inputs
        while len(idx_a)<self.num_apic:
            print "Having to increase - sit A"
            idx_a += idx_a
        while len(idx_d)<self.num_basal:
            print "Having to increase - sit B"
            idx_d += idx_d
        print len(idx_a),len(idx_d)
        print self.num_apic, self.num_basal
        idxzip_e = [otherbaseslabels,np.concatenate((idx_a[:self.num_apic],idx_d[:self.num_basal]))]
        
        idx_a = list(hh._load_input_sites('simdata/sites_%s_apic_mindist%g_shuffle%g'%(celltype,150,shuffle_id)))
        idx_d = list(hh._load_input_sites('simdata/sites_%s_dend_mindist%g_shuffle%g'%(celltype,150,shuffle_id)))
        otherbaseslabels_i = ['idx_i_%g'%i for i in range(self.num_sites_in)]
        othersections_i = ['apic']*self.i_num_apic + ['dend']*self.i_num_basal
        # check that there are enough sites for our desired number of inputs
        while len(idx_a)<self.i_num_apic:
            print "Having to increase - sit C"
            idx_a += idx_a
        while len(idx_d)<self.i_num_basal:
            print "Having to increase - sit D", len(idx_d)
            print idx_d, type(idx_d)
            idx_d = idx_d + idx_d
        print len(idx_a),len(idx_d)
        print self.i_num_apic, self.i_num_basal
        idxzip_i = [otherbaseslabels_i,np.concatenate((idx_a[:self.i_num_apic],idx_d[:self.i_num_basal]))]
            
        #print olocs
        pp['mark_loc'] = {}
        pp['mark_loc']['names'] = ['mysoma']+['tmpA','tmpB'] #otherbaseslabels+otherbaseslabels_i
        pp['mark_loc']['sections'] = ['soma']+['tmpA','tmpB'] #othersections+othersections_i
        pp['mark_loc']['ids'] =  [(0,0.5)] + [{'idx':idxzip_e, 'sections':othersections}] + [{'idx':idxzip_i, 'sections':othersections_i}] 
        

    
        # set up recordings
        # - record voltage from soma
        pp['record_loc'] = {}
        pp['record_loc']['v'] = ['mysoma']
        
        pp['plot'] = { 'voltages': [ ['mysoma', 'v','k-'] ]} 
        
        pp.update({'expname':expbase,
                       'description':'_irr%.3f_factor%.2f_sites%g_ratioEI%.2f_ff%g_rr%g_fb%g_Jex%.1f_Jin%.1f'%(irr,factor,self.num_sites_total,self.ratioEI,freq_ff,freq_rr,freq_fb,Jex,Jin)})
        ##print pp['description']
        es.run_single_experiment(expbase, 'missing', pp)




ff_rates = [10,20,50,60,80]
scan_rates = range(20,101,20) # for fb and rr
irr = 0.001
factor = 0.6
g = 3.0 # set to be the same as balanceEI 
ratioEI = 0.25
EG = AttentionGating()
hh = Helper()
shuffle_id = 5
MAX_SITES = 10000
        
def runt():
    count = 0
    for num_sites in range(50,501,50):
        EG.get_populations(num_sites, ratioEI)            
        for ffr in ff_rates:
            hh.checkSpikeFile(ffr,MAX_SITES,spikefile)
            for rrr in scan_rates:
                hh.checkSpikeFile(rrr,MAX_SITES,spikefile)
                for fbr in scan_rates:
                    # check files exist for ff, rr, and fb rates
                    hh.checkSpikeFile(fbr,MAX_SITES,spikefile)
                    
                    for J in np.arange(0.1,0.31,0.1):
                        EG._run_simulation(irr,factor,ffr,rrr,fbr,[100,-1],[100,-1],J,J*-g,hh,shuffle_id)
                        count += 1
                        #return
    print count
    
def generate_save_locations():
    pass
        
def scan_locations():
    
    # load and run locations sets for cell 
    pass
    
        
runt()





