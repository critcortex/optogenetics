"""
Purpose: rerun several of the in vivo experiments so that n=50 rather than n=1


"""
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import pickle
import sys
import random

class Helper():
    
    def __init__(self):
        self.stg = ntst.StGen()
    
    def get_optdescript(self,irr,factor):
        opt_des = []
        if irr > 0 :
            opt_des.append('whole')
        else:
            opt_des.append('none')
            
        if irr > 0 and factor > 0:
            opt_des.append('whole')
        else:
            opt_des.append('none')
        return opt_des
    
    
    def get_spiketimes(self,rate,num_inputs,tstop):
        spikes = []
        for i in range(num_inputs):
            ss = self.stg.poisson_generator(rate,t_stop=tstop,array=True)
            spikes.append(ss)
        #print 'Spikes = ',spikes
        return spikes
    
    
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

    def get_input_sites(self,mindist=50,shuffles=10):
        """ 
        Extended from version appearing in 140616_replay
        """
        
        __import__('Neuron') 
        neuronclass = sys.modules['Neuron']
        CellClass = getattr(neuronclass,self.celltype)
        self.cell = CellClass()

        for section in self.cellsections:
            idx_a = self.cell.select_idx_posn_bydistance(section,mindist=mindist)[0]
            
            for i in range(shuffles):
                random.shuffle(idx_a)
                self._save_input_sites('shuffledSites_%s_%s_mindist%g_shuffle%g'%(self.celltype,section,mindist,i),idx_a) 

                        

class ExpStellate():
    
    def __init__(self):
        self.expbase = "150616_n50_stellate"
        
        self.helper = Helper()
        self.es = ExpSetup()
        
        self.cell =  ['Neuron', 'SHStellate']
        self.sections = ['dendrite']
     
     
        self.tstop = 2700
        self.freqs = [0.,50.,80.]
        self.freqs = range(10,201,10)
    
        self.light_on = 700
        self.light_dur = 1000
        #factors = [0]+[ 0.0625, 0.125,0.25,0.375,0.5,0.75,1.,2.,4.,8.]
        #factors = [0]+[ 0.0625, 0.125,0.25,0.375,0.5,0.75,1.,1.5,2.]
        self.factors = [ 0.0625, 0.125,0.25,0.375,0.5,0.75,1.,1.5,2.]
        #irrs = np.linspace(0.0025, 0.05, 11)
        #irrs = [0.003,0.007,0.012]
        self.irrs = [0.007,0.012]
        
        #iclamp_amps = np.arange(-1.0,1.7,0.2)
        #iclamp_amps = np.arange(1.8,3.1,0.2)
        self.iclamp_amps = np.arange(-1.,3.1,0.2)
        #iclamp_amps = np.arange(5.,10.1,1.)
        
        #Js = np.arange(1.,5.)
        self.Js = [2.]
        self.Ji = [6.]
        #nsite_range = range(20,51,20)
        # first value in each tuple is for # exc, 2nd value is for # inh
        self.nsite_range = [[40,14]]
        
        self.dist = [100,-1]
        """
        # NOTE - changed this from original setup, where sites were chosen as one or the other
        self.dist_vitro = [100,150]  ## SITE 1
        self.dist_vitro = [100.3,100.4] ## SITE 2
        self.dist_vitro = [304,305]  ## SITE 3
        """
        self.dist_vitro = [[100,150],[100.3,100.4],[304.,305.]]
        self.nsites = 1
        
        
        ##' TEST PARAMS
        self.factors = [0.5]
        self.freqs = [60]
        self.ntrials = 1
        
    
    def scan_locations_optogen_invivo(self):
        #iclamp_amps = [0]
        count = 0
        for freq in self.freqs:
            for factor in self.factors:
                for irr in self.irrs: 
                    for (j,J) in enumerate(self.Js):
                        for nsites in self.nsite_range:
                            for trial in range(self.ntrials):
                            
                                
                                pp = {}
                                # neuron model and params
                                pp['cell'] = ['Neuron', 'SHStellate']
                                pp['cell_params'] = {}
                                
                                # opsin
                                chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
                                hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
                                
                                pp['opsindict'] = {}
                                if irr > 0 :
                                    pp['opsindict']['ChR'] =  {'soma': chrdict,
                                                               'dendrite':chrdict}
                                    pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
                                else:
                                    pp['ChR_areas'] = {'none'      : [None]}
                                    
                                if irr > 0 and factor > 0:
                                    pp['opsindict']['NpHR'] =  {'soma': hrdict,
                                                                'dendrite':hrdict}
                                    pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
                                    
                                else:
                                    pp['NpHR_areas'] = {'none'      : [None]}
                                    
                                # general settings 
                                pp['experiment_type'] = 'opsinonly'
                                pp['savedata'] = True 
                                
                                pp['tstart'] = 0
                                pp['tstop'] = self.tstop
                                
                                mindist = self.dist[0]
                                idx = self.helper._load_input_sites('shuffledSites_%s_dendrite_mindist%g_shuffle%g'%(self.cell[1],mindist,trial))
                                nsites_e = nsites[0]
                                nsites_i = nsites[1]
                                
                                otherbaseslabels = ['idx_%g'%i for i in range(nsites_e)]+['idx_i%g'%i for i in range(nsites_i)]
                                othersections = ['dendrite']*(nsites_e+nsites_i)
                                idxzip = [otherbaseslabels,idx[:nsites_e+nsites_i]]
                                
                                print i, J, self.Ji
                                weights = np.concatenate((np.ones(nsites_e)*J,  np.ones(nsites_i)*self.Ji[j])) 
                                
                                
                                
                                pp['mark_loc'] = {}
                                pp['mark_loc']['names'] = ['mysoma']+otherbaseslabels
                                pp['mark_loc']['sections'] = ['soma']+othersections
                                pp['mark_loc']['ids'] = [(0,0.5)] +  [{'idx':idxzip, 'sections':othersections}]
                                
                                pp['record_loc'] = {}
                                pp['record_loc']['v'] = ['mysoma'] #+labels
                                pp['record_loc']['ina'] = ['mysoma']
                                pp['record_loc']['ik'] = ['mysoma']
                                
                                pp['stim_spiketrains'] = True
                                #TODO: fix spike times to be fixed??
                                print freq,nsites,self.tstop
                                pp['spiketrains'] = [{'tstims': self.helper.get_spiketimes(freq,nsites_e+nsites_i,self.tstop), 'locations': otherbaseslabels, 'weights': weights,  'el': 0.1}]
                                    
                                
                                
                                vplots_soma = [['mysoma','v','k']]
                    
                                #iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
                                pp['plot'] = {1:vplots_soma} 
                                    
                                
                                pp['num_threads'] = 1
                                                           
                                pp.update({'expname':self.expbase,
                                           'description':'irr%.3f_factor%.2f_freq%g_J%g_trial%g'%(irr,factor,freq,J,trial)})
                    
                                #self.es.run_single_experiment(self.expbase, 'missing', pp)
                                count += 1
                                self.es.run_single_experiment(self.expbase, 'local', pp)
                                return 
                                """
                                
                                NE = run_stimulation.NeuronExperiment()
                                ES = ExpSetup()
                                dp = ES.get_default_params()
                                dp.update(pp)
                                NE.main(dp)
                                return pp
                        
                            """
        print '%g jobs submitted'%count
    
    def scan_locations_optogen_invitro(self,loc):
        #iclamp_amps = [0]
        count = 0
    
        for factor in self.factors:
            for irr in self.irrs: 
                for Ia in self.iclamp_amps:
                    
                    pp = {}
                    # neuron model and params
                    pp['cell'] = ['Neuron', 'SHStellate']
                    pp['cell_params'] = {}
                    
                    # opsin
                    chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
                    hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
                    
                    pp['opsindict'] = {}
                    if irr > 0 :
                        pp['opsindict']['ChR'] =  {'soma': chrdict,
                                                   'dendrite':chrdict}
                    pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
                else:
                    pp['ChR_areas'] = {'none'      : [None]}
                    
                if irr > 0 and factor > 0:
                    pp['opsindict']['NpHR'] =  {'soma': hrdict,
                                                'dendrite':hrdict}
                    pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
                    
                else:
                    pp['NpHR_areas'] = {'none'      : [None]}
                    
                # general settings 
                pp['experiment_type'] = 'opsinonly'
                pp['savedata'] = True # False #True
                
                pp['tstart'] = 0
                pp['tstop'] = self.tstop
                
                
                labels = ['stim%g'%i for i in range(self.nsites)]
                
                
                pp['mark_loc'] = {}
                pp['mark_loc']['names'] = ['mysoma']+labels
                pp['mark_loc']['sections'] = ['soma']+['dendrite']*self.nsites
                pp['mark_loc']['ids'] = [(0,0.5)] + [('select_section_posn_bydistance',{'sectionarea':'dendrite','mindist':dists[0],'maxdist':dists[1]}) for (i,dists) in enumerate(self.dist_vitro)] 
                
                
                pp['record_loc'] = {}
                pp['record_loc']['v'] = ['mysoma'] #+labels
                pp['record_loc']['ina'] = ['mysoma']
                pp['record_loc']['ik'] = ['mysoma']
                
                pp['stim_iclamp'] = True
                #print labels[0]
                #pp['spiketrains'] = [{'tstims': get_spiketimes(freq,nsites), 'locations': labels, 'weights':np.ones(nsites)*J,  'el': 0.1}]
                if loc == 'soma':
                    stimloc = 'mysoma'
                else:
                    stimloc = labels[-1]
                pp['iclamp'] = [{'tstim':0,  'location': stimloc, 'amp':Ia, 'duration':self.tstop}]
                       
                
                
                vplots_soma = [['mysoma','v','k']]
    
                #iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
                pp['plot'] = {1:vplots_soma} 
                    
                
                pp['num_threads'] = 1
                                           
                pp.update({'expname':self.expbase,
                           'description':'irr%.3f_factor%.2f_I%.2f_stimloc_%s'%(irr,factor,Ia,stimloc)})
    
                self.es.run_single_experiment(self.expbase, 'missing', pp)
                #es.run_single_experiment(expbase,'names',pp)
                count += 1
                #es.run_single_experiment(expbase, 'local', pp)
                #return 
                """
                
                NE = run_stimulation.NeuronExperiment()
                ES = ExpSetup()
                dp = ES.get_default_params()
                dp.update(pp)
                NE.main(dp)
                return pp
        
            """
        print '%g jobs submitted'%count
       

class ExpL5PC():
    
    
    def __init__(self):
        
        self.cell =  ['Neuron', 'L5PC']
        self.sections = ['apic','dend']
        self.helper = Helper()
        self.es = ExpSetup()
        
        self.expbase = "150616_n50_L5PC"



"""

# Create sites for Stellate cells 
expStellate = ExpStellate()
expStellate.helper.celltype = expStellate.cell[1]
expStellate.helper.cellsections = expStellate.sections

#Create sites for L5PC cells
expL5PC = ExpL5PC()
expL5PC.helper.celltype = expL5PC.cell[1]
expL5PC.helper.cellsections = expL5PC.sections

"""        
expSt = ExpStellate()
expSt.scan_locations_optogen_invivo()
        