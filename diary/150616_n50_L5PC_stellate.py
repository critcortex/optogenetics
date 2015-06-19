"""
Purpose: rerun several of the in vivo experiments so that n=50 rather than n=1


"""
from run_experiments import ExpSetup
import run_analysis
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
                
    def create_opsindict(self,irr,factor,chrdict,hrdict,pp,sections):
        pp['opsindict']={}
        if irr > 0 :
            pp['opsindict']['ChR'] =  {'soma': chrdict}
            for subdomain in sections:
                pp['opsindict']['ChR'][subdomain] = chrdict
            pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
        else:
            pp['ChR_areas'] = {'none'      : [None]}
            
        if irr > 0 and factor > 0:
            pp['opsindict']['NpHR'] =  {'soma': hrdict}
            for subdomain in sections:
                pp['opsindict']['NpHR'][subdomain] = hrdict
            pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
            
        else:
            pp['NpHR_areas'] = {'none'      : [None]}
        return pp

    def plot_experiment(self,expname, expdescript, explabel, variables, exp_params={},extractSpikes=True):
        """
        This assumes that the analysis has already been performed
        
        @param
            expname
            expdescript     array of experiment names to go in one set
            explabel        
            variables       array of the variable we're looping over i.e. factors
            exp_params      dictionary of setup params for that experiment i.e. {"tstart":100,"tstart_bg":50} etc
        """
        af = run_analysis.AnalyticFrame()
        # set tstart,tstop for each experiment - technically don't have to do this, as the analysis has already been performed
        af.update_params(exp_params)
        # set up data structure for each experiment
        af.populate_expset(expname,expdescript,explabel, variables)
        #TODO: extract spikes
        if extractSpikes:
            print "Trying to extract spikes ------"
            af.submenu_extractSpikes()
        # return AnalyticFrame for further plotting, etc.
        #af.submenu_save()
        return af
    
    def plot_trialSet(self,basenames,trialInstances,variables,trials,trialLabels,labels,var_format,exp_params={},extractSpikes=True):
        """
        
        
        """
        af = run_analysis.AnalyticFrame()
        # set tstart,tstop for each experiment - technically don't have to do this, as the analysis has already been performed
        af.update_params(exp_params)
        # set up data structure for each experiment
        af.populate_trialset(basenames, trialInstances, variables, trials, trialLabels, labels, var_format)
        
        if extractSpikes:
            print "Trying to extract spikes ------"
            af.submenu_extractSpikes()
            
        return af
 

class ExpStellate():
    
    def __init__(self):
        self.expbase = "150616_n50_stellate"
        
        self.helper = Helper()
        self.es = ExpSetup()
        
        self.cell =  ['Neuron', 'SHStellate']
        self.sections = ['dendrite']
     
     
        self.tstop = 2700
        self.light_on = 700
        self.light_dur = 1000
        
        self.freqs = [0.,50.,80.]
        self.freqs = range(10,201,10)
    
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
        self.dist_vitro = [100.3,100.4] ## SITE 2 <---------- one analysed in paper
        self.dist_vitro = [304,305]  ## SITE 3
        """
        self.dist_vitro = [[100,150],[100.3,100.4],[304.,305.]]
        self.nsites = 1
        
        
        ##' TEST PARAMS
        #self.factors = [0.5]
        #self.freqs = [60]
        #self.irrs = [0.1,1.]
        self.ntrials =  5
        
    
    def run_invivo(self):
       
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
                                """
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
                                """
                                pp = self.helper.create_opsindict(irr, factor, chrdict, hrdict, pp, self.sections)
                                    
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
                                
                                #print i, J, self.Ji
                                weights = np.concatenate((np.ones(nsites_e)*J,  np.ones(nsites_i)*self.Ji[j])) 
                                
                                #print len(otherbaseslabels)
                                #print len(idx)
                                #print idx
                                
                                pp['mark_loc'] = {}
                                pp['mark_loc']['names'] = ['mysoma']+['dict']
                                pp['mark_loc']['sections'] = ['soma']+['dict']
                                pp['mark_loc']['ids'] = [(0,0.5)] +  [{'idx':idxzip, 'sections':othersections}]
                                
                                pp['record_loc'] = {}
                                pp['record_loc']['v'] = ['mysoma'] #+labels
                                pp['record_loc']['ina'] = ['mysoma']
                                pp['record_loc']['ik'] = ['mysoma']
                                
                                pp['stim_spiketrains'] = True
                                #TODO: fix spike times to be fixed??
                                #print freq,nsites,self.tstop
                                pp['spiketrains'] = [{'tstims': self.helper.get_spiketimes(freq,nsites_e+nsites_i,self.tstop), 'locations': otherbaseslabels, 'weights': weights,  'el': 0.1}]
                                    
                                
                                
                                vplots_soma = [['mysoma','v','k']]
                    
                                #iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
                                #pp['plot'] = {1:vplots_soma} 
                                    
                                
                                pp['num_threads'] = 1
                                                           
                                pp.update({'expname':self.expbase,
                                           'description':'irr%.3f_factor%.2f_freq%g_J%g_trial%g'%(irr,factor,freq,J,trial)})
                    
                                self.es.run_single_experiment(self.expbase, 'missing', pp)
                                count += 1
                                #self.es.run_single_experiment(self.expbase, 'local', pp)
                                
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
    
    def run_invitro(self,loc):
        #iclamp_amps = [0]
        count = 0
    
        for factor in self.factors:
            for irr in self.irrs: 
                for Ia in self.iclamp_amps:
                    for trial in range(0,self.ntrials):
                    
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
                        pp['mark_loc']['names'] = ['mysoma']+['stimloc']
                        pp['mark_loc']['sections'] = ['soma']+['dendrite']
                        pp['mark_loc']['ids'] = [(0,0.5)] + [('select_section_posn_bydistance',{'sectionarea':'dendrite','mindist':dists[0],'maxdist':dists[1]}) for (i,dists) in enumerate(self.dist_vitro)] 
                        #TODO: update
                        
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
                            stimloc = 'stimloc'
                        pp['iclamp'] = [{'tstim':0,  'location': stimloc, 'amp':Ia, 'duration':self.tstop}]
                               
                        
                        
                        vplots_soma = [['mysoma','v','k']]
            
                        #iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
                        pp['plot'] = {1:vplots_soma} 
                            
                        
                        pp['num_threads'] = 1
                                                   
                        pp.update({'expname':self.expbase,
                                   'description':'irr%.3f_factor%.2f_I%.2f_stimloc_%s_trial%g'%(irr,factor,Ia,stimloc,trial)})
            
                        self.es.run_single_experiment(self.expbase, 'missing', pp)
                        #es.run_single_experiment(expbase,'names',pp)
                        count += 1
                        #es.run_single_experiment(expbase, 'local', pp)
                        #return 

            print '%g jobs submitted'%count
       
    def analyse_invivo(self):
        description = ['whole','whole']
        expParams = {'tstart':self.light_on,'tstop': self.light_on+self.light_dur,
                              'tstart_bg': 50,'tstop_bg':self.light_on,
                              'tstart_post':self.light_on+self.light_dur,'tstop_post':self.tstop}

        for irr in self.irrs:
            
            instances = ["irr%.3f"%irr+"_factor%.2f"%factor+"_freq%g_J2_trial%s_NpHR_whole_ChR_whole" for factor in self.factors]
            trials = range(self.ntrials)
            trialLabels = ['%.2f'%f for f in self.factors]
            labels = ["freq%g"%f for f in self.freqs]
            var_format = '%g'

            af = self.helper.plot_trialSet(self.expbase,instances,self.freqs,trials,trialLabels,labels,var_format,exp_params=expParams,extractSpikes=True)
            
            af.submenu_runFI()
            for exp in af.experimentset:
                exp.calculate_responses('FI')
                exp.calculate_responses('FI_bg')
                exp.calculate_responses('FI_post')
            af.submenu_save()
            af.submenu_print()
            for exp in af.experimentset:
                print exp
                exp.load_experiments()
                exp.collate_results()
                print exp.results
            af.submenu_plot(10, self.expbase+'FIbg_gain_irr%.3f_varyFactor_trialled_invivo_exp%s%s_'%(irr,description[0],description[1]))            
            
    
    def analyse_invitro(self):
        pass

class ExpL5PC():
    
    
    def __init__(self):
        
        self.cell =  ['Neuron', 'L5PC']
        self.sections = ['apic','dend']
        self.helper = Helper()
        self.es = ExpSetup()
        
        self.expbase = "150616_n50_L5PC"
        
        self.tstop = 2700
        self.light_on = 700
        self.light_dur = 1000
        
        self.ntrials = 50
        # values for nsites are: [#apical_e, #apical_i, #basal_e, #basal_i]
        self.nsites = [[50,40,10,8]]
        self.mindist = [[100,150]]
        
        self.Js = [2.]
        self.Ji = [1.]
        
        self.irrs = [0.002,1.]
        self.factors = [0.001,0.125,0.25,0.5,0.75,1.]
        self.freqs = range(2,15,2)+range(15,151,5)
        
        
        
        ### DEBUG PARAMS
        self.ntrials = 5
        #self.factors = [0.5]
        #self.freqs = [80]
        self.irrs = [0.002]
        
        
        
    def run_invivo(self):
        count = 0
        for freq in self.freqs:
            for factor in self.factors:
                for irr in self.irrs: 
                    for (j,J) in enumerate(self.Js):
                        for nsites in self.nsites:
                            for dist in self.mindist:
                                for trial in range(self.ntrials):
                                
                                    pp = {}
                                    # neuron model and params
                                    pp['cell'] = ['Neuron', 'L5PC']
                                    pp['cell_params'] = {}
                                    
                                    # opsin
                                    chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
                                    hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
                                    
                                    pp = self.helper.create_opsindict(irr, factor, chrdict, hrdict, pp, self.sections)
                                    print pp
                                        
                                    # general settings 
                                    pp['experiment_type'] = 'opsinonly'
                                    pp['savedata'] = True 
                                    
                                    pp['tstart'] = 0
                                    pp['tstop'] = self.tstop
                                    
                                    #print j, J, self.Ji
                                    Ji = self.Ji[j]
                                    
                                    mindist = dist[0]
                                    mindist_i = dist[1]
                                    
                                    idx_a = self.helper._load_input_sites('shuffledSites_%s_apic_mindist%g_shuffle%g'%(self.cell[1],mindist,trial))[:nsites[0]]
                                    idx_b = self.helper._load_input_sites('shuffledSites_%s_dend_mindist%g_shuffle%g'%(self.cell[1],mindist,trial))[:nsites[1]]
                                    #print idx_a
                                    #print idx_b
                                    idx_ai = self.helper._load_input_sites('shuffledSites_%s_apic_mindist%g_shuffle%g'%(self.cell[1],mindist_i,trial))[:nsites[2]]
                                    idx_bi = self.helper._load_input_sites('shuffledSites_%s_dend_mindist%g_shuffle%g'%(self.cell[1],mindist_i,trial))[:nsites[3]]
                                    #print idx_ai
                                    #print idx_bi
                                    nsites_e = nsites[0]+nsites[2] #excitatory
                                    nsites_i = nsites[1]+nsites[3] #inhibitory
                                    nsites_a = nsites[0]+nsites[2] #apical
                                    nsites_b = nsites[1]+nsites[3] #basal
                                     
                                    
                                    otherbaseslabels = ['apic_%g'%i for i in range(nsites[0])]+['apic_i%g'%i for i in range(nsites[2])]+['basal_%g'%i for i in range(nsites[1])]+['basal_i%g'%i for i in range(nsites[3])]
                                    othersections = ['apic']*(nsites_a) + ['dend']*(nsites_b)
                                    idx =  np.concatenate((idx_a,  idx_ai,  idx_b,  idx_bi)) 
                                    
                                    #return
                                    idxzip = [otherbaseslabels,idx]
                                    #print idxzip
                                    

                                    weights = np.concatenate((np.ones(nsites_e)*J,  np.ones(nsites_i)*self.Ji[j])) 
                                    
                                    pp['mark_loc'] = {}
                                    pp['mark_loc']['names'] = ['mysoma']+['dict']
                                    pp['mark_loc']['sections'] = ['soma']+['dict']
                                    pp['mark_loc']['ids'] = [(0,0.5)] +  [{'idx':idxzip, 'sections':othersections}]
                                    
                                    pp['record_loc'] = {}
                                    pp['record_loc']['v'] = ['mysoma'] #+labels
                                    pp['record_loc']['ina'] = ['mysoma']
                                    pp['record_loc']['ik'] = ['mysoma']
                                    
                                    pp['stim_spiketrains'] = True
                                    pp['spiketrains'] = [{'tstims': self.helper.get_spiketimes(freq,nsites_e+nsites_i,self.tstop), 'locations': otherbaseslabels, 'weights': weights,  'el': 0.1}]
                                        
                                    vplots_soma = [['mysoma','v','k']]
                                    #pp['plot'] = {1:vplots_soma} 
                                    
                                    pp['num_threads'] = 1
                                                           
                                    pp.update({'expname':self.expbase,
                                               'description':'irr%.3f_factor%.2f_freq%g_J%g_Ji%g_nsites%g_%g_basal%g_%g_trial%g'%(irr,factor,freq,J,Ji,nsites[0],nsites[1],nsites[2],nsites[3],trial)})
                                    
                                    self.es.run_single_experiment(self.expbase, 'missing', pp)
                                    count += 1
                                    #self.es.run_single_experiment(self.expbase, 'local', pp)
                                    #return 
                                    
        print '%g jobs submitted'%count

    def analyse_invivo(self):
        description = ['whole','whole']
        expParams = {'tstart':self.light_on,'tstop': self.light_on+self.light_dur,
                              'tstart_bg': 50,'tstop_bg':self.light_on,
                              'tstart_post':self.light_on+self.light_dur,'tstop_post':self.tstop}
        J = self.Js[0]
        Ji = self.Ji[0]
        nsites = self.nsites[0]

        for tt in range(1, self.ntrials+1):
            for irr in self.irrs:
                
                instances = ["irr%.3f"%irr+"_factor%.2f"%factor+"_freq%g"+"_J%g_Ji%g_nsites%g_%g_basal%g_%g"%(J,Ji,nsites[0],nsites[1],nsites[2],nsites[3])+"_trial%s_NpHR_whole_ChR_whole" for factor in self.factors]
                trials = range(tt)
                trialLabels = ['%.2f'%f for f in self.factors]
                labels = ["freq%g"%f for f in self.freqs]
                var_format = '%g'
    
                af = self.helper.plot_trialSet(self.expbase,instances,self.freqs,trials,trialLabels,labels,var_format,exp_params=expParams,extractSpikes=True)
                
                af.submenu_runFI()
                for exp in af.experimentset:
                    exp.calculate_responses('FI')
                    exp.calculate_responses('FI_bg')
                    exp.calculate_responses('FI_post')
                af.submenu_save()
                af.submenu_print()
                for exp in af.experimentset:
                    exp.load_experiments()
                    exp.collate_results()
                af.submenu_print()    
                af.submenu_plot(10, self.expbase+'FIbg_gain_irr%.3f_varyFactor_trialled%g_invivo_exp%s%s_'%(irr,tt,description[0],description[1]))            
                af.submenu_plot(5, self.expbase+'FI_gain_irr%.3f_varyFactor_trialled%g_invivo_exp%s%s_'%(irr,tt,description[0],description[1]))            
            


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
#expSt = ExpStellate()
#expSt.scan_locations_optogen_invivo()

#expPC = ExpL5PC()
#expPC.run_invivo()

if __name__ == '__main__':
    #print sys.argv
    if len(sys.argv)<=1:
        print("No input entered")
    else:
        cell =  sys.argv[1]
        if sys.argv[1] == 'L5PC':
            cell = ExpL5PC()
        else:
            cell = ExpStellate()
            
        proctype = sys.argv[2]
        runtype = sys.argv[3]
        
        getattr(cell,'%s_%s'%(proctype,runtype))()
    
