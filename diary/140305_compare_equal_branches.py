import sys
import run_analysis 
import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import run_stimulation
import pylab
import file_io as fio


class Experiment_140305:
    
    
    
    def __init__(self,whole,stimulus):
        self.es = ExpSetup()
        self.whole = (whole == 'whole')
        if not self.whole:
            self.explist = ['partialSame']
        else:
            self.explist = ['whole']
            
        self.stimulus = stimulus
        self.expbase = '140305_compare_equal_branches'
    
        self.stim_location_fracts = [['peripheral',[0.5,1.0]],
                 ['distributed',[0,1.0]],
                 ['proximal',[0,.5]]]
    
        self.tstop = 3500
    
        self.freq = 50.
        self.dist_frac_sections = 50.

        # NB: collection was incorrect, as (4,30,1) should have been (4,30,2)
        #self.collection_same_total = [(1,1,124),(1,2,7),(2,2,6),(2,7,3),(2,1,62),(2,61,2),(4,5,3),(4,2,5),(4,30,1),(4,1,31)]
        # NB: this collection was used for setting up and testing runs
        #self.collection_same_total = [(2,2,6)]
        # NB: this collection was added later, after additional dendritic forms were identified
        #self.collection_same_total = [(31,1,4),(31,3,2),(62,1,2),(11,1,11),(11,10,2)]
        #self.collection_same_total = [(1,1,124),(1,2,7),(2,2,6),(2,7,3),(2,1,62),(2,61,2),(4,5,3),(4,2,5),(4,30,2),(4,1,31)]
        ### NB. this collection is the correct collection
        self.collection_same_total = [(1,1,124),(1,2,7),(2,2,6),(2,7,3),(2,1,62),(2,61,2),(4,5,3),(4,2,5),(4,30,2),(4,1,31),(31,1,4),(31,3,2),(62,1,2),(11,1,11),(11,10,2),(18,1,7),(18,2,3),(18,6,2)]
        #self.collection_same_total = [(4,30,2)]
        self.collection_same_total = [(124,1,1),(1,1,124)]
        self.collection_same_total = [(2,2,6),(2,7,3),(2,1,62),(2,61,2),(4,5,3),(4,2,5),(4,30,2),(4,1,31)]
        #self.collection_same_total = [(1,1,124),(4,30,2),(124,1,1)]    
        #self.collection_same_total = [(1,2,7) ,(2,2,6),(2,1,62),(2,61,2),(4,5,3),(4,2,5),(4,1,31)]
    
        self.collection_branch_total_large = [(3,7),(2,10),(4,6),(10,4),(32,3)]
        self.collection_branch_total_small = [(5,3),(2,5),(30,2),(1,31)]
        
    
        self.light_on = 1000
        self.light_dur = 1500
        self.factors = [0.25,1.,4.,8.]+[0.,0.0625,0.125]+[0.01,0.03]
        self.irrs = np.arange(0.5,2.6,0.5)
        self.factors = [0,0.125,0.25,1.,4.,8.]
        
        self.clamp_amps = np.arange(0.2,5.,0.2)
        self.clamp_amps = [1.5]
        
        
        
        self.factors = [0.5,2,1.5,1.,0.75,1.25]
        #self.factors = [0,0.5,2,1.5,1.,0.75,1.25]
        
        #self.irrs = np.arange(0.05,.36,0.05)
        self.irrs = [0.05]#,0.02]
        self.freqs = range(10,101,10)
        #self.freqs = range(1,10,1)+[15,20]
        self.Js = [2.]# ,2.5]
        #self.Js = [10.]# ,2.5]
        #self.Js = [1.]# ,2.5]
        #self.Js =[1000,10000]
        #self.Js =[10000]
        #self.Js =[15000]
        
        """
        self.collection_same_total = [(124,1,1)]
        self.Js =[15000]
        
        self.collection_same_total = [(31,1,4),(31,3,2),(62,1,2)]
        self.Js =[10000]
        self.freqs = range(10,101,10)
        """
        
        # MOD 150805:
        self.factors = [0.5,0.75,1.,1.25,1.5,2.]
        self.freqs = range(10,101,10)
        self.irrs = [0.05]
        self.collection_same_total = [(2,2,6),(2,7,3),(4,5,3)]
        
        # MOD 150811:
        self.factors = [0.5,0.75,1.,1.25,1.5,2.]
        #self.freqs = range(10,101,10)
        self.irrs = [ a*b for a in [0.01,0.1,1.,10.] for b in [1.,2.,5.]]
        #self.irrs = [ a*b for a in [0.01] for b in [1.,2.,5.]]
        self.collection_same_total = [(2,2,6)]
        self.factors = [0.0,0.5] #0.12,0.25,0.5,0.75,1.]
        # MOD 150813
        #self.Js =[10000]
        #self.freqs = range(10,101,10)
        #self.freqs = [10,50,100]
        #self.factors = [0.0]
        #self.irrs = [ a*b for a in [0.01,0.1,1.,10.] for b in [1.]]
        self.factors = [0.1*a for a in range(11)] + [2.,5.,10.]
        self.irrs = [ a*b for a in [1.0] for b in [1.]]
        #self.Js =[10,100,1000]
        self.collection_same_total = [(1,1,124),(1,2,7),(2,2,6),(4,5,3),(4,1,31),(31,1,4)]
    
    def get_dendritic_colors(self,num_levels):
        """
        soma = level0 will always be black
        """
        num_levels = num_levels +2 # to account for soma and that we don't want to go to complete white
        return ['#333333','#666666','#888888','#aaaaaa','#dddddd']
        
    
    
    def get_dend0_child_locations(self,num_level):
        labels = []
        locs = []
        for i in range(1,num_level+1):
            tmp = 'dend'+'_'.join(['0']*i)
            labels.append(tmp)
            locs.append((tmp,0.5))
        return labels,locs
    
    def get_spiketrain_locations_dend0(self,nseg=10):
        treebase = []
        labels = []
        locs = []
        for i in range(0,nseg):
            treebase.append('dend0')
            tmp = 'dend0'
            labels.append(tmp)
            loc = 1./20*(i*2 +1) 
            locs.append((tmp,loc))
        return treebase,labels,locs
    
    def get_spiketrain_locations_nloc(self,posloc,nseg=10):
        treebase = []
        labels = []
        locs = []
        
        for i in range(0,nseg):
            treebase.append('dend0')
            tmp = 'dend'+'_'.join(['0']*(posloc+1))
            labels.append('my%s_%g'%(tmp,i))
            loc = 1./20*(i*2 +1) 
            locs.append((tmp,loc))
        return treebase,labels,locs
    
    def get_spiketrain_locations_nloc_soma(self,nseg=10):
        treebase = []
        labels = []
        locs = []
        
        for i in range(0,nseg):
            treebase.append('soma')
            tmp = 'soma'
            labels.append('my%s_%g'%(tmp,i))
            loc = 1./20*(i*2 +1) 
            locs.append((tmp,loc))
        return treebase,labels,locs
    
    
    def run_collection_total(self):
        
        #try:
            getattr(self,'run_collection_%s'%self.stimulus)(self.collection_same_total,self.whole)
        #except:
        #    print 'Cannot act on no stimulus'
        
        
        
    def run_collection_generic(self,pp,tree,irr,factor,whole,other_tag_locations=None,runon=True):
        nb,nc,nl = tree
        # neuron model and params
        pp['cell'] = ['Neuron', 'FractalNeuron']
        pp['cell_params'] = {'num_base': nb,
                      'num_split': nc,
                      'num_levels':nl,
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
        # opsin
        show_levels = np.min([nl,5])
        labels,locs = self.get_dend0_child_locations(show_levels)
        #print locs
        chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
        hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
        #chrdict =  {'exp':5e-2, 'irradiance' :irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
        #hrdict =  {'exp':5e-2, 'irradiance' :factor*irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
        
        
        pp['opsindict'] = {}
        if irr > 0 and whole:
            print 'Case 1'
            pp['opsindict']['ChR'] =  {'soma': chrdict}
            for i in range(nb):
                pp['opsindict']['ChR']['dend%g'%i] = chrdict
                
            pp['opsindict']['NpHR'] =  {'soma': hrdict}
            for i in range(nb):
                pp['opsindict']['NpHR']['dend%g'%i] = hrdict
            pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
            pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
        elif irr > 0 and not whole:
            if nb == 1:
                print "Does not make sense to run partialSame for trees where nb=1"
                return
            print 'Case 2'
            
            
            opsin_index = 0
            
            pp['opsindict']['ChR'] =  {'soma': chrdict}
            pp['opsindict']['NpHR'] =  {'soma': hrdict}
            pp['opsindict']['ChR']['dend%g'%opsin_index] = chrdict
            pp['opsindict']['NpHR']['dend%g'%opsin_index] = hrdict
            pp['NpHR_areas'] = {'partialSame':['dend%g'%opsin_index,'soma']}
            pp['ChR_areas'] = {'partialSame':['dend%g'%opsin_index,'soma']}
        elif irr==0 and whole:
            print 'Case 3'
            pp['NpHR_areas'] = {'none'      : [None]}
            pp['ChR_areas'] = {'none'      : [None]}
        else: # irr = 0 and partial only --> so why bother simulating
            print 'Case 4'
            return
        
        if factor ==0:
            pp['NpHR_areas'] = {'none'      : [None]}
            pp['opsindict']['NpHR'] = {}
        
        
        # general settings 
        pp['experiment_type'] = 'opsinonly'
        pp['savedata'] = True # False #True
        
        pp['tstart'] = 0
        pp['tstop'] = self.tstop
        
        otherbaseslabels,othersections = [],[]
        if nb > 1:
            otherbaseslabels = ['dend1']
            othersections    = ['dend1']
            
        if other_tag_locations is not None:
            obases,olabels,olocs = other_tag_locations
            #print olocs
        else:
            obases,olabels,olocs = [],[],[]
            
        #print olocs
        pp['mark_loc'] = {}
        pp['mark_loc']['names'] = ['mysoma']+labels+otherbaseslabels+olabels
        pp['mark_loc']['sections'] = ['soma']+['dend0']*nl+othersections+obases
        pp['mark_loc']['ids'] = [(0,0.5)] + [('get_section_byname',{'sectioname':loc[0],'subdomain':'dend0'}) for loc in locs]+[('get_section_byname',{'sectioname':loc,'subdomain':'dend1'}) for loc in otherbaseslabels]+[(0,loc[1]) for loc in olocs ]
        #print pp['mark_loc']['ids']
        
        pp['record_loc'] = {}
        pp['record_loc']['v'] = ['mysoma']+labels+otherbaseslabels
        #pp['record_loc']['ina'] = ['mysoma']+labels+otherbaseslabels
        #pp['record_loc']['ik'] = ['mysoma']+labels+otherbaseslabels
        #pp['record_loc']['ik'] = labels+otherbaseslabels
        #pp['record_loc']['ica'] = labels+otherbaseslabels
        
        
        vplots_soma = [['mysoma','v','k']]
        """
        iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
        vplots = [['mysoma','v','k']]
        iplots_k = [['mysoma','ik','b']]
        iplots_na = [['mysoma','ina','g']]
        colors = self.get_dendritic_colors(show_levels)
        for i in range(show_levels):
            vplots.append([labels[i],'v',colors[i]])
            iplots_k.append([labels[i],'ik',colors[i]])
            iplots_na.append([labels[i],'ina',colors[i]])
        #print vplots
        if nb>1:
            vplots.append([otherbaseslabels[0],'v','b-'])
        iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
        """
        pp['plot'] = {1:vplots_soma}
        """ #2:vplots,
                      #3:iplots_soma,
                      #4:iplots_k,
                      #5:iplots_na}
        """
        pp['num_threads'] = 1
        if runon:
            #self.es.run_single_experiment(self.expbase, 'cluster', pp)
            self.es.run_single_experiment(self.expbase, 'missing', pp)
            
        else:
            self.es.run_single_experiment(self.expbase, 'local', pp)
            
            
            
    def run_collection_generic_highexp(self,pp,tree,irr,factor,whole,other_tag_locations=None,runon=True):
        nb,nc,nl = tree
        # neuron model and params
        pp['cell'] = ['Neuron', 'FractalNeuron']
        pp['cell_params'] = {'num_base': nb,
                      'num_split': nc,
                      'num_levels':nl,
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
        # opsin
        show_levels = np.min([nl,5])
        labels,locs = self.get_dend0_child_locations(show_levels)
        #print locs
        #chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
        #hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
        chrdict =  {'exp':5, 'irradiance' :irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
        hrdict =  {'exp':5, 'irradiance' :factor*irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
        
        
        pp['opsindict'] = {}
        if irr > 0 and whole:
            print 'Case 1'
            pp['opsindict']['ChR'] =  {'soma': chrdict}
            for i in range(nb):
                pp['opsindict']['ChR']['dend%g'%i] = chrdict
                
            pp['opsindict']['NpHR'] =  {'soma': hrdict}
            for i in range(nb):
                pp['opsindict']['NpHR']['dend%g'%i] = hrdict
            pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
            pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
        elif irr > 0 and not whole:
            if nb == 1:
                print "Does not make sense to run partialSame for trees where nb=1"
                return
            print 'Case 2'
            
            
            opsin_index = 0
            
            pp['opsindict']['ChR'] =  {'soma': chrdict}
            pp['opsindict']['NpHR'] =  {'soma': hrdict}
            pp['opsindict']['ChR']['dend%g'%opsin_index] = chrdict
            pp['opsindict']['NpHR']['dend%g'%opsin_index] = hrdict
            pp['NpHR_areas'] = {'partialSame':['dend%g'%opsin_index,'soma']}
            pp['ChR_areas'] = {'partialSame':['dend%g'%opsin_index,'soma']}
        elif irr==0 and whole:
            print 'Case 3'
            pp['NpHR_areas'] = {'none'      : [None]}
            pp['ChR_areas'] = {'none'      : [None]}
        else: # irr = 0 and partial only --> so why bother simulating
            print 'Case 4'
            return
        
        if factor ==0:
            pp['NpHR_areas'] = {'none'      : [None]}
            pp['opsindict']['NpHR'] = {}
        
        
        # general settings 
        pp['experiment_type'] = 'opsinonly'
        pp['savedata'] = True # False #True
        
        pp['tstart'] = 0
        pp['tstop'] = self.tstop
        
        otherbaseslabels,othersections = [],[]
        if nb > 1:
            otherbaseslabels = ['dend1']
            othersections    = ['dend1']
            
        if other_tag_locations is not None:
            obases,olabels,olocs = other_tag_locations
            #print olocs
        else:
            obases,olabels,olocs = [],[],[]
            
        #print olocs
        pp['mark_loc'] = {}
        pp['mark_loc']['names'] = ['mysoma']+labels+otherbaseslabels+olabels
        pp['mark_loc']['sections'] = ['soma']+['dend0']*nl+othersections+obases
        pp['mark_loc']['ids'] = [(0,0.5)] + [('get_section_byname',{'sectioname':loc[0],'subdomain':'dend0'}) for loc in locs]+[('get_section_byname',{'sectioname':loc,'subdomain':'dend1'}) for loc in otherbaseslabels]+[(0,loc[1]) for loc in olocs ]
        #print pp['mark_loc']['ids']
        
        pp['record_loc'] = {}
        pp['record_loc']['v'] = ['mysoma']+labels+otherbaseslabels
        #pp['record_loc']['ina'] = ['mysoma']+labels+otherbaseslabels
        #pp['record_loc']['ik'] = ['mysoma']+labels+otherbaseslabels
        #pp['record_loc']['ik'] = labels+otherbaseslabels
        #pp['record_loc']['ica'] = labels+otherbaseslabels
        
        
        vplots_soma = [['mysoma','v','k']]
        """
        iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
        vplots = [['mysoma','v','k']]
        iplots_k = [['mysoma','ik','b']]
        iplots_na = [['mysoma','ina','g']]
        colors = self.get_dendritic_colors(show_levels)
        for i in range(show_levels):
            vplots.append([labels[i],'v',colors[i]])
            iplots_k.append([labels[i],'ik',colors[i]])
            iplots_na.append([labels[i],'ina',colors[i]])
        #print vplots
        if nb>1:
            vplots.append([otherbaseslabels[0],'v','b-'])
        iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
        """
        pp['plot'] = {1:vplots_soma}
        """ #2:vplots,
                      #3:iplots_soma,
                      #4:iplots_k,
                      #5:iplots_na}
        """
        pp['num_threads'] = 1
        if runon:
            #self.es.run_single_experiment(self.expbase, 'cluster', pp)
            self.es.run_single_experiment(self.expbase, 'missing', pp)
            
        else:
            self.es.run_single_experiment(self.expbase, 'local', pp)
            
            
            

            
            
    def run_collection_intrinsic(self,whole=True):
        count = 0
        collection = self.collection_same_total
        Ia = 0.1
        w = 5
        for tree in collection:
            nb,nc,nl = tree
            if nb == 1 and not whole:
                continue
            for factor in self.factors:
                for irr in self.irrs: 
                    
                    #labels,locs = self.get_dend0_child_locations(nl)
                    #clamploc = nl-1
                    pp = {}
                    
                    #pp['stim_iclamp_train'] = True
                    #max_freq = 200
                    #t_interval = 1000./max_freq
                    #pp['iclamp_train'] = [{'tstims': np.arange(0.,self.tstop,t_interval),  'location': labels[-1], 'weights':[w],'amp':Ia, 'dur':1.}]
                    
                    
                    pp.update({'expname':self.expbase,
                               'description':'irr%.2f_factor%.2f_nb%g_ns%g_nl%g_highExp3'%(irr,factor,nb,nc,nl)})
                    print pp['description']
                    self.run_collection_generic_highexp(pp,tree,irr,factor,whole)
                    count += 1
                
                
        print '%g jobs submitted'%count
    
    def run_collection_intrinsic_drive(self,whole=True):
        count = 0
        collection = self.collection_same_total
        Ia = 10.
        w = 1.
        for tree in collection:
            nb,nc,nl = tree
            if nb == 1 and not whole:
                continue
            for factor in self.factors:
                for irr in self.irrs: 
                    
                    labels,locs = self.get_dend0_child_locations(nl-1)
                    clamploc = nl-1
                    pp = {}
                    
                    pp['stim_iclamp_train'] = True
                    max_freq = 200
                    t_interval = 1000./max_freq
                    pp['iclamp_train'] = [{'tstims': np.arange(0.,self.tstop,t_interval),  'location': labels[-1], 'weights':[w],'amp':Ia, 'dur':1.}]
                    
                    
                    pp.update({'expname':self.expbase,
                               'description':'irr%.2f_factor%.2f_nb%g_ns%g_nl%g_illDriven'%(irr,factor,nb,nc,nl)})
                    print pp['description']
                    self.run_collection_generic(pp,tree,irr,factor,whole)#,runon=False)
                    return
                    count += 1
                
                
        print '%g jobs submitted'%count
    
    def run_collection_clamp(self,collection,whole=True):
        
        stim_interval = 50.
        count = 0
        for tree in collection:
            nb,nc,nl = tree
            if nb == 1 and not whole:
                continue
            
            show_levels = np.min([nl,5])
            for factor in self.factors:
                for irr in self.irrs: 
                    for Ia in self.clamp_amps:
                        
                        for clamploc in [nl-1]: #range(show_levels):
                            
                            labels,locs = self.get_dend0_child_locations(nl)
                            pp = {}
                            
                            """
                            # i clamp doesn't work with such abstract models, so we're not going to use this input, and instead
                            # supply a regular spike train.
                            pp['stim_iclamp'] = True
                            pp['iclamp'] = [{'amp':Ia,'tstim':0., 'duration':self.tstop,'location':labels[clamploc]}]
                            
                            pp['stim_spiketrains'] = True
                            pp['spiketrains'] = [{'tstims': [np.arange(0,self.tstop,stim_interval)],  'locations': [labels[clamploc]], 'weights':[Ia],  'el': 0.1}]
                            """
                            pp['stim_iclamp_train'] = True
                            max_freq = 80
                            t_interval = 1000./max_freq
                            pp['iclamp_train'] = [{'tstims': np.arange(0.,self.tstop,t_interval),  'location': labels[clamploc], 'weights':[Ia],'amp':Ia, 'dur':5}]
                            
                            pp.update({'expname':self.expbase,
                                       'description':'irr%.2f_factor%.2f_nb%g_ns%g_nl%g_iclamp%.1f_loc%s'%(irr,factor,nb,nc,nl,Ia,clamploc)})
                            
                            self.run_collection_generic(pp,tree,irr,factor,whole)
                            count += 1
                    
        print '%g jobs submitted'%count


    def get_spiketimes(self,rate,num_inputs):
        spikes = []
        stg = ntst.StGen()
        for i in range(num_inputs):
            ss = stg.poisson_generator(rate,t_stop=self.tstop,array=True)
            spikes.append(ss)
        #print 'Spikes = ',spikes
        return spikes
    
    def add_locations(self,num_loc,upper_id,label='dend0'):
        labels = [label+str(i) for i in range(num_loc)]
        ids = range(upper_id)
        np.random.shuffle(ids)
        locs = ids[:num_loc]
        locs = [(l,0.5) for l in locs]
        return labels,locs
        

    def run_collection_spiketrain(self,whole=True):
        collection = self.collection_same_total
        nsegs = 10
        stim_interval = 50.
        count = 0
        for tree in collection:
            nb,nc,nl = tree
            if nb == 1 and not whole:
                continue
            
            show_levels = np.min([nl,5])
            for factor in self.factors:
                for irr in self.irrs: 
                    for freq in self.freqs:
                        for J in self.Js: 
                            
                            for clamploc in [nl-1]:#0]: #,nl-1]: #
                                ###################################################################
                                
                                
                                pp = {}
                                
                                pp['stim_spiketrains'] = True
                                #locations = ['mysoma','myapic']
                                
                                spbases,splabels,splocs = self.get_spiketrain_locations_nloc(clamploc,nsegs)
                                pp['spiketrains'] = [{'tstims': self.get_spiketimes(freq,nsegs),  'locations': splabels, 'weights':np.ones(nsegs)*J,  'el': 0.02}]
                                pp.update({'expname':self.expbase,
                                           'description':'irr%.2f_factor%.2f_nb%g_ns%g_nl%g_spikes%g_loc%s_J%.1f'%(irr,factor,nb,nc,nl,freq,clamploc,J)})
                                print pp['description']
                                self.run_collection_generic(pp,tree,irr,factor,whole,other_tag_locations=(spbases,splabels,splocs))#,runon=False)
                                #return
                                count += 1
                #return
                    
        print '%g jobs submitted'%count

    
    def run_collection_spiketrain_soma(self,whole=True):
        collection = self.collection_same_total
        nsegs = 10
        count = 0
        for tree in collection:
            nb,nc,nl = tree
            if nb == 1 and not whole:
                continue
            for factor in self.factors:
                for irr in self.irrs: 
                    for freq in self.freqs:
                        for J in self.Js: 
                            
                            ###################################################################
                            
                            
                            pp = {}
                            
                            pp['stim_spiketrains'] = True
                            #locations = ['mysoma','myapic']
                            
                            spbases,splabels,splocs = self.get_spiketrain_locations_nloc_soma(nsegs)
                            pp['spiketrains'] = [{'tstims': self.get_spiketimes(freq,nsegs),  'locations': splabels, 'weights':np.ones(nsegs)*J,  'el': 0.02}]
                            pp.update({'expname':self.expbase,
                                       'description':'irr%.2f_factor%.2f_nb%g_ns%g_nl%g_spikes%s_loc%s_J%.1f'%(irr,factor,nb,nc,nl,freq,'soma',J)})
                            
                            self.run_collection_generic(pp,tree,irr,factor,whole,other_tag_locations=(spbases,splabels,splocs))#,runon=False)
                            #return
                            #print pp['description']
                            count += 1
            #return
                    
        print '%g jobs submitted'%count


    
    def analyse_this(self,etype):
        
        af = run_analysis.AnalyticFrame()
        af.update_params({'tstart':0,'tstop':self.tstop,'label_format':'irr%.2f_factor%.2f_nb%g_ns%g_nl%g_spikes%g_loc%s_J%.1f'})
        
        af.update_params({'tstart':self.light_on,'tstop':self.light_on+self.light_dur,
                                  'tstart_bg': 0,'tstop_bg':self.light_on,
                                  'tstart_post':self.light_on+self.light_dur,'tstop_post':self.tstop})
        
        #irr0.10_factor0.00_nb2_ns2_nl6_spikes30_loc0_J2.0_NpHR_none_ChR_whole
        aa = ['_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2]) for tree in self.collection_same_total] 
        print aa
        if etype == 'proximal':
            exp_comp_list = [['irr0.05_factor%.2f'+'_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'spikes%g'+'_loc%g_J%.1f'%(0,self.Js[0])+'_NpHR_%s_ChR_%s'%(optlog,optlog),'tree=%s, %s'%(tree,optlog)] for tree in self.collection_same_total for optlog in self.explist]
        elif etype == 'soma':
            exp_comp_list = [['irr0.05_factor%.2f'+'_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'spikes%g'+'_loc%s_J%.1f'%('soma',1)+'_NpHR_%s_ChR_%s'%(optlog,optlog),'tree=%s, %s'%(tree,optlog)] for tree in self.collection_same_total for optlog in self.explist]
        elif etype == 'distal':
            exp_comp_list = [['irr0.05_factor%.2f'+'_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'spikes%g'+'_loc%g_J%.1f'%(tree[2]-1,self.Js[0])+'_NpHR_%s_ChR_%s'%(optlog,optlog),'tree=%s, %s'%(tree,optlog)] for tree in self.collection_same_total for optlog in self.explist]
        else:
            print 'Rethink your choices and rethink your life ...'
            return
        
        print exp_comp_list
      
    
        expss = [ec[0] for ec in exp_comp_list]
        explabels = [ec[1] for ec in exp_comp_list]
        print explabels
        
        af.populate_expset(self.expbase,expss,explabels, [self.factors,self.freqs])
    
        af.submenu_extractSpikes()
        
        af.submenu_runFI()
        for exp in af.experimentset:
            exp.calculate_responses('FI')
            exp.calculate_responses('FI_bg')
            exp.calculate_responses('FI_post')

        af.submenu_save()
        
            
    def plot_gain(self):
        self.factors.sort()
        
        
        for exp in self.explist:
            
            for tree in self.collection_same_total:
            
                af = run_analysis.AnalyticFrame()
                af.update_params({'tstart':self.light_on,'tstop':self.light_on+self.light_dur})
                
                
                
                exp_comp_list = [['irr0.05'+'_factor%.2f'%(factor)+'_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'spikes%g'+'_loc%g_J%.1f'%(tree[2]-1,self.Js[0])+'_NpHR_%s_ChR_%s'%(exp,exp),'%g'%factor] for factor in self.factors]
                print exp_comp_list
        
                expss = [ec[0] for ec in exp_comp_list]
                explabels = [ec[1] for ec in exp_comp_list]
                af.populate_expset(self.expbase,expss,explabels,[self.freqs])
        
        
            
                af.submenu_load()
                af.submenu_print()
                af.submenu_plot(5, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp,exp))
                af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp,exp))
        

    def analyse_gain_background(self,etype):
        
        af = run_analysis.AnalyticFrame()
        
        
        af.update_params({'tstart':0,'tstop':self.tstop,'label_format':'irr%.2f_factor%.2f_nb%g_ns%g_nl%g_spikes%g_loc%s_J%.1f'})
        """
        af.update_params({'tstart':self.light_on,'tstop':self.light_on+self.light_dur,
                                  'tstart_bg': 0,'tstop_bg':self.light_on,
                                  'tstart_post':self.light_on+self.light_dur,'tstop_post':self.tstop})
        optlocations = ['whole']
        #irr0.10_factor0.00_nb2_ns2_nl6_spikes30_loc0_J2.0_NpHR_none_ChR_whole
        aa = ['_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2]) for tree in self.collection_same_total] 
        print aa
        exp_comp_list = [['irr0.05_factor%.2f'+'_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'spikes%g_loc0_J2.0'+'_NpHR_%s_ChR_%s'%(optlog,optlog),'tree=%s, %s'%(tree,optlog)] for tree in self.collection_same_total for optlog in optlocations]
        print exp_comp_list
    
        expss = [ec[0] for ec in exp_comp_list]
        explabels = [ec[1] for ec in exp_comp_list]
        print explabels
        
        af.populate_expset(self.expbase,expss,explabels, [self.factors,self.freqs])
    
        af.submenu_extractSpikes()
        af.submenu_runFI()
        for exp in af.experimentset:
            exp.calculate_responses('FI_bg')
            exp.calculate_responses('FI_post')
        af.submenu_save()
        
        """
        self.factors.sort()
        for exp in self.explist:
            
            for tree in self.collection_same_total:
            
                af = run_analysis.AnalyticFrame()
                af.update_params({'tstart':self.light_on,'tstop':self.light_on+self.light_dur,
                                  'tstart_bg': 0,'tstop_bg':self.light_on,
                                  'tstart_post':self.light_on+self.light_dur,'tstop_post':self.tstop})
                
                if etype =='proximal':
                    descript = 'proximal'
                    exp_comp_list = [['irr0.05'+'_factor%.2f'%(factor)+'_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'spikes%g'+'_loc%g_J%.1f'%(0,self.Js[0])+'_NpHR_%s_ChR_%s'%(exp,exp),'%g'%factor] for factor in self.factors]
                elif etype =='soma':
                    descript = 'soma'
                    exp_comp_list = [['irr0.05'+'_factor%.2f'%(factor)+'_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'spikes%g'+'_loc%s_J%.1f'%('soma',self.Js[0])+'_NpHR_%s_ChR_%s'%(exp,exp),'%g'%factor] for factor in self.factors]
                elif etype =='distal':
                    descript = 'distal'
                    exp_comp_list = [['irr0.05'+'_factor%.2f'%(factor)+'_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'spikes%g'+'_loc%g_J%.1f'%(tree[2]-1,self.Js[0])+'_NpHR_%s_ChR_%s'%(exp,exp),'%g'%factor] for factor in self.factors]
                print exp_comp_list
                
                expss = [ec[0] for ec in exp_comp_list]
                explabels = [ec[1] for ec in exp_comp_list]
                af.populate_expset(self.expbase,expss,explabels,[self.freqs])
                af.submenu_extractSpikes()
                """
                
                af.submenu_runFI()
                for exp in af.experimentset:
                    exp.calculate_responses('FI')
                    exp.calculate_responses('FI_bg')
                    exp.calculate_responses('FI_post')
                af.submenu_save()

                """
                af.submenu_load()
                #af.submenu_print()
                af.submenu_plot(5, self.expbase+'_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'_inj%s'%(descript))
                #af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
                af.submenu_plot(10, self.expbase+'bg_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'inj%s'%(descript))



    def analyse_vitro(self,etype):
        
        af = run_analysis.AnalyticFrame()
        af.update_params({'tstart':0,'tstop':self.tstop,
                          'label_format':'irr%.2f_factor%.2f_nb%g_ns%g_nl%g_iclamp%.1f_loc%s'})
        
        af.update_params({'tstart':self.light_on,'tstop':self.light_on+self.light_dur,
                                  'tstart_bg': 0,'tstop_bg':self.light_on,
                                  'tstart_post':self.light_on+self.light_dur,'tstop_post':self.tstop})
        
        #irr0.10_factor0.00_nb2_ns2_nl6_spikes30_loc0_J2.0_NpHR_none_ChR_whole
        aa = ['_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2]) for tree in self.collection_same_total] 
        print aa
        if etype == 'proximal':
            exp_comp_list = [['irr0.05_factor%.2f'+'_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'iclamp%.1f'+'_loc%g'%(0)+'_NpHR_%s_ChR_%s'%(optlog,optlog),'tree=%s, %s'%(tree,optlog)] for tree in self.collection_same_total for optlog in self.explist]
        elif etype == 'soma':
            exp_comp_list = [['irr0.05_factor%.2f'+'_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'iclamp%.1f'+'_loc%s'%('soma')+'_NpHR_%s_ChR_%s'%(optlog,optlog),'tree=%s, %s'%(tree,optlog)] for tree in self.collection_same_total for optlog in self.explist]
        elif etype == 'distal':
            exp_comp_list = [['irr0.05_factor%.2f'+'_nb%g_ns%g_nl%g_'%(tree[0],tree[1],tree[2])+'iclamp%.1f'+'_loc%g'%(tree[2]-1)+'_NpHR_%s_ChR_%s'%(optlog,optlog),'tree=%s, %s'%(tree,optlog)] for tree in self.collection_same_total for optlog in self.explist]
        else:
            print 'Rethink your choices and rethink your life ...'
            return
        
        print exp_comp_list
      
    
        expss = [ec[0] for ec in exp_comp_list]
        explabels = [ec[1] for ec in exp_comp_list]
        print explabels
        
        af.populate_expset(self.expbase,expss,explabels, [self.factors,self.freqs])
    
        af.submenu_extractSpikes()
        
        af.submenu_runFI()
        for exp in af.experimentset:
            exp.calculate_responses('FI')
            exp.calculate_responses('FI_bg')
            exp.calculate_responses('FI_post')

        af.submenu_save()
        

                            
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print 'Need to specify an input mode'
    elif sys.argv[1] == 'soma':
        wholecell,inputtype = sys.argv[2],sys.argv[3]
        Experiment_140305(wholecell,inputtype).run_collection_spiketrain_soma()
    elif sys.argv[1] == 'analyse':
        wholecell,inputtype = sys.argv[2],sys.argv[3]
        Experiment_140305(wholecell,inputtype).analyse_this(inputtype)
    elif sys.argv[1] == 'gain':
        wholecell,inputtype = sys.argv[2],sys.argv[3]
        Experiment_140305(wholecell,inputtype).plot_gain()
    elif sys.argv[1] == 'bg':
        wholecell,inputtype = sys.argv[2],sys.argv[3]
        Experiment_140305(wholecell,inputtype).analyse_gain_background(inputtype)        
    else:
        wholecell,inputtype = sys.argv[2],sys.argv[3]
        print wholecell, inputtype
        #try:
        getattr(Experiment_140305(wholecell,inputtype),'run_collection_%s'%sys.argv[1])()
        #except:
        #    print 'Could not execute for ','Exp(%s,%s).run_collection_%s()'%(sys.argv[2],sys.argv[3],sys.argv[1])
    





