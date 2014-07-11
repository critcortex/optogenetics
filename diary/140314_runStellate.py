import sys
import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import run_stimulation
import pylab
import file_io as fio
import run_analysis

class Experiment_140314:
    
    
    
    def __init__(self,whole,stimulus):
        self.es = ExpSetup()
        self.stg = ntst.StGen()
        self.whole = (whole == 'whole')
        print 'whole = ', self.whole
        if not self.whole:
            self.explist = ['partialSame']
        else:
            self.explist = ['whole']
            
        self.stimulus = stimulus
        self.expbase = '140313_testStellate'
    
        self.stim_location_fracts = [['peripheral',[0.5,1.0]],
                 ['distributed',[0,1.0]],
                 ['proximal',[0,.5]]]
    
        self.tstop = 3500
    
        self.freq = 50.
        self.dist_frac_sections = 50.

           
        self.light_on = 1000
        self.light_dur = 1500
        
        self.factors = [0,0.125,0.25,1.,4.,8.]
        #self.factors = [0.125,0.25,1.,4.,8.]
        
        self.clamp_amps = np.arange(0.2,5.,0.2)
        self.clamp_amps = [1.5]
        
        
        
        #self.factors = [0.5,2,1.5,1.,0.75,1.25]
        #self.factors = [0,0.5,2,1.5,1.,0.75,1.25]
        #self.factors = [10.]
        self.irrs = [0.05]#,0.02]
        self.irrs = [0.0005,0.005]#,0.02]
        self.irrs = [0.001]#,0.02]
        self.freqs = range(20,81,10)#101,10)
        
        self.Js = [30.] # for soma simulation
        
        
        #intrinsic
        
        self.irrs = np.arange(0.001,0.0101,0.001)
        self.irrs = np.arange(0.10,1.001,0.1)
        #self.irrs = [0.01]
        self.factors = [0,0.125,0.25,1.,4.]
        
        

    
    def run_collection_generic(self,pp,irr,factor,whole,other_tag_locations=None,runon=True):
        
        # neuron model and params
        pp['cell'] = ['Neuron', 'Stellate']
        pp['cell_params'] = {}
        # opsin
        
        #print locs
        chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
        hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': self.light_dur,'lightdelay':self.light_on,'interpulse_interval':250,  'n_pulses':1}
        
        
        pp['opsindict'] = {}
        areas = ['soma', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'hill', 'node', 'iseg', 'myelin']
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
            areas = ['soma', 'a1', 'a2', 'a3', 'hill', 'node', 'iseg', 'myelin']
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
        
        otherbaseslabels,othersections = [],[]
        
        if other_tag_locations is not None:
            obases,olabels,olocs = other_tag_locations
            #print olocs
        else:
            obases,olabels,olocs = [],[],[]
            
        #print olocs
        pp['mark_loc'] = {}
        pp['mark_loc']['names'] = ['mysoma']+otherbaseslabels+olabels
        pp['mark_loc']['sections'] = ['soma']+othersections+obases
        pp['mark_loc']['ids'] = [(0,0.5)] + [('get_section_byname',{'sectioname':loc,'subdomain':'dend1'}) for loc in otherbaseslabels]+[(0,loc[1]) for loc in olocs ]
        #print pp['mark_loc']['ids']
        
        pp['record_loc'] = {}
        pp['record_loc']['v'] = ['mysoma']+otherbaseslabels
        pp['record_loc']['ina'] = ['mysoma']+otherbaseslabels
        pp['record_loc']['ik'] = ['mysoma']+otherbaseslabels
        #pp['record_loc']['ik'] = labels+otherbaseslabels
        #pp['record_loc']['ica'] = labels+otherbaseslabels
        
        
        vplots_soma = [['mysoma','v','k']]
        iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
        vplots = [['mysoma','v','k']]
        iplots_k = [['mysoma','ik','b']]
        iplots_na = [['mysoma','ina','g']]
        
        iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
        pp['plot'] = {#1:vplots,
                      2:vplots_soma}#,
                      #3:iplots_soma,
                      #4:iplots_k,
                      #5:iplots_na}
        
        pp['num_threads'] = 1
        if runon:
            #self.es.run_single_experiment(self.expbase, 'cluster', pp)
            self.es.run_single_experiment(self.expbase, 'missing', pp)
            
        else:
            self.es.run_single_experiment(self.expbase, 'local', pp)
            
    
    def get_spiketimes(self,rate,num_inputs):
        spikes = []
        
        for i in range(num_inputs):
            ss = self.stg.poisson_generator(rate,t_stop=self.tstop,array=True)
            spikes.append(ss)
        #print 'Spikes = ',spikes
        return spikes
    
        
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
                        
                        for clamploc in range(show_levels):
                            
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
                                       'description':'irr%.3f_factor%.2f_nb%g_ns%g_nl%g_iclamp%.1f_loc%s'%(irr,factor,nb,nc,nl,Ia,clamploc)})
                            
                            self.run_collection_generic(pp,tree,irr,factor,whole)
                            count += 1
                    
        print '%g jobs submitted'%count

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
    

    def add_locations(self,num_loc,upper_id,label='a1'):
        labels = [label+str(i) for i in range(num_loc)]
        ids = range(upper_id)
        np.random.shuffle(ids)
        locs = ids[:num_loc]
        locs = [(l,0.5) for l in locs]
        return labels,locs
        
        
    def run_collection_intrinsic(self,whole=True):
        count = 0
        
        for factor in self.factors:
            for irr in self.irrs: 
                
                
                pp = {}
                
                pp.update({'expname':self.expbase,
                           'description':'irr%.3f_factor%.2f_'%(irr,factor)})
                
                self.run_collection_generic(pp,irr,factor,whole)
                count += 1
                
            
            
        print '%g jobs submitted'%count

    def run_collection_spiketrain(self,collection,whole=True):
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
                            clamploc = 0 ##################################
                            pp = {}
                            
                            pp['stim_spiketrains'] = True
                            #locations = ['mysoma','myapic']
                            
                            spbases,splabels,splocs = self.get_spiketrain_locations_nloc(clamploc,nsegs)
                            pp['spiketrains'] = [{'tstims': self.get_spiketimes(freq,nsegs),  'locations': splabels, 'weights':np.ones(nsegs)*J,  'el': 0.02}]
                            pp.update({'expname':self.expbase,
                                       'description':'irr%.32f_factor%.2f_nb%g_ns%g_nl%g_spikes%g_loc%s_J%.1f'%(irr,factor,nb,nc,nl,freq,clamploc,J)})
                            
                            self.run_collection_generic(pp,tree,irr,factor,whole,other_tag_locations=(spbases,splabels,splocs))#,runon=False)
                            #return
                            count += 1
            #return
                    
        print '%g jobs submitted'%count


    
    def run_collection_spiketrain_soma(self,whole=True):
        
        nsegs = 10
        stim_interval = 50.
        count = 0
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
                                   'description':'irr%.3f_factor%.2f_spikes%s_loc%s_J%.1f'%(irr,factor,freq,'soma',J)})
                        
                        self.run_collection_generic(pp,irr,factor,whole,other_tag_locations=(spbases,splabels,splocs))#,runon=False)
                        #return
                        count += 1
            #return
                    
        print '%g jobs submitted'%count





    
    def analyse_this(self,etype):
        
        af = run_analysis.AnalyticFrame()
        af.update_params({'tstart':0,'tstop':self.tstop,'label_format':'irr%.2f_factor%.2f_spikes%s_loc%s_J%.1f'})
        
        af.update_params({'tstart':self.light_on,'tstop':self.light_on+self.light_dur,
                                  'tstart_bg': 0,'tstop_bg':self.light_on,
                                  'tstart_post':self.light_on+self.light_dur,'tstop_post':self.tstop})
        
        #irr0.10_factor0.00_nb2_ns2_nl6_spikes30_loc0_J2.0_NpHR_none_ChR_whole
        
        if etype == 'proximal':
            exp_comp_list = [['irr%.3f'%self.irrs[0]+'_factor%.2f'+'_spikes%g'+'_loc%g_J%.1f'%(0,self.Js[0])+'_NpHR_%s_ChR_%s'%(optlog,optlog),'%s'%(optlog)]  for optlog in self.explist]
        elif etype == 'soma':
            exp_comp_list = [['irr%.3f'%self.irrs[0]+'_factor%.2f'+'_spikes%g'+'_loc%s_J%.1f'%('soma',self.Js[0])+'_NpHR_%s_ChR_%s'%(optlog,optlog),' %s'%(optlog)] for optlog in self.explist]
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
        
        
        af.update_params({'label_format':'irr%.2f_factor%.2f_spikes%g_loc%s_J%.1f'})
        self.factors.sort()
        for exp in self.explist:
            
        
            af.update_params({'tstart':self.light_on,'tstop':self.light_on+self.light_dur,
                              'tstart_bg': 0,'tstop_bg':self.light_on,
                              'tstart_post':self.light_on+self.light_dur,'tstop_post':self.tstop})
            
            if etype =='proximal':
                descript = 'proximal'
                exp_comp_list = [['irr0.001'+'_factor%.2f'%(factor)+'_'+'spikes%g'+'_loc%g_J%.1f'%(0,self.Js[0])+'_NpHR_%s_ChR_%s'%(exp,exp),'%g'%factor] for factor in self.factors]
            elif etype =='soma':
                descript = 'somal'
                exp_comp_list = [['irr0.001'+'_factor%.2f'%(factor)+'_'+'spikes%g'+'_loc%s_J%.1f'%('soma',self.Js[0])+'_NpHR_%s_ChR_%s'%(exp,exp),'%g'%factor] for factor in self.factors]
            print exp_comp_list
    
            expss = [ec[0] for ec in exp_comp_list]
            explabels = [ec[1] for ec in exp_comp_list]
            af.populate_expset(self.expbase,expss,explabels,[self.freqs])
            
    
            
            af.submenu_load()
            af.submenu_print()
            #af.submenu_plot(5, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
            #af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
            af.submenu_plot(10, self.expbase+'FI_gain_bg_irr%.3f_stellate_%s_varyFactor_exp%s%s_'%(0.001,descript,exp,exp))




























                            
if __name__ == '__main__':
    print sys.argv
    if len(sys.argv) <= 1:
        print 'Need to specify an input mode'
    elif sys.argv[1] == 'soma':
        print 'A'
        wholecell,inputtype = sys.argv[2],sys.argv[3]
        Experiment_140314(wholecell,inputtype).run_collection_spiketrain_soma()
    elif sys.argv[1] == 'analyse':
        print 'Analyse'
        wholecell,inputtype = sys.argv[2],sys.argv[3]
        Experiment_140314(wholecell,inputtype).analyse_this(inputtype)
    elif sys.argv[1] == 'bg':
        print 'BG'
        wholecell,inputtype = sys.argv[2],sys.argv[3]
        Experiment_140314(wholecell,inputtype).analyse_gain_background(inputtype) 
    else:
        print 'B'
        wholecell,inputtype = sys.argv[2],sys.argv[3]
        print wholecell, inputtype
        #try:
        getattr(Experiment_140314(wholecell,inputtype),'run_collection_%s'%sys.argv[1])()
        #except:
        #    print 'Could not execute for ','Exp(%s,%s).run_collection_%s()'%(sys.argv[2],sys.argv[3],sys.argv[1])
    
    """
    
    elif sys.argv[1] == 'gain':
        wholecell,inputtype = sys.argv[2],sys.argv[3]
        Experiment_140314(wholecell,inputtype).plot_gain()
    
    """   




            