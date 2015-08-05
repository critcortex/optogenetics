import sys
import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import run_stimulation
import pylab
import file_io as fio
import run_analysis
import NeuroTools.stgen as ntst 

"""

Log:

14/07/23: Ran in vitro version, with stim0 being drawn from [100,-1]. This 
          meant that the current injection site was chosen randomly from anywhere on the 
          dendritic tree.
          During analysis, I noticed massive variability in the background rate (unsurprisingly).
          So I updated to use stim1 which was drawn from a distance of [100,150]. 
          -- Update:
          Analysed this latest data, still too variable. 
          Now using dist_vitro = [100.3,100.4] --> selects only one idx
          This set will be known as stim2


"""


expbase = '140408_scan_SHstell'


stg = ntst.StGen()
es = ExpSetup()
 
 
tstop = 2700
freqs = [0.,50.,80.]
freqs = range(10,201,10)

light_on = 700
light_dur = 1000
#factors = [0]+[ 0.0625, 0.125,0.25,0.375,0.5,0.75,1.,2.,4.,8.]
factors = [0]+[ 0.0625, 0.125,0.25,0.375,0.5,0.75,1.,1.5,2.]
factors = [ 0.0625, 0.125,0.25,0.375,0.5,0.75,1.,1.5,2.]
irrs = np.linspace(0.0025, 0.05, 11)
irrs = [0.003,0.007,0.012]
irrs = [0.007,0.012]

iclamp_amps = np.arange(-1.0,1.7,0.2)
iclamp_amps = np.arange(1.8,3.1,0.2)
iclamp_amps = np.arange(-1.,3.1,0.2)
#iclamp_amps = np.arange(5.,10.1,1.)

#Js = np.arange(1.,5.)
Js = [2]

#nsite_range = range(20,51,20)
nsite_range = [40]
dist = [100,-1]

dist_vitro = [100,150]
dist_vitro = [100.3,100.4]
dist_vitro = [304,305]
#irrs = [0]
#factors = [0]
#iclamp_amps = [0]

# MODS: 150622 - 
irrs = [0.007]
#iclamp_amps = np.arange(-1.,3.1,0.1)
#iclamp_amps = np.arange(3.1,5.1,0.1)
iclamp_amps = np.arange(-1.0,5.1,0.1)
#iclamp_amps = [2.9]
dist_vitro = [100.3,100.4]

# /MODS

def get_spiketimes(rate,num_inputs):
    spikes = []
    for i in range(num_inputs):
        ss = stg.poisson_generator(rate,t_stop=tstop,array=True)
        spikes.append(ss)
    #print 'Spikes = ',spikes
    return spikes

def scan_locations_optogen_invivo():
    #iclamp_amps = [0]
    count = 0
    for freq in freqs:
        for factor in factors:
            for irr in irrs: 
                for (i,J) in enumerate(Js):
                    for nsites in nsite_range:
                    
                        if freq == 0 and i > 0:
                            continue
                        
                        pp = {}
                        # neuron model and params
                        pp['cell'] = ['Neuron', 'SHStellate']
                        pp['cell_params'] = {}
                        
                        # opsin
                        chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                        hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                        
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
                        pp['tstop'] = tstop
                        
                        
                        labels = ['stim%g'%i for i in range(nsites)]
                        
                        
                        pp['mark_loc'] = {}
                        pp['mark_loc']['names'] = ['mysoma']+labels
                        pp['mark_loc']['sections'] = ['soma']+['dendrite']*nsites
                        pp['mark_loc']['ids'] = [(0,0.5)] + [('select_section_posn_bydistance',{'sectionarea':'dendrite','mindist':dist[0],'maxdist':dist[1]}) for i in range(nsites)] 
                        
                        
                        pp['record_loc'] = {}
                        pp['record_loc']['v'] = ['mysoma'] #+labels
                        pp['record_loc']['ina'] = ['mysoma']
                        pp['record_loc']['ik'] = ['mysoma']
                        
                        pp['stim_spiketrains'] = True
                        pp['spiketrains'] = [{'tstims': get_spiketimes(freq,nsites), 'locations': labels, 'weights':np.ones(nsites)*J,  'el': 0.1}]
                            
                        
                        
                        vplots_soma = [['mysoma','v','k']]
            
                        #iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
                        pp['plot'] = {1:vplots_soma} 
                            
                        
                        pp['num_threads'] = 1
                                                   
                        pp.update({'expname':expbase,
                                   'description':'irr%.3f_factor%.2f_freq%g_J%g'%(irr,factor,freq,J)})
            
                        es.run_single_experiment(expbase, 'missing', pp)
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

def scan_locations_optogen_invitro(loc):
    #iclamp_amps = [0]
    count = 0

    for factor in factors:
        for irr in irrs: 
            for Ia in iclamp_amps:
                
                pp = {}
                # neuron model and params
                pp['cell'] = ['Neuron', 'SHStellate']
                pp['cell_params'] = {}
                
                # opsin
                chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                
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
                pp['tstop'] = tstop
                
                nsites = 4
                labels = ['stim%g'%i for i in range(nsites)]
                
                
                pp['mark_loc'] = {}
                pp['mark_loc']['names'] = ['mysoma']+labels
                pp['mark_loc']['sections'] = ['soma']+['dendrite']*nsites
                pp['mark_loc']['ids'] = [(0,0.5)] + [('select_section_posn_bydistance',{'sectionarea':'dendrite','mindist':dist_vitro[0],'maxdist':dist_vitro[1]}) for i in range(nsites)] 
                
                
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
                    stimloc = loc
                pp['iclamp'] = [{'tstim':0,  'location': stimloc, 'amp':Ia, 'duration':tstop}]
                       
                
                
                vplots_soma = [['mysoma','v','k']]
    
                #iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
                pp['plot'] = {1:vplots_soma} 
                    
                
                pp['num_threads'] = 1
                                           
                pp.update({'expname':expbase,
                           'description':'irr%.3f_factor%.2f_I%.2f_stimloc_%s'%(irr,factor,Ia,stimloc)})
    
                es.run_single_experiment(expbase, 'missing', pp)
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
    print irrs


def scan_locations_optogen_illumination():
    """
    
    Added 150805 - for illumination driven stellate cell response
    """
    count = 0
    
    irrs = [ a*b for a in [0.01, 0.1] for b in [1,2,5]] + list(np.arange(1.,5.1,0.25)) 
    

    for factor in factors:
        for irr in irrs: 
            pp = {}
            # neuron model and params
            pp['cell'] = ['Neuron', 'SHStellate']
            pp['cell_params'] = {}
            
            # opsin
            chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
            hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
            
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
            pp['tstop'] = tstop
            
            pp['mark_loc'] = {}
            pp['mark_loc']['names'] = ['mysoma']
            pp['mark_loc']['sections'] = ['soma']
            pp['mark_loc']['ids'] = [(0,0.5)]
            
            pp['record_loc'] = {}
            pp['record_loc']['v'] = ['mysoma']
            
            pp['num_threads'] = 1
                                           
            pp.update({'expname':expbase,
                       'description':'irr%.3f_factor%.2f_illOnly'%(irr,factor)})
            
            #es.run_single_experiment(expbase, 'local', pp)
            #return
            count += 1
            
    print '%g jobs submitted'%count
            


def get_optdescript(irr,factor):
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
                                
def analyse_locations_optogen_invivo():
    
    
    for irr in irrs:
        af = run_analysis.AnalyticFrame()
        af.update_params({'tstart':0,'tstop':tstop,'label_format':'irr%.2f_factor%.2f_nb%g_ns%g_nl%g_spikes%g_loc%s_J%.1f'})
    
        af.update_params({'tstart':light_on,'tstop':light_on+light_dur,
                              'tstart_bg': 0,'tstop_bg':light_on,
                              'tstart_post':light_on+light_dur,'tstop_post':tstop})
        optlog = get_optdescript(irr,1)
        exp_comp_list = [['irr%.3f_'%irr+'factor%.2f_'%f+'freq%g'+'_J%g'%(2.)+'_NpHR_%s_ChR_%s'%(optlog[1],optlog[0]),'=%.3f'%f] for f in factors]
    
        print exp_comp_list
  
        expss = [ec[0] for ec in exp_comp_list]
        explabels = [ec[1] for ec in exp_comp_list]
        print explabels
        
        af.populate_expset(expbase,expss,explabels, [freqs])
        
        
        af.submenu_extractSpikes()
        
        af.submenu_runFI()
        for exp in af.experimentset:
            exp.calculate_responses('FI')
            exp.calculate_responses('FI_bg')
            exp.calculate_responses('FI_post')
    
        af.submenu_save()
        af.submenu_print()
        #af.submenu_print()
        #af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
        af.submenu_plot(5, expbase+'FI_gain_invivo_varyFactor_irr%.3f'%irr)
        af.submenu_plot(0, expbase+'FI_gain_invivo_varyFactor_irr%.3f'%irr)
        af.submenu_plot(10, expbase+'FI_gain_invivo_varyFactor_irr%.3f'%irr)
        
    """
    #af.submenu_print()
    #af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
    af.submenu_plot(5, expbase+'FI_gain_varyFactor')
    af.submenu_plot(0, expbase)
    af.submenu_plot(10, expbase)
    """
    
    
                                
def analyse_locations_optogen_invitro(stimloc):
    
    for irr in irrs:
        af = run_analysis.AnalyticFrame()
        af.update_params({'tstart':0,'tstop':tstop,'label_format':'irr%.2f_factor%.2f_nb%g_ns%g_nl%g_spikes%g_loc%s_J%.1f'})
    
        af.update_params({'tstart':light_on,'tstop':light_on+light_dur,
                              'tstart_bg': 0,'tstop_bg':light_on,
                              'tstart_post':light_on+light_dur,'tstop_post':tstop})
        optlog = get_optdescript(irr,1)
        exp_comp_list = [['irr%.3f_'%irr+'factor%.2f_'%f+'I%.2f'+'_stimloc_%s'%(stimloc)+'_NpHR_%s_ChR_%s'%(optlog[1],optlog[0]),'=%.3f'%f] for f in factors]
        #irr%.3f_factor%.2f_I%.2f_stimloc_%s'
        print exp_comp_list
  
        expss = [ec[0] for ec in exp_comp_list]
        explabels = [ec[1] for ec in exp_comp_list]
        print explabels
        
        af.populate_expset(expbase,expss,explabels, [iclamp_amps])
        
        
        af.submenu_extractSpikes()
        
        af.submenu_runFI()
        for exp in af.experimentset:
            exp.calculate_responses('FI')
            exp.calculate_responses('FI_bg')
            exp.calculate_responses('FI_post')
       
        af.submenu_save()
        
        af.submenu_print()
        
        
        af.submenu_plot(5, expbase+'FI_gain_invitro_varyFactor_irr%.3f_stimloc%s'%(irr,stimloc))
        af.submenu_plot(0, expbase+'FI_gain_invitro_varyFactor_irr%.3f_stimloc%s'%(irr,stimloc))
        af.submenu_plot(10, expbase+'FI_gain_invitro_varyFactor_irr%.3f_stimloc%s'%(irr,stimloc))
        
    """
    #af.submenu_print()
    #af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
    af.submenu_plot(5, expbase+'FI_gain_varyFactor')
    af.submenu_plot(0, expbase)
    af.submenu_plot(10, expbase)
    """

""" Added the following three methods from 140612_replay_input
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
    
"""
def rename():
    import shutil
    # hack fnction
    loc = "experiments/140408_scan_SHstell/dat/"
    for i in np.arange(3.1,5.01,0.1):
        for f in factors:
            fn = "140408_scan_SHstellirr0.007_factor%.2f_I%.2f_stimloc_stim%s_NpHR_whole_ChR_whole_v.dat"%(f,i,'%g')
            print loc+fn%3
            shutil.move(loc+fn%3, loc+fn%2)


    
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print 'Need to specify an input mode'
    elif sys.argv[1] == 'scan_vivo':
        scan_locations_optogen_invivo()
    elif sys.argv[1] == 'scan_vitro':
        scan_locations_optogen_invitro(sys.argv[2])        
    elif sys.argv[1] == 'scan_illum':
        scan_locations_optogen_illumination()        
    elif sys.argv[1] == 'analyse':
        analyse_locations_optogen_invivo()
    elif sys.argv[1] == 'analyse_vitro':
        analyse_locations_optogen_invitro(sys.argv[2])
    else:
        print 'Unknown option'
                         
