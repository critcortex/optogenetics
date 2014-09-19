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



expbase = '140723_L5PC_basal_irr_ChRonly'


stg = ntst.StGen()
es = ExpSetup()
 
light_on = 1050
light_dur = 1000
tstop = light_on + light_dur + 50 # as buffer

# Optogenetics parameters
factors = [0,1.]
irrs = [ coeff*power for power in [0.001,0.01,0.1,1.] for coeff in [1,2,5] ]

# In vitro parameters
iclamp_amps = np.arange(0.0,5.51,0.5)

# In vivo parameters
freqs = range(20,150,10)
freqs = range(15,150,10) + [150]
freqs = range(2,15,2)+range(15,151,5)
Js = [2.]
nsite_range = [80]
dist = [100,-1]




L5PC_areas = ['soma', 'axon', 'apic','dend']
#L5PC_areas = ['soma','dend']
SHStellate_areas = ['soma','dendrite','axon']

TOTAL_NUM_AREAS = {'L5PC': 4,
                   'SHStellate':2 }
areas = {'L5PC': L5PC_areas,
         'SHStellate':SHStellate_areas}

celltype = 'L5PC'


#######################
factors = [0.125,0.25,0.5,1.,2.,4.,8.]
factors = [0.125,0.25,0.5,0.75,1.,1.5,2.] # for vitro
factors = [0.001,0.125,0.25,0.5,0.75,1.] # for vivo
#factors = [0.001]



irrs = [0.01,0.02]
irrs = [0.001,0.002]


def get_spiketimes(rate,num_inputs):
    spikes = []
    for i in range(num_inputs):
        # TODO: am guessing there's a more efficient way of implementing this
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
                        pp['cell'] = ['Neuron',celltype]
                        pp['cell_params'] = {}
                        
                        # opsin
                        chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                        hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                        
                        pp['opsindict'] = {}
                        if irr > 0 : 
                            pp['opsindict']['ChR'] =  {}
                            for area in areas[celltype]:
                                pp['opsindict']['ChR'][area] = chrdict
                            if len(areas[celltype])==TOTAL_NUM_AREAS[celltype]:
                                pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
                            else:
                                pp['ChR_areas'] = {'partial'      : pp['opsindict']['ChR'].keys()}
                        else:
                            pp['ChR_areas'] = {'none'      : [None]}
                            
                        if irr > 0 and factor > 0 : 
                            pp['opsindict']['NpHR'] =  {}
                            for area in areas[celltype]:
                                pp['opsindict']['NpHR'][area] = hrdict
                            if len(areas[celltype])==TOTAL_NUM_AREAS[celltype]:
                                pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
                            else:
                                pp['NpHR_areas'] = {'partial'      : pp['opsindict']['NpHR'].keys()}
                            
                        else:
                            pp['NpHR_areas'] = {'none'      : [None]}
                            
                        # general settings 
                        pp['experiment_type'] = 'opsinonly'
                        pp['savedata'] = True # False #True
                        
                        pp['tstart'] = 0
                        pp['tstop'] = tstop
                        
                        
                        labels = ['stim%g'%i for i in range(nsites*2)]
                        
                        
                        pp['mark_loc'] = {}
                        pp['mark_loc']['names'] = ['mysoma']+labels
                        pp['mark_loc']['sections'] = ['soma']+['apic']*nsites+['dend']*nsites
                        pp['mark_loc']['ids'] = [(0,0.5)] + [('select_section_posn_bydistance',{'sectionarea':'apic','mindist':dist[0],'maxdist':dist[1]}) for i in range(nsites)] + [('select_section_posn_bydistance',{'sectionarea':'dend','mindist':dist[0],'maxdist':dist[1]}) for i in range(nsites)] 
                        
                        
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
                                   'description':'irr%.3f_factor%.2f_freq%g_J%g_nsites%g'%(irr,factor,freq,J,nsites)})
            
                        es.run_single_experiment(expbase, 'missing', pp)
                        #es.run_single_experiment(expbase, 'names', pp)
                        count += 1
                        #es.run_single_experiment(expbase, 'local', pp)
                        #return 
                        
            #return 
                        
    print '%g jobs submitted'%count


def scan_locations_optogen_invivo_same_total():
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
                        pp['cell'] = ['Neuron', 'L5PC']
                        pp['cell_params'] = {}
                        
                        # opsin
                        chr_irr = factor/(1+factor)*irr
                        hr_irr = 1./(1+factor)*irr
                        chrdict =  {'exp':5e-4, 'irradiance' :chr_irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                        hrdict =  {'exp':5e-4, 'irradiance' :hr_irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                        
                        pp['opsindict'] = {}
                        if irr > 0 :
                            pp['opsindict']['ChR'] =  {}
                            for area in L5PC_areas:
                                pp['opsindict']['ChR'][area] = chrdict
                            pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
                        else:
                            pp['ChR_areas'] = {'none'      : [None]}
                            
                        if irr > 0 and factor > 0:
                            pp['opsindict']['NpHR'] =  {}
                            for area in L5PC_areas:
                                pp['opsindict']['NpHR'][area] = hrdict
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
                        pp['mark_loc']['sections'] = ['soma']+['apic']*nsites
                        pp['mark_loc']['ids'] = [(0,0.5)] + [('select_section_posn_bydistance',{'sectionarea':'apic','mindist':dist[0],'maxdist':dist[1]}) for i in range(nsites)] 
                        
                        
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
                                   'description':'constant_irr%.3f_factor%.2f_freq%g_J%g_nsites%g'%(irr,factor,freq,J,nsites)})
            
                        #es.run_single_experiment(expbase, 'missing', pp)
                        count += 1
                        es.run_single_experiment(expbase, 'names', pp)
                        #es.run_single_experiment(expbase, 'local', pp)
            #return 
                        
    print '%g jobs submitted'%count



def scan_locations_optogen_invitro(stimloc):
    #iclamp_amps = [0]
    count = 0

    for factor in factors:
        for irr in irrs: 
            for Ia in iclamp_amps:
                
                
                pp = {}
                # neuron model and params
                pp['cell'] = ['Neuron',celltype]
                pp['cell_params'] = {}
                
                # opsin
                chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                
                pp['opsindict'] = {}
                if irr > 0 : 
                    pp['opsindict']['ChR'] =  {}
                    for area in areas[celltype]:
                        pp['opsindict']['ChR'][area] = chrdict
                    if len(areas[celltype])==TOTAL_NUM_AREAS[celltype]:
                        pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
                    else:
                        pp['ChR_areas'] = {'partial'      : pp['opsindict']['ChR'].keys()}
                else:
                    pp['ChR_areas'] = {'none'      : [None]}
                    
                if irr > 0 and factor > 0 : 
                    pp['opsindict']['NpHR'] =  {}
                    for area in areas[celltype]:
                        pp['opsindict']['NpHR'][area] = hrdict
                    if len(areas[celltype])==TOTAL_NUM_AREAS[celltype]:
                        pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
                    else:
                        pp['NpHR_areas'] = {'partial'      : pp['opsindict']['NpHR'].keys()}
                    
                else:
                    pp['NpHR_areas'] = {'none'      : [None]}
                    
                # general settings 
                pp['experiment_type'] = 'opsinonly'
                pp['savedata'] = True # False #True
                
                pp['tstart'] = 0
                pp['tstop'] = tstop
                
                nsites = 1
                labels = ['stim%g'%i for i in range(2*nsites)]
                
                
                pp['mark_loc'] = {}
                pp['mark_loc']['names'] = ['mysoma','distal']+labels
                pp['mark_loc']['sections'] = ['soma','apic']+['apic']*nsites+['dend']*nsites
                pp['mark_loc']['ids'] = [(0,0.5),(0,0.972326)] + [('select_section_posn_bydistance',{'sectionarea':'apic','mindist':dist[0],'maxdist':dist[1]}) for i in range(nsites)] + [('select_section_posn_bydistance',{'sectionarea':'dend','mindist':dist[0],'maxdist':dist[1]}) for i in range(nsites)] 
                
                             
                
                pp['record_loc'] = {}
                pp['record_loc']['v'] = ['mysoma'] #+labels
                pp['record_loc']['ina'] = ['mysoma']
                pp['record_loc']['ik'] = ['mysoma']
                
                pp['stim_iclamp'] = True
                #print labels[0]
                #pp['spiketrains'] = [{'tstims': get_spiketimes(freq,nsites), 'locations': labels, 'weights':np.ones(nsites)*J,  'el': 0.1}]
                
                pp['iclamp'] = [{'tstim': 0,  'location': stimloc, 'amp':Ia, 'duration':tstop}]
                       
                
                
                vplots_soma = [['mysoma','v','k']]
    
                #iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
                pp['plot'] = {1:vplots_soma} 
                    
                
                pp['num_threads'] = 1
                                           
                pp.update({'expname':expbase,
                           'description':'irr%.3f_factor%.2f_I%.2f_stimloc_%s'%(irr,factor,Ia,stimloc)})
    
                es.run_single_experiment(expbase, 'missing', pp)
                #es.run_single_experiment(expbase, 'local', pp)
                #es.run_single_experiment(expbase, 'names', pp)
                count += 1
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


def get_optdescript(irr,factor):
    opt_des = []
    if irr > 0 :
        if len(areas[celltype])==TOTAL_NUM_AREAS[celltype]:
            opt_des.append('whole')
        else:
            opt_des.append('partial')
    else:
        opt_des.append('none')
        
    if irr > 0 and factor > 0:
        if len(areas[celltype])==TOTAL_NUM_AREAS[celltype]:
            opt_des.append('whole')
        else:
            opt_des.append('partial')
    else:
        opt_des.append('none')
    return opt_des

                                
def analyse_locations_optogen_invivo():
    
    
    for f in factors:
        af = run_analysis.AnalyticFrame()
    
        af.update_params({'tstart':light_on,'tstop':light_on+light_dur,
                              'tstart_bg': 50,'tstop_bg':light_on,
                              'tstart_post':light_on+light_dur,'tstop_post':tstop})
        #optlog = get_optdescript(irr,1)
        exp_comp_list = [['irr%.3f_'%irr+'factor%.2f_'%f+'freq%g'+'_J%g_nsites%g'%(Js[0],nsite_range[0])+'_NpHR_%s_ChR_%s'%(get_optdescript(irr,f)[1],get_optdescript(irr,f)[0]),'=%.3f'%irr] for irr in irrs]
    
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
        """
        #af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
        af.submenu_plot(5, expbase+'FI_gain_varyIrr_invivo_factor%g'%f)
        af.submenu_plot(0, expbase+'FI_gain_varyIrr_invivo_factor%g'%f)
        af.submenu_plot(10, expbase+'FI_gain_varyIrr_invivo_factor%g'%f)
        #return
        """
        
    """
    #af.submenu_print()
    #af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
    af.submenu_plot(5, expbase+'FI_gain_varyFactor')
    af.submenu_plot(0, expbase)
    af.submenu_plot(10, expbase)
    """
def analyse_factors_invivo():
    
    
    for irr in irrs:
        af = run_analysis.AnalyticFrame()
    
        af.update_params({'tstart':light_on,'tstop':light_on+light_dur,
                              'tstart_bg': 50,'tstop_bg':light_on,
                              'tstart_post':light_on+light_dur,'tstop_post':tstop})
        #optlog = get_optdescript(irr,1)
        exp_comp_list = [['irr%.3f_'%irr+'factor%.2f_'%f+'freq%g'+'_J%g_nsites%g'%(Js[0],nsite_range[0])+'_NpHR_%s_ChR_%s'%(get_optdescript(irr,f)[1],get_optdescript(irr,f)[0]),'%.2f'%f] for f in factors]
    
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
        
        
        #af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
        af.submenu_plot(5, expbase+'FI_gain_varyFactor_invivo_irr%g'%irr+'_NpHR_%s_ChR_%s'%(get_optdescript(irr,f)[1],get_optdescript(irr,f)[0]))
        af.submenu_plot(0, expbase+'FI_gain_varyFactor_invivo_irr%g'%irr+'_NpHR_%s_ChR_%s'%(get_optdescript(irr,f)[1],get_optdescript(irr,f)[0]))
        af.submenu_plot(10, expbase+'FI_gain_varyFactor_invivo_irr%g'%irr+'_NpHR_%s_ChR_%s'%(get_optdescript(irr,f)[1],get_optdescript(irr,f)[0]))
       
        
    """
    #af.submenu_print()
    #af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
    af.submenu_plot(5, expbase+'FI_gain_varyFactor')
    af.submenu_plot(0, expbase)
    af.submenu_plot(10, expbase)
    """

   
                               
def analyse_locations_optogen_invivo_constant():
    
    
    for irr in irrs:
        af = run_analysis.AnalyticFrame()
    
        af.update_params({'tstart':light_on,'tstop':light_on+light_dur,
                              'tstart_bg': 50,'tstop_bg':light_on,
                              'tstart_post':light_on+light_dur,'tstop_post':tstop})
        #optlog = get_optdescript(irr,1)
        exp_comp_list = [['constant_irr%.3f_'%irr+'factor%.2f_'%f+'freq%g'+'_J%g_nsites%g'%(Js[0],nsite_range[0])+'_NpHR_%s_ChR_%s'%(get_optdescript(irr,f)[1],get_optdescript(irr,f)[0]),'=%.3f'%f] for f in factors]
    
        print exp_comp_list
  
        expss = [ec[0] for ec in exp_comp_list]
        explabels = [ec[1] for ec in exp_comp_list]
        print explabels
        
        af.populate_expset(expbase,expss,explabels, [freqs])
        
        
        af.submenu_extractSpikes()
        
        af.submenu_print()
        """
        af.submenu_runFI()
        for exp in af.experimentset:
            exp.calculate_responses('FI')
            exp.calculate_responses('FI_bg')
            exp.calculate_responses('FI_post')
        
        
        
        af.submenu_save()
        #af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
        af.submenu_plot(5, expbase+'FI_gain_varyFactor_irr%.3f'%irr)
        af.submenu_plot(0, expbase+'FI_gain_varyFactor_irr%.3f'%irr)
        af.submenu_plot(10, expbase+'FI_gain_varyFactor_irr%.3f'%irr)
        """
        
    """
    #af.submenu_print()
    #af.submenu_plot(0, self.expbase+'FI_gain_irr%.2f_tree%s_varyFactor_exp%s%s_'%(0.05,tree,exp[0],exp[1]))
    af.submenu_plot(5, expbase+'FI_gain_varyFactor')
    af.submenu_plot(0, expbase)
    af.submenu_plot(10, expbase)
    """


def analyse_locations_optogen_invitro():
    
    
    for f in [1.]: #factors:
        af = run_analysis.AnalyticFrame()
    
        af.update_params({'tstart':light_on,'tstop':light_on+light_dur,
                              'tstart_bg': 50,'tstop_bg':light_on,
                              'tstart_post':light_on+light_dur,'tstop_post':tstop})
        
        exp_comp_list = [['irr%.3f_'%irr+'factor%.2f_'%f+'I%.2f'+'_stimloc_distal'+'_NpHR_%s_ChR_%s'%(get_optdescript(irr,f)[1],get_optdescript(irr,f)[0]),'%.3f'%irr] for irr in irrs]
        #irr%.3f_factor%.2f_I%.2f_stimloc_%s'

        print exp_comp_list
  
        expss = [ec[0] for ec in exp_comp_list]
        explabels = [ec[1] for ec in exp_comp_list]
        print explabels
        
        af.populate_expset(expbase,expss,explabels, [iclamp_amps])
        
        
        af.submenu_extractSpikes()
        
        af.submenu_print()
        
        """
        af.submenu_runFI()
        for exp in af.experimentset:
            exp.calculate_responses('FI')
            exp.calculate_responses('FI_bg')
            exp.calculate_responses('FI_post')
            #return
        af.submenu_save()
        af.submenu_print()
        
        
        af.submenu_plot(5, expbase+'FI_gain_varyIrr')
        af.submenu_plot(0, expbase+'FI_gain_varyIrr')
        af.submenu_plot(10, expbase+'FI_gain_varyIrr')
        """
        

def analyse_factors_invitro():
    
    
    for irr in irrs:
        af = run_analysis.AnalyticFrame()
    
        af.update_params({'tstart':light_on,'tstop':light_on+light_dur,
                              'tstart_bg': 50,'tstop_bg':light_on,
                              'tstart_post':light_on+light_dur,'tstop_post':tstop})
        
        exp_comp_list = [['irr%.3f_'%irr+'factor%.2f_'%f+'I%.2f'+'_stimloc_distal'+'_NpHR_%s_ChR_%s'%(get_optdescript(irr,f)[1],get_optdescript(irr,f)[0]),'%.2f'%f] for f in factors]
        #irr%.3f_factor%.2f_I%.2f_stimloc_%s'

        print exp_comp_list
  
        expss = [ec[0] for ec in exp_comp_list]
        explabels = [ec[1] for ec in exp_comp_list]
        print explabels
        
        af.populate_expset(expbase,expss,explabels, [iclamp_amps])
        
        
        af.submenu_extractSpikes()
        
        af.submenu_print()
        
        
        af.submenu_runFI()
        for exp in af.experimentset:
            exp.calculate_responses('FI')
            exp.calculate_responses('FI_bg')
            exp.calculate_responses('FI_post')
            #return
        af.submenu_save()
        af.submenu_print()
        
        
        af.submenu_plot(5, expbase+'FI_gain_varyFactor_invitro_irr%g'%irr)
        af.submenu_plot(0, expbase+'FI_gain_varyFactor_invitro_irr%g'%irr)
        af.submenu_plot(10, expbase+'FI_gain_varyFactor_invitro_irr%g'%irr)
       
    
    
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print 'Need to specify an input mode'
    elif sys.argv[1] == 'scan_vivo':
        scan_locations_optogen_invivo()
    elif sys.argv[1] == 'scan_vitro':
        scan_locations_optogen_invitro(sys.argv[2])        
    elif sys.argv[1] == 'constant':
        scan_locations_optogen_invivo_same_total()
    elif sys.argv[1] == 'analyse':
        analyse_locations_optogen_invivo()
    elif sys.argv[1] == 'analyse_irrs':
        analyse_factors_invivo()
    elif sys.argv[1] == 'analyse_constant':
        analyse_locations_optogen_invivo_constant()        
    elif sys.argv[1] == 'analyse_vitro':
        analyse_locations_optogen_invitro()
    elif sys.argv[1] == 'analyse_vitro_irrs':
        analyse_factors_invitro()
                                            
                           
