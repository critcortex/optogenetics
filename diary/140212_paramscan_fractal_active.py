import sys
import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import run_stimulation
import pylab
import file_io as fio


expbase = '140212_fractal_paramscan_active'




length = 50
nseg = 10

es = ExpSetup()

stim_location_fracts = [['peripheral',[0.5,1.0]],
             ['distributed',[0,1.0]],
             ['proximal',[0,.5]]]

tstop = 1000

freq = 50.
dist_frac_sections = 50.


num_base = range(1,5)#7)
num_child = range(1,5)
num_levels = range(1,5)

num_base = [1,2,4]
num_child = [1,2,4]
num_levels = [4]

light_on = 250
light_dur = 500
factors = [0.125,0.25,0.5,1.,2.,4.,8.]
factors = [0.25,1.,4.]
#factors = [0.,0.0625,0.125]
factors = [0.01,0.03]
factors = [0.25,1.,4.]+[0.,0.0625,0.125]+[0.01,0.03]
#irrs = [1.]
irrs = np.arange(1.,5.1,0.5)


num_levels = range(2,4)
num_child = [1,2,4]
num_base = [2]
irrs = [1.]
factors = [0.0]

"""
# debug params
tstop = 500
light_on = 0
light_dur = tstop
num_levels = [4]
num_child = [3]
num_base = [7]
iclamp_amps = [2.0]
clamp_tstim = 200
clamp_tdur = 100
num_levels=[5]
irrs = [1.0]
factors = [1.]
"""
#num_base = range(3,7)

def get_dendritic_colors(num_levels):
    """
    soma = level0 will always be black
    """
    num_levels = num_levels +2 # to account for soma and that we don't want to go to complete white
    return ['#333333','#666666','#888888','#aaaaaa','#dddddd']
    


def get_dend0_child_locations(num_level):
    labels = []
    locs = []
    for i in range(1,num_level+1):
        tmp = 'dend'+'_'.join(['0']*i)
        labels.append(tmp)
        locs.append((tmp,0.5))
    return labels,locs


def scan_locations_optogen():
    #iclamp_amps = [0]
    count = 0
    for nb in num_base:
        #print 'nb = ',nb
        for nc in num_child:
            #print 'nc = ',nc
            for nl in num_levels:
                #print 'nl = ',nl
                labels,locs = get_dend0_child_locations(nl)
        
                for factor in factors:
                    for irr in irrs: 
                        
                        pp = {}
                        # neuron model and params
                        pp['cell'] = ['Neuron', 'FractalNeuron']
                        pp['cell_params'] = {'num_base': nb,
                              'num_split': nc,
                              'num_levels':nl}
                        
                        # opsin
                        chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                        hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                        
                        pp['opsindict'] = {}
                        if irr > 0:
                            pp['opsindict']['ChR'] =  {'soma': chrdict}
                            for i in range(nb):
                                pp['opsindict']['ChR']['dend%g'%i] = chrdict
                                
                            pp['opsindict']['NpHR'] =  {'soma': hrdict}
                            for i in range(nb):
                                pp['opsindict']['NpHR']['dend%g'%i] = hrdict
                            pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
                            pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
                        else:
                            pp['NpHR_areas'] = {'none'      : [None]}
                            pp['ChR_areas'] = {'none'      : [None]}
                            
                        # general settings 
                        pp['experiment_type'] = 'opsinonly'
                        pp['savedata'] = True # False #True
                        
                        pp['tstart'] = 0
                        pp['tstop'] = tstop
                        
                        otherbaseslabels,othersections = [],[]
                        if nb > 1:
                            otherbaseslabels = ['dend1']
                            othersections    = ['dend1']
                        
                        pp['mark_loc'] = {}
                        pp['mark_loc']['names'] = ['mysoma']+labels+otherbaseslabels
                        pp['mark_loc']['sections'] = ['soma']+['dend0']*nl+othersections
                        pp['mark_loc']['ids'] = [(0,0.5)] + [('get_section_byname',{'sectioname':loc[0],'subdomain':'dend0'}) for loc in locs]+[('get_section_byname',{'sectioname':loc,'subdomain':'dend1'}) for loc in otherbaseslabels]
                        
                        
                        pp['record_loc'] = {}
                        pp['record_loc']['v'] = ['mysoma']+labels+otherbaseslabels
                        pp['record_loc']['ina'] = ['mysoma']
                        pp['record_loc']['ik'] = ['mysoma']
                        
                        
                        vplots_soma = [['mysoma','v','k']]
                        vplots = [['mysoma','v','k']]
                        colors = get_dendritic_colors(nl)
                        for i in range(nl):
                            vplots.append([labels[i],'v',colors[i]])
                        #print vplots
                        if nb>1:
                            vplots.append([otherbaseslabels[0],'v','b-'])
                        iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
                        pp['plot'] = {1:vplots,
                                      2:vplots_soma,
                                      3:iplots_soma}
                        
                        
                        pp['num_threads'] = 1
                                                   
                        pp.update({'expname':expbase,
                                   'description':'irr%.1f_factor%.2f_nb%g_ns%g_nl%g'%(irr,factor,nb,nc,nl)})
            
                        #print 'num levels = ', nl
                        #print 'stim-location = ',labels[iclamp_index]
                        
                        #print 'going to run on cluster',pp['description']
                        es.run_single_experiment(expbase, 'cluster', pp)
                        #es.run_single_experiment(expbase, 'missing', pp)
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
                                
                                
def scan_locations_optogen_partial(same_as_stimloc=True):

    #factors = [0.00]
    #num_base = [2] ############################################################################
    count = 0
    for nb in num_base:
        if nb <= 1:
            # no use doing this for when there is only one base
            continue
        #print 'nb = ',nb
        for nc in num_child:
            #print 'nc = ',nc
            for nl in num_levels:
                #print 'nl = ',nl
                if nl == 1 and nc > 1:
                    continue # again, no point
                
                labels,locs = get_dend0_child_locations(nl)
                
                for factor in factors:
                    for irr in irrs: 
                                
                        #print 'ic_index = ',iclamp_index

                        pp = {}
                        # neuron model and params
                        pp['cell'] = ['Neuron', 'FractalNeuron']
                        pp['cell_params'] = {'num_base': nb,
                              'num_split': nc,
                              'num_levels':nl}
                        
                        # opsin
                        chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                        hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                        opsin_index = 0
                        pp['NpHR_areas'] = {'partialSame':['dend%g'%opsin_index,'soma']}
                        pp['ChR_areas'] = {'partialSame':['dend%g'%opsin_index,'soma']}
                        
                        pp['opsindict'] = {}
                        pp['opsindict']['ChR'] =  {'soma': chrdict}
                        pp['opsindict']['ChR']['dend%g'%opsin_index] = chrdict
                            
                        pp['opsindict']['NpHR'] =  {'soma': hrdict}
                        pp['opsindict']['NpHR']['dend%g'%opsin_index] = hrdict
                            
                        # general settings 
                        pp['experiment_type'] = 'opsinonly'
                        pp['savedata'] = True # False #True
                        
                        pp['tstart'] = 0
                        pp['tstop'] = tstop
                        
                        
                        otherbaseslabels,othersections = [],[]
                        if nb > 1:
                            otherbaseslabels = ['dend1']
                            othersections    = ['dend1']
                        if nb > 2 and not same_as_stimloc:
                            otherbaseslabels.append('dend2')
                            othersections.append('dend2')
                        
                        pp['mark_loc'] = {}
                        pp['mark_loc']['names'] = ['mysoma']+labels+otherbaseslabels
                        pp['mark_loc']['sections'] = ['soma']+['dend0']*nl+othersections
                        pp['mark_loc']['ids'] = [(0,0.5)] + [('get_section_byname',{'sectioname':loc[0],'subdomain':'dend0'}) for loc in locs]+[('get_section_byname',{'sectioname':loc,'subdomain':loc}) for loc in otherbaseslabels]
                        
                        pp['record_loc'] = {}
                        pp['record_loc']['v'] = ['mysoma']+labels+otherbaseslabels
                        pp['record_loc']['ina'] = ['mysoma']
                        pp['record_loc']['ik'] = ['mysoma']
                        
                       
                        vplots_soma = [['mysoma','v','k']]
                        vplots = [['mysoma','v','k']]
                        colors = get_dendritic_colors(nl)
                        for i in range(nl):
                            vplots.append([labels[i],'v',colors[i]])
                        #print vplots
                        if nb>1:
                            vplots.append([otherbaseslabels[0],'v','b-'])
                        if nb>2 and not same_as_stimloc:
                            vplots.append([otherbaseslabels[1],'v','g-'])
                        iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
                        pp['plot'] = {1:vplots,
                                      2:vplots_soma,
                                      3:iplots_soma}
                        
                        
                        pp['num_threads'] = 1
                                                   
                        pp.update({'expname':expbase,
                                   'description':'irr%.1f_factor%.2f_nb%g_ns%g_nl%g'%(irr,factor,nb,nc,nl)})
            
                        #print 'num levels = ', nl
                        #print 'stim-location = ',labels[iclamp_index]
                        
                        #print 'going to run on cluster', expbase, pp['description']
                        es.run_single_experiment(expbase, 'cluster', pp)
                        #es.run_single_experiment(expbase, 'missing', pp)
                        count += 1
                        #es.run_single_experiment(expbase, 'local', pp)
                        #return pp 
                        """
                        NE = run_stimulation.NeuronExperiment()
                        ES = ExpSetup()
                        dp = ES.get_default_params()
                        dp.update(pp)
                        NE.main(dp)
                        return pp
                
                        """
    print '%g jobs submitted'%count                                





def analyse():
    import run_analysis 
    af = run_analysis.AnalyticFrame()
    af.update_params({'tstart':250,'tstop':750,'label_format':'irr%.1f_factor%.2f_nb%g_ns%g_nl%g'})
    nbs = [1,2,4]
    nss = [1,2,4]
    nls = [4]
    #num_base = range(3,7)
    optlocations = ['whole','partialSame']
    #optlocations = ['whole']
    #optlocations = ['none']
    optlocations = ['partialSame']
    #exp_comp_list = [['_irr%.1f_factor%.2f_freq%g_J%1.f_nstim%g'%(irr,factor,freq,J,nstim),'ext_freq=%g'%(freq)] for freq in freqs ] #for J in Js for nstim in nstims]
    exp_comp_list = [['irr%.1f_factor%.2f'+'_nb%g'%nb+'_ns%g_nl%g'+'_NpHR_%s_ChR_%s'%(optlog,optlog),'nb=%g, %s'%(nb,optlog)] for nb in nbs for optlog in optlocations]
    #exp_comp_list = [['irr%.1f_factor%.2f'+'_nb%g'%nb+'_ns%g_nl%g_iclamploc%g_I%.1f'+'_NpHR_%s_ChR_%s'%(optlog,optlog),'nb=%g, %s'%(nb,optlog)] for nb in num_base for optlog in optlocations]
    print exp_comp_list
    
    expss = [ec[0] for ec in exp_comp_list]
    explabels = [ec[1] for ec in exp_comp_list]
    print explabels
    #num_locations = range(0,len(num_levels))
    af.populate_expset(expbase,expss,explabels, [irrs,factors,nss,nls])
    
    af.submenu_extractSpikes()
    #af.submenu_plot(6,'')
    #af.submenu_plot(7,'raster_normal',)
    #af.submenu_plot(7,'raster_aligned',align=light_durations,offset=100-pre_irr_duration,tmax=800)
    #af.submenu_plot(8,'voltage_irroff_aligned',align=light_durations,offset=100-pre_irr_duration,tmax=800)
    #af.perform_analysis([])
    af.perform_analysis(['cv_isi','fano_factor_isi','mean_rate','isi','cv_kl'])
    
        







    
def plot_gain():
    import run_analysis
    """
    for irr in [1,5]:
        for factor in factors:
            af = run_analysis.AnalyticFrame()
            af.update_params({'tstart':light_on,'tstop':light_on+light_dur})
            
            exp_comp_list = [['_irr%.1f_factor%.2f'%(irr,factor)+'_Idist%.1f_'+'NpHR_%s_ChR_%s'%(exp[1],exp[0]),'ChR %s, NpHR %s'%(exp[0],exp[1])] for exp in explist ]
            print exp_comp_list
            
            expss = [ec[0] for ec in exp_comp_list]
            explabels = [ec[1] for ec in exp_comp_list]
            af.populate_expset(expbase,expss,explabels,[current_amps])
            
            af.submenu_load()
            af.submenu_print()
            af.submenu_plot(0, expbase+'FI_gain_irr%.1f_factor%.2f'%(irr,factor))
     """
     
    explist = [['whole','whole']] 
    explist = [['partialSame','partialSame']]
    #explist = [['apical','apical'],['whole','whole'] ] 
    print irrs
    factors.sort()
    print factors
    nbs = [1,2,4]
    nss = [1,2,4]
    nls = [4]
    
    for exp in explist:
        for nl in nls:
            for ns in nss:
                for nb in nbs:
                        
                    af = run_analysis.AnalyticFrame()
                    af.update_params({'tstart':250,'tstop':750})
 
                    exp_comp_list = [['irr%.1f'+'_factor%.2f'%(factor)+'_nb%g'%nb+'_ns%g_nl%g_'%(ns,nl)+'NpHR_%s_ChR_%s'%(exp[1],exp[0]),'%g'%factor] for factor in factors]
                    print exp_comp_list
                    
                    expss = [ec[0] for ec in exp_comp_list]
                    explabels = [ec[1] for ec in exp_comp_list]
                    af.populate_expset(expbase,expss,explabels,[irrs])
                    
                    af.submenu_load()
                    af.submenu_runFI()
                    af.submenu_print()
                    af.submenu_plot(5, expbase+'FI_gain_nb%g_ns%g_nl%g_varyFactor_exp%s%s_'%(nb,ns,nl,exp[0],exp[1]))
        













                            
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        pass
    elif sys.argv[1] == 'whole':
        scan_locations_optogen()
    elif sys.argv[1] == 'partial':
        if len(sys.argv)==3:
            scan_locations_optogen_partial(sys.argv[2]=='True')
        else:
            scan_locations_optogen_partial()
    elif sys.argv[1] == 'analyse':
        analyse()
    elif sys.argv[1] == 'plotgain':
        plot_gain()





"""
def get_spiketimes(rate,num_inputs,tstop):
    stg = ntst.StGen()
    spikes = []
    for i in range(num_inputs):
        ss = stg.poisson_generator(rate,t_stop=tstop,array=True)
        spikes.append(ss)
    #print 'Spikes = ',spikes
    return spikes

def add_locations(num_loc,label,upper_id):
    sitelabel='mysite'
    labels = [sitelabel+str(i) for i in range(num_loc)]
    ids = range(upper_id)
    np.random.shuffle(ids)
    locs = ids[:num_loc]
    locs = [(l,0.5) for l in locs]
    return labels,locs


def scan_locations_spiketrain():

    for nb in num_base:
        for nc in num_child:
            for nl in num_levels:
                
                num_locations = range(nl)
                # distance of each segment is 50, therefore we can select by distance 
                iclamp_locations = [ (-1,1.1) ]
                #select_segments_bydistance(self,section,mindist=0,maxdist=-1)
                
                for ic in iclamp_amps:
                    for iclamp_loc in iclamp_locations:
                
                
                        pp = {}
                        # neuron model and params
                        pp['cell'] = ['Neuron', 'FractalNeuron']
                        pp['cell_params'] = {'num_base': nb,
                              'num_split': nc,
                              'num_levels':nl}
                        
                        # opsin
                        chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                        hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                        
                        pp['opsindict'] = {}
                        pp['opsindict']['ChR'] =  {'soma': chrdict}
                        for i in range(nb):
                            pp['opsindict']['ChR']['dend%g'%i] = chrdict
                            
                        pp['opsindict']['NpHR'] =  {'soma': hrdict}
                        for i in range(nb):
                            pp['opsindict']['NpHR']['dend%g'%i] = hrdict
                            
                        # general settings 
                        pp['experiment_type'] = 'opsinonly'
                        pp['savedata'] = True
                        
                        pp['tstart'] = 0
                        pp['tstop'] = tstop
                        
                        
                        
                        pp['stim_spiketrains'] = True          
                        pp['mark_loc'] = {}
                        ##labels,locs = add_locations(num_stim_locations,109)
                        ##pp['spiketrains'] = [{'tstims': get_spiketimes(freq,num_stim_locations,tstop),  'locations': labels, 'weights':np.ones(num_stim_locations)*J,  'el': 0.02}]
                        description = 'spiketrain_loc%g_freq%.2f'
                        
                                    
                        # TODO work out how we want to have different                             
                        pp['mark_loc'] = {}
                        pp['mark_loc']['names'] = ['mysoma','iclamp_loc']
                        pp['mark_loc']['sections'] = ['soma','dend0']
                        pp['mark_loc']['ids'] = [(0,0.5),('get_section_byname',{'sectioname':'dend0_0','subdomain':'dend0'})]
                        
                        pp['record_loc'] = {}
                        pp['record_loc']['v'] = ['mysoma']
                        
                        pp['num_threads'] = 1
                                                   
                        pp.update({'expname':expbase,
                                   'description':'_nb%g_ns%g_nl%g_%s'%(nb,nc,nl,description)})
            
                        print 'going to run on cluster'
                        #es.run_single_experiment(expbase, 'cluster', pp)
                        es.run_single_experiment(expbase, 'local', pp)
                        return
                        
                       
"""
            