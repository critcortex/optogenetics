import sys
import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import run_stimulation
import pylab
import file_io as fio


expbase = '140115_fractal_paramscan'




length = 50
nseg = 10

es = ExpSetup()

stim_location_fracts = [['peripheral',[0.5,1.0]],
             ['distributed',[0,1.0]],
             ['proximal',[0,.5]]]

tstop = 500

freq = 50.
dist_frac_sections = 50.

clamp_tstim = 200
clamp_tdur = 100
iclamp_amps = np.arange(0.5,2.6,0.5)

num_base = range(1,7)
num_child = range(1,5)
num_levels = range(1,5)

light_on = 0
light_dur = tstop
factors = [0.125,0.25,0.5,1.,2.,4.,8.]
factors = [0.25,1.,4.]
irrs = [1.]
irrs = np.arange(2.,6.1,1.0)
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




def scan_locations_clamp():
    #iclamp_amps = [0]
    count = 0
    for nb in num_base:
        #print 'nb = ',nb
        for nc in num_child:
            #print 'nc = ',nc
            for nl in num_levels:
                #print 'nl = ',nl
                labels,locs = get_dend0_child_locations(nl)
                
                for ic in iclamp_amps:
                    #print 'ic = ',ic
                    for iclamp_index in range(nl):
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
                                
                                pp['opsindict'] = {}
                                if irr > 0:
                                    pp['opsindict']['ChR'] =  {'soma': chrdict}
                                    for i in range(nb):
                                        pp['opsindict']['ChR']['dend%g'%i] = chrdict
                                        
                                    pp['opsindict']['NpHR'] =  {'soma': hrdict}
                                    for i in range(nb):
                                        pp['opsindict']['NpHR']['dend%g'%i] = hrdict
                                        
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
                                
                                pp['stim_iclamp'] = True
                                pp['iclamp'] = [{'amp':ic,'tstim':clamp_tstim, 'duration':clamp_tdur,'location':labels[iclamp_index]}]
                                
                                vplots = [['mysoma','v','k']]
                                colors = get_dendritic_colors(nl)
                                for i in range(nl):
                                    if i==iclamp_index:
                                        vplots.append([labels[i],'v','r-'])
                                    else:
                                        vplots.append([labels[i],'v',colors[i]])
                                #print vplots
                                if nb>1:
                                    vplots.append([otherbaseslabels[0],'v','b-'])
                                pp['plot'] = {1:[['mysoma','v','k-'],[labels[iclamp_index],'v','r-']],
                                              2:vplots}
                                
                                
                                pp['num_threads'] = 1
                                                           
                                pp.update({'expname':expbase,
                                           'description':'irr%.1f_factor%.2f_nb%g_ns%g_nl%g_iclamploc%g_I%.1f'%(irr,factor,nb,nc,nl,iclamp_index,ic)})
                    
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
                                
                                
def scan_locations_clamp_partial(same_as_stimloc=True):
    
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
                labels,locs = get_dend0_child_locations(nl)
                
                for ic in [0]: ############################################# iclamp_amps:
                    #print 'ic = ',ic
                    for iclamp_index in [0]: ############################################# range(0,nl): ##############################################################
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
                                
                                if same_as_stimloc:
                                    opsin_index = 0
                                    pp['NpHR_areas'] = {'partialSame':['dend%g'%opsin_index,'soma']}
                                    pp['ChR_areas'] = {'partialSame':['dend%g'%opsin_index,'soma']}
                                else:
                                    opsin_index = 1
                                    pp['NpHR_areas'] = {'partialDiff':['dend%g'%opsin_index,'soma']}
                                    pp['ChR_areas'] = {'partialDiff': ['dend%g'%opsin_index,'soma']}

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
                                
                                pp['stim_iclamp'] = True
                                pp['iclamp'] = [{'amp':ic,'tstim':clamp_tstim, 'duration':clamp_tdur,'location':labels[iclamp_index]}]
                                
                                vplots = [['mysoma','v','k']]
                                colors = get_dendritic_colors(nl)
                                for i in range(nl):
                                    if i==iclamp_index:
                                        vplots.append([labels[i],'v','r-'])
                                    else:
                                        vplots.append([labels[i],'v',colors[i]])
                                #print vplots
                                if nb>1:
                                    vplots.append([otherbaseslabels[0],'v','b-'])
                                if nb>2 and not same_as_stimloc:
                                    vplots.append([otherbaseslabels[1],'v','g-'])
                                pp['plot'] = {1:[['mysoma','v','k-'],[labels[iclamp_index],'v','r-']],
                                              2:vplots}
                                
                                
                                pp['num_threads'] = 1
                                                           
                                pp.update({'expname':expbase,
                                           'description':'irr%.1f_factor%.2f_nb%g_ns%g_nl%g_iclamploc%g_I%.1f'%(irr,factor,nb,nc,nl,iclamp_index,ic)})
                    
                                #print 'num levels = ', nl
                                #print 'stim-location = ',labels[iclamp_index]
                                
                                #print 'going to run on cluster', expbase, pp['description']
                                #es.run_single_experiment(expbase, 'cluster', pp)
                                es.run_single_experiment(expbase, 'missing', pp)
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





def location_vs_branching(illumination='whole',nb=1):
    cols = ['r','y','g','b']
    pylab.figure()
    expfiles = 'irr%.1f_factor%.2f_nb%g_ns%g_nl%g_iclamploc%g_I%.1f_NpHR_%s_ChR_%s'
    irr = 1. 
    factor = 1. 
    
    for (i,nl) in enumerate(num_levels):
        for ns in num_child:
            for clamp_loc in range(nl):
                shiftdx = (i%2)*0.1 - 0.05
                shiftdy = (i/2)*0.1 - 0.05
                for ic in iclamp_amps:
                    #print expfiles, (irr,factor,nb,ns,nl,clamp_loc,ic,illumination,illumination)
                    ss = fio.loadspikes(expfiles%(irr,factor,nb,ns,nl,clamp_loc,ic,illumination,illumination), expbase)
                    #print ss
                    #print type(ss)
                    print ss.size, expfiles%(irr,factor,nb,ns,nl,clamp_loc,ic,illumination,illumination)
                    if ss.size>0: 
                        pylab.scatter(clamp_loc+shiftdx,ns+shiftdy,c=cols[i],s=ss.size**2,marker='o',edgecolors='none',alpha=0.33)
    
    pylab.xlabel('current injection site')
    pylab.ylabel('amount of branches')
    pylab.savefig('test_nb%g_ill%s.png'%(nb,illumination))

def __turn_off_border(ax,turnoffs=['right','top']):
        for loc, spine in ax.spines.iteritems():
            if loc not in turnoffs:
                spine.set_position(('outward',5))
                ax.tick_params(direction='out')
            elif loc in turnoffs:
                #spine.set_color('none') # don't draw spine
                spine.set_visible(False)
                #ax.tick_params([])
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')


def get_clamploc_vs_spikes(ic,ill='whole',nb=1):
    expfiles = '*nb%g_*_iclamploc%g_I%.1f*_NpHR_%s_ChR_%s'
    maxes = []
    clamp_locs = range(0,4)
    for cl in clamp_locs:
        
        # find best for each #nb
        spikes,filenames=fio.loadspikerange(expfiles%(nb,cl,ic,ill,ill),expbase)
        # find max
        sp_sizes = [sp.size for sp in spikes]
        max = np.max(sp_sizes)
        print max, " for file", filenames[sp_sizes.index(max)] 
        maxes.append(max)
    return maxes,clamp_locs


def multi_clamploc_spikes(ills=['whole'],nb=[1]):
    colors = np.linspace(0, len(iclamp_amps)+1, 10)
    mymap = pylab.get_cmap("OrRd")
    # get the colors from the color map
    my_colors = mymap(colors)
    fig = pylab.figure()
    axes = fig.add_subplot(111)
    for (i,ill) in enumerate(ills):
        for (j,ic) in enumerate(iclamp_amps):
            maxes,bases = get_clamploc_vs_spikes(ic,ill,nb[i])
            print j,my_colors[j]
            axes.plot(bases,maxes,c=my_colors[ic],marker='.',ms=3,lw=4)
    pylab.xlabel('Clamp location')
    pylab.ylabel('#spikes')
    __turn_off_border(axes,turnoffs=['right','top'])
    pylab.ylim(ymin=0)
    fig.savefig('spike_clamplocation_%s.png'%'_'.join([ills[i]+str(nb[i]) for i in range(len(ills))]))




def multi_nb_spikes(ills=['whole'],colors=['r'],starts=[1]):
    fig = pylab.figure()
    axes = fig.add_subplot(111)
    for (i,ill) in enumerate(ills):
        maxes,bases = get_nb_vs_spikes(ill,starts[i])
        axes.plot(bases,maxes,c=colors[i],marker='.',ms=3,lw=4)
    pylab.xlabel('Polarity')
    pylab.ylabel('#spikes')
    __turn_off_border(axes,turnoffs=['right','top'])
    pylab.ylim(ymin=0)
    fig.savefig('spike_nb_%s.png'%'_'.join(ills))
    
def get_nb_vs_spikes(ill='whole',start=1):
    expfiles = '*_nb%g_*_NpHR_%s_ChR_%s'
    maxes = []
    num_base = range(start,7)
    for nb in num_base:
        
        # find best for each #nb
        spikes,filenames=fio.loadspikerange(expfiles%(nb,ill,ill),expbase)
        # find max
        sp_sizes = [sp.size for sp in spikes]
        max = np.max(sp_sizes)
        print max, " for file", filenames[sp_sizes.index(max)] 
        maxes.append(max)
    return maxes,num_base



def analyse():
    import run_analysis 
    af = run_analysis.AnalyticFrame()
    af.update_params({'tstart':clamp_tstim,'label_format':'irr%.1f_factor%.2f_ns%g_nl%g_iclamploc%g_I%.1f','tstop':clamp_tstim+clamp_tdur})
    
    #irrs = [0]
    iclamp_amps = [0]
    num_locations = [0]
    #num_base = range(3,7)
    optlocations = ['whole','partialDiff','partialSame']
    #optlocations = ['whole']
    #optlocations = ['none']
    #optlocations = ['partialSame']
    #exp_comp_list = [['_irr%.1f_factor%.2f_freq%g_J%1.f_nstim%g'%(irr,factor,freq,J,nstim),'ext_freq=%g'%(freq)] for freq in freqs ] #for J in Js for nstim in nstims]
    exp_comp_list = [['irr%.1f_factor%.2f'+'_nb%g'%nb+'_ns%g_nl%g_iclamploc%g_I%.1f'+'_NpHR_%s_ChR_%s'%(optlog,optlog),'nb=%g, %s'%(nb,optlog)] for nb in num_base for optlog in optlocations]
    #exp_comp_list = [['irr%.1f_factor%.2f'+'_nb%g'%nb+'_ns%g_nl%g_iclamploc%g_I%.1f'+'_NpHR_%s_ChR_%s'%(optlog,optlog),'nb=%g, %s'%(nb,optlog)] for nb in num_base for optlog in optlocations]
    print exp_comp_list
    
    expss = [ec[0] for ec in exp_comp_list]
    explabels = [ec[1] for ec in exp_comp_list]
    print explabels
    #num_locations = range(0,len(num_levels))
    af.populate_expset(expbase,expss,explabels, [irrs,factors,num_child,num_levels,num_locations,iclamp_amps])
    
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
    explist = [['partialDiff','partialDiff']]
    explist = [['whole','whole'],['partialSame','partialSame'],['partialDiff','partialDiff']]
    #explist = [['apical','apical'],['whole','whole'] ] 
    
    factors.sort()
    print factors
    nbs = [1,2,4]
    nss = [1,2,4]
    nls = [4]
    iclamploc = 3
    irrs = [1]
    iclamps = iclamp_amps
    
    for exp in explist:
        for nl in nls:
            for ns in nss:
                for nb in nbs:
                        
                    af = run_analysis.AnalyticFrame()
                    af.update_params({'tstart':250,'tstop':750})
 
                    exp_comp_list = [['irr%.1f'%irrs[0]+'_factor%.2f'%(factor)+'_nb%g'%nb+'_ns%g_nl%g_iclamploc%g'%(ns,nl,iclamploc)+'_I%.1f_'+'NpHR_%s_ChR_%s'%(exp[1],exp[0]),'%g'%factor] for factor in factors]
                    print exp_comp_list
                    
                    expss = [ec[0] for ec in exp_comp_list]
                    explabels = [ec[1] for ec in exp_comp_list]
                    af.populate_expset(expbase,expss,explabels,[iclamps])
                    
                    af.submenu_load()
                    af.submenu_runFI()
                    af.submenu_print()
                    af.submenu_plot(5, expbase+'FI_gain_nb%g_ns%g_nl%g_varyFactor_exp%s%s_'%(nb,ns,nl,exp[0],exp[1]))
        
















                            
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        pass
    elif sys.argv[1] == 'clamp':
        scan_locations_clamp()
    elif sys.argv[1] == 'partial':
        if len(sys.argv)==3:
            scan_locations_clamp_partial(sys.argv[2]=='True')
        else:
            scan_locations_clamp_partial()
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
            