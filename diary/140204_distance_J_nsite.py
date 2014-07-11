import run_analysis
import sys
import NeuroTools.stgen as ntst 
import run_stimulation
from run_experiments import ExpSetup
import numpy as np

expbase = '140204_distance_J_nsites'


tstop = 1200 #2500
light_on = 0 #500
pw = 1000 # pulse width

es = ExpSetup()
areas = es.get_areas_section()

pd = 1.
Jmax = 5.1
Js = np.arange(1.,Jmax,1.)
freqs = [50]
nstims = range(100,401,100)
factors = [0.25,1.,4.]
factors = [0.25,4.]
#factors = [1.]
irrs = [1.]
dist_choice = [[0,-1],[600,-1]]#,[800,-1]]
trials = range(5)

stg = ntst.StGen()

def get_spiketimes(rate,num_inputs):
    spikes = []
    for i in range(num_inputs):
        ss = stg.poisson_generator(rate,t_stop=tstop,array=True)
        spikes.append(ss)
    #print 'Spikes = ',spikes
    return spikes

def add_locations(num_loc,upper_id,label='apic'):
    labels = [label+str(i) for i in range(num_loc)]
    ids = range(upper_id)
    np.random.shuffle(ids)
    locs = ids[:num_loc]
    locs = [(l,0.5) for l in locs]
    return labels,locs

    
def main(runon):
    count = 0
    params = {}
    #print explist, Js, freqs,nstims,irrs,factors
    for J in Js:
        for freq in freqs:
            for nsites in nstims:
                for irr in irrs:
                    for factor in factors:
                        for dist in dist_choice:
                            for trial in trials:
                            
                                
                                params['num_threads'] = 1
                                
                                params['tstart'] = 0
                                params['tstop'] = tstop
                                
                                params['stim_spiketrains'] = True                            
                                labels = ['stim%g'%i for i in range(nsites)]
                                params['spiketrains'] = [{'tstims': get_spiketimes(freq,nsites),  'locations': labels, 'weights':np.ones(nsites)*J,  'el': 0.02}]
                                
                                
                                params['experiment_type'] = 'opsinonly'
                                params['savedata'] = True
                                
                                params['mark_loc'] = {}
                                params['mark_loc']['names'] = ['mysoma','myapic','proximal','distal']+labels
                                params['mark_loc']['sections'] = ['soma','apic','apic','apic']+['apic']*nsites
                                params['mark_loc']['ids'] = [(0,0.5),(0,0.5),(0,0.0753),(0,0.972326)] + [('select_section_posn_bydistance',{'sectionarea':'apic','mindist':dist[0],'maxdist':dist[1]}) for i in range(nsites)] 
                                
                                
                                params['record_loc'] = {}
                                params['record_loc']['v'] = ['mysoma','myapic','proximal','distal']+labels[:5]
                                
                                chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': pw,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                                hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': pw,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
                                
                                params['opsindict'] = { 'ChR' : {'soma': chrdict,
                                                                 'dend':chrdict,
                                                                 'axon': chrdict,
                                                                 'apic':chrdict},
                                                        'NpHR' : {'soma': hrdict,
                                                                 'dend' : hrdict,
                                                                 'axon':hrdict,
                                                                 'apic':hrdict}}
                                                                
                                
                                params.update({
                                               'expname':expbase,
                                               'description':'_irr%.1f_factor%.2f_freq%g_J%.1f_nstim%g_mindist%g_maxdist%g_trial%g'%(irr,factor,freq,J,nsites,dist[0],dist[1],trial)})
                                params['cell'] = ['Neuron','L5PC']
                                params['cell_description'] = 'test'
                                
                                count +=1
                                
                                if runon =='cluster':
                                    #print 'going to run on cluster'
                                    es.run_single_experiment(expbase, runon, params)
                                elif runon == 'local':
                                    pp = es.get_default_params()
                                    NE = run_stimulation.NeuronExperiment()
                                    pp.update(params)
                                    NE.main(pp)
                                    return
    print 'In total, submitted/ran %g jobs'%count
                            
                            
def analyse():
    factor = 1.
    irr = 0.
    freq = 50.
    af = run_analysis.AnalyticFrame()
    af.update_params({'tstart':200,'tstop':tstop,'label_format':'J=%g, nstim=%g %g %g %g'})
        
    exp_comp_list = [['_irr%.1f_factor%.2f_freq%g_J%s_nstim%s_mindist%s_maxdist%s_trial%s_NpHR_whole_ChR_whole'%(irr,factor,freq,'%.1f','%g','%g','%g','%g'),'ext_freq=%g'%(freq)] for freq in freqs ] 
    print exp_comp_list
    
    expss = [ec[0] for ec in exp_comp_list]
    explabels = [ec[1] for ec in exp_comp_list]
    print explabels
    
    mindists = [dc[0] for dc in dist_choice]
    maxdists = [dc[1] for dc in dist_choice]
    af.populate_expset(expbase,expss,explabels, [Js,nstims,mindists,maxdists,trials])
    
    af.submenu_extractSpikes()
    #af.submenu_plot(6,'psth')
    #af.submenu_plot(7,'raster')
    #af.perform_analysis([])
    af.perform_analysis(['cv_isi','fano_factor_isi','mean_rate','isi','cv_kl'])                            
                            
                            
                            
        
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        pass
    elif sys.argv[1] == 'analyse':
        analyse()
    else:
        runon = sys.argv[1] 
        main(runon)                            
                            
                            
                            
                            
                            
                                