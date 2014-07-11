import run_analysis
import sys
import Neuron 
import NeuroTools.stgen as ntst 
import run_stimulation
from run_experiments import ExpSetup
import numpy as np


expbase = '131120_trial_variability'
expbase = '131122_fixedsubset_vartimes'
expbase = '131122_varsubset_fixedtimes'

fixedsubset = False #True
fixedtimes  = True # False# True #

explist = [['none','none']] 
tstop = 1200

es = ExpSetup()
areas = es.get_areas_section()
params = {}
Js = [3.] 
freqs = [50]
nstims = range(40,81,10)
num_trials = 20

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
    
def main():
    count = 0
    for exparea in explist:
        for J in Js:
            for freq in freqs:
                
                if fixedsubset:
                    subsetlabels,subsetlocs = add_locations(109,109) 
                if fixedtimes:
                    spiketimes = get_spiketimes(freq,109)
                
                for extraloc in nstims:
                    
                    for trial in range(num_trials):
                            
                        params['ChR_areas'] = {exparea[0]:areas[exparea[0]]}
                        params['NpHR_areas'] = {exparea[1]:areas[exparea[1]]}
                        params['num_threads'] = 1
                        
                        params['tstart'] = 0
                        params['tstop'] = tstop
                        params['stim_iclamp'] = False # True
                        params['iclamp'] = [{'amp':.5,'tstim':500., 'duration':500.,'location':'proximal'}]
                        params['stim_epsp'] =  False
                        params['epsp'] =  [{'tstim':200, 'EPSPamp':10.1, 'location':'mysoma','risetau':0.5, 'decaytau':5., 'BACdt':0.}]
                        params['stim_spiketrains'] = True
                        
                        
                        if not fixedsubset:
                            labels,locs = add_locations(extraloc,109)
                        else:
                            labels,locs = subsetlabels[:extraloc],subsetlocs[:extraloc]
                        if not fixedtimes:
                            times = get_spiketimes(freq,extraloc)
                        else:
                            times = spiketimes[:extraloc]
                            
                        params['spiketrains'] = [{'tstims': times,  'locations': labels, 'weights':np.ones(extraloc)*J,  'el': 0.02}]
                        
                        params['experiment_type'] = 'opsinonly'
                        params['savedata'] = True
                        
                        params['mark_loc'] = {}
                        params['mark_loc']['names'] = ['mysoma','myapic','proximal','distal']+labels
                        params['mark_loc']['sections'] = ['soma','apic','apic','apic']+['apic']*extraloc
                        params['mark_loc']['ids'] = [(0,0.5),(0,0.5),(0,0.0753),(0,0.972326)] + locs
                        
                        params['record_loc'] = {}
                        params['record_loc']['v'] = ['mysoma']
                        
                        params['plot'] = { 'voltages': [ ['mysoma', 'v','k-'] ]} 
                        params['opsindict'] = { }
                        
                        params.update({'expname':expbase,
                                       'description':'freq%g_J%.1f_nstim%g_trial%g'%(freq,J,extraloc,trial)})
                        
                        params['cell'] = ['Neuron','L5PC']
                        params['cell_description'] = 'test'
                        
                        count +=1
                        #print 'freq%g_J%.1f_nstim%g_trial%g'%(freq,J,extraloc,trial)
                        #print subsetlocs[:10],times[0][:10]
                        
                        if runon =='cluster':
                            print 'going to run on cluster'
                            es.run_single_experiment(expbase, runon, params)
                        elif runon == 'local':
                            pp = es.get_default_params()
                            NE = run_stimulation.NeuronExperiment()
                            pp.update(params)
                            NE.main(pp)
                            return
    print 'In total, submitted/ran %g jobs'%count


def analyse():
    factor = 1
    irr = 1.
    freq = freqs[0]
    J = Js[0] 
    af = run_analysis.AnalyticFrame()
    af.update_params({'tstart':100,'tstop':tstop,'label_format':'trial=%g'})
        
    exp_comp_list = [['freq%g_J%.1f_nstim%g_trial%s_NpHR_none_ChR_none'%(freq,J,nstim,'%g'),'nstim=%g'%(nstim)] for nstim in nstims ] 
    print exp_comp_list
    
    expss = [ec[0] for ec in exp_comp_list]
    explabels = [ec[1] for ec in exp_comp_list]
    print explabels
    
    af.populate_expset(expbase,expss,explabels, [range(num_trials)])
    """
    af.submenu_extractSpikes()
    af.submenu_plot(6,'psth')
    af.submenu_plot(7,'raster')
    af.submenu_plot(8,'popvoltage')
    """
    #af.perform_analysis(['cv_isi','fano_factor_isi','mean_rate','isi','cv_kl'])
    
    af.submenu_plot(9,'compare',traits=['cv_isi','fano_factor_isi','mean_rate','cv_kl'])

    
        
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        pass
    elif sys.argv[1] == 'analyse':
        analyse()
    else:
        runon = sys.argv[1] 
        main()