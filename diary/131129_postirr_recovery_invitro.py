import run_analysis
import sys
import Neuron 
import NeuroTools.stgen as ntst 
import run_stimulation
from run_experiments import ExpSetup
import numpy as np

expbase = '131129_postirr_recovery_invitro'


explist = [['whole','whole'],['whole','none'],['none','whole']] 
#explist = [['whole','none'],['none','whole']] 
#explist = [['whole','whole']]

es = ExpSetup()
areas = es.get_areas_section()
params = {}
Js = [3.] 
freqs = [50]
nstims = [80]

Is = np.arange(0.5,2.6,1.)
stim_locations = ['mysoma','proximal']
light_durations = np.array([5,10,20,50,100,200,500,1000])
post_irr_duration = 1000
pre_irr_duration = 600

factor = 1
irr = 2.

stg = ntst.StGen()

def get_spiketimes(rate,num_inputs,tstop):
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
        """
        for J in Js:
            for freq in freqs:
                for extraloc in nstims:
        """
        for Ia in Is:
            for stim_location in stim_locations:
                for light_dur in light_durations:
                    
                    tstop = pre_irr_duration + light_dur + post_irr_duration
                    
                    params['ChR_areas'] = {exparea[0]:areas[exparea[0]]}
                    params['NpHR_areas'] = {exparea[1]:areas[exparea[1]]}
                    params['num_threads'] = 1
                    
                    params['tstart'] = 0
                    params['tstop'] = tstop
                    params['stim_iclamp'] = True
                    params['iclamp'] = [{'amp':Ia,'tstim':0., 'duration':tstop,'location':stim_location}]
                    params['stim_epsp'] =  False
                    params['epsp'] =  [{'tstim':200, 'EPSPamp':10.1, 'location':'mysoma','risetau':0.5, 'decaytau':5., 'BACdt':0.}]
                    params['stim_spiketrains'] = False
                    #locations = ['mysoma','myapic']
                    
                    #labels,locs = add_locations(extraloc,109)
                    #params['spiketrains'] = [{'tstims': get_spiketimes(freq,extraloc,tstop),  'locations': labels, 'weights':np.ones(extraloc)*J,  'el': 0.02}]
                    
                    params['experiment_type'] = 'opsinonly'
                    params['savedata'] = True
                    
                    params['mark_loc'] = {}
                    params['mark_loc']['names'] = ['mysoma','myapic','proximal','distal']#+labels
                    params['mark_loc']['sections'] = ['soma','apic','apic','apic']#+['apic']*extraloc
                    params['mark_loc']['ids'] = [(0,0.5),(0,0.5),(0,0.0753),(0,0.972326)] #+ locs
                    
                    params['record_loc'] = {}
                    params['record_loc']['v'] = ['mysoma']
                    
                    params['plot'] = { 'voltages': [ ['mysoma', 'v','k-'] ]} 
                    chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':600,'interpulse_interval':250,  'n_pulses':1}
                    hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':600,'interpulse_interval':250,  'n_pulses':1}
                    
                    params['opsindict'] = {}
                    if exparea[0] != 'none':
                        params['opsindict']['ChR'] =  {'soma': chrdict,
                                                         'dend':chrdict,
                                                         'axon': chrdict,
                                                         'apic':chrdict}
                    if exparea[1] != 'none':
                        params['opsindict']['NpHR'] =  {'soma': hrdict,
                                                         'dend' : hrdict,
                                                         'axon':hrdict,
                                                         'apic':hrdict}
                                                   
                    
                    params.update({'expname':expbase,
                                   'description':'I%.1f_injloc%s_irrdur%g_ChR_%s_NpHR_%s'%(Ia,stim_location,light_dur,exparea[0],exparea[1])})
                    
                    params['cell'] = ['Neuron','L5PC']
                    params['cell_description'] = 'test'
                    
                    count +=1
                    print 'I%.1f_injloc%s_irrdur%g_ChR_%s_NpHR_%s'%(Ia,stim_location,light_dur,exparea[0],exparea[1])
                    
                    
                    if runon =='cluster':
                        print 'going to run on cluster'
                        es.run_single_experiment(expbase, runon, params)
                    elif runon == 'local':
                        pp = es.get_default_params()
                        NE = run_stimulation.NeuronExperiment()
                        pp.update(params)
                        NE.main(pp)
                        return
                    else:
                        print 'Invalid destination for job'
    print 'In total, submitted/ran %g jobs'%count


def analyse():
    af = run_analysis.AnalyticFrame()
    af.update_params({'tstart':100,'label_format':'dur=%g'})
        
    ia = Is[0]
    #exp_comp_list = [['_irr%.1f_factor%.2f_freq%g_J%1.f_nstim%g'%(irr,factor,freq,J,nstim),'ext_freq=%g'%(freq)] for freq in freqs ] #for J in Js for nstim in nstims]
    exp_comp_list = [['I%.1f_injloc%s_irrdur%s_ChR_%s_NpHR_%s_NpHR_%s_ChR_%s'%(ia,stimloc,'%g',exparea[0],exparea[1],exparea[1],exparea[0]),'ChR=%s,NpHR=%s,stimloc=%s'%(exparea[0],exparea[1],stimloc)] for exparea in explist for stimloc in stim_locations] #for J in Js for nstim in nstims]
    print exp_comp_list
    
    expss = [ec[0] for ec in exp_comp_list]
    explabels = [ec[1] for ec in exp_comp_list]
    print explabels
    af.populate_expset(expbase,expss,explabels, [light_durations])
    
    #af.submenu_extractSpikes()
    #af.submenu_plot(6,'')
    af.submenu_plot(7,'raster_normal',)
    af.submenu_plot(7,'raster_aligned',align=light_durations,offset=100-pre_irr_duration,tmax=800)
    af.submenu_plot(8,'voltage_irroff_aligned',align=light_durations,offset=100-pre_irr_duration,tmax=800)
    #af.perform_analysis([])
    #af.perform_analysis(['cv_isi','fano_factor_isi','mean_rate','isi','cv_kl'])
    
        
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        pass
    elif sys.argv[1] == 'analyse':
        analyse()
    else:
        runon = sys.argv[1] 
        main()