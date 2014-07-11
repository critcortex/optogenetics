import run_stimulation
import Neuron
from neuron import h
import pylab
import numpy as np
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 


es = ExpSetup()

stg = ntst.StGen()

tstop = 1000

def get_spiketimes(rate,num_inputs):
    spikes = []
    for i in range(num_inputs):
        ss = stg.poisson_generator(rate,t_stop=tstop,array=True)
        spikes.append(ss)
    #print 'Spikes = ',spikes
    return spikes


def basic():
    h.load_file('stdlib.hoc', 'String') 
    h.load_file('stdrun.hoc')
    h.load_file("import3d.hoc")
    
    cell = Neuron.SHStellate()
    #cell = Neuron.L5PC()
    #cell = Neuron.FractalNeuron()
    
    cell.print_section_info()
    
    h.refrac_SpikeOut = 5
    h.vrefrac_SpikeOut = -60
    h.thresh_SpikeOut = -50
    
    #soma = cell.get_soma()
    #spkout = h.SpikeOut(0.5,sec=soma)
    
    #clamp_site = cell.select_section_posn_bydistance('apic')[0]
    # for dendrites @ hill
    clamp_site = cell.select_section_posn_bydistance('dendrite')[2]
    # for soma
    #clamp_site = cell.select_section_posn_bydistance('soma')[0]
    stim = h.IClamp()
    stim.loc(clamp_site[2], sec=clamp_site[1])
    stim.dur = 400
    setattr(stim, 'del',200)
    stim.amp = 15.
    
    soma = clamp_site[1]
    
    time = h.Vector()
    vmsoma = h.Vector()
    
    time.record(h._ref_t)
    vmsoma.record (soma(0.5)._ref_v)
    
    
    h.init()
    h.tstop = 700
    #h.v_init = -62.4
    h.run()

    
    vms = np.array(vmsoma)
    tt = np.array(time)
    pylab.plot(tt,vms)
    pylab.savefig('test_iclamp%g.png'%stim.amp)


def main(runon):
    
    nsites = 10
    freq = 80
    
    tstop = 1000
    light_on = 200
    pw = 500
    irr = 0.1
    factor = 0.01
    dist = [10,200]
    J = 10.
    trial = 0
    
    
    expbase = '140408_SH_stellate'
    
    params = {}
    params['num_threads'] = 1
    
    params['tstart'] = 0
    params['tstop'] = tstop
    
    params['stim_spiketrains'] = True                            
    labels = ['stim%g'%i for i in range(nsites)]
    params['spiketrains'] = [{'tstims': get_spiketimes(freq,nsites),  'locations': labels, 'weights':np.ones(nsites)*J,  'el': 0.02}]
    
    
    params['experiment_type'] = 'opsinonly'
    params['savedata'] = True
    
    params['mark_loc'] = {}
    params['mark_loc']['names'] = ['mysoma']+labels
    params['mark_loc']['sections'] = ['soma']+['dendrite']*nsites
    params['mark_loc']['ids'] = [(0,0.5),(0,0.5),(0,0.0753),(0,0.972326)] + [('select_section_posn_bydistance',{'sectionarea':'dendrite','mindist':dist[0],'maxdist':dist[1]}) for i in range(nsites)] 
    
    
    params['record_loc'] = {}
    params['record_loc']['v'] = ['mysoma']+labels[:5]
    
    vplots_soma = [['mysoma','v','k']]
    params['plot'] = {1:vplots_soma}
    
    chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': pw,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
    hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': pw,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
    
    params['opsindict'] = { 'ChR' : {'soma': chrdict,
                                     'dendrite':chrdict,
                                     'axon': chrdict},
                            'NpHR' : {'soma': hrdict,
                                     'dendrite' : hrdict,
                                     'axon':hrdict}}
                                    
    params['NpHR_areas'] = {'whole'      : params['opsindict']['NpHR'].keys()}
    params['ChR_areas'] = {'whole'      :  params['opsindict']['ChR'].keys()} 
    params.update({
                   'expname':expbase,
                   'description':'_irr%.1f_factor%.2f_freq%g_J%.1f_nstim%g_mindist%g_maxdist%g_trial%g'%(irr,factor,freq,J,nsites,dist[0],dist[1],trial)})
    
    
    params['cell'] = ['Neuron','SHStellate']
    params['cell_description'] = 'test'
    
    
    
    if runon =='cluster':
        es.run_single_experiment(expbase, runon, params)
    elif runon == 'local':
        pp = es.get_default_params()
        NE = run_stimulation.NeuronExperiment()
        pp.update(params)
        NE.main(pp)
        
    
    
#basic()
main('local')
#main('cluster')