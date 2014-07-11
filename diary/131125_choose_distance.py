import run_analysis
import sys
import Neuron 
import NeuroTools.stgen as ntst 
import run_stimulation
from run_experiments import ExpSetup
import numpy as np
from neuron import h

expbase = '131125_choose_distance'

fixedsubset = False #True
fixedtimes  = False #True # False# True #

explist = [['none','none']] 
tstop = 1200

es = ExpSetup()
areas = es.get_areas_section()
params = {}
Js = [3.] 
freqs = [50]
nstims = [80] #range(40,81,10)
num_trials = 1 #20
mindists = range(0,1200,100)
stg = ntst.StGen()

def select_sections(cell,subsection,distmin=0,distmax=-1):
    """

    
    """
    print 'Dmin and dmax are = ',distmin,distmax
    # for subsection
    sectionlist = cell.get_subdomain_list(subsection)
    distance_list = []
    for (i,section) in enumerate(sectionlist):
        
        d = h.distance(0.5,sec=section)
        #print 'distance of section %g is %g ...'%(i,d),
        if distmax == -1 and d>=distmin:
        #    print 'added'
            distance_list.append(i)
        elif d>=distmin and d<=distmax:
        #    print 'added'
            distance_list.append(i)
        #else:
        #    print 'NOT added'
    print 'Found %g sections that matched criteria (dmin=%g,dmax=%g)'%(len(distance_list),distmin,distmax)
    return distance_list

def get_distances(cell,subsection,nquantals=1):
    
    # for subsection
    sectionlist = cell.get_subdomain_list(subsection)
    quantals = get_quantals(nquantals)
    distance_list = [h.distance(q,sec=section) for section in sectionlist for q in quantals]
    return distance_list

def select_sites(cell,subsection,distmin=0,distmax=-1,nquantals=1):
    """

    
    """
    print 'Dmin and dmax are = ',distmin,distmax
    # for subsection
    sectionlist = cell.get_subdomain_list(subsection)
    quantals = get_quantals(nquantals)
    distance_list = []
    for (i,section) in enumerate(sectionlist):
        for q in quantals:
            d = h.distance(q,sec=section)
            #print 'distance of section %g is %g ...'%(i,d),
            if distmax == -1 and d>=distmin:
                #print 'added'
                distance_list.append((i,q))
            elif d>=distmin and d<=distmax:
                #print 'added'
                distance_list.append((i,q))
            #else:
                #print 'NOT added'
    print 'Found %g sections that matched criteria (dmin=%g,dmax=%g,num_quantals=%g)'%(len(distance_list),distmin,distmax,nquantals)
    return distance_list

    
def get_spiketimes(rate,num_inputs):
    spikes = []
    for i in range(num_inputs):
        ss = stg.poisson_generator(rate,t_stop=tstop,array=True)
        spikes.append(ss)
    #print 'Spikes = ',spikes
    return spikes

def add_locations(num_loc,upper_id,label='apic',quantas=[0.5]):
    labels = [label+str(i) for i in range(num_loc)]
    ids = range(upper_id)
    np.random.shuffle(ids)
    locs = ids[:num_loc]
    locs = [(l,q) for l in locs for q in quantas]
    return labels,locs

def convert_sites(sitelist,num_loc,shuffle=True,label='apic'):
    """
    sitelist    list of sites i.e. [(45,0.2),(47,0.5), ...]
                where each element is a tuple that contains the id and precise location
    """
    if shuffle:
        np.random.shuffle(sitelist)
    if num_loc > len(sitelist):
        print "convert_site:num_loc > number of available sites. Setting num_loc = %g"%(len(sitelist))
        num_loc = len(sitelist)
    
    labels = [label+str(i) for i in range(num_loc)]
    locs = sitelist[:num_loc]
    return labels,locs
    
def get_quantals(num_quantas):
    return [(i+1)*1./(num_quantas+1) for i in range(num_quantas)]    
    
def main():
    
    # set up and get sections
    cell = Neuron.L5PC()
    #subsetlocs = select_sections(cell,'apic',400,600)
    subsetlocs = select_sites(cell,'apic',distmin=600,nquantals=8)
    labels,locs = convert_sites(subsetlocs,80)
    
    
    count = 0
    for exparea in explist:
        for J in Js:
            for freq in freqs:
                for extraloc in nstims:
                    for dm in mindists:
                        subsetlocs = select_sites(cell,'apic',distmin=dm,nquantals=4)
                        labels,locs = convert_sites(subsetlocs,extraloc)
                         
                    #for trial in range(num_trials):
                    
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
                        
                        labels,locs = labels[:extraloc],locs[:extraloc]
                        times = get_spiketimes(freq,extraloc)

                            
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
                                       'description':'freq%g_J%.1f_nstim%g_mindist%g'%(freq,J,extraloc,dm)})
                        
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


def analyse_distances():
    import pylab
    # set up and get sections
    cell = Neuron.L5PC()
    #subsetlocs = select_sections(cell,'apic',400,600)
    for nq in [2**i for i in range(8)]:
        pylab.figure()
        distance_list = get_distances(cell,'apic',nquantals=nq)
        pylab.hist(distance_list, 50)
        pylab.title('Distance distribution, nq = %g'%nq)
        pylab.savefig('distance_apic_nq%g.png'%nq)

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
    elif sys.argv[1] == 'dist':
        analyse_distances()
    else:
        runon = sys.argv[1] 
        main()