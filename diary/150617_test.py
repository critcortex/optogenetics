"""
test.py
"""
import run_experiments
expbase = "150617_test"


def test_cell_vitro(cell, cellsections):
    es = run_experiments.ExpSetup()
    irr = 1.
    factor =0.5
    light_on = 100
    light_dur = 200
    tstop = 500
    nsites = 1
    Ia = 1.5
    dist_vitro = [[101,-1]]
    inj_locations = cellsections[0]
    pp = {}
    # neuron model and params
    pp['cell'] = ['Neuron', cell]
    pp['cell_params'] = {}
    
    # opsin
    chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
    hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':250,  'n_pulses':1}
    
    pp['opsindict'] = {}
    if irr > 0 :
        pp['opsindict']['ChR'] = {'soma': chrdict}
        for subdomain in cellsections:
            pp['opsindict']['ChR'][subdomain] = chrdict
        pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
    else:
        pp['ChR_areas'] = {'none'      : [None]}
        
    if irr > 0 and factor > 0:
        pp['opsindict']['NpHR'] =  {'soma': hrdict}
        for subdomain in cellsections:
            pp['opsindict']['NpHR'][subdomain] = hrdict
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
    pp['mark_loc']['names'] = ['mysoma']+['stimloc']
    pp['mark_loc']['sections'] = ['soma']+[inj_locations]
    pp['mark_loc']['ids'] = [(0,0.5)] + [('select_section_posn_bydistance',{'sectionarea':inj_locations,'mindist':dists[0],'maxdist':dists[1]}) for (i,dists) in enumerate(dist_vitro)] 
    #TODO: update
    
    pp['record_loc'] = {}
    pp['record_loc']['v'] = ['mysoma'] #+labels
    pp['record_loc']['ina'] = ['mysoma']
    pp['record_loc']['ik'] = ['mysoma']
    
    pp['stim_iclamp'] = True
    #print labels[0]
    #pp['spiketrains'] = [{'tstims': get_spiketimes(freq,nsites), 'locations': labels, 'weights':np.ones(nsites)*J,  'el': 0.1}]
    stimloc = 'stimloc'
    pp['iclamp'] = [{'tstim':0,  'location': stimloc, 'amp':Ia, 'duration':tstop}]
           
    
    
    vplots_soma = [['mysoma','v','k']]

    #iplots_soma = [['mysoma','ina','g'],['mysoma','ik','b']]
    pp['plot'] = {1:vplots_soma} 
        
    
    pp['num_threads'] = 1
                               
    pp.update({'expname':expbase,
               'description':'cell%s_irr%.3f_factor%.2f_I%.2f_stimloc_%s'%(cell,irr,factor,Ia,stimloc)})


    es.run_single_experiment(expbase, 'local', pp)
    
    
    
    
def test_stellate_vitro():
    test_cell_vitro('SHStellate',['dendrite'])

def test_L5PC_vitro():
    test_cell_vitro('L5PC',['apic','dend'])
    
    
    
    
    
    
#test_L5PC_vitro()
test_stellate_vitro()