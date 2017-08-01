
from run_experiments import ExpSetup


expbase = '170730_supplFigures'

"""
Purpose:

Supplementary figure for paper
"""

es = ExpSetup()


tstop = 2000
light_dur = 1000
light_on = 500


L5PC_areas = ['soma', 'apic', 'dend'] # NB: we don't include 'axon'
Stellate_areas = ['soma', 'dendrite'] # NB: we don't include 'axon'

TOTAL_NUM_AREAS = {'L5PC': 3 ,
                   'SHStellate': 2}
areas = {'L5PC': L5PC_areas,
         'SHStellate': Stellate_areas}


def _run_simulation(celltype, irr, factor):
    """
    @params
        irr       irradiance
        factor    NpHR factor
    """

    pp = {}
    pp['cell'] = ['Neuron', celltype]
    pp['cell_params'] = {}
    pp['tstart'] = 0
    pp['tstop'] = tstop

    # opsin
    chrdict = {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur, 'lightdelay':light_on, 'interpulse_interval':0, 'n_pulses':1}
    hrdict = {'exp':5e-4, 'irradiance' :factor * irr, 'pulsewidth': light_dur, 'lightdelay':light_on, 'interpulse_interval':0, 'n_pulses':1}

    pp['opsindict'] = {}
    pp['opsindict']['ChR'] = {}
    for area in areas[celltype]:
        pp['opsindict']['ChR'][area] = chrdict
    pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}

    pp['opsindict']['NpHR'] = {}
    for area in areas[celltype]:
        pp['opsindict']['NpHR'][area] = hrdict
    pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}

    # general settings
    pp['experiment_type'] = 'opsinonly'
    pp['savedata'] = True



    # mark locations for PSP throughout dendritic trees, along with soma
    pp['mark_loc'] = {}
    pp['mark_loc']['names'] = ['mysoma'] # + labels
    pp['mark_loc']['sections'] = ['soma'] # + ['apic'] * num_apic + ['dend'] * num_basal
    pp['mark_loc']['ids'] = [(0, 0.5)] # + [('select_section_posn_bydistance', {'sectionarea':'apic', 'mindist':adist[0], 'maxdist':adist[1]}) for i in range(num_apic)] + [('select_section_posn_bydistance', {'sectionarea':'dend', 'mindist':bdist[0], 'maxdist':bdist[1]}) for i in range(num_basal)]


    # set up recordings
    # - record voltage from soma
    pp['record_loc'] = {}
    pp['record_loc']['v'] = ['mysoma']

    pp['plot'] = { 'voltages': [ ['mysoma', 'v', 'k-'] ]}

    pp.update({'expname':expbase,
                   'description':'_cell%s_irr%.3f_factor%.2f_irrOnly' % (celltype, irr, factor)})
    es.run_single_experiment(expbase, 'local', pp)



cells = ['L5PC']# SHStellate'] # L5PC'] #

factors = [0.125, 0.25, 0.5, 0.75, 1., 1.5, 2.] # for vitro
factors = [0., 0.0625, 0.125, 0.25, 0.5, 0.75, 1., 2., 4., 8.] # for vivo

irrs = [0.0001] # [1., 0.1, 0.01, 0.001]



for cell in cells:
    for irr in irrs:
        for factor in factors:
            _run_simulation(cell, irr, factor)







