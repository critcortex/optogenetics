
"""
Purpose: determine whether it's possible to get irradiance for opsin activation to vary with distance (mimicing depth)
 
Approaches:
1) add another pair of parameters to opsin dictionary for dropoff with unit distance, and set irradiance to reflect this modified value at depth z
-- requires: 
----- knowing depth of each segment
----- setting values for illumination in src/opsin.py

2) As above, but extending the opsin mod file so that there is an illumination factor (for irradiance penetration), instead of 
setting the illumination directly.
This will allow us to vary the illumination "master" signal, but still only have to calculate the irradiance seen from each segment once 

 

"""

import Neuron
from run_experiments import ExpSetup
import NeuroTools.stgen as ntst 
import numpy as np
import run_stimulation



expbase = '150316_testirradianceDepth'


stg = ntst.StGen()
es = ExpSetup()


factor = 0.7
irr = 0.02
n_pulses = 1
ipi = 65 # interpulse interval

tstop = 10
light_dur = 200
light_on = 150


L5PC_areas = ['soma', 'axon', 'apic','dend']

TOTAL_NUM_AREAS = {'L5PC': 4,
                   'SHStellate':2 }
areas = {'L5PC': L5PC_areas}

celltype = 'L5PC' 


def setup():
 
    pp = {}
    # neuron model and params
    pp['cell'] = ['Neuron',celltype]
    pp['cell_params'] = {}
    
    dampened_dict = {'irr_gradient':0.0007, 'irr_surface':1.0, 'projection': 'y'}
    # opsin
    chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':ipi,  'n_pulses':n_pulses}
    hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':ipi,  'n_pulses':n_pulses}
    hrdict.update(dampened_dict)
    print hrdict
    
    pp['opsindict'] = {}
    pp['opsindict']['ChR'] =  {}
    for area in areas[celltype]:
        pp['opsindict']['ChR'][area] = chrdict    
    pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
    
    pp['opsindict']['NpHR'] =  {}
    for area in areas[celltype]:
        pp['opsindict']['NpHR'][area] = hrdict
    pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
    
    # general settings 
    pp['experiment_type'] = 'opsinonly'
    pp['savedata'] = True # False #True
    
    pp['tstart'] = 0
    pp['tstop'] = tstop
    
    pp['num_threads'] = 1
                               
    pp.update({'expname':expbase,
               'description':'_testIrr'})

    es.run_single_experiment(expbase, 'local', pp)
    
    


def test_cell_morph():

    NE = run_stimulation.NeuronExperiment()
    NE._use_setup_params()
    NE.setup()
    cell = NE.cell
    
    """
    # following code taken from LFPy example
    from matplotlib.collections import PolyCollection
    import matplotlib.pyplot as plt
    zips = []
    for x, z in cell.get_idx_polygons_xy():
        zips.append(zip(x, z))
    
    polycol = PolyCollection(zips,
                             edgecolors='none',
                             facecolors='gray')
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.add_collection(polycol)
    ax.axis(ax.axis('equal'))
    
    plt.savefig('test.png')
    """
    


setup()
#test_cell_morph()


"""



"""















    