from neuron import h, run, init
import run_stimulation
from run_experiments import ExpSetup

h.load_file("nrn/hoc/neu_tree.hoc")

expname = '150126_electrotonic'

def test_treesneu():
    
    # default Tree
    nb,nc,nl = 4,2,5
    pp = {}
    pp['expname'] = '150126_electrotonic_test_treesneu'
        # neuron model and params
    pp['cell'] = ['Neuron', 'FractalNeuronRall']
    pp['cell_params'] = {'num_base': nb,
                  'num_split': nc,
                  'num_levels':nl,
                  'defaultlevel':{'mechanisms':['pas_dend','hh'],#,'ih','ca_lvast','ca_hva','sk3','ske2','na_tat','ca_dyn','im'],
                                  'cm':2.,
                                  'Ra':100},
                  'soma':{'mechanisms':['pas_soma','hh'],#'ih','ca_lvast','ca_hva','sk3','ske2','ktst','kpst','na_et2','na_tat','ca_dyn'],
                          'cm':1.},
                  'dend0':{'mechanisms':['pas_dend','hh'],
                           'cm':2.,
                           'Ra':100},
                  'dend1':{'mechanisms':['pas_dend','hh'],
                           'cm':2.,
                           'Ra':100},               
                  'mechanisms': {'pas_soma':('pas',{'e':-70}),#,'g_pas':1./20}),
                                 'pas_dend':('pas',{'e':-70})},#,,'g_pas':1./10}),}
                                 'hh':('hh',),
            }
    pp['opsindict'] = {}
    pp['record_loc'] = None
    pp['tstop'] = 0
    
    pp['savedata'] = False
    pp['output'] = {}
    pp['output']['neu'] = 'test2.neu'
    
    es = ExpSetup()
    es.run_single_experiment(expname, 'local', pp)
    
    #NE = run_stimulation.NeuronExperiment()
    #NE.params.update(pp)
    #NE.setup()
    #print("Got here fine")
    #h.neu_tree('test.neu')
    
def save_tree(tree):
    nb,nc,nl = tree
    pp = {}
    pp['expname'] = '150126_electrotonic'
        # neuron model and params
    pp['cell'] = ['Neuron', 'FractalNeuronRall']
    pp['cell_params'] = {'num_base': nb,
                  'num_split': nc,
                  'num_levels':nl,
                  'defaultlevel':{'mechanisms':['pas_dend','hh'],#,'ih','ca_lvast','ca_hva','sk3','ske2','na_tat','ca_dyn','im'],
                                  'cm':2.,
                                  'Ra':100},
                  'soma':{'mechanisms':['pas_soma','hh'],#'ih','ca_lvast','ca_hva','sk3','ske2','ktst','kpst','na_et2','na_tat','ca_dyn'],
                          'cm':1.},
                  'dend0':{'mechanisms':['pas_dend','hh'],
                           'cm':2.,
                           'Ra':100},
                  'dend1':{'mechanisms':['pas_dend','hh'],
                           'cm':2.,
                           'Ra':100},               
                  'mechanisms': {'pas_soma':('pas',{'e':-70}),#,'g_pas':1./20}),
                                 'pas_dend':('pas',{'e':-70})},#,,'g_pas':1./10}),}
                                 'hh':('hh',),
            }
    pp['opsindict'] = {}
    pp['record_loc'] = None
    pp['tstop'] = 0
    pp['description'] = 'tree_nb%g_nc%g_nl%g'%(nb,nc,nl)
    pp['savedata'] = False
    pp['output'] = {}
    pp['output']['neu'] = 'tree_nb%g_nc%g_nl%g.neu'%(nb,nc,nl)
    es = ExpSetup()
    es.run_single_experiment(expname, 'cluster', pp)
    
    NE = run_stimulation.NeuronExperiment()
    NE.main(pp)
    h.neu_tree('tree_nb%g_nc%g_nl%g.neu'%(nb,nc,nl))
    
     
def save_tree_L5PC():
    pp = {}
    pp['expname'] = '150126_electrotonic'
        # neuron model and params
    pp['cell'] = ['Neuron', 'L5PC']
    pp['opsindict'] = {}
    pp['record_loc'] = None
    pp['tstop'] = 0
    pp['description'] = 'tree_L5PC'
    pp['savedata'] = False
    pp['output'] = {}
    pp['output']['neu'] = 'tree_L5PC'
    es = ExpSetup()
    es.run_single_experiment(expname, 'local', pp)

def save_tree_SHStellate():
    pp = {}
    pp['expname'] = '150126_electrotonic'
        # neuron model and params
    pp['cell'] = ['Neuron', 'SHStellate']
    pp['opsindict'] = {}
    pp['record_loc'] = None
    pp['tstop'] = 0
    pp['description'] = 'tree_SHstellate'
    pp['savedata'] = False
    pp['output'] = {}
    pp['output']['neu'] = 'tree_SHstellate'
    es = ExpSetup()
    es.run_single_experiment(expname, 'local', pp)
 


def loop_trees():
    collection_same_total_all = [(1,1,124),(1,2,7), \
                                      (2,2,6),(2,7,3),(2,1,62),(2,61,2), \
                                      (4,5,3),(4,2,5),(4,30,2),(4,1,31), \
                                      (11,1,11),(11,10,2), \
                                      (18,1,7),(18,2,3),(18,6,2), \
                                      (31,1,4),(31,3,2), \
                                      (62,1,2),(124,1,1)]
    tt = {(1, 1, 124): (1, 1, 124),
 (1, 123, 2): (1, 123, 2),
 (2, 1, 62): (2, 1, 62),
 (2, 61, 2): (2, 61, 2),
 (3, 1, 41): (3, 1, 41),
 (3, 40, 2): (3, 40, 2),
 (4, 1, 31): (4, 1, 31),
 (4, 2, 5): (4, 2, 5),
 (4, 5, 3): (4, 5, 3),
 (4, 30, 2): (4, 30, 2),
 (5, 1, 25): (5, 1, 25),
 (5, 24, 2): (5, 24, 2),
 (25, 1, 5): (25, 1, 5),
 (25, 4, 2): (25, 4, 2),
 (31, 1, 4): (31, 1, 4),
 (31, 3, 2): (31, 3, 2),
 (41, 1, 3): (41, 1, 3),
 (41, 2, 2): (41, 2, 2),
 (62, 1, 2): (62, 1, 2),
 (124,1, 1): (124,1, 1)}
    collection_same_total_all = tt.keys()
    
    tt =  [(1, 1, 124),
 (1, 2, 7),
 (1, 3, 5),
 (1, 123, 2),
 (2, 1, 62),
 (2, 2, 6),
 (2, 61, 2),
 (3, 1, 41),
 (3, 1, 42),
 (3, 40, 2),
 (3, 41, 2),
 (4, 1, 31),
 (4, 2, 5),
 (4, 5, 3),
 (4, 30, 2),
 (5, 1, 25),
 (5, 24, 2),
 (6, 1, 21),
 (6, 4, 3),
 (6, 20, 2),
 (7, 1, 18),
 (7, 17, 2),
 (9, 1, 14),
 (9, 13, 2),
 (11, 1, 11),
 (11, 10, 2),
 (14, 1, 9),
 (14, 8, 2),
 (18, 1, 7),
 (18, 2, 3),
 (18, 6, 2),
 (21, 1, 6),
 (21, 5, 2),
 (25, 1, 5),
 (25, 4, 2),
 (31, 1, 4),
 (31, 3, 2),
 (41, 1, 3),
 (41, 2, 2),
 (42, 1, 3),
 (42, 2, 2),
 (62, 1, 2)]
    
    collection_same_total_all = tt
    
    i = 0
    for coll in collection_same_total_all:
        save_tree(coll)
        #i = i+1
        #if i>4: 
        #    break
 
 
def test_trees():
    h.load_file('sample.hoc')
    h.neu_tree('sample_out.neu')
    
   
#test_treesneu()
#loop_trees()
#test_trees()

#save_tree_L5PC()
#save_tree_SHStellate()

loop_trees()
