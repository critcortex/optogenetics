# -*- coding: utf-8 -*-
#import run_stimulation
import os, sys
import time 
import logging
import copy
from subprocess import call
import numpy as np
import collections
from multiprocessing import Process, Queue, current_process, freeze_support

"""
ExpEngine : class for modules to handle the running of multiple jobs



"""


class NeuronExperiment():

    def main(self, experiment_dict):
        # This is where you would set up the experiment in NEURON i.e. 
        # h.init() 
        # Note that you should put this class in a separate file
        print("Running an experiment with parameters: ----------------")
        print(experiment_dict['description'])
        print("-------------------------------------------------------")
        print(experiment_dict)
        print('=======================================================')



class Worker(Process):
    """
    The Worker, which is a job that runs on a single thread. When the job runs, it executes
    the experiment by sending the parameters it received from the experiment to the pyNEURON code
    that actually creates the code i.e. class NeuronExperiment
    """
    def __init__(self, queue):
        super(Worker, self).__init__()
        self.queue= queue

    def run(self):
        print 'Worker started with PID=%g'%os.getpid()
        for data in iter( self.queue.get, None ):
            NE = NeuronExperiment()
            NE.main(data)
            
class ExpEngine():
    """
    The Queue engine that handles the workers. This was taken from online examples on the Python website
    """

    def __init__(self):
        print 'Have created queue'
        self.queue = Queue()
    
    def add_experiment(self,expdict):
        """ NB that we have to use a deep copy else we get errors when passing multiple dicts """
        expcopy = copy.deepcopy(expdict)
        self.queue.put(expcopy)  
       
            
    def start_queue(self,num_threads=4):
        '''start some threads, each one will process one job from the queue'''
        self.num_threads = num_threads
        for i in range(self.num_threads):
            Worker( self.queue ).start()

    def end_engine(self):
        # Tell child processes to stop
        for i in range(self.num_threads):
            self.queue.put(None)


        

class ExpSetup():
    """
    A wrapper class that takes usually handles jobs that are submitted to the cluster as well as my local machine. It 
    also handles my default parameters that I use in my experiments.  
    """

    
    def get_default_params(self):
        defaultp= {  'expname'           : None,
                'description'       : '',
                'experiment_type'   : '',
                # neuron params
                'cell'              : None,
                'cell_description'  : '',
                'cell_params'       : {},
                
                'stim_iclamp'       : False,
                'iclamp'            : [{'amp':1.,
                                        'tstim':200.,
                                        'duration':10.,
                                        'location':'soma'},],
                'stim_epsp'         : False,
                'epsp'              : [{'tstim':200,
                                        'EPSPamp':0.1,
                                        'risetau':0.5,
                                        'decaytau':5.,
                                        'BACdt': 0.,
                                        'location':'soma'}],
                'stim_spiketrains'  : False,
                'spiketrains'       : [{'tstims': [[200, 600]], 
                                       'locations': ['soma'], 
                                       'weights':[1.], 
                                       'el': 0.02}],
                'mark_loc'          : {'names':[],'sections':[], 'distances':[],'ids':[]},
                'record_loc'        : {'v': [],
                                       'i': [],
                                       'na': [],
                                       'k':[],
                                       'g':[],
                                       'i_ChR':[],
                                       'i_NpHR':[]},
                'plot'              : {},
                'output'            : {'output_pkl': False}
                
              }
        
        defaultp['savedata'] =True
        defaultp['logdata'] = True
        defaultp['tstart'] = 0 
        defaultp['tstop'] = 4000
        defaultp['num_threads'] = 1
        return defaultp
        
    
    def run_single_experiment(self,expbase,runon='local',params={},checkdatexists=False):
        """
        
        """
        if checkdatexists:
            if self.__check_dat_file_exists(expbase,params):
                # if it already exists, we don't need to run
                #print '---------------------------- would not run'
                return 
        #print 'would run ', expbase+params['description']
        self.main(expbase,params,runon=runon)                    
    
        
        
    def dict_update(self,d, u):
        """
        Need to be able to merge dictionaries of dictionaries (params)
        
        Source: http://stackoverflow.com/questions/3232943/
        """
        for k, v in u.iteritems():
            if isinstance(v, collections.Mapping):
                r = self.dict_update(d.get(k, {}), v)
                d[k] = r
            else:
                d[k] = u[k]
        return d

    def main_run_local(self,expname,expparams):
        """
        
        """
        self.expengine = ExpEngine()
        print "Running for experiment", expname, 'locally'
        # NB: I had some code in here that generated it from multiple experiments. 
        # As this is a skeleton script, I've put it in the easiest form i.e. a list that contains a single parameter and we know 
        # that we're only going to have one  (because we're calling this from run_single_experiment
        # In reality, you'll probably have multiple parameter sets - so go through and define another function which calls
        # this method and passes multiple parameter sets
        """
        if len(expparams)>params['num_threads']:
            print 'WARNING - you have more jobs than threads. Errors will probably occur'
            print '-------------> we will set the number of workers == number of threads. Eek'
            self.expengine.start_queue(len(expparams))
        else:
            self.expengine.start_queue(params['num_threads']) 
        """
        # we set how many threads in the engine
        self.expengine.start_queue(expparams['num_threads']) 
        # add an experiment to the job
        self.expengine.add_experiment(expparams)
        # and let the engine know we're finished (like closing a file)
        self.expengine.end_engine()

        
        
    def main(self,expname,newparams,runon='local'):
        """
        
        """
        params = self.get_default_params()
        newparams = self.dict_update(params, newparams)
        if runon=='local':
            self.main_run_local(expname,newparams)
        # As you don't have any clusters available, I've taken this code out
        #elif runon=='cluster' :
        #    print 'submitting to cluster'
        #    self.main_run_cluster(expname,newparams,runjob=(runon=='cluster'))
        else:
            raise Exception('Unknown option for runon: %s'%runon)


class SampleExperiment():
    """
    A sample experiment, in which I'd go through and define a new dictionary. Note that I usually submit multiple experiments to a cluster
    that runs on my machine; here, I've only submitted one job to the local queue. 
    """
    def __init__(self):
        self.es = ExpSetup()

    def test_sample_exp(self):
        pp = {  'expbase': '150622_test',
              'description':'testing_workerProcess_singleExp'}
        self.es.run_single_experiment(pp['expbase'], 'local', pp)
        
    def test_multiple_exps(self):
        for i in range(4):
            pp = {  'expbase': '150622_test',
                  'description':'testing_workerProcess_number%g'%i}
            self.es.run_single_experiment(pp['expbase'], 'local', pp)


sampleExp = SampleExperiment()
sampleExp.test_sample_exp()
sampleExp.test_multiple_exps()


