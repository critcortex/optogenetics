"""
run_experiments.py 
@author sjarvis
@date january 2013

@description:
helper class for running jobs locally or on cluster.

"""
from multiprocessing import Process, Queue, current_process, freeze_support
import copy

NUM_THREADS = 6 # number of threads/cores available for this machine

class Worker(Process):
"""
Worker object for running an individual instance of an experiment on a node

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
Class to create the queue for processing jobs
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
Wrapper class for running experiments/trials
"""

    def main_run_local(self,expname,params):
        """
        For a collection of experiments with varying parameters, sets up and 
        queues jobs onto a local queue, so that each job runs on a different node.
        This is highly practical for running NEURON simulations :)

        @params
		expname      string, name of experiment
                params       dict, collection of parameters
        """
        self.expengine = ExpEngine()
        print "Running for experiment", expname, 'locally'
        
        # this is a step I perform to setup my filesystem i.e. create folder for the files that will be generated
        # You may want/need to do some setup - so here's the place to do it
        #exps = self.generate_params(expname,params)
        
        if len(exps)>params['num_threads']:
            #TODO work out a better way of handling this. For moment, just run the jobs that we can
            print 'WARNING - you have more jobs than threads. Errors will probably occur'
            print '-------------> we will set the number of workers == number of threads. Eek'
            self.expengine.start_queue(len(exps))
        else:
            self.expengine.start_queue(params['num_threads'])
        
        for exp in exps:
            self.expengine.add_experiment(exp)
        self.expengine.end_engine()


class NeuronExperiment:
"""
The class that does the actual generation of the h. NEURON 
"""    
    
    def __init__(self):
        h.load_file('stdlib.hoc', 'String') 
        h.load_file('stdrun.hoc')
        h.load_file("import3d.hoc")
        #self.stimuli = {}
        #self.outputparams = self.get_default_outputparams()
        #self.params = self.get_default_params()
        
    def main(params):
	"""
 	The actual bulk of an experiment, where NEURON calls are run.
	
	In my code, I have a generic experiment and pass the values for different possibilities
        via a dictionary (params). 
	"""
        pass #as my code won't work for your setup, so I've taken it out
        """
	# In your example you might have the following example from test_classes.py


	h.celsius = 18.5
	h.tstop = 1e2 # set simulation duration (ms)
	h.dt = 0.0025 #set time step (ms)
	h.finitialize(-65) #initialize voltage state


	ax1 = Unmyelinated('1',100,3000,3,1,200)
	# default values are used
	ax1.channel_init(0.120,0.036,0.0003,50,-77,-54.3)
	stim = Stimulus(0,ax1.axon,0,1e2,10.0, 0.1,0.05)


	### Recording vectors ###

	vrec = [ h.Vector() for i in range(n)]
	irec = [ h.Vector() for i in range(n)]
	ikrec = [ h.Vector() for i in range(n)]
	inarec = [ h.Vector() for i in range(n)]
	trec = h.Vector()
	# record time variable
	trec.record(h._ref_t)
	# record voltage from n segments of the section
	for i in range(n):
	    vrec[i].record(ax1.axon(1.0/n*i+0.1)._ref_v)
	# record current from n segments of the section
	for i in range(n):
	    irec[i].record(ax1.axon(1.0/n*i+0.1)._ref_i_cap)
	    ikrec[i].record(ax1.axon(1.0/n*i+0.1)._ref_ik)
	    inarec[i].record(ax1.axon(1.0/n*i+0.1)._ref_ina)


	### LAUNCH SIMULATION ###
	h.run()

	"""


def TestExperiment():
""" Example way that we would use 
