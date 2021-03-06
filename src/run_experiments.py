# -*- coding: utf-8 -*-
# import run_stimulation
import run_stimulation
import file_io as fio
import os, sys
import time
import logging
import copy
from subprocess import call
import numpy as np
import collections
from multiprocessing import Process, Queue, current_process, freeze_support

from celery import Celery, Task

app = Celery('run_experiments', broker='redis://localhost:6379/0')


"""
import logging
"""
"""
#Bit of code for logging to stdout rather than to default log file
#Taken from stackoverflow qn 14058453
root = logging.getLogger()
root.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
root.addHandler(ch)
# /logging
"""

OPSIN_EXP = 5e-4
OPSIN_IRR = 10
OPSIN_EXPDIST = 600

# for expressions, array is [expression, pd, del_light, t_on, t_off, distance*,num_pulses]
# *NB that distance is only relevant for apical_distal or apical_proximal
sampleexpression = [5e-4, 10, 600, 200, 350, 600, 1]
LOCATION_ALL = ['soma', 'basal', 'axon', 'apical_proximal', 'apical_distal']
samplevalues = [[x, sampleexpression] for x in LOCATION_ALL]

areasNone = {'None' : [None] }
areasSoma = {'soma' : ['soma']}
areasAll = {'whole' : LOCATION_ALL}


areaSections = {    'apical'    : ['apical_distal', 'apical_proximal'],
                    'basal'     : ['basal'],
                    'basalsoma' : ['basal', 'soma'],
                    'soma'      : ['soma'],
                    'axon'      : ['axon'] ,
                    'none'      : [None],
                    'whole'     : LOCATION_ALL }

SECTIONS_DICT = areaSections

areasPaper = {   'apical': ['apical_proximal', 'apical_distal'],
            'whole' : LOCATION_ALL,
            'soma'  : ['soma'],
            'axon'  : ['axon'] }



both_expression = [5e-4, 10, 500, 500, 350, 600, 1]
first_expression = [5e-4, 10, 200, 50, 400, 600, 10]
second_expression = [5e-4, 10, 400, 50, 400, 600, 10]
# both_expression = [5e-4,10,2500,3000,350,600]


areasNp = areasAll # dict(areasPaper.items() + areasNone.items()) #  areasSoma #areasPaper  # areasNone #
areasCh = areasAll # dict(areasPaper.items() + areasNone.items()) # areasPaper + areasNone # areasNone #
ChR_expression = first_expression
NpHR_expression = second_expression

experiment_type = 'opsinonly' # 'BAP' # 'BAC' #

class NeuronExperimentTask(Task):
    ignore_result = True

    def run(self, source, *args, **kwargs):
        self.source = source
        params = kwargs['params']
        print "CHECKING RUN:", params['expname'], params['description'], os.getcwd()
        NE = run_stimulation.NeuronExperiment()
        NE.main(params)



class Worker(Process):
    def __init__(self, queue):
        super(Worker, self).__init__()
        self.queue = queue

    def run(self):
        print 'Worker started with PID=%g' % os.getpid()
        for data in iter(self.queue.get, None):
            # sleeptime = data['sleeptime']
            # print '-'*20,'running ', data['expname'], 'is going to sleep for ', sleeptime, '..............', os.getpid()
            # time.sleep(sleeptime)
            NE = run_stimulation.NeuronExperiment()
            NE.main(data)

class ExpEngine():

    def __init__(self):
        print 'Have created queue'
        self.queue = Queue()

    def add_experiment(self, expdict):
        """ NB that we have to use a deep copy else we get errors when passing multiple dicts """
        expcopy = copy.deepcopy(expdict)
        self.queue.put(expcopy)


    def start_queue(self, num_threads=4):
        '''start some threads, each one will process one job from the queue'''
        self.num_threads = num_threads
        for i in range(self.num_threads):
            Worker(self.queue).start()

    def end_engine(self):
        # Tell child processes to stop
        for i in range(self.num_threads):
            self.queue.put(None)


class ExpPostal:


    def send_to_cluster(self, params):
        # TODO: send to scheduler
        pass

    """
    def processor(self):

        for job in iter(self.queue.get,'STOP'):
            print ">> I'm operating on job item: %s"%(job['expname'])
            sj = SampleJob()
            sj.run_job(job)
            #NE = run_stimulation.NeuronExperiment()
            #NE.main(job)
        
        while True:
            if self.queue.empty() == True:
                print "the Queue is empty!"
                sys.exit(1)
            try:
                job = self.queue.get()
                #print ">> I'm operating on job item: %s"%(job['expname'])
                #NE = run_stimulation.NeuronExperiment()
                #NE.main(job)
                sj = SampleJob()
                sj.run_job(job)
                #print ">> I'm finished on job item: %s"%(job['expname'])
                self.queue.task_done()
                #print 'Number of tasks still to go is -->',self.queue.qsize()
            except:
                print "Failed to operate on job"
        """



"""
#----------------------------------------------
class SampleJob():

    def run_job(self,params):
        #sleeptime = int(np.random.randint(1,6))
        sleeptime = params['sleeptime']
        areasChR = [x[0] for x in params['opdict']['ChR']] 
        areasHR = [x[0] for x in params['opdict']['NpHR']]
        #print areasHR, '\t\t', areasChR
        #print 'running ', params['expname'], 'is going to sleep for ', sleeptime
        time.sleep(sleeptime)
        #return

#----------------------------------------------

"""
class ExpSetup():

    def get_areas_section(self):
        return areaSections

    def get_areas_paper(self):
        return areasPaper


    def get_default_params(self):
        defaultp = {  'expname'           : None,
                'description'       : '',
                'experiment_type'   : experiment_type,
                # neuron params
                'cell'              : None,
                'cell_description'  : '',
                'cell_params'       : {},
                # optogen params
                'opsindict'         : {},
                'NpHR_areas'        : {},# areasNp,
                'ChR_areas'         : {},# areasCh,
                'NpHR_times'        : NpHR_expression,
                'ChR_times'         : ChR_expression,
                'description'       : '',
                # old params that are kept in for compatibility
                'soma_stim_DC'      : 0., # 1.,
                'iclamp_amp'        : 0., # 1.,
                'iclamp_start'      : 0., # 100.,
                'iclamp_duration'   : 0., # 100.,
                'iclamp_dist_amp'        : 0., # 1.,
                'iclamp_dist_start'      : 0., # 100.,
                'iclamp_dist_duration'   : 0., # 100.,
                'EPSP_amp'           : 0.0,
                'EPSP_transient'    : 0., # 200.,
                # and their new, improved counterparts
                'stim_iclamp'       : False,
                'iclamp'            : [{'amp':1.,
                                        'tstim':200.,
                                        'duration':10.,
                                        'location':'soma'}, ],
                'stim_epsp'         : False,
                'epsp'              : [{'tstim':200,
                                        'EPSPamp':0.1,
                                        'risetau':0.5,
                                        'decaytau':5.,
                                        'BACdt': 0.,
                                        'location':'soma'}],
                'stim_spiketrains'  : False,
                'spiketrains'       : [{'tstims': [[200, 600], [400, 600]],
                                       'locations': [['soma', 0.5]],
                                       'weights':[1., 1.],
                                       'el': 0.02}],
                'mark_loc'          : {'names':[], 'sections':[], 'distances':[], 'ids':[]},
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

        defaultp['savedata'] = True
        defaultp['logdata'] = True
        defaultp['tstart'] = 0
        defaultp['tstop'] = 4000
        defaultp['num_threads'] = 4
        return defaultp



    def populate_opsin_dict(self, expdict):
        # TODO: populate
        return expdict


    def generate_params(self, expname, expdict, runon='saveonly'):

        paramlist = []
        expdict['expname_family'] = expname

        if len(expdict['NpHR_areas'].keys()) > 0 and len(expdict['ChR_areas'].keys()) > 0:

            for (i, khr) in enumerate(expdict['NpHR_areas'].keys()): # loop for NpHR expression

                for (j, kch) in enumerate(expdict['ChR_areas'].keys()):


                    ChRlocations = expdict['ChR_areas'][kch]
                    ChR_expression = expdict['ChR_times']
                    ChRvalues = [[x, ChR_expression] for x in ChRlocations]
                    ChR_descript = "ChR_" + kch

                    NpHRlocations = expdict['NpHR_areas'][khr]
                    NpHR_expression = expdict['NpHR_times']
                    NpHRvalues = [[x, NpHR_expression] for x in NpHRlocations]
                    NpHR_descript = "NpHR_" + khr

                    opdict = {  'ChR': ChRvalues ,
                                'NpHR': NpHRvalues}

                    # expdict['expname']=expname+expdict['description']+"_%g_"%i+NpHR_descript+"_%g_"%j+ChR_descript
                    expdict['expname'] = expname + expdict['description'] + '_%s_%s' % (NpHR_descript, ChR_descript)
                    expdict['opdict'] = opdict

                    if expdict['logdata']:
                        fio.log_message('About to run an experiment name:%s\nParams:' % expname)
                        # for p in sorted(expdict.keys(),key=str.lower):
                        #    fio.log_message("Param: %s = \t%s"%(p,expdict[p]),timenow=False)

                    # print expdict['expname']
                    paramlist.append(copy.deepcopy(expdict))

        else:
            expdict['expname'] = expname + expdict['description']
            expdict['opdict'] = {}

            if expdict['logdata']:
                fio.log_message('About to run an experiment name:%s\nParams:' % expname)
                # for p in sorted(expdict.keys(),key=str.lower):
                #    fio.log_message("Param: %s = \t%s"%(p,expdict[p]),timenow=False)

            paramlist.append(copy.deepcopy(expdict))


        return paramlist


    def __check_dat_file_exists(self, expbase, params, fullloc=True):
        """
        
        """
        dirloc = fio.get_exp_dat_location(expbase)
        filename = expbase + params['description'] + '_NpHR_%s_ChR_%s' % (params['NpHR_areas'].keys()[0], params['ChR_areas'].keys()[0]) + '.dat'
        # print dirloc+'/'+filename
        return os.path.isfile(dirloc + '/' + filename)

    def run_single_experiment(self, expbase, runon='local', params={}, checkdatexists=False):
        """
        
        """
        if checkdatexists:
            if self.__check_dat_file_exists(expbase, params):
                # if it already exists, we don't need to run
                # print '---------------------------- would not run'
                return
        # print 'would run ', expbase+params['description']
        self.main(expbase, params, runon=runon)

    def run_frequency_range(self, expbase, freqs=[0.5, 1, 2.5, 5, 10], pulsewidth=10., transient=200., n_pulses=10, runon='local', params={}, offset_phase=0.5):
        """
        
        Params:
            freqs         list of frequencies (Hz)
            pulsewidth    either list (same length as freqs) or scalar, to indicate pulse width (ms)
            transient     transient before pulses start (ms)
        """
        conversion = 1000.
        for f in freqs:
            interstim_interval = 1. / f * conversion # convert freq from /s to /ms
            offset = interstim_interval * offset_phase
            t_off = interstim_interval - pulsewidth
            ChR_expression = [OPSIN_EXP, OPSIN_IRR, transient, pulsewidth, t_off, OPSIN_EXPDIST, n_pulses]
            NpHR_expression = [OPSIN_EXP, OPSIN_IRR, transient + offset, pulsewidth, t_off, OPSIN_EXPDIST, n_pulses]
            params.update({'NpHR_times':NpHR_expression, 'ChR_times':ChR_expression, 'description':'_freq%g_pw%g' % (f, pulsewidth)})
            self.main(expbase, params, runon=runon)


    def run_irraidiance_range(self, expbase, irr_range=[1, 2, 5, 10, 20, 50, 100], pulsewidth=10., transient=200., n_pulses=10, freq=2., runon='local', params={}):
        """
        
        Params:
            irr_range     list of irradiances (mW.mm-2)
            pulsewidth    either list (same length as freqs) or scalar, to indicate pulse width (ms)
            transient     transient before pulses start (ms)
        """
        # TODO: extend to allow different irr values for different opsin
        freq_conversion = 1000.
        for irr in irr_range:
            interstim_interval = 1. / freq * freq_conversion # convert freq from /s to /ms
            offset = interstim_interval / 2
            t_off = interstim_interval - pulsewidth
            ChR_expression = [OPSIN_EXP, irr, transient, pulsewidth, t_off, OPSIN_EXPDIST, n_pulses]
            NpHR_expression = [OPSIN_EXP, irr, transient + offset, pulsewidth, t_off, OPSIN_EXPDIST, n_pulses]
            params.update({'NpHR_times':NpHR_expression, 'ChR_times':ChR_expression, 'description':'_irr%g_freq%g_pw%g' % (irr, freq, pulsewidth)})
            self.main(expbase, params, runon=runon)

    def run_test_exp(self, expbase, runon='local', params={}):
        """
        A simple self-contained test case that can be used without having to set up too many params.
        Should eventually be migrated to a JUnit-like testcase
        """
        # TODO: move run_test_exp to testcase
        transient = 100
        pulsewidth = 50
        t_off = 1
        n_pulses = 1
        offset = 0
        ChR_expression = [OPSIN_EXP, OPSIN_IRR, transient, pulsewidth, t_off, OPSIN_EXPDIST, n_pulses]
        NpHR_expression = [OPSIN_EXP, OPSIN_IRR, transient + offset, pulsewidth, t_off, OPSIN_EXPDIST, n_pulses]

        params.update({'NpHR_times':NpHR_expression,
                       'NpHR_areas':  dict(areasAll.items()),# + areasNone.items()),
                       'ChR_times':ChR_expression,
                       'ChR_areas': dict(areasSoma.items() + areasAll.items() + areasNone.items()), #
                       'description':'_testexperiment'})

        self.main(expbase, params, runon=runon)



    def submit_to_cluster(self, expname):
        print "Submitting to cluster ...",
        jdffile = fio.write_jdf_stimulation(expname)
        call(['qsub', jdffile])
        print "Submitted to cluster"


    def dict_update(self, d, u):
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

    def main_run_local(self, expname, params):
        """
        
        """
        self.expengine = ExpEngine()
        print "Running for experiment", expname, 'locally'

        fio.setup_experiment(expname)
        exps = self.generate_params(expname, params)
        # TODO work out a better way of handling this
        if len(exps) > params['num_threads']:
            print 'WARNING - you have more jobs than threads. Errors will probably occur'
            print '-------------> we will set the number of workers == number of threads. Eek'
            self.expengine.start_queue(len(exps))
        else:
            self.expengine.start_queue(params['num_threads'])

        for exp in exps:
            self.expengine.add_experiment(exp)
        self.expengine.end_engine()
        # fio.finish_experiment(expname)

    def main_run_cluster(self, expname, params, runjob=True):
        """
        
        """
        print "Running for experiment", expname, 'on cluster'

        fio.setup_experiment(expname)
        exps = self.generate_params(expname, params)

        for exp in exps:
            print 'Am here'
            fio.savejob(exp, exp['expname'])

            if runjob:
                self.submit_to_cluster(exp['expname'])


    def main_run_missing(self, expname, newparams):
        # check whether job exists in main location
        exps = self.generate_params(expname, newparams)[0]
        testname = exps['expname']

        if fio.check_mainexp_location(testname + '*.dat'):
            print 'File exists A'
            return
        # check whether it exists in exp home
        if fio.check_expfolder(expname, testname, 'dat'):
            print 'File exists B'
            return
        # if still clear it doesn't exist, then run it
        print 'I want to run for file = ', expname, newparams['description']
        self.main_run_cluster(expname, newparams)


    def main_run_celery(self, expname, params):
        print "Running for experiment", expname, 'on celery'

        fio.setup_experiment(expname)
        exps = self.generate_params(expname, params)
        net = NeuronExperimentTask()
        for exp in exps:
            print "Submitting job to celery: %s" % exp['expname']
            print "pre Celery", exp
            net.delay(source=1, params=exp)



    def main_run_names(self, expname, params, runjob=True):
        """
        Shell method, so that we can see the names of files that would be run.
        Used for checking that loops for batch processing are set up correctly
        """
        print "Dry run: ------------------------------------"
        print "Running for experiment", expname, 'on cluster'
        print "                     :", params['description']


    def main(self, expname, newparams, runon='local'):
        """
        
        """
        params = self.get_default_params()
        newparams = self.dict_update(params, newparams)
        if runon == 'local':
            self.main_run_local(expname, newparams)
        elif runon == 'cluster' or runon == 'saveonly':
            print 'submitting to cluster'
            self.main_run_cluster(expname, newparams, runjob=(runon == 'cluster'))
        elif runon == 'celery':
            print 'submitting to celery'
            self.main_run_celery(expname, newparams)
        elif runon == 'missing':
            self.main_run_missing(expname, newparams)
        elif runon == 'names':
            self.main_run_names(expname, params)
        else:
            raise Exception('Unknown option for runon: %s' % runon)



if __name__ == '__main__':
    try:
        expname = sys.argv[1]
    except:
        expname = raw_input("No experiment name specified. Please enter: ")
    ExpSetup.main(expname)


