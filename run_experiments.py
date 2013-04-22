# -*- coding: utf-8 -*-

import run_stimulation
import file_io as fio
import os, sys
import time 
import logging
import copy
from subprocess import call
import numpy as np

from multiprocessing import Process, Queue, current_process, freeze_support


OPSIN_EXP = 5e-4
OPSIN_IRR = 10
OPSIN_EXPDIST = 600

# for expressions, array is [expression, pd, del_light, t_on, t_off, distance*]
# *NB that distance is only relevant for apical_distal or apical_proximal
sampleexpression = [5e-4,10,600,200,350,600]
LOCATION_ALL = ['soma','basal','axon','apical_proximal','apical_distal']
samplevalues = [[x,sampleexpression] for x in LOCATION_ALL]

areasNone= {'None' : [None] }             
areasSoma = {'soma' : ['soma']}
areasAll = {'whole' : LOCATION_ALL}


areaSections = {    'apical'    : ['apical_distal','apical_proximal'],
                    'basal'     : ['basal'],
                    'basalsoma' : ['basal','soma'],
                    'soma'      : ['soma'],
                    'axon'      : ['axon'] ,
                    'none'      : [None],
                    'whole'     : LOCATION_ALL }

SECTIONS_DICT = areaSections

areasPaper = {   'apical': ['apical_proximal','apical_distal'],
            'whole' : LOCATION_ALL,
            'soma'  : ['soma'],
            'axon'  : ['axon'] }



both_expression = [5e-4,10,500,500,350,600,1]
first_expression= [5e-4,10,200,50,400,600,10]
second_expression = [5e-4,10,400,50,400,600,10]
#both_expression = [5e-4,10,2500,3000,350,600]


areasNp = areasAll # dict(areasPaper.items() + areasNone.items()) #  areasSoma #areasPaper  # areasNone # 
areasCh = areasAll # dict(areasPaper.items() + areasNone.items()) # areasPaper + areasNone # areasNone # 
ChR_expression = first_expression
NpHR_expression = second_expression

experiment_type = 'opsinonly' # 'BAP' # 'BAC' #



class Worker(Process):
    def __init__(self, queue):
        super(Worker, self).__init__()
        self.queue= queue

    def run(self):
        print 'Worker started with PID=%g'%os.getpid()
        for data in iter( self.queue.get, None ):
            #sleeptime = data['sleeptime']
            #print '-'*20,'running ', data['expname'], 'is going to sleep for ', sleeptime, '..............', os.getpid()
            #time.sleep(sleeptime)
            NE = run_stimulation.NeuronExperiment()
            NE.main(data)
            
class ExpEngine():
    
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


class ExpPostal:


    def send_to_cluster(self,params):
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
"""

#----------------------------------------------
class ExpSetup():    
    
    
        
        
    def get_areas_section(self):
        return areaSections

    def get_areas_paper(self):
        return areasPaper

    
    def get_default_params(self):
        defaultp= {  'expname'           : None,
                'experiment_type'   : experiment_type,
                'NpHR_areas'        : areasNp,
                'ChR_areas'         : areasCh,
                'NpHR_times'        : NpHR_expression,
                'ChR_times'         : ChR_expression,
                'description'       : '',
                'soma_stim_DC'      : 1.,
                'iclamp_start'      : 100.,
                'iclamp_duration'   : 100.,
                'EPSPamp'           : 0.5
              }
        
        defaultp['savedata'] =True
        defaultp['logdata'] = True
        defaultp['tstart'] = 0 
        defaultp['tstop'] = 4000
        defaultp['num_threads'] = 4
        return defaultp
        


        
        

    
    def populate_opsin_dict(self,expdict):
        #TODO: populate
        return expdict
    
    
    def generate_params(self,expname,expdict,runon='saveonly'):            
        
        paramlist = []
        expdict['expname_family']=expname   
        
        for (i,khr) in enumerate(expdict['NpHR_areas'].keys()): # loop for NpHR expression
    
            for (j,kch) in enumerate(expdict['ChR_areas'].keys()):
    
                
                ChRlocations = expdict['ChR_areas'][kch]
                ChR_expression = expdict['ChR_times']
                ChRvalues= [[x,ChR_expression] for x in ChRlocations] 
                ChR_descript = "ChR_"+kch
    
                NpHRlocations = expdict['NpHR_areas'][khr]
                NpHR_expression = expdict['NpHR_times']
                NpHRvalues = [[x,NpHR_expression] for x in NpHRlocations]
                NpHR_descript = "NpHR_"+khr
    
                opdict = {  'ChR': ChRvalues , 
                            'NpHR': NpHRvalues} 
                
                #expdict['expname']=expname+expdict['description']+"_%g_"%i+NpHR_descript+"_%g_"%j+ChR_descript
                expdict['expname']=expname+expdict['description']+'_%s_%s'%(NpHR_descript,ChR_descript)
                expdict['opdict']=opdict
    
                if expdict['logdata']:
                    fio.log_message('About to run an experiment name:%s\nParams:'%expname)
                    for p in sorted(expdict.keys(),key=str.lower):
                        fio.log_message("Param: %s = \t%s"%(p,expdict[p]),timenow=False)
                        
                print expdict['expname']
                paramlist.append(copy.deepcopy(expdict))
                
        return paramlist

    
    def run_single_experiment(self,expbase,runon='local',params={}):
        """
        
        """
        self.main(expbase,params,runon=runon)                    
    
    def run_frequency_range(self,expbase,freqs=[0.5,1,2.5,5,10],pulsewidth=10.,transient=200.,n_pulses=10,runon='local',params={}):
        """
        
        Params:
            freqs         list of frequencies (Hz)
            pulsewidth    either list (same length as freqs) or scalar, to indicate pulse width (ms)
            transient     transient before pulses start (ms)
        """
        conversion = 1000.
        for f in freqs:
            interstim_interval = 1./f*conversion #convert freq from /s to /ms
            offset = interstim_interval/2
            t_off = interstim_interval - pulsewidth
            ChR_expression= [OPSIN_EXP,OPSIN_IRR,transient,pulsewidth,t_off,OPSIN_EXPDIST,n_pulses]
            NpHR_expression = [OPSIN_EXP,OPSIN_IRR,transient+offset,pulsewidth,t_off,OPSIN_EXPDIST,n_pulses]
            params.update({'NpHR_times':NpHR_expression,'ChR_times':ChR_expression,'description':'_freq%g_pw%g'%(f,pulsewidth)})
            self.main(expbase,params,runon=runon)

    
    def run_irraidiance_range(self,expbase,irr_range=[1,2,5,10,20,50,100],pulsewidth=10.,transient=200.,n_pulses=10,freq=2.,runon='local',params={}):
        """
        
        Params:
            irr_range     list of irradiances (mW.mm-2)
            pulsewidth    either list (same length as freqs) or scalar, to indicate pulse width (ms)
            transient     transient before pulses start (ms)
        """
        #TODO: extend to allow different irr values for different opsin
        freq_conversion = 1000.
        for irr in irr_range:
            interstim_interval = 1./freq*freq_conversion #convert freq from /s to /ms
            offset = interstim_interval/2
            t_off = interstim_interval - pulsewidth
            ChR_expression= [OPSIN_EXP,irr,transient,pulsewidth,t_off,OPSIN_EXPDIST,n_pulses]
            NpHR_expression = [OPSIN_EXP,irr,transient+offset,pulsewidth,t_off,OPSIN_EXPDIST,n_pulses]
            params.update({'NpHR_times':NpHR_expression,'ChR_times':ChR_expression,'description':'_irr%g_freq%g_pw%g'%(irr,freq,pulsewidth)})
            self.main(expbase,params,runon=runon)
    
    def run_test_exp(self,expbase,runon='local',params={}):
        """
        A simple self-contained test case that can be used without having to set up too many params.
        Should eventually be migrated to a JUnit-like testcase
        """
        transient=100
        pulsewidth = 50
        t_off = 1
        n_pulses = 1
        offset = 0
        ChR_expression= [OPSIN_EXP,OPSIN_IRR,transient,pulsewidth,t_off,OPSIN_EXPDIST,n_pulses]
        NpHR_expression = [OPSIN_EXP,OPSIN_IRR,transient+offset,pulsewidth,t_off,OPSIN_EXPDIST,n_pulses]
        
        params.update({'NpHR_times':NpHR_expression,
                       'NpHR_areas':  dict( areasAll.items()),# + areasNone.items()),
                       'ChR_times':ChR_expression,
                       'ChR_areas': dict(areasSoma.items() + areasAll.items() + areasNone.items()), #
                       'description':'_testexperiment'})
        
        self.main(expbase,params,runon=runon)
        

    
    def submit_to_cluster(self,expname):
        jdffile = fio.write_jdf_stimulation(expname)
        call(['qsub',jdffile])
        
        
    def main_run_local(self,expname,newparams):
        """
        
        """
        self.expengine = ExpEngine()
        print "Running for experiment", expname, 'locally'
        params = self.get_default_params()
        params.update(newparams)
        fio.setup_experiment(expname)
        exps = self.generate_params(expname,params)
        #TODO work out a better way of handling this
        if len(exps)>params['num_threads']:
            print 'WARNING - you have more jobs than threads. Errors will probably occur'
            print '-------------> we will set the number of workers == number of threads. Eek'
            self.expengine.start_queue(len(exps))
        else:
            self.expengine.start_queue(params['num_threads'])
        
        for exp in exps:
            self.expengine.add_experiment(exp)
        self.expengine.end_engine()
        fio.finish_experiment(expname)

    def main_run_cluster(self,expname,newparams,runjob=True):
        """
        
        """
        print "Running for experiment", expname, 'on cluster'
        params = self.get_default_params()
        params.update(newparams)
        fio.setup_experiment(expname)
        exps = self.generate_params(expname,params)
        for exp in exps:
            fio.savejob(exp,exp['expname'])
            if runjob:
                self.submit_to_cluster(exp['expname'])
       

    def main(self,expname,newparams,runon='local'):
        """
        
        """
        if runon=='local':
            self.main_run_local(expname,newparams)
        elif runon=='cluster' or runon=='saveonly':
            self.main_run_cluster(expname,newparams,runjob=(runon=='cluster'))
        else:
            raise Exception('Unknown option for runon: %s'%runon)



if __name__ == '__main__': 
    try:
        expname = sys.argv[1]
    except:
        expname = raw_input("No experiment name specified. Please enter: ")
    ExpSetup.main(expname)


