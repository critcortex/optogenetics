import pickle
import csv

import os, sys
import logging
import time

OWNER = 'sjarvis1'
EXECBASE = '~/exp/optogenetics/'
EXPBASE  = '~/exp/optogenetics/experiments'

DEFAULT_WALLTIME = 3

basedir = 'experiments/%s'
subfolders = ['img','dat','pkl','out','gdf']
filemaps = {    '*png' : 'img',
                '*dat' : 'dat',
                '*pkl' : 'pkl',
                '*jdf' : 'out',
                '*gdf' : 'gdf',
                '*.err': 'out',
                '*.out'  : 'out'}

unpackmap = { '*pkl' : 'pkl',
           '*jdf' : 'out'}

def setup_experiment(expname):
    """
    Makes relevant directories under experiments folder 
    """
    #global logger
    
    for subfolder in subfolders:
        directory = basedir%(expname+'/'+subfolder)
        if not os.path.exists(directory):
            os.makedirs(directory)
            print "Created directory:", basedir%(expname+'/'+subfolder)
        
    # creates log file
    print "Creating expDescript file for experiment %s"%(expname)
    logname = get_logfile(expname)
    logging.basicConfig(filename=logname,level=logging.INFO)


def getexps():
    return os.walk(os.path.expanduser(EXPBASE)).next()[1]

    

def finish_experiment(expname):
    """
    Move files to relevant folders, etc.
    """
    #TODO: log file
    for (k,v) in filemaps.iteritems():
        print '%s%s'%(expname,k)
        print basedir%(expname+'/'+v)
        #shutil.move('./%s%s'%(expname,k),  basedir%(expname+'/'+v))
        os.system('mv %s%s %s'%(expname,k, basedir%(expname+'/'+v)))
        print "Moved %s%s to %s"%(expname,k,basedir%(expname+'/'+v))
        
def unpack(expname):
    for (k,v) in unpackmap.iteritems():
        os.system('mv %s/%s .'%(basedir%(expname+'/'+v),k))
        print 'Unpacked %s/%s to .'%(basedir%(expname+'/'+v),k)

def get_exp_dat_location(expname):
    return basedir%expname+'/'+filemaps['*dat']


def get_logfile(expname):
    return '%s'%(expname)+'_expDescript.txt'
        

def log_message(message,timenow=True,level='info'):
    #TODO make level dynamic
    timestr = ''
    if timenow:
        timestr = '[%s]'%time.strftime('%d/%m/%y %HH:%MM')
    logging.info("%s %s"%(timestr,message))

#--------------------------------------------------------------------------------------
def loadjob(filename):
    try:
        f = open('%s.pkl'%filename, 'r') 
        expdict = pickle.load(f)
        print 'Job loaded: %s.pkl'%filename
        return expdict
    except:
        print 'Could not load job:',filename
        raise Exception('Unloaded job')
    
    
def savejob(expdict,filename):
    f = open('%s.pkl'%filename, 'w')
    pickle.dump(expdict, f) 
    print 'Job saved: %s.pkl'%filename
    

def loadresults(filename):
    try:
        f = open('%s_output.pkl'%filename, 'r') 
        expdict = pickle.load(f)
        print 'Job loaded: %s_output.pkl'%filename
        return expdict
    except:
        print 'Could not load output for job:',filename
        raise Exception('Unloaded job')
    
    
def saveresults(expdict,filename):
    f = open('%s_output.pkl'%filename, 'w')
    pickle.dump(expdict, f) 
    print 'Job saved: %s_output.pkl'%filename
    
    

    
#--------------------------------------------------------------------------------------
   
def get_pbs_header():
    txt = ("#PBS -N %s  \n"   #jobname
          "#PBS -e %s.err \n"  
          "#PBS -o %s.out \n" 
          "#PBS -l ncpus=1 -l nodes=1 \n" 
          "#PBS -l walltime=%g:00:00 \n")

    return txt

def get_pbs_simcode():
    txt = ("cd %s \n"
           "python run_stimulation.py run_from_file %s")
    return txt
   
   
def write_jdf_stimulation(expname):
    f = open('%s_sim.jdf'%expname, 'w')
    f.write(get_pbs_header()%(expname+'_sim',expname,expname,DEFAULT_WALLTIME))
    f.write(get_pbs_simcode()%(EXECBASE,expname))
    f.close()
    return expname+'_sim.jdf'

if __name__ == '__main__':
    """
    Params
    ------
    python file_io.py clean 
        Cleans all experiments 
    python file_io.py clean <expname>
        Cleans files for expname only
    """
    
    if len(sys.argv) > 1:
        if sys.argv[1] == 'clean':
            if len(sys.argv) == 2:
                for exp in getexps():
                    finish_experiment(exp)
            else:
                finish_experiment(sys.argv[2])
        elif sys.argv[1] == 'unpack':
            unpack(sys.argv[2])
        
        
        
        
        