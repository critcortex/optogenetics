import pickle
import csv
import numpy as np
import os, sys
import glob
import logging
import time

OWNER = 'sjarvis1'
EXECBASE = '~/git/optogenetics/'
EXPBASE  = '~/git/optogenetics/experiments'
EXPBASE_ONLY = 'experiments' # equivalent to EXPBASE - EXECBASE

GROUPSERV_LOC = 'smb://bg-thefarm-2012/schultz_group_data/Jarvis_comptogen_simdata/neuron'
EXT_HDD = '/media/Seagate Expansion Drive/ic_desktop/git/optogenetics/experiments/%s'

DEFAULT_WALLTIME = 3

basedir = 'experiments/%s'
subfolders = ['img','dat','pkl','out','gdf']
filemaps = {    '*png' : 'img',
                '*dat' : 'dat',
                '*pkl' : 'pkl',
                '*jdf' : 'out',
                '*gdf' : 'gdf',
                '*.err': 'out',
                '*.out': 'out'}


backupmaps = {  'img' : '*png',
                'pkl' : '*pkl',
                'analysis': '*png',
                'gdf' : '*gdf'}

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

def __remove_dat_files(expname,mv_server=True,del_file=False):
    #TODO implement
    """
    Removes the dat files only when they've been scanned for the spikes 
    """
    gdfloc = basedir%(expname+'/gdf/*gdf')
    gdffiles = glob.glob(gdfloc)
    print 'run for loop ------------'
    for file in gdffiles:
        datfile = file.split('/')
        datfile = '/'.join(datfile[:2]+['dat']+['.'.join(datfile[3].split('.')[:-1])]) + '*.dat'
        print datfile
        gl_dat = glob.glob(datfile)
        if len(gl_dat)==1: # if there is more than one file, we should be careful ... 
            print 'file DOES exist'
            print('gvfs-move %s %s/%s/%s'%(gl_dat[0],GROUPSERV_LOC,expname,'dat/'))
            os.system('gvfs-move %s %s/%s/%s'%(gl_dat[0],GROUPSERV_LOC,expname,'dat/'))
        elif len(gl_dat)>1:
            print 'More than one file match for dat: ',datfile
        else:
            print 'File NOT there'
            


def archive_to_groupserver(expname):
    #TODO implement
    """
    
    
    
    # check that exp dir exists on the farm, else mkdir
    fileloc = ''
    if len(glob.glob(fileloc))>0:
        os.system('gvfs-mkdir  %s/%s'%(GROUPSERV_LOC,expname))
    # check that sub dirs exists on the farm, else mkdir them
    for (k,v) in filemaps.iteritems():
        if clause:
            os.system('gvfs-mkdir  %s/%s/%s'%(GROUPSERV_LOC,expname,v))
    # go through and move all the files
    for (k,v) in filemaps.iteritems():
        for file in file k:
            os.system('gvfs-move %s*.%s %s/%s/%s'%(expname,file,GROUPSERV_LOC,expname,v))
    """
    pass
    
    
def backup(destination='HDD',exp=None):
    """
    Move files to relevant folders, etc.
    """
    if destination == 'HDD':
        destdir = EXT_HDD

    if exp is not None:
        _backup_exp(exp,destdir)
        return

    for exp in getexps():
        _backup_exp(exp,destdir)
        return
    
def _backup_exp(expname,destdir):
    
    spacefreedest = destdir.replace(' ','\ ')
    directory = destdir%(expname)
    if not os.path.exists(directory):
        os.makedirs(directory)
        print "Created directory:", directory
    
    for (k,v) in backupmaps.iteritems():
        print '%s%s'%(expname,k)
        print basedir%(expname+'/'+k)
        
        directory = destdir%(expname+'/'+k)
        if not os.path.exists(directory):
            os.makedirs(directory)
            print "Created directory:", directory
        #shutil.move('./%s%s'%(expname,k),  basedir%(expname+'/'+v))
        
        print('cp %s %s'%(basedir%(expname+'/%s/*'%k),destdir%(expname+'/'+k+'/')))
        os.system('cp -r %s %s'%(basedir%(expname+'/%s/*'%k),spacefreedest%(expname+'/'+k+'/')))
        print "Copied %s to %s"%(basedir%(expname+'/%s/*'%k),destdir%(expname+'/'+k))    
    
    
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

def check_mainexp_location(expname):
    """
    Checks whether a file exists in main location
    """
    print 'checking', expname
    return len(glob.glob(expname))>0

def check_expfolder(expname,searchstr,filetype='dat'):
    fileloc = basedir%expname+'/'+filemaps['*%s'%filetype]+'/%s*dat'%(searchstr)
    print 'checking  = ', fileloc
    # does file exist?
    return len(glob.glob(fileloc))>0

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
    
    
def loadspikes(expname,basename):
     
    #try:
    if True:
        #print basename,expname
        EXP_GDF_LOCATION = EXPBASE_ONLY+'/%s/'+filemaps['*gdf']
        #print EXP_GDF_LOCATION
        loadfile = EXP_GDF_LOCATION%basename +"/%s.gdf"%(basename+expname)
        #print loadfile
        #print glob.glob(loadfile)
        
        #print 'Attempting to load spike file:',loadfile,
        spikes = np.loadtxt(loadfile)#,ndmin=1)
        #print '... loaded successfully',
        return spikes
    
    #except TypeError,IndexError:
    #    self.spikes = np.array([])
    #except IOError:
    #    print 'File does not exist'
    #except:
    #    print "\nUnexpected error:", sys.exc_info()[0]
    #    return 

def loadspikerange(expname,basename):
     
    #try:
    if True:
        #print basename,expname
        EXP_GDF_LOCATION = EXPBASE_ONLY+'/%s/'+filemaps['*gdf']
        #print EXP_GDF_LOCATION
        loadfile = EXP_GDF_LOCATION%basename +"/%s.gdf"%(basename+expname)
        #print loadfile
        file_coll = glob.glob(loadfile)
        spikes = []
        for f in file_coll:
            spikes.append(np.loadtxt(f))
        return spikes,file_coll
    
#--------------------------------------------------------------------------------------
   
def get_pbs_header():
    txt = ("#PBS -N %s  \n"   #jobname
          "#PBS -e %s.err \n"  
          "#PBS -o %s.out \n" 
          "#PBS -k oe \n"
          "#PBS -l ncpus=1 -l nodes=1 \n" 
          "#PBS -l walltime=%g:00:00 \n")

    return txt

def get_pbs_simcode():
    txt = ("cd %s \n"
           "python src/run_stimulation.py run_from_file %s \n\n")
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
        
        elif sys.argv[1] == 'backup':
            if len(sys.argv) == 3:
                for exp in getexps():
                    backup(sys.argv[2],exp)
            else:
                backup(sys.argv[2],sys.argv[3])
        
        
        
        
#__remove_dat_files('140305_compare_equal_branches')        
