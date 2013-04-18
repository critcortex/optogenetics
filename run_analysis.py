

import numpy as np
import pylab 
import glob
import time
import pickle as pkl
#import NeuroTools.signals  

DEFAULT_FREQ = [0.5,1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,20,30,35,40,45,50]
DEFAULT_IRRAD = [1,2,5,10,20,50,100]
DEFAULT_PARAMS = {'FREQ':DEFAULT_FREQ,'IRRAD':DEFAULT_IRRAD}
DEFAULT_ANALYSIS_SETTINGS = {'v_th': -20., 'tstart':200,'tstop':2000}

FUNCTION_LOOKUP = {'FI'     :'calc_firingrate',
                   'iPhoto' : 'calc_peak_iPhoto'}

EXP_DAT_LOCATION = 'experiments/%s/dat/'
EXP_IMG_LOCATION = 'experiments/%s/img/'
EXP_PKL_LOCATION = 'experiments/%s/pkl/'

OPSINS_IPHOTO = {'ChR': False, 
             'NpHR': True,
             'ArchT': True}

RC_SETTINGS = {'paper': {'lw': 3, 'colors':['r','b','k','g','y','r','b','k','g','y']},
               'poster': {'lw': 4, 'colors':['r','b','k','g']}}

class AnalyticFrame:
    
    def __init__(self):
        self.all_experiments = {} #TODO: include all experiments
        self.experimentset = []
        self.expplotter = ExperimentPlotter()
        self.analysis_params = DEFAULT_ANALYSIS_SETTINGS
        
    #def add_experimentset(self,expset):
    #    self.experimentset.append(expset)

    def populate_expset(self,basenames,expsubselects,explabels,exp_mainparams):
        if type(basenames) == 'list':
            if len(basenames)==len(expsubselects):
                expbases = basenames
            else:
                print 'Error!!! basenames should be the same length as expsubselects'
        else:
            expbases = [basenames for i in expsubselects]
        
        for (i,expss) in enumerate(expsubselects):
            ee = ExperimentSet(expbases[i],expss,explabels[i])
            ee.populate_experiments(exp_mainparams,self.analysis_params)
            print 'Created experiment set for : ',expbases[i], ' ',expsubselects[i]
            self.experimentset.append(ee)
            print 'Current expbases: ',[str(exp) for exp in self.experimentset]
        
    def submenu_populate_expset(self):
        try:
            expbase = ''
            while True:
                if expbase =='':
                    expbase = raw_input('Enter basename (press return when finished): ')
                else:
                    tmp = raw_input('Current basename = %s. Press y to continue, else enter basename (press return when finished): '%expbase)
                    if not tmp.upper() == 'y':
                        expbase = tmp
                if expbase == '':
                    break
                
                expsubselect = raw_input('Enter subselect for %s: '%expbase)
                print 'Enter params for subselect:'
                # for the moment, we keep these params as fixed 
                # TODO: extend to allow user-defined
                for (i,kv) in enumerate(DEFAULT_PARAMS.iteritems()):
                    (k,v) = kv
                    print '%g - %s'%(i,k)
                exp_mainparams =  v[int(raw_input('Choose type: '))]
                explabel = raw_input('Enter label: ')
#                ee = ExperimentSet(expbase,expsubselect,explabel)
#                ee.populate_experiments(exp_mainparams,self.analysis_params)
#                print 'Created experiment set for : ',expbase, ' ',expsubselect
#                self.experimentset.append(ee)
#                print 'Current expbases: ',[str(exp) for exp in self.experimentset]
                self.populate_expset(expbase,[expsubselect],[explabel],exp_mainparams)
        except:
            print 'Error: could not '

    def run_analysis_menu(self):
        txt = "Enter type of analysis to run: \n\
                0 - enter experiment bases to compare \n\
                1 - f-I curve comparison\n\
                2 - Calc peak photocurrent \n\
                A - Analyse experiment sets \n\
                [D - run default analysis with default parameters] \n\
                G - plot trends (will be prompted for type of plot) \n\
                P - print all values for all loaded experiments \n\
                Q - Quit \n > "
        while True:
            options = raw_input(txt)
            for opt in options.split():
                # ------------------------------------------------------
                if opt == '0':
                    self.submenu_populate_expset()
                # ------------------------------------------------------
                if opt == '1':
                    print 'Default frequencies are: ', DEFAULT_FREQ
                    #TODO: user input for frequencies
                    for exp in self.experimentset:
                        exp.calculate_responses('FI')
                # ------------------------------------------------------        
                if opt == '2':
                    #TODO: user input for frequencies
                    for exp in self.experimentset:
                        exp.calculate_responses('iPhoto')
                # ------------------------------------------------------
                #if opt == 'L':
                #    print 'Load experiments'
                #    for exp in self.experimentset:
                #        exp.set_experiments()
                # ------------------------------------------------------
                if opt == 'A':
                    print 'Analyse experiment sets'
                    analysis_fns = FUNCTION_LOOKUP.keys() # TODO: allow it to be user driven
                    print 'Types of analysis to be performed:',analysis_fns 
                    for exps in self.experimentset:
                        exps.run_analysis(analysis_fns)
                # ------------------------------------------------------        
                if opt.upper() == 'S':
                    print 'Save values'
                    for exp in self.experimentset:
                        exp.save_experiments()
                # ------------------------------------------------------        
                if opt.upper() == 'G':
                    print "Current plots available:"
                    kk = ['FI'] #TODO: get superset of keys for all exp.results
                    for (i,k) in enumerate(kk):
                        print '%g - %s'%(i,k)
                    #try:
                    if True:
                        exptype = int(raw_input('Select type : '))
                        expplot = self.expplotter
                        figname = raw_input('Enter default fig name: ')
                        if figname=='': figname = '%s'%time.strftime('%H%M%S')
                        getattr(expplot,'plot_%s'%kk[exptype])(self.experimentset,savefig=figname)
                        
                    #except:
                        print 'Errored'
                # ------------------------------------------------------
                if opt.upper() == 'P':
                    print 'Current expbases: '
                    for exp in self.experimentset:
                        exp.print_data()
                # ------------------------------------------------------
                if opt.upper() == 'Q':
                    return
    
        """
        try:
            expname = sys.argv[1]
        except:
            expname = raw_input("No experiment name specified. Please enter: ")
            
        print "Running analysis for experiment", expname
        """
class ExperimentPlotter:
    
    def plot_FI(self,experimentsets,savefig=None,settings='paper'):
        rcsettings = RC_SETTINGS[settings]
        #TODO: set automatic color scale, based on number of experimens in expset
        pylab.figure() 
        for (i,expset) in enumerate(experimentsets):
            pylab.plot(expset.expvariables,expset.results['FI'],lw=rcsettings['lw'],c=rcsettings['colors'][i],label=expset.get_label())
        
        # increase scale so that we have some space to place the legend
        pylab.ylim(ymax=pylab.ylim()[1]*1.2)
        # turn the right and top axes off
        #TODO: turn off right and top axes
          
        # Prettify the legend
        leg = pylab.legend(loc=2,fancybox=True)
        leg.get_frame().set_alpha(0.5) 
        
        if savefig is not None:
            figname = '%sFI_%s.png'%(time.strftime('%y%m%d'),savefig)
            pylab.savefig(figname)
            print 'Saved figure as %s'%figname


class ExperimentSet:

    def __init__(self,basename,subselect,label='no_label'):
        """
        Params
            basename        basename for experiment i.e. '130223_noChR_freq%g_pw10'
            expvariables    dictionary of variables to be compared for range. Values can be list, string or None
                            Valid keys are 'FI', etc.
        
        Also creates
            experiments    list of Experiment objects 
        """
        self.basename = basename
        self.subselect = subselect
        self.label = label
        self.expvariables = {}
        self.results = {} #TODO : get rid of this
        self.experiments = []
    
    def __str__(self):
        return self.basename

    def get_label(self):
        return self.label
    
    def populate_experiments(self,exp_params,analysis_params):
        self.expvariables = exp_params
        for v in self.expvariables:
            print v
            print self.basename , self.subselect%v
            self.experiments.append(Experiment(self.basename ,self.subselect%v,analysis_params))
    
    def save_experiments(self):
        print 'Saving ExperimentSet', self.label
        for exp in self.experiments:
            exp.save_results()
    
    #def set_experiments(self,key,analysis_params):
        #for v in self.expvariables[key]:
        #    print v
        #    print self.basename , self.subselect%v
        #    self.experiments.append(Experiment(self.basename ,self.subselect%v,analysis_params))
    #    pass
        
    def run_analysis(self,analysis_fns,recalc=False):
        print analysis_fns
        for exp in self.experiments:
            print 'About to run analysis for ',str(exp)
            exp.calculate_results(analysis_fns, recalc)
    
    def calculate_responses(self,key,recalc=False):
        print 'calculate responses', key
        # TODO: rewrite with debugger/warning info statements
        #if self.expvariables.has_key(key):
        #try:
        if True:
            if not recalc:
                if self.results.has_key(key):
                    print 'Results for %s already exists'%key 
                    return
            tmp_results = []
            for exp in self.experiments:
                tmp_results.append(getattr(exp, FUNCTION_LOOKUP[key])())
            self.results[key] = tmp_results
        #except:
            print "Could not calculate responses, unknown variable value = ", key
        #print 'Current results ',self.results
    
    def print_data(self,indent=4):
        print self.basename,self.subselect
        for (k,v) in self.results.iteritems():
            print ' '*indent, '%s = %s'%(k,v)

class Experiment:
    
    def __init__(self,basename,expname,analysis_params):
        self.basename = basename
        self.expname = expname
        self.fullname = ''
        self.results = {}
        self.analysis_params = analysis_params
        
    def __str__(self):
        return self.basename+self.expname
    
    
    def save_results(self):
        savefile = EXP_PKL_LOCATION%self.basename +"%s.pkl"%(self.basename+self.expname)
        print 'Attempting to save to file:',savefile
        output = open(savefile, 'wb')
        pkl.dump(self.results, output,-1)
        output.close()
        
        
    def load_results(self):
        loadfile = EXP_PKL_LOCATION%self.basename +"%s.pkl"%(self.basename+self.expname)
        pkl_file = open(loadfile, 'rb')
        self.results = pkl.load(pkl_file)
        pkl_file.close()
        
        
    def calculate_results(self,properties=[],recalc=False):
        """
            Assumes that self.results is already loaded
        """
        for prop in properties:
            print 'Looking at property'
            if not recalc and self.results.has_key(prop):
                continue
            elif recalc:
                print 'Calculating value for %s'%prop
            elif not self.results.has_key(prop):
                print 'Value for %s missing; running %s'%(prop, FUNCTION_LOOKUP[prop])
            #result = getattr(self, FUNCTION_LOOKUP[prop])()
            #self.results[prop] = result
            getattr(self, FUNCTION_LOOKUP[prop])()
            
                

    
    def get_dat_file(self,modifier=''):
        try:
            datfile = EXP_DAT_LOCATION%self.basename +"%s%s.dat"%(self.basename+self.expname,modifier)
        except:
            print 'Error with finding file'
            return 
        gg = glob.glob(datfile)
        if len(gg)==1:
            self.fullname = gg[0][:-4] # name of .dat file without the file extension
            return gg[0]
        elif len(gg)==0:
            #TODO: raise error
            print 'No file found: ',file
        elif len(gg)>1:
            #TODO: raise error - too many files
            print 'Too many files found: ',file

    def calc_peak_iPhoto(self):
        for (opsin,values) in OPSINS_IPHOTO.iteritems():
            try:
                self.calc_peak_iPhoto_opsin(opsin,values)
            except:
                pass

    def calc_peak_iPhoto_opsin(self,opsintype,max_iphoto):
        """
            opsin_type
            max_iphoto
        
        """
        data = np.loadtxt(self.get_dat_file('*_i%s'%opsintype)) 
        i_photo = data[0,:]
        if max_iphoto:
            peak_iphoto =  i_photo.max()
        else:   #no prizes for guessing correctly ... we calculate the min
            peak_iphoto = i_photo.min()
        self.results['iPhoto_%s'%opsintype] = peak_iphoto
        return peak_iphoto

    def calc_firingrate(self,**params):
        """
            expname     experiment name
            t_cutout    list of start time and finish times during which firing rate should be calculated
        """
        #load dat curve, get v_soma
        
        data = np.loadtxt(self.get_dat_file()) 
        t,v_soma = data[0,:],data[1,:]
        spikes = self.__extract_spiketimes(t,v_soma)
        if len(spikes)==0:
            print 'No spikes observed during entire duration'
            return 0
        spikes = np.where((spikes >= self.analysis_params['tstart']) & (spikes <= self.analysis_params['tstop']))[0]
        print 'calc_firing rate'
        spikerate = 1000.*len(spikes)/(self.analysis_params['tstop']-self.analysis_params['tstart'])
        self.results['FI'] = spikerate
        return spikerate
    
    def __extract_spiketimes(self,times,v_soma):
        """
        This section is a rewritten version of NeuroTools
        
        """
        above = np.where(v_soma > self.analysis_params['v_th'])[0]
        if len(above) <= 0:
            return []
        else:
            take = np.where(np.diff(above)>1)[0] + 1
            take = np.append(0,take)
            
            return times[above][take]
        

