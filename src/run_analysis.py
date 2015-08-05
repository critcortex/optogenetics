import sys

from matplotlib import cm, colors
from scipy.optimize import curve_fit
import numpy as np
import pylab 
import glob
import time
import pickle as pkl
import itertools
import NeuroTools.signals  
import math
from xdg.Menu import tmp

DEFAULT_FREQ = [0.5,1,1.5,2,2.5,4,5,6,7,8,9,10,15,20,30,35,40,45,50]
DEFAULT_IRRAD = [1,2,5,10,20,50,100]
DEFAULT_PARAMS = {'FREQ':DEFAULT_FREQ,'IRRAD':DEFAULT_IRRAD}
DEFAULT_ANALYSIS_SETTINGS = {'v_th': -20., 'tstart':200,'tstop':2000,'prebuffer':300,'postbuffer':150,'tstart_bg':0,'tstop_bg':0,'tstart_post':0,'tstop_post':0}

FUNCTION_LOOKUP = {'isi':'calc_isi',
                   'isi_bg':'calc_isi_bg',
                   'isi_post':'calc_isi_post',
                   'cv_isi':'calc_feature',
                   'cv_kl':'calc_feature',
                   'fano_factor_isi':'calc_feature',
                   'mean_rate':'calc_feature',
                   'spiketimes':'process_spiketimes',
                   'FI'     :'calc_firingrate',
                   'FI_bg'  : 'calc_firingrate_bg',
                   'FI_post'  : 'calc_firingrate_post',
                   'iPhoto' : 'calc_peak_iPhoto'}

EXP_DAT_LOCATION = 'experiments/%s/dat/'
EXP_IMG_LOCATION = 'experiments/%s/img/'
EXP_PKL_LOCATION = 'experiments/%s/pkl/'
EXP_GDF_LOCATION = 'experiments/%s/gdf/'

CMAPS = ['jet','YlOrRd','YlGnBu_r','spectral','YlGnBu_r','YlOrRd','GnBu','jet','CMAP_ORANGE_BLUE','CMAP_BLUE_ORANGE','CMAP_FULL_BLUE_ORANGE','r_pink','hot_r']


OPSINS_IPHOTO = {'ChR': False, 
             'NpHR': True,
             'ArchT': True}

RC_SETTINGS = {'paper': {'lw': 3, 
                         'lw_fine':1,
                         'alpha':0.5,
                         'dpi':100,
                         'colors':['r','b','k','g','y','r','b','k','g','y'],
                         'cmaps':CMAPS},
               'poster': {'lw': 4, 
                          'lines.linewidth':4,
                          'axes.linewidth':2,
                          'dpi':300,
                          'figure.dpi':300,
                          'savefig.dpi':300,
                          'lw_fine':2,
                          'alpha':0.5,
                          'font.size':16,
                          'colors':['r','b','k','g'],
                          'cmaps':CMAPS}}

class AnalyticFrame:
    
    def __init__(self,cmap_index=0):
        self.all_experiments = {} #TODO: include all experiments
        self.experimentset = []
        self.expplotter = ExperimentPlotter(cmap_index)
        self.analysis_params = DEFAULT_ANALYSIS_SETTINGS
        
    def update_cmap(self,new_cmap_index):
        self.expplotter.update_cmap_index(new_cmap_index)
    
    def update_params(self,newparams):
        self.analysis_params.update(newparams)
        print 'Updated parameters. New parameters are:',self.analysis_params

    """
    def set_single_exp(self,basename,expfilename):
        ee = ExperimentSet(basename,expfilename)
        ee.experiments.append(Experiment(basename ,expfilename,self.analysis_params))
    """ 
    def populate_trialset(self,basenames,trialInstances,variables,trials,trialLabels,label='no_label',var_format="%g"):
 
        for (i,instance) in enumerate(trialInstances):
            if type(basenames) is list:
                basename = basenames[i]
            else:
                basename = basenames
            if type(trialLabels) is list:
                    
                testlabel = trialLabels[i]
            else:
                testlabel = "%s_%g"%(trialLabels,i)   
            
            tes = TrialExperimentSet(basename,[instance],variables,trials,testlabel,label,var_format,self.analysis_params)
            self.experimentset.append(tes)

    def populate_expset(self,basenames,expsubselects,explabels,exp_mainparams):
        
        if type(basenames) is list:
            if len(basenames)==len(expsubselects):
                expbases = basenames
            else:
                print 'Error!!! basenames should be the same length as expsubselects'
        else:
            expbases = [basenames for i in expsubselects]
            
        print len(expbases),expbases
        print len(expsubselects), expsubselects
        print len(explabels), explabels
        
        for (i,expss) in enumerate(expsubselects):

            ee = ExperimentSet(expbases[i],expss,explabels[i])
            ee.populate_experiments(exp_mainparams,self.analysis_params)
            print 'Created experiment set for : ',expbases[i], ' ',expsubselects[i]
            
            self.experimentset.append(ee)
            # TODO: correct following line so that it works with list of basenames
            #print 'Current expbases: ',[str(exp) for exp in self.experimentset]
        
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
            
            
            
    def submenu_analysis_params(self):
        """
        
        """
        try:
            update = raw_input('Would you like to update params? y/(n):')
            if update.upper() != 'Y':
                return
            else:
                print 'Updating analysis parameters. Press enter to keep current value'
                for (k,v) in self.analysis_params.iteritems():
                    new_value = raw_input('Param[%s]. Current value = %s \t\t. Enter new value? '%(k,v))
                    if new_value != '':
                        self.analysis_params[k] = new_value
        except:
            print 'Error with updating analysis params.'
            return 

    def submenu_extractSpikes(self):
        print 'Default frequencies are: ', DEFAULT_FREQ
        #TODO: user input for frequencies
        for exp in self.experimentset:
            exp.calculate_responses('spiketimes')    
    
    def submenu_runFI(self):
        print 'Default frequencies are: ', DEFAULT_FREQ
        #TODO: user input for frequencies
        for exp in self.experimentset:
            exp.calculate_responses('FI')    
            
    
    def submenu_Iphoto(self):
        #TODO: user input for frequencies
        for exp in self.experimentset:
            exp.calculate_responses('iPhoto')
        
    def submenu_analyse(self):
        analysis_fns = FUNCTION_LOOKUP.keys() # TODO: allow it to be user driven
        print 'Types of analysis to be performed:',analysis_fns 
        for exps in self.experimentset:
            exps.run_analysis(analysis_fns)        
            
    def perform_analysis(self,analysis_features=[],recalc=False):
        #for feature in analysis_features:
        #    print feature,'='*60
        for exps in self.experimentset:
            exps.run_analysis(analysis_features,recalc)
        
    def submenu_save(self):
        print 'Save values'
        for exp in self.experimentset:
            exp.save_experiments()        

    def submenu_load(self):
        print 'Load values'
        for exp in self.experimentset:
            exp.load_experiments()       

    def submenu_plot(self,exptype=None,supplied_figname=None,**kwargs):
        print "Current plots available: ----------------------------------"
    
        kk = ['FI','IClamp','IClamp_iPhoto','IClamp_iPhoto_peak','FI_compare','fit_FI','PSTH_population','population_raster','population_voltage','compare_features','fit_FI_bg','compare_2D_FI'] 
        #TODO: get superset of keys for all exp.results
        for (i,k) in enumerate(kk):
            print '%g - %s'%(i,k)
        
        if supplied_figname is None and exptype is None:
            exptype = int(raw_input('Select type : '))
            
            figname = raw_input('Enter default fig name: ')
            if figname=='': figname = '%s'%time.strftime('%H%M%S')
        else:
            figname = supplied_figname

        #print 'params = ',params
        print 'kwargs = ',kwargs
        print figname, exptype, kk[exptype]
        getattr(self.expplotter,'plot_%s'%kk[exptype])(self.experimentset,savefig=figname,**kwargs)
        
    
        
    def submenu_print(self):
        
        for exp in self.experimentset:
            exp.print_exps()
            exp.print_data()
                        
    def run_analysis_menu(self):
        txt = "Enter type of analysis to run: \n\
                0 - enter experiment bases to compare \n\
                1 - f-I curve comparison\n\
                2 - Calc peak photocurrent \n\
                A - Analyse experiment sets \n\
                [D - run default analysis with default parameters] \n\
                S - save/update experimental analysis \n\
                L - load experimental analysis \n\
                M - print and set analysis parameters \n\
                G - plot trends (will be prompted for type of plot) \n\
                P - print all values for all loaded experiments \n\
                Q - Quit \n > "
        while True:
            options = raw_input(txt)
            for opt in options.split():
                #------------------------------------------------------
                if opt == '0':
                    self.submenu_populate_expset()
                # ------------------------------------------------------
                if opt == '1':
                    self.submenu_runFI()
                # ------------------------------------------------------        
                if opt == '2':
                    self.submenu_Iphoto()
                # ------------------------------------------------------        
                if opt.upper() == 'A':
                    self.submenu_analyse()
                # ------------------------------------------------------        
                if opt.upper() == 'S':
                    self.submenu_save()
                # ------------------------------------------------------        
                if opt.upper() == 'L':
                    self.submenu_load()
                # ------------------------------------------------------  
                if opt.upper() == 'G':
                    self.submenu_plot()
                # ------------------------------------------------------
                if opt.upper() == 'P':
                    print 'Printing current expbases: '
                    self.submenu_print()
                # ------------------------------------------------------
                if opt.upper() == 'M':
                    print 'Current analysis settings are:',self.analysis_params
                    self.submenu_analysis_params()
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
        
        
    def _calculateModulationIndex(self,mitype='MI'):
        """
        
        @param 
            mitype : MI or MI_bg
        """
        
        m_ChR2 = self.experimentset[0].results[mitype][1] #ChR2 max
        m_NpHR = self.experimentset[-1].results[mitype][0] # NpHR max
        return m_ChR2/m_NpHR
        
    def get_analysis_values(self,key):
        """
        
        """
        results = []
        for exp in self.experimentset:
            tmp = exp.get_property(key)
            print tmp
            results.append(tmp)
        return results
        
class ExperimentPlotter:
    
    def __init__(self,cmap_index=0):
        self.cmap_index = cmap_index
        self.__update_cmap()
        
    def update_cmap_index(self,new_cmap_index):
        #try:
        if True:
            self.cmap_index = new_cmap_index
            self.__update_cmap()
        #except:
        #    print("Couldn't change cmap -> leaving it unchanged")
        
    def __update_cmap(self,settings='paper'):
        rcsettings = RC_SETTINGS[settings]
        cmap = rcsettings['cmaps'][self.cmap_index]
        print cmap, type(cmap)
        if not cmap.startswith('CMAP'):
            #cNorm = colors.Normalize(vmin=0,vmax=1) 
            #self.cmap = colors.Colormap(rcsettings['cmaps'][self.cmap_index])
            self.cmap = rcsettings['cmaps'][self.cmap_index]
            print "---> ", self.cmap, type(self.cmap)
        else:
            self.cmap = self.load_cmap(rcsettings['cmaps'][self.cmap_index])
    
    def print_header(self,plottype,expset):
        print '='*40
        print 'Creating %s plot for expset = %s'%(plottype,str(expset))
        print '='*40
        
    def load_cmap(self,filename):
        """
        Loads array from filename.txt -->
        
        Expects that each line is an R,G,B, value i.e.
        
        0,61,196,
        6,73,198,
        12,84,199,
        18,95,200,
        
        
        NOTE: If using a site to generate the colormap i.e. colormap.org , then use 
            sed -i 's/;/,/g' FILENAME.txt
            sed -i 's/]/,/g' FILENAME.txt
            sed -i 's/[//g' FILENAME.txt
              
        to convert from 
            [0,61,196;
            6,73,198;
            0,61,196;
            6,73,198] into the desired format
        
        """
        a = np.genfromtxt(filename+'.txt',delimiter=',')
        a = a[:,:3]
        cm = colors.ListedColormap(a/255)
        return cm
        
    
    
    def __get_color(self,cmap,value=0.8):
        cNorm = colors.Normalize(vmin=0,vmax=1) 
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=self.cmap)
        return scalarMap.to_rgba(value)
    
    def __turn_off_border(self,ax,turnoffs=['right','top']):
        for loc, spine in ax.spines.iteritems():
            if loc not in turnoffs:
                spine.set_position(('outward',5))
                ax.tick_params(direction='out')
            elif loc in turnoffs:
                #spine.set_color('none') # don't draw spine
                spine.set_visible(False)
                #ax.tick_params([])
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
    def plot_FI(self,experimentsets,savefig=None,settings='paper',removeNulls=True):
        pylab.close('all')
        rcsettings = RC_SETTINGS[settings]
        print "-="*20,'cmapindex = ',self.cmap_index
        pylab.figure() 
        intensities = np.linspace(0.1, 0.9, len(experimentsets))

        for (i,expset) in enumerate(experimentsets):
            try:
                print 'i=%g, intensities[i]='%i,intensities[i]
                print 'exp variables len', len(expset.expvariables[0])
                print expset.expvariables[0]
                print 'length FI', len(expset.results['FI'])
                if removeNulls:
                    # find all instances where FI == 0 and remove corresponding expvariables and FI entries
                    ys = expset.results['FI']
                    xs = expset.expvariables[0]
                    obs = np.nonzero(ys)[0][0]
                    print type(ys), type(xs)
                    print "---------------------------", obs, ys
                    print "---------------------------", obs, xs
                    print "---------------------------", obs
                    print "---------------------------", obs
                    print "---------------------------", obs
                    print "---------------------------", obs
                    print "---------------------------", obs
                    pylab.plot(xs[obs:],ys[obs:],lw=rcsettings['lw'],c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),label=expset.get_label())
                else:
                    pylab.plot(expset.expvariables[0],expset.results['FI'],lw=rcsettings['lw'],c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),label=expset.get_label())
            except:
                print 'Error with plotting for expset', expset
                continue
        # turn the right and top axes off
        self.__turn_off_border(pylab.gca())
        
        pylab.ylabel('Output frequency (Hz)')
        pylab.xlabel('Input')
        pylab.ylim(ymax=pylab.ylim()[1]*1.2,ymin=-5)
        
        if savefig is not None:
            figname = '%sFI_%s.png'%(savefig,time.strftime('%y%m%d'))
            pylab.savefig(figname)
            print 'Saved figure as %s'%figname        

        
        
        
    def _sort_xsys(self,xs,ys):
        data = zip(xs,ys)
        data.sort(key=lambda tup: (tup[0],tup[1]))
        xs,ys = zip(*data)
        return xs,ys
            
    def _linear(self,x,a,b,c,d):
        return a*x + b
    
    def _d_linear(self,x,a,b,c,d):
        return a
    
    def _log(self,x,a,b,c,d):
        return a*np.log(x-b)+c
    
    def _d_log(self,x,a,b,c,d):
        return a*1./(x-b)
            
    def _poly(self,x,a,b,c,d):
        return a*x**3 + b*x**2 + c*x + d
    
    def _d_poly(self,x,a,b,c,d):
        return 3*a*x**2 + 2*b*x + c

    def _sigmoid(self,x,x0,y0,c,k):
        return c / (1 + np.exp(-k*(x-x0))) + y0
    
    def _d_sigmoid(self,x,x0,y0,c,k):
        return (c*k*np.exp(-k*(x-x0))) / (1 + np.exp(-k*(x-x0)))**2

    def _sigmoid_log(self,x,a,b,c,d):
        pass
      
    def _log_heaviside(self,x,a,b,c,p):
        y = np.zeros(x.shape)
        for i in range(len(y)):
            y[i]=self._log_heaviside_inner(x[i],a,b,c,p)
        return y  
        
        
    def _log_heaviside_inner(self,x,a,b,c,p):
        if x>p:
            return a*np.log(b*x)+c
        else:
            return 0 
        
    def _power_law(self,x,a,b,c,d):
        y = np.zeros(x.shape)
        for i in range(len(y)):
            #print "i ==", i, "x[i] --> ",x[i],
            y[i]=self._power_law_inner(x[i],a,b,c,d)
            if y[i]<0:
                y[i] = 0 
        return y  
        
        
    def _power_law_inner(self,x,a,b,c,d):
        #print "(x-a), and (x-a)**b ==> ", (x-a), (x-a)**b, 
        if (x-a)>=0:
            return c*(((x-a)**b)/(1+(x-a)**b)) ###-d
        else:
            return 0 
        

    def _get_guesstimate(self,xvals,yvals,polyfn):
        
        if polyfn == '_sigmoid':
            """
            print np.median(xvals),np.median(yvals)
            y0guess =  yvals.min()
            cguess = yvals.max() - yvals.min()
            #x0guess = should happen where y0 mean value is
            return [np.median(xvals),y0guess,cguess,1.]
            """
            return [1.,1.,1.,0.]
        elif polyfn == '_poly':
            return [1.,1.,1.,0.]
        elif polyfn == '_log':
            x0 = np.nonzero(yvals)[0][0]
            print 'First nonzero at index',x0
            print xvals[x0]
            return [np.max(xvals),xvals[x0],0.,0.]
        elif polyfn == '_log_heaviside': 
            threshold = 0.1
            xcrossing = np.where(yvals>threshold)[0][0]
            return [1.,1.,1.,xvals[xcrossing]]
        elif polyfn == '_power_law':
            threshold = 0.1
            xcrossing = np.where(yvals>threshold)[0][0]
            return [xvals[xcrossing],1.,yvals[-1],0]
        
        
        else:
            return [1.,1.,1.,0.]
        
                
    def plot_fit_FI(self,experimentsets,savefig=None,settings='poster',polyfn='_poly'):
        
        pylab.close('all')
       
        self._plot_fit_FI(experimentsets,settings=settings,savefig=savefig+'_poly',polyfn='_poly',p0=[1.,1.,1.,0.])
        pylab.close('all')
        
        self._plot_fit_FI(experimentsets,settings=settings,savefig=savefig+'_sigmoid',polyfn='_sigmoid',p0=None)
        pylab.close('all')
        
        self._plot_fit_FI(experimentsets,settings=settings,savefig=savefig+'_linear',polyfn='_linear',p0=None)
        pylab.close('all')
       
        self._plot_fit_FI(experimentsets,settings=settings,savefig=savefig+'_log',polyfn='_log',p0=None)
        pylab.close('all')
        
        self._plot_fit_FI(experimentsets,settings=settings,savefig=savefig+'_log_heaviside',polyfn='_log_heaviside',p0=None)
        pylab.close('all')
       
        self._plot_fit_FI(experimentsets,settings=settings,savefig=savefig+'_power_law',polyfn='_power_law',p0=None)
        pylab.close('all')
       
        self._plot_fit_FI(experimentsets,settings=settings,savefig=savefig+'_mixed',polyfn=['_poly','_sigmoid','_linear','_log','_log_heaviside','_power_law'],p0=None)
        pylab.close('all')
        
        self._plot_fit_FI(experimentsets,settings=settings,savefig=savefig+'_mixed2',polyfn=['_poly','_sigmoid','_linear','_log'],p0=None)
        pylab.close('all')
        
        
    def _plot_fit_FI(self,experimentsets,savefig=None,settings='poster',polyfn='_poly',p0=None):
        self._plot_fit_FI_general(experimentsets,'expvariables','FI','Input',savefig,settings,polyfn,p0)
    
    
    
    def plot_fit_FI_bg(self,experimentsets,savefig=None,settings='poster',polyfn='_poly'):
        
        pylab.close('all')
        
        self._plot_fit_FI_bg(experimentsets,settings=settings,savefig=savefig+'_poly',polyfn='_poly')
        pylab.close('all')
        
        self._plot_fit_FI_bg(experimentsets,settings=settings,savefig=savefig+'_sigmoid',polyfn='_sigmoid',p0=None)
        pylab.close('all')
        
        self._plot_fit_FI_bg(experimentsets,settings=settings,savefig=savefig+'_linear',polyfn='_linear',p0=None)
        pylab.close('all')
        
        self._plot_fit_FI_bg(experimentsets,settings=settings,savefig=savefig+'_log',polyfn='_log',p0=None)
        pylab.close('all')
        
        self._plot_fit_FI_bg(experimentsets,settings=settings,savefig=savefig+'_log_heaviside',polyfn='_log_heaviside',p0=None)
        pylab.close('all')
        
        self._plot_fit_FI_bg(experimentsets,settings=settings,savefig=savefig+'_power_law',polyfn='_power_law',p0=None)
        pylab.close('all')
        
        self._plot_fit_FI_bg(experimentsets,settings=settings,savefig=savefig+'_mixed',polyfn=['_poly','_sigmoid','_linear','_log','_log_heaviside','_power_law'],p0=None)
        pylab.close('all')
        
        self._plot_fit_FI_bg(experimentsets,settings=settings,savefig=savefig+'_mixed2',polyfn=['_poly','_sigmoid','_linear','_log'],p0=None)
        pylab.close('all')        
        
        
    def _plot_fit_FI_bg(self,experimentsets,savefig=None,settings='poster',polyfn='_poly',p0=None):
        if savefig is not None:
            savefig += '_bg'
        self._plot_fit_FI_general(experimentsets,'FI_bg','FI','Background (Hz)',savefig,settings,polyfn,p0)


    def _plot_xs_ys(self,xs,ys,fittedxs,fittedys,labels,xlabel,savefig=None,settings='poster'):
        """
        Params:
            xs        list of arrays for each exp set
            ys        list of arrays for each exp set, for observed data points. Must be same dimensions as ys
            fittedxs  list of arrays for each expset
            fittedys  list of arrays for each exp set, for fitted data. ""
            xlabel    string, describing text on x-axis
            savefig
            settings
            
        """
        rcsettings = RC_SETTINGS[settings]
        pylab.rcParams.update(rcsettings)
        pylab.figure() 
        print xs, len(xs), '<---------'
        intensities = np.linspace(0.1, 0.9, len(xs))
        m = cm.ScalarMappable(cmap=self.cmap)#cmap=rcsettings['cmaps'][self.cmap_index])
        #m = cm.ScalarMappable(cmap=rcsettings['cmaps'][self.cmap_index])
        #m = cm.ScalarMappable(cmap=rcsettings['cmaps'][1])
        m.set_array(intensities)
        
        for (i,x) in enumerate(xs): 
                pylab.scatter(x,ys[i],lw=rcsettings['lw']/2,c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),s=rcsettings['lw']*10,facecolor='none',edgecolor=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),alpha=0.6)
                pylab.plot(fittedxs[i],fittedys[i],lw=rcsettings['lw'],c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]))
            
        # turn the right and top axes off
        self.__turn_off_border(pylab.gca())
        # set labels
        pylab.ylabel('Output frequency (Hz)')
        pylab.xlabel(xlabel)
        
        
        #pylab.ylim(ymin=ys[0]-np.abs(ys[0])*0.1)
        #pylab.ylim(ymax=ys[-1]+np.abs(ys[-1])*0.1)
        pylab.xlim(xmin=xs[0][0]-np.abs(xs[0][0])*0.1)
        pylab.xlim(xmax=xs[0][-1]+np.abs(xs[0][-1])*0.1)
        
        cbarlabels = labels
        dc = intensities[1]-intensities[0]
        mapint = np.linspace(0.1-dc/2, 0.9+dc/2, len(xs)+1)
        cb = pylab.colorbar(m,ticks=intensities,boundaries=mapint)
        cb.set_label('x NpHR Factor')
        cb.set_ticklabels(cbarlabels)
        
        if savefig is not None:
            figname = '%sFIfit_%s.png'%(savefig,time.strftime('%y%m%d'))
            pylab.savefig(figname,dpi=rcsettings['dpi'])
            print 'Saved figure as %s'%figname       
            figname = '%sFIfit_%s.svg'%(savefig,time.strftime('%y%m%d'))
            pylab.savefig(figname,dpi=rcsettings['dpi'])
            print 'Saved figure as %s'%figname        
        
    def get_derivative(self,fname):
        return getattr(self,'_d'+fname)
        
        
    def _fit_curve(self,xvals,yvals,p0,polyfn):
        """
        
        Returns:
            poly    function itself
            params   best fit params
            polyname name of function
            xfitted  [xstart,xend] where start is where fitting started, and xend where fitting finished
        
        """
        max_interation = 10000
        # make a single polyfn into a list so that the subsequent code works for both single and multiple functions
        if isinstance(polyfn, str):
            polyfn = [polyfn]
        poly = [getattr(self,pfn) for pfn in polyfn]
        min_rmse = float("inf")
        best_poly_index = None
        for (p,pfn) in enumerate(poly):
            if p0 is None:
                pguess = self._get_guesstimate(xvals,yvals,polyfn[p])
            else:
                pguess = p0
                
            # if poly is log_heaviside or power_law, then we have to remove
            # datapoints (x_i, y_i) where y_i=0
            #TODO:
            if polyfn[p] == '_log_heaviside' or polyfn[p] == '_power_law':
                #xs = xvals[yvals.nonzero()]
                #ys = yvals[yvals.nonzero()]
                xs = xvals[np.where(yvals>1)]
                ys = yvals[np.where(yvals>1)]
                print "we had ------------------------------"
                print xvals
                print yvals
                
                print "---- and we now have -----------------------------------------------------"
                print xs
                print ys
            else:
                xs = xvals
                ys = yvals
                
            
            if p0 is None:
                pguess = self._get_guesstimate(xs,ys,polyfn[p])
            else:
                pguess = p0
               
            
                 
            if len(polyfn)>1:
                # we're using a mix of functions
                
                # Calculate the RMSE using info from leastsq
                try:
                    popt_params, pcov,infodict, errmsg, ier = curve_fit(pfn, xdata=xs, ydata=ys,p0=pguess,full_output=True,maxfev=max_interation)
                except:
                    continue
                chisq=(infodict['fvec']**2).sum()
                # dof is degrees of freedom
                dof=len(xs)-len(popt_params)
                rmse=np.sqrt(chisq/dof)
                #print('Poly: %s, RMSE = %g'%(pfn,rmse))
                if rmse < min_rmse:
                    best_poly_index = p
                    min_rmse = rmse
                    p_use = pguess
            else:
                # if there's only one method to try .... 
                print "so we have to fit :::::::::::::::::::::::::::::::::::::::::::"
                print xs 
                print ys
                print pguess
                cparams,cpcov = curve_fit(pfn, xdata=xs, ydata=ys,p0=pguess,maxfev=max_interation)
                print "and we found ::::::::::::::::::::::::::::::::::::::" 
                print cparams
                print "------------------------------------ Got to end of _fit_curve"
                #return pfn, cparams,polyfn[0],[xs[0],xs[-1]]
                return pfn, cparams,polyfn[0],[xvals[0],xvals[-1]]
            
        best_poly = poly[best_poly_index]
        cparams,cpcov = curve_fit(best_poly, xdata=xs, ydata=ys,p0=p_use,maxfev=max_interation)
        print "Best poly was %s with a RMSE = %g, params = "%(polyfn[best_poly_index],min_rmse), cparams
        
        return best_poly, cparams, polyfn[best_poly_index],[xvals[0],xvals[-1]]
        

    def _plot_fit_FI_general(self,experimentsets,xval_key,yval_key,xlabel,savefig=None,settings='poster',polyfn='_poly',p0=None):
        """
        Does the fitting, and hands over to _plot_xs_ys to do the plotting
        
        Params
            experimentsets
            xval_key
            yval_key
            xlabel
            savefig
            settings
            polyfn
            p0
        
        """
        
        
        obsxs = []
        obsys = []
        fittedxs = []
        fittedys = []
        labels = []
        
        for expset in experimentsets:
            
            #if True:
            try:
                ys = expset.get_property(yval_key)
                
                if xval_key == 'expvariables':
                    # need to take the first expvariables value ... sigh. Such a good reminder to plan first, code second
                    # TODO : is there a better way to grab the xvals for this scenario
                    xitem = expset.get_property(xval_key)
                    print xitem
                    if type(xitem) is int or type(xitem) is float: 
                        print "uhhuh - got an int/float"
                    elif type(xitem[0]) is int or type(xitem[0]) is float:
                        print "gotta list - good"
                        xs = xitem
                    else:
                        xs = xitem[0]
                    #xs = expset.get_property(xval_key,0)
                    print xs
                    xvals = np.linspace(xs[0], xs[-1],num=len(xs)) 
                else:
                    xs = expset.get_property(xval_key)
                    xvals = xs
                    
                xs,ys = self._sort_xsys(xs, ys)
                # OLD
                svals = np.asarray(xvals).ravel() 
                # NEW
                #TODO: change this to the correct - which wshould be xs
                #svals = np.asarray(xs).ravel() 
                fis = np.asarray(ys).ravel()
                
                
                poly,popt_curve,pname,xfitted = self._fit_curve(svals, fis,p0, polyfn)
                
                upscale_factor = 4
                #if polyfn == '_power_law':
                #    upscale_factor = 4
                print "Plotting fit from ", xfitted[0]
                fitted_xvals = np.linspace(xfitted[0], xfitted[-1],num=len(xs)*upscale_factor) # add more data points
                fitted_svals = np.asarray(fitted_xvals).ravel() 
                
                fitted = poly(fitted_svals,popt_curve[0],popt_curve[1],popt_curve[2],popt_curve[3])
                print "AAAA---------------------"
                print xs
                print ys
                print popt_curve[0],popt_curve[1],popt_curve[2],popt_curve[3]
                print "---------------------BBB"
                
                
                # OLD
                obsxs.append(xvals)
                # NEW
                #obsxs.append(svals)
                obsys.append(fis)
                fittedxs.append(fitted_svals)
                fittedys.append(fitted)
                labels.append(expset.get_label())
                
            except:
                # Usually a RuntimeError that the optimal parameters has not been found
                print 'Error with plotting for expset', expset
                print "Unexpected error:", sys.exc_info()[0]
                continue
        #print '--> ',obsxs, fitted
        # if we actually have something to plot
        print "-=-=-=-=-=-"
        print obsxs
        print obsys
        print "-=-=-=-=-=-"
        if len(labels)>0:    
            self._plot_xs_ys(obsxs, obsys, fittedxs, fittedys, labels, xlabel, savefig, settings)
            
            
            
    def plot_FI_compare(self,experimentsets,savefig=None,settings='paper'):
        rcsettings = RC_SETTINGS[settings]
        pylab.rcParams.update(rcsettings)
        pylab.close('all')
        
        for (i,expset) in enumerate(experimentsets):
            pylab.figure(figsize=(5,5))
            
            
            data = np.array(expset.results['FI'])
            print data
            print len(data)
            print len(expset.expvariables[0]), len(expset.expvariables[1])
            print len(expset.expvariables[0]) * len(expset.expvariables[1])
            data = data.reshape(len(expset.expvariables[0]),len(expset.expvariables[1]))
             
            pylab.pcolor(data,vmax=30.,cmap=cm.YlGnBu) #@UndefinedVariable
            # turn the right and top axes off
            ax = pylab.gca()
            self.__turn_off_border(ax)
            
            
            plotxrange = range(0,len(expset.expvariables[1]),5)
            plotyrange = range(0,len(expset.expvariables[0]),5)
            ax.set_xticks(plotxrange)
            ax.set_xticklabels([expset.expvariables[1][x] for x in plotxrange])
            ax.set_yticks(plotyrange)
            ax.set_yticklabels([expset.expvariables[0][x] for x in plotyrange])
            
            pylab.ylim(0,len(expset.expvariables[0]))
            pylab.xlim(0,len(expset.expvariables[1]))
            
            #pylab.ylabel('Axis 0')
            #pylab.ylim(ymax=expset.expvariables[0][-1])
            #pylab.xlabel('Axis 1')
            #pylab.xlim(xmax=expset.expvariables[1][-1])
            cb = pylab.colorbar()
                    
            if savefig is not None:
                figname = '%s_%s_FIcompare_%s.png'%(expset,savefig,time.strftime('%y%m%d'))
                pylab.savefig(figname,dpi=rcsettings['dpi'])
                print 'Saved figure as %s'%figname
            pylab.close('all')
            
            
    def plot_FIphoto(self,experimentsets,savefig=None,settings='paper'):
        #TODO: plot
        pass
        """
        rcsettings = RC_SETTINGS[settings]
        pylab.figure() 
        intensities = np.linspace(0.1, 0.9, len(experimentsets))
        for (i,expset) in enumerate(experimentsets):
            print 'i=%g, intensities[i]='%i,intensities[i]
            pylab.plot(expset.expvariables,expset.results['FI'],lw=rcsettings['lw'],c=self.get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),label=expset.get_label())
        pylab.ylabel('Output frequency (Hz)')
        pylab.xlabel('Input frequency (Hz)')
        # increase scale so that we have some space to place the legend
        pylab.ylim(ymax=pylab.ylim()[1]*1.2)
        
        # turn the right and top axes off
        ax = pylab.gca()
        for loc, spine in ax.spines.iteritems():
            if loc in ['left','bottom']:
                spine.set_position(('outward',5))
                ax.tick_params(direction='out')
            elif loc in ['right','top']:
                spine.set_color('none') # don't draw spine
        
        # Prettify the legend
        leg = pylab.legend(loc=2,fancybox=True)
        leg.get_frame().set_alpha(0.5) 
        
        if savefig is not None:
            figname = '%sFI_%s.png'%(time.strftime('%y%m%d'),savefig)
            pylab.savefig(figname)
            print 'Saved figure as %s'%figname
        """

    def plot_IClamp(self,experimentsets,savefig=None,settings='paper'):
        """
        Plots soma voltages of multiple experiments onto same figure
        """
        rcsettings = RC_SETTINGS[settings]
        for (i,expset) in enumerate(experimentsets):
            pylab.close('all')
            pylab.figure() 
            intensities = np.linspace(0.1, 0.9, len(expset.experiments))
            for (j,exp) in enumerate(expset.experiments):
                (ts,vsoma) = exp.get_voltage_trace()
                print 'j=%g, intensities[i]='%j,intensities[j],' .................',
                print exp
                pylab.plot(ts,vsoma,lw=rcsettings['lw_fine'],c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[j]),alpha=rcsettings['alpha'])
            self.__turn_off_border(pylab.gca())
            pylab.xlabel('time (ms)')
            pylab.ylabel('v_soma (mV)')
            
            
            if savefig is not None:
                figname = '%s_voltage_trace_%s_expset%g.png'%(savefig,time.strftime('%y%m%d'),i)
                pylab.savefig(figname)
                print 'Saved figure as %s'%figname

    def plot_IClamp_FR(self,experimentsets,savefig=None,settings='paper'):
        """
        
        """
        #TODO: work out what I was trying to do here
        rcsettings = RC_SETTINGS[settings]
        for (i,expset) in enumerate(experimentsets):
            pylab.close('all')
            pylab.figure() 
            intensities = np.linspace(0.1, 0.9, len(expset.experiments))
            for (j,exp) in enumerate(expset.experiments):
                (ts,vsoma) = exp.get_voltage_trace()
                print 'j=%g, intensities[i]='%j,intensities[j],' .................',
                print exp
                pylab.plot(ts,vsoma,lw=rcsettings['lw_fine'],c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[j]),alpha=rcsettings['alpha'])
            pylab.xlabel('time (ms)')
            pylab.ylabel('v_soma (mV)')
            self.__turn_off_border(pylab.gca())
            
            if savefig is not None:
                figname = '%s_voltage_trace_%s_expset%g.png'%(savefig,time.strftime('%y%m%d'),i)
                pylab.savefig(figname)
                print 'Saved figure as %s'%figname
        
        
        
        
    def plot_IClamp_iPhoto(self,experimentsets,savefig=None,settings='paper'):
        """
        
        
        """
        rcsettings = RC_SETTINGS[settings]
        for (i,expset) in enumerate(experimentsets):
            
            ies = expset.expvariables
            
            for opsintype in OPSINS_IPHOTO:
                pylab.close('all')
                fig = pylab.figure() 
                intensities = np.linspace(0.1, 0.9, len(expset.experiments))

                try:
                    for (j,exp) in enumerate(expset.experiments):
                        
                        
                        
                        (ts,iphoto) = exp.get_iPhoto_opsin(opsintype,entire_trace=False)
                        print ts
                        print 'j=%g, intensities[i]='%j,intensities[j],' .................',
                        print exp
                        ax = fig.gca()
                        ax.plot(ts,iphoto,lw=rcsettings['lw_fine'],c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[j]),alpha=rcsettings['alpha'])
                except:
                    pylab.close(fig)
                    continue
                self.__turn_off_border(ax)
                ax.set_xlabel('time (ms)')
                ax.set_ylim(auto=True)
                ax.set_ylabel('I_%s (nA)'%opsintype)
                
                pylab.gcf()
                inset_ax = pylab.axes([.2, .7, .15, .15])
                self._plot_injected_current(ies, rcsettings['cmaps'][self.cmap_index], inset_ax)
                
                if savefig is not None:
                    figname = '%s_iphoto_trace_%s_%s_expset%g.png'%(savefig,time.strftime('%y%m%d'),opsintype,i)
                    fig.savefig(figname)
                    print 'Saved figure as %s'%figname
                    
                    
                    
    def plot_IClamp_iPhoto_peak(self,experimentsets,savefig=None,settings='paper'):
        """
        
        
        """
        rcsettings = RC_SETTINGS[settings]
        print 'In plot_IClamp_iPhoto_peak'
        for (i,expset) in enumerate(experimentsets):
            
            ies = expset.expvariables
            
            for opsintype in OPSINS_IPHOTO:
                print opsintype,OPSINS_IPHOTO[opsintype]
                #ies = []
                iphotos = []
                intensities = np.linspace(0.1, 0.9, len(expset.experiments))
                try:
                #if opsintype == 'NpHR':
                    for (j,exp) in enumerate(expset.experiments):
                        iphoto = exp.calc_peak_iPhoto_opsin(opsintype,OPSINS_IPHOTO[opsintype])
                        #ie = intensities[j]
                        iphotos.append(iphoto)
                        #ies.append(ie)
                        print 'j=%g, intensities[i]='%j,intensities[j],' .................',
                        print exp
                except:
                    continue
                pylab.close('all')
                fig = pylab.figure() 
                ax = fig.gca()
                ax.plot(ies,iphotos,lw=rcsettings['lw'],c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[j]))
                self.__turn_off_border(ax)
                ax.set_xlabel('I_e (nA)')
                ax.set_ylabel('I_%s (nA)'%opsintype)
                
                
                if savefig is not None:
                    figname = '%s_iphoto_peak_%s_%s_expset%g.png'%(savefig,opsintype,time.strftime('%y%m%d'),i)
                    fig.savefig(figname)
                    print 'Saved figure as %s'%figname
                    
                                        
    def _plot_injected_current(self,amp_range,colormap,axes):
        intensities = np.linspace(0.1, 0.9, len(amp_range)) #TODO rescale so that intensities are between 0.1 and 0.9
        xvalues = [0,0.99,1,3,3.01,4]
        for (j,amp) in enumerate(amp_range):
            axes.plot(xvalues,[0,0,amp,amp,0,0],c=self.__get_color(colormap, intensities[j]))
        
        axes.set_title('Input I_e')
        #axes.set_yticks([amp_range[0],amp_range[-1]])
        self.__turn_off_border(axes,turnoffs=['right','top','bottom'])
        axes.set_xticks([])

    def plot_PSTH_population(self,experimentsets,savefig=None,settings='paper',binwidth=50.,tstop=1250):
        rcsettings = RC_SETTINGS[settings]
        print 'In plot_PSTH_population'
        bins =  np.arange(0,tstop,binwidth)
        if bins[-1] < tstop:
            bins = np.append(bins,tstop)
        
        print "There are %g experiment sets ... "%len(experimentsets)
        
        for (i,expset) in enumerate(experimentsets):
            self.print_header('PSTH', expset)
            ies = expset.expvariables
            spikes = []
            count = 0
            print "<---------------------------===============> There are %g experiments in this set, #%g ... "%(len(expset.experiments),i)
            for exp in expset.experiments:
                
                sp = exp.get_spikes()
                print sp ,'-'*40
                if sp is None or sp.size==0:
                    continue
                spikes += list(sp.flatten())
                count += 1
                
                
            if count==0:
                print "Count =0"
                continue
            if len(spikes)==0:
                print "sp = None =0"
                continue
            print spikes
            print count
            pylab.close('all')
            fig = pylab.figure() 
            ax = fig.gca()
            ax.hist(spikes,bins=bins)
            pylab.title('PSTH for %s (n=%g)'%('experiment',count))
            self.__turn_off_border(ax)
            
            if savefig is not None:
                figname = '%s_%s_PSTH_%s.png'%(savefig,str(expset),time.strftime('%y%m%d'))
                fig.savefig(figname)
                print 'Saved figure as %s'%figname
            
            
    def plot_population_raster(self,experimentsets,savefig=None,settings='paper',marker='|',align=None,offset=100,tmax=1500):
        """
        
        Params:
            experimentsets
            savefig            if not None, then the name of the file
            settings           rcsettings for paper or poster
            marker             marker to be used to denote a spike
            align              list of values by which each experiment should be aligned i.e. when an event happened for each trial 
            offset             value to which the events should be aligned. Only included if align is not None
        """
        rcsettings = RC_SETTINGS[settings]
        print 'In plot_population_raster'
           
        
        for (i,expset) in enumerate(experimentsets):
            self.print_header('population_raster', expset)
            spikes = []
            labels = []
            ys = []
            count = 0
            for exp in expset.experiments:
               
                sp = exp.get_spikes()
                print sp ,'-'*40
                if sp is None or sp.size==0:
                    continue
                spikes.append(sp.flatten())
                labels.append(exp.get_label())
                ys.append(np.ones(sp.size)*count)
                count += 1

            if count==0:
                continue
            if len(spikes)==0:
                continue
            pylab.close('all')
            fig = pylab.figure()
            ax = fig.gca()
            for c in range(count):
                if align is None:
                    ax.scatter(spikes[c],ys[c],marker=marker)
                else:
                    ax.scatter(spikes[c]-align[c]+offset,ys[c],marker=marker)
            
            self.__turn_off_border(ax)
            ax.set_xlabel('time (ms)')
            ax.set_yticks(range(count))
            ax.set_yticklabels(labels)
            pylab.xlim((0,tmax))
            #ax.set_ylabel('I_%s (nA)'%opsintype)
            if savefig is not None:
                figname = '%s_%s_raster_%s.png'%(savefig,str(expset),time.strftime('%y%m%d'))
                fig.savefig(figname)
                print 'Saved figure as %s'%figname
            
            
    def plot_population_voltage(self,experimentsets,savefig=None,settings='paper',align=None,offset=100,tmax=1500,plotMean=False,colors=[],cmap=None,colorMean='k',showLegend=True,useCBar=False):
        """
        
        Params:
            experimentsets
            savefig            if not None, then the name of the file
            settings           rcsettings for paper or poster
            align              list of values by which each experiment should be aligned i.e. when an event happened for each trial 
            offset             value to which the events should be aligned. Only included if align is not None
        """
        rcsettings = RC_SETTINGS[settings]
        for (i,expset) in enumerate(experimentsets):
            
            self.print_header('population_voltage', expset)
            
            fig = pylab.figure()
            ax = fig.gca()
            
            for (c,exp) in enumerate(expset.experiments):
                data = exp.get_voltage_trace()
                ts,vs = data[0],data[1]
                
                if align is None:
                    ax.plot(ts,vs,label=exp.get_label())
                else:        
                    ax.plot(ts-align[c]+offset,vs,label=exp.get_label())
                    
                if plotMean:
                    #TODO implement plotMean
                    pass
                
            self.__turn_off_border(ax)
            ax.set_xlabel('time (ms)')
            ax.set_ylabel('mV')
            pylab.xlim((0,tmax))
            if savefig is not None:
                figname = '%s_%s_voltage_%s.png'%(savefig,str(expset),time.strftime('%y%m%d'))
                fig.savefig(figname)
                print 'Saved figure as %s'%figname
            
    def plot_compare_features(self,experimentsets,savefig=None,settings='paper',traits=[],drawMean=True,marker='o'):
        
        rcsettings = RC_SETTINGS[settings]
        allresults = {}
        labels = []
        #setup
        for t in traits:
            allresults[t] = {}
        
        #get results --> our data structure
        for (i,expset) in enumerate(experimentsets):
            print '='*40
            print expset
            print expset.get_label()
            print '='*40
            labels.append(expset.get_label())
            #setup
            for t in traits:
                allresults[t][expset.get_label()] = []
            
            for (j,exp) in enumerate(expset.experiments):
                results = exp.get_results()
                print results
                if len(results.items())==0:
                    continue
                
                for t in traits:
                    if not results.has_key(t):
                        continue
                    allresults[t][expset.get_label()].append(results[t])
         
        
        #plot results
        for t in traits:
            pylab.close('all')
            fig = pylab.figure()
            ax = fig.gca()
            
            results = allresults[t]
            expsets = results.keys()
            expsets.sort()
            means = []
            stds =  []
            
            for expset in expsets:    
                rr = np.array( results[expset])
                means.append(rr.mean())
                stds.append(rr.std())
            
            if drawMean:
                fmt = '-'
            ax.errorbar(range(len(expsets)),means,yerr=stds,fmt=fmt+marker)
            
            self.__turn_off_border(ax)
            ax.set_xticks(range(len(expsets)))
            ax.set_xticklabels(labels)
            pylab.xlim((-1,len(labels)+1))
            
            if savefig is not None:
                figname = '%s_%s_trait%s_%s.png'%(savefig,str(expset),t,time.strftime('%y%m%d'))
                fig.savefig(figname)
                print 'Saved figure as %s'%figname
    
    
    
    def plot_compare_2D_FI(self,experimentsets,savefig=None,settings='paper',axeslabels=[],axestitles=[],firing_rate='FI',plottitle=None):
        """
        Plot 2D histogram of firing rate (specified as FI) against two axes for two different traits
        
        #TODO set transparent for FI==0
        #TODO set vmax to be mod 5
        """
        rcsettings = RC_SETTINGS[settings]
        pylab.rcParams.update(rcsettings)
        pylab.close('all')
        
        for (i,expset) in enumerate(experimentsets):

            #get responses for 2D data set
            data = np.array(expset.results[firing_rate])
            print data
            print expset.expvariables[0], expset.expvariables[1]
            print len(data)
            print len(expset.expvariables[0]), len(expset.expvariables[1])
            print len(expset.expvariables[0]) * len(expset.expvariables[1])
            data = data.reshape(len(expset.expvariables[0]),len(expset.expvariables[1]))
            print data 
            
            #plot dataset
            pylab.close('all')
            pylab.figure()#figsize=(5,5))
            pylab.imshow(data,interpolation='nearest',cmap=self.cmap)
            
            ax = pylab.gca()
            self.__turn_off_border(ax)
            # axes labels and titles
            ax.set_xlabel(axestitles[0])
            ax.set_ylabel(axestitles[1])
            ax.set_xticks(range(len(expset.expvariables[0])))
            ax.set_yticks(range(len(expset.expvariables[1])))
            if len(axeslabels)==0:
                ax.set_xticklabels(expset.expvariables[0])
                ax.set_yticklabels(expset.expvariables[1])
            else:
                ax.set_xticklabels(axeslabels[0])
                ax.set_yticklabels(axeslabels[1])
            
            cb = pylab.colorbar()
            cb.ax.set_ylabel('%s (Hz)'%firing_rate)
            
            if plottitle is not None:
                pylab.title(plottitle)
            
            #save plot
            if savefig is not None:
                figname = '%s_%s_2DcompareFI_%s.png'%(savefig,expset,time.strftime('%y%m%d'))
                pylab.savefig(figname,dpi=rcsettings['dpi'])
                print 'Saved figure as %s'%figname
            pylab.close('all')
            
        
        
    
    """   
        
    def _old_plot_fit_FI(self,experimentsets,savefig=None,settings='poster',polyfn='_poly',p0=None):
        rcsettings = RC_SETTINGS[settings]
        pylab.rcParams.update(rcsettings)
        poly = getattr(self,polyfn)
        pylab.figure() 
        intensities = np.linspace(0.1, 0.9, len(experimentsets))
        from matplotlib import mpl
        norm = mpl.colors.BoundaryNorm(intensities, len(intensities))
        
        #cNorm = colors.Normalize(vmin=0,vmax=1) 
        #scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)
        
        m = cm.ScalarMappable(cmap=rcsettings['cmaps'][self.cmap_index])
        m.set_array(intensities)

        for (i,expset) in enumerate(experimentsets):
            try:
            #if True:
                print 'i=%g, intensities[i]='%i,intensities[i]
                print 'exp variables len', len(expset.expvariables[0])
                print expset.expvariables[0]
                print expset.results
                print 'length FI', len(expset.results['FI'])
                xvals = np.linspace(expset.expvariables[0][0], expset.expvariables[0][-1],num=len(expset.expvariables[0]))
                svals = np.asarray(xvals).ravel() 
                fis = np.asarray(expset.results['FI']).ravel()
                print fis
                if p0 is None:
                    p0 = self._get_guesstimate(expset.expvariables[0],expset.results['FI'],polyfn)
                popt_curve, pcov_curve = curve_fit(poly, svals, fis,p0)
                fitted = poly(svals,popt_curve[0],popt_curve[1],popt_curve[2],popt_curve[3])
                #scatter real results, plot trend
                pylab.scatter(expset.expvariables[0],expset.results['FI'],lw=rcsettings['lw']/2,c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),s=rcsettings['lw']*10,facecolor='none',edgecolor=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),alpha=0.6)
                pylab.plot(svals,fitted,lw=rcsettings['lw'],c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),label=expset.get_label())
            except:
                print 'Error with plotting for expset', expset
                print "Unexpected error:", sys.exc_info()[0]
                continue
            
        # turn the right and top axes off
        self.__turn_off_border(pylab.gca())
        
        pylab.ylabel('Output frequency (Hz)')
        pylab.xlabel('Input')
        
        # TODO: add vertical line at x = 0
        
        
        # increase scale so that we have some space to place the legend
        #pylab.ylim(ymax=pylab.ylim()[1]*1.2)
        pylab.xlim(expset.expvariables[0][0]-0.2,expset.expvariables[0][-1]+0.2)
        
        cbarlabels = [expset.get_label() for expset in experimentsets]
        dc = intensities[1]-intensities[0]
        mapint = np.linspace(0.1-dc/2, 0.9+dc/2, len(experimentsets)+1)
        cb = pylab.colorbar(m,ticks=intensities,boundaries=mapint)
        cb.set_label('x NpHR Factor')
        cb.set_ticklabels(cbarlabels)
        
        
        if savefig is not None:
            figname = '%sFIfit_%s_%s.png'%(savefig,time.strftime('%y%m%d'),polyfn)
            pylab.savefig(figname,dpi=rcsettings['dpi'])
            print 'Saved figure as %s'%figname        

        # Prettify the legend
        leg = pylab.legend(loc=2,fancybox=True)
        leg.get_frame().set_alpha(0.5) 
        
        if savefig is not None:
            figname = '%sFIfit_%s_%s_legend.png'%(savefig,time.strftime('%y%m%d'),polyfn)
            pylab.savefig(figname,dpi=rcsettings['dpi'])
            print 'Saved figure as %s'%figname        
        
                
    
     
    def _old_plot_fit_FI_bg(self,experimentsets,savefig=None,settings='poster',polyfn='_poly',p0=None):
        
        rcsettings = RC_SETTINGS[settings]
        pylab.rcParams.update(rcsettings)
        poly = getattr(self,polyfn)
        pylab.figure() 
        intensities = np.linspace(0.1, 0.9, len(experimentsets))
        from matplotlib import mpl
        norm = mpl.colors.BoundaryNorm(intensities, len(intensities))
        
        #cNorm = colors.Normalize(vmin=0,vmax=1) 
        #scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)
        
        m = cm.ScalarMappable(cmap=rcsettings['cmaps'][self.cmap_index])
        m.set_array(intensities)

        for (i,expset) in enumerate(experimentsets):
            try:
                print 'i=%g, intensities[i]='%i,intensities[i]
                print 'exp variables len', len(expset.expvariables[0])
                print expset.expvariables[0]
                print 'length FI', len(expset.results['FI'])
                xvals = expset.results['FI_bg']
                svals = np.asarray(xvals).ravel()
                fis = np.asarray(expset.results['FI']).ravel()
                svals.sort()
                fis.sort()
                print fis
                if p0 is None:
                    p0 = self._get_guesstimate(expset.results['FI_bg'],expset.results['FI'],polyfn)
                print "About to fit (old):", svals, fis, p0
                popt_curve, pcov_curve = curve_fit(poly, svals, fis,p0)
                fitted = poly(svals,popt_curve[0],popt_curve[1],popt_curve[2],popt_curve[3])
                #scatter real results, plot trend
                pylab.scatter(expset.results['FI_bg'],expset.results['FI'],lw=rcsettings['lw']/2,c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),s=rcsettings['lw']*10,facecolor='none',edgecolor=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),alpha=0.6)
                pylab.plot(svals,fitted,lw=rcsettings['lw'],c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),label=expset.get_label())
            except RuntimeError:
                pylab.scatter(expset.results['FI_bg'],expset.results['FI'],lw=rcsettings['lw']/2,c=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),s=rcsettings['lw']*10,facecolor='none',edgecolor=self.__get_color(rcsettings['cmaps'][self.cmap_index],intensities[i]),alpha=0.6)
                print 'Error with plotting for expset', expset
                print "Unexpected error:", sys.exc_info()[0]
                continue
            except:
                print 'Error with plotting for expset', expset
                print "Unexpected error:", sys.exc_info()[0]
                continue
        # turn the right and top axes off
        self.__turn_off_border(pylab.gca())
        
        pylab.ylabel('Output frequency (Hz)')
        pylab.xlabel('Input')
        
        # TODO: add vertical line at x = 0
        
        
        # increase scale so that we have some space to place the legend
        #pylab.ylim(ymax=pylab.ylim()[1]*1.2)
        #pylab.xlim(expset.expvariables[0][0]-0.2,expset.expvariables[0][-1]+0.2)
        pylab.xlim(xmin=-2)
        
        cbarlabels = [expset.get_label() for expset in experimentsets]
        dc = intensities[1]-intensities[0]
        mapint = np.linspace(0.1-dc/2, 0.9+dc/2, len(experimentsets)+1)
        cb = pylab.colorbar(m,ticks=intensities,boundaries=mapint)#mapint)
        cb.set_label('x NpHR Factor')
        cb.set_ticklabels(cbarlabels)
        
        
        if savefig is not None:
            figname = '%sFIfit_bg_%s_%s.png'%(savefig,time.strftime('%y%m%d'),polyfn)
            pylab.savefig(figname,dpi=rcsettings['dpi'])
            print 'Saved figure as %s'%figname        

        # Prettify the legend
        leg = pylab.legend(loc=2,fancybox=True)
        leg.get_frame().set_alpha(0.5) 
        
        if savefig is not None:
            figname = '%sFIfit_bg_%s_%s_legend.png'%(savefig,time.strftime('%y%m%d'),polyfn)
            pylab.savefig(figname,dpi=rcsettings['dpi'])
            print 'Saved figure as %s'%figname        
            
    
    """    
        



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
        self.expvariables = []
        self.results = {} 
        self.experiments = []
    
    def __str__(self):
        return self.basename+self.subselect

    def get_label(self):
        return self.label
    
    def get_property(self,key,index=None):
        """ Catch for whether we don't know if it's a label or a dictionary label """
        print self.results.keys()
        # try results dictionary first
        if self.results.has_key(key):
            return self.results[key]
        else: # it might be a class property
            try:
                if index is None:
                    return getattr(self,key)
                elif type(index) is int:
                    return getattr(self,key)[index]
            except:
                print 'No property %s associated with this experiment'%key
                return None
            
    
    def populate_experiments(self,exp_params,analysis_params):
        print '='*40
        print exp_params
        print analysis_params

        if type(exp_params) is not list:
            print 'Error'
            return

        if type(exp_params[0]) is list or type(exp_params[0]) is np.ndarray:
            self.expvariables = [list(exp) for exp in exp_params]
        else:
            self.expvariables = exp_params
        
        # Extended for the case where exp_params is lists of lists by using itertools.product 
        print self.expvariables, '<---------------------'
        try:
            for v in list(itertools.product(*self.expvariables)):
                print v
                print self.subselect
                print self.subselect%v
                self.experiments.append(Experiment(self.basename ,self.subselect%v,analysis_params,vals=v))
        except:
            for v in list(itertools.product(self.expvariables)):
                print v
                print self.subselect
                print self.subselect%v
                self.experiments.append(Experiment(self.basename ,self.subselect%v,analysis_params,vals=v))
    
    def save_experiments(self):
        print 'Saving ExperimentSet', self.label
        for exp in self.experiments:
            exp.save_results()
        
    def load_experiments(self):
        print 'Loading ExperimentSet', self.label
        for exp in self.experiments:
            # will fail if doesn't get results for one experiment ->
            #TODO: implement more robust method
            exp.load_results()
            
            #once data is loaded, go through and repopulate self.
            for (k,v) in exp.results.iteritems():
                print k, v
                try:
                    self.results[k].append(v)
                except: #doesn't have key
                    self.results[k] = [v]
        
        
    def run_analysis(self,analysis_fns,recalc=False):
        print 'in expset = ', analysis_fns
        for exp in self.experiments:
            print 'About to run analysis for ',str(exp)
            exp.calculate_results(analysis_fns, recalc)
            exp.save_results()
            
        # calculate experimentset values
        #self.get_datafit()
        
    def get_datafit(self):
        # Try all type of curve to get best fit for injected and bg
        self.results['MI'] = self._get_modulationIndex('expvariables','FI')
        try:
            self.results['MI_bg'] = self._get_modulationIndex('FI_bg','FI')
        except:
            pass
        
    
    
    def _find_nearest(self,array,value):
        return (np.abs(array-value)).argmin()
        
    
    def _get_modulationIndex(self,fit_key,yval_key,eval_at_x=None):
        
        ep = ExperimentPlotter()
        
        ys = self.get_property(yval_key)
        
        if fit_key == 'expvariables':
            xs = self.expvariables[0]
            xvals = np.linspace(xs[0], xs[-1],num=len(xs)) 
        else:
            xs = self.get_property(fit_key)
            # as this is usually an indirect measure i.e. background, it may not be sorted - 
            # so sort ALONG with the corresponding y-values
            xs,ys = ep._sort_xsys(xs,ys)
            xvals = xs
            
        svals = np.asarray(xvals).ravel()
        fis = np.asarray(ys).ravel()
        
        fitted_xvals = np.linspace(xvals[0], xvals[-1],num=len(xs)*4) # add more data points
        fitted_svals = np.asarray(fitted_xvals).ravel() 
        
        # find the best fit
        poly,popt_curve,pname = ep._fit_curve(svals, fis,p0=None,polyfn=['_sigmoid','_poly'])
        #fittedys = poly(fitted_svals,popt_curve[0],popt_curve[1],popt_curve[2],popt_curve[3])
        # grab the derivative of the best fit
        d_poly = ep.get_derivative(pname)
        fitted_dys = d_poly(fitted_svals,popt_curve[0],popt_curve[1],popt_curve[2],popt_curve[3])
        """
        pylab.figure()
        pylab.plot(fitted_svals,fittedys)
        pylab.scatter(svals,fis)
        pylab.savefig('original_%s_%g_%g_%g_%g.png'%(pname,popt_curve[0],popt_curve[1],popt_curve[2],popt_curve[3]))
        
        pylab.figure()
        pylab.plot(fitted_svals,fitted_dys)
        pylab.savefig('derivative_%s_%g_%g_%g_%g.png'%(pname,popt_curve[0],popt_curve[1],popt_curve[2],popt_curve[3]))
        pylab.close('all')
        
        print type(fitted_dys),'------------------------------- %s when calculated on xs[0], xs[-1]=('%(pname),xvals[0], xvals[-1],') are',[fitted_dys.min(),fitted_dys.max()]
        print "fittedys = ",fitted_dys.min(),fitted_dys.max(),fittedys
        """
        # and find the min, max values
        if eval_at_x is None:
            return [fitted_dys.min(),fitted_dys.max()]
        else:
            # TODO: find x index and corresponding t
            return fitted_dys[self._find_nearest(fitted_svals, eval_at_x)]
    
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
                print 'looking up ',FUNCTION_LOOKUP[key]
                try:
                    tmp_v = getattr(exp, FUNCTION_LOOKUP[key])(key)
                except:
                    tmp_v = getattr(exp, FUNCTION_LOOKUP[key])()
                print tmp_v, '-------------------------------------------------------'
                tmp_results.append(tmp_v)
            self.results[key] = tmp_results
        #except:
        #    print "Could not calculate responses, unknown variable value = ", key
        #print 'Current results ',self.results
    
    def print_exps(self,indent=4):
        print self.basename,self.subselect, 'contains experiments:'
        for exp in self.experiments:
            print ' '*indent, exp 
            print ' '*indent, exp.results
            
    
    def print_data(self,indent=4):
        print self.basename,self.subselect, 'contains data:'
        print 'label = ',self.label
        for val in self.expvariables:
            print ' '*indent, '%s'%(val)
        for (k,v) in self.results.iteritems():
            print ' '*indent, '%s = %s'%(k,v)
            
class RangeTrials:
    
    def __init__(self,basenames,trialInstances,variables,trials,trialLabel,label='no_label',var_format="%g"):
        print "0 trialLabels = ",trialLabel
        self.expRange = []
        for (i,instance) in enumerate(trialInstances):
            if type(basenames) is list:
                basename = basenames[i]
            else:
                basename = basenames
            if type(trialLabel) is list:
                    
                testlabel = trialLabel[i]
            else:
                testlabel = "%s_%g"%(trialLabel,i)   
            print "A trialLabels = ",trialLabel
            self.expRange.append(TrialExperimentSet(basename,[instance],variables,trials,testlabel,label,var_format))
            
    def load_experiments(self):
        for tes in self.expRange:
            tes.load_experiments()
        
    def collate_results(self):
        for tes in self.expRange:
            tes.collate_results()
    
    def calculate_responses(self,key,recalc=False):
        for tes in self.expRange:
            tes.calculate_responses(key,recalc=False)

class TrialExperimentSet:
    
    def __init__(self,basename,trialInstances,variables,trials,trialLabel,label='no_label',var_format="%g",analysis_params={}):
        self.experimentsets = [] 
        print "B trialLabels = ",trialLabel
        self.__createTrialset(basename, trialInstances,variables, trials, trialLabel, label,var_format,analysis_params)
        self.results = {}
        

        
    def __createTrialset(self,basename,trialInstances,variables,trials,trialLabel,labels,var_format,analysis_params):
        print "C trialLabels = ",trialLabel
        self.label = trialLabel
        print labels, type(labels)
        print self.label, "=-======"
        self.expvariables = variables
        for (i,instance) in enumerate(trialInstances):
            print("instance = %s ------------------------------------------"%instance)
            for (j,var) in enumerate(variables):

                if type(labels) is list:
                    
                    label = labels[j]
                else:
                    label = "%s_%g"%(labels,j)   
                print("instance --> "+instance%(var,"%g"))
                expset = ExperimentSet(basename,instance%(var,"%g"),label=label)
                expset.populate_experiments(trials,analysis_params)
                self.experimentsets.append(expset)
    
    def get_label(self):
        return self.label
    
    def save_experiments(self):
        print 'Saving TrialExperimentSet', self.label
        for exp in self.experimentsets:
            exp.save_experiments()
    
    def load_experiments(self):
        
        print 'Loading TrialExperimentSet'
        for expset in self.experimentsets:
            
            expset.load_experiments()
            print expset.results
            #TODO: check that expset.label is not a reserved analysis type
            self.results[expset.label] = expset.results

        #print self.results
    
    def calculate_responses(self,key,recalc=False):
        for expset in self.experimentsets:
            expset.calculate_responses(key,recalc=False)
    
    
    def collate_results(self):
        """
        Use at your own risk! This method does NOT check that all results are present for 
        all instances / trials; it instead just collates all instances so that 
            label_freq10 : { 'FI': [10.5,11.,10.4]}
            label_freq20 : { 'FI': [23.2,21.,24.2]}
            --> 
            {'FI' : [10.5,11.,10.4,23.2,21.,24.2] }
        """
        tmp_dict = []
        for expset in self.experimentsets:
            for restype in expset.results.keys():
                if restype not in tmp_dict:
                    tmp_dict.append(restype)
        
        for restype in tmp_dict:
            self.results[restype] = []
            for expset in self.experimentsets:
                if self.results[expset.label].has_key(restype):
                    self.results[restype] = self.results[restype]+self.results[expset.label][restype]
                
                
                
    def print_exps(self,indent=4):
        for exp in self.experimentsets:
            exp.print_exps(indent)
            
    
    def print_data(self,indent=4):
        for exp in self.experimentsets:
            exp.print_data(indent)

        
        
    def calculate_average_responses(self):
        avgs = ['FI','FI_bg','FI_post']
        for expset in self.experimentsets:
            for restype in expset.results.keys():
                 
                #TODO: implement
                pass
        
    def get_property(self,key,index=None):
        """
        Catch for whether we don't know if it's a label or a dictionary label 
        """
        print self.results.keys()
        # try results dictionary first
        if self.results.has_key(key):
            return self.results[key]
        else: # it might be a class property
            try:
                if index is None:
                    return getattr(self,key)
                elif type(index) is int:
                    return getattr(self,key)[index]
            except:
                print 'No property %s associated with this experiment'%key
                return None
        
        #TODO: reimplement this, so that property for all trials is returned - and that should be enough
        

class Experiment:
    
    NEUROTOOLS_METHODS = ['isi','cv_isi','cv_kl','fano_factor_isi','mean_rate']
    
    def __init__(self,basename,expname,analysis_params={},vals={},tryload=True):
        self.basename = basename
        self.expname = expname
        self.fullname = ''
        self.results = {}
        self.analysis_params = analysis_params
        if not self.analysis_params.has_key('tstart'):
            self.analysis_params['tstart']=None
        if not self.analysis_params.has_key('tstop'):
            self.analysis_params['tstop']=None
        self.spikes = None
        self.label = None
        self.values = vals
        if analysis_params.has_key('label_format'):
            self.label = analysis_params['label_format']
        self.spiketrain = None
        
        if tryload:
            try:
                self.load_results()
            except:
                print 'Error with tryload in init'
                
            try:
                self.load_spikes()
            except:
                print 'Error with second tryload'
        
    def __str__(self):
        return self.basename+self.expname
    
    def get_label(self):
        print self.label
        print self.values
        
        if self.label is not None:
            return self.label%self.values
        return str(self)
    
    
    def save_results(self):
        try:
            savefile = EXP_PKL_LOCATION%self.basename +"%s_results.pkl"%(self.basename+self.expname)
            print 'Attempting to save to file:',savefile,
            output = open(savefile, 'wb')
            pkl.dump(self.results, output,-1)
            output.close()
            print '... saved successfully'
        except:
            print "Unexpected error:", sys.exc_info()[0]
        
    def load_results(self):
        try:
        #if True:
            loadfile = EXP_PKL_LOCATION%self.basename +"%s_results.pkl"%(self.basename+self.expname)
            print 'Attempting to load file:',loadfile,
            pkl_file = open(loadfile, 'rb')
            self.results = pkl.load(pkl_file)
            pkl_file.close()
            print '... loaded successfully'
        except IOError as e:
            print 'File does not exist'
            raise e
        except:
            print "Unexpected error:", sys.exc_info()[0]
    
    def save_spikes(self,spikes):
        #TODO : should check that file doesn't already exist, and if it does ...?
        
        try:    
            savefile = EXP_GDF_LOCATION%self.basename +"%s.gdf"%(self.basename+self.expname)
            print 'Attempting to save to file:',savefile,
            np.savetxt(savefile,spikes)
            print '...saved successfully'
        except:
            print "Unexpected error:", sys.exc_info()[0]
            
            
    def load_spikes(self):
        try:
        #if True:
            loadfile = EXP_GDF_LOCATION%self.basename +"%s.gdf"%(self.basename+self.expname)
            print 'Attempting to load spike file:',loadfile,
            spikes = np.loadtxt(loadfile,ndmin=1)
            print '... loaded successfully',
            self.spikes = spikes
            
            print 'and set spikes as field'
        #except TypeError,IndexError:
        #    self.spikes = np.array([])
        except IOError:
            print 'File does not exist'
        except:
            print "\nUnexpected error:", sys.exc_info()[0]
            return
        self.get_spiketrain()
    
            
    def _check_spikefile_exists(self):
        savefile = EXP_GDF_LOCATION%self.basename +"%s.gdf"%(self.basename+self.expname)
        print 'Testing for savefile', savefile
        gg = glob.glob(savefile)
        if len(gg)==1:
            return True
        else:
            return False
        
    def calculate_results(self,properties=[],recalc=False,save=True):
        """
            Assumes that self.results is already loaded
        """
        for prop in properties:
            print 'Looking at property', prop
            #print FUNCTION_LOOKUP
            if not recalc and self.results.has_key(prop):
                continue
            elif recalc:
                print 'Calculating value for %s'%prop
            elif not self.results.has_key(prop):
                print 'Value for %s missing; running %s'%(prop, FUNCTION_LOOKUP[prop])
            result = getattr(self, FUNCTION_LOOKUP[prop])(feature=prop)
            if result is not None:
                self.results[prop] = result
            #getattr(self, FUNCTION_LOOKUP[prop])()
        print self.results
        if save and self.results != {}:
            self.save_results()
    
    def get_dat_file(self,modifier=''):
        try:
            datfile = EXP_DAT_LOCATION%self.basename +"%s%s.dat"%(self.basename+self.expname,modifier)
            print 'location =',datfile
        except:
            print 'Error with resolving name when finding file'
            return 
        gg = glob.glob(datfile)
        if len(gg)==1:
            print 'Found it'                    
            self.fullname = gg[0][:-4] # name of .dat file without the file extension
            return gg[0]
        elif len(gg)==0:
            #TODO: raise error
            print 'No file found: ',datfile
            try:
                datfile = EXP_DAT_LOCATION%self.basename +"%s%s_v.dat"%(self.basename+self.expname,modifier)
                print 'trying location =',datfile
                dd = glob.glob(datfile)
                print dd
                if len(dd)==1:
                    print 'Found it'
                    self.fullname = dd[0][:-6] # name of .dat file without the file extension
                    return dd[0]
                else:
                    print 'Could not find it = '
            except:
                print 'Error = ',sys.exc_info()[0] 
                pass
        elif len(gg)>1:
            #TODO: raise error - too many files
            print 'Too many files found: ',file

    def get_voltage_trace(self):
        data = np.loadtxt(self.get_dat_file()) 
        t,v_soma = data[0,:],data[1,:]
        print data.size
        return (t,v_soma)

    def calc_peak_iPhoto(self,**kwargs):
        #if self.results['iPhoto'] is None:
        #    self.results['iPhoto'] = {}
        #self.results['iPhoto'][opsintype] = peak_iphoto
        for (opsin,values) in OPSINS_IPHOTO.iteritems():
            try:
                #TODO return correct iphoto
                return self.calc_peak_iPhoto_opsin(opsin,values)
            except:
                pass


    def calc_peak_iPhoto_opsin(self,opsintype,max_iphoto,**kwargs):
        """
            opsin_type
            max_iphoto
        
        """
        (ts,i_photo) = self.get_iPhoto_opsin(opsintype)
        if max_iphoto:
            peak_iphoto =  i_photo.max()
        else:   #no prizes for guessing correctly ... we calculate the min
            peak_iphoto = i_photo.min()
        print 'peak iphoto = ',peak_iphoto
        return peak_iphoto
    
    
    def get_iPhoto_opsin(self,opsintype,entire_trace=True):
        data = np.loadtxt(self.get_dat_file('*_i%s'%opsintype)) 
        ts = data[1:,0]
        i_photo = data[1:,1]
        print 'got here ok'
        if not entire_trace:
            #select so that it's 100 before tstart, and 150 after tstop from analysis_params
            tstart = np.where(ts>(self.analysis_params['tstart']-self.analysis_params['prebuffer']))[0][0] #TODO move 100 and 150 as buffer
            tstop = np.where(ts<(self.analysis_params['tstop']+self.analysis_params['postbuffer']))[0][-1]
            print 'tstart and tstop', tstart, 
            print tstop
            ts = ts[tstart:tstop]
            i_photo = i_photo[tstart:tstop]
            print ts
        return ts, i_photo
        
    def calc_isi(self,**params):
        return self.calc_isi_generic(self.analysis_params['tstart'],self.analysis_params['tstop'],field='isi')

    def calc_isi_bg(self,**params):
        #print 'calc_firing rate_bg',
        return self.calc_isi_generic(self.analysis_params['tstart_bg'],self.analysis_params['tstop_bg'],field='isi_bg')
    
    def calc_isi_post(self,**params):
        #print 'calc_firing rate_post',
        return self.calc_isi_generic(self.analysis_params['tstart_post'],self.analysis_params['tstop_post'],field='isi_post')
        
    def calc_isi_generic(self,tstart,tstop,field):
        try:
            if self.spikes is not None:
                spikes = self.spikes
            elif self._check_spikefile_exists():
                spikes = self.load_spikes()
            else:
                data = np.loadtxt(self.get_dat_file()) 
                t,v_soma = data[0,:],data[1,:]
                spikes = self.extract_spiketimes(t,v_soma)
            if len(spikes)==0:
                print 'No spikes observed during entire duration'
                self.results[field] = 0
                return 0
            spikes = np.where((spikes >= tstart) & (spikes <= tstop))[0]
            
            #isi = 
            self.get_spiketrain(tstart,tstop)
            # set isi for bg/post/stim
            self.results[field] = self.calc_feature('isi')
            print('Generated ISI')
            # and while we're here, also set it for cvisi and fano factor
            self.results['cv_'+field] = self.calc_feature('cv_isi')
            self.results['ff_'+field] = self.calc_feature('fano_factor_isi')
            print('Generated CV and FF for ISI')
        except:
            print "Unexpected error:", sys.exc_info()[0]
            print '-'*40
            return -1



    def calc_firingrate(self,**params):
        """
            expname     experiment name
            t_cutout    list of start time and finish times during which firing rate should be calculated
        """
        
        #print 'calc_firing rate',
        return self.calc_firingrate_generic(self.analysis_params['tstart'],self.analysis_params['tstop'],field='FI')
    
    def calc_firingrate_bg(self,**params):
        #print 'calc_firing rate_bg',
        return self.calc_firingrate_generic(self.analysis_params['tstart_bg'],self.analysis_params['tstop_bg'],field='FI_bg')
    
    def calc_firingrate_post(self,**params):
        #print 'calc_firing rate_post',
        return self.calc_firingrate_generic(self.analysis_params['tstart_post'],self.analysis_params['tstop_post'],field='FI_post')
        
    def calc_firingrate_generic(self,tstart,tstop,field):
        try:
            if self.spikes is not None:
                spikes = self.spikes
            elif self._check_spikefile_exists():
                spikes = self.load_spikes()
            else:
                data = np.loadtxt(self.get_dat_file()) 
                t,v_soma = data[0,:],data[1,:]
                spikes = self.extract_spiketimes(t,v_soma)
            if len(spikes)==0:
                print 'No spikes observed during entire duration'
                self.results[field] = 0
                return 0
            spikes = np.where((spikes >= tstart) & (spikes <= tstop))[0]
            
            spikerate = 1000.*len(spikes)/(tstop-tstart)
            self.results[field] = spikerate
            #print spikerate
            return spikerate
        except:
            print "Unexpected error:", sys.exc_info()[0]
            print '-'*40
            return -1
        
    def check_spikes(self):
        self.get_spikes()
        self.get_spiketrain()

    def calc_feature(self,feature,**kwargs):
        
        self.check_spikes()
        if self.spikes is None:
            print 'No spike train found'
            return
        if self.spikes.size == 0:
            print 'Emptry spike train'
            return
        
        if self.NEUROTOOLS_METHODS.count(feature)>0:
            try:
                value = getattr(self.spiketrain,feature)()
            except:
                print "Unexpected error:", sys.exc_info()[0]
                return None
            
            print '%s = '%feature,value
            return value
        else:
            print "Characteristic '%' is currently unsupported"%feature
            return None


    def get_spiketrain(self,tstart=None,tstop=None):
        if tstart is None:
            tstart = self.analysis_params['tstart']
        if tstop is None:
            tstop = self.analysis_params['tstop']
        
        if self.spiketrain is not None:
            pass
        elif self.spiketrain is None and self.spikes is not None:
            self.spiketrain = NeuroTools.signals.SpikeTrain(self.spikes,tstart,tstop)
        elif self.spikes is None:
            print 'Could not create spiketrain'
        return self.spiketrain
        
    
        
    def get_spikes(self):
        if self.spikes is not None:
            pass
        elif self._check_spikefile_exists():
            self.load_spikes()
        else:
            spikes = self.process_spiketimes()
            self.save_spikes(spikes)
            self.spikes = spikes
        
        return self.spikes
    
    def get_results(self):
        if len(self.results.items())==0:
            pass
        else:
            try:
                self.load_results()
            except:
                print "Couldn't load results"
                print sys.exc_info()
        print 'self.results = ',self.results
        return self.results
        
        
    def process_spiketimes(self,**params):
        
        
        #try:
        if True:
            print 'Am going to check if spike file exists ...'
            if self._check_spikefile_exists():
                print 'File already exists'
                return
            datfile = self.get_dat_file()
            print datfile
            data = np.loadtxt(datfile) 
            t,v_soma = data[0,:],data[1,:]
            spikes = self.extract_spiketimes(t,v_soma)
            # we want to have all spikes, not just the ones in our area of analysis
            self.save_spikes(spikes)
            return spikes
        #except:
            print "Unexpected error:", sys.exc_info()[0]
            #pass #
    
    def extract_spiketimes(self,times,v_soma):
        """
        This section is a rewritten version of NeuroTools
        
        """
        above = np.where(v_soma > self.analysis_params['v_th'])[0]
        print 'spikes = ',above
        if len(above) <= 0:
            return []
        else:
            take = np.where(np.diff(above)>1)[0] + 1
            take = np.append(0,take)
            
            return times[above][take]
        

