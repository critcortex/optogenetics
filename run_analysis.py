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

DEFAULT_FREQ = [0.5,1,1.5,2,2.5,4,5,6,7,8,9,10,15,20,30,35,40,45,50]
DEFAULT_IRRAD = [1,2,5,10,20,50,100]
DEFAULT_PARAMS = {'FREQ':DEFAULT_FREQ,'IRRAD':DEFAULT_IRRAD}
DEFAULT_ANALYSIS_SETTINGS = {'v_th': -20., 'tstart':200,'tstop':2000,'prebuffer':300,'postbuffer':150}

FUNCTION_LOOKUP = {'isi':'calc_feature',
                   'cv_isi':'calc_feature',
                   'cv_kl':'calc_feature',
                   'fano_factor_isi':'calc_feature',
                   'mean_rate':'calc_feature',
                   'spiketimes':'process_spiketimes',
                   'FI'     :'calc_firingrate',
                   'iPhoto' : 'calc_peak_iPhoto'}

EXP_DAT_LOCATION = 'experiments/%s/dat/'
EXP_IMG_LOCATION = 'experiments/%s/img/'
EXP_PKL_LOCATION = 'experiments/%s/pkl/'
EXP_GDF_LOCATION = 'experiments/%s/gdf/'

OPSINS_IPHOTO = {'ChR': False, 
             'NpHR': True,
             'ArchT': True}

RC_SETTINGS = {'paper': {'lw': 3, 
                         'lw_fine':1,
                         'alpha':0.5,
                         'dpi':100,
                         'colors':['r','b','k','g','y','r','b','k','g','y'],
                         'cmaps':['YlOrRd','GnBu','jet']},
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
                          'cmaps':['YlOrRd','GnBu','jet']}}



class AnalyticFrame:
    
    def __init__(self):
        self.all_experiments = {} #TODO: include all experiments
        self.experimentset = []
        self.expplotter = ExperimentPlotter()
        self.analysis_params = DEFAULT_ANALYSIS_SETTINGS
        
        
    #def add_experimentset(self,expset):
    #    self.experimentset.append(expset)
    
    def update_params(self,newparams):
        self.analysis_params.update(newparams)
        print 'Updated parameters. New parameters are:',self.analysis_params

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
            
    def perform_analysis(self,analysis_features=[]):
        #for feature in analysis_features:
        #    print feature,'='*60
        for exps in self.experimentset:
            exps.run_analysis(analysis_features)
        
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
    
        kk = ['FI','IClamp','IClamp_iPhoto','IClamp_iPhoto_peak','FI_compare','fit_FI','PSTH_population','population_raster','population_voltage','compare_features'] 
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
class ExperimentPlotter:
    
    def print_header(self,plottype,expset):
        print '='*40
        print 'Creating %s plot for expset = %s'%(plottype,str(expset))
        print '='*40
        
    
    def __get_color(self,cmap,value=0.8):
        cNorm = colors.Normalize(vmin=0,vmax=1) 
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)
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
        
    def plot_FI(self,experimentsets,savefig=None,settings='paper'):
        rcsettings = RC_SETTINGS[settings]
        pylab.figure() 
        intensities = np.linspace(0.1, 0.9, len(experimentsets))

        for (i,expset) in enumerate(experimentsets):
            try:
                print 'i=%g, intensities[i]='%i,intensities[i]
                print 'exp variables len', len(expset.expvariables[0])
                print expset.expvariables[0]
                print 'length FI', len(expset.results['FI'])
                pylab.plot(expset.expvariables[0],expset.results['FI'],lw=rcsettings['lw'],c=self.__get_color(rcsettings['cmaps'][0],intensities[i]),label=expset.get_label())
            except:
                print 'Error with plotting for expset', expset
                continue
        # turn the right and top axes off
        self.__turn_off_border(pylab.gca())
        
        pylab.ylabel('Output frequency (Hz)')
        pylab.xlabel('Input')
        # increase scale so that we have some space to place the legend
        pylab.ylim(ymax=pylab.ylim()[1]*1.2)
        
        if savefig is not None:
            figname = '%sFI_%s.png'%(savefig,time.strftime('%y%m%d'))
            pylab.savefig(figname)
            print 'Saved figure as %s'%figname        

        # Prettify the legend
        leg = pylab.legend(loc=2,fancybox=True)
        leg.get_frame().set_alpha(0.5) 
        
        if savefig is not None:
            figname = '%sFI_%s_legend.png'%(savefig,time.strftime('%y%m%d'))
            pylab.savefig(figname)
            print 'Saved figure as %s'%figname
            
            
    def __poly(self,x,a,b,c,d):
        return a*x**3 + b*x**2 + c*x + d
                
    def plot_fit_FI(self,experimentsets,savefig=None,settings='poster'):
        
        rcsettings = RC_SETTINGS[settings]
        pylab.rcParams.update(rcsettings)
        pylab.figure() 
        intensities = np.linspace(0.1, 0.9, len(experimentsets))
        from matplotlib import mpl
        norm = mpl.colors.BoundaryNorm(intensities, len(intensities))
        
        
        
        #cNorm = colors.Normalize(vmin=0,vmax=1) 
        #scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)
        
        m = cm.ScalarMappable(cmap=rcsettings['cmaps'][0])
        m.set_array(intensities)

        for (i,expset) in enumerate(experimentsets):
            #try:
            if True:
                print 'i=%g, intensities[i]='%i,intensities[i]
                print 'exp variables len', len(expset.expvariables[0])
                print expset.expvariables[0]
                print 'length FI', len(expset.results['FI'])
                xvals = np.linspace(expset.expvariables[0][0], expset.expvariables[0][-1],num=len(expset.expvariables[0]))
                svals = np.asarray(xvals).ravel() 
                fis = np.asarray(expset.results['FI']).ravel()
                
                popt_curve, pcov_curve = curve_fit(self.__poly, svals, fis,p0=[1.,1.,1.,0.])
                fitted = self.__poly(svals,popt_curve[0],popt_curve[1],popt_curve[2],popt_curve[3])
                #scatter real results, plot trend
                pylab.scatter(expset.expvariables[0],expset.results['FI'],lw=rcsettings['lw']/2,c=self.__get_color(rcsettings['cmaps'][0],intensities[i]),s=rcsettings['lw']*10,facecolor='none',edgecolor=self.__get_color(rcsettings['cmaps'][0],intensities[i]),alpha=0.6)
                pylab.plot(svals,fitted,lw=rcsettings['lw'],c=self.__get_color(rcsettings['cmaps'][0],intensities[i]),label=expset.get_label())
            """except:
                print 'Error with plotting for expset', expset
                print "Unexpected error:", sys.exc_info()[0]
                continue
            """
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
            figname = '%sFIfit_%s.png'%(savefig,time.strftime('%y%m%d'))
            pylab.savefig(figname,dpi=rcsettings['dpi'])
            print 'Saved figure as %s'%figname        

        # Prettify the legend
        leg = pylab.legend(loc=2,fancybox=True)
        leg.get_frame().set_alpha(0.5) 
        
        if savefig is not None:
            figname = '%sFIfit_%s_legend.png'%(savefig,time.strftime('%y%m%d'))
            pylab.savefig(figname,dpi=rcsettings['dpi'])
            print 'Saved figure as %s'%figname        
        
        
        
        
        
        
        
        
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

        
             
        

    def plot_FIphoto(self,experimentsets,savefig=None,settings='paper'):
        #TODO: plot
        pass
        """
        rcsettings = RC_SETTINGS[settings]
        pylab.figure() 
        intensities = np.linspace(0.1, 0.9, len(experimentsets))
        for (i,expset) in enumerate(experimentsets):
            print 'i=%g, intensities[i]='%i,intensities[i]
            pylab.plot(expset.expvariables,expset.results['FI'],lw=rcsettings['lw'],c=self.get_color(rcsettings['cmaps'][0],intensities[i]),label=expset.get_label())
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
            
            pylab.figure() 
            intensities = np.linspace(0.1, 0.9, len(expset.experiments))
            for (j,exp) in enumerate(expset.experiments):
                (ts,vsoma) = exp.get_voltage_trace()
                print 'j=%g, intensities[i]='%j,intensities[j],' .................',
                print exp
                pylab.plot(ts,vsoma,lw=rcsettings['lw_fine'],c=self.__get_color(rcsettings['cmaps'][0],intensities[j]),alpha=rcsettings['alpha'])
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
            
            pylab.figure() 
            intensities = np.linspace(0.1, 0.9, len(expset.experiments))
            for (j,exp) in enumerate(expset.experiments):
                (ts,vsoma) = exp.get_voltage_trace()
                print 'j=%g, intensities[i]='%j,intensities[j],' .................',
                print exp
                pylab.plot(ts,vsoma,lw=rcsettings['lw_fine'],c=self.__get_color(rcsettings['cmaps'][0],intensities[j]),alpha=rcsettings['alpha'])
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
            
                fig = pylab.figure() 
                intensities = np.linspace(0.1, 0.9, len(expset.experiments))

                try:
                    for (j,exp) in enumerate(expset.experiments):
                        
                        
                        
                        (ts,iphoto) = exp.get_iPhoto_opsin(opsintype,entire_trace=False)
                        print ts
                        print 'j=%g, intensities[i]='%j,intensities[j],' .................',
                        print exp
                        ax = fig.gca()
                        ax.plot(ts,iphoto,lw=rcsettings['lw_fine'],c=self.__get_color(rcsettings['cmaps'][0],intensities[j]),alpha=rcsettings['alpha'])
                except:
                    pylab.close(fig)
                    continue
                self.__turn_off_border(ax)
                ax.set_xlabel('time (ms)')
                ax.set_ylim(auto=True)
                ax.set_ylabel('I_%s (nA)'%opsintype)
                
                pylab.gcf()
                inset_ax = pylab.axes([.2, .7, .15, .15])
                self._plot_injected_current(ies, rcsettings['cmaps'][0], inset_ax)
                
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
                fig = pylab.figure() 
                ax = fig.gca()
                ax.plot(ies,iphotos,lw=rcsettings['lw'],c=self.__get_color(rcsettings['cmaps'][0],intensities[j]))
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
            pylab.clf()
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
        self.results = {} #TODO : get rid of this
        self.experiments = []
    
    def __str__(self):
        return self.basename+self.subselect

    def get_label(self):
        return self.label
    
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

        for v in list(itertools.product(*self.expvariables)):
            #print v
            #print self.basename , self.subselect%v
            self.experiments.append(Experiment(self.basename ,self.subselect%v,analysis_params,vals=v))
    
    def save_experiments(self):
        print 'Saving ExperimentSet', self.label
        for exp in self.experiments:
            exp.save_results()
    
    def load_experiments(self):
        print 'Loading ExperimentSet', self.label
        for exp in self.experiments:
            exp.load_results()        
            #once data is loaded, go through and repopulate self.
            for (k,v) in exp.results.iteritems():
                print k, v
                try:
                    self.results[k].append(v)
                except: #doesn't have key
                    self.results[k] = [v]
        print self.results
    #def set_experiments(self,key,analysis_params):
        #for v in self.expvariables[key]:
        #    print v
        #    print self.basename , self.subselect%v
        #    self.experiments.append(Experiment(self.basename ,self.subselect%v,analysis_params))
    #    pass
        
    def run_analysis(self,analysis_fns,recalc=False):
        print 'in expset = ', analysis_fns
        for exp in self.experiments:
            print 'About to run analysis for ',str(exp)
            exp.calculate_results(analysis_fns, recalc)
            exp.save_results()
    
    def calculate_responses(self,key,recalc=False):
        print 'calculate responses', key
        # TODO: rewrite with debugger/warning info statements
        #if self.expvariables.has_key(key):
        try:
        #if True:
            if not recalc:
                if self.results.has_key(key):
                    print 'Results for %s already exists'%key 
                    return
            tmp_results = []
            for exp in self.experiments:
                #print 'looking up ',FUNCTION_LOOKUP[key]
                tmp_v = getattr(exp, FUNCTION_LOOKUP[key])()
                #print tmp_v, '-------------------------------------------------------'
                tmp_results.append(tmp_v)
            self.results[key] = tmp_results
        except:
            print "Could not calculate responses, unknown variable value = ", key
        #print 'Current results ',self.results
    
    def print_exps(self,indent=4):
        print self.basename,self.subselect, 'contains experiments:'
        for exp in self.experiments:
            print ' '*indent, exp, '\t'*4, exp.results
            
    
    def print_data(self,indent=4):
        print self.basename,self.subselect, 'contains data:'
        print 'label = ',self.label
        for val in self.expvariables:
            print ' '*indent, '%s'%(val)
        for (k,v) in self.results.iteritems():
            print ' '*indent, '%s = %s'%(k,v)
            


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
        

    def calc_firingrate(self,**params):
        """
            expname     experiment name
            t_cutout    list of start time and finish times during which firing rate should be calculated
        """
        #load dat curve, get v_soma
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
                self.results['FI'] = 0
                return 0
            spikes = np.where((spikes >= self.analysis_params['tstart']) & (spikes <= self.analysis_params['tstop']))[0]
            print 'calc_firing rate'
            spikerate = 1000.*len(spikes)/(self.analysis_params['tstop']-self.analysis_params['tstart'])
            self.results['FI'] = spikerate
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


    def get_spiketrain(self):
        if self.spiketrain is not None:
            pass
        elif self.spiketrain is None and self.spikes is not None:
            self.spiketrain = NeuroTools.signals.SpikeTrain(self.spikes,self.analysis_params['tstart'],self.analysis_params['tstop'])
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
        
        
        try:
            print 'Am going to check if spike file exists ...'
            if self._check_spikefile_exists():
                print 'File already exists'
                return
            data = np.loadtxt(self.get_dat_file()) 
            t,v_soma = data[0,:],data[1,:]
            spikes = self.extract_spiketimes(t,v_soma)
            # we want to have all spikes, not just the ones in our area of analysis
            self.save_spikes(spikes)
            return spikes
        except:
            pass #print "Unexpected error:", sys.exc_info()[0]
    
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
        

