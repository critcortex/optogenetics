# -*- coding: utf-8 -*-


import file_io as fio

from neuron import h
from nrn    import *
import pylab
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

opsin_dict = {}
opsin_dict['ChR'] = {'ycoord':-5 , 'color':'b'}
opsin_dict['NpHR'] = {'ycoord':-10 , 'color':'#FFA500'} # orange

class NeuronExperiment:
    
    #TODO update so that any cell morphology
    
    def __init__(self):
        h.load_file('stdlib.hoc', 'String') 
        h.load_file('stdrun.hoc')
        h.load_file("import3d.hoc")
        self.stimuli = {}
    
    def get_default_params(self):
        params = { 'tstart': 300,
                   'tstop' : 1500,
                   'proxpt': 400,
                   'distpt': 620,
                   'vinit' :-80,
                   'squareamp': 1.,
                   'EPSP_transient':5,
                   'imax': 0.5}
        return params
    
    def setup(self,params): #tstart=300,tstop=1500,proxpt=400,distpt=620,vinit=-80,squareamp=1.,imax=0.5):
        """
        Sets up and initializes parameters and simulation environment, including cell morphology   
        """
        global squareAmp, proximalpoint, distalpoint
        global risetau, decaytau, BACdt
        
        #====================== General files and tools =====================
        #h.load_file("nrngui.hoc")
    
        #====================== cvode =======================================
        h('objref cvode')
    
        h.cvode = h.CVode()
        h.cvode.active(1)
    
        #=================== creating cell object ===========================
        
        h('objref L5PC')
    
        h('strdef morphology_file')
        h.morphology_file = "morphologies/cell1.asc"
    
        h.load_file("models/L5PCbiophys3.hoc")
        h('load_file("models/L5PCtemplate.hoc")')
        h.L5PC = h.L5PCtemplate(h.morphology_file)
    
        #=================== settings ================================
        #v_init = vinit   # orig val: -80
        
        
            
        #somatic pulse settings
        squareAmp = params['squareamp']  # original value: 1.9
    
        #EPSP settings
        risetau = 0.5
        decaytau = 5
        BACdt = params['EPSP_transient']
        
    
        proximalpoint = params['proxpt']
        distalpoint = params['distpt']
        self.set_proxdist_locations(proximalpoint,distalpoint)
    
        h('objref tstart')
        h.tstart = params['tstart']
        h.tstop = params['tstop']
    
    
    def set_proxdist_locations(self,proximalpoint,distalpoint):
        h('objref sl')
        h('double distSite[2]')
        h('double proxSite[2]')
        h.sl = h.List()
        maxdiam = 0
        
        h.sl = h.L5PC.locateSites("apic",distalpoint)
        axsites = int(h.sl.count())
        j = 0
        for i in range(axsites):
            dd1 = h.sl.object(i).x[1]
            dd = h.L5PC.apic[int(h.sl.object(i).x[0])](dd1).diam
            if (dd > maxdiam) :
                j = i
                maxdiam = dd 
        h.distSite[0] = h.sl.object(j).x[0]
        h.distSite[1] = h.sl.object(j).x[1]
        
        
        h.sl = h.L5PC.locateSites("apic",proximalpoint)
        axsites = int(h.sl.count())
        j = 0
        
        for i in range(axsites):
            dd1 = h.sl.object(i).x[1]
            dd = h.L5PC.apic[int(h.sl.object(i).x[0])](dd1).diam
            if (dd > maxdiam):
                j = i
                maxdiam = dd 
            
        h.proxSite[0] = h.sl.object(j).x[0]
        h.proxSite[1] = h.sl.object(j).x[1]
    
    
    
    def set_experiment_type (self,params):
        """
        Define the type of experiment type to be performed.
        Possible values are:
            BAP (default)      current injected at soma
            CaBurst            EPSP-like event injected at apical location
            BAC                current injected at soma and EPSP events
            opsinonly          no current or EPSP injection; intended to examine the solo contribution of opsins
        
        """
        # TODO: extend this to allow for concurrent stimulation types easily - possibly define as list or dict of stimuli  
        experiment_type = params['experiment_type']
        global somastimamp, EPSPamp
        somastimamp = params['iclamp_amp']
        EPSPamp = params['EPSP_amp']

        print "Setting up for experiment type", experiment_type
        
        if (experiment_type=="BAP"):
            #somastimamp = squareAmp
            EPSPamp = 0
        
        elif (experiment_type=="CaBurst"):
            somastimamp = 0
            #EPSPamp = Imax*3
        
        elif (experiment_type=="BAC"):
            pass
            #somastimamp = squareAmp
            #EPSPamp = Imax
            
        elif (experiment_type=="opsinonly"):
            somastimamp = 0
            EPSPamp = 0
        
        else:
            print "WARNING: no experiment type was set"
            print "Valid options: [BAP | CaBurst | BAC | opsinonly]"
            print "Assuming = BAP"
            
            #somastimamp = squareAmp
            EPSPamp = 0
        
        
    
    
    def add_optogenetics (self,opsindict):
        """
        Loads relevant hoc files for opsins in revelant sections of neuron. 
        Params: 
              optoexpress     = Dictionary of opsins, where keys are opsin names, 
                                while the value is a list of the locations the opsin is to 
                                be expressed in 
        Example: 
        
        Equivalent to [for NpHR in apical distal:
            h.load_file("NpHR-in-apical_distal.hoc")  
            print "loaded NpHR apical successful"
            h('HR_in_apical_distal(5e-4,10,600,200,350,600)')
        """
        for (opsin, loclist) in opsindict.iteritems():
    
            loclist.sort(key=lambda x:x[0],reverse=True)    
            for loc in loclist: #we reverse the sort so that soma is first
                if loc[0] is not None: # we allow users to specify None in order to make looping easier
                    print "Locations specified for opsin %s: "%opsin, loc[0]
                    filename = '%s-in-%s.hoc'%(opsin,loc[0])
                    try:
                        h.load_file(filename)
                        print "loaded %s %s successful"%(opsin,loc[0])
                    except:
                        print "Error: could not load hoc file for opsin %s in location: "%opsin, filename
                    ss = '%s_in_%s(%s)'%(opsin,loc[0],','.join([str(x) for x in loc[1]]))
                    print ss
                    h(ss)
                else:
                    print "No locations specified for opsin : ",opsin
        
    
    def set_stimulus (self,params) :
        """
        Set stimuli objects:
            st1        IClamp object at soma
            st2        IClamp object at distal location on apical dendrite (distSite)
            syn1       EPSP event at distal location on apical dendrite (distSite)
        """
        #======================== stimulus settings ============================
        global vec
        vec = {}
        
        """
        # Determine location of distSite
        h('objref sl')
        h.sl = h.List()
        h('double distSite[2]')
        h.sl = h.L5PC.locateSites("apic",distalpoint)
        axsites = int(h.sl.count())
        maxdiam = 0
        j = 0
        for i in range(axsites):
            dd1 = h.sl.object(i).x[1]
            dd = h.L5PC.apic[int(h.sl.object(i).x[0])](dd1).diam
            if (dd > maxdiam) :
                j = i
                maxdiam = dd 
        h.distSite[0] = h.sl.object(j).x[0]
        h.distSite[1] = h.sl.object(j).x[1]
        """
        # ------------------------ Current pulses
        h('objref st1, st2')
        h.st1 = h.IClamp()
        h.st1.loc(0.5, sec=h.L5PC.soma[0])
        h.st1.amp = somastimamp
        #st1.del = tstart --> error as del is reserved word in python
        st1_delay = params['iclamp_start']
        st1_duration = params['iclamp_duration']
        setattr(h.st1, 'del', st1_delay)
        h.st1.dur = st1_duration
        h('L5PC.soma st1')
        
        
        h('access L5PC.apic[distSite[0]]')
        h.st2 = h.IClamp(h.distSite[1])
        setattr(h.st2, 'del', params['iclamp_dist_start'])
        h.st2.dur = params['iclamp_dist_duration']
        h.st2.amp = params['iclamp_dist_amp']
        h('L5PC.apic[distSite[0]] { st2 } ') ########## TODO: extend and make this an active possibility to use
        
        
        # ------------------------ Dendritic EPSP-like current
        # added to segments on apical dendrite that are distal to <distalpoint>
        h('objref syn1,isyn, tvec')
        
        h.isyn = h.Vector()
        h.tvec = h.Vector()
        
        h.syn1 = h.epsp(h.distSite[1])
        h.syn1.tau0 = risetau       
        h.syn1.tau1 = decaytau   
        h.syn1.onset = h.tstart + BACdt  
        h.syn1.imax = EPSPamp
        h('L5PC.apic[distSite[0]] { syn1 } ')
        
        h.cvode.record(h.syn1._ref_i,h.isyn,h.tvec)
        
        #    
    #    # Locate the electrode at the center of the soma
    #    h('objref stim, hoc_vector, mytvec')
    #    tmppy_tvec = range(int(h.tstop))
    #    #h.mytvec = h.Vector(2)
    #    #h.mytvec.x[0] = 0
    #    #h.mytvec.x[1] = h.tstop
    #    h.mytvec = h.Vector(tmppy_tvec)
    #    VecT = h.Vector([0, 1000])
    #    
    #    stim = h.IClamp()
    #    stim.loc(0.5, sec=h.L5PC.soma[0])
    #    stim.dur = 1e9
    #    
    #    # Setting recording paradigm
    #    #h.stim.delay = 100
    #    #h.stim.amp = 100
    #    #h.stim.dur = 400
    #    
    #    import numpy.random as nprnd
    #    stimvector = nprnd.randint(1000, size=2000)
    #    print stimvector[:10]
    #    print tmppy_tvec[:10]
    #    hoc_vector = h.Vector(stimvector)
    #    #h('hoc_vector.play(&stim.i,mytvec)')
    #    hoc_vector.play(stim._ref_amp, h.tvec, 1)
        print "-----------------------------------------------------------------------==============================================="
        
    def get_location(self,location):
        #TODO: update so it returns proper section
        return h.L5PC.soma[0]
    
        
        
    def __create_iclamp(self,params, location):
        h.st1 = h.IClamp()
        h.st1.loc(0.5, sec=self.get_location(location))
        h.st1.amp = params['amp']
        h.st1.dur = params['duration']
        #st1.del = tstart --> error as del is reserved word in python
        st1_delay = params['start']
        setattr(h.st1, 'del', st1_delay)
        
        
        
    def __create_epsp(self,params, location):
        h.syn1 = h.epsp(self.get_location(location))
        h.syn1.tau0 = params['risetau']       
        h.syn1.tau1 = params['decaytau']
        h.syn1.onset = h.tstart + params['BACdt']
        h.syn1.imax = params['EPSPamp']
        #h('L5PC.apic[distSite[0]] { syn1 } ')
    
    def __create_spiketrain(self,params,location):
        name = 'test'
        """
        self.syn[name] = h.ExpSyn(self.get_location(location))
        self.stim[name] = h.Vector(np.sort(params['tstim'])) # Converting tstim into a NEURON vector (to play in NEURON)
        self.vplay[name] = h.VecStim() # Creating play vectors to interface with NEURON
        self.vplay[name].play(self.stim[name])  # Connecting vector to VecStim object to play them
        self.netcon[name] = h.NetCon(self.vplay[name], self.syn[name]) # Building the netcon object to connect the stims and the synapses
        self.netcon[name].weight[0] = params['w'] # Setting the individual weights
        """
        
        self.syn[name] = h.ExpSyn(self.get_location(location))
        self.stim[name] = h.Vector(np.sort(params['tstim'])) # Converting tstim into a NEURON vector (to play in NEURON)
        self.vplay[name] = h.VecStim() # Creating play vectors to interface with NEURON
        self.vplay[name].play(self.stim[name])  # Connecting vector to VecStim object to play them
        self.netcon[name] = h.NetCon(self.vplay[name], self.syn[name]) # Building the netcon object to connect the stims and the synapses
        self.netcon[name].weight[0] = params['w'] # Setting the individual weights
    
    # TODO: implement and replace for old set_stimulus
    def new_set_stimulus(self,params):
        
        for (stimtype,stimparams) in params.iteritems():
            for location in stimparams['locations']:
                # work out what type of stimulus it is, create corresponding object
                if stimtype == 'IClamp':
                    stim = self.__create_iclamp(stimparams['stimparams'],location)
                elif stimtype == 'EPSP': 
                    stim = self.__create_epsp(stimparams['stimparams'],location)
                elif stimtype == 'spiketrain':
                    stim = self.__create_spiketrain(stimparams['stimparams'],location)
                #TODO: update it so that stim is returned to an appropriate location
    
    
    
    def setup_record(self,opsindict={}):
        """
        Set up recording devices and structures.
        Recording:
            vsoma    voltage from soma
            vdend    voltage from distal location on apical dendrites
            
        
        """
        
        #======================== recording settings ============================
        h('objref vsoma, vdend, proxClamp, vdend2, isoma, gsoma, isoma_k, isoma_na')
    
        h.vsoma = h.Vector()
        h.isoma = h.Vector()
        h('access L5PC.soma')
        h('cvode.record(&v(0.5),vsoma,tvec)')
        h.cvode.record(h.st1._ref_i,h.isoma,h.tvec)
        h.isoma_k = h.Vector()
        h.isoma_na = h.Vector()
        h('cvode.record(&ik(.5),isoma_k,tvec)')
        h('cvode.record(&ina(.5),isoma_na,tvec)')
        
        # TODO: would also like to record conductances from soma, etc.
        #h.gsoma = h.Vector()
        #h.cvode.record(h.st1._ref_g,h.gsoma,h.tvec)
        
        # Record voltage at distal point on apical tree
        h.vdend = h.Vector()
        h('access L5PC.apic[distSite[0]]')
        h.cvode.record(h.L5PC.apic[int(h.distSite[0])](h.distSite[1])._ref_v,h.vdend,h.tvec)
        
    
        # Record voltage at proximal point on apical tree
        h('access L5PC.apic[proxSite[0]]')
        h.vdend2 = h.Vector()
        h.cvode.record(h.L5PC.apic[int(h.proxSite[0])](h.proxSite[1])._ref_v,h.vdend2,h.tvec)
    
        # Following lines are only useful when we're wanting to plot location of proxSite
    #    h('access L5PC.apic[proxSite[0]]')
    #    h.proxClamp = h.IClamp(h.proxSite[1])
    #    h.proxClamp.amp = 0
    #    h('L5PC.apic[proxSite[0]] { proxClamp }')
    
        # for each opsin, include vector for corresponding photocurrent
        h('access L5PC.soma')
        for (opsin,opsinval) in opsindict.iteritems():
            print 'testing for %s'%opsin, opsinval
            if self.is_opsin_present(opsinval,opsin):
                h('objref i_%s'%opsin)
                h('i_%s = new Vector()'%opsin)
                h('cvode.record(&%s[0].i%s,i_%s,tvec)'%(opsin,opsin,opsin)) 
                print 'Recording from i_%s'%opsin
    
    
    def setup_plot (self):
        
        ## ARTEFACT OF NEURON CODE. Note that we're now plotting using python matplot methods instead
        pass
        #h('objref gV, gI, s')
    
        #h.s = h.Shape(shape=h.L5PC.all)
        #h.s.color_list(h.L5PC.axonal,2,2)
        #h.s.color_list(h.L5PC.somatic,5)
        #h.s.color_list(h.L5PC.basal,4)
        #h.s.color_list(h.L5PC.apical,1)
        #h.s.point_mark(h.st2,2)
        ##h.s.point_mark(h.proxClamp,3)
        #h.s.show(0)
        
    
    
    
    def simulate_exp(self):
        
        print "Init experiment ...",
        h.init()
        print "Init successful. About to run experiment ...",
        h.run()
        print "...finished"
        
    def is_opsin_present(self,opsinvalue, opsin):
        """
        opsinvalue = opsindict[opsin]
        """
        print opsin, opsinvalue[0][0],'======================================================= is Opsin present? ', 
        bb = opsinvalue[0][0] is not None
        print bb
        return opsinvalue[0][0] is not None
        
    def save_data(self,expname,savedata=False,opsindict={}):
        if savedata:
            import numpy as np
            mat = np.matrix([h.tvec,h.vsoma,h.vdend,h.vdend2,h.isoma,h.isoma_k,h.isoma_na])
            np.savetxt(expname+".dat",mat)
            
            if len(opsindict.keys())>0:
                
                # A horrible hack to get around python/NEURON conversion fun
                """
                h('objref list_i_opsin')
                for opsin in opsindict.keys():
                    if not self.is_opsin_present(opsindict[opsin],opsin):
                        continue
                    h('list_i_opsin = new List()')
                    # add it to opsin list
                    h('list_i_opsin.append(i_%s)'%opsin)
                    #previously: we were only saving the opsin current
                    # mat = np.matrix(h.list_i_opsin)
                    #but it's better to save both time and value (as we're offgrid) therefore:
                    mat = np.matrix([h.tvec,np.array(h.list_i_opsin.object(0))])
                    # save data to dat file
                    np.savetxt(expname+"_i%s.dat"%opsin,mat)
                """
                for opsin in opsindict.keys():
                    if not self.is_opsin_present(opsindict[opsin],opsin):
                        continue
                    h('objref savdata')
                    h('savdata = new File()')
                    h('savdata.wopen("%s_i%s.dat")'%(expname,opsin))
    
                    #h('savdata.printf("t SThcells[2].soma.v(0.5)\n")')
    
                    h('objref tempmatrix')
                    h('tempmatrix = new Matrix()')
                    h('tempmatrix.resize(tvec.size(),2)')
                    h('tempmatrix.setcol(0, tvec)')
                    h('tempmatrix.setcol(1, i_%s)'%opsin)
                    h('tempmatrix.fprint(savdata, " %g")')
                    h('savdata.close()')
    
    def plot_optogenetic(self,opsin_type,ton,toff,yoffset):
        """ 
        Plots lines at top of figure to indicate when each opsin was active 
        Y-value and color set in static dict 'opsin_dict' 
        """
        yval = opsin_dict[opsin_type]['ycoord']
        color = opsin_dict[opsin_type]['color']
        ax = pylab.gca()
        x = [ton,toff]
        y = np.array([yoffset+yval,yoffset+yval])
        ax.fill_between(x,y-2,y+2,color=color)
    
    def run_plots(self,params):
        """
        Generates plots for voltage and current at soma, distal and proximal locations on apical dendrites
        """
        lw = 2
        
        # Plot voltage at soma and dendrites (apical proximal and distal)
        pylab.figure(1)
        pylab.plot(h.tvec,h.vsoma,lw=lw,c='k',label='v_soma')
        pylab.plot(h.tvec,h.vdend,lw=lw,c='r',label='v_dend')
        pylab.plot(h.tvec,h.vdend2,lw=lw,c='b',label='v_dend2')
        pylab.xlim(h.tstart-20,h.tstop+20)
        pylab.ylim(-120,40)
        # If optogenetics were included, draw blocks for times that illumination occurred in appropriate colours 
        if params.has_key('opdict'):
            for (opsin,opexpressions) in params['opdict'].iteritems():
                for opexp in opexpressions:
                    if opexp[0] is None or opexp[0].lower() == 'none':
                        continue
                    for pulsenum in range(opexp[1][6]):                
                        pulse_start = opexp[1][2]+pulsenum*(opexp[1][3]+opexp[1][4])
                        self.plot_optogenetic(opsin,pulse_start,pulse_start+opexp[1][3],yoffset=40)
                    # once we've plotted an activation for one area, that should be sufficient i.e. we don't need to plot apical *and* soma, only the first 
                    # TODO: think how to extend this to allow for different areas to be indicated i.e. ChR in soma vs ChR in apical dendritic arbor
                    break
        pylab.title('V')
        ax = pylab.gca()
        for loc, spine in ax.spines.iteritems():
            if loc in ['left','bottom']:
                spine.set_position(('outward',5))
                ax.tick_params(direction='out')
            elif loc in ['right','top']:
                spine.set_color('none') 
        pylab.legend()
        pylab.xlabel('time (ms)')
        pylab.ylabel('V (mV)')
        
        # Plot currents at soma and i_syn
        pylab.figure(2)
        pylab.plot(h.tvec,h.isyn,lw=lw,c='g',label='i_syn')
        pylab.plot(h.tvec,h.isoma,lw=lw,c='k',label='i_soma')
        if params.has_key('opdict'):
            for (opsin,opexpressions) in params['opdict'].iteritems():
                for opexp in opexpressions:
                    if opexp[0] is None or opexp[0].lower() == 'none':
                        continue
                    h('objref list_i_opsin')
                    h('list_i_opsin = new List()')
                    h('list_i_opsin.append(i_%s)'%opsin)
                    pylab.plot(h.tvec,h.list_i_opsin.object(0),color=opsin_dict[opsin]['color'],label='i_%s'%opsin)
                    break
        pylab.xlim(h.tstart-20,h.tstop+20)
        #pylab.ylim(-3,6)
        pylab.title('I')
        ax = pylab.gca()
        for loc, spine in ax.spines.iteritems():
            if loc in ['left','bottom']:
                spine.set_position(('outward',5))
                ax.tick_params(direction='out')
            elif loc in ['right','top']:
                spine.set_color('none') 
        pylab.legend()
        pylab.xlabel('time (ms)')
        pylab.ylabel('I (nA)')
        
        
        if params['expname'] is not None:
            savename = params['expname']
            pylab.figure(1)
            pylab.savefig(savename+'_voltage.png')
            pylab.figure(2)
            pylab.savefig(savename+'_current.png')
            print "Saved figures under %s*.png"%savename
            pylab.close('all')
        else:
            pylab.show()
            
            
    
    def debug(self):
        print "tvec range = ",h.tvec.min(),h.tvec.max()
        print "isyn range = ",h.isyn.min(),h.isyn.max()
        print "vdend range = ",h.vdend.min(),h.vdend.max()
        print "vdend2 range = ",h.vdend2.min(),h.vdend2.max()
        print "vsoma range =", h.vsoma.min(), h.vsoma.max()
        print "isoma range = ",h.isoma.min(),h.isoma.max()
       
     
    
    def main(self,exp_keydict):
        print '\n'*2, '='*40, exp_keydict['expname']
        keydict = self.get_default_params()
        keydict.update(exp_keydict)
        self.setup(keydict)#tstart=keydict['tstart'],tstop=keydict['tstop'])
        self.set_experiment_type(keydict)
        self.set_stimulus(keydict)
        self.add_optogenetics(keydict['opdict'])
        self.setup_record(keydict['opdict'])
        self.setup_plot()
        self.simulate_exp()
        self.save_data(keydict['expname'],keydict['savedata'],keydict['opdict'])
        self.run_plots(keydict)
        self.debug()
    
    
    def run_from_file(self,paramfilename):
        
        job = fio.loadjob(paramfilename)
        self.main(job)
        return job['expname_family']
    
if __name__ == '__main__':
    NE = NeuronExperiment()
    if len(sys.argv)>1:
        print sys.argv[1:]
        if sys.argv[1]=='run_from_file': # TODO: make this more generic
            for fn in sys.argv[2:]:
                print 'run from file-----------', fn 
                if True: #try:
                    basename = NE.run_from_file(fn)
                #try:
                #    fio.finish_experiment(basename)
                #except:
                #    print 'Error with running job',fn
    
       
#if __name__ == '__main__':
#    default_dict = {}
#    default_dict['experiment_type'] = None
#    default_dict['opdict'] = {}
#    default_dict['expname'] = "defaulttest"
#    NE = NeuronExperiment()
#    NE.main(default_dict)
