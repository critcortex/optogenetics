# -*- coding: utf-8 -*-


import file_io as fio

from neuron import h, run, init
#from nrn    import *

import Neuron  
import opsin as oplib

import pylab
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import random

opsin_dict = {}
opsin_dict['ChR'] = {'ycoord':-5 , 'color':'b'}
opsin_dict['NpHR'] = {'ycoord':-10 , 'color':'#FFA500'} # orange


# 
# class L5PC:
#     
#     def __init__(self):
#         h.load_file('stdlib.hoc') 
#         h.load_file('stdrun.hoc')
#         h.load_file("import3d.hoc")
#         h.load_file("models/L5PCbiophys3.hoc")
#         h.load_file("models/L5PCtemplate.hoc")
#         
#         self.cell = h.L5PCtemplate("morphologies/cell1.asc")
#         
#     def get_cell(self):
#         return self.cell
#     
#     def get_soma(self):
#         return self.cell.soma
#     
#     def get_basal(self):
#         return self.cell.basal
#     
#     def get_apical(self):
#         return self.cell.apic
# 
#     def get_axon(self):
#         return self.cell.axon
#     
#     
#     
    
    

class NeuronExperiment:
    
    
    def __init__(self):
        h.load_file('stdlib.hoc', 'String') 
        h.load_file('stdrun.hoc')
        h.load_file("import3d.hoc")
        self.stimuli = {}
        self.outputparams = self.get_default_outputparams()
        self.params = self.get_default_params()
        
    
    def get_default_params(self):
        params = { 'tstart': 0,
                   'tstop' : 1500,
                   'proxpt': 400,
                   'distpt': 620,
                   'vinit' :-80,
                   'squareamp': 1.,
                   'EPSP_transient':5,
                   'imax': 0.5}
        return params
    
    def _use_setup_params(self):
        params = {'cell': ['Neuron','L5PC'],
                  'cell_params':{},
                  'tstart': 0,
                  'tstop': 100,
                  }
        self.params = params
    
    def get_default_outputparams(self):
        outputparams = {}
        return outputparams
    
    def setup(self): #tstart=300,tstop=1500,proxpt=400,distpt=620,vinit=-80,squareamp=1.,imax=0.5):

        """ 
        Cell type
        """
        #print "Got to here"
        __import__(self.params['cell'][0]) 
        neuronclass = sys.modules[self.params['cell'][0]]
        CellClass = getattr(neuronclass,self.params['cell'][1])
        self.cell = CellClass(self.params['cell_params'])

        """
        if self.params['cell'] is None:
            self.cell = Neuron.L5PC()
        else:
            self.cell = self.params['cell']
        """
        
        #====================== General files and tools =====================
        #h.load_file("nrngui.hoc")
        """
    
        #====================== cvode =======================================
        h('objref cvode')
    
        h.cvode = h.CVode()
        h.cvode.active(1)
        """
        #=================== creating cell object ===========================
        """
        h('objref L5PC')
    
        h('strdef morphology_file')
        h.morphology_file = "morphologies/cell1.asc"
    
        h.load_file("models/L5PCbiophys3.hoc")
        h('load_file("models/L5PCtemplate.hoc")')
        h.L5PC = h.L5PCtemplate(h.morphology_file)
        """
        
        #proximalpoint = self.params['proxpt']
        #distalpoint = self.params['distpt']
        ######self.set_proxdist_locations(proximalpoint,distalpoint)
    
        h('objref tstart')
        h.tstart = self.params['tstart']
        #h('objref tstop')
        h.tstop = self.params['tstop']
        #h.dt = 0.01
        h.t = 0
        
        self.tag_locations()
        
        
    def set_experiment_type(self):
        pass
    
    def add_optogenetics(self):
        """
        print '....... ADDING TO SOMA ------------------'
        opsintype = 'ChR'
        opsinProtein = getattr(h,opsintype)()
        try:
            opsinProtein.loc(0.5, sec=self.cell.get_soma()[0])
        except:
            opsinProtein.loc(0.5, sec=self.cell.get_soma())
        
        print self.cell.get_num_segments('apical')
        print self.cell.get_num_segments('dend0')
        print '....... ADDED TO SOMA ------------------', '\n'*20
        """
        
        
        opsinExp = oplib.Opsin() 
        oplist = opsinExp.express_opsin(self.cell, self.params)
        self.outputparams['expressed_opsin'] = oplist
        
        
    def add_optogenetics_old (self,opsindict):
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
   


    # TODO: implement and replace for old set_stimulus
    def new_set_stimulus(self):
        """
        
        """
        #print self.params
        stimtypes = [key for key in self.params.keys() if key.startswith('stim_') and self.params[key]]
        print '\nThe following stimulus types were seen:', stimtypes, '\n'
        for k in stimtypes:
            #print self.params[k.replace('stim_','')]
            self.outputparams[k] = {}
            for stimulus in self.params[k.replace('stim_','')]:
                print 'going to run: ',"getattr(self,'_create_%s')(stimulus)"%k.replace('stim_','') #,' where stimulus = ', stimulus
                #print stimulus,'--------------------------------'
                getattr(self,'_create_%s'%k.replace('stim_',''))(stimulus)
        #print self.outputparams
        
    def get_ids_distance(self,section,distances=(100,200),cell=None):
        if cell is None:
            try:
                cell = self.cell
            except:
                print "Function requires cell"
                return
        ids = cell.locateSites(section, distances[1])
    
    
    def set_stimulus(self):
        pass
    
       
    def get_location(self,location):
        #TODO: update so it returns proper section
        return getattr(self.cell, 'get_%s'%location)
    """
    def get_section(self,section,id=0):
        #TODO convert to using section lists
        print section, id
        try:
            print 'A'
            return getattr(self.cell, 'get_%s'%section)()[id]
        except :
            print 'B'
            try:
                return getattr(self.cell, 'get_%s'%section)()
            except:
                print "Fail - could not get_section %s"%section
    """
    def get_section(self,location,id=0):
        #print 'get_section:'
        #print self.cell.get_subdomains()
        #print location, id
        #ls = self.cell.get_subdomain_list('soma')
        #print ls.hname()
        #if True:
        try:
            #print 'A', location
            ls = self.cell.get_subdomain_list(location)
            tmp = ls.object(id)
            #print tmp.hname()
            return tmp
        except:
            print sys.exc_info()[0]
            #print 'B'        
    
    def tag_locations(self):
        """
        Go through and populate Python reference/pointers to given sections in our NEURON morphology
        Warns if ids or distances are invalid
        
        Input:
            names        list of pointer names i.e. ['Site1','Site2']
            sections     list of sections
            ids          list of tuples (ids,sitex) i.e. [(0,0.5),(36,0.97)]. Should be same length as name

        Returns:
            zipped list of valid names and ids
        """
        print 'TAGGING LOCATIONS'
        self.cell.print_section_info()
        try:
            names = self.params['mark_loc']['names']
            sections = self.params['mark_loc']['sections']
            ids = self.params['mark_loc']['ids']
        except:
            print 'No names, sections and ids included to mark for locations'
            return
            
                    
        marked_loc = {}
        """
        i=-1
        try:
            i = sections.index('soma')
        except:
            pass
        if i>-1:
            #process singular: soma
            print sections
            print self.cell.get_soma()[0]
            marked_loc[names[i]] = [sections[i], self.cell.get_soma()[0],0.5]
            # remove soma from names/sections
            sections.pop(i)
            names.pop(i)
            ids.pop(i)
        """    
        print 'we have find for', names, sections, ids
        
        for (i,name) in enumerate(names):
            print '---------------------------------------',i , type(ids[i])
            
            if marked_loc.has_key(name):
                # keep in for when we have excess for dictionary of idx
                print 'Already has entry for ', name
                continue
            if type(ids[i])==dict:
                # take a bit more care with what we do next
                if ids[i].has_key('idx'):
                        
                        sublistnames,idxlist = ids[i]['idx']
                        #section_area = ids[i]['sections'][0]
                        section_segments = self.cell.convert_idx_secids(idxlist,ids[i]['sections'])
                        print "has key idx"
                        print "sublistnames ====> ",sublistnames
                        print len(sublistnames)
                        print len(ids[i]['sections'])
                        
                        for (tmp_id,sbname) in enumerate(sublistnames):
                            marked_loc[sbname] = [ids[i]['sections'][tmp_id], section_segments[tmp_id][1], section_segments[tmp_id][2], section_segments[tmp_id][0]]
                else:
                    print "Error: couldn't work out what kind of dictionary we had"
                        
            elif type(ids[i][0])==int:
                print 'Got id and section'
                marked_loc[name] = [sections[i], self.get_section(sections[i],ids[i][0]), ids[i][1], ids[i][0]]
            elif type(ids[i][0])==str:
                # DANGER WILL ROBINSON!
                #string is a function call to select ... so 
                print 'Got dynamic selection - ', ids[i]
                print "cell.%s(%s)"%(ids[i][0],ids[i][1]),
                print self.cell,
                #secid = getattr(self.cell,ids[i][0])(**ids[i][1])[0]
                ids_list = getattr(self.cell,ids[i][0])(**ids[i][1])
                print "... found %g elements"%len(ids_list)
                random.shuffle(ids_list)
                secid = ids_list[0]
                print secid
                marked_loc[name] = [sections[i], self.get_section(sections[i],secid[0]),  secid[2], secid[0]]
            else:
                print type(ids[i][0]), ids[i][0]
                
            print '---asdasd------------------------------------'
        print 'Have hopefully tagged all locations' #, marked_loc
        self.outputparams['tagged_locations'] = marked_loc
         
    
    def _create_spiketrains(self, stparams):
        #print '\n'*5, 'params = ', stparams
        tstims = stparams['tstims']
        el = stparams['el']
        weights = stparams['weights']
        locations = stparams['locations']
        #print 'tagged locations = ', self.outputparams['tagged_locations']
        
        for i, tstim in enumerate(tstims):
            tagloc = self.outputparams['tagged_locations'][locations[i]]
            #print tagloc, locations[i]
            expsyn = self.add_ExpSyn(section=locations[i], position=tagloc[2], name='Stream'+str(i), tstim=tstim, w=weights[i]*el)
            self.outputparams['stim_spiketrains'][locations[i]] = expsyn
            print 'Added ExpSyn=Stream%g at position %s'%(i,str(locations[i])),'with first 5 times: ',tstim[:5],'...'
        
    
    def add_ExpSyn(self, section='soma', position=0.5, name='default', tstim=[50], w=.001):
        """
        Create/replace an Expsyn synapse on a given section which is active at the time in tstim

        Comments
        --------
        The sort command is here to make sure that tstim are in the right order. This method
        requires the pre-compiling of vecstim.mod by NEURON.

        Function supplied by R. Caze 2013
        """
        expsyn = {}
        expsyn['syn'] = h.ExpSyn()
        location =self.outputparams['tagged_locations'][section][1]
        print 'Adding for name=%s, section=%s'%(name,section)
        expsyn['syn'].loc(position,sec=location)
        expsyn['stim'] = h.Vector(np.sort(tstim)) # Converting tstim into a NEURON vector (to play in NEURON)
        expsyn['vplay'] = h.VecStim() # Creating play vectors to interface with NEURON
        expsyn['vplay'].play(expsyn['stim'])  # Connecting vector to VecStim object to play them
        expsyn['netcon'] = h.NetCon(expsyn['vplay'], expsyn['syn']) # Building the netcon object to connect the stims and the synapses
        expsyn['netcon'].weight[0] = w # Setting the individual weights
        return expsyn

        
    def _create_iclamp(self,params):
        """
        Delivers an current injection using NEURON's IClamp object
        
        """
        print '\n'*3, 'iclamp = params = ', params
        tag_location = params['location']
        location =self.outputparams['tagged_locations'][tag_location][1]
        position =self.outputparams['tagged_locations'][tag_location][2]
        iclamp = h.IClamp()
        iclamp.loc(position, sec=location)
        iclamp.amp = params['amp']
        iclamp.dur = params['duration']
        #st1.del = tstart --> error as del is reserved word in python
        setattr(iclamp, 'del', params['tstim'])
        self.outputparams['stim_iclamp'][tag_location] = iclamp
    
    def _create_iclamp_train(self,params):
        """
        Uses the IPulse3 to supply current to the neuron:
        
        From the IPulse3 documentation:
            Generates a train of current pulses of variable amplitude
            User specifies dur (pulse duration), per (period, i.e. interval 
            between pulse onsets), and num (number of pulses).
            Ensures that period is longer than pulse duration.
        
        """
        print '\n'*5, 'iclamp = params = ', params
        tag_location = params['location']
        
        location =self.outputparams['tagged_locations'][tag_location][1]
        position =self.outputparams['tagged_locations'][tag_location][2]
        clamptrain = {}

        ipulse = h.Ipulse3()
        ipulse.loc(position,sec=location)
        ipulse.dur = params['dur']
        ipulse.amp = params['amp']
        
        clamptrain['stim'] = h.Vector(np.sort(params['tstims']))
        clamptrain['vplay'] = h.VecStim()
        clamptrain['vplay'].play(clamptrain['stim']) 
        clamptrain['netcon'] = h.NetCon(clamptrain['vplay'], ipulse)
        clamptrain['ipulse'] = ipulse
        self.outputparams['stim_iclamp_train'][tag_location] = clamptrain
                
    def _create_epsp(self,params):
        print '\n'*5, 'epsp params = ', params        
        tag_location = params['location']
        location =self.outputparams['tagged_locations'][tag_location][1]
        position =self.outputparams['tagged_locations'][tag_location][2]
        epsp = h.epsp()
        epsp.loc(position, sec=location)
        epsp.tau0 = params['risetau']       
        epsp.tau1 = params['decaytau']
        epsp.onset = params['tstim'] + params['BACdt']
        epsp.imax = params['EPSPamp']
        self.outputparams['stim_epsp'][tag_location] = epsp


    
    def old_create_iclamp(self,params, location=None):
        print '\n'*5, 'iclamp = params = ', params
        pass
        h.st1 = h.IClamp()
        h.st1.loc(0.5, sec=self.get_location(location))
        print '\n'*5, 'params = ', params
        h.st1.amp = params['amp']
        h.st1.dur = params['duration']
        #st1.del = tstart --> error as del is reserved word in python
        st1_delay = params['tstim']
        setattr(h.st1, 'del', st1_delay)
        
        
        
    def old_create_epsp(self,params, location=None):
        print '\n'*5, 'epsp params = ', params        
        pass
        h.syn1 = h.epsp(self.get_location(location))
        h.syn1.tau0 = params['risetau']       
        h.syn1.tau1 = params['decaytau']
        h.syn1.onset = h.tstart + params['BACdt']
        h.syn1.imax = params['EPSPamp']
        #h('L5PC.apic[distSite[0]] { syn1 } ')
    

    def _record_v(self):
        pass
        #h.st1._ref_v
    
    def _record_i(self):
        pass
    
    def _record_g(self):
        pass
        #h.st1._ref_g
    
    def _record_x(self):
        pass
    
        
    def setup_record(self):
        """
        Set up recording devices and structures.

        """
        if self.params['record_loc'] is None or len(self.params['record_loc']) == 0:
            print("No values/locations to record")
            return        
        
        #Record time
        self.rec_t = h.Vector()
        self.rec_t.record(h._ref_t)
        #print self.outputparams.keys()
        # Now go through and start recording voltage/current/etc
        
        
        for (recordable,locations) in self.params['record_loc'].iteritems():
            record_dict = {}
            for loc in locations:
                sec = self.outputparams['tagged_locations'][loc][1]
                posn = self.outputparams['tagged_locations'][loc][2]
                tmprec = h.Vector()
                #if recordable == 'iphoto':
                if recordable == 'iChR' or recordable == 'iNpHR' :
                    print('Attempting to record iphoto:')
                    #fnrec = getattr(sec,'%s'%recordable)
                    
                    """
                    Note that this will only ever record the first section in the Neuron object. 
                    As this is typically the soma, these results therefore only ever hold for the soma
                    
                    """
                    #TODO see if this can be extended for any section 
                        
                    if recordable == 'iChR':
                        print('Recording iChR')
                        #tmprec.record(h.ref(v.o(0).iChR)) # record in soma
                        tmprec.record(h.ChR[0]._ref_iChR)
                    elif recordable == 'iNpHR':
                        print('Recording iNpHR')
                        #tmprec.record(h.ref(v.o(0).iNpHR))
                        tmprec.record(h.NpHR[0]._ref_iNpHR)
                    else:
                        print 'Could not record photocurrent for %s'%recordable
                    
                        #print tmprec
                        #tmprec.record(h.ChR[0]._ref_iChR)
                        #h('tmprec.record(&ChR[0].iChR,i_ChR,rec_t)')
                        #break
                       
                        #tmp = getattr(firstobj,'iNpHR')
                    #tmprec.record(tmp)  
                else:
                    fnrec = getattr(sec(posn),'_ref_%s'%recordable)
                    tmprec.record(fnrec)
                record_dict['%s_%s'%(recordable,loc)] = tmprec
                
    
                print('Added %s_rec for location=%s _____________ as %s_%s'%(recordable,loc,recordable,loc))
            self.outputparams['%s_rec'%recordable] = record_dict
            #print self.outputparams
        
        
        
        """
        # for each opsin, include vector for corresponding photocurrent
        opsindict = self.params['opdict']
        h('access L5PC.soma')
        for (opsin,opsinval) in opsindict.iteritems():
            print 'testing for %s'%opsin, opsinval
            if self.is_opsin_present(opsinval,opsin):
                h('objref i_%s'%opsin)
                h('i_%s = new Vector()'%opsin)
                h('cvode.record(&%s[0].i%s,i_%s,tvec)'%(opsin,opsin,opsin)) 
                print 'Recording from i_%s'%opsin
                
        """ # TODO remove once debugged
        """
        self.rec_t = h.Vector()
        self.rec_t.record(h._ref_t)
        # Record Voltage
        self.rec_v = h.Vector()
        self.rec_v.record(self.cell.soma(0.5)._ref_v)
        """
        
        
        """
            tmpi = h.Vector()
            tmpi.record(sec(0.5)._ref_ina)
            tmpi2 = h.Vector()
            tmpi2.record(sec(0.5)._ref_ik)
            #record_vrefs['g_%s'%loc] = h.Vector()
            #record_vrefs['g_%s'%loc].record(sec(0.5)._ref_v)
        """        
        
        
        
            
#         
#         
#         
#         #======================== recording settings ============================
#         h('objref vsoma, vdend, proxClamp, vdend2, isoma, gsoma, isoma_k, isoma_na')
#     
#         h.vsoma = h.Vector()
#         h.isoma = h.Vector()
#         h('access L5PC.soma')
#         h('cvode.record(&v(0.5),vsoma,tvec)')
#         #h.cvode.record(h.st1._ref_i,h.isoma,h.tvec)
#         h.isoma_k = h.Vector()
#         h.isoma_na = h.Vector()
#         h('cvode.record(&ik(.5),isoma_k,tvec)')
#         h('cvode.record(&ina(.5),isoma_na,tvec)')
#         
#         
#         """
#         # TODO: would also like to record conductances from soma, etc.
#         #h.gsoma = h.Vector()
#         #h.cvode.record(h.st1._ref_g,h.gsoma,h.tvec)
#         
#         # Record voltage at distal point on apical tree
#         h.vdend = h.Vector()
#         h('access L5PC.apic[distSite[0]]')
#         h.cvode.record(h.L5PC.apic[int(h.distSite[0])](h.distSite[1])._ref_v,h.vdend,h.tvec)
#         
#     
#         # Record voltage at proximal point on apical tree
#         h('access L5PC.apic[proxSite[0]]')
#         h.vdend2 = h.Vector()
#         h.cvode.record(h.L5PC.apic[int(h.proxSite[0])](h.proxSite[1])._ref_v,h.vdend2,h.tvec)
#     
#         # Following lines are only useful when we're wanting to plot location of proxSite
#     #    h('access L5PC.apic[proxSite[0]]')
#     #    h.proxClamp = h.IClamp(h.proxSite[1])
#     #    h.proxClamp.amp = 0
#     #    h('L5PC.apic[proxSite[0]] { proxClamp }')
#         """
#         
#         
        
    def _run_simulation(self,cell, variable_dt=False, atol=0.001):
        '''
        Running the actual simulation in NEURON, simulations in NEURON
        is now interruptable.
        '''
        print 'A'
        cell.tstopms = 50
        cell.timeres_NEURON = 0.01
        cell.v_init = -65.
        cell.tstartms = 0.
        print 'B'
        h.dt = cell.timeres_NEURON
        
        cvode = h.CVode()
        
        #don't know if this is the way to do, but needed for variable dt method
        if variable_dt:
            cvode.active(1)
            cvode.atol(atol)
        else:
            cvode.active(0)
        print 'C'
        #initialize state
        h.finitialize(cell.v_init)
        print 'D'
        #initialize current- and record
        if cvode.active():
            cvode.re_init()
        else:
            h.fcurrent()
        h.frecord_init()
        print 'E'
        ##Starting simulation at tstart
        h.t = cell.tstartms
        print 'F'
        ##################################################cell._loadspikes()
        print 'G'
        #print sim.time and realtime factor at intervals
        counter = 0.
        t0 = time.time()
        ti = h.t
        print 'H'
        if cell.tstopms > 1000:
            interval = 1 / cell.timeres_NEURON * 100
        else:
            interval = 1 / cell.timeres_NEURON * 10
        print 'I'
        while h.t < cell.tstopms:
            #print 'J'
            h.fadvance()
            #print 'K'
            counter += 1.
            #print 'L'
            if np.mod(counter, interval) == 0:
                #print 'M'
                rtfactor = (h.t - ti)  * 1E-3 / (time.time() - t0)
                print('t = %.0f, realtime factor: %.3f' % (h.t, rtfactor))
                t0 = time.time()
                ti = h.t
        
        
        
    
    def simulate_exp(self):
        print 'Got here'
        print "Init experiment ..."
        #self.cell.initialise()
        """
        init()
        run(self.params['tstop'])
        
        """
        h.init()
        h.run()
        
        
        print "Init successful. About to run experiment ..."
        try:
            pass
        except:
            print 'eerrrer'
        print "...finished"
        
    def is_opsin_present(self,opsinvalue, opsin):
        """
        opsinvalue = opsindict[opsin]
        """
        #print opsin, opsinvalue[0][0],'======================================================= is Opsin present? ', 
        #bb = opsinvalue[0][0] is not None
        #print bb
        return opsinvalue[0][0] is not None
        
        
    def num_recording_areas(self,typerec):
        print self.params.keys()
        if self.params.has_key('record_loc'):
        
            if self.params['record_loc'] is not None and self.params['record_loc'].has_key(typerec):
                return len(self.params['record_loc'][typerec])
        
        return 0
    """
    def get_recording_prefixs(self,prefix):
        pass
    
    def _save_data(self,prefix):
        prefixs = self.get_recording_prefixs(prefix)
        if len(prefixs) == 0:
            print('No data to save for prefix %s!'%prefix)
            return
        mat = np.zeros((len(prefixs)+1,int(len(self.rec_t))))
        mat[0,:] = self.rec_t
        for (i,rec) in enumerate(prefixs):
            mat[i+1,:] = self.outputparams['%s_rec']['%s_%s'%(prefix,self.params['record_loc']['v'][i])]
    """   
    
    def save_data(self):#,expname,savedata=False,opsindict={}):
        """
        
        
        """
        
        print('==================')
        print self.params['record_loc']
        print self.outputparams
        print('============record_loc======')
        
        
        expname = self.params['expname']
        if "savedata" not in self.params or not self.params['savedata']: 
            print("Not saving data")
            return
        
        import numpy as np
        # Updated: 06/10/15 - 
        # save recorded voltages, currents, etc.
        # new stronger method for saving to file
        # save all voltages in same file
        # same for all currents, in their separate i
        # same for conductances

        # Version 2 06/10/15 - save everything into its own file with extension. Sigh. Much easier
        if self.params.has_key('record_loc') and self.params['record_loc'] is not None:
            for (k,v) in self.params['record_loc'].iteritems():
                if len(v)==0:
                    continue
                # else print
                mat = np.zeros((2,int(len(self.rec_t))))
                mat[0,:] = self.rec_t
                mat[1,:] = self.outputparams['%s_rec'%k]['%s_%s'%(k,self.params['record_loc'][k][0])]
                np.savetxt(expname+"_%s.dat"%k,mat)
                
        
        """
        #mat = np.matrix([self.rec_t])
        mat = np.zeros((self.num_recording_areas('v')+1,int(len(self.rec_t))))
        mat[0,:] = self.rec_t
        for i in range(self.num_recording_areas('v')):
            mat[i+1,:] = self.outputparams['v_rec']['v_%s'%self.params['record_loc']['v'][i]]
        np.savetxt(expname+"_v.dat",mat)
        # TODO also save tag names as header info
        
        
        print self.outputparams.keys()
        try:
            print self.outputparams['iphoto_rec']
        
            # hack to make sure that we're recording iPhoto ok
            mat = np.zeros((2,int(len(self.rec_t))))
            mat[0,:] = self.rec_t
            mat[1,:] = self.outputparams['iphoto_rec']['iphoto_mysoma']
            np.savetxt('test_iphoto.txt',mat)
        except:
            print 'error with saving data ... '
            
        """
        #if len(self.params['opdict'].keys())==0:
        #######    return
        """
        for opsin in self.params['opdict'].keys():
            if not self.is_opsin_present(self.params['opdict'][opsin],opsin):
                continue
            h('objref savdata')
            h('savdata = new File()')
            h('savdata.wopen("%s_i%s.dat")'%(expname,opsin))
            h('objref tempmatrix')
            h('tempmatrix = new Matrix()')
            h('tempmatrix.resize(tvec.size(),2)')
            h('tempmatrix.setcol(0, tvec)')
            h('tempmatrix.setcol(1, i_%s)'%opsin)
            h('tempmatrix.fprint(savdata, " %g")')
            h('savdata.close()')
        """ # TODO opsin remove
    
    def old_save_data(self,expname,savedata=False,opsindict={}):
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
    
    
    def plot_optogenetics(self,max_val=40,min_val=-100):
        
        if self.params.has_key('opsindict'):
            #print self.params['opsindict']
            for (opsin,oplocations) in self.params['opsindict'].iteritems():
                for (oploc,opexp) in oplocations.iteritems():
                    #print opexp
                    if oploc.lower() == 'none':
                        continue
                    for pulsenum in range(opexp['n_pulses']):                
                        pulse_start = opexp['lightdelay']+pulsenum*(opexp['pulsewidth']+opexp['interpulse_interval'])
                        yoffset = (max_val+min_val)/2.
                        self.plot_optogenetic(opsin,pulse_start,pulse_start+opexp['pulsewidth'],yoffset=yoffset)
                    # once we've plotted an activation for one area, that should be sufficient i.e. we don't need to plot apical *and* soma, only the first 
                    # TODO: think how to extend this to allow for different areas to be indicated i.e. ChR in soma vs ChR in apical dendritic arbor
                    break
    
    def plot_optogenetic(self,opsin_type,ton,toff,yoffset):
        """ 
        Plots lines at top of figure to indicate when each opsin was active 
        Y-value and color set in static dict 'opsin_dict' 
        """
        yval = opsin_dict[opsin_type]['ycoord']
        color = opsin_dict[opsin_type]['color']
        print 'plotting opsin --> ', opsin_type, color
        ax = pylab.gca()
        x = [ton,toff]
        y = np.array([yoffset+yval,yoffset+yval])
        ax.fill_between(x,y-2,y+2,color=color)
    
    def _get_vector(self,lineinfo):
        name = lineinfo[0]
        electrode = lineinfo[1]
        print electrode,name
        return self.outputparams['%s_rec'%electrode]['%s_%s'%(electrode,name)]
    
    def run_plots(self):
        lw = 2
        
        if "plot" not in self.params or len(self.params['plot'].keys())==0:
            print('No plots to create')
            return
        
        for figid in self.params['plot'].keys():
            print figid
            pylab.figure()
            
            for line in self.params['plot'][figid]:
                print line
                if len(line)==3: #means we had a line description
                    pylab.plot(self.rec_t,self._get_vector(line),line[2])
                else:
                    pylab.plot(self.rec_t,self._get_vector(line))
            print "\nMax was ", self._get_vector(line).max(),'\n'
            max_val = self._get_vector(line).max()
            min_val = self._get_vector(line).min()
            pylab.xlim(self.rec_t[0]-20,self.rec_t[-1]+20)
            ax = pylab.gca()
            for loc, spine in ax.spines.iteritems():
                if loc in ['left','bottom']:
                    spine.set_position(('outward',5))
                    ax.tick_params(direction='out')
                elif loc in ['right','top']:
                    spine.set_color('none') 
            pylab.xlabel('time (ms)')
            
            self.plot_optogenetics(max_val=max_val,min_val=min_val)
            
            if self.params['expname'] is not None:
                savename = '%s_fig%s.png'%(self.params['expname'],figid)
                pylab.savefig(savename)
                print "Saved figure under %s*.png"%savename
                pylab.close('all')
            else:
                pylab.show()
        
  
            
    
    def old_run_plots(self,params):
        """
        Generates plots for voltage and current at soma, distal and proximal locations on apical dendrites
        """
        lw = 2
        
        
        # Plot voltage at soma and dendrites (apical proximal and distal)
        pylab.figure(1)
        pylab.plot(h.tvec,h.vsoma,lw=lw,c='k',label='v_soma')
        #pylab.plot(h.tvec,h.vdend,lw=lw,c='r',label='v_dend')
        #pylab.plot(h.tvec,h.vdend2,lw=lw,c='b',label='v_dend2')
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
        
        """
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
        """
        
        if params['expname'] is not None:
            savename = params['expname']
            pylab.figure(1)
            pylab.savefig(savename+'_voltage.png')
            #pylab.figure(2)
            #pylab.savefig(savename+'_current.png')
            print "Saved figures under %s*.png"%savename
            pylab.close('all')
        else:
            pylab.show()
            
    def save_outputfile(self):
        
        if self.params['output']['output_pkl']:
            self.outputparams['time_saved'] = time.strftime('%y/%m/%d:%H.%M')
            for (k,v) in self.outputparams.iteritems():
                print k, '-->',v
            fio.saveresults(self.outputparams, self.params['expname'])
            
        if self.params['output'].has_key('neu'):
            h.load_file("nrn/hoc/neu_tree.hoc")

            try:
                name = self.params['output']['neu']
                if not name.endswith('.neu'):
                    name = name+'.neu'
                h.neu_tree(name)
                print('Saved neuron to NEURON .hoc format : %s'%name)
            except:
                print('Error with trying to save supplied name')
                h.neu_tree('%s.neu'%self.params['expname'])
        

    def print_stats(self):
        # TODO: implement
        pass
        #print self.outputparams['v_rec']
    
    def printstatus(self,msg):
        print '='*30
        print msg
        print '='*30
        
    
    def main(self,exp_keydict):
        print '\n'*2, '='*40, exp_keydict['expname']
        
        self.params.update(exp_keydict)
        self.printstatus('Setting up')
        self.setup()
        self.set_experiment_type()
        self.printstatus('Setting stimulus')
        self.set_stimulus()
        self.new_set_stimulus()
        self.printstatus('Setting optogenetics')
        self.add_optogenetics()
        #self.add_optogenetics_old(self.params['opdict'])
        self.printstatus('Setting recording')
        self.setup_record()
        self.printstatus('Simulating')
        self.simulate_exp()
        self.printstatus('Saving data')
        self.save_data()
        self.printstatus('Running plots')
        self.run_plots()
        self.printstatus('Saving output')
        self.save_outputfile()
        self.printstatus('Printing stats')
        self.print_stats()
    
    
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
