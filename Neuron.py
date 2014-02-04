from neuron import nrn, h, run, init
#from nrn    import *
import sys
import numpy as np
import math
from LFPy import Cell

class NeuronInterface(Cell):
    
    def __init__(self):
        self._create_sectionlists()
        self.totnsegs = self._calc_totnsegs()
        self.verbose = True
        self._collect_geometry()
        try:
            self.generate_seg_section()
        except:
            pass
    
    def get_cell(self):
        raise NotImplementedError( "Not implemented")

    def get_soma(self):
        raise NotImplementedError( "Not implemented")
 
    def get_axon(self):
        raise NotImplementedError( "Not implemented")
    
    def get_basal(self):
        raise NotImplementedError( "Not implemented")
    
    def get_apical(self):
        raise NotImplementedError( "Not implemented")
    
    def locateSites(self,section,distance):
        raise NotImplementedError( "Not implemented")
       
    def select_segments_bydistance(self,sectionarea,mindist=0,maxdist=-1):
        """
        Returns segments by distance
        
        """
        somaidx = self.somaidx[0]
        #print somaidx
        idxs = self.get_idx(sectionarea)
        #print idxs
        dists = np.array([self.get_intersegment_distance(somaidx, idx) for idx in idxs])
        #print dists
        if maxdist <=0:
            segments = np.where(dists >= mindist)
        else:
            segments = np.where((dists >= mindist) & (dists<= maxdist))
        
        return segments
    
    def select_section_posn_bydistance(self,sectionarea,mindist=0,maxdist=-1):
        #print self.seg_sec.keys()
        segments = self.select_segments_bydistance(sectionarea, mindist, maxdist)
        #print segments
        offset = self.get_idx(sectionarea)[0]
        secids = [self.seg_sec[sectionarea][seg+offset] for seg in segments[0]]
        #print secids
        return secids
    
    """
    def get_section_from_segment(self,sectionarea,idx):
        secname = self.get_idx_name(int(idx))[1]
    """
        
    def generate_seg_section(self):
        #print 'generating segsection'
        self.seg_sec = {}
        for subdom in self.get_subdomains():
            #print 'subdom = ', subdom
            self.seg_sec[subdom] = {}
            for idx in self.get_idx(subdom):
                idxname = self.get_idx_name(int(idx))[1]
                #print '--', idx, idxname
                sec = [(i,s) for (i,s) in enumerate(self.get_subdomain_list(subdom)) if s.name() == idxname]
                #print sec
                self.seg_sec[subdom][int(idx)] = (sec[0][0],sec[0][1], self.get_idx_name(int(idx))[2])
        #print self.seg_sec
    
    def get_section_byname(self,sectioname,subdomain):
        sec =  [(i,s) for (i,s) in enumerate(self.get_subdomain_list(subdomain)) if s.name() == sectioname]
        if len(sec) == 0:
            return None
        else:
            return [(sec[0][0],sec[0][1],0.5)]
    
    
    def get_max_distance(self,section):
        somaidx = self.somaidx[0]
        idxs = self.get_idx(section)
        
        dists = np.array([self.get_intersegment_distance(somaidx, idx) for idx in idxs])
        return np.max(dists)
    
    def get_subdomains(self):
        # TODO we could actually automate this so that we can pull the names out of this.allsecnames
        return []
    
    def get_subdomain_list(self,name):
        #print 'get_subdomain_list: Was asked for name=',name
        #print 'type = ',self.domainlists[name].hname()
        try:
            return self.domainlists[name]
        except:
            raise NotImplementedError( "Not implemented")

    def print_section_info(self):
        for subdom in self.get_subdomains():
            print 'Subdomain list %s has %g elements'%(subdom,self.domainlists[subdom].count())


    def initialise(self):
        print 'Initialization of cell type not implemented. If segmentation faults exist, could be a reason why ... '
        pass

class L5PC(NeuronInterface):
    
    def get_default_params(self):
        return {'model':'L5PCbiophys3.hoc','cellmorph':'cell1.asc'}
   
    def __init__(self,inparams={}): 
        params = self.get_default_params()
        params.update(inparams)
        h.load_file('stdlib.hoc') 
        h.load_file('stdrun.hoc')
        h.load_file("import3d.hoc")
        h.load_file("models/%s"%params['model'])
        h.load_file("models/L5PCtemplate.hoc")
        
        self.cell = h.L5PCtemplate("morphologies/%s"%params['cellmorph'])
        super(L5PC,self).__init__()
        self.populate_sectionlists()
        # note that we have to regenerate seg_sections
        self.generate_seg_section()
        print 'Created L5PC neuron'
        
    def get_cell(self):
        return self.cell
    
    def get_soma(self):
        return self.cell.soma
    
    def get_basal(self):
        return self.cell.dend
    
    def get_apical(self):
        return self.cell.apic
    
    def get_apic(self):
        return self.cell.apic

    def get_axon(self):
        return self.cell.axon
    
    def locateSites(self, sectionString, distance):
        return self.cell.locateSites(sectionString,distance)
    
    
    def populate_sectionlists(self):
        self.domainlists = {}
        for subd in self.get_subdomains():
            self.domainlists[subd] = h.List() #h.SectionList()
        
        print self.domainlists['soma'].hname()
        
        for (i,sec) in enumerate(self.allseclist):
            try:               
                # this is hacky but we have to do it this way, as a typical name is L5Template[3].apic[34] 
                matched_domain = [subd for subd in self.get_subdomains() if sec.name().find(subd)>=0][0]
                #self.domainlists[matched_domain].append(sec=sec) ### fine only when list is a SectionList
                self.domainlists[matched_domain].append(sec)
                #print 'Added another to subdomain %s'%(matched_domain)
            except:
                print '%s does not seem to belong to any subdomain'%sec.name()
                continue
                
    def get_num_segments(self,section):
        """ SUPERCEDED BY LFPy.Cell.
        """
        try:
            count = getattr(self.cell, 'nSec%s'%(section.capitalize()))
            return int(count)
        except:
            print 'Section %s unknown; returning 0 for section count'%section
            return 0 

    def get_subdomains(self):
        return ['soma', 'axon', 'apic','dend']



# Mechanism and Section are wrappers from Andrew Davidson neuronpy modified and commented by Romain Caze
# They enable to interface NEURON and Python more easily.
class Mechanism(object):
    """
    Create mechanism which will be inserted in the membrane of a cell

    Examples
    --------
    >>> leak = Mechanism('pas', {'e': -65, 'g': 0.0002})
    >>> hh = Mechanism('hh')
    """
    def __init__(self, name, parameters={}):
        """
        Parameters
        ----------
        name: a char
            the label of the mechanism
        parameters: a dictionary
            contains the different parameter of a mechanism
        """
        self.name = name
        self.parameters = parameters

    def insert_into(self, section):
        """
        Method used to insert a mechanism into a section

        Parameters
        ----------
        section: a NEURON section
            the section where the mechanism needs to be inserted

        """
        section.insert(self.name)
        for name, value in self.parameters.items():
            for segment in section:
                mech = getattr(segment, self.name)
                setattr(mech, name, value)



class Section(nrn.Section):
    """
    Create a NEURON section with certain mechanism inserted in it

    Examples
    --------
    >>> soma = Section(L=30, diam=30, mechanisms=[hh, leak])
    >>> apical = Section(L=600, diam=2, nseg=5, mechanisms=[leak],
    ...                  parent=soma, connection_point=0)
    """
    def __init__(self, L, diam, nseg=1, Ra=100, cm=1, mechanisms=[], parent=None, connection_point=1,name=None,startpt=None,endpt=None):        
        """
        Parameters
        ----------
        L: a float
            length in micrometers
        diam: a float
            diameter in micrometers
        nseg: an int
            number of segements
        Ra: a float
            surface axial resistance Ohm/micrometer square
        cm: a float
            capacitance in F/micrometer square
        mechanisms: a list
            mechanisms to be inserted (need to be created beforehand)
        parent: a NEURON section
            section to which this section is coming from
        connection_point: a float between 0 and 1
            where this section is connected to its parent

        """
        
        if name is not None:
            nrn.Section.__init__(self,name=name)
        else:
            nrn.Section.__init__(self)
        # set geometry
        if startpt is None and endpt is None:
            self.L = L
        else:
            h.pt3dadd(startpt[0],startpt[1],startpt[2],diam,sec=self)
            h.pt3dadd(endpt[0],endpt[1],endpt[2],diam,sec=self)
        self.diam = diam
        self.nseg = nseg
        # set cable properties
        self.Ra = Ra
        self.cm = cm
        # connect to parent section
        if parent:
            self.connect(parent, connection_point, 0)
        # add the mechanisms
        for mechanism in mechanisms:
            mechanism.insert_into(self)
        
            
    """def name(self):
        return self.aname
    """
    def record_spikes(self, threshold=-30):
        """
        Record the number of spikes produced in this section,
        which is the number of time a voltage is crossed in
        the middle of a section

        Parameters
        ----------
        threshold: a float
            voltage determining the presence or not of a spike

        Returns
        -------
        nothing, but change the self.spikecount

        """
        self.spiketimes = h.Vector()
        self.spikecount = h.APCount(0.5, sec=self)
        self.spikecount.thresh = threshold
        self.spikecount.record(self.spiketimes)


class SimpleNeuron(NeuronInterface):
    
    def __init__(self):
        super(SimpleNeuron,self).__init__()
        
    
    def create_soma(self,soma_params):
        return Section(**soma_params) # soma_params['x'], soma_params['y'], Ra=soma_params['Ra'], mechanisms=soma_params['mechs'],name='soma')
    
    def create_dendrites(self,dend_params,label):
        #return Section(dend_params['dlength'], dend_params['ddiam'], nseg=dend_params['nseg'], Ra=dend_params['Ra'], parent=dend_params['parent'], mechanisms=dend_params['mechs'],name=label)
        return Section(**dend_params)
    
    def get_soma(self):
        return self.soma

    def get_cell(self):
        return self
    
    def get_dend(self,dend_id):
        return getattr(self,'dend%g'%dend_id)
    

    def initialise(self, vrest=-65):
        '''
        Initialise the model, to launch before each simulation
        '''
        print 'initialized the model'
        for sec in h.allsec():
            h.finitialize(vrest, sec)
            h.fcurrent(sec)
        h.frecord_init()
    
    def get_default_params(self):
        hh = Mechanism('hh')
        pas = Mechanism('pas', {'e':-65,'g':0.0001})
        params = {'soma': {'L':10,
                           'diam':10,
                           'Ra':150,
                           'name':'soma',
                           'startpt':[0,0,0],
                           'endpt': [0,0,0],
                           'mechanisms':[hh,pas]}, 
                  'dend0': {'L':50.,
                            'name':'dend0',
                            'diam'  : 0.4,
                            'nseg':  10,
                            'Ra' : 150,
                            'startpt':[0,0,0],
                            'endpt': [50,0,0],
                            'mechanisms':[pas]},
                  'dend1': {'L':50.,
                            'name':'dend1',
                            'diam'  : 0.4,
                            'nseg':  10,
                            'Ra' : 150,
                            'startpt':[0,0,0],
                            'endpt': [-50,0,0],
                            'mechanisms':[pas]},}
        return params

    
class SimpleBipolar(SimpleNeuron):
    
    
    def __init__(self,inparams={}):
        params = self.get_default_params()
        params.update(inparams)
        self.create_bipolar(params)
        self.subdomains = ['soma',  'dend0','dend1'] 
        self.params = params
        self.populate_sectionlists()
        super(SimpleBipolar,self).__init__()
        # can be written as 
        #super().__init__()
        # when we get to Python 3.0
        
    def create_bipolar(self,params):
        self.soma = self.create_soma(params['soma'])
        for i in range(2):
            label = 'dend%g'%i
            params[label]['parent'] = self.soma
        # TODO: make this dynamic --> SimpleNPolar
        self.dend0 = self.create_dendrites(params['dend0'],label='dend0')
        self.dend1 = self.create_dendrites(params['dend1'],label='dend1')

    def count_dendrites(self):
        return 2        
        
    def get_apical(self):
        return self.dend0
    
    def get_basal(self):
        return self.dend1
    
    def get_num_segments(self,section):
        return self.subdomains.count(section)

    def get_subdomains(self):
        return self.subdomains   
    
    def populate_sectionlists(self):
        self.domainlists = {}
        for subd in self.get_subdomains():
            self.domainlists[subd] = h.List()
        self.domainlists['soma'].append(self.soma)
        self.domainlists['dend0'].append(self.dend0)
        self.domainlists['dend1'].append(self.dend1)
        """        
        self.domainlists['soma'].append(sec=self.soma)
        self.domainlists['dend0'].append(sec=self.dend0)
        self.domainlists['dend1'].append(sec=self.dend1)
        """
        
class SimpleUnipolar(SimpleNeuron):
    

    def __init__(self,inparams={}):
        params = self.get_default_params()
        params.update(inparams)
        self.create_unipolar(params)
        self.subdomains = ['soma',  'dend0']   
        self.params = params
        self.populate_sectionlists()
        super(SimpleUnipolar,self).__init__()
        
        
    def create_unipolar(self,params):
        self.soma = self.create_soma(params['soma'])
        params['dend0']['parent'] = self.soma
        self.dend0 = self.create_dendrites(params['dend0'],label='dend0')
        
    def count_dendrites(self):
        return 1
    
    def get_apical(self):
        return self.dend0

    def get_num_segments(self,section):
        return self.subdomains.count(section)

    def get_subdomains(self):
        return self.subdomains   

    def populate_sectionlists(self):
        self.domainlists = {}
        for subd in self.get_subdomains():
            self.domainlists[subd] = h.List()
        
        self.domainlists['soma'].append(self.soma)
        self.domainlists['dend0'].append(self.dend0)
    
    
#TODO : could extend this to SimpleNPolar (i.e. n number of child dendrites))




class FractalNeuron(SimpleNeuron):
    
    def get_default_params(self):
        #pp = super.get_default_params()
        pp = super(FractalNeuron,self).get_default_params()
        
        pas = Mechanism('pas', {'e':-65,'g':0.0001})
        params = {'defaultlevel': {'L':50.,
                            'diam'  : 0.4,
                            'nseg':  10,
                            'Ra' : 150,
                            'mechanisms':[pas]},}
        pp.update(params)
        pp['num_base']=1
        pp['num_split']=2
        pp['num_levels']=4
        return pp
    
    
    
    def __init__(self,inparams={}):
        
        
        
        params = self.get_default_params()
        params.update(inparams)
        self.params = params
        # sectionlists
        self.domainlists = {}
        # create soma and populate corresponding domainlist
        self.soma = self.create_soma(params['soma'])
        self.domainlists['soma'] = h.List()
        self.domainlists['soma'].append(self.soma)
        # create dendrites
        dendnames = self.__create_fractal_dendrites(params['num_base'],params['num_split'],params['num_levels'])
        # set subdomain names
        self.subdomains = ['soma'] + dendnames
        super(FractalNeuron,self).__init__()
        
        print "Created Fractal neuron"
        
        # NORMALLY WOULDN'T HAVE THIS
        #self.draw_fractal_tree('%s_nb%g_ns%g_nl%g'%('131220_confirm_fractal',params['num_base'],params['num_split'],params['num_levels']))
        
    def __create_fractal_dendrites(self,num_base,num_split,num_levels):
        
        dendnames = []
        for i in range(num_base):
            dn = 'dend%g'%i
            print 'creating ',dn
            # create domainlist
            self.domainlists[dn] = h.List()
            # add to domainnames
            dendnames.append(dn)
            # create section
            if self.params.has_key('level0'):
                pp = self.params['level0']
            else:
                pp = self.params['defaultlevel']
            
            position,max_subtend_theta,theta = self._get_firstchild_fractal_position(i,num_base,pp['L'])
            print 'base_%g = '%i,position,theta
            
            pp['parent'] = self.soma
            pp['startpt'] = [0,0,0]
            pp['endpt'] = position
            pp['name'] = 'dend%g'%i
            
            setattr(self,'dend%g'%i,Section(**pp))
            print 'creating ---'
            # recursive limiting case
            print 'creating ---'
            # recursive step
            thisbranch = getattr(self,'dend%g'%i)
            self.domainlists[dn].append(thisbranch)
            print 'creating children'
            self.__create_child_dendrite(thisbranch,0,num_levels,num_split,name=dn,theta_parent=max_subtend_theta,posn_parent=position,theta_relative=theta)
            
        print 'bbb',dendnames
        return dendnames
        
        
        
        
    def __create_child_dendrite(self,parent_branch,parent_level,max_level,num_split,name,theta_parent,posn_parent,theta_relative):
        #print('child_dendrite: ', parent_level, max_level)
        if parent_level+1 >= max_level:
            #print('oops, should stop here')
            return 
        #pp = self.params['level%g'%(parent_level+1)]
        if self.params.has_key('level%g'%(parent_level+1)):
            pp = self.params['level%g'%(parent_level+1)]
        else:
            pp = self.params['defaultlevel']
        
        for j in range(num_split):
            nn = name+'_%g'%j 
            print '-'*(parent_level+1),'creating section ', parent_level+1, j, '----------------------',nn
            # create section, with this branch as parent
            
            position,max_subtend_theta,theta = self._get_child_fractal_position(j,num_split,length=pp['L'],parent_theta_subtend=theta_parent,parent_position=posn_parent,theta_relative=theta_relative)
            print 'child_%s = '%nn,position,theta
            
            pp['parent'] = parent_branch
            pp['startpt'] = posn_parent
            pp['endpt'] = position
            pp['name'] = nn
            
            dj =  Section(**pp)
            # add to domain list
            parent_family = name.split('_')[0]
            self.domainlists[parent_family].append(dj)
            # then do the fractal step
            self.__create_child_dendrite(dj,parent_level+1,max_level,num_split,name=nn,theta_parent=max_subtend_theta,posn_parent=position,theta_relative=theta)
        
    def _get_firstchild_fractal_position(self,i,num_siblings,length):
        origin = np.zeros(3)
        return self._get_fractal_position(i,num_siblings,parent_theta_subtend=2*math.pi,parent_position=origin,length=length,theta_relative=0,use_theta_offset=False)
    
    def _get_child_fractal_position(self,i,num_siblings,parent_theta_subtend,parent_position,length,theta_relative):
        return self._get_fractal_position(i,num_siblings,parent_theta_subtend,parent_position,length,theta_relative)
    
    def _get_fractal_position(self,i,num_siblings,parent_theta_subtend,parent_position,length,theta_relative,use_theta_offset=True):
        """
        
        """
        print 'in get_fractal_posn:'
        # angle relative to origin. Important for grand+children of soma, as everything else is done relative to this
        # note that if we're at first child, then parent_position = [0,0,0], which gives theta_relative = 0, which is fine.
        #theta_relative = np.arctan2(parent_position[1],parent_position[0]) 
        print 10*' ','theta rel = ', theta_relative
        print 10*' ','parent theta = ', parent_theta_subtend
        
        # theta offset:  
        #theta_offset = 0 # TODO give proper value. Should be max of 90, 
        if use_theta_offset and num_siblings > 1:
            theta_offset = parent_theta_subtend*0.5*(num_siblings-1)/num_siblings
        else:
            theta_offset = 0
        print 10*' ','theta_offset = ', theta_offset
        angle = i*parent_theta_subtend/num_siblings + theta_relative - theta_offset
        #angle = i*parent_theta_subtend/num_siblings + theta_relative + theta_offset
        print 10*' ','base angle = ', i*parent_theta_subtend/num_siblings
        print 10*' ','updated angle = ', angle
        new_position = parent_position + np.array([np.cos(angle)*length, np.sin(angle)*length, 0.])
        print 10*' ','new position = ',new_position
        max_subtend = parent_theta_subtend/(num_siblings+1)
        return new_position,  max_subtend, angle
        
        
    def get_subdomains(self):
        return self.domainlists
    
    def populate_sectionlists(self):
        # don't need to do anything as this performed while elements are being created
        pass
        
    def draw_fractal_tree(self,savename):
        from matplotlib.collections import PolyCollection
        import matplotlib.pyplot as plt
        
        
        
        x,y,z,r = self._collect_pt3d()
        #ss = np.ones(len(x))
        zips = []
        for (i,pt) in enumerate(x):
            zips.append(zip(pt,y[i]))
        print zips
        #z = zip(x, y)
        polycol = PolyCollection(zips,
                                 edgecolors='gray', linewidths=2,
                                 facecolors='gray', closed=False)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax.add_collection(polycol)
        ax.scatter(x,y)
        #for z in polys:
        #    ax.plot(z[0])
        ax.axis(ax.axis('equal'))
        
        #plt.show()
        fig.savefig('%s.png'%savename)

        
        
        
        

class Simulation(object):
    """
    Create and control a simulation

    Example
    -------
    >>> cell = BipolarNeuron()
    >>> sim = Simulation(cell)
    >>> sim.go()
    >>> sim.show()
    """
    def __init__(self, cell, sim_time=100, dt=0.01):
        '''
        Parameters
        ----------
        cell: a BipolarNeuron object
            the cell which will be simulated
        sim_tim: an integer
            the time of the simulation
        dt: a float
            the integration timestep (fix)
        '''
        self.cell = cell
        self.sim_time = sim_time
        self.dt = dt
        self.syn, self.stim, self.vplay, self.netcon = {}, {}, {}, {}

    def add_Exp2Syn(self, section='soma', position=0.5, name='default', tstim=[50], w=.001):
        """
        Create/replace an Exp2syn synapse on a given section which is active at the time in tstim

        Comments
        --------
        The sort command is here to make sure that tstim are in the right order. This method
        requires the pre-compiling of vecstim.mod by NEURON.

        """
        self.syn[name] = h.Exp2Syn(self.cell.__getattribute__(section)(position))
        self.stim[name] = h.Vector(np.sort(tstim)) # Converting tstim into a NEURON vector (to play in NEURON)
        self.vplay[name] = h.VecStim() # Creating play vectors to interface with NEURON
        self.vplay[name].play(self.stim[name])  # Connecting vector to VecStim object to play them
        self.netcon[name] = h.NetCon(self.vplay[name], self.syn[name]) # Building the netcon object to connect the stims and the synapses
        self.netcon[name].weight[0] = w # Setting the individual weights

    def set_IClamp(self, name='IClamp', delay=1, amp=-1, dur=3):
        """
        Create a current clamp point process on the soma

        """
        stim = h.IClamp(self.cell.soma(0.5))
        stim.delay = delay
        stim.amp = amp
        stim.dur = dur
        self.stim[name] = stim
    """
    def show(self):
        x = np.array(self.rec_t)
        y = np.array(self.rec_v)
        plt.plot(x, y)
        plt.xlabel("Time [ms]")
        plt.ylabel("Voltage [mV]")
        #plt.axis(ymin=-120, ymax=-50)
        plt.show()
    """
    def set_recording(self):
        """
        Set a recording vector to the soma
        """
        # Record Time
        self.rec_t = h.Vector()
        self.rec_t.record(h._ref_t)
        # Record Voltage
        self.rec_v = h.Vector()
        self.rec_v.record(self.cell.soma(0.5)._ref_v)

    def get_recording(self):
        """
        Return the recordings of the somatic voltage and a time axis
        """
        time = np.array(self.rec_t)
        voltage = np.array(self.rec_v)
        return time, voltage

    def go(self, sim_time=None):
        """
        Launch a simulation of a given time

        Parameters
        ----------
        sim_time: an integer
            the time in millisecond of the simulation
            it replaces the self.sim_time if defined

        """

        h.t = 0
        #Start recording
        self.set_recording()
        #h.dt = self.dt
        #self.cell.initialise()
        init()
        #while h.t < self.sim_time: #I was using this procedure before do not know which one is better
        #    h.fadvance()
        if sim_time:
            run(sim_time)
        else:
            run(self.sim_time)

        
        
        
        
        
        


