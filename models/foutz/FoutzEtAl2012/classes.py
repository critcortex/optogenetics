from neuron import h
from functions import *

class Hu(object):
    def __init__(self, fiberD=2.0):
        """ Initialize Cell by loading hoc file """
        #-----------------------------------------------------------------------
        # Cell Parameters
        #-----------------------------------------------------------------------
        self.fiberD = fiberD
        #-----------------------------------------------------------------------
        # Optical properties of tissue
        #-----------------------------------------------------------------------
        # Scattering only Kubelka-Munk Model
        #scatter_coefficient      = 4.0 # fit from gradinaru et al
        # 10.3    # 1/mm #S = 10.3 # for rat # Aravanis et al JNE 2007 p s146
        # Kubelka-Munk General model:
        # *Ignoring* light spreading geometry in gradinaru data
        # Fit with matlab fit toolbox, nonlinear least square fit
        #
        # c(X) = (sqrt((1+K/S)^2-1)) / ((1+K/S) * sinh((sqrt((1+K/S)^2-1)) * S * X)
        #                 + (sqrt((1+K/S)^2-1)) * cosh((sqrt((1+K/S)^2-1)) *
        #                S * X))
        #
        # Coefficients (with 95% confidence bounds):
        h('absorbance_coefficient = 0.1249') # (1/mm) # Range: (0.05233, 0.1975)
        h('scatter_coefficient = 7.37')      # (1/mm) # Range: (6.679, 8.062)
        #-----------------------------------------------------------------------
        # Construct Cell
        #-----------------------------------------------------------------------
        self._construct_cell()
        #-----------------------------------------------------------------------
        # Create action potential counters
        #-----------------------------------------------------------------------
        self.apcs = []
        self._ap_counters()
        #-----------------------------------------------------------------------
        # Scaling
        #-----------------------------------------------------------------------
        self.current_scale=1
        #-----------------------------------------------------------------------
        # Visualization
        #-----------------------------------------------------------------------
        self.mlab_cell = None
        self._recordings = {}
    def __str__(self):
        return "HU"
    def _construct_cell(self):
        """ Layer V Pyramidal Cell Hu"""
        h.load_file("cell.hoc")
        self.root = h.soma
        h('number_of_apc = 1')
        h('required_aps = 1')
        h('axonnodes = 14')
    def _ap_counters(self):
        """ Create action potential counters, esp. useful for threshold
        calculation """
        self.apcs = []
        self.apc_times = []
        if h.number_of_apc > 0:
            try:
                sec = h.node[h.axonnodes.__int__()-2]
            except AttributeError:
                print "No node compartments!"
                return 0
            apc = h.APCount(0.5, sec=sec)
            apc.thresh = 0 # mV
            apc_time = h.Vector()
            apc.record(apc_time)
            self.apcs.append(apc)
            self.apc_times.append(apc_time)
        if h.number_of_apc == 2:
            sec = h.node[1]
            apc = h.APCount(0.5, sec=sec)
            apc.thresh = 0 # mV
            apc_time = h.Vector()
            apc.record(apc_time)
            self.apcs.append(apc)
            self.apc_times.append(apc_time)
        else:
            if h.number_of_apc>2:
                raise ValueError,"Too many apc counters; only 1 or 2 allowed"
    def intensity(self,sec):
        return sec.irradiance_chanrhod * sec.Tx_chanrhod
    def photon_flux(self, sec):
        """ Determine the light intensity at a given section (photons/ms cm2)
        """
        return h.photons_chanrhod * sec.Tx_chanrhod #  photons/ms cm2
    def photons(self,sec):
        section_area = h.area(0.5,sec=sec) # um2
        section_intensity = self.intensity(sec) #  photons/ms cm2
    def root_section(self):
        return h.SectionRef().root
    def build_tree(self, func,segfunc=False):
        """ func must act on a neuron section
        """
        from numpy import array
        print "-"*100
        def append_data(sec, xyzdv, parent_id, connections,func,segfunc):
            """ Append data to xyzdv
            """
            if not segfunc: v=func(sec)
            n = int(h.n3d(sec=sec))
            for ii in xrange(1, n):
                x = h.x3d(ii,sec=sec)
                y = h.y3d(ii,sec=sec)
                z = h.z3d(ii,sec=sec)
                d = h.diam3d(ii,sec=sec)
                if segfunc:
                    if n==1:v=func(sec(0.5))
                    else:v = func(sec(ii/float(n-1)))
                xyzdv.append([x,y,z,d,v])
                child_id = len(xyzdv)-1
                if len(xyzdv)>1:
                    connections.append([child_id, parent_id])
                parent_id = child_id
            return xyzdv, connections

        def append_children_data(parent, parent_id, xyzdv, connections, func, segfunc):
            sref = h.SectionRef(sec=parent)
            if sref.child:
                for child in sref.child:
                    xyzdv, connections = append_data(child, xyzdv, parent_id, connections, func, segfunc)
                    xyzdv, connections = append_children_data(parent = child,
                                                              parent_id = len(xyzdv)-1,
                                                              xyzdv = xyzdv,
                                                              connections = connections,
                                                              func = func,
                                                              segfunc = segfunc)
            return xyzdv, connections

        # Find data and connections
        root_section = self.root_section()
        if segfunc:
            if root_section.nseg==1:
                v = func(root_section(0.5))
            else:
                v = func(root_section(0.0))
        else:
            v=func(root_section)
        xyzdv = [[h.x3d(0,sec=root_section),h.y3d(0,sec=root_section),h.z3d(0,sec=root_section),h.diam3d(0,sec=root_section),v]]
        xyzdv, connections = append_data(root_section, xyzdv, 0, [],func,segfunc)
        xyzdv, connections = append_children_data(root_section,len(xyzdv)-1,xyzdv,connections,func,segfunc)
        self.xyzdv = array(xyzdv)
        self.connections = array(connections)
    def move(self, xyz, move_mlab=False):
        """ Move visualization and cell """
        from neuron import h
        if move_mlab:
            if self.mlab_cell:
                self.mlab_cell.mlab_source.x = self.mlab_cell.mlab_source.x + xyz[0]
                self.mlab_cell.mlab_source.y = self.mlab_cell.mlab_source.y + xyz[1]
                self.mlab_cell.mlab_source.z = self.mlab_cell.mlab_source.z + xyz[2]
        tree = h.SectionList()
        tree.wholetree(sec=self.root)
        for sec in tree:
            for ii in xrange(h.n3d(sec=sec).__int__()):
                x=h.x3d(ii,sec=sec)
                y=h.y3d(ii,sec=sec)
                z=h.z3d(ii,sec=sec)
                d=h.diam3d(ii,sec=sec)
                h.pt3dchange(ii,x+float(xyz[0]),y+float(xyz[1]),z+float(xyz[2]),d)
    def retrieve_coordinates(self, sec):
        xyzds = []
        for ii in xrange(int(h.n3d(sec=sec))):
            xyzds.append([h.x3d(ii,sec=sec),
                          h.y3d(ii,sec=sec),
                          h.z3d(ii,sec=sec),
                          h.diam3d(ii,sec=sec)])
        return xyzds
    def display(self, func, segfunc=False, scaling=1, replace=True, clim=None, colormap='jet'):
        ''' Display current cell in mayavi
        '''
        #from neuron import h
        from numpy import array, vstack
        try:
            from enthought.mayavi import mlab
            from enthought.mayavi.mlab import pipeline
        except:
            from mayavi import mlab
            from mayavi.mlab import pipeline
        if replace:
            try:self.mlab_cell.parent.parent.parent.parent.parent.parent.remove()
            except AttributeError:pass
        ### Turn off vtk warnings # # # # # # # # # # # # # # # # # # # # # # #
        from vtk import vtkObject
        o = vtkObject
        o.GetGlobalWarningDisplay()
        o.SetGlobalWarningDisplay(0) # Turn it off.

        self.build_tree(func, segfunc)
        xs = self.xyzdv[:,0]
        ys = self.xyzdv[:,1]
        zs = self.xyzdv[:,2]

        # don't want scaling for soma segments
        diams = self.xyzdv[:,3]
        nonsoma = (diams < 15) # non-somatic
        diams += diams*nonsoma*(scaling-1)
        #diams = self.xyzdv[:,3] * scaling # larger scaling makes neurons more visible
        data = self.xyzdv[:,4]
        edges = self.connections

        # Display in mayavi
        pts = pipeline.scalar_scatter(xs, ys, zs, diams/2.0,
                                      name=str(self))
        dataset = pts.mlab_source.dataset
        dataset.point_data.get_array(0).name = 'diameter'
        dataset.lines = vstack(edges)

        array_id = dataset.point_data.add_array(data.T.ravel())
        dataset.point_data.get_array(array_id).name = 'data'
        dataset.point_data.update()

        #### Create tube with diameter data
        src = pipeline.set_active_attribute(pts,
                                            point_scalars='diameter')
        stripper = pipeline.stripper(src)
        tube = pipeline.tube(stripper,
                             tube_sides = 8,
                             tube_radius = 1)
        tube.filter.capping = True
        tube.filter.use_default_normal = False
        tube.filter.vary_radius = 'vary_radius_by_absolute_scalar'
        #tube.filter.radius_factor = 90.0 # just for making movies
        src2 = pipeline.set_active_attribute(tube, point_scalars='data')

        lines = pipeline.surface(src2,colormap = colormap)
        if clim:
            from numpy import array
            lines.parent.scalar_lut_manager.use_default_range = False
            lines.parent.scalar_lut_manager.data_range = array(clim)
        self.mlab_cell = lines
    def plot(self, func, scaling = 1, segfunc=False, clim=None,cmap=None):
        """ plot cell in matplotlib line plot collection
        """
        from numpy import array, linspace
        from matplotlib.collections import LineCollection
        from matplotlib import pyplot
        self.build_tree(func,segfunc)
        pts   = self.xyzdv[:,:2]
        edges = self.connections
        diam  = self.xyzdv[:,3]
        data  = self.xyzdv[:,4]
        print "DATA RANGE: ",data.min(),data.max()
        # Define colors
        if not cmap:
            from matplotlib.cm import jet as cmap
        if not clim:
            clim=[data.min(),data.max()]
        a = (data - clim[0])/(clim[1]-clim[0])
        # Define line segments
        segments = []
        for edge in edges:
            segments.append([pts[edge[0],:], pts[edge[1],:]])
        # Build Line Collection
        collection = LineCollection(segments = array(segments),
                                    linewidths = diam*scaling,
                                    colors=cmap(a))
        collection.set_array(data)
        collection.set_clim(clim[0], clim[1])
        pyplot.gca().add_collection(collection,autolim=True)
        pyplot.axis('equal')
        return collection
    def channels_in_list(self,seclist):
        channels = 0
        for sec in seclist:
            if h.ismembrane('chanrhod',sec = sec):
                for seg in sec:
                    rho  = seg.channel_density_chanrhod/1e8 # 1/cm2 --> 1/um2
                    area = h.area(seg.x, sec=sec) # um2
                    n = rho * area
                    channels += n
        return channels
    def area_in_list(self,seclist):
        area = 0
        for sec in seclist:
            if h.ismembrane('chanrhod',sec = sec):
                for seg in sec:
                    area += h.area(seg.x, sec=sec) # um2
        return area
    def illuminated_area_in_list(self,seclist,Tx_threshold = 0.001):
        area = 0
        for sec in seclist:
            if h.ismembrane('chanrhod',sec = sec):
                for seg in sec:
                    if seg.Tx_chanrhod>Tx_threshold:
                        area += h.area(seg.x, sec=sec) # um2
        return area
    def open_channels_in_list(self,seclist):
        open_channels = 0
        for sec in seclist:
            if h.ismembrane('chanrhod', sec = sec):
                for seg in sec:
                    rho  = seg.channel_density_chanrhod/1e8 # 1/cm2 --> 1/um2
                    area = h.area(seg.x, sec=sec) # um2
                    try:
                        f_open = seg.o2_chanrhod + seg.o1_chanrhod # open fraction # 4 state model
                    except:
                        f_open = seg.o1_chanrhod # open fraction # 3 state model
                    n = f_open * rho * area
                    open_channels += n
        return open_channels
    def get_axonal(self):
        """ Additional iseg compartment
        """
        secs = [h.hill]
        secs.extend([sec for sec in h.ais])
        secs.extend([h.nakeaxon])
        secs.extend([sec for sec in h.myelin])
        secs.extend([sec for sec in h.node])
        return secs
    def get_axonal_channels(self):
        return self.channels_in_list(self.axonal)
    def get_open_axonal_channels(self):
        return self.open_channels_in_list(self.axonal)
    def get_dendritic(self):
        secs = [sec for sec in h.somatodendritic]
        #secs.extend([h.hill])
        return secs
    def get_dendritic_channels(self):
        return self.channels_in_list(self.dendritic)
    def get_open_dendritic_channels(self):
        return self.open_channels_in_list(self.dendritic)
    def get_apical_tuft(self):
        """ Return a list of all sections which make up the apical tuft, starting
        at the branch point
        """
        secs=[]
        for ii in xrange(23,len(h.dend11)):
            secs.append(h.dend11[ii])
        return secs
    def get_apical_shaft(self):
        """ Return the sections which compose the apical shaft
        """
        secs=[]
        for ii in [0,4,10,16,18,20,22]:
            secs.append(h.dend11[ii])
        return secs
    def get_basilar_tuft(self):
        """ Return the dendritic sections which compose the basilar tuft
        """
        secs=[]
        for ii,dendrite in enumerate((h.dend1,h.dend2,h.dend3,h.dend4,h.dend5,h.dend6,
                                      h.dend7,h.dend8,h.dend9,h.dend10,h.dend11)):
            if ii==10:
                for jj,sec in enumerate(dendrite):
                    if jj < 23: # Apical Tuft
                        if jj not in [0,4,10,16,18,20,22]: # Apical Shaft
                            secs.append(sec)
            else:
                for sec in dendrite:
                    secs.append(sec)
        return secs
    def get_somatic(self):
        """ Return the sections which compose the Soma
        """
        return [h.soma]
    def set_density_distribution(self, distribution=0.5, n_channels = 1e7):
        """ Set density in dendritic compartments
        distribution: 0.0 - Higher Somatic density
                      0.5 - Uniform distribution
                      1.0 - Higher Apical density
        """
        # Find the maximal distance between the soma and all dendrities
        max_distance = 0
        for sec in self.dendritic:
            for seg in sec:
                max_distance = max([max_distance, self.seg_section_distance(seg)])

        for sec in self.dendritic:
            for seg in sec:
                distance=self.seg_section_distance(seg)
                s = (max_distance-distance)/(max_distance) # soma centric weighting
                a = (distance)/(max_distance)              # apical centric weighting
                W = distribution
                seg.channel_density_chanrhod  =  s*(1-W) + a*W
        scale = n_channels/self.dendritic_channels
        for sec in self.dendritic:
            for seg in sec:
                seg.channel_density_chanrhod = scale*seg.channel_density_chanrhod
        assert 0.001 > (n_channels - self.dendritic_channels)
    def get_open_channels(self):
        return self.open_channels_in_list(h.allsec())
    def get_total_channels(self):
        return self.channels_in_list(h.allsec())
    def get_icat(self):
        """ Determine the total amount of channelrhodopsin current in the cell
        """
        icat = 0
        for sec in h.allsec():
            if h.ismembrane('chanrhod',
                            sec = sec):
                for seg in sec:
                    i = seg.icat_chanrhod # (mA/cm2)
                    area = h.area(seg.x, sec=sec)/1e8      # cm2
                    icat += area * i # mA
        return icat
    def set_required_aps(self,stimulator,additional_aps=1):
        """ Determine the number of action potentials normally occuring, so that
        we can set a goal number of additional aps
        Also, we check to make sure that the stimulator's amplitude is set to 0
        """
        from neuron import h
        initial_amplitude=stimulator.amplitude
        stimulator.amplitude=0
        h.run()
        h.required_aps=0
        assert self.apc_times,'No action potential counters'
        for apct in self.apc_times:
            h.required_aps=max((h.required_aps,len(apct)))
        print "** NO STIM APS: %d; Goal APS: %d **" % (h.required_aps,
                                                       h.required_aps+additional_aps)
        h.required_aps += additional_aps # usually require one additional ap
        stimulator.amplitude=initial_amplitude
    def set_tstop(self,tstop,stimulator,additional_aps=1):
        h.tstop=tstop
        self.set_required_aps(stimulator,additional_aps)
    def get_response(self):
        """ Determine if an action potential has occurred
        """
        # Determine if any of the counters saw less than the requisite number of action potentials
        for apct in self.apc_times:
            if len(apct) < h.required_aps:
                return False
        # If all counters saw the requisite number, than return True
        return True
    def seg_section_distance(self,seg,root_section=None):
        """ Returns the distance between each segment of section, and the
        root_section
        """
        if not root_section:root_section=self.root
        h.distance(0, root_section(0.5).x, sec=root_section)
        return h.distance(seg.x, sec=seg.sec)
    response = property(get_response)
    open_channels = property(get_open_channels)
    total_channels = property(get_total_channels)
    icat = property(get_icat)
    apical_tuft  = property(get_apical_tuft)
    apical_shaft = property(get_apical_shaft)
    basilar_tuft = property(get_basilar_tuft)
    dendritic = property(get_dendritic)
    dendritic_channels = property(get_dendritic_channels)
    open_dendritic_channels = property(get_open_dendritic_channels)
    somatic = property(get_somatic)
    axonal = property(get_axonal)
    axonal_channels = property(get_axonal_channels)
    open_axonal_channels = property(get_open_axonal_channels)
class Optrode(object):
    """ xyz0 is the base, and xyz1 is the tip """
    def __init__(self,origin,delay=1,duration=5,initial_amplitude=38.0,distance=1000,pulses=1,frequency=1):
        from numpy import pi
        self.origin=h.secname(sec=origin)
        self.sec=h.Section(name=str(self))
        self.sec.L=1000
        self.sec.diam=200 # um # Aravanis: 200 um # Gradinaru: 400 um
        self.stim=h.ostim(0.5,sec=self.sec)
        self.delay=delay
        self.pulses=pulses
        self.frequency=frequency
        self.duration=duration
        self.amplitude=initial_amplitude
        h.setpointer(h._ref_source_irradiance_chanrhod, 'irradiance',self.stim)
        h.setpointer(h._ref_source_photons_chanrhod,   'photons',self.stim)
        h.setpointer(h._ref_source_flux_chanrhod,      'flux',self.stim)
        h.setpointer(h._ref_tstimon_chanrhod,      'tstimon',self.stim)
        h.setpointer(h._ref_tstimoff_chanrhod,      'tstimoff',self.stim)
        self.stim.radius=self.sec.diam/2.0
        self.stim.pulses=self.pulses
        self.stim.isi = 1000 / self.frequency - self.duration #in ms
        self.stim.amp=initial_amplitude
        self.absorbance_coefficient = 0.1249 # (1/mm) # Range: (0.05233, 0.1975)
        self.scatter_coefficient = 7.37      # (1/mm) # Range: (6.679, 8.062)
        self.n = 1.36        # index of refraction of gray matter
        self.NA = 0.37        # numerical aperture of the optical fiber
        #self.NA = 0.48
        self.set_distance(origin, distance)
    def __str__(self):
        return "Fiber Optic"
    def __info__(self):
        info=str(self)
        info+="\n-Section: %s" % self.origin
        xyz0,xyz1=self.xyz
        info+="\n-Base: %g,%g,%g" % (xyz0[0],xyz0[1],xyz0[2])
        info+="\n-Tip: %g,%g,%g" % (xyz1[0],xyz1[1],xyz1[2])
        info+="\n-Delay: %s" % self.delay
        info+="\n-Duration: %s" % self.duration
        info+="\n-Amplitude: %s" % self.amplitude
        info+="\n-Length: %g" % self.length
        info+="\n-Diameter: %g" % self.diameter
        info+="\n-closest_section: %s" % h.secname(sec=self.closest_section)
        return info
    def _find_axon_trajectory(self, center):
        """ Find the normalized vector which describes axon trajectory """
        from numpy import sqrt
        most_distant_node = find_mean_section_coordinates(h.node[-2])
        trajectory = most_distant_node - center
        trajectory /= sqrt(sum(trajectory**2))  # Normalize
        return trajectory
    def interpxyz(self):
        """ interpolated data, spaced at regular intervals
        """
                # First, need to interpolate centers unto all compartments; from interpxyz.hoc
        for sec in h.allsec():
            #if h.ismembrane('chanrhod',sec=sec):
            if h.ismembrane('chanrhod',sec=sec):
                nn = h.n3d(sec=sec).__int__()
                xx = h.Vector(nn)
                yy = h.Vector(nn)
                zz = h.Vector(nn)
                length = h.Vector(nn)
                for ii in xrange(nn):
                    xx.x[ii] = h.x3d(ii,sec=sec)
                    yy.x[ii] = h.y3d(ii,sec=sec)
                    zz.x[ii] = h.z3d(ii,sec=sec)
                    length.x[ii] = h.arc3d(ii,sec=sec)
                # to use Vector class's .interpolate() must first scale the
                # independent variable i.e. normalize length along centroid
                length.div(length.x[nn-1])
                # initialize the destination "independent" vector
                rr = h.Vector(sec.nseg+2)
                rr.indgen(1./sec.nseg)
                rr.sub(1./(2.*sec.nseg))
                rr.x[0]=0.
                rr.x[sec.nseg+1]=1.
                # length contains the normalized distances of the pt3d points
                # along the centroid of the section.  These are spaced at
                # irregular intervals.
                # range contains the normalized distances of the nodes along the
                # centroid of the section.  These are spaced at regular intervals.
                # Ready to interpolate.
                xint = h.Vector(sec.nseg+2)
                yint = h.Vector(sec.nseg+2)
                zint = h.Vector(sec.nseg+2)
                xint.interpolate(rr, length, xx)
                yint.interpolate(rr, length, yy)
                zint.interpolate(rr, length, zz)
                # for each node, assign the xyz values to x_xtra, y_xtra, z_xtra
                # don't bother computing coords of the 0 and 1 ends
                # also avoid writing coords of the 1 end into the last internal node's coords
                for ii in range(1,sec.nseg+1):
                    xr = rr.x[ii]
                    #sec(xr).x_chanrhod = xint.x[ii]
                    #sec(xr).y_chanrhod = yint.x[ii]
                    #sec(xr).z_chanrhod = zint.x[ii]
                    sec(xr).x_chanrhod = xint.x[ii]
                    sec(xr).y_chanrhod = yint.x[ii]
                    sec(xr).z_chanrhod = zint.x[ii]
    def find_illumination(self,X,Y,Z,spreading=True,scattering=True):
        from numpy import sqrt
        def gaussian(r,radius):
            """ r is displacement from center
            95.4 % of light is within the radius (2 standard deviations)
            constant energy in distribution
            """
            from numpy import array, pi,sqrt, exp
            r = 2*array(r)/array(radius)
            dist = (1/sqrt(2*pi)) * exp((r**2)/(-2))
            return dist/0.4
        def kubelka_munk(distance):
            """
            distance to center of optrode, approximates mean distance to all points along surface of optrode
            distance in um
            """
            from numpy import sqrt,sinh,cosh
            #K = 0.1248e-3 # 1/um
            #S = 7.37e-3   # 1/um
            K = self.absorbance_coefficient * 1e-3 # (1/um) # Range: (0.05233, 0.1975)
            S = self.scatter_coefficient * 1e-3    # (1/um) # Range: (6.679, 8.062)
            a = 1 + K / S                # unitless
            b = sqrt(a ** 2 - 1)         # unitless
            Tx = b / (a * sinh(b * S * distance) + b * cosh(b * S * distance)) # distance in um - losses due to absorption and scattering through the tissue on top of losses due to beam quality?
            Tx[distance<0]=0 # negative values set to zero
            return Tx
        def apparent_radius(z,radius):
            """ Find the apparent radius at a distance z
            """
            from numpy import tan
            return radius + z*tan(self.theta_div)
        def spread(z):
            """ irradiance loss due to spreading
            """
            from numpy import sqrt,pi
            rho = self.radius * sqrt(((self.n/self.NA)**2) - 1)
            return rho**2 / ((z + rho)**2)
        r,z = find_cylindrical_coords(X,Y,Z,self.xyz)
        if scattering: # kubelka-munk scattering
            Kx = kubelka_munk(sqrt(r**2+z**2))
        else:
            Kx = 1
        if spreading: # conservation of energy spreading
            Sx = spread(z)
            radius = apparent_radius(z,self.radius)
        else:
            Sx = 1
            radius = self.radius
        return Sx * Kx * gaussian(r,radius)
    def calc_tx(self):
        """ Set the fractional illumionation for all sections which have chanrhod
        density mechanisms """
        from numpy import sqrt, sinh, cosh, arccos, tan, array, dot, isnan, pi
        self.interpxyz()

        # !All units should be in um
        #---------------------------------------------------------------------------
        # Location of each segment
        #---------------------------------------------------------------------------
        seg_xyz = []
        for sec in h.allsec():
            #if h.ismembrane('chanrhod', sec=sec):
            if h.ismembrane('chanrhod', sec=sec):
                for seg in sec:
                    xyz = seg.x_chanrhod, seg.y_chanrhod, seg.z_chanrhod
                    seg_xyz.append(xyz)
        seg_xyz = array(seg_xyz)
        Tx = self.find_illumination(seg_xyz[:,0],seg_xyz[:,1],seg_xyz[:,2])
        #---------------------------------------------------------------------------
        # Set Tx_chanrhod
        #---------------------------------------------------------------------------
        ii = 0
        #for sec in h.allsec():
            #if h.ismembrane('chanrhod', sec=sec):
                #for seg in sec:
                    #seg.Tx_chanrhod = Tx[ii]
                    #ii+=1
        for sec in h.allsec():
            if h.ismembrane('chanrhod', sec=sec):
                for seg in sec:
                    seg.Tx_chanrhod = Tx[ii]
                    ii+=1
    def rotate_optrode(self, angle, distance, sec, axis_defining_plane='z'):
        """ rotate optrode around a given sec,
        to a given angle,
        at a given distance
        in the plane defined by the axon and the given vectro
        """
        from numpy import cross, sqrt, array

        # Determine center of rotation
        center = find_mean_section_coordinates(sec)

        # Determine starting position
        xyz0 = self._find_axon_trajectory(center) * (distance+1000)  # end optrode along axon
        xyz1 = self._find_axon_trajectory(center) * distance         # start optrode along axon

        # Find vector perpendicular to axon and given axis
        if axis_defining_plane=='x':
            uvw = cross(xyz1, array([1,0,0]))
        elif axis_defining_plane=='y':
            uvw = cross(xyz1, array([0,1,0]))
        elif axis_defining_plane=='z':
            uvw = cross(xyz1, array([0,0,1]))
        else:
            raise ValueError,'No such plane: %s' % axis_defining_plane
        uvw /= sqrt(sum(uvw**2))  # Normalize

        # Rotate the optrode around the vector
        new_xyz0 = rotate_point_around_vector(xyz0,
                                              uvw,
                                              angle)
        new_xyz1 = rotate_point_around_vector(xyz1,
                                              uvw,
                                              angle)

        # Redefine optrode position
        self.set_position([float(new_xyz0[0]),float(new_xyz1[0])],
                         [float(new_xyz0[1]),float(new_xyz1[1])],
                         [float(new_xyz0[2]),float(new_xyz1[2])])
    def draw(self):
        """ Draw a visual representation of the optrode """
        from neuron import gui
        if not hasattr(self,'gOpt'):
            from numpy import pi
            self.gOpt=h.Shape(0)
            self.gOpt.view(-2100, -2100, 4200, 4200, 230, 450, 200.64, 200.32)
            self.gOpt.rotate(0, 0, 0, pi/2, 0, 0)
        h.pt3dclear(sec = self.sec)
        h.pt3dadd(self.x[0],
                  self.y[0],
                  self.z[0],
                  self.diameter,
                  sec = self.sec)
        h.pt3dadd(self.x[1],
                  self.y[1],
                  self.z[1],
                  self.diameter,
                  sec = self.sec)
        self.pOpt0 = h.IClamp(0,
                              sec = self.sec)
        self.pOpt1 = h.IClamp(1,
                              sec = self.sec)

        self.gOpt.point_mark(self.pOpt0, 1)  # make start black
        self.gOpt.point_mark(self.pOpt1, 3)  # make output blue
        self.gOpt.exec_menu("Show Diam")
        self.gOpt.exec_menu("3D Rotate")
    def display(self,show_illum=True,dvol=50,cladding=None,scattering=True,spreading=True,bounds = None):
        """ Display optrode in mayavi
        """
        if not bounds:
            if hasattr(self,"bounds"):
                bounds = self.bounds
            else:
                bounds = [[-500,500],[-500,500],[-500,500]]
        self.bounds = bounds
        if cladding==None:
            if hasattr(self,'mlab_cladding'):
                cladding = True
            else:
                cladding = False
        if hasattr(self,'mlab_tube'):
            self.mlab_tube.parent.parent.remove()
        if hasattr(self,'mlab_illum'):
            self.mlab_illum.parent.parent.remove()
        if hasattr(self,'mlab_cladding'):
            self.mlab_cladding.parent.parent.remove()
        try:
            from enthought.mayavi import mlab
        except:
            from mayavi import mlab
        cyl = mlab.plot3d(self.x,
                          self.y,
                          self.z,
                          name='optrode',
                          color=(.9,.9,.9))
        self.mlab_tube = cyl.parent.parent
        self.mlab_tube.filter.capping = True
        self.mlab_tube.filter.number_of_sides = 20
        self.mlab_tube.filter.radius = self.radius
        if cladding:
            from numpy import array
            clad= mlab.plot3d(self.x,
                              self.y,
                              self.z-array([0.1,0.1]),
                              name='cladding',
                              color=(0.5,0.5,0.5),
                              opacity = 0.5)
            self.mlab_cladding = clad.parent.parent
            self.mlab_cladding.filter.capping = True
            self.mlab_cladding.filter.number_of_sides = 20
            self.mlab_cladding.filter.radius = self.radius*2
            self.mlab_cladding.children[0].children[0].actor.property.backface_culling = True
        if show_illum:
            from numpy import diff, mgrid,array,matrix,c_,cos,sin,arctan,ones
            from numpy.linalg import norm
            x = self.xyz[1,0]
            y = self.xyz[1,1]
            z = self.xyz[1,2]
            X,Y,Z = mgrid[x+bounds[0][0]:x+bounds[0][1]:dvol*1j,
                          y+bounds[1][0]:y+bounds[1][1]:dvol*1j,
                          z+bounds[2][0]:z+bounds[2][1]:dvol*1j]

            Tx = self.find_illumination(X,Y,Z,spreading,scattering)
            self.mlab_illum = mlab.contour3d(X,Y,Z,Tx,
                                             #opacity=0.8,
                                             transparent=True,
                                             vmin=0.001,
                                             vmax=0.1,
                                             contours=[t for t in [0.1,0.01,0.001]])
            self.mlab_illum.parent.scalar_lut_manager.use_default_range = False
            self.mlab_illum.parent.scalar_lut_manager.data_range = array([ 0.001,  0.1   ])
            self.mlab_illum.parent.scalar_lut_manager.lut.scale='log10'
            #self.mlab_illum.actor.property.backface_culling = True
            self.mlab_illum.actor.property.frontface_culling = True
    def record(self,args=('i',)):
        self._recording={'t':h.Vector()}
        self._recording['t'].record(h._ref_t)
        for k in args:
            self._recording[k]=h.Vector()
            exec('self._recording[k].record(self.stim._ref_%s)' % k)
    def plot(self,show=False):
        from matplotlib.pyplot import plot,legend
        recordings=self.recordings
        t=recordings['t']
        for k,v in recordings.items():
            if k != 't':
                plot(t,v,label=k)
        legend()
        if show:
            from matplotlib.pyplot import show as shw
            shw()
    def get_duration(self):
        return self.stim.dur
    def set_duration(self,duration):
        self.stim.dur=duration
    def get_amplitude(self):
        return self.stim.amp
    def set_amplitude(self,amplitude):
        self.stim.amp=amplitude
    def get_delay(self):
        return self.stim.delay
    def set_delay(self,delay):
        self.stim.delay=delay
    def get_recordings(self):
        from numpy import array
        recordings={}
        for k in self._recording.keys():
            recordings[k]=array(self._recording[k].to_python())
        return recordings
    def set_position(self, x, y, z):
        """ Move electrode to new coordinates:
            x0, y0, z0: optrode input
            x1, x1, z1: optrode output
        """
        from numpy.linalg import norm
        self._x=x
        self._y=y
        self._z=z
        xyz0,xyz1=self.xyz
        #self.sec.L=dist(xyz0,xyz1)
        self.sec.L = norm(xyz1-xyz0)
        h.pt3dclear(sec=self.sec)
        h.pt3dadd(float(x[0]), float(y[0]), float(z[0]),
                  self.radius * 2, sec=self.sec)
        h.pt3dadd(float(x[1]), float(y[1]), float(z[1]),
                  self.radius * 2, sec=self.sec)
        self.calc_tx()
        if hasattr(self, "gOpt"): self.draw()
        if hasattr(self, "mlab_tube"): self.display()
    def set_longitudinal_radial(self, longitudinal_percent, radial_distance, center, terminal_node=None):
        from numpy import array, sign
        from numpy.linalg import norm

        if not terminal_node:
            try:
                terminal_node = h.node[-1]
            except:
                raise StandardError('No h.node[-1] compartment!')

        # Place electrode along length of neuron
        fraction_length = longitudinal_percent / 100.0
        xyz_soma = find_mean_section_coordinates(center)
        xyz_terminal_node = find_mean_section_coordinates(terminal_node)
        xyz_optrode = xyz_soma + \
                    fraction_length * (xyz_terminal_node - xyz_soma) + \
                    array([0, 0, -1 * radial_distance])

        # Find displacement from underneath center (for analysis)
        # essentially, when plotting, it is more useful to demonstrate the
        # displacement in mm rather than the  fractional longitudinal distance
        longitudinal_displacement = sign(longitudinal_percent) * norm(fraction_length * (xyz_terminal_node - xyz_soma))

        # Set optrode at location
        self.set_position([xyz_optrode[0], xyz_optrode[0]],
                         [xyz_optrode[1], xyz_optrode[1]],
                         [xyz_optrode[2] - self.length, xyz_optrode[2]])
        return xyz_optrode, longitudinal_displacement
    def get_distance(self, sec):
        """ Helper method to find the distance between the optrode tip and a
        given compartment (section)
        """
        from numpy.linalg import norm
        from numpy import array
        from neuron import nrn
        optrode_output = self.xyz[1]
        if isinstance(sec,nrn.Section):
            return norm(find_mean_section_coordinates(sec) -
                        optrode_output)
        elif isinstance(sec,nrn.Segment):
            seg = sec
            print find_section_coordinates(seg.sec),seg.x
            raise StandardError,"Not yet implemented"
        else:
            raise TypeError, "Wrong type: %s" % type(sec)
    def set_distance(self, sec, z_distance):
        """ Set optrode a certain distance in z-direction below the given section
        directed upwards """
        x, y, z = find_mean_section_coordinates(sec)
        self.set_position([x,x],
                          [y,y],
                          [z - z_distance - self.length,
                           z - z_distance])
        assert (z_distance == self.get_distance(sec))
    def get_x(self):
        return self._x
    def set_x(self, val):
        from numpy import array
        assert len(val)==2
        self._x=array([float(val[0]),
                       float(val[1])])
        from numpy.linalg import norm
        xyz0,xyz1=self.xyz
        self.length = norm(xyz1-xyz0)
    def get_y(self):
        return self._y
    def set_y(self, val):
        from numpy import array
        assert len(val)==2
        self._y=array([float(val[0]),
                       float(val[1])])
        from numpy.linalg import norm
        xyz0,xyz1=self.xyz
        self.length = norm(xyz1-xyz0)
    def get_z(self):
        return self._z
    def set_z(self, val):
        from numpy import array
        assert len(val)==2
        self._z=array([float(val[0]),
                       float(val[1])])
        from numpy.linalg import norm
        xyz0,xyz1=self.xyz
        self.length = norm(xyz1-xyz0)
    def get_xyz(self):
        from numpy import array
        x=self.x
        y=self.y
        z=self.z
        return array([[x[0],y[0],z[0]],[x[1],y[1],z[1]]])
    def set_xyz(self,val):
        from numpy import array
        val = array(val).astype(float)
        #print "SHAPE: ",val.shape
        if val.shape == (3,2):
            self.set_position(val[0,:],val[1,:],val[2,:])
        elif val.shape == (2,3):
            self.set_position(val[:,0],val[:,1],val[:,2])
        else:
            raise ValueError, "Trying to set xyz with inappropriate array %s" % str(val)
    def get_length(self):
        return self.sec.L
    def set_length(self,val):
        self.sec.L = val
        x0, x1 = self.x
        y0, y1 = self.y
        z0, z1 = self.z
        new_x0, new_y0, new_z0 = point_along_vector([x1, y1, z1],
                                                    [x0-x1, y0-y1, z0-z1],
                                                    val)
        self.set_position([new_x0, x1],
                         [new_y0, y1],
                         [new_z0, z1])
        from numpy import array
        from numpy.linalg import norm
        xyz0,xyz1=self.xyz
        dist = norm(xyz1-xyz0)
        assert(approx_equal(self.length, dist))
    def get_diameter(self):
        """ Diameter in um
        """
        assert self.sec.diam==self.stim.radius*2
        return self.sec.diam
    def set_diameter(self,diameter):
        """ Diameter in um
        """
        self.radius=diameter/2.0
    def get_radius(self):
        return self.stim.radius
    def set_radius(self,radius):
        self.sec.diam=radius * 2.0
        self.stim.radius=radius
    def get_closestsection(self):
        """ Find the closest chanrhod+ section """
        closest_sec=None
        closest_sec_distance=1e9
        for sec in h.allsec():
            #if h.ismembrane('chanrhod', sec = sec):
            if h.ismembrane('chanrhod', sec = sec):
                sec_distance = self.get_distance(sec)
                if sec_distance<closest_sec_distance:
                    closest_sec=sec
                    closest_sec_distance=sec_distance
        return closest_sec
    def get_intensity(self):
        return self.stim.intensity
    def radiant_power(self,seclist):
        """ Find the radiant power integrated over all sections in seclist
        """
        p = 0
        for sec in seclist:
            for seg in sec:
                area = h.area(seg.x, sec=sec)*1e-8 # um2 --> cm2
                #p += (area * self.amplitude * sec.Tx_chanrhod) # Add section's radiant power
                p += (area * self.amplitude * sec.Tx_chanrhod)
        return p
    def get_sec(self):
        return self._sec
    def set_sec(self,sec):
        self._sec = sec
    def get_theta_div(self):
        from numpy import arcsin
        return arcsin(self.NA/self.n)
    duration=property(get_duration,set_duration)
    amplitude=property(get_amplitude,set_amplitude)
    delay=property(get_delay,set_delay)
    recordings=property(get_recordings)
    theta_div = property(get_theta_div)
    sec = property(get_sec,set_sec)
    x = property(get_x, set_x)
    y = property(get_y, set_y)
    z = property(get_z, set_z)
    xyz = property(get_xyz,set_xyz)
    length = property(get_length, set_length)
    closest_section = property(get_closestsection)
    diameter=property(get_diameter,set_diameter)
    radius=property(get_radius,set_radius)
    intensity=property(get_intensity)
    info=property(__info__)
class Sim(object):
    """ Serial simulation object
    setup up environment with the setup_func, and run a trial function with each
    parameter in a parameter set
    """
    def __init__(self,cell, stimulator, output_filename):
        self.cell = cell
        self.stimulator = stimulator
        self.output_filename=output_filename
        open(output_filename,'w') # Erase output file
    def distance_threshold(self,param_set):
        """ Find the threshold at various distances
        """
        self.data=[]
        for param in param_set:
            distance=param['Distance (um)']
            self.stimulator.diameter=param['Fiber Optic Diameter (mm)']*1e3 # mm --> um
            self.stimulator.set_distance(self.cell.root,
                                               distance)
            param['Threshold (W/cm2)']=self.find_threshold()
            self.data.append(param)
    def main(self,param_set):
        h.tstop = 60
        self.distance_threshold(param_set)
        self.flush()
    def flush(self):
        # Write data results
        f=open(self.output_filename,'a')
        keys=self.data[0].keys()
        f.write(','.join(keys)+'\n')
        for data in self.data:
            vals=[]
            for key in keys:
                vals.append(str(data[key]))
            f.write(','.join(vals)+'\n')
        f.close()
    def find_threshold(self,upper_limit=1e4,error_threshold=1e-3,verbose=True,additional_aps=1):
        from numpy import log10
        if not self.stimulator.amplitude:
            self.stimulator.amplitude=1
        self.cell.set_tstop(h.tstop,self.stimulator,additional_aps)
        supra=upper_limit/0.9
        sub=0
        if verbose:print "+ _SubT____, Amplitude, _SupraT__, _Error___ +"
        while (supra-sub)/self.stimulator.amplitude > error_threshold: # while error larger than threshold
            h.run()
            if self.cell.response:response='+'
            else:response='-'
            if verbose:print "%s %0.3e, %0.3e, %0.3e, %0.3e : %s" % (response,
                                                                     sub,
                                                                     self.stimulator.amplitude,
                                                                     supra,
                                                                     (supra-sub)/self.stimulator.amplitude,
                                                                     (4+int(log10((supra-sub)/self.stimulator.amplitude)))*'*')
            if self.cell.response:
                supra=self.stimulator.amplitude
                self.stimulator.amplitude = (supra+sub)/2.0
            else:
                if self.stimulator.amplitude >= upper_limit:
                    # Probably not going to reach upper
                    if verbose:print "** Upper threshold reached **"
                    self.stimulator.amplitude=0
                    h.run()
                    return float('nan')
                sub=self.stimulator.amplitude
                self.stimulator.amplitude = min([(supra+self.stimulator.amplitude)/2.0,self.stimulator.amplitude*2,upper_limit])
        h.run()
        if not self.cell.response:
            self.stimulator.amplitude=supra
            h.run()
        if verbose:print "** THRESHOLD: %g **" % self.stimulator.amplitude
        return self.stimulator.amplitude
class Data(object):
    def __init__(self, filename, seperator=','):
        from numpy import array
        data={}
        for line in open(filename,'r'):
            data_list = [t.strip() for t in line.strip().split(seperator)]
            line_type = self.data_or_header(data_list)
            if line_type =='data':
                for ii,x in enumerate(data_list):
                    if not data.has_key(ii):data[ii]=[]
                    try:
                        data[ii].append(float(x))
                    except:
                        data[ii].append(x.strip())
            elif line_type == 'header':
                headers=data_list
        if not data:
            raise ValueError,"Empty data file, or not seperated by %s" % seperator
        for k in data.keys():
            data[k] = array(data[k])
        self.filename=filename
        self.data={}
        for ii,header in enumerate(headers):
            d = data[ii]
            try:
                f = d.astype(float)
                i = f.astype(int)
                if (f==i).all():
                    self.data[header]=i
                else:
                    self.data[header]=f
            except ValueError:
                self.data[header]=d
    def data_or_header(self,data_list):
        for x in data_list:
            try:
                float(x)
                return 'data'
            except ValueError:
                pass
        return 'header'
    def sort(self, cname):
        """ Sort data by column with header name <<cname>>
        """
        from numpy import array, argsort
        sorted_indices = argsort(self.data[cname])
        for k in self.data.keys():
            self.data[k] = array(self.data[k])[sorted_indices]
    def set_slice(self, mask):
        self.slice={}
        for k in self.headers:
            self.slice[k] = self.data[k][mask]
    def __str__(self):
        return str(self.data)
    def get_headers(self):
        return self.data.keys()
    headers = property(get_headers)
