from neuron import nrn, h, hclass, run, init
#from nrn    import *
import sys


class Neuron:
    
    def __init__(self):
        raise NotImplementedError( "Should have implemented init method to load neuron" )
    
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
    
    def get_domains(self):
        raise NotImplementedError( "Not implemented")
        


class L5PC(Neuron):
    
    def get_default_params(self):
        return {'model':'L5PCbiophys3.hoc','cellmorph':'cell1.asc'}
   
    def __init__(self,inparams): #TODO make model and cellmorph star params
        params = self.get_default_params()
        params.update(inparams)
        h.load_file('stdlib.hoc') 
        h.load_file('stdrun.hoc')
        h.load_file("import3d.hoc")
        h.load_file("models/%s"%params['model'])
        h.load_file("models/L5PCtemplate.hoc")
        
        self.cell = h.L5PCtemplate("morphologies/%s"%params['cellmorph'])
        
    def get_cell(self):
        return self.cell
    
    def get_soma(self):
        return self.cell.soma
    
    def get_basal(self):
        return self.cell.basal
    
    def get_apical(self):
        return self.cell.apic
    
    def get_apic(self):
        return self.cell.apic

    def get_axon(self):
        return self.cell.axon
    
    def locateSites(self, sectionString, distance):
        return self.cell.locateSites(sectionString,distance)


# Wrappers from Andrew Davidson neuronpy modified and commented by Romain Caze
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
    def __init__(self, L, diam, nseg=1, Ra=100, cm=1, mechanisms=[], parent=None, connection_point=1):
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

        nrn.Section.__init__(self)
        # set geometry
        self.L = L
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


class SimpleNeuron(Neuron):
    
    def create_soma(self,soma_params):
        return Section(soma_params['x'], soma_params['y'], Ra=soma_params['Ra'], mechanisms=soma_params['mechs'])
    
    def create_dendrites(self,dend_params):
        return Section(dend_params['dlength'], dend_params['ddiam'], nseg=dend_params['nseg'], Ra=dend_params['Ra'], parent=dend_params['parent'], mechanisms=dend_params['mechs'])
    
    def get_soma(self):
        return self.soma
    
    
class SimpleBipolar(SimpleNeuron):
    
    def get_default_params(self):
        hh = Mechanism('hh')
        pas = Mechanism('pas', {'e':-65,'g':0.0001})
        params = {'soma': {'x':10,
                           'y':10,
                           'Ra':150,
                           'mechs':[hh,pas]}, 
                  'dend0': {'dlength':50.,
                            'ddiam'  : 0.4,
                            'nseg':  10,
                            'Ra' : 150,
                            'mechs':[pas]},
                  'dend1': {'dlength':50.,
                            'ddiam'  : 0.4,
                            'nseg':  10,
                            'Ra' : 150,
                            'mechs':[pas]},}
        return params
    
    def __init__(self,inparams={}):
        params = self.get_default_params()
        params.update(inparams)
        self.cell = self.create_bipolar(params)
        self.params = params
        
    def create_bipolar(self,params):
        self.soma = self.create_soma(params['soma'])
        params['dend0']['parent'] = self.soma
        params['dend1']['parent'] = self.soma
        self.dend0 = self.create_dendrites(params['dend0'])
        self.dend1 = self.create_dendrites(params['dend1'])

    def count_dendrites(self):
        return 2        
        
    def get_apical(self):
        return self.dend0
    
    def get_basal(self):
        return self.dend1
        
        
class SimpleUnipolar(SimpleNeuron):
    
    def get_default_params(self):
        hh = Mechanism('hh')
        pas = Mechanism('pas', {'e':-65,'g':0.0001})
        params = {'soma': {'x':10,
                           'y':10,
                           'Ra':150,
                           'mechs':[hh,pas]}, 
                  'dend0': {'dlength':50.,
                            'ddiam'  : 0.4,
                            'nseg':  10,
                            'Ra' : 150,
                            'mechs':[pas]},
                  'dend1': {'dlength':50.,
                            'ddiam'  : 0.4,
                            'nseg':  10,
                            'Ra' : 150,
                            'mechs':[pas]},}
        return params
    
    def __init__(self,inparams={}):
        params = self.get_default_params()
        params.update(inparams)
        self.cell = self.create_unipolar(params)
        self.params = params
        
    def create_unipolar(self,params):
        self.soma = self.create_soma(params['soma'])
        params['dend0']['parent'] = self.soma
        self.dend0 = self.create_dendrites(params['dend0'])
        
    def count_dendrites(self):
        return 1
    
    def get_apical(self):
        return self.dend0

#TODO : could extend this to SimpleNPolar (i.e. n number of child dendrites))
