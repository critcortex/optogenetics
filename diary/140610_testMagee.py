from neuron import h
import Neuron
"""
###############################################
# Version 1: worked fine
h.load_file('stdlib.hoc') 
h.load_file('stdrun.hoc')
h.load_file("import3d.hoc")
h.xopen("comptogen_load_magree2000.hoc")

cell = h.magee_ca1

"""
###############################################
# Version 2: writing as a NeuronInterface

class CA1Magee(Neuron.NeuronInterface):
    
    def get_default_params(self):
        return {'v_init': -62.66}
    
    def __init__(self,inparams={}):
        params = self.get_default_params()
        params = self.dict_update(params,inparams)
        h.load_file('stdlib.hoc') 
        h.load_file('stdrun.hoc')
        h.load_file("import3d.hoc")
        h.xopen("comptogen_load_magree2000.hoc")
        
        self.cell = h.magee_ca1()
        self.soma = self.get_soma()
        
        super(CA1Magee,self).__init__()
        
        #self.populate_sectionlists()
        print "domain lists: ",self.domainlists
        ss = self.domainlists['soma']
        for sec in ss.allsec():
            print sec.name()
        
        
        # note that we have to regenerate seg_sections
        self.generate_seg_section()
        
        print 'Created CA1 (MAgee) neuron'



        
    def get_subdomains(self):
        return ['soma','d'] 
    
    def get_soma(self):
        return self.cell.soma
    
    #def get_subdomain_list(self,name):
    #    return self.get_subdomains()
    
    def populate_sectionlists(self):
        """
        self.domainlists = {}
        for subd in self.get_subdomains():
            self.domainlists[subd] = h.List() #h.SectionList()
        
        for (i,sec) in enumerate(self.allseclist):
            try:               
                # this is hacky but we have to do it this way, as a typical name is L5Template[3].apic[34] 
                matched_domain = [subd for subd in self.get_subdomains() if sec.name().find(subd)>=0][0]
                #self.domainlists[matched_domain].append(sec=sec) ### fine only when list is a SectionList
                self.domainlists[matched_domain].append(sec)
            except:
                print '%s does not seem to belong to any subdomain'%sec.name()
                continue
           """
        pass 
            
tmp = CA1Magee()


###############################################
# Version 3: 
def test_CA1():
    from neuron import h
    
    cell = CA1Magee()
    print cell.somapos
    print cell._calc_totnsegs()
    #clamp_site = cell.select_section_posn_bydistance('d')[0]
    #iclamp = h.IClamp()
    #iclamp.loc(clamp_site[2], sec=clamp_site[1])
    
    return cell

#tmp = test_CA1() 
