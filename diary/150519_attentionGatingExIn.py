
import Neuron
from run_experiments import ExpSetup
import run_analysis
import NeuroTools.stgen as ntst 
import numpy as np

expbase = '150519_attentionGatingExIn'

"""
Purpose:

Model for having input of three different types to Pyramidal cell (L5PC) neuron:
- feedforward (input)
- recurrent 
- feedback (attention)
and for:
- keeping input rate constant
- varying recurrent and feedback rates
- with and without optogenetics

Intent: that the change in gain modulation will be feed back into changes 
in recurrent population -- which will in turn modulate this cell

To consider:
-- how many inputs for basal vs apic - correct for 
-- synchrony between input signals?

"""

stg = ntst.StGen()
es = ExpSetup()


tstop = 2500# 2200
light_dur = 1000
light_on = 1200


L5PC_areas = ['soma', 'apic','dend'] # NB: we don't include 'axon'

TOTAL_NUM_AREAS = {'L5PC': 3 }
areas = {'L5PC': L5PC_areas}

celltype = 'L5PC' 


gradients = np.arange(0.0,0.00101,0.0001)
"""
balanceEI = 4. # balance of E:I is 4:1
num_sites_ex = 800 # total number of connections
ratioEI = 1./(balanceEI + 1.)
num_sites_in = int(num_sites_ex*ratioEI)
num_sites_total = num_sites_ex + num_sites_in

############ FOR Excitatory
num_ff = int(0.1*num_sites_ex) # which will all be in the apical dendrites
num_fb = int(0.1*num_sites_ex) # which will all be in the basal dendrites
num_rr = int(0.8*num_sites_ex) # which will be split between the apical and dendrites
# how we split our recurrent connections between apical and basal
rr_apic = int(0.5*num_rr)
rr_basal = num_rr - rr_apic
# so then the total number of synaptic sites on the apical and basal dendrites are:  
num_apic = num_ff + rr_apic
num_basal = num_fb + rr_basal
############ FOR Inhibitory
i_num_ff = int(0.1*num_sites_in) # which will all be in the apical dendrites
i_num_fb = int(0.1*num_sites_in) # which will all be in the basal dendrites
i_num_rr = int(0.8*num_sites_in) # which will be split between the apical and dendrites
# how we split our recurrent connections between apical and basal
i_rr_apic = int(0.5*i_num_rr)
i_rr_basal = i_num_rr - i_rr_apic
# so then the total number of synaptic sites on the apical and basal dendrites are:  
i_num_apic = i_num_ff + i_rr_apic
i_num_basal = i_num_fb + i_rr_basal

Jex = 0.1
g = 3.0 # set to be the same as balanceEI 
Jin = -g*Jex
"""
class AttentionGating():

    def get_populations(self,num_sites_ex,ratioEI):
        self.num_sites_ex = num_sites_ex
        self.ratioEI = ratioEI
        self.num_sites_in = int(self.num_sites_ex*self.ratioEI)
        self.num_sites_total = self.num_sites_ex + self.num_sites_in
    
        ############ FOR Excitatory
        self.num_ff = int(0.1*num_sites_ex) # which will all be in the apical dendrites
        self.num_fb = int(0.1*num_sites_ex) # which will all be in the basal dendrites
        self.num_rr = int(0.8*num_sites_ex) # which will be split between the apical and dendrites
        # how we split our recurrent connections between apical and basal
        self.rr_apic = int(0.5*self.num_rr)
        self.rr_basal = self.num_rr - self.rr_apic
        # so then the total number of synaptic sites on the apical and basal dendrites are:  
        self.num_apic = self.num_ff + self.rr_apic
        self.num_basal = self.num_fb + self.rr_basal
        ############ FOR Inhibitory
        self.i_num_ff = int(0.1*self.num_sites_in) # which will all be in the apical dendrites
        self.i_num_fb = int(0.1*self.num_sites_in) # which will all be in the basal dendrites
        self.i_num_rr = int(0.8*self.num_sites_in) # which will be split between the apical and dendrites
        # how we split our recurrent connections between apical and basal
        self.i_rr_apic = int(0.5*self.i_num_rr)
        self.i_rr_basal = self.i_num_rr - self.i_rr_apic
        # so then the total number of synaptic sites on the apical and basal dendrites are:  
        self.i_num_apic = self.i_num_ff + self.i_rr_apic
        self.i_num_basal = self.i_num_fb + self.i_rr_basal
    
    
    
    
    def get_spiketimes(self,rate,num_inputs):
        spikes = []
        for i in range(num_inputs):
            ss = stg.poisson_generator(rate,t_stop=tstop,array=True)
            spikes.append(ss)
        #print 'Spikes = ',spikes
        return spikes
    
    
    def _run_simulation(self,irr,factor,freq_ff,freq_rr,freq_fb,adist,bdist,Jex,Jin,num_sites):
        """
        @params
            irr       irradiance
            factor    NpHR factor
        """
        
        pp = {}
        pp['cell'] = ['Neuron',celltype]
        pp['cell_params'] = {}
        pp['tstart'] = 0
        pp['tstop'] = tstop
        
        # opsin
        chrdict =  {'exp':5e-4, 'irradiance' :irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':0,  'n_pulses':1}
        hrdict =  {'exp':5e-4, 'irradiance' :factor*irr, 'pulsewidth': light_dur,'lightdelay':light_on,'interpulse_interval':0,  'n_pulses':1}
        
        pp['opsindict'] = {}
        pp['opsindict']['ChR'] =  {}
        for area in areas[celltype]:
            pp['opsindict']['ChR'][area] = chrdict    
        pp['ChR_areas'] = {'whole'      : pp['opsindict']['ChR'].keys()}
        
        pp['opsindict']['NpHR'] =  {}
        for area in areas[celltype]:
            pp['opsindict']['NpHR'][area] = hrdict
        pp['NpHR_areas'] = {'whole'      : pp['opsindict']['NpHR'].keys()}
        
        # general settings 
        pp['experiment_type'] = 'opsinonly' #TODO change
        pp['savedata'] = True 
        pp['logdata'] = False
        
        
        # mark locations for PSP throughout dendritic trees, along with soma
        pp['mark_loc'] = {}
        #labels = ['estim%g'%i for i in range(num_sites_ex)] + ['istim%g'%i for i in range(num_sites_in)]
        labels = ['estim_apic','estim_basal','istim_apic','istim_basal']
        full_labels = ['estim_apic%g'%i for i in range(self.num_apic)] +['estim_basal%g'%i for i in range(self.num_basal)] +['istim_apic%g'%i for i in range(self.i_num_apic)]+ ['istim_basal%g'%i for i in range(self.i_num_basal)]
        pp['mark_loc']['names'] = ['mysoma']+labels # both ex and inh sites
        pp['mark_loc']['sections'] = ['soma']+['apic']*self.num_apic+['dend']*self.num_basal+['apic']*self.i_num_apic+['dend']*self.i_num_basal
        # going to separate this as Python is complaining about        
        ids_ex_apic = [('select_section_posn_bydistance',{'sectionarea':'apic','mindist':adist[0],'maxdist':adist[1],'nsites':self.num_apic}) ] 
        ids_ex_basal = [('select_section_posn_bydistance',{'sectionarea':'dend','mindist':bdist[0],'maxdist':bdist[1],'nsites':self.num_basal})] 
        ids_inh_apic = [('select_section_posn_bydistance',{'sectionarea':'apic','mindist':adist[0],'maxdist':adist[1],'nsites':self.i_num_apic})] 
        ids_inh_basal = [('select_section_posn_bydistance',{'sectionarea':'dend','mindist':bdist[0],'maxdist':bdist[1],'nsites':self.i_num_basal})]        
        pp['mark_loc']['ids'] = [(0,0.5)] + ids_ex_apic + ids_ex_basal + ids_inh_apic + ids_inh_basal
        
        pp['stim_spiketrains'] = True                            
        tstims_ex = self.get_spiketimes(freq_ff,self.num_ff) + self.get_spiketimes(freq_rr,self.num_rr) + self.get_spiketimes(freq_fb,self.num_fb)
        tstims_inh = self.get_spiketimes(freq_ff,self.i_num_ff) + self.get_spiketimes(freq_rr,self.i_num_rr) + self.get_spiketimes(freq_fb,self.i_num_fb)
        pp['spiketrains'] = [{'tstims': tstims_ex+tstims_inh,  'locations': full_labels, 'weights':np.concatenate([np.ones(self.num_sites_ex)*Jex,np.ones(self.num_sites_in)*Jin]), 'el': 0.02}]
    
    
        # set up recordings
        # - record voltage from soma
        pp['record_loc'] = {}
        pp['record_loc']['v'] = ['mysoma']
        
        pp['plot'] = { 'voltages': [ ['mysoma', 'v','k-'] ]} 
        
        pp.update({'expname':expbase,
                       'description':'_irr%.3f_factor%.2f_sites%g_ff%g_rr%g_fb%g_Jex%.1f_Jin%.1f'%(irr,factor,self.num_sites_total,freq_ff,freq_rr,freq_fb,Jex,Jin)})
        ##print pp['description']
        es.run_single_experiment(expbase, 'missing', pp)




ff_rates = [10,20,50,60,80]
scan_rates = range(20,101,20) # for fb and rr
irr = 0.001
factor = 0.6
g = 3.0 # set to be the same as balanceEI 
ratioEI = 0.25
EG = AttentionGating()

def runt():
    count = 0
    for num_sites in range(100,501,100):
        EG.get_populations(num_sites, ratioEI)            
        for ffr in ff_rates:
            for rrr in scan_rates:
                for fbr in scan_rates:
                    for J in np.arange(0.1,0.31,0.1):
                        EG._run_simulation(irr,factor,ffr,rrr,fbr,[100,-1],[100,-1],J,J*-g,num_sites)
                        count += 1
                        #return
                
    print count
    
    
runt()





