import sys
from neuron import h
import Neuron
import numpy as np
import pylab


h.load_file('stdlib.hoc', 'String') 
h.load_file('stdrun.hoc')
h.load_file("import3d.hoc")

# <codecell>

pp = {'num_split': 2, 
      'dend0': {'mechanisms': ['pas_dend', 'hh'],'diam':1.,'L':600}, 
       'dend1': {'mechanisms': ['pas_dend', 'hh'],'diam':1.,'L':600}, 
       'soma': {'mechanisms': ['pas_soma', 'hh'], 'diam':30,'L':30, 'Ra': 200, 'cm':1}, 
       'mechanisms': {'pas_dend': ('pas', {'e': -45,'g':0.0005}), 'pas_soma': ('pas', {'e': -65,'g':0.0005}),'hh' : ('hh',{}) } , 
       'num_base': 2, 
       'hh': ('hh',), 
       'defaultlevel': {'mechanisms': ['pas_dend', 'hh'],'diam':1.,'L':600}, 
       'num_levels': 3}


# <codecell>

def generate_membrane_mechanisms(params):
        # generate Mechansim objects
        mechanisms = {}
        for (mechkey,val) in params['mechanisms'].iteritems():
            mechanisms[mechkey] = Neuron.Mechanism(val[0],val[1])
        # replace listings in params dict with generated newly Mechanism objects
        # TODO: improve this method, as it is name dependent and therefore frail
        for (key,val) in params.iteritems():
            if type(val) is dict and val.has_key('mechanisms'):
                params[key]['mechanisms'] = [ mechanisms[mk] for mk in val['mechanisms']]
        return params

# <codecell>

def test_single():
    pp = generate_membrane_mechanisms(pp)
    
    # <codecell>
    
    soma = Neuron.Section(**pp['soma'])
    #pp['dend0']['parent'] = soma
    #dend0 = Neuron.Section(**pp['dend0'])
    
    # <codecell>
    
    
    stim = h.IClamp()
    stim.loc(0.5, sec=soma)
    stim.dur = 400
    setattr(stim, 'del',200)
    stim.amp = 0.509
    
    # <codecell>
    
    
    time = h.Vector()
    vmsoma = h.Vector()
    
    time.record(h._ref_t)
    vmsoma.record (soma(0.5)._ref_v)
    
    
    h.init()
    h.tstop = 700
    h.run()
    
    vms = np.array(vmsoma)
    tt = np.array(time)
    pylab.plot(tt,vms)
    
    pylab.savefig('test.png')


def scan(Ra,g_pas,diam,Ia):
    pp = {'num_split': 2, 
      'dend0': {'mechanisms': ['pas_dend', 'hh'],'diam':1.,'L':600}, 
       'dend1': {'mechanisms': ['pas_dend', 'hh'],'diam':1.,'L':600}, 
       'soma': {'mechanisms': ['pas_soma', 'hh'], 'diam':30,'L':30, 'Ra': 200, 'cm':1}, 
       'mechanisms': {'pas_dend': ('pas', {'e': -45,'g':0.0005}), 'pas_soma': ('pas', {'e': -65,'g':0.0005}),'hh' : ('hh',{}) } , 
       'num_base': 2, 
       'hh': ('hh',), 
       'defaultlevel': {'mechanisms': ['pas_dend', 'hh'],'diam':1.,'L':600}, 
       'num_levels': 3}
    
    pp['soma']['Ra'] = float(Ra)
    pp['soma']['diam'] = float(diam)
    pp['soma']['L'] = float(diam)
    pp['mechanisms']['pas_soma'] = ('pas', {'e': -65,'g':float(g_pas)})

    pp = generate_membrane_mechanisms(pp)
    
    # <codecell>
    
    soma = Neuron.Section(**pp['soma'])
    #pp['dend0']['parent'] = soma
    #dend0 = Neuron.Section(**pp['dend0'])
    
    # <codecell>
    
    
    stim = h.IClamp()
    stim.loc(0.5, sec=soma)
    stim.dur = 500
    setattr(stim, 'del',200)
    stim.amp = float(Ia)
    
    # <codecell>
    
    
    time = h.Vector()
    vmsoma = h.Vector()
    
    time.record(h._ref_t)
    vmsoma.record (soma(0.5)._ref_v)
    
    
    h.init()
    
    
    
    h.tstop = 900
    h.run()
    
    vms = np.array(vmsoma)
    tt = np.array(time)
    pylab.plot(tt,vms)
    
    pylab.savefig('140307_scan_Ra%s_gpas%s_diam%s_Ia%s.png'%(Ra,g_pas,diam,Ia))
    
def create_jdf(Ra,g_pas,diam,Ia):
    import file_io as fio
    from subprocess import call
    expname = '140307_scan_Ra%g_gpas%g_diam%g_Ia%g'%(Ra,g_pas,diam,Ia)
    f = open('%s_sim.jdf'%expname, 'w')
    f.write(fio.get_pbs_header()%(expname+'_sim',expname,expname,24))
    f.write('cd git/optogenetics/ \n python diary/140307_test_soma_only.py run_exp %g %g %g %g'%(Ra,g_pas,diam,Ia))
    f.close()
    print 'Created ','%s_sim.jdf'%expname
    call(['qsub','%s_sim.jdf'%expname])

def scan_range():
    Ras = np.arange(100,301,50)
    g_pases = [0.2,0.002,0.0002,0.00002]
    diams = np.arange(10,201,10)
    Ias = np.arange(0.5,10.1,0.5)
    for Ra in Ras:
        for g_pas in g_pases:
            for diam in diams:
                for Ia in Ias:
                    create_jdf(Ra,g_pas,diam,Ia)
    
    
       
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        pass
    elif sys.argv[1] == 'run_exp':
        scan(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
                                
                            
                            
  