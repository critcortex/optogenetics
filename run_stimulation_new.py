# -*- coding: utf-8 -*-
#


from neuron import *
from nrn    import *



def setup(tstart=300,tstop=1500,proxpt=400,distpt=620,vinit=-80,squareamp=1.,imax=0.5):
    
    
    global squareAmp, Imax, proximalpoint, distalpoint
    global risetau, decaytau, BACdt
    
    #====================== General files and tools =====================
    #h.load_file("nrngui.hoc")

    #====================== cvode =======================================
    h('objref cvode')

    h.cvode = h.CVode()
    h.cvode.active(1)

    #=================== creating cell object ===========================
    h.load_file('stdlib.hoc', 'String') 
    h.load_file('stdrun.hoc')
    h.load_file("import3d.hoc")
    h('objref L5PC')

    h('strdef morphology_file')
    h.morphology_file = "morphologies/cell1.asc"

    h.load_file("models/L5PCbiophys3.hoc")
    h('load_file("models/L5PCtemplate.hoc")')
    h.L5PC = h.L5PCtemplate(h.morphology_file)

    #=================== settings ================================
    v_init = vinit   # orig val: -80
    BACdt = 5
    
        
    #somatic pulse settings
    squareAmp = squareamp  # original value: 1.9

    #EPSP settings
    risetau = 0.5
    decaytau = 5
    Imax = imax

    proximalpoint = proxpt
    distalpoint = distpt

    h('objref tstart')
    h.tstart = tstart
    h.tstop = tstop




def set_experiment_type (experiment_type):
    global somastimamp, EPSPamp
    # Type = [BAP | CaBurst | BAC]
    
    ##experiment_type = $1
    print "Setting up for experiment type", experiment_type
    
    if (experiment_type=="BAP"):
        somastimamp = squareAmp
        EPSPamp = 0
    
    elif (experiment_type=="CaBurst"):
        somastimamp = 0
        EPSPamp = Imax*3
    
    elif (experiment_type=="BAC"):
        somastimamp = squareAmp
        EPSPamp = Imax
        
    elif (experiment_type=="opsinonly"):
        somastimamp = 0
        EPSPamp = 0
    
    else:
        print "WARNING: no experiment type was set"
        print "Assuming = somastimulation"
        
        somastimamp = squareAmp
        EPSPamp = 0
    
    


def add_optogenetics (opsindict):
    # Params: 
    #      optoexpress     = Dictionary of opsins, where keys are opsin names, 
    #                        while the value is a list of the locations the opsin is to 
    #                        be expressed in 
    # Example: 
    #h.load_file("NpHR-in-apical_distal.hoc")  
    #print "loaded NpHR apical successful"
    #h.HR_in_apic(5e-4,10,600,200,350,600)
    
    for (opsin, loclist) in opsindict.iteritems():
        for loc in loclist:
            
            if loc[0] is not None: # we allow users to specify None in order to make looping easier
                filename = '%s-in-%s.hoc'%(opsin,loc[0])
                print "filename =",filename
                try:
                    h.load_file(filename)
                    print "loaded %s %s successful"%(opsin,loc[0])
                except:
                    print "Error: could not load hoc file for opsin/location: ", filename
                ss = '%s_in_%s(%s)'%(opsin,loc[0],','.join([str(x) for x in loc[1]]))
                print ss
                h(ss)
            else:
                print "No locations specified for this opsin"
    

def set_distal_proximal ():
    h('objref sl')

    h.sl = h.List()
    h('double siteVec[2]')
    h('double site_dend_distal[2]')
    
    # find distal points
    h.sl = h.L5PC.locateSites("apic",distalpoint)
    axsites = h.sl.count()
    maxdiam = 0
    j= 0
    for i in range(axsites):
        dd1 = h.sl.object(i).x[1]
        dd = h.L5PC.apic[int(h.sl.object(i).x[0])](dd1).diam
        if (dd > maxdiam) :
            j = i
            maxdiam = dd 
    h.site_dend_distal[0] = h.sl.object(j).x[0]
    h.site_dend_distal[1] = h.sl.object(j).x[1]
    
    #h('double site_dend_distal[2]')
    #h.site_dend_distal = h.siteVec
    
    # find proximal points
    h('double site_dend_proximal[2]')
    h.sl = h.List()
    h.sl = h.L5PC.locateSites("apic",proximalpoint)
    maxdiam = 0
    j = 0
    axsites = int(h.sl.count())
    
    for i in range(axsites):
        dd1 = h.sl.object(i).x[1]
        dd = h.L5PC.apic[int(h.sl.object(i).x[0])](dd1).diam
        if (dd > maxdiam):
            j = i
            maxdiam = dd 
        

    h.site_dend_proximal[0] = h.sl.object(j).x[0]
    h.site_dend_proximal[1] = h.sl.object(j).x[1]    



def set_stimulus () :
    
    #======================== stimulus settings ============================
    #global siteVec, st1, st2, sl, syn1, ns, con1, isyn, tvec
    global vec
    vec = {}
    
    #Somatic pulse
    h('objref st1')
    
    #st1 = new IClamp(0.5) # from .hoc
    h.st1 = h.IClamp()
    h.st1.loc(0.5, sec=h.L5PC.soma[0])
    h.st1.amp = somastimamp
    #st1.del = tstart --> error as del is reserved word in python
    setattr(h.st1, 'del', h.tstart)
    h.st1.dur = h.tstop-h.tstart-100
    
    h('L5PC.soma st1')
    
    #Dendritic EPSP-like current
    h('objref sl,st2,ns,syn1,con1,isyn, tvec')

    h.isyn = h.Vector()
    h.tvec = h.Vector()
    h('double siteVec[2]')
    
    ## find distal points
    #h.sl = h.L5PC.locateSites("apic",distalpoint)
    #axsites = h.sl.count()
    #maxdiam = 0
    #j= 0
    #for i in range(axsites):
        #dd1 = h.sl.object(i).x[1]
        #dd = h.L5PC.apic[int(h.sl.object(i).x[0])](dd1).diam
        #if (dd > maxdiam) :
            #j = i
            #maxdiam = dd 
    #h.siteVec[0] = h.sl.object(j).x[0]
    #h.siteVec[1] = h.sl.object(j).x[1]
    
    #h('double site_dend_distal[2] = siteVec')
    
    ## find proximal points
    #h.sl = h.List()
    #h.sl = h.L5PC.locateSites("apic",proximalpoint)
    #maxdiam = 0
    #j = 0
    #axsites = int(h.sl.count())
    
    #for i in range(axsites):
        #dd1 = h.sl.object(i).x[1]
        #dd = h.L5PC.apic[int(h.sl.object(i).x[0])](dd1).diam
        #if (dd > maxdiam):
            #j = i
            #maxdiam = dd 
        

    #h.siteVec[0] = h.sl.object(j).x[0]
    #h.siteVec[1] = h.sl.object(j).x[1]
    #h('double site_dend_proximal[2] = siteVec')
        
    
    
    h('access L5PC.apic[site_dend_distal[0]]')

    h.st2 = h.IClamp(h.site_dend_distal[1])
    h.st2.amp = 0

    h('L5PC.apic[site_dend_distal[0]] { st2 } ')
    
    h.syn1 = h.epsp(h.siteVec[1])
    h.syn1.tau0 = risetau       
    h.syn1.tau1 = decaytau   
    h.syn1.onset = h.tstart + BACdt  
    h.syn1.imax = EPSPamp
    h('L5PC.apic[site_dend_distal[0]] { syn1 } ')
    
    h.cvode.record(h.syn1._ref_i,h.isyn,h.tvec)
    


def setup_record():
    
    #======================== recording settings ============================
    h('objref vsoma, vdend, recSite, vdend2, isoma,idend,idend2')

    ###########################################
    # Record voltages
    # - vsoma
    # - vdend
    # - vdend2
    ###########################################


    h.vsoma = h.Vector()
    h('access L5PC.soma')
    #tmpsoma = h.L5PC.soma
    #h.cvode.record(tmpsoma(0.5)._ref_v,h.vsoma,h.tvec)
    h('cvode.record(&v(0.5),vsoma,tvec)')
    
    h('access L5PC.apic[site_dend_distal[0]]')
    h.vdend = h.Vector()
    h.cvode.record(h.L5PC.apic[int(h.site_dend_distal[0])](h.site_dend_distal[1])._ref_v,h.vdend,h.tvec)

    h('access L5PC.apic[site_dend_proximal[0]]')
    h.vdend2 = h.Vector()
    h.cvode.record(h.L5PC.apic[int(h.site_dend_proximal[0])](h.site_dend_proximal[1])._ref_v,h.vdend2,h.tvec)

    
    h('access L5PC.apic[site_dend_proximal[0]]')
    h.recSite = h.IClamp(h.site_dend_proximal[1])
    h.recSite.amp = 0
    h('L5PC.apic[site_dend_proximal[0]] { recSite }')

    
    ###########################################
    # Record currents
    # - isoma
    # - idend
    # - idend2
    ###########################################
    
    h('access L5PC.soma')
    h.isoma = h.Vector()
    h.cvode.record(h.st1._ref_i,h.isoma,h.tvec)

    #h('access L5PC.apic[site_dend_distal[0]]')
    #h.idend = h.Vector()
    #h.cvode.record(h.L5PC.apic[int(h.site_dend_distal[0])](h.site_dend_distal[1])._ref_i,h.idend,h.tvec)

    #h('access L5PC.apic[site_dend_proximal[0]]')
    #h.idend2 = h.Vector()
    #h.cvode.record(h.L5PC.apic[int(h.site_dend_proximal[0])](h.site_dend_proximal[1])._ref_i,h.idend2,h.tvec)




def setup_plot ():
    
    
    pass
    #h('objref gV, gI, s')

    #h.s = h.Shape(shape=h.L5PC.all)
    #h.s.color_list(h.L5PC.axonal,2,2)
    #h.s.color_list(h.L5PC.somatic,5)
    #h.s.color_list(h.L5PC.basal,4)
    #h.s.color_list(h.L5PC.apical,1)
    #h.s.point_mark(h.st2,2)
    ##h.s.point_mark(h.recSite,3)
    #h.s.show(0)
    



def simulate_exp():
    
    print "Init and run experiment ...",
    h.init()
    print "...init'ed...",
    h.run()
    print "...finished"
    
    
    
def save_data(expname,savedata=False):
    if savedata:
        import numpy as np
        mat = np.matrix([h.tvec,h.vsoma,h.vdend,h.vdend2])
        mat.dump(expname+".dat")
        
    

def run_plots(savename=None):
    lw = 2
    import pylab
    pylab.figure(1)
    pylab.plot(h.tvec,h.vsoma,lw=lw,c='k')
    pylab.plot(h.tvec,h.vdend,lw=lw,c='r')
    pylab.plot(h.tvec,h.vdend2,lw=lw,c='b')
    pylab.xlim(h.tstart-20,h.tstop+20)
    pylab.ylim(-120,40)
    pylab.title('V')
    
    pylab.figure(2)
    pylab.plot(h.tvec,h.isyn,lw=lw)
    pylab.plot(h.tvec,h.isoma)
    pylab.xlim(h.tstart-20,h.tstop+20)
    pylab.ylim(-3,6)
    pylab.title('I')
    
    if savename is not None:
        pylab.figure(1)
        pylab.savefig(savename+'_voltage.png')
        pylab.figure(2)
        pylab.savefig(savename+'_current.png')
        print "Saved figures under %s*.png"%savename
        pylab.close('all')
    else:
        pylab.show()
    

def debug():
    print "tvec range = ",h.tvec.min(),h.tvec.max()
    print "isyn range = ",h.isyn.min(),h.isyn.max()
    print "vdend range = ",h.vdend.min(),h.vdend.max()
    print "vdend2 range = ",h.vdend2.min(),h.vdend2.max()
    print "vsoma range =", h.vsoma.min(), h.vsoma.max()
    print "isoma range = ",h.isoma.min(),h.isoma.max()
   

def main(keydict):
    setup(tstart=keydict['tstart'],tstop=keydict['tstop'])
    set_experiment_type(keydict['experiment_type'])
    set_distal_proximal ()
    set_stimulus()
    add_optogenetics(keydict['opdict'])
    setup_record()
    setup_plot()
    simulate_exp()
    save_data(keydict['expname'],keydict['savedata'])
    run_plots(keydict['expname'])
    debug()
    
if __name__ == '__main__':
    default_dict = {}
    default_dict['experiment_type'] = None
    default_dict['opdict'] = {}
    default_dict['expname'] = "defaulttest"
    main(default_dict)
