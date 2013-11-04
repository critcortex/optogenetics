: NpHR 3-state model 2012-06-11

NEURON  {
	POINT_PROCESS NpHR
	ELECTRODE_CURRENT i
	RANGE irradiance, flux, g, gph, gv, iNpHR 
    	RANGE gbar, gpump                        
    	RANGE e1, k20, k3dark, k3light 
    	RANGE lightdelay, pulsewidth, interpulse_interval, n_pulses   
}

UNITS {
	(nA) = (nanoamp)
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

: illumination  (it is changed via hoc file)
   irradiance = 1	            :(mW/mm2)	  
   lightdelay = 30                (ms)
   pulsewidth = 100               (ms) <0, 1e9>  : duration of ON phase
   interpulse_interval = 50               (ms) <0, 1e9>  : duration of OFF phase
   n_pulses = 2                                : how many to deliver
  
    
: rate parameters  (/ms)
    e1 = 0.01         (/ms)     
    k20 = 1          (/ms)
    k3dark  = 0.0002  (/ms)	
    k3light = 0.0135  (/ms)

: conductance    
    gpump = 0.02 : fitting parameter for the peak current (in nA) 
                : recorded at saturation light @ -70mv voltage clamp (e.g. =1)   
	
    gbar = 1 : defined in the main prog as HR_expression*Area

: light constants   
    Er0 = 0.1   : (mW/mm2)
:    lambda = 470 	   :(nm)
:    flux0 = 2e16    :(photons/(sec*cm2))


: reversal potential for NpHr is about -400mV    
    e_rev = -400        :(mV)  
    
}

ASSIGNED {

    flux
    k1   (/ms)
    k2   (/ms)  
    k3   (/ms)
    
    gph
    gv
    g    
    i    (nA)
    v    (mV) 
    
    on
    tally         : how many more to deliver
   	 
    iNpHR
	
}

STATE { cl op de }

BREAKPOINT {
    SOLVE kin METHOD sparse
    gph = op*gpump     : light dependency op=[0:1] 
    gv =  e_rev-v      : voltage dependency 
    g = gph*gv
    i = gbar*g         : -ve ELECTRODE_CURRENT i hiperpolarizes the cell 
    iNpHR=-i	       : for plotting
}

INITIAL { 
    	cl=1
  	op=0
  	de=0
	
	flux = 0
  	i = 0
  	tally = n_pulses
  	if (tally > 0) {
	      net_send(lightdelay, 1)
      	      on = 0
	      tally = tally - 1
  	}

}

 
KINETIC kin {
    rates()
    ~ cl <-> op (k1, 0)
    ~ op <-> de (k2, 0)
    ~ de <-> cl (k3, 0)
    CONSERVE cl + op + de = 1
}


PROCEDURE rates() {

    k1 = e1*flux/Er0
    k2 = k20
    k3 = k3dark + k3light*log(1+flux/Er0)  

}  

NET_RECEIVE (w) {
       : ignore any but self-events with flag == 1
       
       if (flag == 1) {
          if (on == 0) {
          
             : turn it on   
             flux = irradiance 
             on = 1
             
             : prepare to turn it off
             net_send(pulsewidth, 1)
          } else {
             : turn it off
             flux = 0   
             on = 0
             if (tally > 0) {
                : prepare to turn it on again
                net_send(interpulse_interval, 1)
                tally = tally - 1
             }
          }
       }
    }