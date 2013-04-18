: ChR 6-state model 2012-06-11

NEURON  {
	POINT_PROCESS ChR
	ELECTRODE_CURRENT i
	RANGE Er, flux, g, gph, gv, iChR 
    	RANGE gbar                        
    	RANGE e1, e3, b10, a2dark, a2light, b2dark, b2light, a30, a40, gama
    	RANGE delay, pulsewidth, interpulse_interval, n_pulses   
}

UNITS {
	(nA) = (nanoamp)
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

: illumination  (it is changed via hoc file)
   Er = 1	            :(mW/mm2)	 irradaince 
   delay = 30                (ms)
   pulsewidth = 100               (ms) <0, 1e9>  : duration of ON phase
   interpulse_interval = 50               (ms) <0, 1e9>  : duration of OFF phase
   n_pulses = 1                                : how many to deliver
  
    
: rate parameters  (/ms)
    e1 = 0.004        (/ms)     
    b10 = 0.08          (/ms)  : Gd1 O1->C1
    a2dark  = 0.012  (/ms)	
    a2light = 0.010  (/ms)
    b2dark  = 0.006  (/ms)	
    b2light = 0.003  (/ms)

    a30 = 0.08   : Gd2 O2->C2
    e3  = 0.003
    a40 = 0.00033  : 1/taudark, taudark=3s

: conductance    
    gama = 0.05
    gbar = 1 : defined in the main prog as HR_expression*Area

: light constants   
    Er0 = 0.1   : (mW/mm2)
:    lambda = 470 	   :(nm)
:    flux0 = 2e16    :(photons/(sec*cm2))


: reversal potential for NpHr is about -400mV        e_rev = -400        :(mV)  
: inward rectifier ChR2 conductance    
   v0 = 43   		  :(mV)
   v1= -4.1       

}

ASSIGNED {

    flux
    
    a11   (/ms)
    a12   (/ms)  
    b1   (/ms)
    
    a2
    b2

    a3
    b31
    b32
    a4

    gph
    gv
    g    
    i    (nA)
    v    (mV) 
    
    on
    tally         : how many more to deliver
   	 
    iChR
}

STATE { c1 i1 o1 o2 i2 c2 }

BREAKPOINT {
    SOLVE kin METHOD sparse
    gph = o1+gama*o2     : light dependency op=[0:1] 
    gv =  (1-exp(-v/v0))/v1       : voltage dependency (equal 1 at Vm=-70mV)
    g = gph*gv
    i = gbar*g       : +ve ELECTRODE_CURRENT derpolarise the cell 
    iChR=-i	       : for plotting
}

INITIAL { 
    	c1=1
  	o1=0
	o2=0
	c2=0 	
	i1=0
	i2=0

	flux = 0
  	i = 0
  	tally = n_pulses
  	if (tally > 0) {
	      net_send(delay, 1)
      	      on = 0
	      tally = tally - 1
  	}

}

 
KINETIC kin {
    rates()
    ~ c1 <-> i1 (a11, 0)
    ~ i1 <-> o1 (a12, 0)
    ~ o1 <-> c1 (b1, 0)      
    ~ o1 <-> o2 (a2, b2)    
    ~ o2 <-> i2 (0, b31)
    ~ i2 <-> c2 (0, b32)
    ~ o2 <-> c2 (a3, 0)    
    ~ c2 <-> c1 (a4, 0)
    CONSERVE c1 + i1 + o1 + o2 + i2 + c2 = 1

}


PROCEDURE rates() {

    a11 = e1*flux/Er0
    a12 = 1
    b1  = b10

    a2 = a2dark + a2light*log(1+flux/Er0)  
    b2 = b2dark + b2light*log(1+flux/Er0)  

    a3  = a30
    b31 = 1
    b32 = e3*flux/Er0
    a4  = a40
}  

NET_RECEIVE (w) {
       : ignore any but self-events with flag == 1
       if (flag == 1) {
          if (on == 0) {
          
             : turn it on   
             flux = Er 
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