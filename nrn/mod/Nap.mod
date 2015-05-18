COMMENT
-----------------------------------------------------------------------------
from ModelDB:

Uebachs M, Opitz T, Royeck M, Dickhof G, Horstmann MT, Isom LL, Beck H (2010) 
Efficacy Loss of the Anticonvulsant Carbamazepine in Mice Lacking Sodium Channel 
beta Subunits via Paradoxical Effects on Persistent Sodium Currents. J Neurosci 30:8489-501 

as used by
Welday AC, Shlifer IG, Bloom ML, Zhang K, Blair HT (2011) 
Cosine Directional Tuning of Theta Cell Burst Frequencies:
Evidence for Spatial Coding by Oscillatory Interference.
J. Neurosci. 31:16157-16176
-----------------------------------------------------------------------------
ENDCOMMENT

TITLE nap

NEURON {
	SUFFIX nap
	USEION na READ ena WRITE ina
	RANGE  gbar, thegna, sh, scalerate
	: GLOBAL minf, mtau 
	:, hinf, mtau, htau
}

PARAMETER {
	gbar = .0052085   	(mho/cm2)
	sh = 0  (mV)
	eNa = 55 	(mV)
	scalerate = 1	(ms)	
	ena		(mV)            
	celsius (degC)
	v 		(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	thegna		(mho/cm2)
	minf 	
	mtau (ms)
}
 

STATE { m }

UNITSOFF

BREAKPOINT {
    SOLVE states METHOD cnexp
	        		
	thegna =gbar*m       
	ina = thegna * (v - eNa)
	} 

INITIAL {
	mtau = scalerate
	minf = (1/(1+exp(-(v+52.3-sh)/6.8)))      	
	m=minf  
	
}

DERIVATIVE states {   
    
	mtau = scalerate
	minf = (1/(1+exp(-(v+52.3-sh)/6.8))) 	         	
	m' = (minf-m)/mtau
}



UNITSON
