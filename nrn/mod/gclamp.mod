COMMENT
conductance clamp defined by
    i = g * (v - e)      i(nanoamps), g(microsiemens);
ENDCOMMENT
					       
NEURON {
	POINT_PROCESS GClamp
	RANGE gamp, e, i, g
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	gamp=0 	(uS)	<0,1e9>
	e=0	(mV)
}

ASSIGNED { v (mV) i (nA)  g (uS)}

BREAKPOINT {
	g = gamp
	i = g*(v - e)
}
