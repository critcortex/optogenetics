NEURON {
	POINT_PROCESS SpikeOut
	GLOBAL thresh, refrac, vrefrac, grefrac
	NONSPECIFIC_CURRENT i
}

PARAMETER {
	thresh = 1 (millivolt)
	refrac = 5 (ms)
	vrefrac = 0 (millivolt)
	grefrac = 100 (microsiemens) :clamp to vrefrac
}

ASSIGNED {
	i (nanoamp)
	v (millivolt)
	g (microsiemens)
}

INITIAL {
	net_send(0, 3)
	g = 0
}

BREAKPOINT {
	i = g*(v - vrefrac)
}

NET_RECEIVE(w) {
	if (flag == 1) {
		net_event(t)
		net_send(refrac, 2)
		v = vrefrac
		g = grefrac
	}else if (flag == 2) {
		g = 0
	}else if (flag == 3) {
		WATCH (v > thresh) 1
	}	
}
