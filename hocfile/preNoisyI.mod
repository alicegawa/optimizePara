TITLE Fluctuating current pre load

NEURON{
    POINT_PROCESS preNoisyI
    RANGE setEta
    RANGE iave,sigma
    RANGE del,dur
    ELECTRODE_CURRENT i
}

UNITS{
    (nA) = (nanoamp) 
    (mV) = (millivolt)
    (umho) = (micromho)
}

PARAMETER{
    dt			(ms)
    iave		(nA)
    sigma		(nA)
    del			(ms)
    dur			(ms)
}

ASSIGNED {
    v	(mV)
    i 	(nA)
    eta (1)
}

BREAKPOINT {
    at_time(del)
    at_time(del+dur)
    
    if (t < del + dur && t >= del) {
	i = iave + sigma*eta
    }else{
	i = 0
    }
}

PROCEDURE setEta(inputEta) {
    eta = inputEta
}