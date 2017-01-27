TITLE Sustained potassium current from kloppenburg_horner_1998

UNITS {
      (mA) = (milliamp)
      (mV) = (millivolt)
}
 
NEURON {
       SUFFIX ks_kloppenburg
       USEION k READ ek WRITE ik
       RANGE  gkbar
} 
 
PARAMETER {
	  gkbar = 5.9145 (mho/cm2)
}
 
ASSIGNED {
	 v (mV)
	 ek (mV)
	 ik (mA/cm2)
	 gk (mho/cm2)
	 minf
	 mtau (ms)
	 hinf
	 htau (ms)
}
 
STATE {
      m
}
 
BREAKPOINT {
	   SOLVE states METHOD cnexp
	   gk = gkbar * m
	   ik = gk * (v - ek)
}
 
INITIAL {
	rates(v)
	m = minf
}

DERIVATIVE states {
	  rates(v)
	  m' = (minf - m) / mtau
}

UNITSOFF
 
PROCEDURE rates(v (mV)) {
        LOCAL tmp
        TABLE minf, mtau DEPEND dt, celsius FROM -120 TO 120 WITH 700

	tmp = exp(  -(v + 15.2) / 33.4  )
        minf = 1 / ( (1 + tmp) ^ 3 )

	tmp = exp(  (v - 10) / 17  )
	mtau = 1.6 + ( 10 / (1 + tmp) )
}

UNITSON