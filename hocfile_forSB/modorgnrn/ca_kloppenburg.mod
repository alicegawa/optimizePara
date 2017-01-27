TITLE Calcium current from kloppenburg_horner_1998

UNITS {
      (mA) = (milliamp)
      (mV) = (millivolt)
}
 
NEURON {
       SUFFIX ca_kloppenburg
       USEION ca READ eca WRITE ica
       RANGE  gcabar
} 
 
PARAMETER {
	  gcabar = 2.54777 (mho/cm2)
}
 
ASSIGNED {
	 v (mV)
	 eca (mV)
	 ica (mA/cm2)
	 gca (mho/cm2)
	 minf
	 mtau (ms)
	 hinf
	 htau (ms)
}
 
STATE {
      m h
}
 
BREAKPOINT {
	   SOLVE states METHOD cnexp
	   gca = gcabar * m * m * m * h
	   ica = gca * (v - eca)
}
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

DERIVATIVE states {
	  rates(v)
	  m' = (minf - m) / mtau
	  h' = (hinf - h) / htau
}

UNITSOFF

PROCEDURE rates(v (mV)) {
        LOCAL tmp
        TABLE minf, mtau, hinf, htau DEPEND dt, celsius FROM -120 TO 120 WITH 700

	tmp = exp(  -(v + 17) / 4  )
        minf = 1 / (1 + tmp)
	tmp = exp(  (v + 22) / 14  )
	mtau = ( 32 / (1 + tmp) ) + 1.8

	tmp = exp(  (v + 30) / 12  )
        hinf = 1 / (1 + tmp)
	htau = 3 * mtau
}

UNITSON