TITLE A type potassium current from kloppenburg_horner_1998

UNITS {
      (mA) = (milliamp)
      (mV) = (millivolt)
}
 
NEURON {
       SUFFIX ka_kloppenburg
       USEION k READ ek WRITE ik
       RANGE  gkbar
} 
 
PARAMETER {
	  gkbar = 1.228389 (mho/cm2)
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
      m h
}
 
BREAKPOINT {
	   SOLVE states METHOD cnexp
	   gk = gkbar * m * h
	   ik = gk * (v - ek)
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

	tmp = exp(  -(v + 20.4) / 17.4  )
        minf = 1 / ( (1 + tmp) ^ 3 )

	tmp = exp(  (v + 40.8) / 12.3  )
        hinf = 1 / (1 + tmp)

	tmp = exp(  (v + 9) / 10.5  )
	mtau = 5 + ( 60 / (1 + tmp) )

	htau = 7 * mtau
}

UNITSON