TITLE Sodium current from kloppenburg_horner_1998

UNITS {
      (mA) = (milliamp)
      (mV) = (millivolt)
}
 
NEURON {
       SUFFIX na_kloppenburg
       USEION na READ ena WRITE ina
       RANGE  gnabar
} 
 
PARAMETER {
	  gnabar = 19.997 (mho/cm2)
}
 
ASSIGNED {
	 v (mV)
	 ena (mV)
	 ina (mA/cm2)
	 gna (mho/cm2)
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
	   gna = gnabar * m * m * m * h
	   ina = gna * (v - ena)
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

	tmp = exp(  -(v + 40) / 5  )
        minf = 1 / (1 + tmp)
	tmp = exp(  (v + 23) / 6.5  )
	mtau = (( 6.5 / (1 + tmp) ) + 0.9) * 5

	tmp = exp(  (v + 45) / 15  )
        hinf = 1 / (1 + tmp)
	tmp = exp( (v + 20) / 8 )
	htau = (( 1.3 / (1 + tmp) ) + 0.49) * 5
}

UNITSON