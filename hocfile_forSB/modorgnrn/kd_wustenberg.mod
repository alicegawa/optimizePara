TITLE Delayed potassium current from wu¨stenberg_baxter_2004

UNITS {
      (mA) = (milliamp)
      (mV) = (millivolt)
}
 
NEURON {
       SUFFIX kd_wustenberg
       USEION k READ ek WRITE ik
       RANGE  gkbar
} 
 
PARAMETER {
	  gkbar = 8.11 (mho/cm2)
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
        LOCAL  tmp_1
        TABLE minf, mtau, hinf, htau DEPEND dt, celsius FROM -120 TO 120 WITH 700

	tmp_1 = exp( (-20.1 - v) / 16.1 )
        minf = 1 / ( (1 + tmp_1) ^ 3 )
	tmp_1 = exp( (v - 20) / 20 )
	mtau = ( (5.0 - 0.5) / (1 + tmp_1) ) + 0.5

	tmp_1 = exp( (v + 74.7) / 7 )
        hinf = 1 / (1 + tmp_1)
	tmp_1 = exp( (v - 52) / 15 )
	htau = ( (200 - 150) / (1 + tmp_1) ) + 150
}

UNITSON