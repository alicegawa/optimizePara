TITLE  Belmabrouk et al. 2011
:GPeA is changed based on Belmanbrouk et al. 2011
:
: Na+, K, Ca_T, Leakage, Ca diffusion, SK, and A(K)
: 



NEURON {
    SUFFIX GPeA
    NONSPECIFIC_CURRENT ilk
    USEION ca READ cai, cao WRITE ica, cai
    USEION k READ ki, ko WRITE ik
    USEION na READ nai, nao WRITE ina
    RANGE ina, ik, ica, iA
    RANGE gnabar, ena, m_inf, a_m, v_m, ena_fixed       : fast sodium
    RANGE gkdrbar, ek, W_inf,tau_W, ikD, lambda, sigma, a_W, v_W, ek_fixed               : delayed K rectifier
    RANGE gl, el, ilk                                    : leak
    RANGE gcatbar, eca, icaT, X, X_inf, tau_X, a_X, v_X : T-type ca current
    RANGE gkcabar, ek, ikAHP, Csk, q_inf, C_gamma, tau_sk, a_sk                      : ca dependent SK current
    RANGE kca, vol, caGain                               : ca dynamics
    RANGE A_inf, B, B_inf, tau_B, gA, a_A, a_B, v_A, v_B :A(K) leak current
    RANGE offset
    RANGE coef_log, coef_scale_C, coef_intercept
    RANGE power_q, qinf_scale
}


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
    (molar) = (1/liter)
    (mM)	= (millimolar)
    FARADAY = (faraday) (coulomb)  :units are really coulombs/mole
    PI	= (pi) (1)
}

PARAMETER {
    R = 8.31441 (Gas constant)
    T 		(Absolute temp)
    celsius		(degC)
    
    :Fast Na channel
    gnabar   = 120e-3 (S/cm2):49e-3 (S/cm2) 
    a_m = 0.065
    v_m = -31 (mV)
    ena_fixed = 55 (mV)
    : delayed K rectifier 
    gkdrbar  = 5e-3 (S/cm2):57e-3	(S/cm2)
    sigma = 1 (1)
    a_W = 0.055
    v_W = -46 (mV)
    ek_fixed = -70 (mV)
    lambda = 0.05 :the parameter of time constant here0.3 and init +10 ha antei
    :Leakage current
    gl	= 0.3e-3 (S/cm2):0.35e-3	(S/cm2)
    el	= -60	(mV)
    
    :Ca dynamics
    kca   = 2        (1/ms)
    area
    vol = 3.355e-11  (L) :~20um radius sphere
    caGain = .1
    
    :T-type ca current
    gcatbar   = 0.1e-3 (S/cm2):5e-3 (S/cm2)  
    tau_X = 50 (ms)
    a_X = 0.071
    v_X = -10 (mV)
    :AHP current (Ca dependent K current)
    gkcabar   = 2e-3 (S/cm2) :1e-3 (S/cm2) 
    C_gamma = 113e-6(mM)
    tau_sk =  250 (ms)
    a_sk = 2.5e-4 (1):0.08 (1)
    :A(K)
    gA =  220e-3(S/cm2):22e-3 (S/cm2)
    tau_B = 20 (ms)
    a_A = 0.02
    a_B = -0.1
    v_A = -20 (mV)
    v_B = -70 (mV)
    counter = 0 (1)
    offset = 0 (mM)
    coef_log = 2.508 (1)
    coef_scale_C = 0.05
    coef_intercept = 1.120
    power_q = 2
    q_inf_scale = 1.15
}

ASSIGNED {
    v	(mV)
    ina	(mA/cm2)
    ik	(mA/cm2) 
    ikD	(mA/cm2)   
    ikAHP	(mA/cm2)  
    ica	(mA/cm2) 
    icaT	(mA/cm2) 
    ilk	(mA/cm2)
    iA      (mA/cm2)
    
    :Fast Na
    m_inf
    ena           (mV)   := 60
    
    :K rectifier
    W_inf
    tau_W	(ms)
    ek         (mV) := -90
    
    :T-type ca current
    eca           (mV)   :calc from Nernst
    X_inf
    :AHP (Ca dependent K current)
    q_inf (mM)
    :k
    A_inf
    B_inf
}

STATE {
    W 
    cai (mM) <1e-10>
    cao (mM) <1e-10>
    nai (mM) <1e-10>
    nao (mM) <1e-10>
    ki (mM) <1e-10>
    ko (mM) <1e-10>
    
    X
    Csk (mM) <1e-10>
    B
}


BREAKPOINT {
    SOLVE states METHOD cnexp
    
    T = 273 + celsius - 9.5
    :ena = -(R*T)/FARADAY*log(nai/nao)*1000
    :ek = (R*T)/FARADAY*log(ko/ki)*1000
    :eca = -(R*T)/FARADAY*log(cai/cao)*1000/2
    
    ina = gnabar * m_inf*m_inf*m_inf*(1-W) * (v - ena)
   
    ikD = gkdrbar * (W/sigma)^4 * (v - ek)
    ikAHP = gkcabar * (v - ek)*q_inf^power_q:*q_inf
    iA = gA*A_inf*B*(v-ek)
    ik=ikD+ikAHP+iA
    ilk = gl * (v - el)
    ica = gcatbar*X*(v-eca)    
}

DERIVATIVE states {   
    evaluate_fct(v)
    W' = (W_inf - W)/tau_W
    X' = (X_inf - X)/tau_X
    :(Ica mA/cm2)*(area um2)*(1e-8 cm2/um2)*(1e-3 A/mA)*(1/(2*F) mol/C)*(1e-3 sec/msec)*(1e3 mMol/mol)(1/volume 1/L)=(mM/msec)
    cai' = caGain*(-ica*area*1e-11/(2*FARADAY*vol) - kca*cai)
    Csk' = -a_sk*ica-(Csk-C_gamma)/tau_sk
    B' = (B_inf - B)/tau_B
}

UNITSOFF

INITIAL {
    evaluate_fct(v)
    W = W_inf   
    X = X_inf
    Csk = C_gamma 
    B = B_inf
}

PROCEDURE evaluate_fct(v(mV)) {     
    m_inf = 1/(1+exp(-2*a_m*(v-v_m)))
    W_inf = 1/(1+exp(-2*a_W*(v-v_W)))
    tau_W = 1/(lambda*(exp(a_W*(v-v_W))+exp(-a_W*(v-v_W))))
    X_inf = 1/(1+exp(-2*a_X*(v-v_X)))
    
    q_inf = q_inf_scale/(1+exp(-coef_intercept - coef_log*llog((Csk-C_gamma-offset)/coef_scale_C)))
    :q_inf = 1.2/(1+exp(-1.120-2.508*llog((Csk-C_gamma)/0.00005)))
    
    A_inf = 1/(1+exp(-2*a_A*(v-v_A)))
    B_inf = 1/(1+exp(-2*a_B*(v-v_B)))
    counter = counter + 1
    : if(counter == 1 && t < 10){
    : 	printf("#time\tA_inf\tB_inf\tCsk\tCr\tq_inf\tcai\tica\tina\tik\tX\n")
    : 	printf("#setting of parmas\n")
    : 	printf("#a_sk = %g, gkdrbar = %g, gA = %g\n",a_sk, gkdrbar, gA)
    : }
    if(counter > 40 && t >= 100 && t<= 2000){
	:print every 1 milisecond
	counter = 0
	:printf("llog((Csk-C_gamma-offset)/0.050) = %g\n", llog((Csk-C_gamma-offset)/0.050))
	:printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", t, A_inf, B_inf, Csk, C_gamma, q_inf, cai, ica, ina, ik, X)
    }
}

FUNCTION llog(x) {
    if(x<1e-3){
	llog = -7
    }else{
	llog = log(x)
    }
}
    

UNITSON