INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX hh1
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
    NONSPECIFIC_CURRENT il
    RANGE gnabar, gkbar, gl, el, gna, gk
    GLOBAL minf, hinf, ninf, mtau, htau, ntau
    POINTER dMemPot
}

UNITS {
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gnabar = .12 (S/cm2)
    gkbar = .036 (S/cm2)
    gl = .0003 (S/cm2)
    el = -54.3 (mV)
    ena = 50    (mV)
    ek = -90   (mV)
    gsoma = 100 (nS)
    celsius = 36    (degC)
    dt              (ms)
    v               (mV)
    vtraub  = -63   (mV)
}

STATE {
    m h n
    sMemPot (mV)
}

ASSIGNED {
    ina     (mA/cm2)
    ik      (mA/cm2)
    il      (mA/cm2)
    gna (S/cm2)
    gk (S/cm2)
    minf
    hinf
    ninf
    mtau (ms)
    htau (ms)
    ntau (ms)
    tadj	
    dMemPot (mV)    
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = gnabar*m*m*m*h
    ina = gna*(v - ena)
    gk = gkbar*n*n*n*n
    ik = gk*(v - ek)
    il = gl*(v - el)
}

INITIAL {
    rates(v, celsius)
    m = minf
    h = hinf
    n = ninf
    sMemPot = -78.8361 (mV)
}

DERIVATIVE states {
    rates(v, celsius)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau
    n' = (ninf-n)/ntau
    sMemPot' = -(ina+ik+il) + gsoma * (dMemPot - sMemPot)
}

PROCEDURE rates(v, celsius)
{
    LOCAL  alpha, beta, sum, q10
    q10 = 3^((celsius - 6.3)/10)
    :"m" sodium activation system
    alpha = .1 * vtrap(-(v+40),10)
    beta =  4 * exp(-(v+65)/18)
    sum = alpha + beta
    mtau = 1/(q10*sum)
    minf = alpha/sum
    :"h" sodium inactivation system
    alpha = .07 * exp(-(v+65)/20)
    beta = 1 / (exp(-(v+35)/10) + 1)
    sum = alpha + beta
    htau = 1/(q10*sum)
    hinf = alpha/sum
    :"n" potassium activation system
    alpha = .01*vtrap(-(v+55),10)
    beta = .125*exp(-(v+65)/80)
    sum = alpha + beta
    ntau = 1/(q10*sum)
    ninf = alpha/sum
}

:FUNCTION vtrap(x,y) {
    : use built in exprelr(z) = z/(exp(z)-1), which handles the z=0 case correctly
:    vtrap = y*exprelr(x/y)
:}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(Exp(x/y)-1)
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
} 