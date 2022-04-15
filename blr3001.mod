TITLE BLR rate

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON { 
	SUFFIX blr3001
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(pA) = (picoampere)
	(molar) = (1/liter)
        (uM) = (micromolar)
	(pC) = (picocoulomb)
	(nS) = (nanosiemens)
	(nF) = (nanofarad)
}

PARAMETER {
	capInv = 18.727 (/nF)
	cc1lin = 11.4109 (/s)
	cc2 = 10.5600 (/s)
	ck1lin = 0.7781 (/s)
	ck2 = 0.3430 (/s)
	clmax = 2.3131 (nS)
	cnmax = 1.1266 (nS)
	cx1lin = 0.8088 (/s)
	cx2 = 9.2763 (/s)
	ef = 1.5789 (/s)
	gl = 0.8056 (nS)
	hmc1 = 10.5250 (uM)
	hmc2 = 19.4752 (uM)
	inf = 1.0162 (uM/pC)
	inhmax = 1.5105
	k1 = 0.0011 (/uM-s)
	k2 = 12.1312 (/s)
	kI = 1.8815 (uM)
	kinh = 1.0074 (uM)
	kinhcng = 2.3202 (uM)	
	n1 = 2.9406 
	n2 = 1.2006
	nI = 1.5708
	ninh = 4.8356
	ninhcng = 2.8548
	pd = 0.6752 (/s)
	r1 = 1.1133 (/s)
	r2 = 5.2182 (/s)
	smax = 189.0821 (uM/s)
	vcl = -0.0246 (mV)
	vcng = -2.9351e-7 (mV)
	vl = -78.8361 (mV)
:from code of Multiscale, as values are not given in paper
	Rtot = 1
	Gtot = 1
:parameters for odor stimulus
:reduced sharpness to 0.01 for avoiding out of range error
	sharpness = 0.01 (ms)
	ostim = 300.0 (uM)
	ipulse
	hv1
	hv2
	od (uM)
	Icng (pA)
	inhcng (1)
	JNCX (/s)
	cc1 (/s)
	ck1 (/s)
	cx1 (/s)
	IClCa (pA)
	IL (pA)
	syncAMP (uM/s)
	kG (uM/s)
	rG (uM/s)
	synth (/s)
	degrad (/s)
}

STATE {
	blr (uM)
	aG (uM)
	aCaMK (uM)
	cAMP (uM)
	CaCaM (uM)
	Ca (uM)
	IX (uM)	
	memPot (mV)
}

ASSIGNED {
	intTerm1
	intTerm2
:terms added for calculating powers separately using intermediate terms and pow()
	powTerm1
	powTerm2
	powTerm3
	powTerm4
	powTerm5
	powTerm6
	powTerm7
	powTerm8
	itot (pA)
	grad
}

BREAKPOINT {
:Heaviside function for single pulse
	hv1 = 1/(1+(exp(-(t-0.2)/sharpness)))
	hv2 = 1/(1+(exp(-(t-60)/sharpness)))
	ipulse = hv1 - hv2 
	:if ( t > 0.2 && t <= 1.2 ) {
	:	ipulse = 1
	:}
	:else {
	:	ipulse = 0
	:}
:Odor stimulus
	od = ostim * ipulse
	:capInv = 1/capCilia
	SOLVE blrrate METHOD cnexp

	:printf("t: %g, blr: %g, aG: %g, synth: %g, degrad: %g, syncAMP: %g, Icng: %g, IX: %g, v: %g \n", t, blr, aG, cAMP, Ca, CaCaM, aCaMK, IX, v)	
}

DERIVATIVE blrrate {
	:calcParameters(blr, aG, cAMP, aCaMK, Vcilia, CaCaM, Ca, IX)
	:printf("t: %g, blr: %g, aG: %g, cAMP: %g, Ca: %g, CaCaM: %g, aCaMK: %g, IX: %g, Vcilia: %g \n", t, blr, aG, cAMP, Ca, CaCaM, aCaMK, IX, Vcilia)
 	:printf("t: %g, syncAMP: %g, Icng: %g, inhcng: %g, JNCX: %g, cc1: %g, ck1: %g, cx1: %g, IClCa: %g, IL: %g, Itotal: %g grad: %g \n", t, syncAMP, Icng, inhcng, JNCX, cc1, ck1, cx1, IClCa, IL, Itotal, grad)
	:calcParameters()
	blr' = k1 * od * (Rtot - blr) - r1 * blr
	aG' = kG * (Gtot - aG) - rG
	cAMP' = syncAMP - pd * cAMP
	Ca' = inf * Icng - JNCX - (cc1 - cc2 * CaCaM)
	CaCaM' = cc1 - cc2 * CaCaM
	aCaMK' = ck1 - ck2 * aCaMK
	IX' = cx1 - cx2 * IX
	memPot' = capInv * (Icng+IClCa+IL) * 0.001
	calcParameters()
	:printf("t: %g, Icng: %g, IClCa: %g IL: %g \n", t, Icng, IClCa ,IL)
}

INITIAL {
	:calcParameters(blr, aG, cAMP, aCaMK, Vcilia, CaCaM, Ca, IX)
	blr=1.e-8 (uM)
	aG=1.e-8 (uM)
	cAMP=1.e-8 (uM)
	Ca=1.e-8 (uM)
	CaCaM=1.e-8 (uM)
	aCaMK=1.e-8 (uM)
	IX=1.e-8 (uM)
	memPot=-78.8361 (mV)
	inhcng = 0 
	calcParameters()
}

:PROCEDURE calcParameters(blr, aG, cAMP, aCaMK, Vcilia, CaCaM, Ca, IX) {
:calculation of intermediate terms for DE
PROCEDURE calcParameters() {
	kG = k2 * blr
	rG = r2 * aG
	intTerm1 = aCaMK/kinh
	:powTerm1 = calcExp(intFrac(aCaMK,kinh), ninh)
	powTerm1 = calcExp(intTerm1, ninh)
	syncAMP = aG * (smax/(1+powTerm1))
	:printf("t: %g, aG: %g", t , aG)
	powTerm2 = calcExp(cAMP, n1)
	intTerm2 = inhcng * hmc1
	:powTerm3 = calcExp(intProd(inhcng, hmc1), n1)
	powTerm3 = calcExp(intTerm2, n1)
	powTerm4 = calcExp(CaCaM, ninhcng)
	powTerm5 = calcExp(kinhcng, ninhcng)
	inhcng = 1 + (((inhmax-1) * powTerm4)/(powTerm4 + powTerm5))
	Icng = cnmax * (powTerm2 /(powTerm2  + powTerm3)) * (vcng - memPot)
	:intTerm3 = IX/kI
	powTerm6 = calcExp(intFrac(IX, kI),nI)
	JNCX = (ef * Ca)/(1+powTerm6)
	:JNCX = ef * Ca
	cc1 = intProd(cc1lin, Ca)
	ck1 = intProd(ck1lin, CaCaM)
	cx1 = intProd(cx1lin, Ca)
	powTerm7 = calcExp(Ca, n2)
	powTerm8 = calcExp(hmc2, n2)
	IClCa = clmax * (powTerm7 / (powTerm7 + powTerm8)) * (vcl - memPot)
	IL = intProd(gl, (vl - memPot))
	itot = (Icng+IClCa+IL)
	:printf("t: %g, Vcilia: %g \n", t, Vcilia)
	:printf("t: %g, powterm1: %g, powterm2: %g, powterm3: %g, powterm4: %g, powterm5: %g, powterm6: %g, powterm7: %g, powterm8: %g \n", powTerm1, powTerm2, powTerm3, powTerm4, powTerm5, powTerm6, powTerm7, powTerm8)
	:grad = Itotal * capInv
}

:to calculate power terms
FUNCTION calcExp(term1, term2) {
	calcExp = pow(term1, term2)
}

:to evaluate fractions
FUNCTION intFrac(term3, term4) {
	intFrac = term3/term4
}

:to evaluate product
FUNCTION intProd(term5, term6) {
	intProd = term5 * term6
}