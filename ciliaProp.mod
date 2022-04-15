TITLE CILIAProps

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ciliaProp
	POINTER cilMemPot1, cilMemPot2, cilMemPot3, cilMemPot4
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
}
 
STATE {
	ciliaMemPot (mV) : Vcilia
}

ASSIGNED {
	cilMemPot1 (mV)
	cilMemPot2 (mV)
	cilMemPot3 (mV)
	cilMemPot4 (mV)
	cilMemPot5 (mV)
	cilMemPot6 (mV)
}

BREAKPOINT {
	:SOLVE ciliaEq METHOD cnexp
	ciliaMemPot = (cilMemPot1 + cilMemPot2 + cilMemPot3 + cilMemPot4) / 4
}

:DERIVATIVE ciliaEq {
	:ciliaMemPot' = cilMemPot1' + cilMemPot2'
:}

INITIAL {
	ciliaMemPot = -78.8361 (mV)
}
	

