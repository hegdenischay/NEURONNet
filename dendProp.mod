TITLE DEND prop

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX dendProp
	POINTER ciliaMemPoten
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
	dcapInv = 1 (/nF)
	gdend = 10 (nS) 
}
 
STATE {
	dmemPot (mV) : Vdendrite
}

ASSIGNED {
	ciliaMemPoten (mV)
}

BREAKPOINT {
	SOLVE dendEq METHOD cnexp
}

DERIVATIVE dendEq {
	dmemPot' = dcapInv * gdend * (ciliaMemPoten - dmemPot) 
}

INITIAL {
	dmemPot=-78.8361 (mV)
}
	

