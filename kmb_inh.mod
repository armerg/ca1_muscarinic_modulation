TITLE CA1 KM channel from Mala Shah
: M. Migliore June 2006

NEURON {
	SUFFIX kmb_inh
	POINTER kppt, kcpt
	USEION kcnq READ kcnqi CHARGE 0
	USEION pip2_kcnq READ pip2_kcnqi CHARGE 0
	USEION k READ ek WRITE ik CHARGE 1
	RANGE  gbar, ik, sh, perc_test, perc_i, m, base_perc
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 		(mV)
	ek
	celsius 	(degC)
	gbar=.0001 	(mho/cm2)
    vhalfl=-40   	(mV)
	kl=-10
    vhalft=-42   	(mV)
    a0t=0.003      	(/ms)
    zetat=7    	(1)
    gmt=.4   	(1)
	q10=5
	b0=60
	st=1
	sh =0

	base_perc = 0.3816

	perc_test = 1.0
}


STATE {
	m
}

ASSIGNED {
	ik (mA/cm2)
	inf
	tau
	taua
	taub
	kppt
	kcpt
	perc_i
	kcnqi
	pip2_kcnqi
}

INITIAL {
	rate(v)
	m=inf
}

BREAKPOINT {
	SOLVE state METHOD cnexp

	perc_i = (pip2_kcnqi/(kcnqi + pip2_kcnqi))^2 / base_perc
	ik = perc_test*perc_i*gbar*m^st*(v-ek)
}

FUNCTION alpt(v(mV)) {
  	alpt = exp(0.0378*zetat*(v-vhalft-sh))
}

FUNCTION bett(v(mV)) {
  	bett = exp(0.0378*zetat*gmt*(v-vhalft-sh))
}

DERIVATIVE state {
	rate(v)
	m' = (inf - m)/tau
}

PROCEDURE rate(v (mV)) { :callable from hoc
	LOCAL a,qt
	qt=q10^((celsius-35)/10)
	inf = (1/(1 + exp((v-vhalfl-sh)/kl)))
	a = alpt(v)
	tau = b0 + bett(v)/(a0t*(1+a))
}














