TITLE Mod file Km deactivation by M1 AChR activation

COMMENT

	Model of the Kv7 (KCNQ) Potassium Channel and its deactivation due to 
	activation of m1 Muscarinic Acetylcholine Receptor (mOxomR) by acetylcholine
	(Oxom). This mod file implements a kinetic model of the G protein cascade
	with equations and parameters from [1] and [2]. 	

	[1]	B. H. Falkenburger, J. B. Jensen, and B. Hille, “Kinetics of PIP2 
	metabolism and KCNQ2/3 channel regulation studied with a 
	voltage-sensitive phosphatase in living cells.,” J. Gen. Physiol., 
	vol. 135, no. 2, pp. 99–114, 2010.

	[2] B. H. Falkenburger, J. B. Jensen, and B. Hille, “Kinetics of M1 
	muscarinic receptor and G protein signaling to phospholipase C in living 
	cells,” J. Gen. Physiol., vol. 135, no. 2, pp. 99–114, 2010.

	[3] E. K. Pissadaki, K. Sidiropoulou, M. Reczko, and P. Poirazi, “Encoding 
    of spatio-temporal input characteristics by a CA1 pyramidal neuron model,” 
    PLoS Comput. Biol., vol. 6, no. 12, 2010.

	[4] B. H. Falkenburger, E. J. Dickson, and B. Hille, “Quantitative 
	properties and receptor reserve of the DAG and PKC branch of Gq-coupled 
	receptor signaling,” J Gen Physiol, vol. 141, no. 5, pp. 537–555, 2013.

	[5] M. Kruse, O. Vivas, A. Traynor-Kaplan, and B. Hille, “Dynamics of 
	Phosphoinositide-Dependent Signaling in Sympathetic Neurons,” J. Neurosci., 
	vol. 36, no. 4, pp. 1386–1400, 2016.

	Km channel mechanism from kmb.mod M. Migliore June 2006

    Kinetics have been altered from [5] in order to speed up PIP2 formation and
    hydrolysis to match CA1 pyramidal cell behavior. This was done by Adam 
    Mergenthal on 4/18/2017.

ENDCOMMENT

NEURON {
	THREADSAFE SUFFIX m1_kmb
	USEION k READ ek WRITE ik
	USEION ACh READ AChi VALENCE 0
	RANGE rip3
    RANGE gbar, ik, perc_i
	RANGE R, G, Gbg, RG, RGbg, RLA, RLGbgA
	RANGE DAG, PIP2, PIP2_KCNQ, Ga_GTP, PLC, Ga_GTP_PLC
	RANGE perc_open
}

UNITS {
	(mA)	= (milliamp)
	(mV) 	= (millivolt)
}

PARAMETER {

	R_init = 15.87  : (um^-2)
    G_init = 40  : (um^-2) determined as a ratio see [2]
	PLC_init = 3.12  : (um^-2) 

	: Parameters from [2]

	alpha = 100
	KL1A = 0.2

	k_oxom_diff = 0.05  : (msec^-1)

	kf_L1 = 2.78  : (mM^-1 msec^-1) kr_L1/KL1
	kr_L1 = 0.00555  : (msec^-1)
	
	kr_L1A = 0.00215  : (msec^-1)
	
	kf_L2 = 2.78  : (mM^-1 msec^-1)
	kr_L2 = 0.0000556  : (msec^-1) kf_L1*KL2

	kf_L2A = 2.78  : (mM^-1 ms^-1)

	kf_G1 = 0.00000000027  : (um^2 ms^-1)
	kr_G1 = 0.0068  : (ms^-1)

	kf_G2 = 0.0000027  : (um^2 ms^-1)
	kr_G2 = 0.00068  : (ms^-1)

	kf_G2A = 0.0000027  : (um^2 ms^-1)
	kr_G2A = 0.00068  : (ms^-1)

	k_NX_RLG = 0.00065  : (ms^-1)
	k_NX_RLGA = 0.00065  : (ms^-1)
	k_NX_G = 0.000000015  : (ms^-1)
	k_NX_P = 0.0047  : (ms^-1)

	k_GTPase1 = 0.000026  : (ms^-1)
	k_GTPase2 = 0.015  : (ms^-1)

	k_PLCassoc = 0.01  : (um^2 ms^-1)
	k_PLCdiss = 0.0071  : (ms^-1) 
	PLC_basal = 0
	PLC_efficiency = 1

	k_reconst = 0.001  : (um^2 s^-1)

	: Parameters from [1]
	k_PLC = 0.03  : 0.0003 (um^2 ms^-1)
	k_PLConPIP = 0.000014  : (um^2 ms^-1)

	k_4K = 0.00008  : (ms^-1)
	k_4P = 0.012  : (ms^-1)
	k_5K = 0.002  : (ms^-1)
	k_5P = 0.0028  : (ms^-1)

	speed_PIP2_buffer = 1  : (ms^-1) [4]
	fold_PIP2 = 3  : [4,5]

	speed_KCNQ_PIP2 =  0.00005 : 0.00005  : (um^2 ms^-1)
	KD_KCNQ_PIP2 = 2000  : (um^-2)

	k_IP3ase = 1.292 : 0.646  : 0.00013  : (ms^-1)
	k_DAGase = 0.02  : 0.0002  : (ms^-1)

    t_inh = 10
    perc_test = 1.0
	perc_open_ss = 0.3812

	gbar=0.0001		(mho/cm2)
	vhalfl=-40		(mV)
	kl=-10
	vhalft=-42		(mV)
	a0t=0.003		(/ms)
	zetat=7
	gmt=0.4
	q10=5
	b0=60
	st=1
	sh=0
}

ASSIGNED {

    celsius (degC)
	area
	v 		(mV)
    ik 		(mA/cm2)
	inf
	tau
	taua
	taub
	ek		(mV)

	kf_L1A
	kr_L2A
	KL2A

	oxom_diff

	L1
	L1A
	L2
	L2A
	G1
	G2
	G2A
	L2b
	L2bA
	G1b
	G2b
	G2bA	
	NX_G
	NX_RG
	NX_RLG
	NX_RLGA
	NX_P
	PLC_assoc
	PLC_diss
	G_reconst
	GTPase1
	GTPase2
	K_4
	P_4
	K_5
	P_5
	VSP
	rate_PLC
	rip3
	rate_PLConPI4P

	DAGase

	PIP2_buffer
    tot_KCNQ
	KCNQ_bind_PIP2
	perc_open
}

STATE {
	: M1 mOxomR model states

	R
	RL
	RLA
	RG
	RLG
	RLGA

	AChi
	Oxom

	RGbg
	RLGbg
	RLGbgA

	G
	Gbg

	Ga_GTP
	Ga_GDP
	Ga_GTP_PLC
	Ga_GDP_PLC
		
	PLC
	
	: PIP2 model states
	PI
	PI4P
	PIP2
	PIP2_bound
	DAG
	
	: KV7 (KCNQ) model states
	KCNQ
	PIP2_KCNQ

	m
}

BREAKPOINT {
	SOLVE states METHOD cnexp

	tot_KCNQ = PIP2_KCNQ + KCNQ
	perc_open = (PIP2_KCNQ/tot_KCNQ)^2

 	ik = perc_test*perc_open*gbar*m*(v - ek)
}

DERIVATIVE states {
	rates(v)
    m' = (inf-m)/tau

	R' = -L1 - L1A - G1 - G1b
	RL' = L1 - G2 - G2b
	RG' = G1 - L2 - NX_RG
	RLG' = L2 + G2 - NX_RLG

	RLA' = L1A - G2A - G2bA
	RLGA' = L2A + G2A - NX_RLGA
	RLGbgA' = L2bA + G2bA + NX_RLGA 

	RGbg' = G1b - L2b + NX_RG
	RLGbg' = L2b + G2b + NX_RLG

	G' = -G1 - NX_G - G2 + G_reconst
	Gbg' = NX_G - G1b - G2b - G_reconst
	
	Ga_GTP' = NX_G + NX_RG + NX_RLG + NX_RLGA - GTPase1 - PLC_assoc

	Ga_GDP' = GTPase1 + PLC_diss - G_reconst

	Ga_GTP_PLC' = PLC_assoc + NX_P - GTPase2
	Ga_GDP_PLC' = GTPase2 - NX_P - PLC_diss

	PLC' = -PLC_assoc + PLC_diss

	: PI' = P_4 - K_4 Clamped in [1]
	PI4P' = K_4 - P_4 + VSP + P_5 - K_5 - rate_PLConPI4P
	PIP2' = K_5 - P_5 - VSP - rate_PLC - KCNQ_bind_PIP2 - PIP2_buffer
	PIP2_bound' = PIP2_buffer

	DAG' = -DAGase + rate_PLC + rate_PLConPI4P

	: Kv7 State Derivatives
	PIP2_KCNQ' = KCNQ_bind_PIP2
    KCNQ' = -KCNQ_bind_PIP2

}

INITIAL {
	: M1 mOxomR model states
	
	kf_L1A = kr_L1A/KL1A  : (uM^-1 ms^-1)
	KL2A = KL1A/alpha 
	kr_L2A = kf_L1A*KL2A  : (s^-1) kf_L1A*KL2A

	Oxom = 0

	R = R_init
	RL = 0
	RG = 0
	RLG = 0

	RLA = 0
	RLGA = 0

	RGbg = 0
	RLGbg = 0
	RLGbgA = 0

	G = G_init
	Gbg = 0

	Ga_GTP = 0
	Ga_GDP = 0 
	Ga_GTP_PLC = 0
	Ga_GDP_PLC = 0
		
	PLC = PLC_init

	Ga_GTP = 0
	Ga_GDP = 0
	Ga_GTP_PLC = 0
	Ga_GDP_PLC = 0

	PI = 226975
	PI4P = 4540
	PIP2 = 3232
	PIP2_bound = 6464

	DAG = 13  : (uM^-2)
	
	: KV7 (KCNQ) model states
	KCNQ = 1.53 : 4
	PIP2_KCNQ = 2.47  : 0

	m = inf

    rates(v)
}

PROCEDURE rates(v) {
	: Equations from Kv7 model used in [3]
    LOCAL a, b
    a = exp(0.0378*zetat*(v - vhalft))
    b = exp(0.0378*zetat*gmt*(v - vhalft))
    inf = 1/(1 + exp((v-vhalfl)/kl))
    tau = b0 + b/(a0t*(1 + a))

	: M1 Model Rates [2]

	L1 = kf_L1*R*Oxom - kr_L1*RL
	L1A = kf_L1A*R*AChi - kr_L1A*RLA
	L2 = kf_L2*Oxom*RG - kr_L2*RLG
	L2A = kf_L2A*AChi*RG - kr_L2A*RLGA
	G1 = kf_G1*G*R - kr_G1*RG
	G2 = kf_G2*G*RL - kr_G2*RLG
	G2A = kf_G2*G*RLA - kr_G2*RLGA

	G1b = kf_G1*Gbg*R - kr_G1*RGbg
	G2b = kf_G2*Gbg*RL - kr_G2*RLGbg
	G2bA = kf_G2A*Gbg*RLA - kr_G2A*RLGbgA
	L2b = kf_L2*Oxom*RGbg - kr_L2*RLGbg

	NX_RLG = k_NX_RLG*RLG
	NX_RLGA = k_NX_RLGA*RLGA
	NX_G = k_NX_G*G
	NX_RG = k_NX_G*RG
	NX_P = k_NX_P*Ga_GDP_PLC

	GTPase1 = k_GTPase1*Ga_GTP
	GTPase2 = k_GTPase2*Ga_GTP_PLC

	PLC_assoc = k_PLCassoc*PLC*Ga_GTP
	PLC_diss = k_PLCdiss*Ga_GDP_PLC

	G_reconst = k_reconst*Gbg*Ga_GDP

	: G Protein Cascade Rates [1]

	K_4 = fold_PIP2*k_4K*PI
	P_4 = k_4P*PI4P

	K_5 = fold_PIP2*k_5K*PI4P
	P_5 = k_5P*fold_PIP2*PIP2

	rate_PLC = PIP2*(PLC_basal + (k_PLC*fold_PIP2*Ga_GTP_PLC))
	rip3 = rate_PLC*area
	rate_PLConPI4P = PI4P*PLC_efficiency*(PLC_basal + Ga_GTP_PLC*k_PLC)

    PIP2_buffer = (fold_PIP2 - 1)*speed_PIP2_buffer*PIP2 - speed_PIP2_buffer*PIP2_bound

	: Kv7 (KCNQ) Rates [1]

	KCNQ_bind_PIP2 = speed_KCNQ_PIP2*PIP2*KCNQ - speed_KCNQ_PIP2*KD_KCNQ_PIP2*PIP2_KCNQ
}


