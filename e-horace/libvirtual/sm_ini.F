* sm_ini.F
* all parameters for the Standard Model are initialized here
* this file is part of FormCalc
* last modified 19 Jul 01 th
	subroutine model_ini(Ecms, reset, *)
	implicit none
	double precision Ecms
	logical reset

	double precision ALPHAS2
	external ALPHAS2
#include "model.h"

* Warning: With the PDG values used here, the CKM matrix is
* not strictly unitary. In cases where GIM cancellations are
* important, this parametrization may give wrong results!

	data CKM /
     +    1d0,   0D0, 0D0,
     +    0D0,   1D0, 0D0,
     +    0D0,   0D0, 1D0 /

c	data MLE  / ME,  MM,  ML  /
c	data MLE2 / ME2, MM2, ML2 /
c	data MQU  / MU,  MC,  MT  /
c	data MQU2 / MU2, MC2, MT2 /
c	data MQD  / MD,  MS,  MB  /
c	data MQD2 / MD2, MS2, MB2 /


	MLE(1) = ME
	MLE(2) = MM
	MLE(3) = ML
	MLE2(1) = ME2
	MLE2(2) = MM2
	MLE2(3) = ML2


	MQU(1) = MU
	MQU(2) = MC
	MQU(3) = MT
	MQU2(1) = MU2
	MQU2(2) = MC2
	MQU2(3) = MT2


	MQD(1) = MD
	MQD(2) = MS
	MQD(3) = MB
	MQD2(1) = MD2
	MQD2(2) = MS2
	MQD2(3) = MB2


	


	Alfas = ALPHAS2(Ecms)
	Alfas2 = Alfas**2
	GS = sqrt(4*pi*Alfas)

	if(.not. reset) return

	EL = sqrt(4*pi*Alfa)

* boson masses

c	MH2 = MH**2
	SW = sqrt(SW2)

#ifdef DEBUG

	print *, "============================"
	print *, "SM parameters"
	print *, "----------------------------"
	print *, "1/Alfa =", 1/Alfa
	print *, "Alfas(s) =", Alfas
	print *, "GF =", GF
	print *, "SW2 =", SW2
	print *, "----------------------------"
	print *, "Gauge-boson masses:"
	print *, "MW =", MW
	print *, "MZ =", MZ
	print *, "----------------------------"
	print *, "Fermion masses:"
	print *, "ME =", ME
	print *, "MM =", MM
	print *, "ML =", ML
	print *, "MU =", MU
	print *, "MC =", MC
	print *, "MT =", MT
	print *, "MD =", MD
	print *, "MS =", MS
	print *, "MB =", MB
	print *, "----------------------------"
	print *, "Higgs mass:"
	print *, "MH =", MH
	print *, "============================"

#endif

	end

	double precision function alphas2(S)
	double precision S, alphasmz

	alphas2 = 0.118d0



	end
