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

	double complex mwcmplx

#include "model.h"


	common /resonant/ mwcmplx

* Warning: With the PDG values used here, the CKM matrix is
* not strictly unitary. In cases where GIM cancellations are
* important, this parametrization may give wrong results!

*	data CKM /
*     +    .97495D0, .2225D0,  .0035D0,
*     +    .222D0,   .97415D0, .040D0,
*     +    .009D0,   .039D0,   .99915D0 /
	data CKM /
     +    1D0, 0d0, 0d0,
     +    0d0, 1d0, 0D0,
     +    0D0, 0D0, 1D0 /

	data MLE  / ME,  MM,  ML  /
	data MLE2 / ME2, MM2, ML2 /

	data MQU  / MU,  MC,  MT  /
	data MQU2 / MU2, MC2, MT2 /

	data MQD  / MD,  MS,  MB  /
	data MQD2 / MD2, MS2, MB2 /

	Alfas = ALPHAS2(Ecms)
	Alfas2 = Alfas**2
	GS = sqrt(4*pi*Alfas)

	if(.not. reset) return

	EL = sqrt(4*pi*Alfa)


* boson masses

	mwcmplx = dcmplx(mz,-gz/2d0)

	MH2 = MH**2
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

	function alphas2(ecm)
	implicit none
	double precision alphas2,ecm,b0,mz,asmz,pi
	integer nf

	pi = 3.14159265d0
	nf=5
	asmz = 0.118d0
	b0 = 11d0-2d0/3d0*nf
	mz = 91.187
	alphas2 = asmz/(1d0+b0*asmz/2d0/pi*log(ecm/mz))

	end
