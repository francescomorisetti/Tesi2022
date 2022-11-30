* model.h
* common blocks for the model parameters
* this file is part of FormCalc
* last modified 30 Oct 01 th
	double precision pi, degree, sqrt2, hbar_c2
	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (degree = pi/180D0)
	parameter (sqrt2 = 1.41421356237309504880168872421D0)
	parameter (hbar_c2 = 3.8937966D8)
*         = hbar c^2 in picobarn
	double complex cI
	parameter (cI = (0D0, 1D0))

#ifndef polar
#define polar(r, theta) r*exp(theta*degree*cI)
#endif
* SM parameters
	double precision EL, Alfa, Alfa2, GF, GS, Alfas, Alfas2
        double precision MW, MW2, MZ, MZ2, GW, GZ
	double precision SW, SW2, CW, CW2
	double precision MH, MH2, MG0, MG02, MGp, MGp2
	double precision ME, ME2, MM, MM2, ML, ML2, MLE(3), MLE2(3)
	double precision MU, MU2, MC, MC2, MT, MT2, MQU(3), MQU2(3)
	double precision MD, MD2, MS, MS2, MB, MB2, MQD(3), MQD2(3)
	double complex CKM(3, 3)
	common/feynartshorace/Alfa, Alfa2, GF,
     . MW, MW2, MZ, MZ2, GW, GZ, SW, SW2, CW, CW2,
     . MH, MH2, 
     . ME, ME2, MM, MM2, ML, ML2, 
     . MU, MU2, MC, MC2, MT, MT2, 
     . MD, MD2, MS, MS2, MB, MB2
c$$$	parameter (Alfa = 1/137.03599911D0, Alfa2 = Alfa**2)
c$$$	parameter (GF = 1.16637D-5)
c$$$	parameter (MZ = 91.1876D0, MZ2 = MZ**2)
c$$$	parameter (MW = 80.425D0, MW2 = MW**2)
c$$$	parameter (GW = 2.124D0)
c$$$	parameter (SW2 = (MZ2 - MW2)/MZ2)
c$$$	parameter (CW = MW/MZ, CW2 = CW**2)
c$$$	parameter (MM = 0.51099892D-3, MM2 = MM**2)
c$$$	parameter (ME = 105.658369D-3, ME2 = ME**2)
c$$$	parameter (ML = 1776.99D-3, ML2 = ML**2)
c$$$	parameter (MU = 66.D-3, MU2 = MU**2)
c$$$	parameter (MC = 1.2D0, MC2 = MC**2)
c$$$	parameter (MT = 178D0, MT2 = MT**2)
c$$$	parameter (MD = 66.D-3, MD2 = MD**2)
c$$$	parameter (MS = 150D-3, MS2 = MS**2)
c$$$	parameter (MB = 4.3D0, MB2 = MB**2)
c$$$	parameter (MH = 115D0, MH2 = MH**2 )
	common /sm_para/
     +    CKM, MLE, MQU, MQD, MLE2, MQU2, MQD2,
     +    EL, GS, Alfas, Alfas2,
     +    MG0, MG02, MGp, MGp2
* MSSM parameters

	double complex USf(2, 2, 4, 3)
	double complex UCha(2, 2), VCha(2, 2), ZNeu(4, 4)
	double complex Af(4, 3), Au, Ad, MUE
	double precision M_1, M_2
	double precision MSf(2, 4, 3), MSf2(2, 4, 3)
	double precision MCha(2), MNeu(4), MCha2(2), MNeu2(4)
	double precision MSNE(3), MSLE1(3), MSQU1(3), MSQD1(3)
	double precision MSLE2(3), MSQU2(3), MSQD2(3)
	double precision Mh0, MHH, MA0, MHp, MGl, MSusy
	double precision Mh02, MHH2, MA02, MHp2, MGl2
	double precision CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	double precision CA, SA, CA2, SA2, C2A, S2A
	double precision CAB, SAB, CBA, SBA

	common /mssm_para/
     +    USf, UCha, VCha, ZNeu,
     +    Af, Au, Ad, MUE,
     +    M_1, M_2,
     +    MSf, MSf2,
     +    MCha, MNeu, MCha2, MNeu2,
     +    MSNE, MSLE1, MSQU1, MSQD1,
     +    MSLE2, MSQU2, MSQD2,
     +    Mh0, MHH, MA0, MHp, MGl, MSusy,
     +    Mh02, MHH2, MA02, MHp2, MGl2,
     +    CB, SB, TB, CB2, SB2, TB2, C2B, S2B,
     +    CA, SA, CA2, SA2, C2A, S2A,
     +    CAB, SAB, CBA, SBA

#ifndef CKMC
#define CKMC(a, b) dconjg(CKM(a, b))
#define USfC(a, b, t, g) dconjg(USf(a, b, t, g))
#define VChaC(a, b) dconjg(VCha(a, b))
#define UChaC(a, b) dconjg(UCha(a, b))
#define ZNeuC(a, b) dconjg(ZNeu(a, b))
#endif
