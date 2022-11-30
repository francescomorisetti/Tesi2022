      subroutine getkappamuno(dkl,dku,dkd)
      include 'shared.inc'
      double precision kappamuno
      double complex fdg,aggfer
      include 'libzvirtual/uubar/renconst.h'

      deltauv = getdelta()
      call setdelta(0.d0)
      call uucalc_renconst  ! for Z

      c2    = 1.d0 - s2tw
      drho1 = dreal(dMWsq1-c2*dMZsq1)/mw/mw
      dk1   = -sqrt(c2/s2tw) * (-dZAZ1*mz*mz/2.d0-dZZA1*mz*mz/2.d0) / (mz*mz - aggfer(mz*mz))
      dk1   = dk1 + alpha/pi * c2/s2tw * log(c2)
      dk1   = dk1 - alpha/4.d0/pi/s2tw * (0.5d0*fdg(1.d0/c2) + 
     .        4.d0*c2*gdg(1.d0/c2) - 
     .     (1.d0 - 6.d0*s2tw+ 8.d0*s2tw*s2tw )/4.d0/c2 * fdg(1.d0))
      dkl = 8.d0*mw*mw*gf/sqrt(2.d0)*(dk1 + c2/s2tw*drho1)


      dk1   = -sqrt(c2/s2tw) * (-dZAZ1*mz*mz/2.d0-dZZA1*mz*mz/2.d0) / (mz*mz - aggfer(mz*mz))
      dk1   = dk1 + alpha/pi * c2/s2tw * log(c2)
      dk1   = dk1 - alpha/4.d0/pi/s2tw * ((0.5d0-s2tw/3.d0)*fdg(1.d0/c2) + 
     .        4.d0*c2*gdg(1.d0/c2) - 
     .     (1.d0 - 4.d0*s2tw + 32.d0/9.d0*s2tw*s2tw )/4.d0/c2 * fdg(1.d0))
      dku = 8.d0*mw*mw*gf/sqrt(2.d0)*(dk1 + c2/s2tw*drho1)


      dk1   = -sqrt(c2/s2tw) * (-dZAZ1*mz*mz/2.d0-dZZA1*mz*mz/2.d0) / (mz*mz - aggfer(mz*mz))
      dk1   = dk1 + alpha/pi * c2/s2tw * log(c2)
      dk1   = dk1 - alpha/4.d0/pi/s2tw * ((0.5d0-2.d0/3.d0*s2tw)*fdg(1.d0/c2) + 
     .        4.d0*c2*gdg(1.d0/c2) - 
     .     (1.d0 - 2.d0*s2tw + 8.d0/9.d0*s2tw*s2tw )/4.d0/c2 * fdg(1.d0))
      dkd = 8.d0*mw*mw*gf/sqrt(2.d0)*(dk1 + c2/s2tw*drho1)

      call setdelta(deltauv)
      call uucalc_renconst   ! for Z
      return
      end

      function getdrhoir(int)
      include 'shared.inc'
* int = 1 --> 1     loop
* int = 2 --> 1 + 2 loop
      xt = gf*mt*mt/sqrt(2.d0)/8.d0/pi/pi
      getdrhoir = 3.d0*xt
      if (int.eq.2) getdrhoir = getdrhoir-3.d0*(2.d0*pi*pi-19.d0)*xt*xt
      return
      end

      function getrhofi(q2,itype)
      include 'shared.inc'
      double precision iaqa
      complex*16 fdg
! itype = 1 --> up type quark
! itype = 2 --> down type quark
      csi = mh*mh/mz/mz
      c2  = 1.d0 - s2tw
      if (itype.eq.1) then
         iaqa = 5.d0/3.d0
         qaqa = 13.d0/9.d0
      else
         iaqa = 4.d0/3.d0
         qaqa = 10.d0/9.d0
      endif
      xw = q2/mw/mw
      xz = q2/mz/mz

      vfi = (1.d0-s2tw*(2.d0-iaqa))*fdg(xw) + 8.d0*c2*gdg(xw) 
     .     - fdg(xz)/2.d0/c2*(1.d0-3.d0*s2tw*iaqa+6.d0*s2tw*s2tw * qaqa)

      anum = 1.d0 + alpha/4.d0/pi/s2tw*
     .     (3.d0/4.d0/s2tw * log(c2) - 7.d0/4.d0 +
     .  3.d0/4.d0*csi*(log(c2/csi)/(c2-csi)+1.d0/c2*log(csi)/(1.d0-csi))
     .     )
      anum = anum + alpha/4.d0/pi/s2tw * vfi + AZZbosratio(q2)
      den = 1.d0 - AZZferratio(q2)

      getrhofi = anum/den
      return
      end

      function AZZbosratio(q2)
      include 'shared.inc'
      double complex AZZbos,dZZZbos,dZZZfer
      include 'libzvirtual/uubar/renconst.h'
      call getdZZZ(dZZZfer,dZZZbos)
      delta = 0.5d0
      if (abs(q2-mz*mz).lt.delta) then
         AZZbosratio = (AZZbos(0.d0)-AZZbos(mz*mz))/mz/mz + dZZZbos
      else
         AZZbosratio = AZZbos(q2)-q2/mz/mz*AZZbos(mz*mz) 
     .        + (q2-mz*mz)/mz/mz*AZZbos(0.d0)
         AZZbosratio = AZZbosratio/(q2-mz*mz)
      endif
      return
      end

      function AZZferratio(q2)
      include 'shared.inc'
      double complex AZZfer,AZZfer0
      double complex dZZZbos,dZZZfer
      include 'libzvirtual/uubar/renconst.h'
      call getdZZZ(dZZZfer,dZZZbos)
      delta = 0.5d0
      if (abs(q2-mz*mz).lt.delta) then
         AZZferratio = (AZZfer0(0.d0)-AZZfer(mz*mz))/mz/mz + dZZZfer
      else
         AZZferratio = AZZfer(q2)-q2/mz/mz*AZZfer(mz*mz) 
     .        + (q2-mz*mz)/mz/mz*AZZfer0(0.d0)
         AZZferratio = AZZferratio/(q2-mz*mz)
      endif
      return
      end

      subroutine getdZZZ(fer,bos)
      implicit character (a-s,u-z)
      implicit double complex (t)
      double complex sgzmzren, sggmzren, sggmz, sgz0, sgzmz,sgg0,szzmz
      double complex Agg, AgZ, AZZ, gzmz,fer,bos,dZZZ
      double precision dmztmp
      integer ite
      external Agg, AgZ, AZZ
#include "libzvirtual/uubar/vars.h"
      gzmz = dcmplx(0d0,gz*mz)
      
      dZZZ = (Alfa*CW2)/(6.D0*Pi*SW2)
!     -   +(5*Alfa*CW2)/(4.D0*Pi*SW2)*dble(B0(MZ2,MW2,MW2))
     -     +(3*Alfa)/(8.D0*CW2*Pi*SW2)*dble(B1(MZ2,0.D0,0.D0)) 
     -     +(Alfa*(9 - 12*SW2 + 8*SW2**2))/(24.D0*CW2*Pi*SW2)*
     -     dble(B1(MZ2,MB2,MB2)) + 
     -     (Alfa*(9 - 24*SW2 + 32*SW2**2))/(24.D0*CW2*Pi*SW2)*
     -     dble(B1(MZ2,MC2,MC2)) + 
     -     (Alfa*(9 - 12*SW2 + 8*SW2**2))/(24.D0*CW2*Pi*SW2)*
     -     dble(B1(MZ2,MD2,MD2)) + 
     -     (Alfa*(1 - 4*SW2 + 8*SW2**2))/(8.D0*CW2*Pi*SW2)*
     -     dble(B1(MZ2,ME2,ME2)) + 
     -     (Alfa*(1 - 4*SW2 + 8*SW2**2))/(8.D0*CW2*Pi*SW2)*
     -     dble(B1(MZ2,ML2,ML2)) + 
     -     (Alfa*(1 - 4*SW2 + 8*SW2**2))/(8.D0*CW2*Pi*SW2)*
     -     dble(B1(MZ2,MM2,MM2)) + 
     -     (Alfa*(9 - 12*SW2 + 8*SW2**2))/(24.D0*CW2*Pi*SW2)*
     -     dble(B1(MZ2,MS2,MS2)) + 
     -     (Alfa*(9 - 24*SW2 + 32*SW2**2))/(24.D0*CW2*Pi*SW2)*
     -     dble(B1(MZ2,MT2,MT2)) + 
     -     (Alfa*(9 - 24*SW2 + 32*SW2**2))/(24.D0*CW2*Pi*SW2)*
     -     dble(B1(MZ2,MU2,MU2))
!     -   +(Alfa*CW2)/(2.D0*Pi*SW2)*dble(B1(MZ2,MW2,MW2))
     -     +(3*Alfa*MB2)/(8.D0*CW2*Pi*SW2)*dble(DB0(MZ2,MB2,MB2)) + 
     -     (3*Alfa*MC2)/(8.D0*CW2*Pi*SW2)*dble(DB0(MZ2,MC2,MC2)) + 
     -     (3*Alfa*MD2)/(8.D0*CW2*Pi*SW2)*dble(DB0(MZ2,MD2,MD2)) + 
     -     (Alfa*ME2)/(8.D0*CW2*Pi*SW2)*dble(DB0(MZ2,ME2,ME2))
!     -   -(Alfa*MW2)/(4.D0*CW2**2*Pi*SW2)*dble(DB0(MZ2,MH2,MZ2)) 
     -     +(Alfa*ML2)/(8.D0*CW2*Pi*SW2)*dble(DB0(MZ2,ML2,ML2)) + 
     -     (Alfa*MM2)/(8.D0*CW2*Pi*SW2)*dble(DB0(MZ2,MM2,MM2)) + 
     -     (3*Alfa*MS2)/(8.D0*CW2*Pi*SW2)*dble(DB0(MZ2,MS2,MS2)) + 
     -     (3*Alfa*MT2)/(8.D0*CW2*Pi*SW2)*dble(DB0(MZ2,MT2,MT2)) + 
     -     (3*Alfa*MU2)/(8.D0*CW2*Pi*SW2)*dble(DB0(MZ2,MU2,MU2))
!     -   +(Alfa*(CW2**2*(2*MW2 + 5*MZ2) - 2*MW2*SW2**2))/
!     -     (4.D0*CW2*Pi*SW2)*dble(DB0(MZ2,MW2,MW2))
     -     -(3*Alfa)/(4.D0*CW2*Pi*SW2)*dble(DB00(MZ2,0.D0,0.D0)) - 
     -     (Alfa*(9 - 12*SW2 + 8*SW2**2))/(12.D0*CW2*Pi*SW2)*
     -     dble(DB00(MZ2,MB2,MB2)) - 
     -     (Alfa*(9 - 24*SW2 + 32*SW2**2))/(12.D0*CW2*Pi*SW2)*
     -     dble(DB00(MZ2,MC2,MC2)) - 
     -     (Alfa*(9 - 12*SW2 + 8*SW2**2))/(12.D0*CW2*Pi*SW2)*
     -     dble(DB00(MZ2,MD2,MD2))
      dZZZ = dZZZ - (Alfa*(1 - 4*SW2 + 8*SW2**2))/
     -     (4.D0*CW2*Pi*SW2)*dble(DB00(MZ2,ME2,ME2))
!     -   +Alfa/(4.D0*CW2*Pi*SW2)*dble(DB00(MZ2,MH2,MZ2)) 
     -     -(Alfa*(1 - 4*SW2 + 8*SW2**2))/(4.D0*CW2*Pi*SW2)*
     -     dble(DB00(MZ2,ML2,ML2)) - 
     -     (Alfa*(1 - 4*SW2 + 8*SW2**2))/(4.D0*CW2*Pi*SW2)*
     -     dble(DB00(MZ2,MM2,MM2)) - 
     -     (Alfa*(9 - 12*SW2 + 8*SW2**2))/(12.D0*CW2*Pi*SW2)*
     -     dble(DB00(MZ2,MS2,MS2)) - 
     -     (Alfa*(9 - 24*SW2 + 32*SW2**2))/(12.D0*CW2*Pi*SW2)*
     -     dble(DB00(MZ2,MT2,MT2)) - 
     -     (Alfa*(9 - 24*SW2 + 32*SW2**2))/(12.D0*CW2*Pi*SW2)*
     -     dble(DB00(MZ2,MU2,MU2))
!     -   +(Alfa*(9*CW2**2 - 2*CW2*SW2 + SW2**2))/(4.D0*CW2*Pi*SW2)*
!     -    dble(DB00(MZ2,MW2,MW2))
     -     +(3*Alfa*MZ2)/(8.D0*CW2*Pi*SW2)*dble(DB1(MZ2,0.D0,0.D0)) + 
     -     (Alfa*MZ2*(9 - 12*SW2 + 8*SW2**2))/(24.D0*CW2*Pi*SW2)*
     -     dble(DB1(MZ2,MB2,MB2)) + 
     -     (Alfa*MZ2*(9 - 24*SW2 + 32*SW2**2))/(24.D0*CW2*Pi*SW2)*
     -     dble(DB1(MZ2,MC2,MC2)) + 
     -     (Alfa*MZ2*(9 - 12*SW2 + 8*SW2**2))/(24.D0*CW2*Pi*SW2)*
     -     dble(DB1(MZ2,MD2,MD2)) + 
     -     (Alfa*MZ2*(1 - 4*SW2 + 8*SW2**2))/(8.D0*CW2*Pi*SW2)*
     -     dble(DB1(MZ2,ME2,ME2)) + 
     -     (Alfa*MZ2*(1 - 4*SW2 + 8*SW2**2))/(8.D0*CW2*Pi*SW2)*
     -     dble(DB1(MZ2,ML2,ML2)) + 
     -     (Alfa*MZ2*(1 - 4*SW2 + 8*SW2**2))/(8.D0*CW2*Pi*SW2)*
     -     dble(DB1(MZ2,MM2,MM2)) + 
     -     (Alfa*MZ2*(9 - 12*SW2 + 8*SW2**2))/(24.D0*CW2*Pi*SW2)*
     -     dble(DB1(MZ2,MS2,MS2)) + 
     -     (Alfa*MZ2*(9 - 24*SW2 + 32*SW2**2))/(24.D0*CW2*Pi*SW2)*
     -    dble(DB1(MZ2,MT2,MT2)) + 
     -     (Alfa*MZ2*(9 - 24*SW2 + 32*SW2**2))/(24.D0*CW2*Pi*SW2)*
     -     dble(DB1(MZ2,MU2,MU2))
!     -   +(Alfa*CW2*MZ2)/(2.D0*Pi*SW2)*dble(DB1(MZ2,MW2,MW2))
      
      fer = dZZZ
      dZZZ =
     -     +(5*Alfa*CW2)/(4.D0*Pi*SW2)*dble(B0(MZ2,MW2,MW2))
     -     +(Alfa*CW2)/(2.D0*Pi*SW2)*dble(B1(MZ2,MW2,MW2))
     -     -(Alfa*MW2)/(4.D0*CW2**2*Pi*SW2)*dble(DB0(MZ2,MH2,MZ2)) 
     -     +(Alfa*(CW2**2*(2*MW2 + 5*MZ2) - 2*MW2*SW2**2))/
     -     (4.D0*CW2*Pi*SW2)*dble(DB0(MZ2,MW2,MW2))
      dZZZ = dZZZ
     -     +Alfa/(4.D0*CW2*Pi*SW2)*dble(DB00(MZ2,MH2,MZ2)) 
     -     +(Alfa*(9*CW2**2 - 2*CW2*SW2 + SW2**2))/(4.D0*CW2*Pi*SW2)*
     -     dble(DB00(MZ2,MW2,MW2))
     -     +(Alfa*CW2*MZ2)/(2.D0*Pi*SW2)*dble(DB1(MZ2,MW2,MW2))
        bos = dZZZ
        return
	end







      function fdg(x)
      implicit double precision (a-h,o-z)
      double complex fdg,im
      im = (0.d0,1.d0)
      pi = 4.d0*atan(1.d0)
      fdg = im*pi*(2.d0/x + 3.d0 - 2.d0*(1.d0-1.d0/x)**2*log(1.d0+x)) +
     .     2.d0/x + 3.5d0 - (3.d0 + 2.d0/x)*log(x) +
     .     (1.d0+1.d0/x)**2 * (2.d0 * ddiloghere(1.d0/(1.d0 + x))
     .     - pi**2/3.d0 + log(1.d0+x)**2)
      return
      end

      function gdg(x)
      implicit double precision (a-h,o-z)

      gdg = (sqrt((4.d0-x)/x)*atan(sqrt(x/(4.d0-x))) - 1.d0  )* (1.d0/x + 0.5d0)
     .     + 9.d0/8.d0+1.d0/2.d0/x - (1.d0+1.d0/2.d0/x)*4.d0/x*(atan(sqrt(x/(4.d0-x))))**2 
      return
      end

      function getdeltar()
      include 'shared.inc'
      include 'libvirtual/renconst.h'
      double complex wselfnonresonant,cdeltar
      external wselfnonresonant
      double complex wselfatzero
      external wselfatzero
      double complex wselfp
      external wselfp
c      common /renconst/ dMWsq1, dMZsq1, dSW1, dZAA1, dZe1, dZW1
c      common /renconst/ dZZA1, dZfL1

      c2tw = 1.d0 - s2tw
      stw  = sqrt(s2tw)
      ctw  = sqrt(c2tw)
      qmw  = mw**2
      qmz  = mz**2

      cdeltar = -dzaa1 - c2tw/s2tw*(dmzsq1/qmz - dmwsq1/qmw) +
     .     (wselfatzero()-wselfnonresonant(qmw))/qmw + 
     .     ctw/stw*dzza1 +
     .     alpha/4.d0/pi/s2tw * (6.d0 + (7.d0-4.d0*s2tw)/2.d0/s2tw *
     .     log(c2tw))

      getdeltar = cdeltar
      return
      end

      double complex function Aggfer(SS)
      implicit character (a-s,u-z)
      implicit double complex (t)
      double precision SS
#include "libzvirtual/uubar/vars.h"

      Aggfer= (Alfa*SS)/(6.*Pi) + (Alfa*A0(MB2))/(3.*Pi) + (4*Alfa*A0(MC2))/(3.*Pi) + 
     -     (Alfa*A0(MD2))/(3.*Pi) + (Alfa*A0(ME2))/Pi + (Alfa*A0(ML2))/Pi + 
     -     (Alfa*A0(MM2))/Pi + (Alfa*A0(MS2))/(3.*Pi) + (4*Alfa*A0(MT2))/(3.*Pi) + 
     -     (4*Alfa*A0(MU2))/(3.*Pi)
     -     - (2*Alfa*B00(SS,MB2,MB2))/(3.*Pi) - 
     -     (8*Alfa*B00(SS,MC2,MC2))/(3.*Pi) - (2*Alfa*B00(SS,MD2,MD2))/(3.*Pi) - 
     -     (2*Alfa*B00(SS,ME2,ME2))/Pi - (2*Alfa*B00(SS,ML2,ML2))/Pi - 
     -     (2*Alfa*B00(SS,MM2,MM2))/Pi - (2*Alfa*B00(SS,MS2,MS2))/(3.*Pi) - 
     -     (8*Alfa*B00(SS,MT2,MT2))/(3.*Pi) - (8*Alfa*B00(SS,MU2,MU2))/(3.*Pi) + 
     -     (Alfa*SS*B1(SS,MB2,MB2))/(3.*Pi) + 
     -     (4*Alfa*SS*B1(SS,MC2,MC2))/(3.*Pi) + (Alfa*SS*B1(SS,MD2,MD2))/(3.*Pi) + 
     -     (Alfa*SS*B1(SS,ME2,ME2))/Pi + (Alfa*SS*B1(SS,ML2,ML2))/Pi + 
     -     (Alfa*SS*B1(SS,MM2,MM2))/Pi + (Alfa*SS*B1(SS,MS2,MS2))/(3.*Pi) + 
     -     (4*Alfa*SS*B1(SS,MT2,MT2))/(3.*Pi) + 
     -     (4*Alfa*SS*B1(SS,MU2,MU2))/(3.*Pi)
      end function Aggfer

        double complex function AZZfer(SS)
	implicit character (a-s,u-z)
	implicit double complex (t)
        double complex AZZ
        double precision SS
#include "libzvirtual/uubar/vars.h"
        AZZ= (Alfa*CW2*SS)/(6.*Pi*SW2) - (Alfa*A0(MB2))/(2.*CW2*Pi) + 
     -  (3*Alfa*A0(MB2))/(8.*CW2*Pi*SW2) + (Alfa*SW2*A0(MB2))/(3.*CW2*Pi) - 
     -  (Alfa*A0(MC2))/(CW2*Pi) + (3*Alfa*A0(MC2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SW2*A0(MC2))/(3.*CW2*Pi) - (Alfa*A0(MD2))/(2.*CW2*Pi) + 
     -  (3*Alfa*A0(MD2))/(8.*CW2*Pi*SW2) + (Alfa*SW2*A0(MD2))/(3.*CW2*Pi) - 
     -  (Alfa*A0(ME2))/(2.*CW2*Pi) + (Alfa*A0(ME2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SW2*A0(ME2))/(CW2*Pi) 
!- (Alfa*A0(MH2))/(16.*CW2*Pi*SW2)
     -  -(Alfa*A0(ML2))/(2.*CW2*Pi) + (Alfa*A0(ML2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SW2*A0(ML2))/(CW2*Pi) - (Alfa*A0(MM2))/(2.*CW2*Pi) + 
     -  (Alfa*A0(MM2))/(8.*CW2*Pi*SW2) + (Alfa*SW2*A0(MM2))/(CW2*Pi) - 
     -  (Alfa*A0(MS2))/(2.*CW2*Pi) + (3*Alfa*A0(MS2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SW2*A0(MS2))/(3.*CW2*Pi) - (Alfa*A0(MT2))/(CW2*Pi) + 
     -  (3*Alfa*A0(MT2))/(8.*CW2*Pi*SW2) + (4*Alfa*SW2*A0(MT2))/(3.*CW2*Pi) - 
     -  (Alfa*A0(MU2))/(CW2*Pi) + (3*Alfa*A0(MU2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SW2*A0(MU2))/(3.*CW2*Pi) 
!+ (Alfa*A0(MW2))/(4.*Pi)
!     -  -(9*Alfa*CW2*A0(MW2))/(8.*Pi*SW2)
! - (Alfa*SW2*A0(MW2))/(8.*CW2*Pi) 
!     -  -(Alfa*A0(MZ2))/(16.*CW2*Pi*SW2)
     -  +(3*Alfa*MB2*B0(SS,MB2,MB2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MC2*B0(SS,MC2,MC2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MD2*B0(SS,MD2,MD2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*ME2*B0(SS,ME2,ME2))/(8.*CW2*Pi*SW2) 
!     -  -(Alfa*MW2*B0(SS,MH2,MZ2))/(4.*CW2**2*Pi*SW2)
     -  +(Alfa*ML2*B0(SS,ML2,ML2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*MM2*B0(SS,MM2,MM2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MS2*B0(SS,MS2,MS2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MT2*B0(SS,MT2,MT2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MU2*B0(SS,MU2,MU2))/(8.*CW2*Pi*SW2)
!     -  +(Alfa*CW2*MW2*B0(SS,MW2,MW2))/(2.*Pi*SW2) 

        AZZ = AZZ
!     -  +(5*Alfa*CW2*SS*B0(SS,MW2,MW2))/(4.*Pi*SW2)
!     -  -(Alfa*MW2*SW2*B0(SS,MW2,MW2))/(2.*CW2*Pi)
     -  -(3*Alfa*B00(SS,0D0,0D0))/(4.*CW2*Pi*SW2) +
     -  (Alfa*B00(SS,MB2,MB2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MB2,MB2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,MB2,MB2))/(3.*CW2*Pi) + 
     -  (2*Alfa*B00(SS,MC2,MC2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MC2,MC2))/(4.*CW2*Pi*SW2) - 
     -  (8*Alfa*SW2*B00(SS,MC2,MC2))/(3.*CW2*Pi) + 
     -  (Alfa*B00(SS,MD2,MD2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MD2,MD2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,MD2,MD2))/(3.*CW2*Pi) + 
     -  (Alfa*B00(SS,ME2,ME2))/(CW2*Pi) - 
     -  (Alfa*B00(SS,ME2,ME2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,ME2,ME2))/(CW2*Pi)
!     -  +(Alfa*B00(SS,MH2,MZ2))/(4.*CW2*Pi*SW2)
     -  +(Alfa*B00(SS,ML2,ML2))/(CW2*Pi) - 
     -  (Alfa*B00(SS,ML2,ML2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,ML2,ML2))/(CW2*Pi) + 
     -  (Alfa*B00(SS,MM2,MM2))/(CW2*Pi) - 
     -  (Alfa*B00(SS,MM2,MM2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,MM2,MM2))/(CW2*Pi) + 
     -  (Alfa*B00(SS,MS2,MS2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MS2,MS2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,MS2,MS2))/(3.*CW2*Pi) + 
     -  (2*Alfa*B00(SS,MT2,MT2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MT2,MT2))/(4.*CW2*Pi*SW2) - 
     -  (8*Alfa*SW2*B00(SS,MT2,MT2))/(3.*CW2*Pi) + 
     -  (2*Alfa*B00(SS,MU2,MU2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MU2,MU2))/(4.*CW2*Pi*SW2) - 
     -  (8*Alfa*SW2*B00(SS,MU2,MU2))/(3.*CW2*Pi)
!     -  -(Alfa*B00(SS,MW2,MW2))/(2.*Pi) + 
!     -  (9*Alfa*CW2*B00(SS,MW2,MW2))/(4.*Pi*SW2) + 
!     -  (Alfa*SW2*B00(SS,MW2,MW2))/(4.*CW2*Pi)
     -  +(3*Alfa*SS*B1(SS,0D0,0D0))/(8.*CW2*Pi*SW2) - 
     -  (Alfa*SS*B1(SS,MB2,MB2))/(2.*CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MB2,MB2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,MB2,MB2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MC2,MC2))/(CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MC2,MC2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SS*SW2*B1(SS,MC2,MC2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MD2,MD2))/(2.*CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MD2,MD2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,MD2,MD2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,ME2,ME2))/(2.*CW2*Pi) + 
     -  (Alfa*SS*B1(SS,ME2,ME2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,ME2,ME2))/(CW2*Pi) - 
     -  (Alfa*SS*B1(SS,ML2,ML2))/(2.*CW2*Pi) + 
     -  (Alfa*SS*B1(SS,ML2,ML2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,ML2,ML2))/(CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MM2,MM2))/(2.*CW2*Pi) + 
     -  (Alfa*SS*B1(SS,MM2,MM2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,MM2,MM2))/(CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MS2,MS2))/(2.*CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MS2,MS2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,MS2,MS2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MT2,MT2))/(CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MT2,MT2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SS*SW2*B1(SS,MT2,MT2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MU2,MU2))/(CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MU2,MU2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SS*SW2*B1(SS,MU2,MU2))/(3.*CW2*Pi)
!     -  +(Alfa*CW2*SS*B1(SS,MW2,MW2))/(2.*Pi*SW2)

        AZZfer = AZZ

        end function AZZfer



        double complex function AZZfer0(SS)
	implicit character (a-s,u-z)
	implicit double complex (t)
        double complex AZZ
        double precision SS
#include "libzvirtual/uubar/vars.h"
        AZZ= (Alfa*CW2*SS)/(6.*Pi*SW2) - (Alfa*A0(MB2))/(2.*CW2*Pi) + 
     -  (3*Alfa*A0(MB2))/(8.*CW2*Pi*SW2) + (Alfa*SW2*A0(MB2))/(3.*CW2*Pi) - 
     -  (Alfa*A0(MC2))/(CW2*Pi) + (3*Alfa*A0(MC2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SW2*A0(MC2))/(3.*CW2*Pi) - (Alfa*A0(MD2))/(2.*CW2*Pi) + 
     -  (3*Alfa*A0(MD2))/(8.*CW2*Pi*SW2) + (Alfa*SW2*A0(MD2))/(3.*CW2*Pi) - 
     -  (Alfa*A0(ME2))/(2.*CW2*Pi) + (Alfa*A0(ME2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SW2*A0(ME2))/(CW2*Pi) 
!- (Alfa*A0(MH2))/(16.*CW2*Pi*SW2)
     -  -(Alfa*A0(ML2))/(2.*CW2*Pi) + (Alfa*A0(ML2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SW2*A0(ML2))/(CW2*Pi) - (Alfa*A0(MM2))/(2.*CW2*Pi) + 
     -  (Alfa*A0(MM2))/(8.*CW2*Pi*SW2) + (Alfa*SW2*A0(MM2))/(CW2*Pi) - 
     -  (Alfa*A0(MS2))/(2.*CW2*Pi) + (3*Alfa*A0(MS2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SW2*A0(MS2))/(3.*CW2*Pi) - (Alfa*A0(MT2))/(CW2*Pi) + 
     -  (3*Alfa*A0(MT2))/(8.*CW2*Pi*SW2) + (4*Alfa*SW2*A0(MT2))/(3.*CW2*Pi) - 
     -  (Alfa*A0(MU2))/(CW2*Pi) + (3*Alfa*A0(MU2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SW2*A0(MU2))/(3.*CW2*Pi) 
!+ (Alfa*A0(MW2))/(4.*Pi)
!     -  -(9*Alfa*CW2*A0(MW2))/(8.*Pi*SW2)
! - (Alfa*SW2*A0(MW2))/(8.*CW2*Pi) 
!     -  -(Alfa*A0(MZ2))/(16.*CW2*Pi*SW2)
     -  +(3*Alfa*MB2*B0(SS,MB2,MB2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MC2*B0(SS,MC2,MC2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MD2*B0(SS,MD2,MD2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*ME2*B0(SS,ME2,ME2))/(8.*CW2*Pi*SW2) 
!     -  -(Alfa*MW2*B0(SS,MH2,MZ2))/(4.*CW2**2*Pi*SW2)
     -  +(Alfa*ML2*B0(SS,ML2,ML2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*MM2*B0(SS,MM2,MM2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MS2*B0(SS,MS2,MS2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MT2*B0(SS,MT2,MT2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MU2*B0(SS,MU2,MU2))/(8.*CW2*Pi*SW2)
!     -  +(Alfa*CW2*MW2*B0(SS,MW2,MW2))/(2.*Pi*SW2) 

        AZZ = AZZ
!     -  +(5*Alfa*CW2*SS*B0(SS,MW2,MW2))/(4.*Pi*SW2)
!     -  -(Alfa*MW2*SW2*B0(SS,MW2,MW2))/(2.*CW2*Pi)
     -  +(Alfa*B00(SS,MB2,MB2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MB2,MB2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,MB2,MB2))/(3.*CW2*Pi) + 
     -  (2*Alfa*B00(SS,MC2,MC2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MC2,MC2))/(4.*CW2*Pi*SW2) - 
     -  (8*Alfa*SW2*B00(SS,MC2,MC2))/(3.*CW2*Pi) + 
     -  (Alfa*B00(SS,MD2,MD2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MD2,MD2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,MD2,MD2))/(3.*CW2*Pi) + 
     -  (Alfa*B00(SS,ME2,ME2))/(CW2*Pi) - 
     -  (Alfa*B00(SS,ME2,ME2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,ME2,ME2))/(CW2*Pi)
!     -  +(Alfa*B00(SS,MH2,MZ2))/(4.*CW2*Pi*SW2)
     -  +(Alfa*B00(SS,ML2,ML2))/(CW2*Pi) - 
     -  (Alfa*B00(SS,ML2,ML2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,ML2,ML2))/(CW2*Pi) + 
     -  (Alfa*B00(SS,MM2,MM2))/(CW2*Pi) - 
     -  (Alfa*B00(SS,MM2,MM2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,MM2,MM2))/(CW2*Pi) + 
     -  (Alfa*B00(SS,MS2,MS2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MS2,MS2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,MS2,MS2))/(3.*CW2*Pi) + 
     -  (2*Alfa*B00(SS,MT2,MT2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MT2,MT2))/(4.*CW2*Pi*SW2) - 
     -  (8*Alfa*SW2*B00(SS,MT2,MT2))/(3.*CW2*Pi) + 
     -  (2*Alfa*B00(SS,MU2,MU2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MU2,MU2))/(4.*CW2*Pi*SW2) - 
     -  (8*Alfa*SW2*B00(SS,MU2,MU2))/(3.*CW2*Pi)
!     -  -(Alfa*B00(SS,MW2,MW2))/(2.*Pi) + 
!     -  (9*Alfa*CW2*B00(SS,MW2,MW2))/(4.*Pi*SW2) + 
!     -  (Alfa*SW2*B00(SS,MW2,MW2))/(4.*CW2*Pi)
     -  -(Alfa*SS*B1(SS,MB2,MB2))/(2.*CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MB2,MB2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,MB2,MB2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MC2,MC2))/(CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MC2,MC2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SS*SW2*B1(SS,MC2,MC2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MD2,MD2))/(2.*CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MD2,MD2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,MD2,MD2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,ME2,ME2))/(2.*CW2*Pi) + 
     -  (Alfa*SS*B1(SS,ME2,ME2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,ME2,ME2))/(CW2*Pi) - 
     -  (Alfa*SS*B1(SS,ML2,ML2))/(2.*CW2*Pi) + 
     -  (Alfa*SS*B1(SS,ML2,ML2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,ML2,ML2))/(CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MM2,MM2))/(2.*CW2*Pi) + 
     -  (Alfa*SS*B1(SS,MM2,MM2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,MM2,MM2))/(CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MS2,MS2))/(2.*CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MS2,MS2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,MS2,MS2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MT2,MT2))/(CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MT2,MT2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SS*SW2*B1(SS,MT2,MT2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MU2,MU2))/(CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MU2,MU2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SS*SW2*B1(SS,MU2,MU2))/(3.*CW2*Pi)
!     -  +(Alfa*CW2*SS*B1(SS,MW2,MW2))/(2.*Pi*SW2)

        AZZfer0 = AZZ

        end function AZZfer0





        double complex function AZZbos(SS)
	implicit character (a-s,u-z)
	implicit double complex (t)
        double complex AZZ
        double precision SS
#include "libzvirtual/uubar/vars.h"
        AZZ= 
     .       -(Alfa*A0(MH2))/(16.*CW2*Pi*SW2)
     .       + (Alfa*A0(MW2))/(4.*Pi)
     .       -(9*Alfa*CW2*A0(MW2))/(8.*Pi*SW2)
     .       - (Alfa*SW2*A0(MW2))/(8.*CW2*Pi) 
     ,       -(Alfa*A0(MZ2))/(16.*CW2*Pi*SW2)
     ,       -(Alfa*MW2*B0(SS,MH2,MZ2))/(4.*CW2**2*Pi*SW2)
     ,       +(Alfa*CW2*MW2*B0(SS,MW2,MW2))/(2.*Pi*SW2) 

        AZZ = AZZ
     ,       +(5*Alfa*CW2*SS*B0(SS,MW2,MW2))/(4.*Pi*SW2)
     ,       -(Alfa*MW2*SW2*B0(SS,MW2,MW2))/(2.*CW2*Pi)
     ,       +(Alfa*B00(SS,MH2,MZ2))/(4.*CW2*Pi*SW2)
     ,       -(Alfa*B00(SS,MW2,MW2))/(2.*Pi) + 
     ,       (9*Alfa*CW2*B00(SS,MW2,MW2))/(4.*Pi*SW2) + 
     ,       (Alfa*SW2*B00(SS,MW2,MW2))/(4.*CW2*Pi)
     ,       +(Alfa*CW2*SS*B1(SS,MW2,MW2))/(2.*Pi*SW2)

        AZZbos = AZZ

        end function AZZbos



        double complex function AgZ(SS)
	implicit character (a-s,u-z)
	implicit double complex (t)
        double precision SS
#include "libzvirtual/uubar/vars.h"
      AgZ= -(Alfa*CW*SS)/(6.*Pi*SW) - (Alfa*A0(MB2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SW*A0(MB2))/(6.*CW*Pi) + (Alfa*SW2*A0(MB2))/(6.*CW*Pi*SW) - 
     -  (Alfa*A0(MC2))/(2.*CW*Pi*SW) + (2*Alfa*SW*A0(MC2))/(3.*CW*Pi) + 
     -  (2*Alfa*SW2*A0(MC2))/(3.*CW*Pi*SW) - (Alfa*A0(MD2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SW*A0(MD2))/(6.*CW*Pi) + (Alfa*SW2*A0(MD2))/(6.*CW*Pi*SW) - 
     -  (Alfa*A0(ME2))/(4.*CW*Pi*SW) + (Alfa*SW*A0(ME2))/(2.*CW*Pi) + 
     -  (Alfa*SW2*A0(ME2))/(2.*CW*Pi*SW) - (Alfa*A0(ML2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SW*A0(ML2))/(2.*CW*Pi) + (Alfa*SW2*A0(ML2))/(2.*CW*Pi*SW) - 
     -  (Alfa*A0(MM2))/(4.*CW*Pi*SW) + (Alfa*SW*A0(MM2))/(2.*CW*Pi) + 
     -  (Alfa*SW2*A0(MM2))/(2.*CW*Pi*SW) - (Alfa*A0(MS2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SW*A0(MS2))/(6.*CW*Pi) + (Alfa*SW2*A0(MS2))/(6.*CW*Pi*SW) - 
     -  (Alfa*A0(MT2))/(2.*CW*Pi*SW) + (2*Alfa*SW*A0(MT2))/(3.*CW*Pi) + 
     -  (2*Alfa*SW2*A0(MT2))/(3.*CW*Pi*SW) - (Alfa*A0(MU2))/(2.*CW*Pi*SW) + 
     -  (2*Alfa*SW*A0(MU2))/(3.*CW*Pi) + (2*Alfa*SW2*A0(MU2))/(3.*CW*Pi*SW) + 
     -  (Alfa*CW*A0(MW2))/(Pi*SW) + (Alfa*CW2*A0(MW2))/(4.*CW*Pi*SW) - 
     -  (Alfa*SW2*A0(MW2))/(4.*CW*Pi*SW) - 
     -  (Alfa*CW*MW2*B0(SS,MW2,MW2))/(2.*Pi*SW) - 
     -  (5*Alfa*CW*SS*B0(SS,MW2,MW2))/(4.*Pi*SW) - 
     -  (Alfa*MW2*SW*B0(SS,MW2,MW2))/(2.*CW*Pi) + 
     -  (Alfa*B00(SS,MB2,MB2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SW*B00(SS,MB2,MB2))/(3.*CW*Pi) - 
     -  (Alfa*SW2*B00(SS,MB2,MB2))/(3.*CW*Pi*SW) + 
     -  (Alfa*B00(SS,MC2,MC2))/(CW*Pi*SW) - 
     -  (4*Alfa*SW*B00(SS,MC2,MC2))/(3.*CW*Pi) - 
     -  (4*Alfa*SW2*B00(SS,MC2,MC2))/(3.*CW*Pi*SW) + 
     -  (Alfa*B00(SS,MD2,MD2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SW*B00(SS,MD2,MD2))/(3.*CW*Pi) - 
     -  (Alfa*SW2*B00(SS,MD2,MD2))/(3.*CW*Pi*SW) + 
     -  (Alfa*B00(SS,ME2,ME2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SW*B00(SS,ME2,ME2))/(CW*Pi) - 
     -  (Alfa*SW2*B00(SS,ME2,ME2))/(CW*Pi*SW) + 
     -  (Alfa*B00(SS,ML2,ML2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SW*B00(SS,ML2,ML2))/(CW*Pi) - 
     -  (Alfa*SW2*B00(SS,ML2,ML2))/(CW*Pi*SW) + 
     -  (Alfa*B00(SS,MM2,MM2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SW*B00(SS,MM2,MM2))/(CW*Pi) - 
     -  (Alfa*SW2*B00(SS,MM2,MM2))/(CW*Pi*SW) + 
     -  (Alfa*B00(SS,MS2,MS2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SW*B00(SS,MS2,MS2))/(3.*CW*Pi) - 
     -  (Alfa*SW2*B00(SS,MS2,MS2))/(3.*CW*Pi*SW) + 
     -  (Alfa*B00(SS,MT2,MT2))/(CW*Pi*SW) - 
     -  (4*Alfa*SW*B00(SS,MT2,MT2))/(3.*CW*Pi) - 
     -  (4*Alfa*SW2*B00(SS,MT2,MT2))/(3.*CW*Pi*SW) + 
     -  (Alfa*B00(SS,MU2,MU2))/(CW*Pi*SW) - 
     -  (4*Alfa*SW*B00(SS,MU2,MU2))/(3.*CW*Pi) - 
     -  (4*Alfa*SW2*B00(SS,MU2,MU2))/(3.*CW*Pi*SW) - 
     -  (2*Alfa*CW*B00(SS,MW2,MW2))/(Pi*SW) - 
     -  (Alfa*CW2*B00(SS,MW2,MW2))/(2.*CW*Pi*SW) + 
     -  (Alfa*SW2*B00(SS,MW2,MW2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MB2,MB2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SS*SW*B1(SS,MB2,MB2))/(6.*CW*Pi) + 
     -  (Alfa*SS*SW2*B1(SS,MB2,MB2))/(6.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MC2,MC2))/(2.*CW*Pi*SW) + 
     -  (2*Alfa*SS*SW*B1(SS,MC2,MC2))/(3.*CW*Pi) + 
     -  (2*Alfa*SS*SW2*B1(SS,MC2,MC2))/(3.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MD2,MD2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SS*SW*B1(SS,MD2,MD2))/(6.*CW*Pi) + 
     -  (Alfa*SS*SW2*B1(SS,MD2,MD2))/(6.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,ME2,ME2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SS*SW*B1(SS,ME2,ME2))/(2.*CW*Pi) + 
     -  (Alfa*SS*SW2*B1(SS,ME2,ME2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,ML2,ML2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SS*SW*B1(SS,ML2,ML2))/(2.*CW*Pi) + 
     -  (Alfa*SS*SW2*B1(SS,ML2,ML2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MM2,MM2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SS*SW*B1(SS,MM2,MM2))/(2.*CW*Pi) + 
     -  (Alfa*SS*SW2*B1(SS,MM2,MM2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MS2,MS2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SS*SW*B1(SS,MS2,MS2))/(6.*CW*Pi) + 
     -  (Alfa*SS*SW2*B1(SS,MS2,MS2))/(6.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MT2,MT2))/(2.*CW*Pi*SW) + 
     -  (2*Alfa*SS*SW*B1(SS,MT2,MT2))/(3.*CW*Pi) + 
     -  (2*Alfa*SS*SW2*B1(SS,MT2,MT2))/(3.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MU2,MU2))/(2.*CW*Pi*SW) + 
     -  (2*Alfa*SS*SW*B1(SS,MU2,MU2))/(3.*CW*Pi) + 
     -  (2*Alfa*SS*SW2*B1(SS,MU2,MU2))/(3.*CW*Pi*SW) - 
     -  (Alfa*CW*SS*B1(SS,MW2,MW2))/(2.*Pi*SW)

      end function AgZ














      function ctsnew()
      include 'shared.inc'
      double precision deltar
      double complex ctsnew,wselfp
      double complex dze1new,dswnew,dzwnew,dzf1,dzf2,dzf3,dzf4,dmw2
      external wselfp,getdelta
      common/newct/dze1new,dswnew,dzwnew,dzf1,dzf2,dzf3,dzf4,dmw2
      common/deltaerre/deltar

      denom = horaces - mw**2

      ctsnew =0.d0 * dmw2/denom + ! if included in dww
     >        2.d0*dze1new -
     >        2.d0*dswnew/sqrt(s2tw) +
     >        0.5d0*dzf1 +  
     >        0.5d0*dzf2 +  
     >        0.5d0*dzf3 +  
     >        0.5d0*dzf4
     >        -deltar !....+ 0.00435918367346938775d0/2.d0
c     >        - wselfp(mw*mw)
c      print*,dze1new,dswnew/sqrt(s2tw),dzf1,dzf2,dzf3,dzf4

c      print*,dze1new
c      stop

c      print*,2.d0*dze1new-deltaalpha(mz**2)
c      stop
      return
      end
******************************
      double complex function wselfp(s)
      include 'shared.inc'
      double precision mlep(3)
      double precision lm,mw2,mh2,mz2,mg2
      double precision mquu(3),mqud(3)
      double complex B0,wlep,wbos,wbos1,wbos2,wbos3,wquark,dB0
      external B0,dB0

      mg2 = getlambda()

      s2 = s2tw
      c2 = 1.d0 - s2

      mw2 = mw*mw
      mz2 = mz*mz
      mh2 = mh*mh

      mlep(1) = me**2
      mlep(2) = mmu**2
      mlep(3) = mtau**2

      mquu(1) = mu**2
      mquu(2) = mc**2
      mquu(3) = mt**2 

      mqud(1) = md**2 
      mqud(2) = ms**2
      mqud(3) = mb**2 

      wlep   = 0.d0
      wquark = 0.d0
      wbos   = 0.d0
      wbos1   = 0.d0
      wbos2   = 0.d0
      wbos3   = 0.d0

      do k = 1,3
         lm = mlep(k)
         wlep = wlep
     >        -B0(s,0.d0,lm) - s*dB0(s,0.d0,lm) +lm/2.d0*dB0(s,0.d0,lm)

         wlep = wlep + 1.d0/3.d0
     >        -lm**2/2.d0/s/s* (B0(s,0.d0,lm) - B0(0.d0,0.d0,lm))

        wlep = wlep +lm**2/2.d0/s*dB0(s,0.d0,lm)
      enddo
      wlep = wlep * 2.d0/3.d0/2.d0/s2
      
      do k = 1,3
         um = mquu(k)
         dm = mqud(k)
         wquark = wquark
     >        -B0(s,um,dm)
     >        -(s-(um+dm)/2.d0)*dB0(s,um,dm) 
     >        + 1.d0/3.d0 
     >        -(um - dm)**2/2.d0/s/s* (B0(s,um,dm) - B0(0.d0,um,dm))
     >        +(um - dm)**2/2.d0/s*dB0(s,um,dm)
      enddo
      wquark = wquark * 2.d0/3.d0/2.d0/s2 * 3.d0

      wbos1 = (2.d0*mw2+5.d0*s)*dB0(s,mw2,mg2) + 5.d0*B0(s,mw2,mg2)
c   wbos1 = (2.d0*mw2+5.d0*s)*dB0(s,mw2,mg2)+5.d0*B0resc(s,mg2,mw2,gw)
     >        + 1.d0/3.d0
     >        + mw2**2/s/s*(B0(s,mw2,mg2)-B0(0.d0,mw2,mg2))
c     >        + mw2**2/s/s*(B0resc(s,mg2,mw2,gw)-B0(0.d0,mw2,mg2))
     >        - mw2**2/s*dB0(s,mw2,mg2)
      wbos1 = wbos1*2.d0/3.d0

      wbos2 = ((40.d0*c2-1.d0)*s+(16.d0*c2+54.d0-10.d0/c2)*mw2)
     >        *dB0(s,mw2,mz2)+
     >        (40.d0*c2-1.d0)*B0(s,mw2,mz2)
     >        + (4.d0*c2 - 1.d0) * 2.d0/3.d0
     >        + (8.d0*c2 + 1.d0) * (mw2-mz2)**2/s/s
     >          *(B0(s,mw2,mz2)-B0(0.d0,mw2,mz2))
     >        - (8.d0*c2 + 1.d0) * (mw2-mz2)**2/s
     >          *dB0(s,mw2,mz2)

      wbos2 = wbos2/12.d0/s2

      wbos3 = -B0(s,mw2,mh2) + 
     >        (2.d0*mh2 - 10.d0*mw2 - s)*dB0(s,mw2,mh2) +
     >        (mw2-mh2)**2/s/s*(B0(s,mw2,mh2)-B0(0.d0,mw2,mh2)) -
     >        (mw2-mh2)**2/s*dB0(s,mw2,mh2) 
     >         -2.d0/3.d0

      wbos3 = wbos3/12.d0/s2
      
      wbos = wbos1 + wbos2 + wbos3

      wlep   = -alpha/4.d0/pi * wlep
      wbos   = -alpha/4.d0/pi * wbos
      wquark = -alpha/4.d0/pi * wquark

      wselfp = wlep + wbos + wquark
      return
      end
**
      double complex function Pair(i, j)
      implicit none
      integer i, j
#include "2to2.h"
      Pair =
     +     vec(0, i)*vec(0, j) -
     +     vec(1, i)*vec(1, j) -
     +     vec(2, i)*vec(2, j) -
     +     vec(3, i)*vec(3, j)
      end
*     * Eps calculates -I times the Levi-Civita tensor contracted with
*     * four four-vectors.
      double complex function Eps(i, j, k, l)
      implicit none
      integer i, j, k, l
#include "2to2.h"
      double complex x01, x02, x03, x12, x13, x23
      double complex Cross
      integer mu_, nu_
      Cross(mu_, nu_) =
     +     vec(mu_, k)*vec(nu_, l) - vec(nu_, k)*vec(mu_, l)
      x01 = Cross(0, 1)
      x02 = Cross(0, 2)
      x03 = Cross(0, 3)
      x12 = Cross(1, 2)
      x13 = Cross(1, 3)
      x23 = Cross(2, 3)
      Eps =
     +     vec(0, i)*(vec(1, j)*x23 - vec(2, j)*x13 + vec(3, j)*x12) -
     +     vec(1, i)*(vec(0, j)*x23 - vec(2, j)*x03 + vec(3, j)*x02) +
     +     vec(2, i)*(vec(0, j)*x13 - vec(1, j)*x03 + vec(3, j)*x01) -
     +     vec(3, i)*(vec(0, j)*x12 - vec(1, j)*x02 + vec(2, j)*x01)
*     Note: although Eps is defined as -I*(Levi-Civita tensor) in FormCalc,
*     we have to multiply with I (not -I) here because we've neglected
*     the metric in the summation above.
      Eps = cI*Eps
      end

      double complex function wselfatzero()
      include 'shared.inc'
      double precision mlep(3)
      double precision lm,mw2,mh2,mz2,mg2
      double precision mquu(3),mqud(3)
      double complex B0,wlep,wbos,wbos1,wbos2,wbos3,wquark,dB0
      external B0,dB0
      s = 0.d0
      mg2 = getlambda()
      s2 = s2tw
      c2 = 1.d0 - s2
      mw2 = mw*mw
      mz2 = mz*mz
      mh2 = mh*mh
      mlep(1) = me**2
      mlep(2) = mmu**2
      mlep(3) = mtau**2
      mquu(1) = mu**2
      mquu(2) = mc**2
      mquu(3) = mt**2 

      mqud(1) = md**2 
      mqud(2) = ms**2
      mqud(3) = mb**2 

      wlep    = 0.d0
      wquark  = 0.d0
      wbos    = 0.d0
      wbos1   = 0.d0
      wbos2   = 0.d0
      wbos3   = 0.d0

      do k = 1,3
         lm = mlep(k)
         wlep = wlep 
     >        -(s - lm/2.d0)*B0(0.d0,0.d0,lm) +
     >        lm*B0(0d0,lm,lm) + 
     >        lm**2/2.d0*dB0(0.d0,0.d0,lm)
      enddo

      wlep = wlep * 2.d0/3.d0/2.d0/s2

      do k = 1,3
         um = mquu(k)
         dm = mqud(k)
         wquark = wquark
     >        -(s - (um+dm)/2.d0)*B0(0.d0,um,dm) +
     >        um*B0(0d0,um,um) + 
     >        dm*B0(0d0,dm,dm) + 
     >        (um - dm)**2/2.d0* (dB0(0.d0,um,dm))
      enddo

      wquark = wquark * 2.d0/3.d0/2.d0/s2 * 3.d0

      wbos1 = (2.d0*mw2)*B0(0.d0,mw2,mg2)-2.d0*mw2*B0(0.d0,mw2,mw2)
     >        - mw2**2*(dB0(0.d0,mw2,mg2))

      wbos1 = wbos1*2.d0/3.d0

      wbos2 = ((16.d0*c2+54.d0-10.d0/c2)*mw2)
     >        *B0(0.d0,mw2,mz2)
     >        - (16.d0*c2 + 2.d0)
     >          *(mw2*B0(0.d0,mw2,mw2)+mz2*B0(0.d0,mz2,mz2)) 
     >        + (4.d0*c2 - 1.d0) * 2.d0/3.d0*0.d0
     >        - (8.d0*c2 + 1.d0) * (mw2-mz2)**2
     >          *(dB0(0.d0,mw2,mz2))
      wbos2 = wbos2/12.d0/s2

      wbos3 = (2.d0*mh2 - 10.d0*mw2)*B0(0.d0,mw2,mh2) -
     >        2.d0*mw2 * B0(0.d0,mw2,mw2) - 2.d0*mh2*B0(0.d0,mh2,mh2) -
     >        (mw2-mh2)**2*(dB0(0.d0,mw2,mh2)) 

      wbos3 = wbos3/12.d0/s2
      
      wbos = wbos1 + wbos2  + wbos3

      wlep   = -alpha/4.d0/pi * wlep
      wbos   = -alpha/4.d0/pi * wbos
      wquark = -alpha/4.d0/pi * wquark

      wselfatzero = wlep + wbos + wquark
      return
      end

      double complex function wselfnonresonant(s)
      include 'shared.inc'
      double precision mlep(3)
      double precision lm,mw2,mh2,mz2,mg2
      double precision mquu(3),mqud(3)
      double complex B0,wlep,wbos,wbos1,wbos2,wbos3,wquark
      external B0

      mg2 = getlambda()

      s2 = s2tw
      c2 = 1.d0 - s2

      mw2 = mw*mw
      mz2 = mz*mz
      mh2 = mh*mh

      mlep(1) = me**2
      mlep(2) = mmu**2
      mlep(3) = mtau**2

      mquu(1) = mu**2
      mquu(2) = mc**2
      mquu(3) = mt**2 

      mqud(1) = md**2 
      mqud(2) = ms**2
      mqud(3) = mb**2 

      wlep    = 0.d0
      wquark  = 0.d0
      wbos    = 0.d0
      wbos1   = 0.d0
      wbos2   = 0.d0
      wbos3   = 0.d0

      do k = 1,3
         lm = mlep(k)
         wlep = wlep
     >        -(s - lm/2.d0)*B0(s,0.d0,lm) + s/3.d0 +
     >        lm*B0(0d0,lm,lm) + 
     >        lm**2/2.d0/s* (B0(s,0.d0,lm) - B0(0.d0,0.d0,lm))
      enddo

      wlep = wlep * 2.d0/3.d0/2.d0/s2

      do k = 1,3
         um = mquu(k)
         dm = mqud(k)
         wquark = wquark 
     >        -(s - (um+dm)/2.d0)*B0(s,um,dm) + s/3.d0 +
     >        um*B0(0d0,um,um) + 
     >        dm*B0(0d0,dm,dm) + 
     >        (um - dm)**2/2.d0/s* (B0(s,um,dm) - B0(0.d0,um,dm))
      enddo

      wquark = wquark * 2.d0/3.d0/2.d0/s2 * 3.d0

      wbos1 = (2.d0*mw2+5.d0*s)*B0(s,mw2,mg2)-2.d0*mw2*B0(0.d0,mw2,mw2)
     >        - mw2**2/s*(B0(s,mw2,mg2)-B0(0.d0,mw2,mg2)) + s/3.d0
      wbos1 = wbos1*2.d0/3.d0

      wbos2 = ((40.d0*c2-1.d0)*s+(16.d0*c2+54.d0-10.d0/c2)*mw2)
     >        *B0(s,mw2,mz2)
     >        - (16.d0*c2 + 2.d0)
     >          *(mw2*B0(0.d0,mw2,mw2)+mz2*B0(0.d0,mz2,mz2)) 
     >        + (4.d0*c2 - 1.d0) * 2.d0/3.d0*s
     >        - (8.d0*c2 + 1.d0) * (mw2-mz2)**2/s
     >          *(B0(s,mw2,mz2)-B0(0.d0,mw2,mz2))
      wbos2 = wbos2/12.d0/s2

      wbos3 = (2.d0*mh2 - 10.d0*mw2 - s)*B0(s,mw2,mh2) -
     >        2.d0*mw2 * B0(0.d0,mw2,mw2) - 2.d0*mh2*B0(0.d0,mh2,mh2) -
     >        (mw2-mh2)**2/s*(B0(s,mw2,mh2)-B0(0.d0,mw2,mh2)) 
     >         -2.d0/3.d0*s
      wbos3 = wbos3/12.d0/s2
      wbos = wbos1 + wbos2  + wbos3
      wlep   = -alpha/4.d0/pi * wlep
      wbos   = -alpha/4.d0/pi * wbos
      wquark = -alpha/4.d0/pi * wquark
      wselfnonresonant = wlep + wbos + wquark
      return
      end

c$$$      double complex function zzself(s)
c$$$      include 'shared.inc'
c$$$      double precision mlep(3)
c$$$      double precision lm,mw2,mh2,mz2,mg2
c$$$      double precision mquu(3),mqud(3)
c$$$      double complex B0,wlep,wbos,wbos1,wbos2,wbos3,wquark
c$$$      external B0
c$$$
c$$$      mg2 = getlambda()
c$$$
c$$$      s2 = s2tw
c$$$      c2 = 1.d0 - s2
c$$$
c$$$      mw2 = mw*mw
c$$$      mz2 = mz*mz
c$$$      mh2 = mh*mh
c$$$
c$$$      mlep(1) = me**2
c$$$      mlep(2) = mmu**2
c$$$      mlep(3) = mtau**2
c$$$
c$$$      mquu(1) = mu**2
c$$$      mquu(2) = mc**2
c$$$      mquu(3) = mt**2 
c$$$
c$$$      mqud(1) = md**2 
c$$$      mqud(2) = ms**2
c$$$      mqud(3) = mb**2 
c$$$
c$$$      wlep    = 0.d0
c$$$      wquark  = 0.d0
c$$$      wbos    = 0.d0
c$$$      wbos1   = 0.d0
c$$$      wbos2   = 0.d0
c$$$      wbos3   = 0.d0
c$$$
c$$$      do k = 1,3
c$$$         lm = mlep(k)
c$$$         wlep = wlep +
c$$$     >        -(s - lm/2.d0)*B0(s,0.d0,lm) + s/3.d0 +
c$$$     >        lm*B0(0d0,lm,lm) + 
c$$$     >        lm**2/2.d0/s* (B0(s,0.d0,lm) - B0(0.d0,0.d0,lm))
c$$$      enddo
c$$$
c$$$      wlep = wlep * 2.d0/3.d0/2.d0/s2
c$$$
c$$$      do k = 1,3
c$$$         um = mquu(k)
c$$$         dm = mqud(k)
c$$$         wquark = wquark +
c$$$     >        -(s - (um+dm)/2.d0)*B0(s,um,dm) + s/3.d0 +
c$$$     >        um*B0(0d0,um,um) + 
c$$$     >        dm*B0(0d0,dm,dm) + 
c$$$     >        (um - dm)**2/2.d0/s* (B0(s,um,dm) - B0(0.d0,um,dm))
c$$$      enddo
c$$$
c$$$      wquark = wquark * 2.d0/3.d0/2.d0/s2 * 3.d0
c$$$
c$$$      wbos1 = (2.d0*mw2+5.d0*s)*B0(s,mw2,mg2)-2.d0*mw2*B0(0.d0,mw2,mw2)
c$$$     >        - mw2**2/s*(B0(s,mw2,mg2)-B0(0.d0,mw2,mg2)) + s/3.d0
c$$$      wbos1 = wbos1*2.d0/3.d0
c$$$
c$$$      wbos2 = ((40.d0*c2-1.d0)*s+(16.d0*c2+54.d0-10.d0/c2)*mw2)
c$$$     >        *B0(s,mw2,mz2)
c$$$     >        - (16.d0*c2 + 2.d0)
c$$$     >          *(mw2*B0(0.d0,mw2,mw2)+mz2*B0(0.d0,mz2,mz2)) 
c$$$     >        + (4.d0*c2 - 1.d0) * 2.d0/3.d0*s
c$$$     >        - (8.d0*c2 + 1.d0) * (mw2-mz2)**2/s
c$$$     >          *(B0(s,mw2,mz2)-B0(0.d0,mw2,mz2))
c$$$      wbos2 = wbos2/12.d0/s2
c$$$
c$$$      wbos3 = (2.d0*mh2 - 10.d0*mw2 - s)*B0(s,mw2,mh2) -
c$$$     >        2.d0*mw2 * B0(0.d0,mw2,mw2) - 2.d0*mh2*B0(0.d0,mh2,mh2) -
c$$$     >        (mw2-mh2)**2/s*(B0(s,mw2,mh2)-B0(0.d0,mw2,mh2)) 
c$$$     >         -2.d0/3.d0*s
c$$$      wbos3 = wbos3/12.d0/s2
c$$$      wbos = wbos1 + wbos2  + wbos3
c$$$      wlep   = -alpha/4.d0/pi * wlep
c$$$      wbos   = -alpha/4.d0/pi * wbos
c$$$      wquark = -alpha/4.d0/pi * wquark
c$$$      wselfnonresonant = wlep + wbos + wquark
c$$$      return
c$$$      end
