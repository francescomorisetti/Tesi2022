      function elmat2form(p1,p2,p3,p4,k,Qu,Qd,Qw,Qe)
      implicit none
      common/fordebugging/idebugging
!! u dbar --> W+ --> e+ nue gamma
!! p1 p2             p4 p3  k
!! attention: it is different from my usual convention!!
      integer i,idebugging

      double precision me,m1,m2,ch1,ch2,chfs,mfs
      integer lepton,imirror

      double precision dot,elmat2form,q2,mup,tridot
      complex*16 elmat2,ecomplex_,e_
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),k(0:3),
     > ptmp(0:3), q(0:3)
      
      double precision kp1,kp2,kp3,kp4,kp12,kp22,kp32,kp42
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      double precision p1p22,p1p32,p1p42,p2p32,p2p42,p3p42
      double precision p1p1,p2p2,p3p3,p4p4
      double precision me2,mup2,md2
      double precision me4,mup4,md4
*      double precision mwm2,mwm4,mwm6,mwm8
      complex*16 mw2m1,mw2m2,mw2m3,mw2m4,mwc2,mwc2d
      complex*16 mwc2m1,mwc2m2,mwc2m3,mwc2m4
      complex*16 mwc2dm1,mwc2dm2,mwc2dm3,mwc2dm4
      double precision kp1m1,kp2m1,kp3m1,kp4m1,kp1m2,kp2m2,kp3m2,kp4m2

      complex*16 im,s,sf,prop12,prop12d,prop34,prop34d,i_
      double precision prop12m2,prop34m2
      double precision absSf2,propdo,absS2,propel
      complex*16 propSc,propSf,propSfc,propSSfc
      double precision absSf
      double precision propup

      double precision meorig,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt
      double precision mw2,charges(4),amasses(4)
!      double precision eps1234
      complex*16 eps1234
      double precision gw,gz,el2,Qu,Qd,Qe,Qu2,Qd2,Qe2,gwrun,gamw
!      complex*16 Qw,Qw2       ! due to Baur-Zeppenfeld recipe
      double precision Qw,Qw2,bet,xx,ombet
      double precision alpha,alpha_gf,s2tw,gf,pi,convfac
      common/t_widths/gw,gz
      double precision dum1(-5:5),dum2(-5:5),dum3(-5:5)
      common/masseshor/meorig,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     dum1,dum2,dum3
      common/processhor/m1,m2,ch1,ch2,chfs,mfs,lepton,imirror
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/chargesmasses_cm/charges,amasses
      data im /(0.d0,1.d0)/

      i_ = im

      gamw = gw/mw

!      Qu =  2.d0/3.d0
!      Qd = -1.d0/3.d0
!      Qw =  1.d0                    ! Lopez-Castro recipe
!      Qw =  1.d0 * (1.d0 + im*gamw) ! Baur - Zeppenfeld recipe
!      Qe = -1.d0

      Qe2 = Qe*Qe
      Qd2 = Qd*Qd
      Qw2 = Qw*Qw
      Qu2 = Qu*Qu

      el2 = alpha*4.d0*pi
      
      kp1 = dot(k,p1)
      kp2 = dot(k,p2)
      kp3 = dot(k,p3)
      kp4 = dot(k,p4)

      mup  = amasses(1)
      me   = amasses(3)
      me2  = me*me
      mup2 = mup*mup
      md2  = amasses(2)*amasses(2)

      mup4 = mup2*mup2
      md4  = md2*md2
      me4  = me2*me2

cc ! trying to avoid NaNs at extremely high energies      
      if (kp4.eq.0.d0) then
c         bet   = sqrt(1.d0 - me2/p4(0)/p4(0))
c         ombet = 1.d0 - bet
         xx    = -me2/p4(0)/p4(0)
         ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .        - 0.078125d0 * xx*xx*xx)
         kp4 = p4(0)*k(0)*ombet
      endif
cc
      
      p1p1 = dot(p1,p1)
      p2p2 = dot(p2,p2)
      p3p3 = dot(p3,p3)
      p4p4 = dot(p4,p4)

      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4)
      p3p4 = dot(p3,p4)

      kp12 = kp1*kp1
      kp22 = kp2*kp2
      kp32 = kp3*kp3
      kp42 = kp4*kp4
      
      p1p22 = p1p2*p1p2
      p1p32 = p1p3*p1p3
      p1p42 = p1p4*p1p4
      p2p32 = p2p3*p2p3
      p2p42 = p2p4*p2p4
      p3p42 = p3p4*p3p4


*      mwm2 = mw**(-2)
*      mwm4 = mw**(-4)
*      mwm6 = mw**(-6)
*      mwm8 = mw**(-8)

      mw2  = mw*mw
      mwc2  = mw2 - im*gw*mw
      mwc2d = mw2 + im*gw*mw

      do i = 0,3
         ptmp(i) = p1(i) + p2(i)
      enddo      
      q2 = dot(ptmp,ptmp)
      gwrun = mw2*gamw           ! Lopez-Castro recipe
!      gwrun = q2*gamw            ! Baur - Zeppenfeld recipe
      s  = 1.d0/(q2 - mwc2)

cc      if (idebugging.eq.1) print*,'sq.f 1  q2 =',q2

!      s  = 1.d0/(q2 - mw2) ! null width!!

      do i = 0,3
         ptmp(i) = p4(i) + p3(i)
      enddo
      q2 = dot(ptmp,ptmp)
      gwrun = mw2*gamw           ! Lopez-Castro recipe
!      gwrun = q2*gamw            ! Baur - Zeppenfeld recipe
      sf  = 1.d0/(q2 - mwc2)
!      sf  = 1.d0/(q2 - mw2) ! null width!!

cc      if (idebugging.eq.1) print*,'sq.f 2  q2 =',q2

      prop12  = s
      prop12d = conjg(s)
      prop34  = sf
      prop34d = conjg(sf)

cc      if (idebugging.eq.1) then
cc         print*,'sq.f 3  prop1 =',prop12*prop12d
cc         print*,'sq.f 4  prop2 =',prop34*prop34d
cc      endif

      mw2m1 = 1.d0/mw2 ! Ward identities-violating or null width
      mw2m1 = 1.d0/(mw2 - im*gw*mw)     ! Lopez-Castro recipe

!      mw2m1 = 1.d0/mw2 * (1.d0+im*gamw) ! Baur - Zeppenfeld recipe
      mw2m2 = mw2m1*mw2m1
      mw2m3 = mw2m2*mw2m1
      mw2m4 = mw2m3*mw2m1

      mwc2m1  = 1.d0/mwc2
      mwc2dm1 = 1.d0/mwc2d

!!!!!!!!!!!!!!!!!!!!!!!!!!
!      mwc2m1  = 1.d0/mw2
!      mwc2dm1 = 1.d0/mw2
!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mwc2m2  = mwc2m1*mwc2m1
      mwc2m3  = mwc2m2*mwc2m1
      mwc2m4  = mwc2m3*mwc2m1
      mwc2dm2 = mwc2dm1*mwc2dm1
      mwc2dm3 = mwc2dm2*mwc2dm1
      mwc2dm4 = mwc2dm3*mwc2dm1

      kp1m1 = 1.d0/kp1
      kp2m1 = 1.d0/kp2
      kp3m1 = 1.d0/kp3
      kp4m1 = 1.d0/kp4

      kp1m2 = kp1m1*kp1m1
      kp2m2 = kp2m1*kp2m1
      kp3m2 = kp3m1*kp3m1
      kp4m2 = kp4m1*kp4m1

      propel =  0.5d0 * kp4m1
      propup = -0.5d0 * kp1m1
      propdo = -0.5d0 * kp2m1

      prop12m2 = prop12*prop12d
      prop34m2 = prop34*prop34d
      do i = 0,3
         q(i) = k(i)
      enddo
      eps1234 = e_(p1,p2,p3,p4)
      include 'form/formme.f'
!      elmat2form = RealPart(elmat2) * el2**3/s2tw**2 / 1024d0
      elmat2form = elmat2 * el2**3/s2tw**2 / 1024d0

      return
      end

***********************************
      function elmat2born(p1,p2,p3,p4)
      implicit none
!! u dbar --> W+ --> e+ nue
!! p1 p2             p4 p3
!! attention: it is different from my usual convention!!
      integer i

      double precision dot,elmat2form,M2,q2,mup,elmat2born
      complex*16 elmat2,ecomplex_,e_
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),k(0:3),
     > ptmp(0:3), q(0:3)
      
      double precision kp1,kp2,kp3,kp4,kp12,kp22,kp32,kp42
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      double precision p1p22,p1p32,p1p42,p2p32,p2p42,p3p42
      double precision p1p1,p2p2,p3p3,p4p4
      double precision me2,mup2,md2,meorig
      double precision me4,mup4,md4
*      double precision mwm2,mwm4,mwm6,mwm8
      complex*16 mw2m1,mw2m2,mw2m3,mw2m4,mwc2,mwc2d
      complex*16 mwc2m1,mwc2m2,mwc2m3,mwc2m4
      complex*16 mwc2dm1,mwc2dm2,mwc2dm3,mwc2dm4
      double precision kp1m1,kp2m1,kp3m1,kp4m1,kp1m2,kp2m2,kp3m2,kp4m2

      complex*16 im,s,sf,prop12,prop12d,prop34,prop34d,i_
      double precision prop12m2,prop34m2
      double precision absSf2,propdo,absS2,propel
      complex*16 propSc,propSf,propSfc,propSSfc
      double precision absSf
      double precision propup

      double precision me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt
      double precision mw2,propmodulo2,charges(4),amasses(4)
!      double precision eps1234
      complex*16 eps1234
      double precision gw,gz,el2,Qu,Qd,Qe,Qu2,Qd2,Qe2,gwrun,gamw
!      complex*16 Qw,Qw2       ! due to Baur-Zeppenfeld recipe
      double precision Qw,Qw2
      double precision alpha,alpha_gf,s2tw,gf,pi,convfac
      common/t_widths/gw,gz
      double precision dum1(-5:5),dum2(-5:5),dum3(-5:5)
      common/masseshor/meorig,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     dum1,dum2,dum3
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/chargesmasses_cm/charges,amasses
      data im /(0.d0,1.d0)/

      i_ = im

      gamw = gw/mw

      el2 = alpha*4.d0*pi
            
      p1p1 = dot(p1,p1)
      p2p2 = dot(p2,p2)
      p3p3 = dot(p3,p3)
      p4p4 = dot(p4,p4)

      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4)
      p3p4 = dot(p3,p4)
      
      p1p22 = p1p2*p1p2
      p1p32 = p1p3*p1p3
      p1p42 = p1p4*p1p4
      p2p32 = p2p3*p2p3
      p2p42 = p2p4*p2p4
      p3p42 = p3p4*p3p4

      me   = amasses(3)
      mup2 = amasses(1)**2
      md2  = amasses(2)**2
      me2  = me*me

      mup4 = mup2*mup2
      md4  = md2*md2
      me4  = me2*me2

      mw2   = mw*mw
      mwc2  = mw2 - im*gw*mw
      mwc2d = mw2 + im*gw*mw

      do i = 0,3
         ptmp(i) = p1(i) + p2(i)
      enddo      
      q2 = dot(ptmp,ptmp)
      gwrun = mw2*gamw           ! Lopez-Castro recipe
!      gwrun = q2*gamw            ! Baur - Zeppenfeld recipe
      s  = 1.d0/(q2 - mwc2)
!      s  = 1.d0/(q2 - mw2) ! null width!!

      prop12  = s
      prop12d = conjg(s)

      propmodulo2 = prop12*prop12d


      mw2m1 = 1.d0/mw2            ! Ward identities-violating or null width
cccc      mw2m1 = 1.d0/(mw2 - im*gw*mw)     ! Lopez-Castro recipe

!      mw2m1 = 1.d0/mw2 * (1.d0+im*gamw) ! Baur - Zeppenfeld recipe
      mw2m2 = mw2m1*mw2m1
      mw2m3 = mw2m2*mw2m1
      mw2m4 = mw2m3*mw2m1

      mwc2m1  = 1.d0/mwc2
      mwc2dm1 = 1.d0/mwc2d

CCC no width in kmuknu
      mwc2m1  = 1.d0/mw2
      mwc2dm1 = 1.d0/mw2
CCC
      mwc2m2  = mwc2m1*mwc2m1
      mwc2m3  = mwc2m2*mwc2m1
      mwc2m4  = mwc2m3*mwc2m1
      mwc2dm2 = mwc2dm1*mwc2dm1
      mwc2dm3 = mwc2dm2*mwc2dm1
      mwc2dm4 = mwc2dm3*mwc2dm1
      include 'form/bornme.f'
!      elmat2born = RealPart(elmat2)*el2**2/s2tw**2/64d0 * propmodulo2
      elmat2born = elmat2*el2**2/s2tw**2/64d0 * propmodulo2
      return
      end
*
      function bornm2ale(p1,p2,p3,p4)
! 1--> u, 2--> dbar, 3--> neutrino, 4--> positron
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),ptmp(0:3)
      double precision me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,mw2
      double precision gw,gz
      double precision alpha,alpha_gf,s2tw,gf,pi,convfac


      complex*16 M2,eps1234,mwc2m1,mwc2dm1,im,e_


      common/t_widths/gw,gz
      double precision dum1(-5:5),dum2(-5:5),dum3(-5:5)
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     dum1,dum2,dum3
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac

      data im /(0.d0,1.d0)/


      mw2 = mw*mw
      mup2 = mu*mu
      me2 = me*me
      md2 = md*md

      mup4 = mup2*mup2
      me4 = me2*me2
      md4 = md2*md2

      mwc2m1  =  1.d0/(mw2 - im*gw*mw)
      mwc2dm1 =  1.d0/(mw2 + im*gw*mw)


      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4)
      p3p4 = dot(p3,p4)

      do k = 0,3
         ptmp(k) = p1(k) + p2(k)
      enddo
      s = dot(ptmp,ptmp)

      propag2 = 1.d0/( (s-mw2)**2 + gw**2*mw2  )

!      eps1234 = e_(p1,p2,p3,p4)
      eps1234 = 0.d0 ! per conservazione!!!!

      M2=        (4.*(-2.*md2*me2*mw2*p1p3 -
     -     1.*(2.*me2*mup2*mw2 - 4.*(Gw**2*mw2 + mw2**2)*p1p4)*
     -     p2p3 + me2*(2.*md2*mup2 + (md2 + mup2)*p1p2)*p3p4))/
     -     ((Gw**2*mw2 + mw2**2)*
     -     (Gw**2*mw2 + (md2 + mup2 - 1.*mw2 + 2.*p1p2)**2))
      bornm2ale = M2
      return
      end
