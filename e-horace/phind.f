      subroutine phasespacephindz_OLD(p1,p2,p3,p4,p5,qph,m1,m2,m3,
     .     w,phsp,ie)
!for Z
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
      dimension dir3(0:3),dir4(0:3),ptmp(0:3),vers(0:3),pcoll(0:3)
      double precision kone(0:3),k1(0:3),k2(0:3),pp(0:3),pw(0:3)
      double precision lambda
      double precision mq(-5:5),chq(-5:5),itre(-5:5)
      character*1 boson
      integer emissionpattern(1)
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (imaxph = 40)
      double precision qph(imaxph,0:3),oph(imaxph)
      common/nrandom/ir
      real*4 rnd(2),csi(1)
      common/radpattern/nph(4)      
      common/phindchannel/bbc,pw0,qmodc,c1qc,ich
! from shared.inc      
      common/exp_cuts/etamax,ptmin,cmin,cmax,ptminmiss,
     .     tmmin,mass2,ismear,irec,icutmt,ischemealpha
      common/vectorboson/boson
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     mq,chq,itre
      common/t_widths/gw,gz
      common/propstest/testprop1,testprop2
      common/testphspz/phspmin,phspmax
      data phspmin,phspmax /1d12,-1d12/
      common/icountphspphindz/icount
      data icount /0/
      icount = icount + 1

! using the same routine as for W
c      call phasespacephind(p1,p2,p3,p4,p5,qph,m1,m2,m3,
c     .     w,phsp,ie)
c      return
      
      
      nphot = 1

      ie       = 0
      w        = 1.d0
      phsp     = 1.d0

      ei = p1(0)+p2(0)
      if (ei.lt.(m1+m2+m3)) then 
         phsp = 0.d0
         w    = 0.d0
         ie   = 1
         return
      endif

*****************
      beta    = (1.d0+m3/p1(0))*(1.d0-m3/p1(0))
      beta    = sqrt(beta)
      beta5   = beta
      beta3   = sqrt(1.d0-m1**2/max(p1(0),2.d0*m1)**2)

      dir3(0) =  p1(0)
      dir3(1) =  0.d0
      dir3(2) =  0.d0
      dir3(3) =  p1(0)*beta
      dir4(0) =  dir3(0)
      ptmp(0) =  p1(0)
      ptmp(1) =  0.d0
      ptmp(2) =  0.d0
      ptmp(3) = -p1(0)*beta3
      do k = 1,3
         dir4(k) = -dir3(k)
         dir3(k) = dir3(k)/beta*beta3
      enddo               
      beta1 = sqrt(tridot(p1,p1))/p1(0)
************************
      vm = mw
      vg = gw
      if (boson.eq.'Z') then
         vm = mz
         vg = gz
      endif

      q2max = (ei-m3)**2
*************
!modified 5 12 2007, CARLO
      q2minmass = (m1+m2)**2
      q2min = (max(tmmin,ptmin+ptminmiss))**2
      q2min = max(q2min,q2minmass)
      if (q2max.lt.q2min) then 
         phsp = 0.d0
         w    = 0.d0
         ie   = 1
         return
      endif
************

      gm = vm*vg
      A     =  atan(q2min/gm - vm/vg)
      anbw  = (atan(q2max/gm - vm/vg) - A)/gm

      pflatc = 0.5d0

      pich1 = 0.5d0
      pich2 = 1.d0 - pich1

      call wraprng(csi,1)
      ich = 2
      if (csi(1).le.pich1) ich = 1

      if (ich.eq.1) then

         pbw = 0.5d0
         pphot = 1.d0 -pbw
         call wraprng(rnd,2)

         if (rnd(1).lt.pbw) then
            q2 = gm*tan(gm*anbw*rnd(2)+A) + vm*vm
            w = w * anbw * ((q2 - vm*vm)**2 + gm*gm)
         else
            an = log(q2max/q2min)
            q2 = q2min * exp(an*rnd(2))
            w  = w * an * q2
         endif

         sqla = sqrt(lambda(ei**2,q2,m3**2))
         pmod = sqla/2.d0/ei
         p5(0) = sqrt(pmod**2 + m3**2)

         call wraprng(csi,1)
!! NEVER FLAT!!
         if (csi(1).lt.1.d0) then
            dir4(0) = p5(0)
            dir4(1) = 0.d0
            dir4(2) = 0.d0
            dir4(3) = -pmod
            call wraprng(rnd,2)
            call collinear(0.d0,rnd(1),rnd(2),dir4,vers,wverg)
            w = w * wverg
            c = vers(3)
            s = sqrt(abs(1.d0 - c*c))
            if (s.gt.0.d0) then
               sphi = vers(1)/s
               cphi = vers(2)/s
            else
               sphi = 0.d0
               cphi = 1.d0
            endif
         else
            call wraprng(rnd,2)
            c = 2.d0*rnd(1) - 1.d0
            phi = 2.d0*pi*rnd(2)
            w   = w * 4.d0 * pi
            sphi = sin(phi)
            cphi = cos(phi)
            s = sqrt(abs(1.d0 - c*c))
         endif

         p5(1) = pmod * sphi * s
         p5(2) = pmod * cphi * s
         p5(3) = pmod * c
         pw(0) = ei - p5(0)
         pw(1) = - p5(1)
         pw(2) = - p5(2)
         pw(3) = - p5(3)
         phsp1 = sqla/8.d0/ei**2
         
         call wraprng(rnd,2)
         c = 2.d0*rnd(1) - 1.d0
         s = sqrt(abs(1.d0-c*c))
         phi = 2.d0*pi*rnd(2)
         w   = w *4.d0*pi     
         sphi = sin(phi)
         cphi = cos(phi)
         
         sqla = sqrt(lambda(q2,m1**2,m2**2))
         p3mod = sqla/2.d0/sqrt(q2)
         p3(0) = sqrt(p3mod**2 + m1**2)
         p3(1) = p3mod * sphi * s
         p3(2) = p3mod * cphi * s
         p3(3) = p3mod * c
         p4(0) = sqrt(q2) - p3(0)
         p4(1) = - p3(1)
         p4(2) = - p3(2)
         p4(3) = - p3(3)
         
         phsp2 = sqla/8.d0/q2
         
         call new_boost(pw,p3,p3,-1)
         call new_boost(pw,p4,p4,-1)

         phsp = phsp1*phsp2/(2.d0*pi)**5
ch1 OK
      else
*-----------
         call wraprng(csi,1)
         q2 = (q2max-q2min)*csi(1) + q2min
         w = w * (q2max-q2min)

         sqla = sqrt(lambda(ei**2,q2,m3**2))
         pmod = sqla/2.d0/ei
         p5(0) = sqrt(pmod**2 + m3**2)

         call wraprng(rnd,2)
!         if (rnd(1).lt.0.5d0) then
         if (rnd(1).lt.0.d0) then
            ab = 2.d0*p1(0)*p5(0)*beta1*pmod/p5(0)
            bb = -2.d0*p1(0)*p5(0) - vm*vm + dot(p1,p1) + m3*m3
            cb = gm
            ak = atan((-ab+bb)/cb)
            ancbwt =  (atan((ab+bb)/cb) - ak)/ab/cb
            c = (cb*tan(ab*cb*ancbwt*rnd(2)+ak) - bb)/ab
            w = w * ancbwt*((ab*c+bb)**2 + cb*cb)
         else
            ag = 2.d0*p1(0)*p5(0) - dot(p1,p1) - m3*m3
            bg = 2.d0*p1(0)*p5(0)*beta1*pmod/p5(0)
            agmbg = ag-bg

!      trucco a-b = (a^2-b^2)/(a+b)
            am12 = dot(p1,p1)
            am52 = m3**2
            sam2 = am12+am52
            e1   = p1(0)
            e5   = p5(0)
            agmbgv2 = 4.d0*e1**2*am52+4.d0*e5**2*am12-4.d0*e1*e5*sam2
     .           +sam2**2-4.d0*am12*am52
            agmbgv2 = agmbgv2/(ag+bg)
*****
            agmbg = agmbgv2

            if (agmbg.lt.0.d0) then
               print*,'phind.f line 225 ---!!',agmbg
               bg = ag - 1.d-13
               agmbg = 1.d-13
            endif

            an = -1.d0/bg*log((agmbg)/(ag+bg))
            c  = ag - (ag+bg)*exp(-an*bg*rnd(2))
            c  = c/bg
            w = w * an * (ag - bg * c)
         endif
         call wraprng(csi,1)
         phi = 2.d0*pi*csi(1)
         w   = w * 2.d0*pi

         s = sqrt(abs(1.d0 - c*c))
         sphi = sin(phi)
         cphi = cos(phi)
         
         p5(1) = pmod * sphi * s
         p5(2) = pmod * cphi * s
         p5(3) = pmod * c
         pw(0) = ei - p5(0)
         pw(1) = - p5(1)
         pw(2) = - p5(2)
         pw(3) = - p5(3)
         qmod = sqrt(tridot(pw,pw))

         phsp1 = sqla/8.d0/ei/ei

         call wraprng(rnd,2)
         call collinear(0.d0,rnd(1),rnd(2),ptmp,vers,wverg)
         w = w * wverg

         c1q  = tridot(vers,pw)/qmod
         
         b2 = qmod*qmod
         ak = 2.d0*qmod*c1q

         U  = pw(0)**2 + m1**2 - m2**2 - b2
         aa = 4.d0*pw(0)**2 - ak*ak
         bb = -2.d0*U*ak
         cc = 4.d0*pw(0)**2*m1*m1 - U*U

         arg = bb*bb - 4.d0*aa*cc
         if (arg.lt.0.d0) then
            w = 0.d0
            phsp = 0.d0
            ie = 1
            return
         endif

         p3mod1 = (-bb + sqrt(arg))/2.d0/aa
         p3mod2 = (-bb - sqrt(arg))/2.d0/aa

         p3jac = 1.d0
         if (p3mod1.ge.0.d0.and.p3mod2.lt.0.d0) then
            p3mod = p3mod1
         elseif(p3mod1.lt.0.d0.and.p3mod2.ge.0.d0) then
            p3mod = p3mod2
         elseif(p3mod1.ge.0.d0.and.p3mod2.ge.0.d0) then
            p3jac = 2.d0
            p3mod = p3mod1
            call wraprng(csi,1)
            if (csi(1).le.0.5d0) then
               p3mod = p3mod2
            endif
         else
            w = 0.d0
            phsp = 0.d0
            ie = 1
            return
         endif
         w = w * p3jac

         p3(0) = sqrt(p3mod**2+m1*m1)
         p3(1) = p3mod*vers(1)
         p3(2) = p3mod*vers(2)
         p3(3) = p3mod*vers(3)
         p4(0) = pw(0) - p3(0)
         p4(1) = pw(1) - p3(1)
         p4(2) = pw(2) - p3(2)
         p4(3) = pw(3) - p3(3)

         bbb = p3mod/p3(0)
         phsp2 = bbb*p3mod/abs(bbb*pw(0)-qmod*c1q)/4.d0

         bbc = bbb
         pw0 = pw(0)
         qmodc = qmod
         c1qc = c1q

         phsp = phsp1*phsp2/(2.d0*pi)**5
!! mu+ mu- symmetrization
         call wraprng(csi,1)
         if (csi(1).lt.0.5d0) call exchange_mom(p3,p4)

cc         w = w *2.d0 !! THIS WAS A BUG !!!
      endif

      c5 = p5(3)/sqrt(tridot(p5,p5))

      ab = 2.d0*p1(0)*p5(0)*beta1*pmod/p5(0)
      bb = -2.d0*p1(0)*p5(0) - vm*vm + dot(p1,p1) + m3*m3
      cb = gm
      ak = atan((-ab+bb)/cb)
      anbwt =  (atan((ab+bb)/cb) - ak)/ab/cb

      anq2g = log(q2max/q2min)


      ag  = 2.d0*p1(0)*p5(0) - dot(p1,p1) - m3*m3
      bg  = 2.d0*p1(0)*p5(0)*beta1*pmod/p5(0)

      u1 = 2.d0*p1(0)*p5(0)
      u2 = dot(p1,p1) + m3*m3
      u3 = beta1*pmod/p5(0)
      agmbg = u1*(1.d0 - u2/u1-u3)
c      b1t = 1.d0 - beta1
c      b5t = 1.d0 - pmod/p5(0)
c      agmbgv2 = b1t+b5t - b1t*b5t - u2/u1
c      agmbgv2 = u1 * agmbgv2

***
!      trucco a-b = (a^2-b^2)/(a+b)
      am12 = dot(p1,p1)
      am52 = m3**2
      sam2 = am12+am52
      e1   = p1(0)
      e5   = p5(0)
      agmbgv2 = 4.d0*e1**2*am52+4.d0*e5**2*am12-4.d0*e1*e5*sam2
     .     +sam2**2-4.d0*am12*am52
      agmbgv2 = agmbgv2/(ag+bg)
***
      agmbg = agmbgv2

      if (agmbg.lt.0.d0) then
         print*,'phind.f line 361 ---!!',agmbg
         bg = ag - 1.d-13
         agmbg = 1.d-13
      endif

      angt = -1.d0/bg*log((agmbg)/(ag+bg))

      do k = 0,3
         ptmp(k) = p3(k)+p4(k)
      enddo
      q2 = dot(ptmp,ptmp)

      b5    = sqrt(tridot(p5,p5))/p5(0)
      b3    = sqrt(tridot(p3,p3))/p3(0)
      b4    = sqrt(tridot(p4,p4))/p4(0)
      anc5  = 1.d0/b5*log((1.d0+b5)/(1.d0-b5))
      anc3  = 1.d0/b3*log((1.d0+b3)/(1.d0-b3))
      anc4  = 1.d0/b4*log((1.d0+b4)/(1.d0-b4))

      rbw   = 1.d0/anbw/((q2-vm*vm)**2+gm*gm)    !  * anbw 
      rq2f  = 1.d0/(q2max - q2min)               !  * (q2max - q2min)
      rq2g  = 1.d0/q2/anq2g                      !  * anq2g
      rcmuf = 1.d0/2.d0                          !  * 2.d0
      rcuf  = 1.d0/2.d0                          !  * 2.d0

! never flat!
      rcuf  = 0.d0

      rcug  = 1.d0/dot(p2,p5)*p5(0)*p2(0)/anc5   !  * anc5
      rcmumg = 1.d0/dot(p2,p3)*p3(0)*p2(0)/anc3  !  * anc3
      rcmupg = 1.d0/dot(p2,p4)*p4(0)*p2(0)/anc4  !  * anc4
      rbwt  = 1.d0/((ab*c5+bb)**2 + gm*gm)/anbwt !  * anbwt 

      rbwt  = 0.d0

      ddd = (ag - bg * c5)
      if (ddd.le.0.d0) then
         print*,'phind.f line ~394!!!'
         print*,'EVENT REJECTED,printing p1,p2,p3,p4,p5'
         print*,p1
         print*,p2
         print*,p3
         print*,p4
         print*,p5
         ie = 1
         phsp = 0.d0
         w    = 0.d0
         return
      endif

      rgt   = 1.d0/angt/ddd           !  * angt

      reg1 = (rbw+rq2g) * rcmuf * (rcuf+rcug)
      reg2 = rq2f * (rbwt+rgt) * (rcmumg+rcmupg)

      if (ich.eq.1) w = w * reg1/(pich1*reg1 + pich2*reg2)
      if (ich.eq.2) w = w * reg2/(pich1*reg1 + pich2*reg2)

      qph(1,0) = p5(0)
      qph(1,1) = p5(1)
      qph(1,2) = p5(2)
      qph(1,3) = p5(3)

      return
      end
**********************************************************************
      subroutine phasespacephind(p1,p2,p3,p4,p5,qph,m1,m2,m3,
     .     w,phsp,ie)
! for W
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3),p34(0:3)
      dimension p1boo(0:3),p2boo(0:3)
      dimension dir3(0:3),dir4(0:3),ptmp(0:3),vers(0:3),pcoll(0:3)
      double precision kone(0:3),k1(0:3),k2(0:3),pp(0:3),pw(0:3)
      double precision lambda
      double precision mq(-5:5),chq(-5:5),itre(-5:5)
      integer emissionpattern(1)
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (imaxph = 40)
      double precision qph(imaxph,0:3),oph(imaxph)
      common/nrandom/ir
      real*4 rnd(2),csi(1)
      common/radpattern/nph(4)      
      common/phindchannel/bbc,pw0,qmodc,c1qc,ich
! from shared.inc      
      character*1 boson
      common/vectorboson/boson
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     mq,chq,itre
      common/t_widths/gw,gz

      common/phspindcmn/ifirst
      data ifirst /0/
      
      common/exp_cuts/etamax,ptmin,cmin,cmax,ptminmiss,
     .     tmmin,mass2,ismear,irec,icutmt,ischemealpha

      external lambda
      double precision ml1,ml2,ml3,ml4,ml5,ml34,ml342
!! trying again

      common/mvgvbwaMCphind/ambwxs,agbwxs
      
      
      if (ifirst.eq.0) then
         print*,'Print in phind phasespace:'
         do j = 1,10
            print*,'tunare pch1 e pch2'
c            print*,'pt lepton cut hard wired'
         enddo
         ifirst = 1
      endif
      
      nphot = 1
      npart = nphot + 2
      nvar  = 3*nphot + 2

      ncharged = 4
      if (boson.eq.'W') ncharged = 3
      ie       = 0
      w        = 1.d0
      duepigrechi = (2.d0*pi)**(-nvar)

      phsp     = duepigrechi
      
      ei = p1(0)+p2(0)
      if (ei.lt.(m1+m2+m3)) then 
         phsp = 0.d0
         w    = 0.d0
         ie   = 1
         return
      endif
      
      ml1 = sqrt(dot(p1,p1))
      ml2 = 0.d0
      ml3 = m1
      ml4 = m2
      ml5 = m3

c      ml1 = 10.d0
c      ml2 = 10.d0
c      ml3 = 10.d0
c      ml4 = 20.d0
c      ml5 = 30.d0
      
*     * following the notes
*     * first sampling p3 + p4 as bw
      sl    = dot(p1+p2,p1+p2)
      ssl   = sqrt(sl)
      q2min = max((ml3+ml4)**2,tmmin*tmmin)
      q2max = (ssl - ml5)**2

      if (q2max.lt.q2min) then
         phsp = 0.d0
         w    = 0.d0
         ie   = 1
         return
      endif
      ml12 = ml1*ml1
      ml22 = ml2*ml2
      ml32 = ml3*ml3
      ml42 = ml4*ml4
      ml52 = ml5*ml5

      if (boson.eq.'W') then
         mv = mw
         gv = gw 
      else
         mv = mz
         gv = gz 
      endif
      ambwxs = mv
      agbwxs = gv

      mv2    = mv*mv
      gv2    = gv*gv
      mv2gv2 = mv2*gv2
      
** **
*     * **
      if (boson.eq.'W') then
c         call get_invm2(q2max,q2min,ml342,wm342)
         call bwxssample(q2max,q2min,ml342,wm342)

c     call bwovermcomplex(q2max,q2min,ml342,wm342)  ! not needed the phase space Q-> p3 p4 doesn't go
                                                       ! as 1/M!!! it was a mistake... aargh...
         w = w * wm342
      else                      ! if boson eq Z
         print*,'should not enter here!'
      endif
      
      ml34 = sqrt(ml342)

      qmod  = sqrt(lambda(sl,ml52,ml342)/sl)*0.5d0
      eq    = 0.5d0*(sl - ml52 + ml342)/ssl
      e5    = ssl - eq
      
      p1mod = sqrt(p1(0)*p1(0)-ml12) ! sqrt(tridot(p1,p1))
      p2mod = p2(0)                  !sqrt(tridot(p2,p2))      

** generating t1Q or t2Q to sample the propagator in diagram 3 or 2
c          get_tlimits(M2,ma2, mb2, m12,  m22, tp,tm)
*     I generate |t1Q| according to 1/(|t1Q| + ml52) in ich = 1
*     I generate  t2Q  according to spacelike bw in ich = 2
      tp1 = (p1(0)-eq)**2 - (p1mod - qmod)**2
      tm1 = (p1(0)-eq)**2 - (p1mod + qmod)**2
      tp2 = (p2(0)-eq)**2 - (p2mod - qmod)**2
      tm2 = (p2(0)-eq)**2 - (p2mod + qmod)**2         
            
      tmax = abs(tm1)
      tmin = abs(tp1)
      an1  = log(tmax+ml52) - log(tmin+ml52)
      tmax2 = abs(tm2)
      tmin2 = abs(tp2)
c      an2  = bwnorm(tp2,tm2)
      an2  = bwxsnorm(tmax2,tmin2)

c      pbwt = 0.5d0
      pbwt = anbwt/(anbwt+anoot)
      pgt  = 1.d0 - pbwt

c      pch1 = 0.25d0 ! seems to be the optimum
      pch1 = an1/(an1+an2)
      pch2 = 1.d0 - pch1
      ich  = 2
      call wraprng(csi,1)
      if (csi(1).le.pch1) ich = 1
      
      if (ich.eq.1) then
c         call get_tlimits(sl,ml12,ml22,ml342,ml52,tp,tm)
         tp = tp1
         tm = tm1
      else
c         call get_tlimits(sl,ml22,ml12,ml342,ml52,tp,tm)
         tp = tp2
         tm = tm2
      endif      

      
      if (ich.eq.1) then
         an   = an1
         call wraprng(csi,1)
         abst = exp(an*csi(1))*(tmin+ml52) - ml52
         w    = w * an * (abst + ml52)
         t    = -abst
      else
         if (boson.eq.'W') then
c            call get_invm2(tp,tm,t,wt)
            call bwxssample(tmax2,tmin2,t,wt)
            t = -t
            w = w * wt
         else ! if boson eq Z
            
         endif
      endif
            
      if (ich.eq.1) then
c         cq = (t - ml12 - ml342 + 2.d0*p1(0)*eq)/p1mod/qmod*0.5d0
         cq = 0.5d0*(t-tm)/p1mod/qmod - 1.d0
      else
c         cq = (t - ml22 - ml342 + 2.d0*p2(0)*eq)/p2mod/qmod*0.5d0
         cq = 0.5d0*(t-tm)/p2mod/qmod - 1.d0
         cq = -cq
      endif
      if (( abs(cq)-1.d0 ).ge.0.d0) then
         cq = cq/abs(cq)
      endif
      sq = sqrt(1.d0 - cq*cq)

      call wraprng(csi,1)
      phi  = 2.d0*pi*csi(1)
      sphi = sin(phi)
      cphi = cos(phi)
      w    = w * 2.d0*pi
      
      pw(0) = eq
      pw(1) = qmod*sq*sphi
      pw(2) = qmod*sq*cphi
      pw(3) = qmod*cq

cc      if (ich.eq.2) call rot(-1,p2,pw,pw)
      
      p5(0)   =  e5
      p5(1:3) = -pw(1:3)
      
      phsp = phsp * 0.25d0 / sqrt(lambda(sl,ml12,ml22)) ! for both channels

***   regulator for the two channels
      if (boson.eq.'W') then
         t1q = dot(p1-pw,p1-pw)
         t2q = dot(p2-pw,p2-pw)
         fch1 = 1.d0/abs(t1q - ml52)/an1
c         fch2 = 1.d0/( (t2q - mv2)**2 + mv2gv2 )/an2
         fch2 = abs(t2q)/( (t2q - mv2)**2 + mv2gv2 )/an2

         den = pch1*fch1 + pch2*fch2
         if (ich.eq.1) then
            anum = fch1
         else
            anum = fch2
         endif
         w = anum/den * w
      else

      endif      
****  end regulator      
      
c      ml152 = dot(p1-p5,p1-p5)
! now I have to decay p3 and p4, where pw is at rest
!     uhummm

      p3mod = 0.5d0*sqrt(lambda(ml342,ml32,ml42)/ml342)
      
      call new_boost(pw,p2,p2boo,1)

      p3(0) = sqrt(p3mod*p3mod + ml32)
      p2mod = p2boo(0) !sqrt(tridot(p2boo,p2boo))


      b3 = sqrt(1.d0 - ml32/p3(0)/p3(0))
c      b3  = p3mod/p3(0)
      ombet = 1.d0 - b3
      if (b3.eq.1.d0) then
         xx    =  -ml32/p3(0)/p3(0)
         ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .        - 0.078125d0 * xx*xx*xx)
         b3 = 1.d0 - ombet
      endif
cc
      anc = 1.d0/b3 * log((1.d0+b3)/ombet)

c      nc3 = 0
c 111  continue
c      nc3 = nc3 + 1
      
      call wraprng(csi,1)
      phi  = 2.d0*pi*csi(1)
      sphi = sin(phi)
      cphi = cos(phi)
      
      call wraprng(csi,1)
      c23 = 1.d0 - (1.d0 + b3)*exp(-b3*anc*csi(1))
      c23 = c23/b3
      s23 = sqrt(1.d0 - c23*c23)

      w = w * anc * (1.d0 - b3*c23) * 2.d0 * pi !/ nc3
      phsp = phsp * p3mod/ml34*0.25d0      
      
      p3(0) = sqrt(p3mod*p3mod + ml32)
      p3(1) = p3mod*s23*sphi
      p3(2) = p3mod*s23*cphi
      p3(3) = p3mod*c23

      call rot(-1,p2boo,p3,p3)

      p4(0)   =  ml34 - p3(0)
      p4(1:3) = -p3(1:3)            
      
      call new_boost(pw,p3,p3,-1)

! trying to implement pt cut here... !! mmhhh doesn't work...
c      pt32 = ptmp(1)*ptmp(1)+ptmp(2)*ptmp(2)
c      if (pt32.lt.ptmin*ptmin) goto 111

      call new_boost(pw,p4,p4,-1)
      
      qph(1,0) = p5(0)
      qph(1,1) = p5(1)
      qph(1,2) = p5(2)
      qph(1,3) = p5(3)
      return
      end
******************************************************
**********************************************************************
      subroutine phasespacephindz(p1,p2,p3,p4,p5,qph,m1,m2,m3,
     .     w,phsp,ie)
! for W
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3),p34(0:3)
      dimension p1boo(0:3),p2boo(0:3)
      dimension dir3(0:3),dir4(0:3),ptmp(0:3),vers(0:3),pcoll(0:3)
      double precision kone(0:3),k1(0:3),k2(0:3),pp(0:3),pw(0:3)
      double precision lambda
      double precision mq(-5:5),chq(-5:5),itre(-5:5)
      integer emissionpattern(1)
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (imaxph = 40)
      double precision qph(imaxph,0:3),oph(imaxph)
      common/nrandom/ir
      real*4 rnd(2),csi(1)
      common/radpattern/nph(4)      
      common/phindchannel/bbc,pw0,qmodc,c1qc,ich
! from shared.inc      
      character*1 boson
      common/vectorboson/boson
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     mq,chq,itre
      common/t_widths/gw,gz

      common/phspindcmn/ifirst
      data ifirst /0/
      
      common/exp_cuts/etamax,ptmin,cmin,cmax,ptminmiss,
     .     tmmin,mass2,ismear,irec,icutmt,ischemealpha

      external lambda
      double precision ml1,ml2,ml3,ml4,ml5,ml34,ml342
!! trying again

      common/mvgvbwaMCphind/ambwxs,agbwxs
      
      
      if (ifirst.eq.0) then
         print*,'Print in phind phasespace:'
         do j = 1,10
            print*,'tunare pch1 e pch2'
c            print*,'pt lepton cut hard wired'
         enddo
         ifirst = 1
      endif
      
      nphot = 1
      npart = nphot + 2
      nvar  = 3*nphot + 2

      ncharged = 4
      ie       = 0
      w        = 1.d0
      duepigrechi = (2.d0*pi)**(-nvar)

      phsp     = duepigrechi
      
      ei = p1(0)+p2(0)
      if (ei.lt.(m1+m2+m3)) then 
         phsp = 0.d0
         w    = 0.d0
         ie   = 1
         return
      endif
      
      ml1 = sqrt(dot(p1,p1))
      ml2 = 0.d0
      ml3 = m1
      ml4 = m2
      ml5 = m3

c      ml1 = 10.d0
c      ml2 = 10.d0
c      ml3 = 10.d0
c      ml4 = 20.d0
c      ml5 = 30.d0
      
*     * following the notes
*     * first sampling p3 + p4 as bw
      sl    = dot(p1+p2,p1+p2)
      ssl   = sqrt(sl)
      q2min = max((ml3+ml4)**2,tmmin*tmmin)
      q2max = (ssl - ml5)**2

      if (q2max.lt.q2min) then
         phsp = 0.d0
         w    = 0.d0
         ie   = 1
         return
      endif
      ml12 = ml1*ml1
      ml22 = ml2*ml2
      ml32 = ml3*ml3
      ml42 = ml4*ml4
      ml52 = ml5*ml5

      if (boson.eq.'W') then
         mv = mw
         gv = gw 
      else
         mv = mz
         gv = gz 
      endif
      ambwxs = mv
      agbwxs = gv

      mv2    = mv*mv
      gv2    = gv*gv
      mv2gv2 = mv2*gv2
      
** **
*     * **
      if (boson.eq.'W') then
         print*,'should not enter here!'
      else ! if boson eq Z
c         anbws1 = bwnorm(q2max,q2min)
         anbws1  = bwxsnorm(q2max,q2min)

         anooq2 = log(q2max/q2min)  
c         anooq2 = 1.d0/q2min - 1.d0/q2max
         
c         pbws = 0.5d0
         pbws = anbws1/(anbws1+anooq2)
         pgs  = 1.d0 - pbws
         
         call wraprng(csi,1)
         if (csi(1).le.pbws) then
c            call get_invm2(q2max,q2min,ml342,wm342)
            call bwxssample(q2max,q2min,ml342,wm342)
         else
            call wraprng(csi,1)
            ml342 = q2min * exp(anooq2*csi(1))
            wm342 = anooq2 * ml342
            
c            ml342 = q2min / (1.d0 - q2min*anooq2*csi(1))
         endif
c         den = pbws / anbws1/((ml342-mv2)**2+mv2gv2)
c     .       + pgs / anooq2/ml342
         den = pbws / anbws1/((ml342-mv2)**2+mv2gv2)*ml342
     .       + pgs / anooq2/ml342
c         den = pbws / anbws1/((ml342-mv2)**2+mv2gv2)
c     .       + pgs / anooq2/ml342/ml342

         w = w / den 
      endif
         
      ml34 = sqrt(ml342)

      qmod  = sqrt(lambda(sl,ml52,ml342)/sl)*0.5d0
      eq    = 0.5d0*(sl - ml52 + ml342)/ssl
      e5    = ssl - eq

c      eq = sqrt(qmod*qmod + ml342)
c      e5 = ssl - eq
      
      p1mod = sqrt(p1(0)*p1(0)-ml12) ! sqrt(tridot(p1,p1))
      p2mod = p2(0)                  !sqrt(tridot(p2,p2))      

c          get_tlimits(M2,ma2, mb2, m12,  m22, tp,tm)
*     I generate |t1Q| according to 1/(|t1Q| + ml52) in ich = 1 for diagrams 7 & 8
*     I generate  t2Q  according to spacelike bw in ich = 2 for diagrams 3 4 5 6
      tp1 = (p1(0)-eq)**2 - (p1mod - qmod)**2
      tm1 = (p1(0)-eq)**2 - (p1mod + qmod)**2
      tp2 = (p2(0)-eq)**2 - (p2mod - qmod)**2
      tm2 = (p2(0)-eq)**2 - (p2mod + qmod)**2         
            
      tmax = abs(tm1)
      tmin = abs(tp1)
      an1  = log(tmax+ml52) - log(tmin+ml52)
      tmax2 = abs(tm2)
      tmin2 = abs(tp2)
c      an2  = bwnorm(tp2,tm2)
      an2  = bwxsnorm(tmax2,tmin2)
      anoot = log(tmax2/tmin2)
      anbwt = an2

c      pbwt = 0.5d0
      pbwt = anbwt/(anbwt+anoot)
      pgt  = 1.d0 - pbwt

c      pch1 = 0.25d0 ! seems to be the optimum
c      pch1 = an1/(an1+an2)
c      pch2 = 1.d0 - pch1
      pch1 = 1.d0 / 3.d0 !2.d0/6.d0
      pch2 = 1.d0 - pch1
      ich  = 2
      call wraprng(csi,1)
      if (csi(1).le.pch1) ich = 1
      
      if (ich.eq.1) then
c         call get_tlimits(sl,ml12,ml22,ml342,ml52,tp,tm)
         tp = tp1
         tm = tm1
      else
         call get_tlimits(sl,ml22,ml12,ml342,ml52,tpa,tma)
         tp = tp2
         tm = tm2
      endif      

c      if (ich.eq.2.and.abs(tp).lt.1d-16) then
c         tpt = (eq - qmod) * (eq+qmod - 2.d0*p2(0))        
c         print*,tp,tpa,tpt,p2(0),eq,qmod
c      endif
      
      if (ich.eq.1) then
         an   = an1
         call wraprng(csi,1)
         abst = exp(an*csi(1))*(tmin+ml52) - ml52
         w    = w * an * (abst + ml52)
         t    = -abst
      else
         if (boson.eq.'W') then
         else ! if boson eq Z
            
            call wraprng(csi,1)
            if (csi(1).le.pbwt) then
               icht = 1
c               call get_invm2(tp,tm,t,wt)
               call bwxssample(tmax2,tmin2,t,wt)
               t = -t
               w = w * wt
            else
               icht = 2
               call wraprng(csi,1)
               at = tmin2 * exp(anoot*csi(1))
               wt = anoot*at
               
               t  = -at
               w = w * wt               
            endif
c            den = pbwt / anbwt/((t-mv2)**2+mv2gv2)*abs(t)
c     .          + pgt / anoot/abs(t)
c            w = w /den
         endif
      endif
            
      if (ich.eq.1) then
c         cq = (t - ml12 - ml342 + 2.d0*p1(0)*eq)/p1mod/qmod*0.5d0
         cq = 0.5d0*(t-tm)/p1mod/qmod - 1.d0
      else
c         cq = (t - ml22 - ml342 + 2.d0*p2(0)*eq)/p2mod/qmod*0.5d0
         cq = 0.5d0*(t-tm)/p2mod/qmod - 1.d0
         cq = -cq
      endif
      if (( abs(cq)-1.d0 ).ge.0.d0) then
         cq = cq/abs(cq)
      endif
      sq = sqrt(1.d0 - cq*cq)

      call wraprng(csi,1)
      phi  = 2.d0*pi*csi(1)
      sphi = sin(phi)
      cphi = cos(phi)
      w    = w * 2.d0*pi
      
      pw(0) = eq
      pw(1) = qmod*sq*sphi
      pw(2) = qmod*sq*cphi
      pw(3) = qmod*cq

cc      if (ich.eq.2) call rot(-1,p2,pw,pw)
      
      p5(0)   =  e5
      p5(1:3) = -pw(1:3)
      
      phsp = phsp * 0.25d0 / sqrt(lambda(sl,ml12,ml22)) ! for both channels

***   regulator for the two channels
      if (boson.eq.'W') then
      else
         t1q = dot(p1-pw,p1-pw)
         t2q = dot(p2-pw,p2-pw)
         fch1 = 1.d0/abs(t1q - ml52)/an1
         
         fsch21 = 1.d0/( (t2q - mv2)**2 + mv2gv2 )/anbwt * abs(t2q)
         fsch22 = 1.d0/abs(t2q)/anoot
c         fsch22 = 1.d0/t2q/t2q/anoot         
         
         den = pch1*fch1 + pch2* ( pbwt* fsch21 + pgt * fsch22 )
         if (ich.eq.1) then
            anum = fch1
         else
            if (icht.eq.1) anum = fsch21
            if (icht.eq.2) anum = fsch22
c            anum = ( pbwt* fsch21 + pgt * fsch22 )
         endif
         w = anum/den * w

      endif      
****  end regulator      
      
c      ml152 = dot(p1-p5,p1-p5)
! now I have to decay p3 and p4, where pw is at rest
!     uhummm

      p3mod = 0.5d0*sqrt(lambda(ml342,ml32,ml42)/ml342)
      
      call new_boost(pw,p2,p2boo,1)

      p3(0) = sqrt(p3mod*p3mod + ml32)
      p2mod = p2boo(0) !sqrt(tridot(p2boo,p2boo))


      b3 = sqrt(1.d0 - ml32/p3(0)/p3(0))
c      b3  = p3mod/p3(0)
      ombet = 1.d0 - b3
      if (b3.eq.1.d0) then
         xx    =  -ml32/p3(0)/p3(0)
         ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .        - 0.078125d0 * xx*xx*xx)
         b3 = 1.d0 - ombet
      endif
cc
      anc = 1.d0/b3 * log((1.d0+b3)/ombet)

c      nc3 = 0
c 111  continue
c      nc3 = nc3 + 1
      
      call wraprng(csi,1)
      phi  = 2.d0*pi*csi(1)
      sphi = sin(phi)
      cphi = cos(phi)
      
      call wraprng(csi,1)
      c23 = 1.d0 - (1.d0 + b3)*exp(-b3*anc*csi(1))
      c23 = c23/b3
      s23 = sqrt(1.d0 - c23*c23)

      w = w * anc * (1.d0 - b3*c23) * 2.d0 * pi !/ nc3
      phsp = phsp * p3mod/ml34*0.25d0      
      
      p3(0) = sqrt(p3mod*p3mod + ml32)
      p3(1) = p3mod*s23*sphi
      p3(2) = p3mod*s23*cphi
      p3(3) = p3mod*c23

!! double channel for mu+ or mu-      
      call wraprng(csi,1)
      if (csi(1).lt.0.5d0) then
         p3(1:3) = -p3(1:3)
      endif
      w = w * 2.d0/
     .  (1.d0 - b3*c23)/(1.d0/(1.d0 - b3*c23) + 1.d0/(1.d0 + b3*c23))
      
      call rot(-1,p2boo,p3,p3)

      p4(0)   =  ml34 - p3(0)
      p4(1:3) = -p3(1:3)            
      
      call new_boost(pw,p3,p3,-1)

! trying to implement pt cut here... !! mmhhh doesn't work...
c      pt32 = ptmp(1)*ptmp(1)+ptmp(2)*ptmp(2)
c      if (pt32.lt.ptmin*ptmin) goto 111

      call new_boost(pw,p4,p4,-1)

!     ! exchanging p3 and p3 to implement a double channel
c      call wraprng(csi,1)
c      w1 = 1.d0/dot(p2,p3)
c      w2 = 1.d0/dot(p2,p4)      
c      anum = w1
c      if (csi(1).lt.0.5d0) then
c         ptmp = p3
c         p3   = p4
c         p4   = ptmp
c         anum = w2
c      endif
c      w = w * anum/(w1+w2)*2.d0


    
      
      qph(1,0) = p5(0)
      qph(1,1) = p5(1)
      qph(1,2) = p5(2)
      qph(1,3) = p5(3)
      return
      end
******************************************************
      subroutine get_tlimits(M2,ma2,mb2,m12,m22,tp,tm)
      implicit double precision (a-h,m,l,o-z)
      external lambda
      tav = ma2 + m12 - (M2 + ma2 - mb2)*(M2 + m12 - m22)*0.5d0/M2
      ts  = sqrt(lambda(M2,m12,m22)*lambda(M2,ma2,mb2))*0.5d0/M2
      tp  = tav + ts
      tm  = tav - ts
      return
      end
******************************************************
      subroutine decay2b(rnd1,rnd2,mass,m1,m2,p1,p2,
     >     phsp,w,ie)
      implicit double precision (a-h,m,l,o-z)
      double precision rnd1,rnd2
      dimension p1(0:3),p2(0:3)
      parameter (pi = 3.1415926535897932384626433832795029D0)

      if (mass.lt.(m1+m2)) then
         ie   = 1
         w    = 0.d0
         phsp = 0.d0
         return
      endif

      ie   = 0
      w    = 1.d0
      phsp = 1.d0

      cth = 2.d0*rnd1 - 1.d0
      phi = 2.d0*pi*rnd2
      w   = 4.d0*pi*w
      
      sth  = sqrt(1.d0 - cth**2)
      cphi = cos(phi)
      sphi = sin(phi)
      
      sqla = sqrt(lambda(mass**2,m1**2,m2**2))
      phsp = phsp * sqla/8.d0/mass/mass

      pmod = sqla/2.d0/mass

      p1(0) = sqrt(pmod**2 + m1**2)
      p1(1) = pmod * sphi * sth
      p1(2) = pmod * cphi * sth
      p1(3) = pmod * cth

! to avoid round-offs when pmod is really small!!
      p2(0) = sqrt(pmod**2+m2**2)
      p2(1) = - p1(1)
      p2(2) = - p1(2)
      p2(3) = - p1(3)
      return
      end
***********************************************************
      function bwxsnorm(q2max,q2min)
      implicit double precision (a-h,o-z)
      common/mvgvbwaMCphind/am,ag
      common/bwxssamplecmn/aa,am2,amag,ooamag,ooaa,ifirst
      common/bwxssamplecmn2/emax,emin,F1min,F2min,d1,d2,ds
      data ifirst /0/
      if (ifirst.eq.0) then
         am2    = am*am
         amag   = am*ag
         ooamag = 1.d0/amag
         aa     = am2/amag
         ooaa   = 1.d0/aa
         ifirst = 1
      endif
         
      emin   = (q2min - am2)*ooamag
      emax   = (q2max - am2)*ooamag
         
      F1min = aa*atan(emin)
      F2min = 0.5d0 * log(emin*emin+1.d0)
         
      d1   = aa * atan(emax) - F1min
      d2   = 0.5d0 * log(emax*emax+1.d0) - F2min
      ds   = d1 + d2         
      
      bwxsnorm = ds
      
      return
      end
      
      subroutine bwxssample(q2max,q2min,s,w)
      implicit double precision (a-h,o-z)
      real*4 csi(1)
      common/mvgvbwaMCphind/am,ag
      common/bwxssamplecmn/aa,am2,amag,ooamag,ooaa,ifirst
      common/bwxssamplecmn2/emax,emin,F1min,F2min,d1,d2,ds
      data ifirst /0/
! from AcerMC !
      
      if (ifirst.eq.0) then
         am2    = am*am
         amag   = am*ag
         ooamag = 1.d0/amag
         aa     = am2/amag
         ooaa   = 1.d0/aa
         ifirst = 1
      endif
      
      emin   = (q2min - am2)*ooamag
      emax   = (q2max - am2)*ooamag
         
      F1min = aa*atan(emin)
      F2min = 0.5d0 * log(emin*emin+1.d0)
         
      d1   = aa * atan(emax) - F1min
      d2   = 0.5d0 * log(emax*emax+1.d0) - F2min
      ds   = d1 + d2
      
 333  call wraprng(csi,1)
      if (csi(1).lt.d2/ds) then
         call wraprng(csi,1) !! puttana eva!
         x =  d2*csi(1) + F2min
         e =  sqrt(exp(2.d0*x)-1.d0)         
      else
         call wraprng(csi,1)
         x = d1 * csi(1) + F1min
         e = tan(x*ooaa)
         if (e.lt.0.d0) then
            r = (aa+e)/aa
            call wraprng(csi,1)
            if (csi(1).ge.r) e = -e
            if ((e.gt.emax).or.(e.lt.emin)) goto 333
         endif
      endif
      s = amag * e + am2
      w = ds * ((s-am2)*(s-am2) + amag*amag )/s
      return
      end
********************
      subroutine get_invm2(q2max,q2min,minv2,w)
      implicit double precision (a-h,m,o-z)
      double precision minv2
      real*4 rnd(1)

      character*1 boson
      double precision mq(-5:5),chq(-5:5),itre(-5:5)
      
      common/vectorboson/boson
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     mq,chq,itre
      common/t_widths/gw,gz
      
      common/getinvmcommon2/am,ag,am2,ag2,amag,amsag,am2ag2,ifirst
      data ifirst/0/

      if (ifirst.eq.0) then
         if (boson.eq.'W') then
            am     = mw
            ag     = gw
         else
            am     = mz
            ag     = gz
         endif
         am2    = am*am
         ag2    = ag*ag
         amag   = am*ag
         amsag  = am/ag
         am2ag2 = am2*ag2

         ifirst = 1
      endif
      
      call wraprng(rnd,1)
      csi = 1.d0*rnd(1)
      
      iflat = 0
      if (iflat.eq.1) then
         minv2 = (q2max - q2min)*csi + q2min
         w = q2max - q2min
      else
         anbw  = bwnorm(q2max,q2min)         
         tmp   = atan(q2min/amag-amsag)
         minv2 = amag*( tan(amag*anbw*csi+tmp)+amsag )
         w     = anbw*((minv2-am2)**2+am2ag2)
      endif

      return
      end
********************************
      function bwnorm(q2max,q2min)
      implicit double precision (a-h,m,o-z)
      common/getinvmcommon2/am,ag,am2,ag2,amag,amsag,am2ag2,ifirst
      common/bwnormcommon/uozmzw,zmszw,ifirst2
      data ifirst2/0/

      character*1 boson
      double precision mq(-5:5),chq(-5:5),itre(-5:5)
      
      common/vectorboson/boson
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     mq,chq,itre
      common/t_widths/gw,gz
      
      if (ifirst2.eq.0) then

         if (boson.eq.'W') then
            am     = mw
            ag     = gw
         else
            am     = mz
            ag     = gz
         endif
         am2    = am*am
         ag2    = ag*ag
         amag   = am*ag
         amsag  = am/ag
         am2ag2 = am2*ag2
         
         uozmzw = 1.d0/am/ag
         zmszw  = am/ag
         ifirst2 = 1
      endif
      bwnorm=uozmzw*(atan(q2max*uozmzw-zmszw)
     >        - atan( q2min*uozmzw -zmszw) )
      return
      end
********************************
      subroutine bwovermcomplex(q2max,q2min,m2,w)
      implicit double precision (a-b,d-h,m,o-z)
      implicit double complex (c)
      common/getinvmcommon2/am,ag,am2,ag2,amag,amsag,am2ag2,ifirst
      double complex i_
      common/bwnormcommon/uozmzw,zmszw,ifirst2
      common/bwomcommoncomplex/i_,cz,csz,czs,cszs,ak,ifirst3
      common/bwomvarious/error
      data ifirst3/0/
      real*4 csi(1)
      
      if (ifirst3.eq.0) then
! to set-up the common!         
         call get_invm2(q2max,q2min,test,test2)
c         
         error = 1.d-5
         
         i_ = (0.d0,1.d0)
         cz  = am*am + i_*ag*am
         csz = sqrt(cz)
         czs  = dconjg(cz)
         cszs = dconjg(csz)
         ak = 2.d0*i_/(czs-cz)
         ifirst3 = 1
      endif

      qmax = sqrt(q2max)
      qmin = sqrt(q2min)
c      qmax = am + 3.d0
c      qmin = am - 3.d0

      fmin = funq(qmin)
      an   = ak * (funq(qmax) - fmin)

      call wraprng(csi,1)
      sleft = an * csi(1)/ak + fmin
      sleft = -sleft ! always negative  abs(sleft)
      
      xup   = qmax
      xdown = qmin
      x     = (xup+xdown)*0.5d0
      delta = 1.d0
      k = 0
      do while (delta.gt.error)
         k = k + 1
c         deltaf = abs(funq(x)) - sleft
         deltaf = -funq(x) - sleft ! funq always negative
         if (deltaf.gt.0.d0) then
            xup   = x
         else
            xdown = x
         endif
         x     = (xup+xdown)*0.5d0
         delta = abs(x-xup)/xdown ! relative error
      enddo

      m2 = x*x
      w  = an * x * ((m2-am2)**2 + am2ag2)
c      npoints = 10000
c      delta   = (qmax - qmin)/npoints
c      open(34,file='complex.txt',status='unknown')
c      x = qmin - delta
c      do k = 1,npoints
c         x = x+delta
c         cf = 1.d0/csz*log((csz+x)/(csz-x))
c         cf = 2.d0*i_*cf/(czs - cz)
c         write(34,*)x,Real(cf),ak*funq(x),dimag(cf)
c      enddo
c      close(34)
c      stop
      return
      end
**********************************************************************
      function funq(x)
      implicit double precision (a-b,d-h,o-z)
      implicit double complex (c)
      double complex i_
      common/bwomcommoncomplex/i_,cz,csz,czs,cszs,ak,ifirst3
      funq  = dImag(1.d0/csz*log((csz+x)/(csz-x)))
      return
      end
**********************************************************************
      subroutine phasespacephind_OLD(p1,p2,p3,p4,p5,qph,m1,m2,m3,
     .     w,phsp,ie)
! for W
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
      dimension dir3(0:3),dir4(0:3),ptmp(0:3),vers(0:3),pcoll(0:3)
      double precision kone(0:3),k1(0:3),k2(0:3),pp(0:3),pw(0:3)
      double precision lambda
      double precision mq(-5:5),chq(-5:5),itre(-5:5)
      character*1 boson
      integer emissionpattern(1)
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (imaxph = 40)
      double precision qph(imaxph,0:3),oph(imaxph)
      common/nrandom/ir
      real*4 rnd(2),csi(1)
      common/radpattern/nph(4)      
      common/phindchannel/bbc,pw0,qmodc,c1qc,ich
! from shared.inc      
! from shared.inc      
      common/exp_cuts/etamax,ptmin,cmin,cmax,ptminmiss,
     .     tmmin,mass2,ismear,irec,icutmt,ischemealpha
      common/vectorboson/boson
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     mq,chq,itre
      common/t_widths/gw,gz
      nphot = 1
      npart = nphot + 2
      nvar  = 3*nphot + 2

      ncharged = 4
      if (boson.eq.'W') ncharged = 3
      ie       = 0
      w        = 1.d0
      phsp     = 1.d0
      duepigrechi = (2.d0*pi)**(-nvar)

      ebper2 = p1(0)+p2(0)
      if (ebper2.lt.(m1+m2+m3)) then 
         phsp = 0.d0
         w    = 0.d0
         ie   = 1
         return
      endif

*****************
c      print*,''
c      print*,p1
c      print*,p2
c      print*,p1+p2
      
      ei   = p1(0)+p2(0)
      beta = (1.d0+m3/p1(0))*(1.d0-m3/p1(0))
      beta = sqrt(beta)

      dir3(0) = p1(0)
      dir3(1) = 0.d0
      dir3(2) = 0.d0
      dir3(3) = p1(0)*beta
      dir4(0) = dir3(0)
      ptmp(0) = dir3(0)
      pp(0)   = p1(0)+p2(0)
      beta5   = beta
      beta3   = sqrt(1.d0-m1**2/max(p1(0),2.d0*m1)**2)
      do k = 1,3
         pp(k)   = p1(k)+p2(k)
         dir4(k) = -dir3(k)
         ptmp(k) = -dir3(k)/beta*beta3
         dir3(k) = dir3(k)/beta*beta3
      enddo               
      beta1 = sqrt(tridot(p1,p1))/p1(0)
************************
      vm = mw
      vg = gw
      if (boson.eq.'Z') then
         vm = mz
         vg = gz
      endif

      ssl = sqrt(dot(p1+p2,p1+p2))
      
      q2max = (ei-m3)**2
      q2min = (m1+m2)**2

!!!      
      q2minmass = (m1+m2)**2
      q2min = (max(tmmin,ptmin+ptminmiss))**2
      q2min = max(q2min,q2minmass)
!!!!  

      if (q2max.lt.q2min) then
         phsp = 0.d0
         w    = 0.d0
         ie   = 1
         return
      endif
      
      gm = vm*vg
      A     =  atan(q2min/gm - vm/vg)
      anbw  = (atan(q2max/gm - vm/vg) - A)/gm
      
      pflatc = 0.5d0
      anc3 = 1.d0/beta3*log((1.d0+beta3)/(1.d0-beta3))
      anc5 = 1.d0/beta5*log((1.d0+beta5)/(1.d0-beta5))
      
      e5min = m3
      e5max = ei/2.d0

      gam = 0.91d0
      angam =  1.d0/(1.d0-gam)*(e5max-e5min)**(1.d0-gam)
      pgam = 0.2d0

      pich1 = 0.5d0 
      pich2 = 1.d0 - pich1
      pich3 = 0.d0 !!1.d0 - pich1 - pich2

      call wraprng(csi,1)
      ich = 2
      if (csi(1).le.pich1) ich = 1
      if (csi(1).gt.pich1.and.csi(1).le.(pich1+pich2)) ich = 2
!!! ONLY 2 CHANNELS !!!

      if (ich.eq.1) then
!     phase space splitting
         call wraprng(csi,1)
         q2 = gm*tan(gm*anbw*csi(1)+A) + vm*vm
         w = w * anbw * ((q2 - vm*vm)**2 + gm*gm)
         
         sqla = sqrt(lambda(ei**2,q2,m3**2))
         pmod = sqla/2.d0/ei
         p5(0) = sqrt(pmod**2 + m3**2)

         dir4(0) = p5(0)
         dir4(1) = 0.d0
         dir4(2) = 0.d0
         dir4(3) = -pmod

         call wraprng(rnd,2)
         call collinear(pflatc,rnd(1),rnd(2),dir4,vers,wverg)
         w = w * wverg
         c = vers(3)
         s = sqrt(abs(1.d0 - c*c))
         if (s.gt.0.d0) then
            sphi = vers(1)/s
            cphi = vers(2)/s
         else
            sphi = 0.d0
            cphi = 1.d0
         endif
         p5(1) = pmod * sphi * s
         p5(2) = pmod * cphi * s
         p5(3) = pmod * c
         pw(0) = ei - p5(0)
         pw(1) = - p5(1)
         pw(2) = - p5(2)
         pw(3) = - p5(3)
         phsp1 = sqla/8.d0/ei**2
         
         call wraprng(rnd,2)
         c = 2.d0*rnd(1) - 1.d0
         s = sqrt(abs(1.d0-c**2))
         w = w * 2.d0
         phi = 2.d0*pi*rnd(2)
         w   = w *2.d0*pi     
         sphi = sin(phi)
         cphi = cos(phi)
         
         sqla = sqrt(lambda(q2,m1**2,m2**2))
         p3mod = sqla/2.d0/sqrt(q2)
         p3(0) = sqrt(p3mod**2 + m1**2)
         p3(1) = p3mod * sphi * s
         p3(2) = p3mod * cphi * s
         p3(3) = p3mod * c
         p4(0) = sqrt(q2) - p3(0)
         p4(1) = - p3(1)
         p4(2) = - p3(2)
         p4(3) = - p3(3)
         
         phsp2 = sqla/8.d0/q2
         
         call new_boost(pw,p3,p3,-1)
         call new_boost(pw,p4,p4,-1)

         phsp = phsp1*phsp2/(2.d0*pi)**5

      elseif (ich.eq.2) then
*-----------
         pflq2 = 0.5d0
         call wraprng(csi,1)
         if (csi(1).le.pflq2) then
            ich2 = 1
            call wraprng(csi,1)
            q2 = (q2max-q2min)*csi(1) + q2min
            w = w * (q2max-q2min)
         else
            ich2 = 2
            call wraprng(csi,1)
            q2 = gm*tan(gm*anbw*csi(1)+A) + vm*vm
            w = w * anbw * ((q2 - vm*vm)**2 + gm*gm)
         endif

! consider ich2....

         sqla = sqrt(lambda(ei**2,q2,m3**2))
         pmod = sqla/2.d0/ei
         p5(0) = sqrt(pmod**2 + m3**2)

         ab = 2.d0*p1(0)*p5(0)*beta1*pmod/p5(0)
         bb = -2.d0*p1(0)*p5(0) - vm*vm + dot(p1,p1) + m3*m3
         cb = gm
         ak = atan((-ab+bb)/cb)
         ancbwt =  (atan((ab+bb)/cb) - ak)/ab/cb
         call wraprng(csi,1)
         c = (cb*tan(ab*cb*ancbwt*csi(1)+ak) - bb)/ab
         w = w * ancbwt*((ab*c+bb)**2 + cb*cb)
         call wraprng(csi,1)
         phi = 2.d0*pi*csi(1)
         w   = w * 2.d0*pi

         s = sqrt(abs(1.d0 - c*c))
         sphi = sin(phi)
         cphi = cos(phi)
         
         p5(1) = pmod * sphi * s
         p5(2) = pmod * cphi * s
         p5(3) = pmod * c
         pw(0) = ei - p5(0)
         pw(1) = - p5(1)
         pw(2) = - p5(2)
         pw(3) = - p5(3)
         qmod = sqrt(tridot(pw,pw))

         phsp1 = sqla/8.d0/ei/ei

         call wraprng(rnd,2)
         if (ich2.eq.2) then
            call collinear(1.d0,rnd(1),rnd(2),ptmp,vers,wverg)
            w = w * wverg
         else
            call collinear(0.d0,rnd(1),rnd(2),ptmp,vers,wverg)
            w = w * wverg
         endif

         c1q  = tridot(vers,pw)/qmod
         
         b2 = qmod*qmod
         ak = 2.d0*qmod*c1q

         U  = pw(0)**2 + m1**2 - m2**2 - b2
         aa = 4.d0*pw(0)**2 - ak*ak
         bb = -2.d0*U*ak
         cc = 4.d0*pw(0)**2*m1*m1 - U*U

         arg = bb*bb - 4.d0*aa*cc
         if (arg.lt.0.d0) then
            w = 0.d0
            phsp = 0.d0
            ie = 1
            return
         endif

         p3mod1 = (-bb + sqrt(arg))/2.d0/aa
         p3mod2 = (-bb - sqrt(arg))/2.d0/aa

         p3jac = 1.d0
         if (p3mod1.ge.0.d0.and.p3mod2.lt.0.d0) then
            p3mod = p3mod1
         elseif(p3mod1.lt.0.d0.and.p3mod2.ge.0.d0) then
            p3mod = p3mod2
         elseif(p3mod1.ge.0.d0.and.p3mod2.ge.0.d0) then
            p3jac = 2.d0
            p3mod = p3mod1
            call wraprng(csi,1)
            if (csi(1).le.0.5d0) then
               p3mod = p3mod2
            endif
         else
            w = 0.d0
            phsp = 0.d0
            ie = 1
            return
         endif
         w = w * p3jac

         p3(0) = sqrt(p3mod**2+m1*m1)
         p3(1) = p3mod*vers(1)
         p3(2) = p3mod*vers(2)
         p3(3) = p3mod*vers(3)
         p4(0) = pw(0) - p3(0)
         p4(1) = pw(1) - p3(1)
         p4(2) = pw(2) - p3(2)
         p4(3) = pw(3) - p3(3)

         bbb = p3mod/p3(0)
         phsp2 = bbb*p3mod/abs(bbb*pw(0)-qmod*c1q)/4.d0

         bbc = bbb
         pw0 = pw(0)
         qmodc = qmod
         c1qc = c1q

         phsp = phsp1*phsp2/(2.d0*pi)**5
      else
! THIS CHANNEL IS NO MORE ACTIVE.....
         call wraprng(csi,1)
         q2 = gm*tan(gm*anbw*csi(1)+A) + vm*vm
         w = w * anbw * ((q2 - vm*vm)**2 + gm*gm)

         sqla = sqrt(lambda(ei**2,q2,m3**2))
         pmod = sqla/2.d0/ei
         p5(0) = sqrt(pmod**2 + m3**2)

         ab = 2.d0*p1(0)*p5(0)*beta1*pmod/p5(0)
         bb = -2.d0*p1(0)*p5(0) - vm*vm + dot(p1,p1) + m3*m3
         cb = gm
         ak = atan((-ab+bb)/cb)
         ancbwt =  (atan((ab+bb)/cb) - ak)/ab/cb

         call wraprng(csi,1)
         c = (cb*tan(ab*cb*ancbwt*csi(1)+ak) - bb)/ab
         w = w * ancbwt*((ab*c+bb)**2 + cb*cb)

         call wraprng(csi,1)
         phi = 2.d0*pi*csi(1)
         w   = w * 2.d0*pi

         s = sqrt(abs(1.d0 - c*c))
         sphi = sin(phi)
         cphi = cos(phi)
         
         p5(1) = pmod * sphi * s
         p5(2) = pmod * cphi * s
         p5(3) = pmod * c
         pw(0) = ei - p5(0)
         pw(1) = - p5(1)
         pw(2) = - p5(2)
         pw(3) = - p5(3)
         phsp1 = sqla/8.d0/ei**2
         
         call wraprng(rnd,2)
         c = 2.d0*rnd(1) - 1.d0
         s = sqrt(abs(1.d0-c**2))
         w = w * 2.d0
         phi = 2.d0*pi*rnd(2)
         w   = w *2.d0*pi     
         sphi = sin(phi)
         cphi = cos(phi)
         
         sqla = sqrt(lambda(q2,m1**2,m2**2))
         p3mod = sqla/2.d0/sqrt(q2)
         p3(0) = sqrt(p3mod**2 + m1**2)
         p3(1) = p3mod * sphi * s
         p3(2) = p3mod * cphi * s
         p3(3) = p3mod * c
         p4(0) = sqrt(q2) - p3(0)
         p4(1) = - p3(1)
         p4(2) = - p3(2)
         p4(3) = - p3(3)
         
         phsp2 = sqla/8.d0/q2
         
         call new_boost(pw,p3,p3,-1)
         call new_boost(pw,p4,p4,-1)

         phsp = phsp1*phsp2/(2.d0*pi)**5

      endif

      c5 = p5(3)/sqrt(tridot(p5,p5))

      ab = 2.d0*p1(0)*(1.d0-beta1*beta5*c5)
      bb = vm**2-dot(p1,p1)-m3*m3
      cb = vm*vg
      ak = atan((ab*e5min+bb)/cb)
      anbwt = (atan((ab*e5max+bb)/cb) - ak)/ab/cb

      do k = 0,3
         ptmp(k) = p3(k)+p4(k)
      enddo
      q2 = dot(ptmp,ptmp)

      cg3 = tridot(p3,p2)/sqrt(tridot(p3,p3))/p2(0)
      cg5 = tridot(p5,p2)/sqrt(tridot(p5,p5))/p2(0)

      rcd  = 1.d0/anc5/(1.d0-beta5*cg5)*(1.d0-pflatc)+
     .     1.d0/2.d0*pflatc
      rcmu = 1.d0/anc3/(1.d0-beta3*cg3)*(1.d0-pflatc)+
     .     1.d0/2.d0*pflatc

      rcqw = 1.d0/dot(p3,ptmp)*p3(0)*ptmp(0)

      rcd  = 1.d0/(1.d0-beta5*cg5)*(1.d0-pflatc)+
     .     1.d0*pflatc


      rcmu = 1.d0/(1.d0-beta3*cg3)*(1.d0-pflatc)+
     .     1.d0*pflatc

      rcmu = 1.d0/(1.d0-beta3*cg3)


c      rcd = rcd * anc5
c      rcmu = rcmu * anc3

      rbw = 1.d0/anbw/((q2-vm*vm)**2+gm*gm)
      rbwt = 1.d0/anbwt/((ab*p5(0)+bb)**2 + cb*cb) 
      rbw = rbw*anbw
      rbwt = rbwt*anbwt

      ab = 2.d0*p1(0)*p5(0)*beta1*sqrt(tridot(p5,p5))/p5(0)
      bb = -2.d0*p1(0)*p5(0) - vm*vm + dot(p1,p1) + m3*m3
      cb = gm
      ak = atan((-ab+bb)/cb)
      ancbwt =  (atan((ab+bb)/cb) - ak)/ab/cb
      c5 = p5(3)/sqrt(tridot(p5,p5))
      
      rcbwt = 1.d0/ancbwt/((ab*c5+bb)**2+cb*cb)
      rcbwt = rcbwt * ancbwt
      
      reg1 = rbw*rcd
      reg2 = rcbwt*(rbw+rcmu)

      if (ich.eq.1) w = w * reg1/(pich1*reg1 + pich2*reg2 + pich3*reg3)
      if (ich.eq.2) w = w * reg2/(pich1*reg1 + pich2*reg2 + pich3*reg3)

ccc      print*,p5(0)-(ei**2+m3**2-q2)/2.d0/ei

      
      qph(1,0) = p5(0)
      qph(1,1) = p5(1)
      qph(1,2) = p5(2)
      qph(1,3) = p5(3)

      return
      end

      function elmat2phindw(iii,p1o,p2o,p3o,p4o,ko,Qu,Qd,Qw,Qe)
!! u dbar --> W+ --> e+ nue gamma
!! p1 p2             p4 p3  k
!! attention: it is different from my usual convention!!
      implicit none
      common/fordebugging/idebugging
      integer i,idebugging,iii
      double precision me,m1,m2,ch1,ch2,chfs,mfs
      integer lepton,imirror
      double precision dot,elmat2phindw,q2,mup
      double precision pr1,pr2,pr3,pr4,pr5,pr6,pr7,pr8,pr9,pr10
      complex*16 elmat2,ecomplex_,e_
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),k(0:3),
     > ptmp(0:3), q(0:3)
      double precision p1o(0:3),p2o(0:3),p3o(0:3),p4o(0:3),ko(0:3)
      double precision kp1,kp2,kp3,kp4,kp12,kp22,kp32,kp42
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      double precision p1p22,p1p32,p1p42,p2p32,p2p42,p3p42
      double precision p1p1,p2p2,p3p3,p4p4
      double precision me2,mup2,md2
      double precision me4,mup4,md4
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
      complex*16 eps1234
      double precision gw,gz,el2,Qu,Qd,Qe,Qu2,Qd2,Qe2,gwrun,gamw
!      complex*16 Qw,Qw2       ! due to Baur-Zeppenfeld recipe
      double precision Qw,Qw2
      double precision alpha,alpha_gf,s2tw,gf,pi,convfac
      common/t_widths/gw,gz
      double precision dum1(-5:5),dum2(-5:5),dum3(-5:5)
      common/masseshor/meorig,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     dum1,dum2,dum3
      common/processhor/m1,m2,ch1,ch2,chfs,mfs,lepton,imirror
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/chargesmasses_cm/charges,amasses
      common/propagateurs/pr1,pr2,pr3,pr4,pr5,pr6,pr7,pr8,pr9,pr10
      data pr1,pr2,pr3,pr4,pr5,pr6,pr7,pr8,pr9,pr10 /10*0.d0/
      data im /(0.d0,1.d0)/

      do i = 0,3
         if (iii.eq.1) then
            p1(i) =  p1o(i)
            p2(i) = -p2o(i)
         else
            p1(i) = -p1o(i)
            p2(i) =  p2o(i)
         endif
         p3(i) =  p3o(i)
         p4(i) =  p4o(i)
         k(i)  = -ko(i)
      enddo

      i_ = im
      gamw = gw/mw

!      Qw =  1.d0                    ! Lopez-Castro recipe
!      Qw =  1.d0 * (1.d0 + im*gamw) ! Baur - Zeppenfeld recipe

c      Qe = 0.d0
c      qd = -1.d0
c      qu = 1.d0
c      Qd = 0.d0
c      Qw  = 0.d0

      Qe2 = Qe*Qe
      Qd2 = Qd*Qd
      Qw2 = Qw*Qw
      Qu2 = Qu*Qu

      el2 = alpha*4.d0*pi
      
      kp1 = dot(k,p1)
      kp2 = dot(k,p2)
      kp3 = dot(k,p3)
      kp4 = dot(k,p4)
      
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

      mup = amasses(1)
      me  = amasses(3)
      mup2 = mup*mup
      md2  = dot(p2,p2)
      me2  = me*me

      mup4 = mup2*mup2
      md4  = md2*md2
      me4  = me2*me2

*      mwm2 = mw**(-2)
*      mwm4 = mw**(-4)
*      mwm6 = mw**(-6)
*      mwm8 = mw**(-8)

      mw2  = mw*mw
      mwc2  = mw2 - im*gw*mw !* 0.d0
      mwc2d = mw2 + im*gw*mw !* 0.d0

      do i = 0,3
         ptmp(i) = p1(i) + p2(i)
      enddo      
      q2 = dot(ptmp,ptmp)
      gwrun = mw2*gamw           ! Lopez-Castro recipe
!      gwrun = q2*gamw            ! Baur - Zeppenfeld recipe
      s  = 1.d0/(q2 - mwc2)

c      print*,''
c      print*,q2,mwc2
      
cc      if (idebugging.eq.1) print*,'sq.f 1  q2 =',q2

      do i = 0,3
         ptmp(i) = p4(i) + p3(i)
      enddo
      q2 = dot(ptmp,ptmp)
      gwrun = mw2*gamw           ! Lopez-Castro recipe
!      gwrun = q2*gamw            ! Baur - Zeppenfeld recipe
      sf  = 1.d0/(q2 - mwc2)

c      print*,q2

      
      prop12  = s
      prop12d = conjg(s)
      prop34  = sf
      prop34d = conjg(sf)

      mw2m1 = 1.d0/mw2 ! Ward identities-violating or null width
      mw2m1 = 1.d0/(mw2 - im*gw*mw)     ! Lopez-Castro recipe
!      mw2m1 = 1.d0/mw2 * (1.d0+im*gamw) ! Baur - Zeppenfeld recipe

      mw2m2 = mw2m1*mw2m1
      mw2m3 = mw2m2*mw2m1
      mw2m4 = mw2m3*mw2m1

      mwc2m1  = 1.d0/mwc2
      mwc2dm1 = 1.d0/mwc2d

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

      elmat2 = -elmat2 !! 5.68 Peskin-Schroeder

      elmat2phindw = elmat2 * el2**3/s2tw**2 / 1024d0
      return
      end


      function elmat2phindz(iii,p1o,p2o,p3o,p4o,ko,Qq,Qf)
      implicit none
      common/fordebugging/idebugging
      integer i,idebugging,iii
      double precision dot,elmat2phindz,q2
      complex*16 elmat2,ecomplex_,e_
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),k(0:3),
     > ptmp(0:3), q(0:3)
      double precision p1o(0:3),p2o(0:3),p3o(0:3),p4o(0:3),ko(0:3)
      double precision kp1,kp2,kp3,kp4,kp12,kp22,kp32,kp42
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      double precision p1p22,p1p32,p1p42,p2p32,p2p42,p3p42
      double precision p1p1,p2p2,p3p3,p4p4
      double precision mf,mf2,mf4,mq,mq2,mq4,s,mzm2,mzm4
      double precision Qq,Qf,Qq2,Qf2,Qq4,Qf4
      double precision kp1m1,kp2m1,kp3m1,kp4m1,kp1m2,kp2m2,kp3m2,kp4m2
      complex*16 im
      double precision charges(4),amasses(4)
      complex*16 eps1234
      double precision gw,gz,s12mmz2,s34mmz2,dg12,dg122,dg34,dg342
      double precision dz12m2,dz34m2
      double precision e2,e4,gs2sctw2,gs2sctw4,g2,gzmz,gzmz2
      double precision gaf,gvf,gaq,gvq,gaf2,gvf2,gaq2,gvq2,gfp
      double precision alpha,alpha_gf,s2tw,gf,pi,convfac
      double precision m1,m2,ch1,ch2,chfs,mfs
      double precision meorig,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt
      integer lepton,imirror
      common/t_widths/gw,gz
      double precision dum1(-5:5),dum2(-5:5),dum3(-5:5)
      common/masseshor/meorig,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     dum1,dum2,dum3
      common/processhor/m1,m2,ch1,ch2,chfs,mfs,lepton,imirror
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/chargesmasses_cm/charges,amasses
      common/imaginaryzphind/im
      data im /(0.d0,1.d0)/

      do i = 0,3
         if (iii.eq.1) then
            p1(i) =  p1o(i)
            p2(i) =  -p2o(i)
         else
            p1(i) = -p1o(i)
            p2(i) =  p2o(i)
         endif
         p3(i) =  p3o(i)
         p4(i) =  p4o(i)
         k(i)  =  -ko(i)
      enddo


      e2 = alpha*4.d0*pi
      e4 = e2*e2
      g2 = e2/s2tw
c      gfp = sqrt(2.d0)*e2/s2tw/8.d0/mw/mw
c      g2 = gfp/sqrt(2.d0)*mw*mw*8.d0
      gs2sctw2 = g2/4.d0/(1.d0-s2tw)
      gs2sctw4 = gs2sctw2*gs2sctw2

      gaf  = -0.5d0 ! electron weak isospin!
      gvf  =  gaf - 2.d0*Qf*s2tw
      gaf2 =  gaf*gaf
      gvf2 =  gvf*gvf
      if (Qq.gt.0) then
         gaq = 0.5d0
         gvq = gaq - 2.d0*Qq*s2tw
      else
         gaq = -0.5d0
         gvq =  gaq - 2.d0*Qq*s2tw
      endif
      gaq2 = gaq*gaq
      gvq2 = gvq*gvq

      Qf2 =Qf*Qf
      Qf4 =Qf2*Qf2
      Qq2 =Qq*Qq
      Qq4 =Qq2*Qq2

      kp1 = dot(k,p1)
      kp2 = dot(k,p2)
      kp3 = dot(k,p3)
      kp4 = dot(k,p4)
      
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

      mf = amasses(3)
      mf2 = mf*mf
      mf4 = mf2*mf2
      mq = amasses(1)
      mq2 = mq*mq
      mq4 = mq2*mq2

      mzm2 = 1.d0/mz/mz
      mzm4 = mzm2*mzm2
      gzmz = gz*mz
      gzmz2= gzmz*gzmz

      do i = 0,3
         ptmp(i) = p1(i) + p2(i)
      enddo      
      q2      = dot(ptmp,ptmp)
      dg12    = 1.d0/q2
      dg122   = dg12*dg12
      s12mmz2 = q2 - mz*mz
      dz12m2  = 1.d0/((q2-mz**2)**2 + gzmz2)

      do i = 0,3
         ptmp(i) = p4(i) + p3(i)
      enddo
      q2      = dot(ptmp,ptmp)

      dg34    = 1.d0/q2
      dg342   = dg34*dg34
      s34mmz2 = q2 - mz*mz
      dz34m2  = 1.d0/((q2-mz**2)**2 + gzmz2)

      kp1m1 = 1.d0/kp1
      kp2m1 = 1.d0/kp2
      kp3m1 = 1.d0/kp3
      kp4m1 = 1.d0/kp4

      kp1m2 = kp1m1*kp1m1
      kp2m2 = kp2m1*kp2m1
      kp3m2 = kp3m1*kp3m1
      kp4m2 = kp4m1*kp4m1

****      eps1234 = e_(p1,p2,p3,p4)
      eps1234 = (0.d0,0.d0)
      elmat2 = (0.d0,0.d0)
      include 'form/zrad.f'
ccc      include 'form/phindz.f' ! used to select diagrams....
      elmat2 = -elmat2 !! 5.68 Peskin-Schroeder
      elmat2phindz = elmat2
      return
      end
