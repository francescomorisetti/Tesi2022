**** SAMPLING
      subroutine multiplicityv2(eps,ecms,n,w)
      implicit double precision (a-h,m,l,o-z)
      integer lepton
      character*1 boson
      real*4 csi(1)
      character*6 ord
      double precision mq(-5:5),chq(-5:5),itre(-5:5)
      common/vectorboson/boson
      common/qedORDER/ord
      common/nphot_mode/nphotmode
!! SIGNALS  n = -1003  ---> 3 body
!! SIGNALS  n = -1002  ---> 2 body
*from shared.inc
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt
     .     ,mq,chq,itre
      common/processhor/m1,m2,ch1,ch2,chfs,mfs,lepton,imirror
      common/partons/ipart1,ipart2

      if (ipart1.eq.0.and.ipart2.eq.0) then
         w = 1.d0
         n = 0
         return
      endif

      if (nphotmode.ge.0) then
         w = 1.d0
         n = nphotmode
         return
      endif
      if (ord.eq.'born') then
         w = 1.d0
         n = 0
         return
      endif
      if (ord.eq.'alpha') then
         sffa = 0.3d0
c         if (ecms.gt.500d0)  sffa = 0.05d0
c         if (ecms.gt.1000d0) sffa = 0.01d0
         p0 = sffa
         p1 = 1.d0 - sffa      

         p1003 = 0d0
         p1002 = 0d0         

         p0 = 4.d0
         p1 = 6.d0

         ptot  = p0+p1+p1003+p1002
         p0    = p0/ptot
         p1    = p1/ptot
         p1002 = p1002/ptot
         p1003 = p1003/ptot

         call wraprng(csi,1)
         x = csi(1)
         if (x.le.p0) then 
            n = 0
            w = 1.d0 / p0
         elseif (x.le.p0+p1) then
            n = 1
            w = 1.d0 / p1
         elseif (x.le.p0+p1+p1002) then
            n = -1002
            w = 1.d0/p1002
         else
            n = -1003
            w = 1.d0/p1003
         endif
         return
      endif
!-------------------------------------
! faster (and more elegant!) poissonian distribution
      if (ord.eq.'exp') then
         pusual = 1.d0
         p1002  = 0.d0
         p1003  = 0.d0

         ptot   = pusual+p1002+p1003
         pusual = pusual/ptot
         p1002  = p1002/ptot
         p1003  = p1003/ptot
         call wraprng(csi,1)
         x = csi(1)
         if (x.le.pusual) then
            scale  = ecms**2
            aieps  = -2.d0*log(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
            lcoll1 = log(scale/m1**2)  - 1.d0
            lcoll2 = log(scale/m2**2)  - 1.d0
            lcoll3 = log(scale/mfs**2) - 1.d0
            lcoll4 = 0.d0
            if (boson.eq.'Z') lcoll4 = log(scale/mfs**2) - 1.d0
            lcoll  = 4.d0/9.d0*lcoll1+1.d0/9.d0*lcoll2+lcoll3+lcoll4
            arg    = alpha/2.d0/pi * lcoll * aieps

            ffs  = exp(-arg)
            call wraprng(csi,1)
            cs   = 1.d0 * csi(1)
            ptot = ffs
            pn   = ptot
            n    = 0
            w    = 1.d0/pn
            do while(cs.gt.ptot)
               n    = n + 1
               pn   = pn*arg/n
               ptot = ptot + pn
            enddo
            w = 1.d0/pn/pusual
         elseif (x.le.pusual+p1002) then
            n = -1002
            w = 1.d0/p1002
         else
            n = -1003
            w = 1.d0/p1003
         endif
      endif
      return
      end
**
      subroutine multiplicity(eps,ecms,n,w)
      implicit double precision (a-h,m,l,o-z)
      double precision mq(-5:5),chq(-5:5),itre(-5:5)
      integer lepton
      real*4 csi(1)
      character*6 ord
      common/qedORDER/ord
      common/nphot_mode/nphotmode
*from shared.inc
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt
     .     ,mq,chq,itre
      common/processhor/m1,m2,ch1,ch2,chfs,mfs,lepton,imirror
      if (nphotmode.ge.0) then
         w = 1.d0
         n = nphotmode
         return
      endif

      if (ord.eq.'born') then
         w = 1.d0
         n = 0
         return
      endif

      if (ord.eq.'alpha') then
         sffa = 0.3d0
         if (ecms.gt.500d0) sffa = 0.05d0
         if (ecms.gt.1000d0) sffa = 0.01d0
         p0 = sffa
         p1 = 1.d0 - sffa      
         ptot = 1.d0
         call wraprng(csi,1)
         if (csi(1).lt.p0) then 
            n = 0
            w = ptot / p0
         else
            n = 1
            w = ptot / p1
         endif
         return
      endif

      scale = ecms**2
      aieps = -2.d0*log(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
      lcoll1 = log(scale/m1**2) - 1
      lcoll2 = log(scale/m2**2) - 1
      lcoll3 = log(scale/mfs**2) - 1

      lcoll = lcoll1+lcoll2+lcoll3

      arg = alpha/2.d0/pi * lcoll * aieps
! faster (and more elegant!) poissonian distribution
      if (ord.eq.'exp') then
         ffs = exp(-arg)
         call wraprng(csi,1)
         cs = 1.d0 * csi(1)
         ptot = ffs
         pn   = ptot
         n = 0
         w = 1.d0/pn
         do while(cs.gt.ptot)
            n = n + 1
            pn = pn*arg/n
            ptot = ptot + pn
         enddo
         w = 1.d0/pn
      endif
      return
      end
*
      subroutine get_pattern(ncharged,n,ep)
      integer n,ep(n),ncharged
      real*4 csi(1)
      do k = 1,n
         call wraprng(csi,1)
         ep(k) = 1.d0*ncharged*csi(1) + 1
         if (ep(k).gt.ncharged) ep(k) = ncharged ! to avoid round-offs...
      enddo
      return
      end
*
      subroutine get_pattern_v2(ncharged,n,ep,w)
      implicit double precision (w)
      integer n,ep(n),ncharged,i
      dimension wi(ncharged)
      real*4 csi(1)
      w = 1.d0
      wi(1) = 4.d0/9.d0
      wi(2) = 1.d0/9.d0
      wi(3) = 1.d0
      wtot  = wi(1)+wi(2)+wi(3)
      wi(1) = wi(1)/wtot
      wi(2) = wi(2)/wtot
      wi(3) = wi(3)/wtot
      do k = 1,n
         i = 1
         wc = wi(1)
         call wraprng(csi,1)
         do while(csi(1).gt.wc)
            i = i + 1
            wc = wc + wi(i)
         enddo
         w = w / wi(i)
         ep(k) = i
         if (ep(k).gt.ncharged) ep(k) = ncharged ! to avoid round-offs...
      enddo
      return
      end
*** PHASE SPACE
      subroutine phasespace(p1,p2,p3,p4,qph,nphot,m1,m2,esoft,
     .     w,phsp,ie)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension dir3(0:3),dir4(0:3),ptmp(0:3),vers(0:3),pcoll(0:3)
      double precision kone(0:3),k1(0:3),k2(0:3),pp(0:3)
      double precision lambda
      double precision mq(-5:5),chq(-5:5),itre(-5:5)
      character*1 boson
      integer emissionpattern(nphot)
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (imaxph = 40)
      double precision qph(imaxph,0:3),oph(imaxph)
      common/cosgp3/cgp3
      common/nrandom/ir
      real*4 rnd(2),csi(1)
      common/radpattern/nph(4)      
      common/masseinitial/muno,mdue
      common/partialweights/wlepang,wphoten,wphotang,reg
!     from shared.inc
      common/partons/ipart1,ipart2      
      common/vectorboson/boson
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     mq,chq,itre
      common/t_widths/gw,gz
      npart = nphot + 2
      nvar  = 3*nphot + 2
      ncharged = 4
      if (boson.eq.'W') ncharged = 3
      ie       = 0
      w        = 1.d0
      wlepang  = 1.d0
      wphoten  = 1.d0
      wphotang = 1.d0
      reg      = 1.d0
      phsp     = 1.d0
      duepigrechi = (2.d0*pi)**(-nvar)
      iwider = 0
      ntot   = nphot
      egmin  = esoft
      do j = 0,3
         do k = 1,imaxph
            qph(k,j) = 0.d0
         enddo
      enddo

      eb     = p1(0)
      ebper2 = 2.d0*eb

      if (ebper2.lt.(m1+m2+nphot*esoft)) then 
         phsp = 0.d0
         w    = 0.d0
         ie   = 1
         return
      endif
*****************
      ir = 1
      call wraprng(rnd,2)

      if (nphot.eq.0) then
         mass = p1(0)+p2(0)
         sqla = sqrt(lambda(mass**2,m1**2,m2**2))
         phsp = phsp * sqla/8.d0/mass**2
         call get_cos_fer(rnd(1),nphot,c,wcos)
         w = w * wcos
         s = sqrt(1.d0 - c**2)
         call get_phi(rnd(2),phi,wph)
         w = w*wph
         wlepang = wcos*wph
         sphi = sin(phi)
         cphi = cos(phi)
         p1mod = sqla/2.d0/mass
         p3(0) = sqrt(p1mod**2 + m1**2)
         p3(1) = p1mod * sphi * s
         p3(2) = p1mod * cphi * s
         p3(3) = p1mod * c
         p4(0) = mass - p3(0)
         p4(1) = - p3(1)
         p4(2) = - p3(2)
         p4(3) = - p3(3)
         phsp = phsp * duepigrechi
         return
      endif
*******************
******************
* the following for ng >= 1
      pflatc = 0.d0 ! for collinear sampling, flat fraction
      call get_cos_fer(rnd(ir),nphot,cth,wcos)
      w = w * wcos
      sth = sqrt(1.d0 - cth**2)
      call get_phi(rnd(ir+1),phi,wph)
      w = w*wph
      wlepang = wcos*wph
      sphi = sin(phi)
      cphi = cos(phi)
      ir = ir + 2
      do k = 1,4
         nph(k) = 0
      enddo
      wpattern = 1.d0
      call get_pattern(ncharged,nphot,emissionpattern)
c      call get_pattern_v2(ncharged,nphot,emissionpattern,wpattern)
      w = w*wpattern
      do k = 1,nphot
         nph(emissionpattern(k)) = nph(emissionpattern(k)) + 1
      enddo
      n1 = nph(1)
      n2 = nph(2)
      n3 = nph(3)
      n4 = nph(4)

      ei   = p1(0)+p2(0)
      beta = (1.d0+m1/p1(0))*(1.d0-m1/p1(0))

*** it can appens, for example for muon-neutrino-gamma final state
*** that ei > mmu + mnu + esoft, but p1(0)<mmu...!!!
      if (beta.le.0.d0) beta = 0.1d0


c      betaold = sqrt(1.d0 - m1**2/p1(0)**2)
      beta = sqrt(beta)

      dir3(0) = p1(0)
      dir3(1) = p1(0)*beta*cphi*sth
      dir3(2) = p1(0)*beta*sphi*sth
      dir3(3) = p1(0)*beta*cth         
      dir4(0) = dir3(0)
      pp(0)   = p1(0)+p2(0)
      do k = 1,3
         pp(k)    = p1(k)+p2(k)
         dir4(k) = -dir3(k)
      enddo               
** PHOTONS GENERATION...
      nis = n1 + n2
      nfs = n3 + n4
** initial state photons - 1
      isharebw = 0
      if (n1.ge.1) isharebw = isharebw + 1 
      if (n2.ge.1) isharebw = isharebw + 1 
      if (n3.ge.1) isharebw = isharebw + 1 
      if (n4.ge.1) isharebw = isharebw + 1 

      if (boson.eq.'W') then
         amv = mw
         agv = gw
      else
         amv = mz
         agv = gz
      endif
      
      if (nis.gt.0) then
         if (n1.gt.0) then
            do k = 1,max(n1,1)
               pflat1 = pflatc
c               pflat1 = 0.1d0
               c1 = 2.d0/3.d0
               c2 = 1.d0/3.d0
               call wraprng(rnd,2)
               call collinearm2(c1,c2,muno,pflat1,rnd(1),rnd(2),
     .              p1,vers,wverg)
               cgp3 = tridot(vers,dir3)/beta/dir3(0)
               ir = ir + 2
               w = w * wverg
               wphotang = wphotang * wverg

               pbwbase = 0.2d0/isharebw
               if (boson.eq.'W') then
                  if (pp(0).gt.90.d0)  pbwbase = 0.7d0/isharebw
                  if (pp(0).gt.120.d0) pbwbase = 0.8d0/isharebw
               else
                  if (pp(0).gt.101.d0)  pbwbase = 0.7d0/isharebw
                  if (pp(0).gt.131.d0) pbwbase = 0.8d0/isharebw
               endif
               pbw = pbwbase

               q2min = (m1+m2)**2
               call wraprng(csi,1)
               call photon_energy(pbw,csi(1),pp,vers,q2min,esoft,
     .              p1(0),amv,agv,om,wom)
               oph(k) = om
               ir = ir + 1
               w = w * wom
               wphoten = wphoten * wom

               qph(k,0) = oph(k)
               qph(k,1) = oph(k)*vers(1)
               qph(k,2) = oph(k)*vers(2)
               qph(k,3) = oph(k)*vers(3)

               pp(0)    = pp(0) - qph(k,0)
               pp(1)    = pp(1) - qph(k,1)
               pp(2)    = pp(2) - qph(k,2)
               pp(3)    = pp(3) - qph(k,3)
            enddo
         endif
         if (n2.gt.0) then
            do k = 1,max(n2,1)
               pflat2 = pflatc
c               pflat2 = 0.1d0
               c1 = 1.d0/3.d0
               c2 = 2.d0/3.d0
               call wraprng(rnd,2)
               call collinearm2(c1,c2,mdue,pflat2,rnd(1),rnd(2),p2,
     .              vers,wverg)
               cgp3 = tridot(vers,dir3)/beta/dir3(0)

               ir = ir + 2
               w = w * wverg
               wphotang = wphotang * wverg

               pbwbase = 0.2d0/isharebw
               if (boson.eq.'W') then
                  if (pp(0).gt.90.d0)  pbwbase = 0.7d0/isharebw
                  if (pp(0).gt.120.d0) pbwbase = 0.8d0/isharebw
               else
                  if (pp(0).gt.101.d0)  pbwbase = 0.7d0/isharebw
                  if (pp(0).gt.131.d0) pbwbase = 0.8d0/isharebw
               endif
               pbw = pbwbase

               q2min = (m1+m2)**2
               call wraprng(csi,1)
               call photon_energy(pbw,csi(1),pp,vers,q2min,esoft,
     .              p1(0),amv,agv,om,wom)
               oph(k+n1) = om
               ir = ir + 1
               w = w * wom
               wphoten = wphoten * wom

               qph(k+n1,0) = oph(k+n1)
               qph(k+n1,1) = oph(k+n1)*vers(1)
               qph(k+n1,2) = oph(k+n1)*vers(2)
               qph(k+n1,3) = oph(k+n1)*vers(3)

               pp(0)    = pp(0) - qph(k+n1,0)
               pp(1)    = pp(1) - qph(k+n1,1)
               pp(2)    = pp(2) - qph(k+n1,2)
               pp(3)    = pp(3) - qph(k+n1,3)
            enddo
         endif
      endif

      if (dot(pp,pp).lt.(m1+m2)**2.or.pp(0).lt.0.d0) then
         ie = 1
         phsp = 0.d0
         w  = 0.d0
         return
      endif
** final state photons

      iref = 3
      if (boson.eq.'Z') then
         call wraprng(csi,1)
         if (csi(1).gt.0.5) iref = 4
         wref = 1.d0
         if (nfs.gt.0) then
            if (n4.eq.0) iref = 3
            if (n3.eq.0) iref = 4
         endif
         w = w*wref
      endif

***********
      if (nfs.gt.0) then
         if (iref.eq.3) then
            if (n3.gt.0) then
               do k = 1,max(n3,1)
                  c1 = 1.d0
                  c2 = 1.d0
                  call wraprng(rnd,2)

                  call collinearm2(c1,c2,m1,pflatc,rnd(1),rnd(2),dir3,
     .                 vers,wverg)
                     
c     pflatt = 0.05d0
c     call collinear(pflatt,rnd(1),rnd(2),dir3,vers,wverg)
                  cgp3 = tridot(vers,dir3)/beta/dir3(0)
                  
                  ir = ir + 2
                  w = w * wverg
                  wphotang = wphotang * wverg
                  
                  pbwbase = 0.2d0/isharebw
                  if (boson.eq.'W') then
                     if (pp(0).gt.90.d0)  pbwbase = 0.7d0/isharebw
                     if (pp(0).gt.120.d0) pbwbase = 0.8d0/isharebw
                  else
                     if (pp(0).gt.101.d0)  pbwbase = 0.7d0/isharebw
                     if (pp(0).gt.131.d0) pbwbase = 0.8d0/isharebw
                  endif
                  pbw = pbwbase
                  
                  q2min = (m1+m2)**2
                  call wraprng(csi,1)
                  call photon_energy(pbw,csi(1),pp,vers,q2min,esoft,
     .                 p1(0),amv,agv,om,wom)
                  oph(k+n1+n2) = om
                  ir = ir + 1
                  w = w * wom
                  wphoten = wphoten * wom
                  
                  qph(k+n1+n2,0) = oph(k+n1+n2)
                  qph(k+n1+n2,1) = oph(k+n1+n2)*vers(1)
                  qph(k+n1+n2,2) = oph(k+n1+n2)*vers(2)
                  qph(k+n1+n2,3) = oph(k+n1+n2)*vers(3)
                  
                  pp(0)    = pp(0) - qph(k+n1+n2,0)
                  pp(1)    = pp(1) - qph(k+n1+n2,1)
                  pp(2)    = pp(2) - qph(k+n1+n2,2)
                  pp(3)    = pp(3) - qph(k+n1+n2,3)
               enddo
            endif
*********
            if (n4.gt.0) then
               if (dot(pp,pp).lt.(m1+m2)**2.or.pp(0).lt.0.d0) then
                  ie = 1
                  phsp = 0.d0
                  w  = 0.d0
                  return
               endif
               do k = 0,3
                  pcoll(k) = pp(k) - dir3(k)
               enddo
               pcmodmu = 1.d0/sqrt(tridot(pcoll,pcoll))            
               pcoll(0) = dir4(0)
               bcoll    = sqrt(1.d0 - m1**2/pcoll(0)**2)
               pcoll(1) = pcoll(1)*pcmodmu*bcoll*pcoll(0)
               pcoll(2) = pcoll(2)*pcmodmu*bcoll*pcoll(0)
               pcoll(3) = pcoll(3)*pcmodmu*bcoll*pcoll(0)

               do k = 1,n4
                  c1 = 1.d0
                  c2 = 1.d0
                  call wraprng(rnd,2)
                  call collinearm2(c1,c2,m2,pflatc,rnd(1),rnd(2),pcoll,
     .                 vers,wverg)
                  cgp3 = tridot(vers,dir3)/beta/dir3(0)
                  
                  ir = ir + 2
                  w = w * wverg
                  wphotang = wphotang * wverg
                  
                  pbwbase = 0.2d0/isharebw
                  if (boson.eq.'W') then
                     if (pp(0).gt.90.d0)  pbwbase = 0.7d0/isharebw
                     if (pp(0).gt.120.d0) pbwbase = 0.8d0/isharebw
                  else
                     if (pp(0).gt.101.d0)  pbwbase = 0.7d0/isharebw
                     if (pp(0).gt.131.d0) pbwbase = 0.8d0/isharebw
                  endif
                  pbw = pbwbase
                  
                  q2min = (m1+m2)**2
                  call wraprng(csi,1)
                  call photon_energy(pbw,csi(1),pp,vers,q2min,esoft,
     .                 p1(0),amv,agv,om,wom)
                  oph(k+n1+n2+n3) = om
                  ir = ir + 1
                  w = w * wom
                  wphoten = wphoten * wom
                  
                  qph(k+n1+n2+n3,0) = oph(k+n1+n2+n3)
                  qph(k+n1+n2+n3,1) = oph(k+n1+n2+n3)*vers(1)
                  qph(k+n1+n2+n3,2) = oph(k+n1+n2+n3)*vers(2)
                  qph(k+n1+n2+n3,3) = oph(k+n1+n2+n3)*vers(3)
                  
                  pp(0)    = pp(0) - qph(k+n1+n2+n3,0)
                  pp(1)    = pp(1) - qph(k+n1+n2+n3,1)
                  pp(2)    = pp(2) - qph(k+n1+n2+n3,2)
                  pp(3)    = pp(3) - qph(k+n1+n2+n3,3)
               enddo
            endif
         else ! if (iref.eq.4)
            if (n4.gt.0) then
               do k = 1,max(n4,1)
                  c1 = 1.d0
                  c2 = 1.d0
                  call wraprng(rnd,2)
                  call collinearm2(c1,c2,m2,pflatc,rnd(1),rnd(2),dir4,
     .                 vers,wverg)
                  cgp3 = tridot(vers,dir3)/beta/dir3(0)
                  
                  ir = ir + 2
                  w = w * wverg
                  wphotang = wphotang * wverg
                  
                  pbwbase = 0.2d0/isharebw
                  if (boson.eq.'W') then
                     if (pp(0).gt.90.d0)  pbwbase = 0.7d0/isharebw
                     if (pp(0).gt.120.d0) pbwbase = 0.8d0/isharebw
                  else
                     if (pp(0).gt.101.d0)  pbwbase = 0.7d0/isharebw
                     if (pp(0).gt.131.d0) pbwbase = 0.8d0/isharebw
                  endif
                  pbw = pbwbase
                  
                  q2min = (m1+m2)**2
                  call wraprng(csi,1)
                  call photon_energy(pbw,csi(1),pp,vers,q2min,esoft,
     .                 p1(0),amv,agv,om,wom)
                  oph(k+n1+n2+n3) = om
                  ir = ir + 1
                  w = w * wom
                  wphoten = wphoten * wom
                  
                  qph(k+n1+n2+n3,0) = oph(k+n1+n2+n3)
                  qph(k+n1+n2+n3,1) = oph(k+n1+n2+n3)*vers(1)
                  qph(k+n1+n2+n3,2) = oph(k+n1+n2+n3)*vers(2)
                  qph(k+n1+n2+n3,3) = oph(k+n1+n2+n3)*vers(3)
                  
                  pp(0)    = pp(0) - qph(k+n1+n2+n3,0)
                  pp(1)    = pp(1) - qph(k+n1+n2+n3,1)
                  pp(2)    = pp(2) - qph(k+n1+n2+n3,2)
                  pp(3)    = pp(3) - qph(k+n1+n2+n3,3)
               enddo
            endif

            if (n3.gt.0) then
               if (dot(pp,pp).lt.(m1+m2)**2.or.pp(0).lt.0.d0) then
                  ie = 1
                  phsp = 0.d0
                  w  = 0.d0
                  return
               endif
               do k = 0,3
                  pcoll(k) = pp(k) - dir4(k)
               enddo

               pcmodmu = 1.d0/sqrt(tridot(pcoll,pcoll))            
               pcoll(0) = dir3(0)
               bcoll    = sqrt(1.d0 - m2**2/pcoll(0)**2)
               pcoll(1) = pcoll(1)*pcmodmu*bcoll*pcoll(0)
               pcoll(2) = pcoll(2)*pcmodmu*bcoll*pcoll(0)
               pcoll(3) = pcoll(3)*pcmodmu*bcoll*pcoll(0)

               do k = 1,n3
                  c1 = 1.d0
                  c2 = 1.d0
                  call wraprng(rnd,2)
                  call collinearm2(c1,c2,m1,pflatc,rnd(1),rnd(2),pcoll,
     .                 vers,wverg)
                  cgp3 = tridot(vers,dir3)/beta/dir3(0)

                  ir = ir + 2
                  w = w * wverg
                  wphotang = wphotang * wverg
                  
                  pbwbase = 0.2d0/isharebw
                  if (boson.eq.'W') then
                     if (pp(0).gt.90.d0)  pbwbase = 0.7d0/isharebw
                     if (pp(0).gt.120.d0) pbwbase = 0.8d0/isharebw
                  else
                     if (pp(0).gt.101.d0)  pbwbase = 0.7d0/isharebw
                     if (pp(0).gt.131.d0) pbwbase = 0.8d0/isharebw
                  endif
                  pbw = pbwbase
                  
                  q2min = (m1+m2)**2
                  call wraprng(csi,1)
                  call photon_energy(pbw,csi(1),pp,vers,q2min,esoft,
     .                 p1(0),amv,agv,om,wom)
                  oph(k+n1+n2) = om
                  ir = ir + 1
                  w = w * wom
                  wphoten = wphoten * wom
                  
                  qph(k+n1+n2,0) = oph(k+n1+n2)
                  qph(k+n1+n2,1) = oph(k+n1+n2)*vers(1)
                  qph(k+n1+n2,2) = oph(k+n1+n2)*vers(2)
                  qph(k+n1+n2,3) = oph(k+n1+n2)*vers(3)
                  
                  pp(0)    = pp(0) - qph(k+n1+n2,0)
                  pp(1)    = pp(1) - qph(k+n1+n2,1)
                  pp(2)    = pp(2) - qph(k+n1+n2,2)
                  pp(3)    = pp(3) - qph(k+n1+n2,3)
               enddo
            endif
         endif
      endif

***********************************
      do k = 0,3
         kone(k) = 0.d0
      enddo
      do j = 0,3
         do k = 1,nphot
            kone(j) = kone(j) + qph(k,j)
         enddo
         ptmp(j) = p1(j) + p2(j) - kone(j)
      enddo
      if (dot(ptmp,ptmp).lt.(m1+m2)**2.or.ptmp(0).lt.0.d0) then
         ie = 1
         phsp = 0.d0
         w  = 0.d0
         return
      endif
      amkone2 = dot(kone,kone)
***
c old      call wraprng(csi,1)
c old      ich = nphot * csi(1) + 1
c old      ich = min(ich,nphot) ! to avoid round-offs
c old      iref = 3
c old      if (ich.gt.(n1+n2+n3)) iref = 4
***         
      if (iref.eq.3) then
         vx  = sth*cphi
         vy  = sth*sphi
         vz  = cth
      else
         vx  = -sth*cphi
         vy  = -sth*sphi
         vz  = -cth
      endif

      o   = kone(0)
      ox  = kone(1)
      oy  = kone(2)
      oz  = kone(3)
      odv = ox*vx + oy*vy + oz*vz

      abig = ei**2 - 2.d0 * ei * o + amkone2 + m1**2 - m2**2
      bbig = -2.d0 * (ei - o)
      cbig = -2.d0 * odv
      arg  = abig**2*bbig**2-m1**2*bbig**2*(bbig**2-cbig**2)
      if (arg.lt.0.d0) then
         phsp = 0.d0
         w    = 0.d0
         ie = 2
         return
      endif
      p1mod_1 = abig*cbig + sqrt(arg)
      p1mod_1 = p1mod_1 / (bbig**2 - cbig**2)
      p1mod_2 = abig*cbig - sqrt(arg)
      p1mod_2 = p1mod_2 / (bbig**2 - cbig**2)
      p1jac = 1.d0
      if (p1mod_1.gt.0.d0.and.p1mod_2.gt.0.d0) then
         call wraprng(csi,1)
         if (csi(1).lt.0.5) p1mod = p1mod_1
         if (csi(1).ge.0.5) p1mod = p1mod_2
         p1jac = 2.d0
      endif
      if (p1mod_1.gt.0.d0.and.p1mod_2.lt.0.d0) then
         p1mod = p1mod_1
      endif
      if (p1mod_1.lt.0.d0.and.p1mod_2.gt.0.d0) then
         p1mod = p1mod_2
      endif
      if (p1mod_1.lt.0.d0.and.p1mod_2.lt.0.d0) then
         phsp = 0.d0
         w    = 0.d0
         ie   = 3
         return
      endif

      w       = w * p1jac
      wphoten = wphoten * p1jac

      p3(0) = sqrt(p1mod**2 + m1**2)
      p3(1) = p1mod*vx
      p3(2) = p1mod*vy
      p3(3) = p1mod*vz      
      p4(0) =  ei - p3(0) - o
      p4(1) = -ox - p3(1)
      p4(2) = -oy - p3(2)
      p4(3) = -oz - p3(3)

!      phasespace!
      prodomega = 1.d0
      do k = 1,nphot
         prodomega = prodomega * qph(k,0)
      enddo
      den = dabs( p4(0) + p3(0)*(1.d0 + odv/p1mod) )
      phsp = prodomega * p1mod / den
      phsp = phsp*duepigrechi/2.d0**(npart)
      if (iref.eq.4) call exchange_mom(p3,p4)
***********************      
! regolatore
!!! only on the final state -- as in BabaYaga... 
      den = 0.d0
      if (n3.gt.0.and.n4.gt.0) then ! this implies only Z!!!
         prodp3 = 1.d0
         prodp4 = 1.d0
         do k = nis+1,nphot
            do i = 0,3
               ptmp(i) = qph(k,i)
            enddo
            p3q = dot(p3,ptmp)
            p4q = dot(p4,ptmp)

ccc trying to avoid NaNs at extremely high energies            
            if (p3q.eq.0.d0) then
               xx  = -m1*m1/p3(0)/p3(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               p3q = p3(0)*ptmp(0)*ombet
            endif
            if (p4q.eq.0.d0) then
               xx  = -m2*m2/p4(0)/p4(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               p4q = p4(0)*ptmp(0)*ombet
            endif
cccc

            if (k.le.(n1+n2+n3)) then
               prodp3=prodp3*
     .              (1.d0/p3q*p3(0)*ptmp(0))
            endif
            if (k.gt.(n1+n2+n3)) then
               prodp4=prodp4*
     .              (1.d0/p4q*p4(0)*ptmp(0))
            endif
         enddo
         if (iref.eq.3) anum = prodp3
         if (iref.eq.4) anum = prodp4
         regulator1 = anum/(prodp3+prodp4) * 2.d0
         w = w * regulator1
         reg = reg * regulator1
      endif
*******************
* secondo regolatore
      ifour = 1
      if (boson.eq.'W') ifour = 0
      clog1 = 1.d0
      clog2 = 1.d0
      clog3 = 1.d0
      clog4 = 1.d0
      den  = 1.d0
      anum = 1.d0
      if (n1.gt.0) then
         do k = 1,n1
            do i = 0,3
               ptmp(i) = qph(k,i)
            enddo            

            p1q = dot(p1,ptmp)
            p2q = dot(p2,ptmp)
            p3q = dot(p3,ptmp)
            p4q = dot(p4,ptmp)

ccc trying to avoid NaNs at extremely high energies            
            if (p3q.eq.0.d0) then
               xx  = -m1*m1/p3(0)/p3(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               p3q = p3(0)*ptmp(0)*ombet
            endif
            if (p4q.eq.0.d0) then
               xx  = -m2*m2/p4(0)/p4(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               p4q = p4(0)*ptmp(0)*ombet
            endif
cccc
            
            den = den * (
     .           1.d0/p1q*p1(0)*ptmp(0)/clog1 +
     .           1.d0/p2q*p2(0)*ptmp(0)/clog2 +
     .           1.d0/p3q*p3(0)*ptmp(0)/clog3 +
     .     ifour*1.d0/p4q*p4(0)*ptmp(0)/clog4 )
            anum = anum/p1q*p1(0)*ptmp(0)/clog1
         enddo
      endif
      if (n2.gt.0) then
         do k = n1+1,n1+n2
            do i = 0,3
               ptmp(i) = qph(k,i)
            enddo

            p1q = dot(p1,ptmp)
            p2q = dot(p2,ptmp)
            p3q = dot(p3,ptmp)
            p4q = dot(p4,ptmp)

ccc trying to avoid NaNs at extremely high energies            
            if (p3q.eq.0.d0) then
               xx  = -m1*m1/p3(0)/p3(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               p3q = p3(0)*ptmp(0)*ombet
            endif
            if (p4q.eq.0.d0) then
               xx  = -m2*m2/p4(0)/p4(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               p4q = p4(0)*ptmp(0)*ombet
            endif
cccc
            
            den = den * (
     .           1.d0/p1q*p1(0)*ptmp(0)/clog1 +
     .           1.d0/p2q*p2(0)*ptmp(0)/clog2 +
     .           1.d0/p3q*p3(0)*ptmp(0)/clog3 +
     .     ifour*1.d0/p4q*p4(0)*ptmp(0)/clog4 )
            anum = anum/p2q*p2(0)*ptmp(0)/clog2
         enddo
      endif
      if (n3.gt.0) then
         do k = n1+n2+1,n1+n2+n3
            do i = 0,3
               ptmp(i) = qph(k,i)
            enddo

            p1q = dot(p1,ptmp)
            p2q = dot(p2,ptmp)
            p3q = dot(p3,ptmp)
            p4q = dot(p4,ptmp)

ccc   trying to avoid NaNs at extremely high energies            
            if (p3q.eq.0.d0) then
               xx  = -m1*m1/p3(0)/p3(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               p3q = p3(0)*ptmp(0)*ombet
            endif
            if (p4q.eq.0.d0) then
               xx  = -m2*m2/p4(0)/p4(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               p4q = p4(0)*ptmp(0)*ombet
            endif
cccc
                        
            den = den * (
     .           1.d0/p1q*p1(0)*ptmp(0)/clog1 +
     .           1.d0/p2q*p2(0)*ptmp(0)/clog2 +
     .           1.d0/p3q*p3(0)*ptmp(0)/clog3 +
     .     ifour*1.d0/p4q*p4(0)*ptmp(0)/clog4 )
            anum = anum/p3q*p3(0)*ptmp(0)/clog3
         enddo
      endif
      if (n4.gt.0) then
         do k = n1+n2+n3+1,nphot
            do i = 0,3
               ptmp(i) = qph(k,i)
            enddo

            p1q = dot(p1,ptmp)
            p2q = dot(p2,ptmp)
            p3q = dot(p3,ptmp)
            p4q = dot(p4,ptmp)

ccc   trying to avoid NaNs at extremely high energies            
            if (p3q.eq.0.d0) then
               xx  = -m1*m1/p3(0)/p3(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               p3q = p3(0)*ptmp(0)*ombet
            endif
            if (p4q.eq.0.d0) then
               xx  = -m2*m2/p4(0)/p4(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               p4q = p4(0)*ptmp(0)*ombet
            endif
cccc

            den = den * (
     .           1.d0/p1q*p1(0)*ptmp(0)/clog1 +
     .           1.d0/p2q*p2(0)*ptmp(0)/clog2 +
     .           1.d0/p3q*p3(0)*ptmp(0)/clog3 +
     .     ifour*1.d0/p4q*p4(0)*ptmp(0)/clog4 )
            anum = anum/p4q*p4(0)*ptmp(0)/clog4
         enddo
      endif
      regulator2 = 1.d0* ncharged**(nphot) * anum/den
!!!   ATTENZIONE CON get_pattern usare la linea precedente!!!
! con get_pattern_v2 quella successiva!!
cc      regulator2 = 1.d0 * anum/den
      w   = w * regulator2
      reg = reg * regulator2
***********************
      return
      end
*
*** MATRIX ELEMENTS
      subroutine real_photons(model,ng,ecms,p1,p2,qph,ie,imtx,summt2,
     .     iemap,emtxsub,iphind)
      implicit double precision (a-h,m,o-z)
      character*10 model
      character*6 eikng ! multip single
      parameter (imaxph = 40)
      dimension iclose(imaxph)
      real*4 csi(1)
      character*1 boson
      double precision lcoll,mass,lambda
      double precision mq(-5:5),chq(-5:5),itre(-5:5)
      dimension xap(imaxph),apx(imaxph)
      dimension p1(0:3),p2(0:3),pin1(0:3),pin2(0:3),q(0:3)
      dimension p1r(0:3),p2r(0:3),p1b(0:3),p2b(0:3)
      dimension q1(0:3),q2(0:3),q3(0:3),q4(0:3),qq(0:3),qr(0:3)
      dimension p1tmp(0:3),p2tmp(0:3),p3tmp(0:3),p4tmp(0:3)
      dimension ptmp(0:3),di(imaxph),fc(imaxph)
      dimension pin1bl(0:3),pin2bl(0:3),p1bl(0:3),p2bl(0:3),vboost(0:3)
      dimension qph(imaxph,0:3),qphtmp(imaxph,0:3),randvec(3)
      dimension pin1sub(0:3),pin2sub(0:3),p1sub(0:3),p2sub(0:3)
      common/momentainitial/pin1,pin2
      common/reducedtoborn/p1b,p2b,iref
      common/momentainitialred/pin1b(0:3),pin2b(0:3),pin1r(0:3)
     >     ,pin2r(0:3)
      common/forborncrosssection/phsp2b,flux2b,bornme,bornmeclean
      common/radpattern/nph(4)
      common/chargesmasses_cm/charges(4),amasses(4)
      common/iclosest/iclose
      common/closesubtraction/iclosesub(imaxph)
      common/scalapersottrazione/sqrtscalesub
      common/kforsub/ksub
      common/singleseikonals/singlephotons(40)
      character*6 ord
      common/qedORDER/ord
!from shared.inc
      common/vectorboson/boson
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt
     .     ,mq,chq,itre
      common/processhor/m1,m2,ch1,ch2,chfs,mfs,lepton,imirror
      common/iplsubtraction/iplsub
      common/gen_event/ebeam,ebeam1,ebeam2,xmin,xmax,shat,x1pdf,x2pdf
      common/partons/ipart1,ipart2
      common/ipartsmirrrored/ipart1mirr,ipart2mirr
! imtx is to say that the matrix el. has been calculated..
! summt2 is the squared ME summed over all spins
! NB the 16.d0 in front of the oal ME are necessary...
      common/sampling_weight/wsampling,spaziofasi,flusso
      common/counter_realphotons/icountrf,ii1,ii2,ii3
      data icountrf,ii1,ii2,ii3 /0,0,0,0/
      icountrf = icountrf + 1
      emtxsub = 0.d0
      iemap = 0
      Qu = charges(1)
      Qd = charges(2)
      Qe = charges(3)
      Qw = -Qe

      Qqu = charges(1)
      Qel = charges(3)

      imtx  = 0
      summt2 = 0.d0
      if (ie.gt.0) return
      imtx = 1

      if (iphind.eq.1) then
         i_i_i = 2
c         if (ipart1mirr.eq.2.or.ipart1mirr.eq.4) i_i_i = 1
         if (abs(ipart1mirr).eq.2.or.abs(ipart1mirr).eq.4) i_i_i = 1
         do i=0,3
            q(i) = qph(1,i)
         enddo            
         if (boson.eq.'W') then
            if (i_i_i.eq.1) then
               summt2=16.d0*
     .              elmat2phindw(i_i_i,pin1,q,p2,p1,pin2,Qu,Qd,Qw,Qe)
            else
               summt2=16.d0*
     .              elmat2phindw(i_i_i,q,pin1,p2,p1,pin2,Qu,Qd,Qw,Qe)
            endif
         else
            ! ipart2mirr = 0 sempre
            if (ipart1mirr.gt.0) then
               i_i_i = 1
               summt2 = elmat2phindz(i_i_i,pin1,q,p1,p2,pin2,Qqu,Qel)
            else
               i_i_i = 2
               summt2 = elmat2phindz(i_i_i,q,pin1,p1,p2,pin2,Qqu,Qel)
            endif
         endif
         return
      endif

      if (ng.eq.0) then
         do k = 0,3
            pin1b(k) = pin1(k)
            pin2b(k) = pin2(k)
            p1b(k) = p1(k)
            p2b(k) = p2(k)
         enddo
      else
         call closest('full',ng,qph,pin1,pin2,p1,p2,q1,q2,q3,q4,iclose)
         do k = 0,3
            p1tmp(k) = pin1(k)  - q1(k)
            p2tmp(k) = pin2(k)  - q2(k)
            qq(k)    = p1tmp(k) + p2tmp(k)
         enddo
         call maptozero(p1tmp,p2tmp,p1,p2,amasses,pin1b,pin2b,p1b,p2b,
     $        iemap)
!         call maptozerov2(ng,pin1,pin2,p1,p2,qph,iclose,amasses,
!     .        pin1b,pin2b,p1b,p2b,iemap)
*[ safest choice
         if (iemap.eq.1.and.ng.ge.1.and.ord.ne.'alpha') then
ccc         if (iemap.eq.1.and.ng.gt.1) then
*]
c$$$            print*,'ng ',ng
c$$$            print*,'pin1 ',pin1
c$$$            print*,'pin2 ',pin2
c$$$            print*,'p1 ',p1
c$$$            print*,'p2 ',p2
c$$$            print*,'pin1b ',pin1b
c$$$            print*,'pin2b ',pin2b
c$$$            print*,'p1b ',p1b
c$$$            print*,'p2b ',p2b
            return
         endif
         iemap = 0
      endif
****
! to be passed in a common!
      flux2b = 8.d0*pin1b(0)*pin2b(0)
      mass = pin1b(0)+pin2b(0)
      sqla = sqrt(lambda(mass**2,0.d0,0.d0))
      phsp2b = sqla/8.d0/mass**2 /(2.d0*pi) / 4.d0 ! spins...
****
** doing as if only the k-th photon is present
      if (model.eq.'matchedps') then
         if (ng.eq.0) then
            iref = 1
            if (boson.eq.'W') then 
               bornme = elmat2born(pin1b,pin2b,p2b,p1b)
            else
               if (ipart1.eq.0.and.ipart2.eq.0) then
                  bornme = elmat2borngg(pin1b,pin2b,p1b,p2b)
                  bornmeclean = bornme
               else
                  bornme = elmat2bornz(pin1b,pin2b,p1b,p2b,1)
                  bornmeclean = elmat2bornz(pin1b,pin2b,p1b,p2b,0)
               endif
            endif
            summt2 = bornme
            imtx = 1
            return
         endif
         if (ng.eq.1) then
            do i=0,3
               q(i) = qph(1,i)
            enddo            
            if (boson.eq.'W') then
               summt2 = 16.d0*elmat2form(pin1,pin2,p2,p1,q,Qu,Qd,Qw,Qe)
               bornme = elmat2born(pin1b,pin2b,p2b,p1b)
            else
               summt2 = elmat2formz(pin1,pin2,p1,p2,q,Qqu,Qel)
               bornme = elmat2bornz(pin1b,pin2b,p1b,p2b,1)               
               bornmeclean = elmat2bornz(pin1b,pin2b,p1b,p2b,0)

               if (ord.eq.'exp') summt2 = summt2 * bornme/bornmeclean

            endif
            imtx = 1
            
*** also for ng=1 I need to calculate the subtraction term...
            ei1r = pin1(0)
            ei2r = pin2(0)
            ef3r = p1(0)
            ef4r = p2(0)
            if (iclose(1).eq.1) then
               ef1r  = ei1r - q(0)
               xjaci = 1.d0/ei1r/2.d0
               xapi  = q(0)/ei1r
               apxi  = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi  = apxi/(1.d0 - xapi)
            endif
            if (iclose(1).eq.2) then
               ef2r  = ei2r - q(0)
               xjaci = 1.d0/ei2r/2.d0
               xapi  = q(0)/ei2r
               apxi  = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi  = apxi/(1.d0 - xapi)
            endif
            if (iclose(1).eq.3) then 
               ei3r = ef3r + q(0)
               xjaci = 1.d0/ei3r/2.d0
               xapi = q(0)/ei3r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi = apxi/(1.d0 - xapi) ! attenzione !!
            endif
            if (iclose(1).eq.4) then 
               ei4r = ef4r + q(0)
               xjaci = 1.d0/ei4r/2.d0
               xapi = q(0)/ei4r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi = apxi/(1.d0 - xapi) ! attenzione !!
            endif

*[[[[[
** PS approximation
c            ksub = -1000
c            call eikonal_factor('single',ng,pin1,pin2,p1,p2,q,
c     >           qph,eiki)
c            summt2approx = bornme*apxi*eiki*q(0)*xjaci
c            if (q(0).lt.0.5d0) then 
c               print*,'----- matching.f line 909'
cc               print*,summt2approx,summt2
c               print*,summt2approx/summt2,summt2/summt2approx
c            endif
*]]]]]

            sub     = 0.d0
            emtxsub = 0.d0
            riscala = 1.d0
            if (boson.eq.'Z') riscala = bornmeclean/bornme
            k      = 1
            scale  = sqrtscalesub
            angsub = subtractionnew(k,scale,amasses(1),amasses(2),
     .           charges(1),charges(2),q)
            sub = bornme*apxi*q(0)*xjaci*angsub*riscala
            sub = sub*iplsub

            summt2  = summt2 - sub
            emtxsub = sub

c            if (isnan(summt2)) then
c               print*,'matching line 1140',sub
c            endif

            return
         endif

         if (boson.eq.'W') then
            bornme = elmat2born(pin1b,pin2b,p2b,p1b)
         else
            bornme = elmat2bornz(pin1b,pin2b,p1b,p2b,1)
            bornmeclean = elmat2bornz(pin1b,pin2b,p1b,p2b,0)
         endif

         eikng = 'multip'
         ksub = -1000
         call eikonal_factor(eikng,ng,pin1,pin2,p1,p2,q,qph,eikfull)

         ei1 = pin1(0)
         ei2 = pin2(0)
         ef1 = ei1
         ef2 = ei2
         ef3 = p1(0)
         ef4 = p2(0)
         prodomega = 1.d0
         do k = 1,ng
            prodomega = prodomega*qph(k,0)
            if (iclose(k).eq.3) ef3 = ef3 + qph(k,0)
            if (iclose(k).eq.4) ef4 = ef4 + qph(k,0)
         enddo
         ei3 = ef3
         ei4 = ef4

         prodapx = 1.d0
         xjac    = 1.d0
         delta   = 0.d0
         prod    = 1.d0
         delta_eik   = 0.d0
         prod_eik    = 1.d0
         do k = 1,ng
            do i=0,3
               q(i) = qph(k,i)
            enddo            
**** PS QUANTITIES
            if (iclose(k).eq.1) then
               xjac = xjac/ef1/2.d0
               xap(k) = qph(k,0)/ef1
               apx(k) = (1.d0 + (1.d0 - xap(k))**2)/xap(k)
               apx(k) = apx(k)/(1.d0 - xap(k)) ! attenzione !!
               ef1    = ef1 - qph(k,0)
            endif
            if (iclose(k).eq.2) then
               xjac = xjac/ef2/2.d0
               xap(k) = qph(k,0)/ef2
               apx(k) = (1.d0 + (1.d0 - xap(k))**2)/xap(k)
               apx(k) = apx(k)/(1.d0 - xap(k)) ! attenzione !!
               ef2    = ef2 - qph(k,0)
            endif
            if (iclose(k).eq.3) then
               xjac = xjac/ef3/2.d0
               xap(k) = qph(k,0)/ef3
               apx(k) = (1.d0 + (1.d0 - xap(k))**2)/xap(k)
               apx(k) = apx(k)/(1.d0 - xap(k)) ! attenzione !!
               ef3    = ef3 - qph(k,0)
            endif
            if (iclose(k).eq.4) then
               xjac = xjac/ef4/2.d0
               xap(k) = qph(k,0)/ef4
               apx(k) = (1.d0 + (1.d0 - xap(k))**2)/xap(k)
               apx(k) = apx(k)/(1.d0 - xap(k)) ! attenzione !!
               ef4    = ef4 - qph(k,0)
            endif
            prodapx = prodapx * apx(k)
*****
** doing as if only the k-th photon is present !!!
            i1 = 0
            i2 = 0
            i3 = 0
            i4 = 0
            if (iclose(k).eq.1) i1 = 1
            if (iclose(k).eq.2) i2 = 1
            if (iclose(k).eq.3) i3 = 1
            if (iclose(k).eq.4) i4 = 1
            do ki = 0,3
               p1tmp(ki) = pin1(ki) - i1 * q(ki)
               p2tmp(ki) = pin2(ki) - i2 * q(ki)
               qq(ki)    = p1tmp(ki) + p2tmp(ki)
            enddo
            call maptozero(p1tmp,p2tmp,p1,p2,amasses,pin1bl,pin2bl,
     .           p1bl,p2bl,iemap)
            iemap = 0
            if (boson.eq.'W') then
               bornmel = elmat2born(pin1bl,pin2bl,p2bl,p1bl)
            else
               bornmel = elmat2bornz(pin1bl,pin2bl,p1bl,p2bl,0)
            endif
*******
            is = 0
            if (iclose(k).eq.1.or.iclose(k).eq.2) is = 1
            call maptoone(is,qq,p1,p2,q,amasses,pin1r,pin2r,p1r,p2r,qr)
!            call maptoonev2(ng,k,pin1b,pin2b,p1b,p2b,q,iclose,amasses
!     .           ,pin1r,pin2r,p1r,p2r,qr,iemap)
*??            if (iemap.gt.0) then
*??               return
*??            endif

            ksub = k
            call eikonal_factor('single',ng,pin1r,pin2r,p1r,p2r,qr,
     >           qph,eiki)
      
            if (boson.eq.'W') then
               oal=16.d0*elmat2form(pin1r,pin2r,p2r,p1r,qr,Qu,Qd,Qw,Qe)
            else
               oal=elmat2formz(pin1r,pin2r,p1r,p2r,qr,Qqu,Qel)
            endif
            
            ei1r = pin1r(0)
            ei2r = pin2r(0)
            ef3r = p1r(0)
            ef4r = p2r(0)
            if (iclose(k).eq.1) then 
               ef1r = ei1r - qr(0)
               xjaci = 1.d0/ei1r/2.d0
               xapi = qr(0)/ei1r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi = apxi/(1.d0 - xapi) ! attenzione !!
            endif
            if (iclose(k).eq.2) then 
               ef2r = ei2r - qr(0)
               xjaci = 1.d0/ei2r/2.d0
               xapi = qr(0)/ei2r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi = apxi/(1.d0 - xapi) ! attenzione !!
            endif
            if (iclose(k).eq.3) then 
               ei3r = ef3r + qr(0)
               xjaci = 1.d0/ei3r/2.d0
               xapi = qr(0)/ei3r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi = apxi/(1.d0 - xapi) ! attenzione !!
            endif
            if (iclose(k).eq.4) then 
               ei4r = ef4r + qr(0)
               xjaci = 1.d0/ei4r/2.d0
               xapi = qr(0)/ei4r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi = apxi/(1.d0 - xapi) ! attenzione !!
            endif
            approxi =  apxi*eiki*qr(0)*xjaci
            
c     if (oal.lt.0.d0) then
c               print*,'WARNING: O(alpha) < 0 \a \a!',ng,k,oal
c               oal = bornmel * approxi
c            endif
! subtraction!
            scale = sqrtscalesub
            sub = bornmel*approxi/eiki*
     .           subtractionnew(k,scale,amasses(1),amasses(2),
     .           charges(1),charges(2),qr)
            sub = sub*iplsub
            oal = oal - sub
            rescaletomatcheikfull = 1.d0 * eiki/singlephotons(k)
! rescaling by "rescaletomatcheikfull" is useful when subtraction is on,
! but when it is off the rescaling can make factcorr negative!!
            deltai = (oal - approxi*bornmel)/approxi *
     .                                       rescaletomatcheikfull
            delta  = delta + deltai
            factcorr = (1.d0 + deltai/bornmel)
            di(k) = deltai
            fc(k) = factcorr
            prod = prod * factcorr
            
ccc            if (factcorr.lt.0.d0) print*,k,ng,rescaletomatcheikfull

*** matched with eikonal
            deltai_eik = (oal/eiki - bornmel)
            delta_eik  = delta_eik + deltai_eik
            prod_eik = prod_eik * (1 + deltai_eik/bornmel)
***
         enddo

         summt2 = bornme*prodapx*eikfull*prodomega*xjac
         summt2_eik = eikfull*bornme*prod_eik

!! MATCHED PARTON SHOWER !
! EUREKA !!!
! Sun Aug 14 01:42:08 CEST 2005
         summt2 = summt2 * prod
!!!!!!!!!!!!!!!!!!!!!!!!
         imtx = 1
         return
      endif
*******************
      if (model.eq.'eikonal') then
         eikng = 'multip'
         ksub = -1000
         call eikonal_factor(eikng,ng,pin1,pin2,p1,p2,q,qph,eikonal)
         bornme = elmat2born(pin1b,pin2b,p2b,p1b)
                      !     the convention on the momenta is different! 
         summt2 = eikonal * bornme
         imtx = 1
         return
      endif
*
      if (model.eq.'ps') then
         if (ng.eq.0) then
            iref = 1
            bornme = elmat2born(pin1b,pin2b,p2b,p1b)
            summt2 = bornme
            imtx = 1
            return
         endif
*
         bornme = elmat2born(pin1b,pin2b,p2b,p1b)

         ei1 = pin1(0)
         ei2 = pin2(0)
         ef1 = ei1
         ef2 = ei2
         ef3 = p1(0)
         ef4 = p2(0)

         prodomega = 1.d0
         do k = 1,ng
            prodomega = prodomega*qph(k,0)
            if (iclose(k).eq.3) ef3 = ef3 + qph(k,0)
            if (iclose(k).eq.4) ef4 = ef4 + qph(k,0)
         enddo
         ei3 = ef3
         ei4 = ef4
         prodapx = 1.d0
         xjac    = 1.d0
         do k = 1,ng
** PS QUANTITIES
            if (iclose(k).eq.1) then
               xjac = xjac/ef1/2.d0
               xap(k) = qph(k,0)/ef1
               apx(k) = (1.d0 + (1.d0 - xap(k))**2)/xap(k)
               apx(k) = apx(k)/(1.d0 - xap(k)) ! attenzione !!
               ef1    = ef1 - qph(k,0)
            endif
            if (iclose(k).eq.2) then
               xjac = xjac/ef2/2.d0
               xap(k) = qph(k,0)/ef2
               apx(k) = (1.d0 + (1.d0 - xap(k))**2)/xap(k)
               apx(k) = apx(k)/(1.d0 - xap(k)) ! attenzione !!
               ef2    = ef2 - qph(k,0)
            endif
            if (iclose(k).eq.3) then
               xjac = xjac/ef3/2.d0
               xap(k) = qph(k,0)/ef3
               apx(k) = (1.d0 + (1.d0 - xap(k))**2)/xap(k)
               apx(k) = apx(k)/(1.d0 - xap(k)) ! attenzione !!
               ef3    = ef3 - qph(k,0)
            endif
            if (iclose(k).eq.4) then
               xjac = xjac/ef4/2.d0
               xap(k) = qph(k,0)/ef4
               apx(k) = (1.d0 + (1.d0 - xap(k))**2)/xap(k)
               apx(k) = apx(k)/(1.d0 - xap(k)) ! attenzione !!
               ef4    = ef4 - qph(k,0)
            endif
            prodapx = prodapx * apx(k)
****
         enddo
         if (ng.gt.1) then
            eikng = 'multip'
            ksub = -1000
            call eikonal_factor(eikng,ng,pin1,pin2,p1,p2,q,qph,eikfull)
         else
            do i=0,3
               q(i) = qph(1,i)
            enddo
            eikng = 'single'
            ksub  = 1
            call eikonal_factor(eikng,ng,pin1,pin2,p1,p2,q,qph,eikfull)
         endif
         summt2 = bornme*prodapx*eikfull*prodomega*xjac ! questo funziona !...
         imtx = 1
         return
      endif
*
      if (model.eq.'matchedps_RR') then
         if (ng.eq.0) then
            iref = 1
            bornme = elmat2born(pin1b,pin2b,p2b,p1b)
            summt2 = bornme
            imtx = 1
            return
         endif
         if (ng.eq.1) then
            do i=0,3
               q(i) = qph(1,i)
            enddo            
            summt2 = 16.d0*elmat2form(pin1,pin2,p2,p1,q,Qu,Qd,Qw,Qe)
            bornme = elmat2born(pin1b,pin2b,p2b,p1b)
            imtx = 1
            return            
         endif
*
         bornme = elmat2born(pin1b,pin2b,p2b,p1b)
         eikng = 'multip'
         ksub  = -1000
         call eikonal_factor(eikng,ng,pin1,pin2,p1,p2,q,qph,eikfull)

         ei1 = pin1(0)
         ei2 = pin2(0)
         ef1 = ei1
         ef2 = ei2
         ef3 = p1(0)
         ef4 = p2(0)
         prodomega = 1.d0
         do k = 1,ng
            prodomega = prodomega*qph(k,0)
            if (iclose(k).eq.3) ef3 = ef3 + qph(k,0)
            if (iclose(k).eq.4) ef4 = ef4 + qph(k,0)
         enddo
         ei3 = ef3
         ei4 = ef4

         prodapx = 1.d0
         xjac    = 1.d0
         delta   = 0.d0
         prod    = 1.d0
         delta_eik   = 0.d0
         prod_eik    = 1.d0
         do k = 1,ng
            do i=0,3
               q(i) = qph(k,i)
            enddo            
**** PS QUANTITIES
            if (iclose(k).eq.1) then
               xjac = xjac/ef1/2.d0
               xap(k) = qph(k,0)/ef1
               apx(k) = (1.d0 + (1.d0 - xap(k))**2)/xap(k)
               apx(k) = apx(k)/(1.d0 - xap(k)) ! attenzione !!
               ef1    = ef1 - qph(k,0)
            endif
            if (iclose(k).eq.2) then
               xjac = xjac/ef2/2.d0
               xap(k) = qph(k,0)/ef2
               apx(k) = (1.d0 + (1.d0 - xap(k))**2)/xap(k)
               apx(k) = apx(k)/(1.d0 - xap(k)) ! attenzione !!
               ef2    = ef2 - qph(k,0)
            endif
            if (iclose(k).eq.3) then
               xjac = xjac/ef3/2.d0
               xap(k) = qph(k,0)/ef3
               apx(k) = (1.d0 + (1.d0 - xap(k))**2)/xap(k)
               apx(k) = apx(k)/(1.d0 - xap(k)) ! attenzione !!
               ef3    = ef3 - qph(k,0)
            endif
            if (iclose(k).eq.4) then
               xjac = xjac/ef4/2.d0
               xap(k) = qph(k,0)/ef4
               apx(k) = (1.d0 + (1.d0 - xap(k))**2)/xap(k)
               apx(k) = apx(k)/(1.d0 - xap(k)) ! attenzione !!
               ef4    = ef4 - qph(k,0)
            endif
            prodapx = prodapx * apx(k)
*****
            is = 0
            if (iclose(k).eq.1.or.iclose(k).eq.2) is = 1           


            iemap = 0
            call maptoone(is,qq,p1,p2,q,amasses,pin1r,pin2r,p1r,p2r,qr)
!            call maptoonev2(ng,k,pin1b,pin2b,p1b,p2b,q,iclose,amasses
!     .           ,pin1r,pin2r,p1r,p2r,qr,iemap)

            if (iemap.gt.0) then
               return
            endif
            ksub = k
            call eikonal_factor('single',ng,pin1r,pin2r,p1r,p2r,qr,
     >           qph,eiki)
      
            oal = 16.d0*elmat2form(pin1r,pin2r,p2r,p1r,qr,Qu,Qd,Qw,Qe)

            ei1r = pin1r(0)
            ei2r = pin2r(0)
            ef3r = p1r(0)
            ef4r = p2r(0)
            if (iclose(k).eq.1) then 
               ef1r = ei1r - qr(0)
               xjaci = 1.d0/ei1r/2.d0
               xapi = qr(0)/ei1r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi = apxi/(1.d0 - xapi) ! attenzione !!
            endif
            if (iclose(k).eq.2) then 
               ef2r = ei2r - qr(0)
               xjaci = 1.d0/ei2r/2.d0
               xapi = qr(0)/ei2r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi = apxi/(1.d0 - xapi) ! attenzione !!
            endif
            if (iclose(k).eq.3) then 
               ei3r = ef3r + qr(0)
               xjaci = 1.d0/ei3r/2.d0
               xapi = qr(0)/ei3r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi = apxi/(1.d0 - xapi) ! attenzione !!
            endif
            if (iclose(k).eq.4) then 
               ei4r = ef4r + qr(0)
               xjaci = 1.d0/ei4r/2.d0
               xapi = qr(0)/ei4r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi = apxi/(1.d0 - xapi) ! attenzione !!
            endif

            approxi =  apxi*eiki*qr(0)*xjaci

            if (oal.lt.0.d0) then
               print*,'WARNING: O(alpha) < 0 \a \a!',ng,k,oal
               oal = bornme * approxi
            endif

!            deltail = (oal/approxi - bornmel) 
            deltai = (oal/approxi - bornme)
            delta  = delta + deltai
            factcorr = (1.d0 + deltai/bornme)
            di(k) = deltai
            fc(k) = factcorr
            prod = prod * factcorr
*** matched with eikonal
            deltai_eik = (oal/eiki - bornme) 
            delta_eik  = delta_eik + deltai_eik
            prod_eik = prod_eik * (1 + deltai_eik/bornme)
***
         enddo            
         summt2 = bornme*prodapx*eikfull*prodomega*xjac
         summt2_eik = eikfull*bornme*prod_eik

!! MATCHED PARTON SHOWER !!!!!!!!!
! EUREKA !!! 
! Sun Aug 14 01:42:08 CEST 2005
!
         summt2 = summt2 * prod
!
!!!!!!!!!!!!!!!!!!!!!!!!
         imtx = 1
         return
      endif
************************
*** ALPHA
c$$$      if (model.eq.'alpha') then
c$$$         if (ng.eq.0) then
c$$$            iref = 1
c$$$            bornme = elmat2born(pin1b,pin2b,p2b,p1b)
c$$$            summt2 = bornme
c$$$            imtx = 1
c$$$            return
c$$$         endif
c$$$         if (ng.eq.1) then
c$$$            do i=0,3
c$$$               q(i) = qph(1,i)
c$$$            enddo            
c$$$            summt2 = 16.d0*elmat2form(pin1,pin2,p2,p1,q,Qu,Qd,Qw,Qe)
c$$$!            call alpha_me('R',ng,p2,p1,qph,imtx,emtxalpha)
c$$$!            summt2 = emtxalpha
c$$$            bornme = elmat2born(pin1b,pin2b,p2b,p1b)
c$$$            imtx = 1
c$$$            return            
c$$$         endif
c$$$         call alpha_me('R',ng,p2,p1,qph,ialpha,emtxalpha)
c$$$         imtx   = ialpha
c$$$         summt2 = emtxalpha
c$$$         if (imtx.eq.0) iemap = 0
c$$$         return
c$$$      endif
**********************
      return
      end
*
*** SOFT + VIRTUAL
      subroutine svfactor(model,ng,ecms,p1,p2,eps,sv)
      implicit double precision (a-h,l,m,o-z)
      dimension p1(0:3),p2(0:3),pin1(0:3),pin2(0:3),p1b(0:3),p2b(0:3)
      dimension pin1p(0:3),pin2p(0:3),p1p(0:3),p2p(0:3)
      dimension qphempty(40,0:3)
      double complex duedv,factor,duedvres,wselfp,Vfffp,dboxfp,selfenfp
      double complex dww,dbox,vff,vffnonresonant,ctsnew
      double complex dwwnonresonant,dboxnonresonant
      external dww,dbox,vff,vffnonresonant,ctsnew,Vfffp,dboxfp,selfenfp
      external dwwnonresonant,dboxnonresonant,wselfp
      character*10 model	
      character*1 boson
      double precision mq(-5:5),chq(-5:5),itre(-5:5)
      logical varlog
      dimension locmasses(4),loccharges(4)
      double precision masses(4),ptmp(0:3)
      character*6 ord
      common/qedORDER/ord
      common/svfordebug/svdebug,emtx,emtxsub,wnphot
      common/momentainitial/pin1,pin2
      common/reducedtoborn/p1b,p2b,iref
      common/momentainitialred/pin1b(0:3),pin2b(0:3),pin1r(0:3)
     >     ,pin2r(0:3)
      common/forborncrosssection/phsp2b,flux2b,bornme,bornmeclean
      common /regulators/ESOFTMAX
      common/chargesmasses_cm/charges(4),amasses(4)
      common/counter_svfactor/icountsv
      data icountsv /0/
!from shared.inc
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/t_widths/gw,gz
      common/vectorboson/boson
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt
     .     ,mq,chq,itre
      common/iplsubtraction/iplsub
      common/subtractedsv/subtractedoalsv
      common/scalapersottrazione/sqrtscalesub
** for virtual libraries
      common/kinvars/horaces,horacet,horaceu
      common/ddkinvars/dds,ddt,ddu
      common/uukinvars/uus,uut,uuu
      character*5 zIS
      common/selectvirtualz/zIS
      common/partons/ipart1,ipart2

***
      common/debugging/lcoll1,lcoll2
      common/ifirstsvf/ifirst
      data ifirst /0/
      
      data qphempty /160*0.d0/
      icountsv = icountsv + 1

      subtractedoalsv = 0.d0
      if (ipart1.eq.0.and.ipart2.eq.0) then
         sv = 1.d0
         return
      endif
      if (ord.eq.'born'.or.(ord.eq.'alpha'.and.ng.eq.1)) then
         sv = 1.d0
         return
      endif
*.+.+.+.+.+.+.+.+
      do k = 0,3
         ptmp(k) = pin1b(k) + pin2b(k)
      enddo
      s = dot(ptmp,ptmp)
      horaces = s
      do k = 0,3
         ptmp(k) = pin1b(k) - p2b(k)
      enddo
      horacet = dot(ptmp,ptmp)
      do k = 0,3
         ptmp(k) = pin1b(k) - p1b(k)
      enddo
      horaceu = dot(ptmp,ptmp)

      dds = horaces
      uus = horaces
      uut = horaceu
      uuu = horacet
      ddt = horaceu
      ddu = horacet
      bcks  = s
      bckhs = horaces
      bckht = horacet
      bckhu = horaceu
cccccccc Fixed  in common    ESOFTMAX = eps * sqrt(s)/2.d0
      varlog        = .TRUE.
      loccharges(1) = charges(1)
      loccharges(2) = charges(2)
      loccharges(3) = charges(4)
      loccharges(4) = charges(3)
      locmasses(1)  = amasses(1)
      locmasses(2)  = amasses(2)
      locmasses(3)  = amasses(4)
      locmasses(4)  = amasses(3)

      if (ord.eq.'alpha') then
         if (ng.ge.1) then
            sv = 1.d0
            return
         endif
! here only if ng > 0
         
         if (model.eq.'ps'.or.model.eq.'eikonal') then 
            epsa = 2.d0*esoftmax/sqrt(s)
            aieps  = -2.d0*log(epsa)-1.5d0+2.d0*epsa-0.5d0*epsa**2
            lcoll1 = eikonalintegral(pin1,pin2,p1,p2,charges,amasses)
! TEST SUBTRACTION alpha
            subscale = sqrtscalesub**2
            lcoll = loccharges(1)**2*(log(subscale/amasses(1)**2)- 1.d0)
     .           + loccharges(2)**2*(log(subscale/amasses(2)**2)- 1.d0)
            lcoll1 = lcoll1 - iplsub*lcoll
! END TEST SUBTRACTION alpha
            argffs1= alpha/2.d0/pi * lcoll1 * aieps               
            svapp = 1.d0 - argffs1
            sv = svapp
            return
         endif            

         varlog = .TRUE.
         virtual = 0.d0
         loop = 0.d0
         if (boson.eq.'W') then
            call squared_me(tree,loop,s,0,0,0,0,0,0,0,0,varlog)
*[[[[
* PS APPROX
*            epsa = 2.d0*esoftmax/sqrt(s)
*            aieps  = -2.d0*log(epsa)-1.5d0+2.d0*epsa-0.5d0*epsa**2
*            lcoll1 = eikonalintegral(pin1,pin2,p1,p2,charges,amasses)
*            argffs1= alpha/2.d0/pi * lcoll1 * aieps               
*            svapp = 1.d0 - argffs1
*            sv = svapp
*]]]]
         else                   !! Z

*[[[[[
* PS APPROX with subtraction
c$$$            epsa = eps
c$$$            aieps  = -2.d0*log(epsa)-1.5d0+2.d0*epsa-0.5d0*epsa**2
c$$$            lcoll1 = eikonalintegral(pin1,pin2,p1,p2,charges,amasses)
c$$$            argffs1= alpha/2.d0/pi * lcoll1 * aieps               
c$$$            svapp = 1.d0 - argffs1
c$$$            sv = svapp
c$$$            
c$$$            aieps  = -2.d0*log(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
c$$$            subscale = sqrtscalesub**2
c$$$            lcoll = loccharges(1)**2*(log(subscale/amasses(1)**2)- 1.d0)
c$$$     .           + loccharges(2)**2*(log(subscale/amasses(2)**2)- 1.d0)
c$$$            argffs1= alpha/2.d0/pi * lcoll * aieps            
c$$$            svsub = -argffs1
c$$$            subtractedoalsv = -iplsub*svsub ! in common
c$$$            sv    = sv - iplsub * svsub
c$$$            svdebug = sv
c$$$            return
*]]]]]
            if (zIS.eq.'uubar')
     .           call uusquared_me(tree,loop,s,0,0,0,0,0,0,0,0,varlog)
            if (zIS.eq.'ddbar')
     .           call ddsquared_me(tree,loop,s,0,0,0,0,0,0,0,0,varlog)
         endif

c         print*,loop
c         if (icountsv.eq.10) stop
         
         riscala = 1.d0
         if (boson.eq.'Z') riscala = bornmeclean/bornme
         
         virtual = loop*16.d0   ! 16.d0 comes from tree/loop !!
         phmass = sqrt(getlambda())
         softff = soft_integral(phmass,esoftmax,pin1,pin2,p2,p1,
     >        loccharges,locmasses)*alpha/2.d0/pi
         sv = (bornme*(1.d0 + riscala*softff) + virtual)/bornme

         aieps  = -2.d0*log(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
         subscale = sqrtscalesub**2
         lcoll = loccharges(1)**2*(log(subscale/amasses(1)**2)- 1.d0)
     .        + loccharges(2)**2*(log(subscale/amasses(2)**2)- 1.d0)
         argffs1= alpha/2.d0/pi * lcoll * aieps            
         svsub = -argffs1
         subtractedoalsv = -iplsub*svsub ! in common
         sv    = sv - iplsub * svsub * riscala
         svdebug = sv

c         print*,'>>>>>>>>>>>>',sv,svsub,softff

         return
      endif

cc WITHOUT MASSES
c           softdk =  softffDK(esoftmax,phmass,locmasses(1),locmasses(2),
c     .           locmasses(4),
c     .           horaces,horacet,horaceu,2.d0/3.d0,-1.d0/3.d0,-1.d0,
c     .           alpha)
c            softff = softdk
cc WITHOUT MASSES

      if (ord.eq.'exp') then
*.+.+.+.+.+.+.+.+
         do k = 0,3
! nominal s !
c            pin1p(k) = pin1(k)
c            pin2p(k) = pin2(k)
! end
! going at reduced s!
            pin1p(k) = pin1b(k)
            pin2p(k) = pin2b(k)
            p1p(k)   = p1b(k)
            p2p(k)   = p2b(k)
! end
            ptmp(k)  = pin1p(k) + pin2p(k)
         enddo
         s = dot(ptmp,ptmp)
         horaces = s
         ei = sqrt(s)
C below only for "nominal s"
c$$$         ea = (ei**2 + amasses(3)**2 -amasses(4)**2)/2.d0/ei
c$$$         eb = ei - ea
c$$$         beta1 = sqrt(1.d0 - amasses(3)**2/ea**2)
c$$$         beta2 = sqrt(1.d0 - amasses(4)**2/eb**2)
c$$$         p1p(0) = ea
c$$$         p1p(1) = ea*beta1*p1b(1)/sqrt(tridot(p1b,p1b))
c$$$         p1p(2) = ea*beta1*p1b(2)/sqrt(tridot(p1b,p1b))
c$$$         p1p(3) = ea*beta1*p1b(3)/sqrt(tridot(p1b,p1b))
c$$$         p2p(0) = eb
c$$$         p2p(1) = -p1p(1)
c$$$         p2p(2) = -p1p(2)
c$$$         p2p(3) = -p1p(3)
         do k = 0,3
            ptmp(k) = pin1p(k) - p2p(k)
         enddo
         horacet = dot(ptmp,ptmp)
         do k = 0,3
            ptmp(k) = pin1p(k) - p1p(k)
         enddo
         horaceu = dot(ptmp,ptmp)

         dds = horaces
         uus = horaces
         uut = horaceu
         uuu = horacet
         ddt = horaceu
         ddu = horacet
         aieps = -2.d0*log(eps)-1.5d0+2.d0*eps-0.5d0*eps**2

         lcoll1  = eikonalintegral(pin1,pin2,p1,p2,charges,amasses)
         lcoll2  = eikonalintegral(pin1p,pin2p,p1p,p2p,charges,amasses)

         if (boson.eq.'W') then
            bornmel = elmat2born(pin1p,pin2p,p2p,p1p)
         else
            bornmel = elmat2bornz(pin1p,pin2p,p1p,p2p,0)
         endif
         
         argffs1 = alpha/2.d0/pi * lcoll1 * aieps
         argffs2 = alpha/2.d0/pi * lcoll2 * aieps 
** PS APPROX
c         if (boson.eq.'Z') then
c            sv = exp(-argffs1)
c            return
c         endif
***
         if (model.eq.'ps'.or.model.eq.'eikonal') then
            subscale = sqrtscalesub**2
            lcoll = loccharges(1)**2*(log(subscale/amasses(1)**2)- 1.d0)
     .            + loccharges(2)**2*(log(subscale/amasses(2)**2)- 1.d0)
            lcoll1 = lcoll1 - lcoll*iplsub
            argffs1 = alpha/2.d0/pi * lcoll1 * aieps
            sv = exp(-argffs1)
            svdebug = sv
            return
         endif
         subscale = sqrtscalesub**2
         lcollsub = loccharges(1)**2*(log(subscale/amasses(1)**2)- 1.d0)
     .        + loccharges(2)**2*(log(subscale/amasses(2)**2)- 1.d0)
         lcollsub = lcollsub*iplsub
         tbsub = alpha/2.d0/pi * lcollsub * aieps * iplsub
         subtractedoalsv = iplsub*tbsub

*#######################################
c         goto 111 ! PS APPROX FOR SV!!!
*#######################################

         varlog = .TRUE.
         if (boson.eq.'W') then
            call squared_me(tree,loop,s,0,0,0,0,0,0,0,0,varlog)
         else
            if (zIS.eq.'uubar')
     .           call uusquared_me(tree,loop,s,0,0,0,0,0,0,0,0,varlog)
            if (zIS.eq.'ddbar')
     .           call ddsquared_me(tree,loop,s,0,0,0,0,0,0,0,0,varlog)
         endif
         virtual = loop*16.d0   ! 16.d0 comes from tree/bornme !!

         phmass = sqrt(getlambda())
         softff = soft_integral(phmass,esoftmax,pin1p,pin2p,p2p,p1p,
     >        loccharges,locmasses)*alpha/2.d0/pi
         vplussoft = bornmel*softff + virtual
         argex    = vplussoft/bornmel
         deltasv  = argex + argffs2

*#######################################
c 111     deltasv = 0.d0 ! PS APPROX FOR SV!!!
*#######################################

         svcorr = 1.d0+deltasv

ccc   svcorr = exp(deltasv) !!!
         
         svdebug = svcorr
         if (svcorr.lt.0.d0.or.svcorr.gt.2.d0) then
c$$$            print*,'SV WARNING WARNING:'
c$$$            print*,'nphot= ',ng
c$$$            print*,'deltasv=',deltasv
c$$$            print*,'virtual=',virtual
c$$$            print*,'bornme= ',bornmel,bornmel-16.d0*tree
c$$$            print*,'esoftmax= ',esoftmax
c$$$            print*,'sqrt(shat)= ',sqrt(horaces)
c$$$            print*,'sqrt(-that) = ',sqrt(-horacet)
c$$$            print*,'sqrt(-uhat) =',sqrt(-horaceu)
c$$$            print*,'cos =',(2.d0*horacet + horaces)/horaces
c$$$            print*,pin1p
c$$$            print*,pin2p
c$$$            print*,p1p
c$$$            print*,p2p
c$$$            print*,'s+t+u =',horaces+horacet+horaceu
            if (ifirst.lt.10) then
               print*,'matching.f line 2035'               
               print*,'<<>> svcorr negative or large<<>>',
     .              svcorr,sqrt(horaces)
               ifirst = ifirst + 1
               print*,'NOT RESETTING TO 1 as in the default!'
            endif
c            svcorr = 1.d0
            call dump_for_debug(3,p1,p2,qphempty)
         endif
         expon = exp(-(argffs1-tbsub))
         sv = expon * svcorr

c         if (isnan(sv)) then
c            print*,'line 1928 matching.f'
c            print*,sv,deltasv,softff
c         endif
         
         svdebug = svdebug * expon
         return
         endif
      return
      end
***** 
      function eikonalintegral(p1,p2,p3,p4,ch,masses)
      implicit double precision (a-h,l,m,o-z)
      character*1 boson
      common/vectorboson/boson
      parameter (npart = 4)
      dimension eta(npart),p(0:3),q(0:3),pmat(npart,0:3)
      dimension masses(npart),ch(npart)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
! ch are the charges of the field (not anti-field!!)
! Use this convention: the integer factor in front of the charge
! must be
! -1 --> for incoming particle
! -1 --> for outgoing anti-particle
! +1 --> for outgoing particle
! +1 --> for incoming anti-particle
      if (boson.eq.'W') then
!! It is correct for W+ !!!!!!!!!!
         eta(1) = -1.d0 * ch(1)
         eta(2) =  1.d0 * ch(2)
         eta(3) = -1.d0 * ch(3)
         eta(4) =  1.d0 * ch(4)
      else
         eta(1) = -1.d0 * ch(1)
         eta(2) =  1.d0 * ch(2)
         eta(3) =  1.d0 * ch(3)
         eta(4) = -1.d0 * ch(4)
      endif
      do k = 0,3
         pmat(1,k) = p1(k)
         pmat(2,k) = p2(k)
         pmat(3,k) = p3(k)
         pmat(4,k) = p4(k)
      enddo
! off-diagonal contributions
      softint = 0.d0
      do i = 1,npart-1
         do j = i+1,npart
            etaij = eta(i) * eta(j)
            if (abs(etaij).gt.1.d-3) then

ccc OLD: gives NaNs               
c               call rescale_momenta(npart,i,j,pmat,masses,p,q)
c               q2 = dot(q,q)
c               vl = dot(p,p) - q2
c               vl = vl/2.d0
c               v  = vl/(p(0) - q(0))
c               arglog = 1.d0+2.d0*vl/q2
c               terminfra = log(arglog)
c               tot = terminfra  ! + termfinite(p,q,v)
c               tot = - tot *dot(p,q)/vl* etaij
c               softint = softint + tot
c               totbck = tot
ccc   
********************** new: from Hto4l
               call rescale_momenta_new(npart,i,j,pmat,masses,p,q,rho)
               q2 = masses(j)*masses(j)
               vl = 0.5d0*(rho*rho*masses(i)*masses(i) - q2)
               v  = vl/(p(0) - q(0)) 
               arglog = 1.d0+2.d0*vl/q2               
               if (arglog.gt.0.d0) then 
                  tot = log(arglog)
               else
                  tot = 0.d0
               endif

               tot = -tot *dot(p,q)/vl* etaij

c               print*,tot,totbck,tot/totbck
               
               softint  = softint + tot
**************************
               
            endif
         enddo
      enddo
      soft_integral = softint

      
!     diagonal contributions
c$$$      softint = 0.d0
c$$$      do  i = 1,npart
c$$$         etaii = eta(i)**2
c$$$         if (etaii.gt.1.d-3) then
c$$$            betai=(1.d0-masses(i)/pmat(i,0))*(1.d0+masses(i)/pmat(i,0))
c$$$            betai=sqrt(betai)
c$$$!     m2se2 = masses(i)**2/pmat(i,0)**2
c$$$!     umb2  = (1.d0 - betai)*(1.d0 + betai)
c$$$!     m2se2/umb2 is always = 1 !!!!
c$$$            term = log((1.d0-betai)/(1.d0+betai))/betai 
c$$$            term = -2.d0 * term * etaii
c$$$            softint = softint + term
c$$$         endif
c$$$      enddo
      softint = -1.d0 * (eta(1)**2+eta(2)**2+eta(3)**2+eta(4)**2)
      soft_integral = soft_integral + softint
      eikonalintegral = soft_integral

      return
      end
      
*** ROUTINES
      subroutine closest(check,ng,qph,p1,p2,p3,p4,q1,q2,q3,q4,iclose)
      implicit double precision (a-h,o-z)
      parameter (maxph = 40)
      character*1 boson
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension q1(0:3),q2(0:3),q3(0:3),q4(0:3)
      dimension qph(maxph,0:3)
      character*4 check,check2
      dimension iclose(maxph)
      integer emission(maxph)
      common/radpattern/nph(4)
      common/vectorboson/boson
* FOR NEW SUBTRACTION (superseded...)
      common/cosinesforsubtraction/cosu,cosd
      common/closesubtraction/iclosesub(maxph)
****
      if (check.ne.'full'.and.check.ne.'inst') then
         print*,'In subroutine "closest", wrong value of check!!'
         stop
      endif
      do k = 1,maxph
         iclose(k)    = 0
         iclosesub(k) = 0
      enddo
      check2 = check ! to avoid segmentation fault...

      do k = 0,3
         q1(k) = 0.d0
         q2(k) = 0.d0
         q3(k) = 0.d0
         q4(k) = 0.d0
      enddo

      p1mod = sqrt(tridot(p1,p1))
      p2mod = sqrt(tridot(p2,p2))
      p3mod = sqrt(tridot(p3,p3))
      p4mod = sqrt(tridot(p4,p4))
      
      do k = 1,ng
         do i = 0,3
            q(i) = qph(k,i)
         enddo
ccc
         if (k.le.nph(1)) emission(k) = 1
         if (k.gt.nph(1).and.k.le.(nph(1)+nph(2))) emission(k) = 2
         if (k.gt.(nph(1)+nph(2)).and.k.le.(nph(1)+nph(2)+nph(3)))
     .        emission(k) = 3
         if (k.gt.(nph(1)+nph(2)+nph(3)))
     .        emission(k) = 4
ccc
         c1 = tridot(p1,q)/p1mod/q(0)
         c2 = tridot(p2,q)/p2mod/q(0)
         c3 = tridot(p3,q)/p3mod/q(0)
         c4 = tridot(p4,q)/p4mod/q(0)

         if (check2.eq.'inst') then
            if (c1.gt.c2) iclose(k) = 1
            if (c2.gt.c1) iclose(k) = 2
         else
            if (boson.eq.'W') then
               if (c1.gt.c2.and.c1.gt.c3) iclose(k) = 1
               if (c2.gt.c1.and.c2.gt.c3) iclose(k) = 2
               if (c3.gt.c2.and.c3.gt.c1) iclose(k) = 3
            else
               if (c1.gt.c2.and.c1.gt.c3.and.c1.gt.c4) iclose(k) = 1
               if (c2.gt.c1.and.c2.gt.c3.and.c2.gt.c4) iclose(k) = 2
               if (c3.gt.c2.and.c3.gt.c1.and.c3.gt.c4) iclose(k) = 3
               if (c4.gt.c1.and.c4.gt.c2.and.c4.gt.c3) iclose(k) = 4
            endif
            if (c1.gt.cosu) iclosesub(k) = 1
            if (c1.lt.cosd) iclosesub(k) = 2
         endif
ccc         iclose(k) = emission(k)
         ic = iclose(k)
         if (ic.eq.1) then
            do j = 0,3
               q1(j) = q1(j) + q(j)
            enddo
         endif            
         if (ic.eq.2) then
            do j = 0,3
               q2(j) = q2(j) + q(j)
            enddo
         endif
         if (ic.eq.3) then
            do j = 0,3
               q3(j) = q3(j) + q(j)
            enddo
         endif
         if (ic.eq.4) then
            do j = 0,3
               q4(j) = q4(j) + q(j)
            enddo
         endif
      enddo
      return
      end
***
      subroutine maptozero(p1,p2,p3,p4,masses,p1r,p2r,p3r,p4r,iemap)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension p1r(0:3),p2r(0:3),p3r(0:3),p4r(0:3)
      dimension q(0:3),p(0:3)
      double precision masses(4)

      iemap = 0

      do k=0,3
         q(k) = p1(k) + p2(k)
         p(k) = p3(k) + p4(k)
      enddo

      call new_boost(p,p3,p3r,1)
      call new_boost(p,p4,p4r,1)
      
      sshat = sqrt(dot(q,q)) ! invariant

      e = sshat/2.d0

      if (e.lt.(masses(1)+masses(2))) then
         e = max(masses(1),masses(3))+max(masses(2),masses(3)) + e/1d3
         iemap = 1
      endif

      ei   = 2.d0*e

      ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
      eb = ei - ea

      if (ea.le.masses(1).or.eb.le.masses(2)) then
! here it should never pass
         iemap = 1
         if (ea.le.masses(1)) ea = max(masses(1),masses(3))
         if (eb.le.masses(2)) eb = max(masses(2),masses(3))
         ei = ea + eb
      endif

      beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
      beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

      p1r(0) = ea
      p1r(1) = 0.d0
      p1r(2) = 0.d0
      p1r(3) = ea*beta1

      p2r(0) =  eb
      p2r(1) =  0.d0
      p2r(2) =  0.d0
      p2r(3) = -eb*beta2

      abig = ei**2 + masses(3)**2 - masses(4)**2
      bbig = -2.d0 * ei 
      cbig = 0.d0
      arg  = abig**2*bbig**2-masses(3)**2*bbig**2*(bbig**2-cbig**2)

      p1mod = abig*cbig + sqrt(arg)
      p1mod = p1mod / (bbig**2 - cbig**2)

      e1    = sqrt(masses(3)**2 + p1mod**2)
      beta3 = p1mod/e1  

      p3mod = sqrt(tridot(p3r,p3r))
      p4mod = sqrt(tridot(p4r,p4r))

      p3r(0) = e1
      p3r(1) = e1*beta3 * p3r(1)/p3mod
      p3r(2) = e1*beta3 * p3r(2)/p3mod
      p3r(3) = e1*beta3 * p3r(3)/p3mod

      p4r(0) = ei - p3r(0)
      p4r(1) = -p3r(1)
      p4r(2) = -p3r(2)
      p4r(3) = -p3r(3)

      return
      end
*******
      subroutine maptoone(is,bigq,p3,p4,q,masses,p1r,p2r,p3r,p4r,qr)
      implicit double precision (a-h,o-z)
      dimension bigq(0:3),p3(0:3),p4(0:3)
      dimension p1r(0:3),p2r(0:3),p3r(0:3),p4r(0:3)
      dimension q(0:3),p(0:3),qr(0:3),ptmp(0:3)
      double precision masses(4)
      common/momentainitial/pin1(0:3),pin2(0:3)
      common/counter_mapptoone/icountmto
      data icountmto /0/
      icountmto = icountmto + 1
! is = 1 q is initial state photon
! is = 0 q is final state photon
      shat = dot(bigq,bigq)
      if (is.eq.1) then
         do k = 0,3
            p(k) = p3(k)+p4(k)
         enddo
         call new_boost(p,p3,p3r,1)
         call new_boost(p,p4,p4r,1)
         e = sqrt(shat)/2.d0

         ei   = 2.d0*e
         abig = ei**2 + masses(3)**2 - masses(4)**2
         bbig = -2.d0 * ei
         cbig = 0.d0
         arg  = abig**2*bbig**2-masses(3)**2*bbig**2*(bbig**2-cbig**2)

         p1mod = abig*cbig + sqrt(arg)
         p1mod = p1mod / (bbig**2 - cbig**2)

         e1    = sqrt(masses(3)**2 + p1mod**2)
         beta3 = p1mod/e1

         p3mod = sqrt(tridot(p3r,p3r))
         p4mod = sqrt(tridot(p4r,p4r))

         p3r(0) = e1
         p3r(1) = e1*beta3 * p3r(1)/p3mod
         p3r(2) = e1*beta3 * p3r(2)/p3mod
         p3r(3) = e1*beta3 * p3r(3)/p3mod

         p4r(0) = ei - p3r(0)
         p4r(1) = -p3r(1)
         p4r(2) = -p3r(2)
         p4r(3) = -p3r(3)
         ! p3r p4r to be intended in the bigq frame

         do k = 0,3
            p(k)    = bigq(k)+q(k)
            ptmp(k) = pin1(k)+pin2(k)-q(k)
         enddo

         scale = sqrt(dot(bigq,bigq)/dot(ptmp,ptmp))

         do k = 0,3
            qr(k) = scale*q(k)
            ptmp(k) = pin1(k)+pin2(k)
         enddo

         e = sqrt(dot(ptmp,ptmp))/2.d0*scale

         ei = 2.d0*e
         ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
         eb = ei - ea
         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

         p1r(0) = ea
         p1r(1) = 0.d0
         p1r(2) = 0.d0
         p1r(3) = ea*beta1

         p2r(0) = eb
         p2r(1) = 0.d0
         p2r(2) = 0.d0
         p2r(3) = -eb*beta2

         do k = 0,3
            ptmp(k) = p1r(k)+p2r(k)-qr(k)
         enddo

 !        call new_boost(ptmp,qr,qr,1)
         call new_boost(ptmp,p3r,p3r,-1)
         call new_boost(ptmp,p4r,p4r,-1)

      else
         e=sqrt(shat)/2.d0

         ei = 2.d0*e
         ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
         eb = ei - ea
         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

         p1r(0) = ea
         p1r(1) = 0.d0
         p1r(2) = 0.d0
         p1r(3) = ea*beta1

         p2r(0) =  eb
         p2r(1) =  0.d0
         p2r(2) =  0.d0
         p2r(3) = -eb*beta2

         do k = 0,3
            ptmp(k) = p3(k)+p4(k)
         enddo

         call new_boost(ptmp,p3,p3r,1)
         call new_boost(ptmp,p4,p4r,1)
         call new_boost(ptmp,q,qr,1)
         
         do k = 0,3
            ptmp(k) = p3(k)+p4(k)+q(k)
         enddo

         scale = sqrt(dot(bigq,bigq)/dot(ptmp,ptmp))

         do k = 0,3
            qr(k) = scale*qr(k)
            ptmp(k) = (p4(k)+p3(k))*scale
         enddo

         e = sqrt(dot(ptmp,ptmp))/2.d0

         ei   = 2.d0*e
         abig = ei**2 + masses(3)**2 - masses(4)**2
         bbig = -2.d0 * ei
         cbig = 0.d0
         arg  = abig**2*bbig**2-masses(3)**2*bbig**2*(bbig**2-cbig**2)

         p1mod = abig*cbig + sqrt(arg)
         p1mod = p1mod / (bbig**2 - cbig**2)

         e1    = sqrt(masses(3)**2 + p1mod**2)
         beta3 = p1mod/e1

         p3mod = sqrt(tridot(p3r,p3r))
         p4mod = sqrt(tridot(p4r,p4r))

         p3r(0) = e1
         p3r(1) = e1*beta3 * p3r(1)/p3mod
         p3r(2) = e1*beta3 * p3r(2)/p3mod
         p3r(3) = e1*beta3 * p3r(3)/p3mod

         p4r(0) = ei - p3r(0) 
         p4r(1) = -p3r(1)
         p4r(2) = -p3r(2)
         p4r(3) = -p3r(3)
         ! p3r p4r to be intended in the current p frame

         do k = 0,3
            ptmp(k) = p3r(k)+p4r(k)+qr(k)
         enddo
         call new_boost(ptmp,p3r,p3r,1)
         call new_boost(ptmp,p4r,p4r,1)
         call new_boost(ptmp,qr,qr,1)
      endif
      return
      end
*******
      subroutine eikonal_factor(eikng,ng,p1,p2,p3,p4,q,qph,eikonal)
      implicit double precision (a-h,o-z)
      double precision me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     mq(-5:5),chq(-5:5),itre(-5:5)
      character*6 eikng
      character*1 boson
      parameter (imaxph = 40)
      dimension iclose(imaxph)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension qph(imaxph,0:3),q(0:3),peik(0:3),eta(4)
      common/chargesmasses_cm/charges(4),amasses(4)
      common/kforsub/ksub
      common/scalapersottrazione/sqrtscalesub
      common/singleseikonals/singlephotons(40)
! from shared.inc
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,mq,chq,itre
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/vectorboson/boson
      common/iplsubtraction/iplsub
! Use this convention: the integer factor in front of the charge
! must be
! -1 --> for incoming particle
! -1 --> for outgoing anti-particle
! +1 --> for outgoing particle
! +1 --> for incoming anti-particle
      if (boson.eq.'W') then
!! this is for W+ !!!!!!!!!!!!
         eta(1) = -1.d0 * charges(1)
         eta(2) =  1.d0 * charges(2)
         eta(3) = -1.d0 * charges(3)
         eta(4) =  1.d0 * charges(4)
      else
         eta(1) = -1.d0 * charges(1)
         eta(2) =  1.d0 * charges(2)
         eta(3) =  1.d0 * charges(3)
         eta(4) = -1.d0 * charges(4)
      endif

      eikonal = 1.d0
      if (ng.eq.0) return

      nglocal = ng
      if (eikng.eq.'single') nglocal = 1

      ec2 = 4.d0*pi*alpha ! electron charge^2
      do k = 1,nglocal
         if (eikng.eq.'multip') then
            do i = 0,3
               q(i) = qph(k,i)
            enddo
         endif
         oop1dk = 1.d0/dot(p1,q) * eta(1)
         oop2dk = 1.d0/dot(p2,q) * eta(2)

         if (boson.eq.'W') then
            oop4dk = 0.d0

! trying to avoid NaNs at extremely high energies
            oop3dk = dot(p3,q) / eta(3)
            if (oop3dk.eq.0.d0) then
               xx  = -amasses(3)*amasses(3)/p3(0)/p3(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               oop3dk = p3(0)*q(0)*ombet/eta(3)
            endif
!!!
            
            oop3dk = 1.d0/oop3dk
         else
! trying to avoid NaNs at extremely high energies
            oop3dk = dot(p3,q) / eta(3)
            oop4dk = dot(p4,q) / eta(4)
            if (oop3dk.eq.0.d0) then
               xx  = -amasses(3)*amasses(3)/p3(0)/p3(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               oop3dk = p3(0)*q(0)*ombet/eta(3)
            endif
            if (oop4dk.eq.0.d0) then
               xx  = -amasses(4)*amasses(4)/p4(0)/p4(0)
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
               oop4dk = p4(0)*q(0)*ombet/eta(4)
            endif
!!!

            oop3dk = 1.d0/oop3dk
            oop4dk = 1.d0/oop4dk
         endif

         do i = 0,3
            peik(i)=p1(i)*oop1dk+p2(i)*oop2dk+p3(i)*oop3dk+p4(i)*oop4dk
         enddo         
         singlephoton = dot(peik,peik)
         sub = 0.d0
         if (eikng.eq.'multip') ksub = k
         scale = sqrtscalesub
         sub = subtractionnew(ksub,scale,amasses(1),amasses(2),
     -        2.d0/3.d0,1.d0/3.d0,q)
         sub = sub*iplsub
         singlephoton = -ec2 * singlephoton - sub
         eikonal = eikonal * singlephoton
         if (eikng.eq.'multip') singlephotons(k) = singlephoton
      enddo
      eikonal = eikonal/nfactorial(nglocal)

      return
      end
