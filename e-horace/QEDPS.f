**************** subroutines for QED radiation in Parton Shower
*                approximation
      subroutine qed_stuff(E1,E2,q2,ccm,phi,p3,p4,qph,x1,x2,x3,x4,ie)      
      include 'shared.inc'
      dimension x(40),xvect(10)
      ie = 0
      if (qedorder.eq.'exp') then
*      
         qminis1 = m1
         qminis2 = m2

         tmp = charge_factor(1)*ch1
         call dqed(q2,qminis1,tmp,xvect,x1)
         do k = 1, 10
            x(k) = xvect(k)
         enddo
*      
         tmp = charge_factor(2)*ch2
         call dqed(q2,qminis2,tmp,xvect,x2)  
         do k = 1, 10
            x(k+10) = xvect(k)
         enddo
*     
         tmp = charge_factor(3)*chfs
         if (boson.eq.'Z') tmp = -1.d0*charge_factor(3)

         call dqed(q2,mfs,tmp,xvect,x3)  
         do k = 1, 10
            x(k+20) = xvect(k)
         enddo
*     
         tmp = charge_factor(4)*0.d0
         if (boson.eq.'Z') tmp = 1.d0*charge_factor(4)

         call dqed(q2,mfs,tmp,xvect,x4)  
         do k = 1, 10
            x(k+30) = xvect(k)
         enddo      
*     
      else                      ! if (qedorder.eq.'alpha')
*
         if (boson.eq.'Z') then
            chfs1 = -1.d0
            chfs2 = 1.d0
         else
            chfs1 = chfs
            chfs2 = 0.d0
         endif
         call dqed_oalpha(q2,m1,m2,mfs,mfs,ch1,ch2,chfs1,
     >        chfs2,x,x1,x2,x3,x4)
      endif
*
      if (boson.eq.'Z') then
         am1 = mfs
         am2 = mfs
      else
         am1 = mfs
         am2 = 0.d0
      endif

      if (new_kine.eq.1) then 
         call kinematics(E1,E2,m1,m2,am1,am2,x1,x2,x3,x4,x,
     >                ccm,phi,p3,p4,qph,ie)
      else
         call kinematics_cfrWINHAC(E1,E2,m1,m2,am1,am2,x1,x2,x3,x4,x,
     >        ccm,phi,p3,p4,qph,ie)
      endif
*
      return
      end
*
*************************************************************************
* QED Parton Shower at order alpha.....
      subroutine dqed_oalpha(q2,ama,amb,amc,amd,car1,car2,car3,car4,x,
     >                       x1,x2,x3,x4)
      include 'shared.inc'
      dimension x(40)
      double precision mqua(4),qqua(4)
      real*4 rvec(2)
      
      do k=1,40
         x(k) = 0.d0
      enddo       
      
      x1 = 1.d0
      x2 = 1.d0
      x3 = 1.d0
      x4 = 1.d0

      qqua(1) = (charge_factor(1)*car1)**2
      qqua(2) = (charge_factor(2)*car2)**2
      qqua(3) = (charge_factor(3)*car3)**2
      qqua(4) = (charge_factor(4)*car4)**2
     
      qsum = qqua(1) + qqua(2) + qqua(3) + qqua(4)
      if (qsum.lt.1.d-1) return

      mqua(1) = ama**2*exp(1.d0)
      mqua(2) = amb**2*exp(1.d0)

      mqua(3) = amc**2*exp(1.d0)
!!!!!!!!!!!! SF SCALE !!!!
      if (boson.eq.'W') mqua(3) = amc**2*exp(2.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!
      mqua(4) = amd**2*exp(1.d0)

!! for gg - induced 
c      mqua(1) = mqua(3)
c      mqua(2) = mqua(4)
!!      
      pnb = ffsoa(q2,mqua,qqua)

      if (pnb.lt.0.d0) then
         print*,' WARNING: Sudakov Form Factor < 0 !',pnb
      endif
      
      call wraprng(rvec,2)      
      r1 = 1.d0*rvec(1)
      r2 = 1.d0*rvec(2)
      
      if (r1.lt.pnb) return
      
      leg   = 1
      pnorm = (1.d0-pnb)*2.d0*pi/alpha/iplus
      pirr  = qqua(leg)*log(q2/mqua(leg))/pnorm
      do while(r2.gt.pirr)
         leg  = leg  + 1
         pirr = pirr + qqua(leg)*log(q2/mqua(leg))/pnorm
      enddo
      
      call ap_vertex(y)
      y = 1.d0 - y 
      
      index    = 10*(leg-1)+1
      x(index) = y
      
      x1 = 1.d0 - x(1)
      x2 = 1.d0 - x(11)
      x3 = 1.d0 - x(21)
      x4 = 1.d0 - x(31)
      return
      end
*******************************************************      
      function ffsoa(q2,mqa,qqa)
      include 'shared.inc'
      double precision mqa(4),qqa(4)
      ffsoa = 0.d0
      do i = 1,4
         ffsoa = ffsoa + qqa(i)*log(q2/mqa(i))
      enddo      
      ffsoa = ffsoa*alpha/2.d0/pi*iplus      
      ffsoa = 1.d0 - ffsoa      
      return
      end
*************************************************************************
      subroutine dqed(q2,m,ch,xv,x)
*
* q2 is the SF Q^2 scale, m and ch are the mass and of the emitting particle,
* xv(10) are the photons energy fractions extracted by the PS and x is
* the energy fraction left to the emitting particle
*
      include 'shared.inc'
      common/n_phot/nphot
      real*4 r(1)
      double precision m
      dimension xv(10)
*      
      do k = 1, 10
         xv(k) = 0.d0
      enddo
      x     = 1.d0
      nphot = 0
*
      if (abs(ch).lt.1.d-1) return
*
      q2min = m**2*exp(1.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!! SF SCALE !!!
      if (boson.eq.'W') q2min = m**2*exp(2.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ak2   = q2min
*      
      irigenera = 1

      do while(irigenera.eq.1)
         call wraprng(r,1)
*
         rand = 1.d0*r(1)
         ffsn = ffs(q2,ak2,ch)
         if (rand.le.ffsn) return
         nphot = nphot + 1
*     
         akp2 = -2.d0*pi/(ch**2*alpha)/iplus*log(rand)
         akp2 = ak2*exp(akp2)                  
*     
         call ap_vertex(y)
*     
         if(nphot.le.10) then 
            xv(nphot) = 1.d0-y
            x         = x*y
         endif     
         ak2 = akp2
      enddo 
      return
      end
*
************************************************************************
      function ffs(q2,qp2,c)
* no-branching probability (Sudakov form factor) 
      include 'shared.inc'
      ffs = -c**2*alpha/2.d0/pi*log(q2/qp2)*iplus
      ffs = exp(ffs)
      return
      end
************************************************************************
      subroutine ap_vertex(x)
*  x generation according to ap splitting function  
*  (1+x^2)/(1-x), 0 <= x <= 1-eps
      include 'shared.inc'
      real*4 r(2)
*      
      alne = log(eps)
      irigenera = 1
      do while(irigenera.eq.1)
         call wraprng(r,2)
         cx = r(1)
         omx = exp(cx*alne)
         x = 1.d0 - omx
         rx = r(2)*2.d0/omx
         px = (1.d0 + x**2)/omx
         if (rx.lt.px) irigenera = 0
      enddo
      return
      end
***********************************************************************
      subroutine kinematics(E1,E2,ama,amb,amc,amd,xa,xb,xc,xd,x,
     >                      ccm,phi,p3,p4,qph,ie)
      include 'shared.inc'
      dimension x(40),vboost(0:3),xf(4),d3(0:3),d4(0:3)
      dimension ad3(0:3),ad4(0:3)
      double precision mf(4),mtotfs,tmp2(0:3),kone(0:3),ver(0:3)
      real*4 rvec(2),csi(1)
      
      ie = 0
      do k=0,3
         vboost(k) = 0.d0
      enddo      
      
      do k = 1,40
         qph(k,0) = 0.d0
         qph(k,1) = 0.d0
         qph(k,2) = 0.d0
         qph(k,3) = 0.d0
      enddo      
      
      xf(1) = xa
      xf(2) = xb
      xf(3) = xc
      xf(4) = xd
      
      mf(1) = ama
      mf(2) = amb
      mf(3) = amc
      mf(4) = amd

      scm  = sqrt(1.d0 - ccm*ccm)
      sphi = sin(phi)
      cphi = cos(phi)
      
      mtotfs = mf(3)+mf(4)
*
* the beta and gamma of the boost due to QED radiation are calculated
* 
      bqed = (xf(1)-xf(2))/(xf(1)+xf(2))
! I want for the moment only FS-radiation: check that bqed is always = 0
      if (bqed.gt.1.d-8) call messages(6)

      gqed = 1.d0/sqrt(1.d0-bqed*bqed)
*
* The momenta of the FS particles in the CM are reconstructed and boosted back
*
      E3 = xf(1)*E1*gqed*(1.d0-bqed)  ! E3 and E4 must be equal...
      E4 = xf(2)*E2*gqed*(1.d0+bqed)

! I try to mantain fully massive kinematics for final state....      
      Etot = E3 + E4      
      if (Etot.lt.mtotfs) then
         ie = 1
         return
      endif
      
      E3 = (Etot**2+mf(3)**2-mf(4)**2)/2.d0/Etot
      E4 = Etot-E3
      
      b3 = (1.d0+mf(3)/E3)*(1.d0-mf(3)/E3)
      b3 = sqrt(b3)
      
      b4 = (1.d0+mf(4)/E4)*(1.d0-mf(4)/E4)
      b4 = sqrt(b4)
            
      p3(1) = b3*E3*sphi*scm
      p3(2) = b3*E3*cphi*scm
      p3(3) = b3*E3*ccm
      p3(0) = E3

      p4(1) = -p3(1)   
      p4(2) = -p3(2)   
      p4(3) = -p3(3)   
      p4(0) =  E4   
      
      vboost(3) = -bqed            
      call boost(gqed,vboost,p3,p3)
      call boost(gqed,vboost,p4,p4)

      xprod = xf(1)*xf(2)*xf(3)*xf(4)

      if (dabs(xprod-1.d0).lt.1.d-10) return ! no radiation

      E3after = p3(0)*xf(3)
      E4after = p4(0)*xf(4)
      if ((E3after.lt.mf(3)).or.(E4after.lt.mf(4))) then
         ie = 2
         return
      endif

      b3 = (1.d0+mf(3)/E3after)*(1.d0-mf(3)/E3after)
      b3 = sqrt(b3)      
      b4 = (1.d0+mf(4)/E4after)*(1.d0-mf(4)/E4after)
      b4 = sqrt(b4)
      p3m = sqrt(tridot(p3,p3))
      p4m = sqrt(tridot(p4,p4))
      
      d3(0) = E3after
      d4(0) = E4after
      do k=1,3
         d3(k) = d3(0)*b3*p3(k)/p3m
         d4(k) = d4(0)*b4*p4(k)/p4m
      enddo

! IREF e' settato indip. dai fotoni!!
      call wraprng(csi,1)
      iref = 3
      ic4 = charge_factor(4)
      if (csi(1).lt.0.5d0) iref = (1-ic4)*3+ic4*4
*
* Once boosted back, p3 and p4 are degradated for QED radiation....
*      
      E3before = p3(0)
      E4before = p4(0)
      
********** PHOTONS are reconstructed!!! ***********
      o1 = E1
      o2 = E2
      o3 = E3before
      o4 = E4before
         
      do k = 1,10         
         qph(k,0)    = o1*x(k)
         qph(k+10,0) = o2*x(k+10)
         qph(k+20,0) = o3*x(k+20)
         qph(k+30,0) = o4*x(k+30)
         
         o1 = o1 - qph(k,0)
         o2 = o2 - qph(k+10,0)
         o3 = o3 - qph(k+20,0)
         o4 = o4 - qph(k+30,0)
      enddo

      numerodifotoni = 0
      do k=0,3
         kone(k)   = 0.d0
      enddo      
      
      do k = 1,40
         eg = qph(k,0)
         if (eg.gt.0.d0) then

            if (k.le.10) ileg = 1
            if (k.gt.10.and.k.le.20) ileg = 2
            if (k.gt.20.and.k.le.30) ileg = 3
            if (k.gt.30) ileg = 4

            numerodifotoni=numerodifotoni+1
            
            if (ptqed.eq.'no') then   
***   Collinear photons....            
!     qph(k,3) = qph(k,0)
               ver(0) = 1.d0
               if (ileg.eq.1) then
                  ver(1) = 0.d0
                  ver(2) = 0.d0
                  ver(3) = ver(0)
               endif
               if (ileg.eq.2) then
                  ver(1) = 0.d0
                  ver(2) = 0.d0
                  ver(3) = -ver(0)
               endif
               if (ileg.eq.3) then
                  d3mod = sqrt(tridot(d3,d3))
                  ver(1) = d3(1)/d3mod
                  ver(2) = d3(2)/d3mod
                  ver(3) = d3(3)/d3mod
               endif
               if (ileg.eq.4) then
                  d4mod = sqrt(tridot(d4,d4))
                  ver(1) = d4(1)/d4mod
                  ver(2) = d4(2)/d4mod
                  ver(3) = d4(3)/d4mod
               endif
***               
               do i = 0,3
                  qph(k,i) = eg*ver(i)
                  kone(i)  = kone(i) + qph(k,i)
               enddo
            else ! fotoni non collineari
c$$$  !!!!!!!!!!! COSTHG NEW CHOICE !!!!!!!!!!!!!
c$$$  call getcthg(cth,b3)
c$$$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               
!               if (ileg.eq.1) call alongp(pm,ver,a1,w1)
!               if (ileg.eq.2) call alongp(pp,ver,a2,w2)

               if (ileg.le.2) call messages(1)
               if (ileg.gt.2) call get_vers(d3,d4,ver)

               do i = 0,3
                  qph(k,i) = qph(k,0)*ver(i)
                  kone(i)  = kone(i) + qph(k,i)
               enddo
            endif
         endif
      enddo
      
*     I calculate p3 e p4
      
      ei = E1 + E2
      amkone2 = dot(kone,kone)
      
      qu2 = ei**2 + amkone2 - 2.d0*ei*kone(0) ! it is the mass^2 for p3 p4
      
      if (qu2.lt.(mf(3)+mf(4))**2) then
         ie = 1
         return
      endif
      
***   
      if (iref.eq.3) then
         do k=0,3
            ver(k)=d3(k)
         enddo
         am3 = mf(3)
         am4 = mf(4)
      endif      
      if (iref.eq.4) then
         do k=0,3
            ver(k)=d4(k)
         enddo
         am3 = mf(4)
         am4 = mf(3)
      endif      
      
      verm = sqrt(tridot(ver,ver))
      vx = ver(1)/verm
      vy = ver(2)/verm
      vz = ver(3)/verm
      
      o  = kone(0)
      ox = kone(1) 
      oy = kone(2)
      oz = kone(3)
      
      odv = ox*vx + oy*vy + oz*vz
      
      abig = ei**2 - 2.d0 * ei * o + amkone2 + am3**2 - am4**2
      bbig = -2.d0 * (ei - o)
      cbig = -2.d0 * odv
      
      arg = abig**2*bbig**2-am3**2*bbig**2*(bbig**2-cbig**2)
      
      if (arg.lt.0.d0) then
         ie = 1
         return
      endif
      
      p3mod_1 = abig*cbig + sqrt(arg)
      p3mod_1 = p3mod_1 / (bbig**2 - cbig**2)
      
      p3mod_2 = abig*cbig - sqrt(arg)
      p3mod_2 = p3mod_2 / (bbig**2 - cbig**2)
      
      p3jac = 1.d0
      if (p3mod_1.gt.0.d0.and.p3mod_2.gt.0.d0) then
         call wraprng(csi,1)
         if (csi(1).lt.0.5) p3mod = p3mod_1
         if (csi(1).gt.0.5) p3mod = p3mod_2
         p3jac = 2.d0
      endif
      if (p3mod_1.gt.0.d0.and.p3mod_2.lt.0.d0) then
         p3mod = p3mod_1
      endif
      if (p3mod_1.lt.0.d0.and.p3mod_2.gt.0.d0) then
         p3mod = p3mod_2
      endif
      if (p3mod_1.lt.0.d0.and.p3mod_2.lt.0.d0) then
         ie = 1
         return
      endif
      
      p3(0) = sqrt(p3mod**2 + am3**2)
      p3(1) = p3mod * vx
      p3(2) = p3mod * vy
      p3(3) = p3mod * vz
      
      p4(0) =  ei - p3(0) - o
      p4(1) = -ox - p3(1)
      p4(2) = -oy - p3(2)
      p4(3) = -oz - p3(3)
      
      if (iref.eq.4) call exchange_mom(p3,p4)
***   
      return
      end
***************************************************
      subroutine get_vers(d3,d4,ver)
      implicit double precision (a-h,o-z)
      dimension d3(0:3),d4(0:3),ver(0:3)
      real*4 csi(1)
      character*1 boson
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/vectorboson/boson

      if (boson.eq.'W') then
         b = sqrt(tridot(d3,d3))/d3(0)
         call wraprng(csi,1)
         phi = 2.d0 * pi * csi(1)
         cp = cos(phi)
         sp = sin(phi)
         call getcthg_W(cth,b)
         sth = sqrt(1.d0-cth*cth)
         ver(0) = 1.d0
         ver(1) = sp * sth
         ver(2) = cp * sth
         ver(3) = cth
         call rot(-1,d3,ver,ver)
      else ! if boson = Z
         istop = 0
         icount = 0
         DO WHILE(istop.eq.0)
            icount = icount + 1
            call wraprng(csi,1)
            if (csi(1).lt.0.5d0) then 
               call alongp(d3,ver,a3,w3)
            else
               call alongp(d4,ver,a4,w4)
            endif
            call alongp_nandw(d3,ver,a3,w3)
            call alongp_nandw(d4,ver,a4,w4)
            call PSangular_fact(d3,d4,ver,aok,wok)
            wwrong = (w3+w4) !/2.d0            
            call wraprng(csi,1)
            tetto = wwrong * 2.1d0
            cfr = tetto * csi(1)           
!           if (tetto.lt.wok) print*,wok/tetto,' <-- tetto basso!'
            if (cfr.lt.wok) istop = 1
         ENDDO
      endif
      return
      end
********************************************************
      subroutine getcthg(c,b)
      include 'shared.inc'
! this subroutine is needed by the old kinematics routine
! (kinematics_cfrWINHAC)....
      if (boson.eq.'W') then
         call getcthg_W(c,b)
      else
         call dummy_subroutine
      endif
      return
      end
*****
      subroutine getcthg_W(c,b)
      implicit double precision (a-h,o-z)
      double precision m2oE2
      real*4 r(2)

      m2oE2 = 1.d0 - b**2
      istop = 0
      icount = 0

      do while(istop.eq.0)
         icount = icount + 1

         call wraprng(r,2)
         an = -1.d0/b*log((1.d0-b)/(1.d0+b))
         c  = (1.d0 - (1.d0+b)*exp(-b*an*r(1)))/b
         
         umbc = 1.d0 - b*c

         w_sampling = 1.d0 / umbc
         r_sampling = w_sampling - 0.5d0 - 0.5d0*m2oE2/umbc/umbc
         
         ratio = r_sampling / w_sampling
        
         csi = 1.d0*r(2)

         if (ratio.gt.1.d0.or.ratio.lt.0.d0) then 
            print*,'WARNING!! in subroutine getcthg_W, ',ratio
         endif
         if (csi.lt.ratio) istop = 1
      enddo
      return
      end
*****
      subroutine dummy_subroutine
      implicit double precision (a-h,o-z)
      ! dummy_subroutine
      call messages(2)
      return
      end
********************************+
      subroutine alongp(p,vers,anorm,w)
      implicit double precision (a-h,o-z)
      dimension p(0:3),vers(0:3),ptmp(0:3)
      real*4 rv(1)
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac

      pmod = sqrt(tridot(p,p))
      b1 = pmod/p(0)
      vers(0) = 1.d0

      call wraprng(rv,1)
      phph = 2.d0*pi*rv(1)
      call wraprng(rv,1)

      an  = -1.d0/b1*log((1.d0-b1)/(1.d0+b1))
      cth = (1.d0 - (1.d0+b1)*
     #     exp(-b1*an*rv(1)))/b1
      
      sth = sqrt(1.d0-cth*cth)
      cph = cos(phph)
      sph = sin(phph)
      
      vers(1) = vers(0)*sph*sth
      vers(2) = vers(0)*cph*sth
      vers(3) = vers(0)*cth
                  
      call rot(-1,p,vers,vers)
      
      anorm = an
      w     = 1.d0
!!!      call alongp_nandw(p,vers,anorm,w) ! e' inutile
      return
      end

      subroutine alongp_nandw(p,vers,anorm,w)
      implicit double precision (a-h,o-z)
      dimension p(0:3),vers(0:3)

      pmod = sqrt(tridot(p,p))
      b1 = pmod/p(0)

      anorm  = -1.d0/b1*log((1.d0-b1)/(1.d0+b1))

      w = p(0)*vers(0)/dot(p,vers)/anorm
! w e' normalizzato a 1!

      return
      end
************
      subroutine PSangular_fact(p3,p4,q,anorm,w)
      implicit double precision (A-H,O-Z)
      dimension p3(0:3),p4(0:3),ptmp(0:3)
      dimension q(0:3),p(4,0:3),tmp1(0:3),tmp2(0:3)
      dimension ieta(4)

! w is the normalized function!!

      ieta(1) =  1
      ieta(2) = -1
      ieta(3) = -1
      ieta(4) =  1
      
      do k = 0,3
         p(3,k) = p3(k)
         p(4,k) = p4(k)
      enddo

      sum = 0.d0
      anorm = 0.d0
      ame2 = dot(p3,p3)
!      ame2 = dot(p4,p4)
      ame  = sqrt(ame2)

      do k=3,4
         do j=3,4
            do i = 0,3
               tmp1(i) = p(j,i)
               tmp2(i) = p(k,i)
            enddo
            etaietaj = -ieta(k)*ieta(j)
            term  = etaietaj
            p1dp2 = dot(tmp1,tmp2)
            term  = term * p1dp2
            term  = term/dot(tmp1,q)/dot(tmp2,q)*q(0)**2
            sum = sum + term
            if (k.ne.j) then
               f = f_for_coherent_norm(p1dp2,ame)
               anorm=anorm + etaietaj*f
            endif
         enddo
      enddo
      w = sum 
      anorm = (anorm - 2.d0)
      w = w/anorm
      return
      end
********************************************************
      function f_for_coherent_norm(pidpj,am)
      implicit double precision (a-h,o-z)
* formula 2.21 tesi di dottorato * k_0^2 / (4*pi)
      x = pidpj
      am2 = am**2
      radice = sqrt(x**2-am2**2)
      arg = (radice + x - am2)/(radice - (x-am2))
      f = x / radice * log(arg)
      f_for_coherent_norm = f
      return
      end

***************************************************************

      subroutine kinematics_cfrWINHAC(E1,E2,ama,amb,amc,amd,xa,xb,
     >     xc,xd,x,ccm,phi,p3,p4,qph,ie)
      
      include 'shared.inc'
      dimension x(40),vboost(0:3),xf(4)
      double precision mf(4),mtotfs,tmp2(0:3),kone(0:3),ver(0:3)
      real*4 rvec(2),csi(1)
      
      ie = 0
      do k=0,3
         vboost(k) = 0.d0
         kone(k)   = 0.d0
      enddo      
      
      do k = 1,40
         qph(k,0) = 0.d0
         qph(k,1) = 0.d0
         qph(k,2) = 0.d0
         qph(k,3) = 0.d0
      enddo      
      
      xf(1) = xa
      xf(2) = xb
      xf(3) = xc
      xf(4) = xd
      
      mf(1) = ama
      mf(2) = amb
      mf(3) = amc
      mf(4) = amd

      scm  = sqrt(1.d0 - ccm*ccm)
      sphi = sin(phi)
      cphi = cos(phi)
      
      mtotfs = mf(3)+mf(4)
*
* the beta and gamma of the boost due to QED radiation are calculated
* 
      bqed = (xf(1)-xf(2))/(xf(1)+xf(2))
      gqed = 1.d0/sqrt(1.d0-bqed*bqed)
*
* The momenta of the FS particles in the CM are reconstructed and boosted back
*
      E3 = xf(1)*E1*gqed*(1.d0-bqed)  ! E3 and E4 must be equal...
      E4 = xf(2)*E2*gqed*(1.d0+bqed)

! I try to mantain fully massive kinematics for final state....
      
      Etot = E3 + E4
      
      if (Etot.lt.mtotfs) then
         ie = 1
         return
      endif
      
      E3 = (Etot**2+mf(3)**2-mf(4)**2)/2.d0/Etot
      E4 = Etot-E3
      
      b3 = (1.d0+mf(3)/E3)*(1.d0-mf(3)/E3)
      b3 = sqrt(b3)
      
      b4 = (1.d0+mf(4)/E4)*(1.d0-mf(4)/E4)
      b4 = sqrt(b4)
            
      p3(1) = b3*E3*sphi*scm
      p3(2) = b3*E3*cphi*scm
      p3(3) = b3*E3*ccm
      p3(0) = E3

      p4(1) = -p3(1)   
      p4(2) = -p3(2)   
      p4(3) = -p3(3)   
      p4(0) =  E4   
      
      vboost(3) = -bqed
            
      call boost(gqed,vboost,p3,p3)
      call boost(gqed,vboost,p4,p4)
*
* Once boosted back, p3 and p4 are degradated for QED radiation....
*      
      p3modmo = 1.d0/sqrt(tridot(p3,p3))
      p4modmo = 1.d0/sqrt(tridot(p4,p4))
      
      E3before = p3(0)
      E4before = p4(0)
      
      p3(0) = p3(0)*xf(3)
      p4(0) = p4(0)*xf(4)

      if ((p3(0).lt.mf(3)).or.(p4(0).lt.mf(4))) then
         ie = 2
         return
      endif
       
      b3 = (1.d0+mf(3)/p3(0))*(1.d0-mf(3)/p3(0))
      b3 = sqrt(b3)
      b4 = (1.d0+mf(4)/p4(0))*(1.d0-mf(4)/p4(0))
      b4 = sqrt(b4)

      do i = 1,3
         p3(i)=p3(i)*p3modmo*b3*p3(0)
         p4(i)=p4(i)*p4modmo*b4*p4(0)
      enddo

********** PHOTONS are reconstructed!!! ***********
      o1 = E1
      o2 = E2
      o3 = E3before
      o4 = E4before
      do k = 1,10
         qph(k,0)    = o1*x(k)
         qph(k+10,0) = o2*x(k+10)
         qph(k+20,0) = o3*x(k+20)
         qph(k+30,0) = o4*x(k+30)

         o1 = o1 - qph(k,0)
         o2 = o2 - qph(k+10,0)
         o3 = o3 - qph(k+20,0)
         o4 = o4 - qph(k+30,0)
     
         if (ptqed.eq.'no') then   
*** Collinear photons....
            qph(k,3)    = qph(k,0)
            qph(k+10,3) = qph(k+10,0)
            qph(k+20,3) = qph(k+20,0)
            qph(k+30,3) = qph(k+30,0)
         else
***   Photons with p_t !      
            b1 = (1+mf(1)/E1)*(1-mf(1)/E1)
            if (b1.gt.0.d0) b1 = sqrt(b1)
            b2 = (1+mf(2)/E2)*(1-mf(2)/E2)
            if (b2.gt.0.d0) b2 = sqrt(b2)
            b3 = (1+mf(3)/E3before)*(1-mf(3)/E3before)
            b3 = sqrt(b3)
            b4 = (1+mf(4)/E4before)*(1-mf(4)/E4before)
            b4 = sqrt(b4)
**       
            if (qph(k,0).gt.0.d0) then
               call wraprng(rvec,2)
               phph = 2.d0*pi*rvec(1)
               
               an = -1.d0/b1*log((1.d0-b1)/(1.d0+b1))
               cth = (1.d0 - (1.d0+b1)*
     #              exp(-b1*an*rvec(2)))/b1
               
               sth = sqrt(1.d0-cth*cth)
               cph = cos(phph)
               sph = sin(phph)
               
               qph(k,1) = qph(k,0)*sph*sth
               qph(k,2) = qph(k,0)*cph*sth
               qph(k,3) = qph(k,0)*cth

            endif
**       
            if (qph(k+10,0).gt.0.d0) then 
               call wraprng(rvec,2)
               phph = 2.d0*pi*rvec(1)
               
               an = -1.d0/b2*log((1.d0-b2)/(1.d0+b2))
               cth = (1.d0 - (1.d0+b2)*
     #              exp(-b2*an*rvec(2)))/b2
               
               sth = sqrt(1.d0-cth*cth)
               cph = cos(phph)
               sph = sin(phph)
               
               qph(k+10,1) = qph(k+10,0)*sph*sth
               qph(k+10,2) = qph(k+10,0)*cph*sth
               qph(k+10,3) = qph(k+10,0)*cth

            endif
**
            if (qph(k+20,0).gt.0.d0) then                
               call wraprng(rvec,2)
               phph = 2.d0*pi*rvec(1)
               
               an = -1.d0/b3*log((1.d0-b3)/(1.d0+b3))
               cth = (1.d0 - (1.d0+b3)*
     #              exp(-b3*an*rvec(2)))/b3

!!!!!!!!!!! COSTHG NEW CHOICE !!!!!!!!!!!!!
               call getcthg(cth,b3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               
               sth = sqrt(1.d0-cth*cth)
               cph = cos(phph)
               sph = sin(phph)
               
               qph(k+20,1) = qph(k+20,0)*sph*sth
               qph(k+20,2) = qph(k+20,0)*cph*sth
               qph(k+20,3) = qph(k+20,0)*cth
               
            endif
**
            if (qph(k+30,0).gt.0.d0) then 
               
               call wraprng(rvec,2)
               phph = 2.d0*pi*rvec(1)
               
               an = -1.d0/b4*log((1.d0-b4)/(1.d0+b4))
               cth = (1.d0 - (1.d0+b4)*
     #              exp(-b4*an*rvec(2)))/b4
               
               sth = sqrt(1.d0-cth*cth)
               cph = cos(phph)
               sph = sin(phph)
               
               qph(k+30,1) = qph(k+30,0)*sph*sth
               qph(k+30,2) = qph(k+30,0)*cph*sth
               qph(k+30,3) = qph(k+30,0)*cth
       
            endif
**
         endif      
      enddo

**
* reporting photons in the lab system....

***************************************************
      tmp2(0) =  1.d0
      tmp2(1) =  0.d0
      tmp2(2) =  0.d0
      tmp2(3) = -1.d0

      do k = 1,10
**       
         if (qph(k+10,0).gt.0.d0) then
            
            ptmp(0) = qph(k+10,0)
            ptmp(1) = qph(k+10,1)
            ptmp(2) = qph(k+10,2)
            ptmp(3) = qph(k+10,3)
            
            call rot(-1,tmp2,ptmp,ptmp)
            
            qph(k+10,0) = ptmp(0) 
            qph(k+10,1) = ptmp(1) 
            qph(k+10,2) = ptmp(2)
            qph(k+10,3) = ptmp(3)
            
         endif
**       
         if (qph(k+20,0).gt.0.d0) then
            
            ptmp(0) = qph(k+20,0)
            ptmp(1) = qph(k+20,1)
            ptmp(2) = qph(k+20,2)
            ptmp(3) = qph(k+20,3)
            
            call rot(-1,p3,ptmp,ptmp)
            
            qph(k+20,0) = ptmp(0) 
            qph(k+20,1) = ptmp(1) 
            qph(k+20,2) = ptmp(2)
            qph(k+20,3) = ptmp(3)
       
         endif
**      
         if (qph(k+30,0).gt.0.d0) then

            ptmp(0) = qph(k+30,0)
            ptmp(1) = qph(k+30,1)
            ptmp(2) = qph(k+30,2)
            ptmp(3) = qph(k+30,3)
            
            call rot(-1,p4,ptmp,ptmp)
            
            qph(k+30,0) = ptmp(0) 
            qph(k+30,1) = ptmp(1) 
            qph(k+30,2) = ptmp(2)
            qph(k+30,3) = ptmp(3)

         endif
**
      enddo
*
* I recalculate p3 and p4 in order to have exact 4-momentum conservation...
* The strategy is analogous to the one used in BABAYAGA
*
      do k=1,40
         do j=0,3
            kone(j)=kone(j)+qph(k,j)
         enddo
      enddo        
      
      if (kone(0).lt.1.d-10) return !!!!!!!!!!!!!
      
      ptmp(0) = E1+E2-kone(0)
      ptmp(1) = -kone(1)
      ptmp(2) = -kone(2)
      ptmp(3) = E1-E2-kone(3)
      
      if (dot(ptmp,ptmp).lt.mtotfs**2) then
         ie = 3 
         return
      endif 
      
      vboost(0) = 0.d0
      vboost(1) = ptmp(1)/ptmp(0)
      vboost(2) = ptmp(2)/ptmp(0)
      vboost(3) = ptmp(3)/ptmp(0)
      
      bcm=sqrt(vboost(1)**2+vboost(2)**2+vboost(3)**2)
      gcm=1.d0/sqrt(1.d0-bcm**2)
 
      call boost(gcm,vboost,ptmp,ptmp)

      Etot = ptmp(0)

      E3 = (Etot**2+mf(3)**2-mf(4)**2)/2.d0/Etot
      E4 = Etot-E3
      
      b3 = (1.d0+mf(3)/E3)*(1.d0-mf(3)/E3)
      b3 = sqrt(b3)
      
      b4 = (1.d0+mf(4)/E4)*(1.d0-mf(4)/E4)
      b4 = sqrt(b4)
*            
      call wraprng(csi,1)

      if (csi(1).lt.0.5) then
         
         p3modmo = 1.d0/sqrt(tridot(p3,p3))
     
         ver(0) = 1.d0          ! I forgot masses....
         ver(1) = p3(1)*p3modmo
         ver(2) = p3(2)*p3modmo
         ver(3) = p3(3)*p3modmo

      else

         p4modmo = 1.d0/sqrt(tridot(p4,p4))
         ver(0) = 1.d0          ! I forgot masses....
         ver(1) = p4(1)*p4modmo
         ver(2) = p4(2)*p4modmo
         ver(3) = p4(3)*p4modmo

      endif

      call boost(gcm,vboost,ver,ver)

      if (csi(1).lt.0.5) then
 
         vermo = 1.d0/sqrt(tridot(ver,ver))

         p3(1) = b3*E3*ver(1)*vermo
         p3(2) = b3*E3*ver(2)*vermo
         p3(3) = b3*E3*ver(3)*vermo
         p3(0) = E3

         p4(1) = -p3(1)   
         p4(2) = -p3(2)   
         p4(3) = -p3(3)   
         p4(0) =  E4

      else

         vermo = 1.d0/sqrt(tridot(ver,ver))

         p4(1) = b4*E4*ver(1)*vermo
         p4(2) = b4*E4*ver(2)*vermo
         p4(3) = b4*E4*ver(3)*vermo
         p4(0) = E4

         p3(1) = -p4(1)   
         p3(2) = -p4(2)   
         p3(3) = -p4(3)   
         p3(0) =  E3      
      endif
      do i=0,3
       vboost(i)=-vboost(i)
      enddo  
      call boost(gcm,vboost,p3,p3)
      call boost(gcm,vboost,p4,p4)
      return
      end
