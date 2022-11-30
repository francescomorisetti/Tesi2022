      subroutine maptozerov2(ng,p1,p2,p3,p4,qph,iclose,masses,
     .     p1r,p2r,p3r,p4r,ierror)
      implicit double precision (a-h,o-z)
      parameter (imaxph=40)
      dimension iclose(imaxph),qph(imaxph,0:3)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension p1r(0:3),p2r(0:3),p3r(0:3),p4r(0:3)
      dimension q(0:3),p(0:3)
      double precision masses(4)
      common/counter_map0b2/icountmap0v2
      data icountmap0v2 /0/
      common/energiesforreduction/ei1,ei2,ei3,ei4,ef1,ef2,ef3,ef4

      ierror = 0
      icountmap0v2 = icountmap0v2 + 1
      ei1 = p1(0)
      ei2 = p2(0)
      ef1 = ei1
      ef2 = ei2
      ef3 = p3(0)
      ef4 = p4(0)
      ei3 = ef3
      ei4 = ef4
      do k = 1,ng
         if (iclose(k).eq.1) ef1 = ef1 - qph(k,0)
         if (iclose(k).eq.2) ef2 = ef2 - qph(k,0)
         if (iclose(k).eq.3) ei3 = ei3 + qph(k,0)
         if (iclose(k).eq.4) ei4 = ei4 + qph(k,0)
      enddo

      if (ef1.lt.masses(1)) then 
         ef1 = masses(1)
         ierror = 1
         return
      endif
      if (ef2.lt.masses(2)) then 
         ef2 = masses(2)
         ierror = 1
         return
      endif
      p1r(0) = ef1
      p1r(1) = 0.d0
      p1r(2) = 0.d0
      p1r(3) = ef1      

      p2r(0) =  ef2
      p2r(1) =  0.d0
      p2r(2) =  0.d0
      p2r(3) = -ef2      

      do k = 0,3
         p(k) = p1r(k)+p2r(k)
      enddo
      eicm = sqrt(dot(p,p)) ! invariant
      ei   = eicm

      ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
      eb = ei - ea
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

      do k=0,3
         q(k) = p1r(k) + p2r(k)
         p(k) = p3(k) + p4(k)
      enddo
      call new_boost(p,p3,p3r,1)
      call new_boost(p,p4,p4r,1)

      sshat = sqrt(dot(q,q)) ! invariant

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
      subroutine maptoonev2(ng,nk,p1b,p2b,p3b,p4b,q,iclose,masses,
     .     p1r,p2r,p3r,p4r,qr,ierror)
      implicit double precision (a-h,o-z)
      parameter (imaxph=40)
      dimension iclose(imaxph)
      dimension p1b(0:3),p2b(0:3),p3b(0:3),p4b(0:3)
      dimension p1r(0:3),p2r(0:3),p3r(0:3),p4r(0:3)
      dimension q(0:3),p(0:3),qr(0:3),ptmp(0:3)
      double precision masses(4)
      common/energiesforreduction/ei1,ei2,ei3,ei4,ef1,ef2,ef3,ef4
      common/counter_map1b2/icountmap1v2
      data icountmap1v2 /0/

      ierror = 0

      icountmap1v2 = icountmap1v2 + 1

      ef1l = ef1
      ef2l = ef2      
      ef3l = ef3
      ef4l = ef4      

      is = 0
      if (iclose(nk).lt.3) is = 1
      
      if (is.eq.1) then
         if (iclose(nk).eq.1) ef1l = ef1 + q(0)
         if (iclose(nk).eq.2) ef2l = ef2 + q(0)        
         p1r(0) = ef1l
         p1r(1) = 0.d0
         p1r(2) = 0.d0
         p1r(3) = ef1l         
         p2r(0) =  ef2l
         p2r(1) =  0.d0
         p2r(2) =  0.d0
         p2r(3) = -ef2l              
         do k = 0,3
            p(k) = p1r(k)+p2r(k)
         enddo         
         call new_boost(p,p1r,p1r,1)
         call new_boost(p,p2r,p2r,1)
         call new_boost(p,q,qr,1)                  

         sshat = sqrt(dot(p,p)) ! invariant
         e = sshat/2.d0
         ei   = 2.d0*e
         
         ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
         eb = ei - ea
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
            p(k)    = p3b(k)+p4b(k)
            ptmp(k) = p1r(k)+p2r(k)-qr(k)
         enddo

         if (dot(ptmp,ptmp).lt.0.d0) then
            ierror = 1
            return
         endif

         ef  = sqrt(dot(p,p)) 
         efn = sqrt(dot(ptmp,ptmp)) 
                  
         ei = efn
         abig = ei**2 + masses(3)**2 - masses(4)**2
         bbig = -2.d0 * ei 
         cbig = 0.d0
         arg  = abig**2*bbig**2-masses(3)**2*bbig**2*(bbig**2-cbig**2)
         p1mod = abig*cbig + sqrt(arg)
         p1mod = p1mod / (bbig**2 - cbig**2)
         e1    = sqrt(masses(3)**2 + p1mod**2)
         beta3 = p1mod/e1  
         p3mod = sqrt(tridot(p3b,p3b))
         p4mod = sqrt(tridot(p4b,p4b))
         p3r(0) = e1
         p3r(1) = e1*beta3 * p3b(1)/p3mod
         p3r(2) = e1*beta3 * p3b(2)/p3mod
         p3r(3) = e1*beta3 * p3b(3)/p3mod
         p4r(0) = ei - p3r(0)
         p4r(1) = -p3r(1)
         p4r(2) = -p3r(2)
         p4r(3) = -p3r(3)         
         call new_boost(ptmp,p3r,p3r,-1)
         call new_boost(ptmp,p4r,p4r,-1)         
      else
         p1r(0) = ef1l
         p1r(1) = 0.d0
         p1r(2) = 0.d0
         p1r(3) = ef1l         
         p2r(0) =  ef2l
         p2r(1) =  0.d0
         p2r(2) =  0.d0
         p2r(3) = -ef2l              
         do k = 0,3
            p(k) = p1r(k)+p2r(k)
         enddo         
         call new_boost(p,p1r,p1r,1)
         call new_boost(p,p2r,p2r,1)
         call new_boost(p,q,qr,1)                  

         sshat = sqrt(dot(p,p)) ! invariant
         e = sshat/2.d0
         ei   = 2.d0*e
         
         ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
         eb = ei - ea
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
            p(k)    = p3b(k)+p4b(k)
            ptmp(k) = p1r(k)+p2r(k)-qr(k)
         enddo
         if (dot(ptmp,ptmp).lt.0.d0) then
            ierror = 1
            return
         endif


         ef  = sqrt(dot(p,p)) 
         efn = sqrt(dot(ptmp,ptmp)) 
                  
         ei = efn
         abig = ei**2 + masses(3)**2 - masses(4)**2
         bbig = -2.d0 * ei 
         cbig = 0.d0
         arg  = abig**2*bbig**2-masses(3)**2*bbig**2*(bbig**2-cbig**2)
         p1mod = abig*cbig + sqrt(arg)
         p1mod = p1mod / (bbig**2 - cbig**2)
         e1    = sqrt(masses(3)**2 + p1mod**2)
         beta3 = p1mod/e1  
         p3mod = sqrt(tridot(p3b,p3b))
         p4mod = sqrt(tridot(p4b,p4b))
         p3r(0) = e1
         p3r(1) = e1*beta3 * p3b(1)/p3mod
         p3r(2) = e1*beta3 * p3b(2)/p3mod
         p3r(3) = e1*beta3 * p3b(3)/p3mod
         p4r(0) = ei - p3r(0)
         p4r(1) = -p3r(1)
         p4r(2) = -p3r(2)
         p4r(3) = -p3r(3)
         
         call new_boost(ptmp,p3r,p3r,-1)
         call new_boost(ptmp,p4r,p4r,-1)
      endif
      return
      end

************
************
************
      subroutine sampleomega(p2,e,bp,cpg,wm,wg,omin,omax,om,w)
      implicit double precision (a-h,o-z)
      common/fordebugging/idebugging
      common/peakom/ompeak
      common/parssampleomega/dir,wir,dbw,wbw,ddef,dp2,wp2
      common/icountso/icount
      data icount /0/
      icount = icount + 1

      w = 1.d0
      dbw = wg/10.d0
      r   = 9d0/10.d0
      dap = omin/r
c      dap = 0.1d0 * omin
*
      delta = min(dbw,dap)
      delta = 0.01d0
      range = omax-omin
      n = range/delta
*
      delta = range/(1.d0*n)
      ecap = 2.d0*e*(1.d0-bp*cpg)
      ompeak = max(omin,(p2 - wm**2)/ecap)

      omp1 = (p2-wg*wm - wm**2)/ecap
      omp2 = (p2+wg*wm - wm**2)/ecap

      dir  = 0.01d0
      wir  = 5.d0
      dbw  = (omp2-omp1)/20.d0
c      dbw  = 0.01d0
      wbw  = 8.d0
      dp2 = dbw/2.d0
      wp2 = wbw*2.d0

      ddef = 1.5d0
      n = nint((omax-omin)/ddef)+nint(wir/dir)+nint(wbw/dbw)
      n = n + nint(wp2/dp2)
      n = n + 3
      if (n.ge.250000) then
         print*,'WARNING: very big n ',n,omin,omax,delta,ompeak
      endif
      call nowsample(n,delta,ecap,p2,wm,wg,omax,omin,om,w)
      return
      end
c
      subroutine nowsample(n,delta,ecap,p2,wm,wg,omax,omin,om,w)
      implicit double precision (a-h,o-z)
      real*4 csi(2)
      dimension x(n+1),f(n+1),cumul(n+1),fint(n+1)
      common/fordebugging/idebugging
      common/icountns/icount
      common/parssampleomega/dir,wir,dbw,wbw,ddef,dp2,wp2
      common/cosgp3/cgp3
      common/peakom/ompeak
      data icount /0/
      icount = icount + 1
      w = 1.d0
c$$$      goto 777
c$$$      x(1) = omin
c$$$      q2   = p2 - ecap*x(1)
c$$$      q21  = q2
c$$$      xf   = 1.d0
c$$$      x1oomax   = x(1)/omax
c$$$      f(1) = bw(q2,wm,wg)*apv(x(1)/omax)
c$$$     .      +xf*bw(q21,wm,wg)*apv(x(1)/omax)
c$$$      sum  = f(1)
c$$$      do k = 2,n
c$$$         x(k) = x(k-1)+delta
c$$$         q2   = p2 - ecap*x(k)
c$$$         f(k) = bw(q2,wm,wg)*apv(x(k)/omax)
c$$$     .         +xf*bw(q21,wm,wg)*apv(x(k)/omax)
c$$$         sum  = sum + f(k)
c$$$      enddo
c$$$      x(n+1) = x(n)+delta
c$$$*----------------------------------------------------------
c$$$ 777  continue
      x(1)  = omin
      q21   = p2 - ecap*x(1)
      bwq1  = bw(q21,wm,wg)
      bwpeak= bw(p2-ecap*ompeak,wm,wg)
      fpeak = bwpeak*apv(ompeak/omax)*(p2-ecap*ompeak)
      xf    = 1.d0
      sum   = 0.d0
      istop = 0
      k     = 0
      en    = sqrt(p2)/2.d0
      omp2  = -2.d0*en* cgp3 / (1.d0 - cgp3)

      do while(istop.eq.0)
         k = k + 1
         dcur = ddef
         if (abs(x(k)-ompeak).lt.wbw/2.d0) then 
            dcur = dbw
         endif
         if (abs(x(k)-omp2).lt.wp2/2.d0) then 
            dcur = dp2
         endif
         if (x(k).lt.wir) dcur = dir
         if (k.gt.n) 
     .        print*,'warning in nowsample',n,k,omax,x(k),dcur,delta
         x(k+1) = x(k)+dcur
         q2   = p2 - ecap*x(k)
         bwq  = bw(q2,wm,wg)
         f(k) = bwq*apv(x(k)/omax)*q2
     .         +xf*bwq1*apv(x(k)/omax)*q21
         f(k) = f(k) * (en-x(k))/(2.d0*en - x(k)*(1.d0-cgp3))**2
         sum  = sum + f(k)
         if (x(k+1).ge.(omax-1.d-6)) then 
            x(k+1) = omax
            istop = 1
         endif
      enddo
      f(k+1) = f(k)/2.d0
      sum = 0.d0
      do i = 1,k
         dx = x(i+1)-x(i)
         fint(i) = f(i)*dx + (f(i+1)-f(i))*0.5d0*dx
         sum = sum + fint(i)
      enddo
*----------------------------------------------------------
      call wraprng(csi,2)
      summ1 = 1.d0/sum
c      cumul(1) = f(1)*summ1
      cumul(1) = fint(1)*summ1
      k = 1
      do while(csi(1).gt.cumul(k))
c         cumul(k+1) = cumul(k) + f(k+1)*summ1
         cumul(k+1) = cumul(k) + fint(k+1)*summ1
         k = k + 1
      enddo
c      om = x(k) + csi(2)*(x(k+1)-x(k))
c      w = sum/f(k)*(x(k+1)-x(k))
c      return
c      w = sum/fint(k)*(x(k+1)-x(k))
*---
      dx  = (x(k+1)-x(k))
      df  = (f(k+1)-f(k))
      der = df/dx
      a   = 0.5d0*der
      b   = f(k)-der*x(k)
      c   = -fint(k)*csi(2) - f(k)*x(k)+a*x(k)**2
      discr = b**2-4.d0*a*c
      om1 = -b + sqrt(discr)
      om1 = om1/2.d0/a
      om2 = -b - sqrt(discr)
      om2 = om2/2.d0/a
      if (om1.ge.x(k).and.om1.le.x(k+1)) om = om1
      if (om2.ge.x(k).and.om2.le.x(k+1)) om = om2
      fx = f(k) + der*(om-x(k))
      fx = fx/fint(k)
      w = sum/fint(k)/fx
c      idebugging = 0
c      if (abs(om-95.6316407d0).lt.0.001d0) then
c         idebugging = 1
c         q2=p2 - ecap*om
c         print*,'sampling ',om,om1,om2,q2
c         print*,'sampling prop',bw(q2,wm,wg)
c      endif
*---
      return
      end
      function bw(q2,wm,wg)
      implicit double precision (a-h,o-z)
      bw = 1.d0/((q2-wm**2)**2+wm**2*wg**2)
      return
      end
      function apv(y)
      implicit double precision (a-h,o-z)
      apv = 1.d0 + (1.d0-y)**2
      apv = apv/y
      return
      end
************
************
************
      subroutine photon_energy(pbw,rnd,pp,vers,q2min,
     .     omin,omax,mw,gw,om,w)
      implicit double precision (a-h,m,o-z)
      real*4 rnd,csi(1)
      dimension pp(0:3),vers(0:3)
!from shared.inc
      common/gen_event/ebeam,ebeam1,ebeam2,xmin,xmax,shat,x1pdf,x2pdf
      common/yminimo/ymincut ! for x1,x2pdf sampling
      w = 1.d0
      bppc1pp = tridot(pp,vers)/pp(0)
      umbc1p  = 1.d0 - bppc1pp
      pp2   = dot(pp,pp)
      q2max = pp2 - 2.d0*omin*umbc1p*pp(0)
      ymax  = q2max/pp2
      ymin  = q2min/pp2
      ymin  = max(ymin,0.d0)
      yjac  =  pp2/2.d0/pp(0)/umbc1p

      gwl  = gw
      gwl2 = 2.d0*gwl

      anbw  = 1.d0/(pp2*mw*gwl)*( atan( ymax*pp2/mw/gwl-mw/gwl )
     >     - atan( ymin*pp2/mw/gwl-mw/gwl ) )
      anbw2 = 1.d0/(pp2*mw*gwl2)*( atan( ymax*pp2/mw/gwl2-mw/gwl2 )
     >     - atan( ymin*pp2/mw/gwl2-mw/gwl2 ) )

      eps  = omin/omax
      anir = -2.d0*log(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
c      anir = log(omax/omin)
ccc TEST
c$$$      goto 555
c$$$      if (pp2.gt.0.d0) then
c$$$         bp  = sqrt(tridot(pp,pp))/pp(0)
c$$$         e   = sqrt(pp2)
c$$$         cpg = 0.d0
c$$$         if (bp.gt.0.d0) cpg = tridot(pp,vers)/sqrt(tridot(pp,pp))
c$$$         call sampleomega(pp2,e,bp,cpg,mw,gw,omin,omax,om,w)
c$$$         return
c$$$      endif
c$$$ 555  continue
ccc TEST
** bi-channel: q2 or om
      if (pp2.lt.0.d0.or.ymin.gt.ymax.or.pp(0).lt.0.d0) then 
          pbw = 0.d0
      endif
      pbw0 = pbw
      pbw  = pbw0*1.d0
      pbw2 = pbw0*0.d0
      pir  = 1.d0 - pbw - pbw2
      extrabw = 1.d0 /yjac
      extrair = 1.d0 / omax
      call wraprng(csi,1)  
      if (csi(1).le.pbw) then
         y    = 1.d0*rnd
         tmp  = atan(ymin*pp2/mw/gwl-mw/gwl)
         y    = gwl*mw/pp2*( tan(pp2*mw*gwl*anbw*y+tmp)+mw/gwl )
         om   = (1-y)*pp2/2.d0/pp(0)/umbc1p
         w     = w * anbw * ((y*pp2 - mw**2)**2+gwl**2*mw**2)
         regi  = 1.d0/( (y*pp2-mw**2)**2+gwl**2*mw**2 )/anbw * yjac
         regi  = regi * extrabw
         xap   = om/omax
      elseif(csi(1).le.(pbw+pbw2)) then
         y    = 1.d0*rnd
         tmp  = atan(ymin*pp2/mw/gwl2-mw/gwl2)
         y    = gwl2*mw/pp2*( tan(pp2*mw*gwl2*anbw2*y+tmp)+mw/gwl2 )
         om   = (1-y)*pp2/2.d0/pp(0)/umbc1p
         w     = w * anbw2 * ((y*pp2 - mw**2)**2+gwl2**2*mw**2)
         regi  = 1.d0/( (y*pp2-mw**2)**2+gwl2**2*mw**2 )/anbw2 * yjac
         regi  = regi * extrabw
         xap   = om/omax
      else
         call ap_vertex_sample(omin,omax,y,x)
         om   = omax*x
         q2   = pp2 - 2.d0 * om * umbc1p * pp(0)
         y    = q2/pp2
         xap  = x
         regi = (1.d0+(1.d0-xap)**2)/xap/anir
         w    = w / regi
      endif
      reg = regi/
     .     (pbw/((y*pp2-mw**2)**2+gwl**2*mw**2)/anbw*extrabw
     .     +pbw2/((y*pp2-mw**2)**2+gwl2**2*mw**2)/anbw2*extrabw
     .     +pir*(1.d0+(1.d0-xap)**2)/xap/anir*extrair)
      w   = w*reg
      return
      end
      subroutine ap_vertex_sample(omin,omax,x,omx)
      implicit double precision (a-h,o-z)
! written by CMCC, essentially from old BABAYAGA
*  x generation according to ap splitting function
*  (1+x^2)/(1-x), 0 <= x <= 1-eps
* added 15/10/2005, modified from babayaga, older releases by CMCC
      real*4 r(2)
      alne = log(omin/omax)
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
************
      subroutine energy_forphind(pbw,rnd,pp,vers,q2min,
     .     omin,omax,mw,gw,om,w)
      implicit double precision (a-h,m,o-z)
      real*4 rnd,csi(1)
      dimension pp(0:3),vers(0:3)
!from shared.inc
      common/gen_event/ebeam,ebeam1,ebeam2,xmin,xmax,shat,x1pdf,x2pdf
      common/yminimo/ymincut ! for x1,x2pdf sampling
      w = 1.d0
      bppc1pp = tridot(pp,vers)/pp(0)
      umbc1p  = 1.d0 - bppc1pp
      pp2   = dot(pp,pp)
      q2max = pp2 - 2.d0*omin*umbc1p*pp(0)
      ymax  = q2max/pp2
      ymin  = q2min/pp2
      ymin  = max(ymin,0.d0)
      yjac  =  pp2/2.d0/pp(0)/umbc1p

      gwl  = gw

      anbw  = 1.d0/(pp2*mw*gwl)*( atan( ymax*pp2/mw/gwl-mw/gwl )
     >     - atan( ymin*pp2/mw/gwl-mw/gwl ) )

** bi-channel: q2 or om
      if (pp2.lt.0.d0.or.ymin.gt.ymax.or.pp(0).lt.0.d0) then 
          pbw = 0.d0
      endif

      pir  = 1.d0 - pbw
      extrabw = 1.d0 /yjac

      call wraprng(csi,1)  
      if (csi(1).le.pbw) then
         y    = 1.d0*rnd
         tmp  = atan(ymin*pp2/mw/gwl-mw/gwl)
         y    = gwl*mw/pp2*( tan(pp2*mw*gwl*anbw*y+tmp)+mw/gwl )
         om   = (1-y)*pp2/2.d0/pp(0)/umbc1p
         w     = w * anbw * ((y*pp2 - mw**2)**2+gwl**2*mw**2)
         regi  = 1.d0/( (y*pp2-mw**2)**2+gwl**2*mw**2 )/anbw * yjac
         regi  = regi * extrabw
         xap   = om/omax
      else
         om = (omax - omin)*rnd + omin

         q2   = pp2 - 2.d0 * om * umbc1p * pp(0)
         y    = q2/pp2
         xap  = x
         regi = 1.d0/(omax-omin)
         w    = w / regi
      endif
      reg = regi/
     .     (pbw/((y*pp2-mw**2)**2+gwl**2*mw**2)/anbw*extrabw
     .     +pir/(omax-omin))
      w   = w*reg
      return
      end
***************
      subroutine photon_energy_notAP_withpow(pbw,rnd,pp,vers,q2min,
c      subroutine photon_energy(pbw,rnd,pp,vers,q2min,
     .     omin,omax,mw,gw,om,w)
      implicit double precision (a-h,m,o-z)
      real*4 rnd,csi(1)
      dimension pp(0:3),vers(0:3)
!from shared.inc
      common/gen_event/ebeam,ebeam1,ebeam2,xmin,xmax,shat,x1pdf,x2pdf
      common/yminimo/ymincut ! for x1,x2pdf sampling
      w = 1.d0
      bppc1pp = tridot(pp,vers)/pp(0)
      umbc1p  = 1.d0 - bppc1pp
      pp2 = dot(pp,pp)
      q2max = pp2 - 2.d0*omin*umbc1p*pp(0)
      ymax = q2max/pp2
      ymin = q2min/pp2
      ymin = max(ymin,0.d0)
      yjac =  pp2/2.d0/pp(0)/umbc1p

      gwl = gw
      anbw = 1.d0/(pp2*mw*gwl)*( atan( ymax*pp2/mw/gwl-mw/gwl )
     >     - atan( ymin*pp2/mw/gwl-mw/gwl ) )
      anir = log(omax/omin)
ccc TEST
c$$$      goto 555
c$$$      if (pp2.gt.0.d0) then
c$$$         bp  = sqrt(tridot(pp,pp))/pp(0)
c$$$         e   = sqrt(pp2)
c$$$         cpg = 0.d0
c$$$         if (bp.gt.0.d0) cpg = tridot(pp,vers)/sqrt(tridot(pp,pp))
c$$$         call sampleomega(pp2,e,bp,cpg,mw,gw,omin,omax,om,w)
c$$$         return
c$$$      endif
c$$$ 555  continue
ccc TEST
** bi-channel: q2 or om
c      ompeak = (pp2-mw**2)/2.d0/pp(0)/umbc1p
c      ompeak = max(ompeak,omin)
c      pa = 1.d0/ompeak*anbw
c      pb = 1.d0/((ymax*pp2-mw**2)**2+gwl**2*mw**2)*anir
c      pb = 2.d0*pb/(pa+2.d0*pb)
c      pa = 1.d0-pb
c      pbw = pa
      if (pp2.lt.0.d0.or.ymin.gt.ymax.or.pp(0).lt.0.d0) then 
          pbw = 0.d0
      endif

      a = 0.d0
      b = -a
      upb   = 1.d0 + b
      anpow = 1.d0/upb*(omax - omin)**upb

c      suman = anpow + anir + anbw
c      prob1 = anir/suman
c      prob2 = anbw/suman
c      prob3 = anpow/suman

      pir  = 1.d0 - pbw
      ppow = 0.d0


      extrabw = 1.d0 /yjac
      extrair = 1.d0
      call wraprng(csi,1)  
      if (csi(1).le.pbw) then
         y    = 1.d0*rnd
         tmp  = atan(ymin*pp2/mw/gwl-mw/gwl)
         y    = gwl*mw/pp2*( tan(pp2*mw*gwl*anbw*y+tmp)+mw/gwl )
         om   = (1-y)*pp2/2.d0/pp(0)/umbc1p
         w     = w * anbw * ((y*pp2 - mw**2)**2+gwl**2*mw**2)
         regi  = 1.d0/( (y*pp2-mw**2)**2+gwl**2*mw**2 )/anbw * yjac
         regi  = regi * extrabw
      elseif (csi(1).le.pbw+pir) then
         om   = omin * exp(1.d0*rnd*anir)
         w    = w * anir * om
         q2   = pp2 - 2.d0 * om * umbc1p * pp(0)
         y    = q2/pp2
         regi = 1.d0/om/anir
         regi = regi*extrair
      elseif (csi(1).le.pbw+pir+ppow) then
         om   = omax - (omax-omin)*(1.d0-rnd)**(1.d0/upb)
         q2   = pp2 - 2.d0 * om * umbc1p * pp(0)
         y    = q2/pp2
         regi = 1.d0!/(omax-om)**a/anpow
         w    = w *1.d0
      endif 
c$$$      reg = regi/(pbw/((y*pp2-mw**2)**2+gwl**2*mw**2)/anbw*extrabw
c$$$     .                +pir/om/anir*extrair
c$$$     .                +ppow/(omax-om)**a/anpow)
      aaa = (omax-om)**a
      reg = regi/(pbw/((y*pp2-mw**2)**2+gwl**2*mw**2)/anbw*extrabw*aaa
     .                +pir/om/anir*extrair*aaa
     .                +ppow/anpow)
      reg = reg * aaa
      w   = w*reg
      return
      end
*
      subroutine get_cos_fer(rnd,ng,c,w)
      implicit double precision (a-h,o-z)
      real*4 rnd,csi(1)
      character*1 boson
! from shared.inc
      double precision m1,m2,ch1,ch2,chfs,mfs
      common/vectorboson/boson
      common/partons/ipart1,ipart2
      common/gen_event/ebeam,ebeam1,ebeam2,xmin,xmax,shat,x1pdf,x2pdf
      common/processhor/m1,m2,ch1,ch2,chfs,mfs,lepton,imirror
      common/cosinesforcuts/cosenomax,cosenomin
      common/ifirstcos2gborn/anorm,al1pb1mb,ifirst
      data ifirst /0/
      iflat = 0

c      sqrts = 2.d0*sqrt(x1pdf*x2pdf)*ebeam
c      if (ng.ge.1.and.sqrts.gt.500d0) iflat = 1
c      iflat = 0

      if (iflat.eq.1) then
         c = 2.d0*rnd - 1.d0
c         cmax = cosenomax
c         cmin = cosenomin
c         c = (cmax-cmin)*rnd + cmin
         w = 2.d0
c         w = cmax-cmin
      else
         if (boson.eq.'W') then
            pflat = 0.03d0
            psing = 1.d0 - pflat
            csi(1) = pflat + 0.1
            if (pflat.gt.0.d0) call wraprng(csi,1)
            if (csi(1).lt.pflat) then
               c  = 2.d0 * rnd - 1.d0
            else
               c = 1.d0-2.d0*(1.d0-rnd)**(1.d0/3.d0)
            endif
            den = psing/(8.d0/3.d0/(1.d0-c)/(1.d0-c))
     .           +pflat/2.d0
            w = 1.d0/den
c - orig            c = 1.d0-2.d0*(1.d0-rnd)**(1.d0/3.d0)
c - orig            w = 8.d0/3.d0/(1.d0-c)/(1.d0-c)
         else  ! if (boson.eq.'Z')
            
            if (ipart1.eq.0.and.ipart2.eq.0) then ! gg -> mu+mu-!!!

* better than flat at high energies               
               e = sqrt(shat)/2.d0
               beta = 1.d0-mfs*mfs/e/e
               beta = sqrt(beta)
               ombet = 1.d0 - beta

               if (beta.eq.1.d0) then
                  xx    = -mfs*mfs/e/e
                  ombet = - 0.5d0*xx*(1.d0-xx*0.250d0+xx*xx*0.125d0
     .                 - 0.078125d0 * xx*xx*xx)
                  beta = 1.d0 - ombet
               endif
               
               al1pb1mb = log((1.d0+beta)/ombet)
               pippo    = log((1.d0+beta)/ombet)
               anorm    = 0.5d0/beta * (al1pb1mb+pippo)
               abig = 2.d0*beta*anorm*rnd - al1pb1mb
               c    = 1.d0/beta * (exp(abig)-1.d0)/(exp(abig)+1.d0)

               if (abs(c).ge.1.d0) then
                  c = c/abs(c)
               endif
               
               w    = anorm*(1.d0-beta*c)*(1.d0+beta*c)
*********************************************
c               c = 2.d0*rnd - 1.d0
c               w = 2.d0
               
            else
               c = ccmgen(w)
            endif

c         c = 2.d0*rnd - 1.d0
c         w = 2.d0
         endif
      endif
      return
      end
*
      subroutine collinearm(m,pflat,rnd1,rnd2,p,vers,w)
      implicit double precision (a-h,o-z)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      dimension p(0:3),vers(0:3),ptmp(0:3)
      real*4 rnd1,rnd2,csi(1)
      double precision m

      vers(0) = 1.d0
      phph = 2.d0*pi*rnd1

cc      pflat = 0.99d0

      if (pflat.gt.0.d0) then 
         call wraprng(csi,1)
      else
         csi(1) = 0.5
      endif

      psing = 1.d0 - pflat

      b1a = (1.d0-m/p(0))*(1.d0+m/p(0))
      b1a = sqrt(b1a)

c      pmod = sqrt(tridot(p,p))
c      b1 = pmod/p(0)         

      b1 = b1a

      an  = -1.d0/b1*log((1.d0-b1)/(1.d0+b1))

      if (csi(1).le.psing) then
         cth = (1.d0 - (1.d0+b1)*
     #        exp(-b1*an*rnd2))/b1

         umcth2 = abs((1.d0 + cth)*(1.d0 - cth)) ! abs to avoid num. problems..
         sth = sqrt(umcth2)
         cph = cos(phph)
         sph = sin(phph)

         vers(1) = vers(0)*sph*sth
         vers(2) = vers(0)*cph*sth
         vers(3) = vers(0)*cth
         call rot(-1,p,vers,vers)
         wi = dot(p,vers)/p(0) *an
      else
         cth = 2.d0*rnd2 - 1.d0
         sth = sqrt(1.d0-cth*cth)
         cph = cos(phph)
         sph = sin(phph)

         vers(1) = vers(0)*sph*sth
         vers(2) = vers(0)*cph*sth
         vers(3) = vers(0)*cth
         
         wi = 2.d0
      endif
c      den = psing*1.d0/dot(p,vers)*p(0)/an + pflat/2.d0
      den = psing*1.d0/(1.d0-b1*cth)/an + pflat/2.d0
      w = 2.d0*pi/den
      return
      end
*
      subroutine collinearm2(c1,c2,m,pflat,rnd1,rnd2,p,vers,w)
      implicit double precision (a-h,o-z)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      dimension p(0:3),vers(0:3),ptmp(0:3)
      real*4 rnd1,rnd2,csi1(1),csi2(1)
      double precision m

      vers(0) = 1.d0
      phph = 2.d0*pi*rnd1

      b = (1.d0-m*m/p(0)/p(0))
c      b = (1.d0-m/p(0))*(1.d0+m/p(0)) ! this is less accurate!
      
*** it can happen, for example for muon-neutrino-gamma final state
*** that ei > mmu + mnu + esoft, but p1(0)<mmu...!!!
      if (b.le.0.d0) b = 0.1d0
      b = sqrt(b)


*###########################################      
      if (b.ge.1.d0) then
         print*,'WARNING: sampling.f line 866'
      endif
*###########################################

      c12  = c1**2
      c12  = 0.d0
      c1c2 = c1*c2
      an  = 1.d0/b*log((1.d0+b)/(1.d0-b))
      an1 = an
      an  = an - 1.d0*c12
      
      istop = 0
      ic    = 0
      am2p02 = m*m/p(0)/p(0)
      do while(istop.eq.0)
         ic = ic + 1
         csi1(1) = rnd2
         if (ic.gt.1) call wraprng(csi1,1)
         call wraprng(csi2,1)
         cth = (1.d0 - (1.d0+b)*exp(-b*an1*csi1(1)))/b
         fcw  = 1.d0/(1.d0-b*cth)
         fcok = fcw*(1.d0 - 0.5d0*am2p02 * fcw * c12)
         if (fcok.lt.0.d0) print*,'fcok lt 0 ',fcok
         r = fcok/fcw
         if (csi2(1).le.r) istop = 1
      enddo

      umcth2 = abs((1.d0 + cth)*(1.d0 - cth)) ! abs to avoid num. problems..
      sth = sqrt(umcth2)
      cph = cos(phph)
      sph = sin(phph)
      
      vers(1) = vers(0)*sph*sth
      vers(2) = vers(0)*cph*sth
      vers(3) = vers(0)*cth
      call rot(-1,p,vers,vers)
      den = fcok/an
      w = 2.d0*pi/den
      return
      end
*
      subroutine collinear(pflat,rnd1,rnd2,p,vers,w)
      implicit double precision (a-h,o-z)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      dimension p(0:3),vers(0:3),ptmp(0:3)
      real*4 rnd1,rnd2,csi(1)
      double precision m
      vers(0) = 1.d0
      phph = 2.d0*pi*rnd1
      if (pflat.gt.0.d0) then 
         call wraprng(csi,1)
      else
         csi(1) = 0.5
      endif
      psing = 1.d0 - pflat
      pmod = sqrt(tridot(p,p))
      b1 = pmod/p(0)         
      an  = -1.d0/b1*log((1.d0-b1)/(1.d0+b1))
      if (csi(1).le.psing) then
         cth = (1.d0 - (1.d0+b1)*
     #        exp(-b1*an*rnd2))/b1
         umcth2 = abs((1.d0 + cth)*(1.d0 - cth)) ! abs to avoid num. problems..
         sth = sqrt(umcth2)
         cph = cos(phph)
         sph = sin(phph)
         vers(1) = vers(0)*sph*sth
         vers(2) = vers(0)*cph*sth
         vers(3) = vers(0)*cth
         call rot(-1,p,vers,vers)
         wi = dot(p,vers)/p(0) *an
      else
         cth = 2.d0*rnd2 - 1.d0
         sth = sqrt(1.d0-cth*cth)
         cph = cos(phph)
         sph = sin(phph)
         vers(1) = vers(0)*sph*sth
         vers(2) = vers(0)*cph*sth
         vers(3) = vers(0)*cth
         wi = 2.d0
      endif
      den = psing*1.d0/dot(p,vers)*p(0)/an + pflat/2.d0
      w = 2.d0*pi/den
      return
      end
*
      subroutine get_phi(rnd,ph,w)
      implicit double precision (a-h,o-z)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      real*4 rnd
      w  = 2.d0*pi
      ph = w*rnd
      return
      end

      function ccmgen(w)
* generate -1 < x < 1 according to 1+x^2 con un bicanale
      implicit double precision (a-h,o-z)
      real*4 random(2)      
      call wraprng(random,2)
c      ccmgen = 2.d0*random(1)-1.d0
c      w      = 2.d0
c      return
      an1   = 2.d0
      an2   = 2.d0/3.d0
      anorm = an1 + an2
      csi   = random(1)*anorm
      if (csi.lt.an1) then         
         x = an1*random(2)-1.d0
      else         
         arg  = 3.d0*(random(2)*an2-1.d0/3.d0)
         aarg = dabs(arg)
	 if (arg.gt.0.d0) then
	    segno =  1.d0
	 else
	    segno = -1.d0
	 endif      
         x = segno*aarg**(1.d0/3.d0)
      endif
      ccmgen = x
      w = anorm/(1.d0+x**2)
      return
      end
