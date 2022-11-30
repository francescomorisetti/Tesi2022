      function soft_integral(l,egmax,p1,p2,p3,p4,ch,masses)
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
         eta(1) = -1.d0 * ch(1)
         eta(2) =  1.d0 * ch(2)
         eta(3) =  1.d0 * ch(3)
         eta(4) = -1.d0 * ch(4)
      else
         eta(1) = -1.d0 * ch(1)
         eta(2) =  1.d0 * ch(2)
         eta(3) = -1.d0 * ch(3)
         eta(4) =  1.d0 * ch(4)
      endif
*
      linfra    = log(egmax/l)
      linfrapl2 = linfra + log(2.d0)

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
c            if ((i.eq.1.or.j.eq.2).and.(i.eq.3.or.j.eq.4)) then 
c               etaij = 0.d0
cc               continue
c            else
cc               etaij = 0.d0
c               continue
c            endif

********************
c            iok = 1
c            if (i.eq.1.and.j.eq.3) iok = 0
c            if (i.eq.1.and.j.eq.4) iok = 0
c            if (i.eq.2.and.j.eq.3) iok = 0
c            if (i.eq.2.and.j.eq.4) iok = 0
c            etaij = etaij * iok
*******************


            if (abs(etaij).gt.1.d-3) then

c OLD: gives NaNs                              
c               call rescale_momenta(npart,i,j,pmat,masses,p,q)
c               q2 = dot(q,q)
c               vl = dot(p,p) - q2
c               vl = vl/2.d0
c               v  = vl/(p(0) - q(0)) 
c               arglog = 1.d0+2.d0*vl/q2
c               terminfra = linfrapl2*log(arglog)         
c               tot = terminfra + termfinite(p,q,v)
c               tot = - tot *2.d0*dot(p,q)/vl* etaij
c               tot = tot * 2.d0 ! this is the double product when
c                                ! squaring the eikonal
c               softint  = softint + tot
c

c NEW: from Hto4l
c               totbck = tot
               call rescale_momenta_new(npart,i,j,pmat,masses,p,q,rho)
               q2 = masses(j)*masses(j)
               vl = 0.5d0*(rho*rho*masses(i)*masses(i) - q2)
               v  = vl/(p(0) - q(0)) 
               arglog = 1.d0+2.d0*vl/q2               
               terminfra = 0.d0
               if (arglog.gt.0.d0) then 
                  terminfra = linfrapl2*log(arglog)
                  tot = terminfra + termfinite_new(masses,p,q,v,rho,i,j)
               else
                  tot = 0.d0
               endif

               tot = -tot *2.d0*dot(p,q)/vl* etaij
               tot =  tot * 2.d0 ! this is the double product when
                                ! squaring the eikonal

               softint  = softint + tot
c               print*,'Z<<<',totbck,tot,totbck/tot
               
               
            endif
         enddo
      enddo      
      soft_integral = softint

!  diagonal contributions
      softint = 0.d0
      do  i = 1,npart
         etaii = eta(i)*eta(i)
         if (etaii.gt.1.d-3) then

c  betai=(1.d0-masses(i)/pmat(i,0))*(1.d0+masses(i)/pmat(i,0))            
            xx =  masses(i)/pmat(i,0)
            xx = -xx*xx

            betai = 1.d0 + xx
            betai = sqrt(betai)
            ombet = 1.d0 - betai

cc ! trying to avoid NaNs at extremely high energies      
            if (betai.eq.1.d0) then
               ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .              - 0.078125d0 * xx*xx*xx)
            endif
cc

            
            term = 2.d0*linfrapl2 + log(ombet/(1.d0+betai))/betai 
            term = -2.d0 * term * etaii
            softint = softint + term
         endif
      enddo
      soft_integral = soft_integral + softint
! divide by 2 to obtain exactly what is inside {} in Dittmaier & Kraemer
! paper (hep-ph/0109062)
      soft_integral = soft_integral*0.5d0

      return
      end
*
***********************************************************************************************
      function termfinite_new(masses,p,q,v,rho,i,j)
      implicit double precision (a-h,m,o-z)
      dimension p(0:3),q(0:3)
      double complex cspen,carg

      parameter (npart=4)
      dimension masses(npart),eta(npart)

      u0   = p(0)
      umod = sqrt(tridot(p,p))
      pp = u0 + umod
      pm = u0 - umod
      ppopm = pp*pp/rho/rho/masses(i)/masses(i) ! the same as pp/pm
      
      u0   = q(0)
      umod = sqrt(tridot(q,q))
      qp = u0 + umod
      qm = u0 - umod
      qpoqm = qp*qp/masses(j)/masses(j)  ! the same as qp/qm
      
      termfinite_new =
     .          log(ppopm)**2/4.d0 + ddiloghere((v-pm)/v) 
     .                             + ddiloghere((v-pp)/v)
     .        - log(qpoqm)**2/4.d0 - ddiloghere((v-qm)/v)
     .                             - ddiloghere((v-qp)/v)

      return
      end
*****************************************************************
***********************************************************************************************
      subroutine rescale_momenta_new(npart,i,j,pmat,masses,p,q,rho)
      implicit double precision (a-h,m,o-z)
      dimension p(0:3),q(0:3),pmat(npart,0:3),masses(npart)
      dimension p1(0:3),p2(0:3)
      
      do k = 0,3
         p1(k) = pmat(i,k)
         p2(k) = pmat(j,k)
      enddo

      m12 = masses(i)*masses(i)
      m22 = masses(j)*masses(j)
      
      p1p2 = dot(p1,p2)

      rho1 = p1p2 + sqrt(p1p2**2 - m12*m22)
      rho1 = rho1 / m12

c      rho2 = p1p2 - sqrt(p1p2**2 - m12*m22)
c      rho2 = rho2 / m12
!  better numerical solution !
      rho2 = m22/m12/rho1

      if ( (rho1*p1(0)-p2(0)) .gt. 0.d0) rho = rho1
      if ( (rho2*p1(0)-p2(0)) .gt. 0.d0) rho = rho2

      p = rho * p1
      q = p2

      return
      end
***********************************************************************************************
      
      subroutine rescale_momenta(npart,i,j,pmat,masses,p,q)
      implicit double precision (a-h,m,o-z)
      dimension p(0:3),q(0:3),pmat(npart,0:3),masses(npart)
      dimension p1(0:3),p2(0:3)
      
      do k = 0,3
         p1(k) = pmat(i,k)
         p2(k) = pmat(j,k)
      enddo

      m12 = masses(i)**2
      m22 = masses(j)**2
      
      p1p2 = dot(p1,p2)

      rho1 = p1p2 + sqrt(p1p2**2 - m12*m22)
      rho1 = rho1 / m12
      rho2 = p1p2 - sqrt(p1p2**2 - m12*m22)
      rho2 = rho2 / m12

cc BARBATRUCCO 
      rho1 =( m22 +rho1**2*m12)/2.d0/p1p2
      rho2 =( m22 +rho2**2*m12)/2.d0/p1p2
ccc
      if ( (rho1*p1(0)-p2(0)) .gt. 0.d0) rho = rho1
      if ( (rho2*p1(0)-p2(0)) .gt. 0.d0) rho = rho2
      do k = 0,3
         p(k) = rho * p1(k)
         q(k) = p2(k)
      enddo
      return
      end
c
      function termfinite(p,q,v)
      implicit double precision (a-h,o-z)
      dimension p(0:3),q(0:3)
      u0   = p(0)
      umod = sqrt(tridot(p,p))
      pp = u0 + umod
      pm = u0 - umod
      u0   = q(0)
      umod = sqrt(tridot(q,q))
      qp = u0 + umod
      qm = u0 - umod
c      termfin = log(pp/pm)**2/4.d0 + ddilog((v-pm)/v) + ddilog((v-pp)/v)
c     >        - log(qp/qm)**2/4.d0 - ddilog((v-qm)/v) - ddilog((v-qp)/v)
      termfin = log(pp/pm)**2/4.d0 + ddiloghere((v-pm)/v) 
     .                             + ddiloghere((v-pp)/v)
     >        - log(qp/qm)**2/4.d0 - ddiloghere((v-qm)/v) 
     .                             - ddiloghere((v-qp)/v)
      termfinite = termfin
      return
      end
*
      function softffDK(de,mg,mu,md,ml,s,t,u,qu,qd,ql,alpha)
c copied from Dittmaier & Kraemer paper...
      implicit double precision (a-h,m,o-z)
      pi = 4.d0 * atan(1.d0)
      ml2 = ml**2
      md2 = md**2
      mu2 = mu**2
      uno = ql**2*(2.d0*log(2.d0*de/mg)*log(ml2/s)+2.d0*log(2.d0*de/mg)
     .    + 0.5d0*log(ml2/s)**2 + log(ml2/s) + pi**2/3.d0 )
      due = qd**2*(2.d0*log(2.d0*de/mg)*log(md2/s)+2.d0*log(2.d0*de/mg)
     .    + 0.5d0*log(md2/s)**2 + log(md2/s) + pi**2/3.d0 )
      tre = qu**2*(2.d0*log(2.d0*de/mg)*log(mu2/s)+2.d0*log(2.d0*de/mg)
     .    + 0.5d0*log(mu2/s)**2 + log(mu2/s) + pi**2/3.d0 )
c      qua =-2.d0*ql*qd*(2.d0*log(2.d0*de/mg)*log(-t/s) - ddilog(-u/t))
c      cin = 2.d0*ql*qu*(2.d0*log(2.d0*de/mg)*log(-u/s) - ddilog(-t/u))
      qua=-2.d0*ql*qd*(2.d0*log(2.d0*de/mg)*log(-t/s)-ddiloghere(-u/t))
      cin= 2.d0*ql*qu*(2.d0*log(2.d0*de/mg)*log(-u/s)-ddiloghere(-t/u))
      softffDK = uno + due + tre + qua + cin
      softffDK = -softffDK * alpha/2.d0/pi
      return
      end
