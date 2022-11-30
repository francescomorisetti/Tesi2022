      function subtractionnew(k,scale,m1,m2,ch1,ch2,q)
      implicit double precision (a-h,m,o-z)
      dimension q(0:3)
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/subnewfirst/ifirst
* FOR NEW SUBTRACTION
      common/cosinesforsubtraction/cosu,cosd
      common/closesubtraction/iclosesub(40)
****   
      data ifirst/1/
      subtractionnew = 0.d0

      if (k.lt.0) return

      subtract1 = 0.d0
      subtract2 = 0.d0
      ik = iclosesub(k)
      cg = q(3)/q(0)

! ATTENZIONE
      q2 = scale**2/4.d0
      beta   = sqrt(1.d0 - m1**2/q2)
      
      ensub  = sqrt(q2)
      singis = ch1**2*(2.d0/(1.d0-beta*cg) -
     -     m1**2/ensub**2/(1.d0-beta*cg)**2)/q(0)**2
      
      subtractionnew = 4.d0*pi*alpha*singis
      subtract1 = subtractionnew
      
! ATTENZIONE
      q2 = scale**2/4.d0
      beta   = sqrt(1.d0 - m2**2/q2)
      ensub  = sqrt(q2)
      singis = ch2**2*(2.d0/(1.d0+beta*cg) -
     -     m2**2/ensub**2/(1.d0+beta*cg)**2)/q(0)**2
      
      subtractionnew = 4.d0*pi*alpha*singis
      subtract2 = subtractionnew
      subtractionnew = subtract1+subtract2
      ifirst = 0
      return
      end
*
      subroutine init_deltaq(scale)
      include 'shared.inc'
      double precision m1bck,m2bck
      common/toreadsubtraction/iread
      parameter (nxes = 5000)
      common/deltaqs/deltaq1(-ifl:ifl,nxes),deltaq2(-ifl:ifl,nxes)

      ipart1bck = ipart1
      ipart2bck = ipart2
      m1bck = m1
      m2bck = m2
      ch1bck = ch1
      ch2bck = ch2

      xpdfminl = log(xpdfmin)
      xpdfmaxl = log(xpdfmax)

      dstep = (xpdfmax - xpdfmin)/nxes
      dstepy = (xpdfmaxl - xpdfminl)/nxes

      x1 = xpdfmin-dstep
      x2 = xpdfmin-dstep
      y1 = xpdfminl-dstepy
      y2 = xpdfminl-dstepy

      if (boson.eq.'W')
     .     open(50,file=subtractionfile,status='unknown')
      if (boson.eq.'Z')
     .     open(50,file=subtractionfile,status='unknown')
      print*,'****************************'
      if (iread.eq.0) then
         print*,'Calculating subtraction array and'
         print*,'writing it to the file ',subtractionfile
         write(50,*)'***** SUBTRACTION FILE *********' 
         write(50,*)'iwhich ebeam pdfscale photon-induced'
         write(50,*)iwhich,ebeam,scale,iphinduced
         write(50,*)'********************************' 
      else
         print*,'Reading subtraction array from the file ',
     .        subtractionfile
         read(50,*)
         read(50,*)
         read(50,*)iw,eb,sc,ipind
         read(50,*)
         if (iw.ne.iwhich.or.
     .        abs(eb-ebeam).gt.1.d-4.or.
     .        abs(sc-scale).gt.1.d-4.or.ipind.ne.iphinduced) then
            print*,'Please, recompute the subtraction file!'
            print*,'variables in the subtraction file & for the run'
            print*,' [they must be the same...]'
            print*,'iwhich            = ',iw,iwhich
            print*,'Ebeam             = ',eb,ebeam
            print*,'pdf scale         = ',sc,scale
            print*,'photon-ind procs  = ',ipind,iphinduced
            print*,'Stopping Horace'
            stop
         endif
      endif
      print*,'n. points ',nxes
!      do k = 1,nxes-1
      do k = 1,nxes
         if (mod(k,500).eq.0) then
            if (iread.eq.0)
     .           print*,'Step n.',k,' of ',nxes,' scale',scale
         endif
         y1 = y1 + dstepy
         y2 = y2 + dstepy
         x1 = exp(y1)
         x2 = exp(y2)
         do ipart1 = -ifl,ifl
            ipart2 = ipart1
cc            if (ipart1.ne.0) then
               m1 = mq(ipart1)
               m2 = mq(ipart2)
               ch1 = chq(ipart1)
               ch2 = chq(ipart2)
               if (iread.eq.0) then ! calculates !
                  call tobesubtractedcuba(scale,x1,x2,q1,q2,bigq1,bigq2)
                  deltaq1(ipart1,k) = bigq1
                  deltaq2(ipart2,k) = bigq2
                  write(50,*)k,ipart1,bigq1,bigq2
               else
                  read(50,*)kkkkk,ipart1dddddd,bigq1,bigq2
                  deltaq1(ipart1,k) = bigq1
                  deltaq2(ipart2,k) = bigq2
               endif
cc            endif
         enddo
      enddo
      print*,'Done!'
      print*,'****************************'
      close(50)
      ipart1 = ipart1bck
      ipart2 = ipart2bck
      m1 = m1bck
      m2 = m2bck
      ch1 = ch1bck
      ch2 = ch2bck
      return
      end
*
      subroutine get_subtraction(x1,x2,dq1,dq2) 
      include 'shared.inc'
      parameter (nxes = 5000)
      common/deltaqs/deltaq1(-ifl:ifl,nxes),deltaq2(-ifl:ifl,nxes)
      common/getsubtra/icountsub
      common/fordebugging/idebugging
      data icountsub /0/
      icountsub = icountsub + 1

      y1 = log(x1)
      y2 = log(x2)

      xpdfminl = log(xpdfmin)
      xpdfmaxl = log(xpdfmax)

      k1 = (y1-xpdfminl)/(xpdfmaxl-xpdfminl) * nxes + 1
      k2 = (y2-xpdfminl)/(xpdfmaxl-xpdfminl) * nxes + 1

      yinf1 = xpdfminl + (xpdfmaxl - xpdfminl)/nxes*(k1 - 1)
      yinf2 = xpdfminl + (xpdfmaxl - xpdfminl)/nxes*(k2 - 1)

      ysup1 = xpdfminl + (xpdfmaxl - xpdfminl)/nxes*(k1)
      ysup2 = xpdfminl + (xpdfmaxl - xpdfminl)/nxes*(k2)

      xinf1 = exp(yinf1)
      xsup1 = exp(ysup1)

      xinf2 = exp(yinf2)
      xsup2 = exp(ysup2)
      
      delta1 = xsup1 - xinf1
      delta2 = xsup2 - xinf2

      a = deltaq1(ipart1,k1)
      if (k1.eq.nxes) then
         b = 0.d0
      else
         b = deltaq1(ipart1,k1+1)
      endif
      if (a.gt.-1.d-14) then 
         a = -abs(a) -1.d-14 
      endif
      if (b.gt.-1.d-14) then
         b = -abs(b) -1.d-14 
      endif
      dq1inf = log(-a)
      dq1sup = log(-b)

      a = deltaq2(ipart2,k2)
      if (k2.eq.nxes) then
         b = 0.d0
      else
         b = deltaq2(ipart2,k2+1)
      endif
      if (a.gt.-1.d-14) then
         a = -abs(a) -1.d-14 
      endif
      if (b.gt.-1.d-14) then
         b = -abs(b) -1.d-14 
      endif
      dq2inf = log(-a)
      dq2sup = log(-b)

      dq1l = dq1inf + (dq1sup - dq1inf)/delta1*(x1 - xinf1)
      dq2l = dq2inf + (dq2sup - dq2inf)/delta2*(x2 - xinf2)

      dq1 = -exp(dq1l)
      dq2 = -exp(dq2l)
      return
      end
******************************************************************************************
      subroutine tobesubtractedcuba(s,x1,x2,q1,q2,bigq1,bigq2)
      include 'shared.inc'
      integer dim
      integer*8 ispin
!      parameter (dim = 1024)
      parameter (dim = 256)
      double precision y1(dim),w1(dim),y2(dim),w2(dim)
      double precision result(2),error(2),prob(2)
      double precision iuser
      character*(*) cuhrefile
      parameter (cuhrefile = "")
      integer integrandcuba
      external integrandcuba
      common/pdfbackup/pdf1bck(-ifl:ifl),pdf2bck(-ifl:ifl)
      common/storewandy/y1,y2,w1,w2
      common/tbsfirstcubaoalpha/ifirst,icount
      data ifirst /1/
      data icount /0/
! anche shared.inc!!
      common/cubacommon/x1jac,x2jac,x1min,x2min,sc,epssub,ai1,ai2,
     .     x1c,x2c
      common/forphotonPDF/sumai1xpdf,sumai2xpdf

      if (ifirst.eq.1) then
         nthr = 0
         nthrmax = 1
c         call cubacores(nthr,nthrmax)
         ifirst = 0
      endif
      
      icount = icount + 1
      sc  = s
      x1c = x1
      x2c = x2
!cuba parameters
      ndim  = 1
      ncomp = 2
      epsrel = 1.d-4
      epsabs = 1.d-4
      iflags = 0
      iminev = 150
      imaxev = 500000
      ikey   = 0
!
      do k = -ifl,ifl
         pdf1bck(k) = pdf1(k)
         pdf2bck(k) = pdf2(k)
      enddo

      epssub = 1.d-4
      xplus  = 1.d0-epssub

      x1min = log(epssub)
      x2min = log(epssub)
      x1max = log(1.d0-x1)
      x2max = log(1.d0-x2)
      x1jac = x1max - x1min
      x2jac = x2max - x2min

      sum1 = 0d0
      sum2 = 0d0

      if (i_pdf.eq.1) then
         call get_pdf(x1,x2,s)
         q1 = max(0.d0,pdf1(ipart1))
         q2 = max(0.d0,pdf2(ipart2))
      else
         q1 = 1.d0
         q2 = 1.d0
      endif

! call to Cuba 1.2
c      call cuhre(ndim,ncomp,integrandcuba,epsrel,epsabs,iflags,
c     .     iminev,imaxev,ikey,nregions,neval,ierr,result,error,prob)

!  call to Cuba 4.2!
      nvecto = 1
      ispin = -1
      iuser = 0.d0
      ikey = 0
      call cuhre(ndim,ncomp,integrandcuba,iuser,nvecto,epsrel,epsabs,
     .     iflags,iminev,imaxev,ikey,cuhrefile,ispin,nregions,neval,
     .     ierr,result,error,prob)
      
*===  END CUBA ====

**SCOMMENTARE
      if (ipart1.ne.0.and.ipart2.ne.0) then
         bigq1 = result(1) - 1.d0*ai1*q1
         bigq2 = result(2) - 1.d0*ai2*q2
      else

         sumai1xpdf = 0.d0
         sumai2xpdf = 0.d0
         do k = -ifl,ifl
            if (k.ne.0) then
               ai1 = 4.5d0+pi**2/3.d0+1.5d0*log(epssub)-log(epssub)**2
!i.e., in Mathematica:
!Series[Integrate[(1+x^2)/(1-x)*(Log[(1-x)/x] - 3/4) + ((9+5x)/4), {x,0,1-eps}, Assumptions -> Im[eps]==0 && eps < 1 && eps > 0], {eps,0,1}]
               ai1 = -ai1 * alpha/2.d0/pi * chq(k)**2
               ai2 =  ai1
               sumai1xpdf = sumai1xpdf + ai1 * max(0.d0,pdf1bck(k)) ! mortacci
               sumai2xpdf = sumai2xpdf + ai2 * max(0.d0,pdf2bck(k))               
            endif
         enddo

         bigq1 = result(1) !- sumai1xpdf ! this is to be activated when term[34] are active below!!
         bigq2 = result(2) !- sumai2xpdf

      endif
**SCOMMENTARE

**COMMENTARE
c      bigq1 = result(1) - 0.d0*ai1*q1
c      bigq2 = result(2) - 0.d0*ai2*q2
**COMMENTARE

      do k = -ifl,ifl
         pdf1(k) = pdf1bck(k)
         pdf2(k) = pdf2bck(k)
      enddo
      ifirst = 0
      return
      end
**
c cuba1.2      subroutine integrandcuba(ndim,xx,ncomp,fun)
cuba 1.2 and 4.2
      integer function integrandcuba(ndim,xx,ncomp,fun)
      include 'shared.inc'
      integer iuser,nvec,core,ndim,ncomp
      common/cubacommon/x1jac,x2jac,x1min,x2min,sc,epssub,ai1,ai2,
     .     x1c,x2c
      common/forphotonPDF/sumai1xpdf,sumai2xpdf
      double precision xx(ndim),fun(ncomp)
      common/dismsbarss/lfc
******* SUBRACTION SCHEME !!!!
      lfc = 1  ! 0 -> MSbar, 1 -> DIS scheme
******* SUBRACTION SCHEME !!!!
      integrandcuba = 0
      
      fun(1) = 0.d0
      fun(2) = 0.d0
      
      ya = x1jac * xx(1) + x1min
      yb = x2jac * xx(1) + x2min
! trucco di Oreste
      z1 = 1.d0-exp(ya)
      z2 = 1.d0-exp(yb)
      if (ipart1.ne.0.and.ipart2.ne.0) then
         call tobeconvoluted(alpha,pi,z1,m1,ch1,sc**2,epssub,apv1,ai1)
         call tobeconvoluted(alpha,pi,z2,m2,ch2,sc**2,epssub,apv2,ai2)
      else
c         for photon PDF these are set below
         ai1 = 0.d0
         ai2 = 0.d0
         sumai1xpdf = 0.d0
         sumai2xpdf = 0.d0
c         ai1 =  4.5d0 + pi**2/3.d0 + 1.5d0*log(epssub) - log(epssub)**2
c         ai1 = -ai1
c         ai2 = ai1
         
      endif
      if (i_pdf.eq.1) then
         call get_pdf(x1c/z1,x2c/z2,sc)
         if (ipart1.ne.0.and.ipart2.ne.0) then
            pdfint1 = max(0.d0,pdf1(ipart1))/z1*(1.d0-z1)*x1jac
            pdfint2 = max(0.d0,pdf2(ipart2))/z2*(1.d0-z2)*x2jac
         endif
      else
         pdfint1 = 1.d0
         pdfint2 = 1.d0
      endif

!!! for photon induced!!

! color factor!!!
      nc = 3
      if (ipart1.ne.0.and.ipart2.ne.0) then
         as2p = alpha/2.d0/pi
         coll1 = log(sc**2/m1**2)
         coll2 = log(sc**2/m2**2)
         disms1=((1.d0-z1)**2+z1**2)*log((1.d0-z1)/z1)-8*z1**2+8*z1-1.d0
         disms2=((1.d0-z2)**2+z2**2)*log((1.d0-z2)/z2)-8*z2**2+8*z2-1.d0

*** SCOMMENTARE
         fun(1) = pdfint1 * apv1 ! only first line originally
! for photon PDF corrections and scheme, refer to Diener, Dittmaier,
! Hollik hep-ph/0509084
     .        + nc * iphinduced * max(0.d0,pdf1(0))/z1*(1.d0-z1)*x1jac
     .        * (( z1**2 + (1.d0-z1)**2 ) * coll1 + lfc*disms1)
     .        * as2p *ch1**2
********
         fun(2) = pdfint2 * apv2 ! only first line originally
     .        + nc * iphinduced * max(0.d0,pdf2(0))/z2*(1.d0-z2)*x2jac
     .        * (( z2**2 + (1.d0-z2)**2 ) * coll2 + lfc*disms2)
     .        * as2p *ch2**2
***SCOMMENTARE

***COMMENTARE
c         fun(1) = 0.d0 * pdfint1 * apv1
c     .        + nc * iphinduced * max(0.d0,pdf1(0))/z1*(1.d0-z1)*x1jac
c     .        * (( z1**2 + (1.d0-z1)**2 ) * coll1 + lfc*disms1)
c     .        * as2p *ch1**2
c         fun(2) = 0.d0 * pdfint2 * apv2
c     .        + nc * iphinduced * max(0.d0,pdf2(0))/z2*(1.d0-z2)*x2jac
c     .        * (( z2**2 + (1.d0-z2)**2 ) * coll2 + lfc*disms2)
c     .        * as2p *ch2**2
***COMMENTARE
      else
         sum1 = 0.d0
         sum2 = 0.d0
         sumai1xpdf = 0.d0
         sumai2xpdf = 0.d0
         do k = -ifl,ifl
            if (k.ne.0) then
                              
! no terms to switch DIS-MSBAR scheme here!
               coll1 = log(sc**2/mq(k)**2) - 1.d0 -2.d0*log(z1)
               term1 = coll1 * max(0.d0,pdf1(k))/z1*(1.d0-z1)*x1jac *
     .              (1.d0+(1.d0-z1)**2)/z1 * alpha/2.d0/pi * chq(k)**2

               term3 = (1.d0+z1**2)/(1.d0-z1)
     .              * (log((1.d0-z1)/z1) - 0.75d0) + (9.d0+5.d0*z1)/4.d0
               term3 = -term3 *  max(0.d0,pdf1(k))/z1*(1.d0-z1)*x1jac
     .                 * alpha/2.d0/pi * chq(k)**2

               sum1 = sum1 + term1 + 0.d0*term3

! no terms to switch DIS-MSBAR scheme here!
               coll2 = log(sc**2/mq(k)**2) - 1.d0 -2.d0*log(z2)
               term2 = coll2 * max(0.d0,pdf2(k))/z2*(1.d0-z2)*x2jac *
     .              (1.d0+(1.d0-z2)**2)/z2 * alpha/2.d0/pi * chq(k)**2

               term4 = (1.d0+z2**2)/(1.d0-z2)
     .              * (log((1.d0-z2)/z2) - 0.75d0) + (9.d0+5.d0*z2)/4.d0
               term4 = -term4 *  max(0.d0,pdf2(k))/z2*(1.d0-z2)*x2jac
     .                 * alpha/2.d0/pi * chq(k)**2

               sum2 = sum2 + term2 + 0.d0 * term4

            endif
         enddo
         extra = 1.d0
         if (boson.eq.'W') extra = 0.d0
         fun(1) = sum1 * iphinduced*extra
         fun(2) = sum2 * iphinduced*extra
      endif
      return
      end
c
      subroutine tobeconvoluted(alp,pi,z,m,ch,s,eps,apz,aiplus)
! - this is the (complete) function (including alpha,pi,L...) to be
!   convoluted with q(x/z,M^2)/z (refer to 3.2 by Dittmaier & Kraemer) 
! - aiplus is its integral between 0 and 1-eps
      implicit double precision (a-h,o-z)
      double precision m
      common/dismsbarss/lfc
      common/ifirsttobeconv/ifirst
      data ifirst/1/

      coll = log(s/m**2) - 1.d0

! 1 --- PURE AP VERTEX --------------------------
      apv    = (1.d0 + z**2)/(1.d0 - z)
      apvint = -2.d0*log(eps)-1.5d0+2.d0*eps-0.5d0*eps**2

      apz    = ch**2 * alp/2.d0/pi* coll * apv
      aiplus = ch**2 * alp/2.d0/pi* coll * apvint
      apzap  = apz
      aipap  = aiplus
c      aiplus = 0.d0
c      goto 444
      goto 443
!------------------------------------------------ 
      if (ifirst.eq.1) then
         do k = 1,10
            print*,'PURE AP IN THE SUBTRACTION!!!'
         enddo
         ifirst = 0
      endif
      goto 444
! 2 --- 3.2 Dittmaier and Kreamer --------------------------
      aieps =  -2.d0*log(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
      aiplog = -log(eps)**2 -2.d0*eps+2.d0+2.d0*eps*log(eps) - 0.25d0
     >        -eps**2*log(eps)-eps**2/2.d0
!!! con mathematica     >        -eps**2*log(eps)/2-eps**2/2.d0/

      extra  = coll - 2.d0*log(1.d0-z)
      apz    = ch**2 * alp/2.d0/pi*(1.d0 + z**2)/(1.d0 - z)*extra
      aiplus = ch**2 * alp/2.d0/pi*(coll*aieps - 2.d0*aiplog)
      apzdk    = apz
      aiplusdk = aiplus

ccc      goto 444
!------------------------------------------------ 
! 3 --- 23 and 24 Baur Keller Wackeroth --------------------------
! implementing BKW (formulae 23, 24) in order to have both DIS and MSbar
! schemes. For lcf = 0 this should be identical to DK 3.2, up to terms
! O(eps)...

 443  continue
      if (ifirst.eq.1.and.lfc.eq.0) then
         do k = 1,10
            print*,'MSbar scheme in subtraction of IS singularities...'
         enddo
         ifirst = 0
      endif

      fvps = 9.d0 + 2.d0*pi**2/3.d0 + 3.d0*log(eps) - 2.d0*log(eps)**2
      fc   = (1.d0 + z**2)/(1.d0 - z)*log((1.d0-z)/z) -
     .     3.d0/2.d0*1.d0/(1.d0 - z) + 2.d0 * z + 3.d0

      coll = log(s/m**2)

      aiplus = 1.d0 - log(eps)-log(eps)**2+(log(eps)+3.d0/4.d0)*coll
     .     - 1.d0/4.d0 * lfc * fvps

      aiplus = -aiplus * ch**2*alp/pi

      apz    = (1.d0 + z**2)/(1.d0 - z)*(coll - 2.d0*log(1.d0-z) - 1.d0)
     .     + lfc*fc
      apz    = alp/2.d0/pi*apz*ch**2
 444  return
      end
***


***************************************************
***************************************************
***************************************************
      subroutine gauleg(x1,x2,x,w,n)
      implicit none
      integer n
      double precision x1,x2,x(n),w(n),pi
      double precision EPS
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (EPS=3.d-14)

      integer i,j,m
      double precision p1,p2,p3,pp,xl,xm,z,z1

      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
*      print*,'m xm xl ',m,xm,xl

      do i=1,m
         z=cos(pi*(i-0.25d0)/(n+0.5d0))
 1       continue
         p1=1.d0
         p2=0.d0
         do j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
         end do
         pp=n*(z*p1-p2)/(z*z-1.d0)
         z1=z
         z=z1-p1/pp
         if (abs(z-z1).gt.EPS) goto 1
         x(i)=xm-xl*z
         x(n+1-i)=xm+xl*z
         w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
         w(n+1-i)=w(i)
      end do
      return
      end
      subroutine tobesubtractedcubapartonic(s,bigq1,bigq2)
      include 'shared.inc'
      integer dim
!      parameter (dim = 1024)
      parameter (dim = 256)
      double precision y1(dim),w1(dim),y2(dim),w2(dim)
      double precision result(2),error(2),prob(2)
      external integrandcubapartonic
      common/pdfbackup/pdf1bck(-ifl:ifl),pdf2bck(-ifl:ifl)
      common/storewandy/y1,y2,w1,w2
      common/tbsfirstcubapartonic/ifirst
      data ifirst /1/

! anche shared.inc!!
      common/cubacommon/x1jac,x2jac,x1min,x2min,sc,epssub,ai1,ai2,
     .     x1c,x2c

      call get_x1x2pdf(s,pdfjac)      
      call get_is(x1pdf,x2pdf,s,wpdfpdf)

      sc  = s
      x1c = x1
      x2c = x2

!cuba parameters
      ndim  = 1
      ncomp = 2
      epsrel = 1.d-10
      epsabs = 1.d-10
      iflags = 0
      iminev = 100000
      imaxev = 1000000
      ikey   = 0
!
      epssub = 1.d-4
      xplus  = 1.d0-epssub

      x1min = log(epssub)
      x2min = log(epssub)
      x1max = log(1.d0)
      x2max = log(1.d0)
      x1jac = x1max - x1min
      x2jac = x2max - x2min

*=== STARTING CUBA
      print*,'for partonic, should not pass here... cuhre commented out'
c 222  call cuhre(ndim,ncomp,integrandcubapartonic,epsrel,epsabs,iflags,
c     .     iminev,imaxev,ikey,nregions,neval,ierr,result,error,prob)
*=== END CUBA ====

      print*,'===== CUBA STATS ====='
      print*,error
      print*,prob
      print*,neval
      print*,'======================'

      saaa = 4.d0*ebeam**2

      s0 = 1.d0/4.d0/3.d0*convfac
      s0 = s0 * 2.d0*pi
      s0 = s0 * gf**2*mw**4/2.d0/pi/saaa/pi
      s0 = s0 * (saaa)**2/4.d0 *8.d0/3.d0
      s0 = s0 * 1.d0/((saaa-mw**2)**2 + gw**2*mw**2)

      bigq1 = (1.d0*result(1) - (ai1)*s0)
      bigq2 = (1.d0*result(2) - (ai2)*s0)
c      bigq1 = (result(1) - ai1*1.d0)
c      bigq2 = (result(2) - ai2*1.d0)
      return
      end
**
      subroutine integrandcubapartonic(ndim,xx,ncomp,fun)
      include 'shared.inc'
      common/cubacommon/x1jac,x2jac,x1min,x2min,sc,epssub,ai1,ai2,
     .     x1c,x2c
      double precision xx(*),fun(*)

      ya = x1jac * xx(1) + x1min
      yb = x2jac * xx(1) + x2min

      z1 = 1.d0-exp(ya)
      z2 = 1.d0-exp(yb)

      call tobeconvoluted(alpha,pi,z1,m1,ch1,sc**2,epssub,apv1,ai1)
      call tobeconvoluted(alpha,pi,z2,m2,ch2,sc**2,epssub,apv2,ai2)

      s = 4.d0*ebeam**2

      s01 = 1.d0/4.d0/3.d0*convfac
      s01 = s01 * 2d0*pi
      s01 = s01 * gf**2*mw**4/2.d0/pi/s/z1/pi
      s01 = s01 * (z1*s)**2/4.d0 *8.d0/3.d0
      s01 = s01 * 1.d0/((z1*s-mw**2)**2 + gw**2*mw**2)

      s02 = 1.d0/4.d0/3.d0*convfac
      s02 = s02 * 2d0*pi
      s02 = s02 * gf**2*mw**4/2.d0/pi/s/z2/pi
      s02 = s02 * (z2*s)**2/4.d0 *8.d0/3.d0
      s02 = s02 * 1.d0/((z2*s-mw**2)**2 + gw**2*mw**2)

      pdfint1 = (s01)*(1.d0-z1)*x1jac
      pdfint2 = (s02)*(1.d0-z2)*x2jac

c      pdfint1 = 1.d0*(1.d0-z1)*x1jac
c      pdfint2 = 1.d0*(1.d0-z2)*x2jac

      fun(1) = pdfint1 * (apv1)
      fun(2) = pdfint2 * (apv2)
      return
      end
