      subroutine setfilesnames
      include 'shared.inc'
      include 'nbin.inc'
      character*150 filesnames(ndistr),afbfile,afbtotfile,tmpname
      character*30  distname(ndistr)
      character*8   distorder,accel,bins,recom
      common/aawrite/ifirst_write,filesnames
      common/write_distributions_extra/afbfile,afbtotfile

      call resetname(tmpname)
      do k = 1,ndistr
         filesnames(k) = tmpname
         distname(k) = tmpname
      enddo
      afbfile = tmpname
      afbtotfile = tmpname
      distorder = tmpname
      accel = tmpname
      bins = tmpname
      recom = tmpname

      if (irec.eq.1) recom='rec_'
      if (irec.eq.0) recom='norec_'

** filling name of distribution files...
      if (charge_factor(3).lt.0.5d0) then 
         distorder='born.dat'
      else
         if (qedorder.eq.'exp')   distorder='best.dat'
         if (qedorder.eq.'alpha') distorder='oal1.dat'
      endif
      
      if (hadr1.eq.hadr2) accel='lhc_'
      if (hadr1.ne.hadr2) accel='tev_'
      if (i_pdf.eq.0)     accel='pl_'
      
      do k = 1,100
         if (path(k:k).ne.' ') lpath=k
      enddo
      do k = 1,ndistr
         call merge_strings(path,'empty_distr.dat',filesnames(k))
      enddo
      call merge_strings(path,'afb_',afbfile)
      call merge_strings(afbfile,distorder,afbfile)
      call merge_strings(path,'afbtot_',afbtotfile)
      call merge_strings(afbtotfile,distorder,afbtotfile)
***
      if (boson.eq.'OLD W DISTRIBUTIONS') then         
      else  ! if boson eq Z or W (new DISTRIBUTIONS!)
         call resetname(distorder)
         call resetname(afbfile)
         call resetname(afbtotfile)

         if (charge_factor(3).lt.0.5d0) then 
            distorder='born_'
         else
            if (qedorder.eq.'exp')   distorder='best_'
            if (qedorder.eq.'alpha') distorder='oal_'
         endif
         call itoa(nbin,bins)

**** DISTRIBUTIONS' NAMES
         distname(1)  = 'ptl1_'
         distname(2)  = 'ptl2_'
         distname(3)  = 'ptn1_'
         distname(4)  = 'mtv1_'
         distname(5)  = 'minv_'
         distname(6)  = 'etal_'
         distname(7)  = 'etan_'
         distname(8)  = 'yv_'
         distname(9)  = 'ptv1_'
         distname(10) = 'ptv2_'
         distname(11) = 'phistar_'
         distname(12) = 'minvext_'
         distname(13) = 'x1pdf_'
         distname(14) = 'logy_'
         distname(15) = 'costh_'
         distname(16) = 'minvfine_'
         distname(17) = 'ptg_'
         distname(18) = 'eneg_'
         distname(19) = 'xmt1_'
         distname(20) = 'xmt2_'
         distname(21) = 'xptl1_'
         distname(22) = 'xptl2_'
         distname(23) = 'minf_'
         distname(24) = 'minb_'
         distname(25) = 'himf_'
         distname(26) = 'himb_'
         
***********************
         do k = 1,ndistr
            call resetname(tmpname)
            call merge_strings(path,tmpname,tmpname)
            call merge_strings(tmpname,accel,tmpname)
            if (boson.eq.'Z') call merge_strings(tmpname,'z_',tmpname)
            if (boson.eq.'W') call merge_strings(tmpname,'w_',tmpname)
            if (lepton.eq.1) call merge_strings(tmpname,'el_',tmpname)
            if (lepton.eq.2) call merge_strings(tmpname,'mu_',tmpname)
            if (lepton.eq.3) call merge_strings(tmpname,'ta_',tmpname)
            call merge_strings(tmpname,distname(k),tmpname)
            call merge_strings(tmpname,distorder,tmpname)
            call merge_strings(tmpname,recom,tmpname)
            call merge_strings(tmpname,bins,tmpname)
            call merge_strings(tmpname,'.dat',tmpname)
            filesnames(k) = tmpname
         enddo
         call merge_strings(path,afbfile,afbfile)
         call merge_strings(path,afbtotfile,afbtotfile)
         call merge_strings(afbfile,accel,afbfile)
         call merge_strings(afbtotfile,accel,afbtotfile)
         call merge_strings(afbfile,'afb_',afbfile)
         call merge_strings(afbtotfile,'afbtot_',afbtotfile)
         if (lepton.eq.1) call merge_strings(afbfile,'el_',afbfile)
         if (lepton.eq.2) call merge_strings(afbfile,'mu_',afbfile)
         if (lepton.eq.3) call merge_strings(afbfile,'ta_',afbfile)
         if (lepton.eq.1) 
     >        call merge_strings(afbtotfile,'el_',afbtotfile)
         if (lepton.eq.2) 
     >        call merge_strings(afbtotfile,'mu_',afbtotfile)
         if (lepton.eq.3) 
     >        call merge_strings(afbtotfile,'ta_',afbtotfile)
         call merge_strings(afbfile,distorder,afbfile)
         call merge_strings(afbtotfile,distorder,afbtotfile)
         call merge_strings(afbfile,recom,afbfile)
         call merge_strings(afbtotfile,recom,afbtotfile)
         call merge_strings(afbfile,bins,afbfile)
         call merge_strings(afbtotfile,bins,afbtotfile)
         call merge_strings(afbfile,'.dat',afbfile)
         call merge_strings(afbtotfile,'.dat',afbtotfile)
*****
      endif
c      print*,'setfilenames'
      return
      end
*****************************************************
      subroutine get_functions(p3,p4,qph,fun)
* Cross section is binned in nbin bins differentiated in these functions
      include 'shared.inc'
      include 'nbin.inc'
      common/yminimo/ymin
      dimension bmi(ndistr),bma(ndistr),fun(ndistr)
      double precision qa(0:3),qb(0:3),qw(0:3)
      double precision pippo(0:3),pluto(0:3),topolino(0:3)
      dimension p3cm(0:3),p4cm(0:3),qphcm(0:3),p1cm(0:3),p2cm(0:3)
      common/eventinCOM/p3cm,p4cm,qphcm
      common/momentainitial/p1cm,p2cm
      common/fordebugging/idebugging
      common/bmibma/bmi,bma
      common/getfunfirst/ifirst
CARLO
      dimension vboost(0:3),p3cmnew(0:3),p4cmnew(0:3)
      common/labp1p2/p1,p2
CARLO


      data ifirst /1/
      if (ifirst.eq.1) then
         do k = 1,ndistr
            bmi(k) = 0.d0
            bma(k) = 1.d0
            fun(k) = 0.5d0
         enddo
      endif
c-----------------------------------
c  with 200 bins we have to set ranges larger than what will be plotted
c------------------------------------
      if (ifirst.eq.1) then
         bmi(1) =   20.d0
         bma(1) =   60.d0
         bmi(2) =  50.d0
         bma(2) = 10050.d0
         bmi(3) =   25.d0
         bma(3) =   75.d0
         bmi(4) =   40.d0
         bma(4) =  100.d0
         bmi(5) =   80.d0
         bma(5) =  220.d0
         bmi(6) =   -2.d0
         bma(6) =   2.d0
         bmi(7) =   -2.d0
         bma(7) =   2.d0
         bmi(8) =   -0.5d0
         bma(8) =   0.5d0
         bmi(9) =    0.d0
         bma(9) =   50.d0
         bmi(10)=    0.d0
         bma(10)=  200.d0
         bmi(11)=    0.d0
         bma(11)=    2.d0
         bmi(12)= 100.d0
         bma(12)= 20100.d0
c         bmi(12)=   -6.d0
c         bma(12)=   14.d0
         bmi(13)=    0.d0
         bma(13)=    1.d-2
         bmi(14)=   -6.d0
         bma(14)=   14.d0
         bmi(15)=   -1.d0
         bma(15)=    3.d0
         bmi(16)=   86.d0
         bma(16)=   96.d0
         bmi(17)=    0.d0
         bma(17)=  100.d0
         bmi(18)=    0.d0
         bma(18)=  100.d0
         bmi(19)=    0.6d0
         bma(19)=    1.8d0
         bmi(20) =   0.6d0
         bma(20) =   1.8d0
         bmi(21) =   0.4d0
         bma(21) =   2.0d0
         bmi(22) =   0.4d0
         bma(22) =   2.0d0
         bmi(23) =   60d0    ! not used
         bma(23) =   180d0    ! not used
         bmi(24) =   60d0    ! not used
         bma(24) =   180d0    ! not used
         bmi(25) =   150d0    ! not used
         bma(25) =   10150d0    ! not used
         bmi(26) =   150d0    ! not used
         bma(26) =   10150d0    ! not used


      endif
c      print*,'inizia'
      fun(1) = sqrt(p3(1)**2 + p3(2)**2)
      fun(2) = fun(1)
      fun(3) = sqrt(p4(1)**2 + p4(2)**2)

      pippo(0) = fun(1)
      pluto(0) = fun(3)
      do k = 1,2
         pippo(k) = p3(k)
         pluto(k) = p4(k)
      enddo
      pippo(3) = 0.d0
      pluto(3) = 0.d0
      do k = 0,3
         topolino(k) = pippo(k) + pluto(k)
      enddo
      tmass2 = abs(dot(topolino,topolino))
      fun(4) = sqrt(tmass2)      

      do k = 0,3
         topolino(k) = p3(k)+p4(k)
      enddo
      fun(5) = sqrt(dot(topolino,topolino))

      fun(6) = rapidity(p3)
      fun(7) = rapidity(p4)
      fun(8) = rapidity(topolino)

      fun(9) = topolino(1)**2 + topolino(2)**2
      fun(9) = sqrt(fun(9)) !!!!!
      fun(10) = fun(9)


      ctseta = tanh(0.5d0*(fun(6)-fun(7)))
      stseta = sqrt(1d0-ctseta**2)
      cosdphi = tridot(p3,p4)/sqrt(tridot(p3,p3))/sqrt(tridot(p4,p4))
      dphi=acos(cosdphi)
      fun(11)=tan((pi-dphi)*0.5d0)*stseta

c      print*,ctseta,cstar(p3,p4),cosdphi,dphi,pi,fun(11)



c      print*,'prime dieci'
      l10=2.302585092994045d0   !!!  log(10)    log_10(x)=log(x)/log(10)
c      call fotonepiuenergetico(qph,qa)

      qa(0)=0d0
      qa(1)=0d0
      qa(2)=0d0
      qa(3)=0d0

      do iga = 1,10
         pippo(0)=qph(iga,0)
         pippo(1)=qph(iga,1)
         pippo(2)=qph(iga,2)
         pippo(3)=qph(iga,3)

         pippoy = rapidity(pippo)

         if ((pippo(0).gt.0.0004d0).and.(pippoy.lt.3d0)) then
            qa(0)=qa(0)+pippo(0)
            qa(1)=qa(1)+pippo(1)
            qa(2)=qa(2)+pippo(2)
            qa(3)=qa(3)+pippo(3)
         end if
      end do

         fun(12) = fun(5)

      egmin = 0.d0
      etag = rapidity(qa)
      if (qa(0).gt.egmin.and.abs(etag).lt.10.d0) then
         deta2 = (fun(6)-etag)**2

         xminnie = sqrt(qa(1)**2+qa(2)**2+qa(3)**2)*
     -            sqrt(p3(1)**2+p3(2)**2+p3(3)**2)
         cosphi = -(dot(p3,qa)-p3(0)*qa(0))/xminnie
 	 dphi2 = acos(cosphi)**2

         rlg = sqrt(deta2+dphi2)

c         print*,'rlg=',rlg
c         fun(11) = sqrt(rlg)
c         fun(12) = log(rlg)/l10

         yexp    = qa(0)/(qa(0)+p3(0))
c         print*,'yexp=',yexp
         fun(13) = yexp**(1.d0/3.d0)
         fun(14) = log(yexp)/l10

      endif

      fun(13)=x1pdf

CARLO

      bis = (p1(3)+p2(3))/(p1(0)+p2(0))
      gis =  1.d0/sqrt(1.d0-bis*bis)
c      if (imirror.eq.0) then
         vboost(0) = 0.d0
         vboost(1) = 0.d0
         vboost(2) = 0.d0
         vboost(3) = bis
c      else
         continue
c         vboost(0) = 0.d0
c         vboost(1) = 0.d0
c         vboost(2) = 0.d0
c         vboost(3) = bis
c      end if

*-- oppure vboost(3) = -bis !??????

      call boost(gis,vboost,p3,p3cmnew)
      call boost(gis,vboost,p4,p4cmnew)

c      print*,p4cmnew(0),p3cmnew(0)
c      print*,p4cmnew(1),p3cmnew(1)
c      print*,p4cmnew(2),p3cmnew(2)
c      print*,p4cmnew(3),p3cmnew(3)
c      print*,' '

      fun(15) = p3cmnew(3)/sqrt(tridot(p3cmnew,p3cmnew))
CARLO

c      if (imirror.eq.0) then
c         xminnie = sqrt(p1cm(1)**2+p1cm(2)**2+p1cm(3)**2)*
c     -             sqrt(p3cm(1)**2+p3cm(2)**2+p3cm(3)**2)
c         fun(15) = -(dot(p3cm,p1cm)-p3cm(0)*p1cm(0))/xminnie
c         fun(15) = tridot(p3cm,p1cm)/xminnie
c         fun(15)=p3cm(3)/tridot(p3cm,p3cm)
c      else
c         continue
c         xminnie = sqrt(p2cm(1)**2+p2cm(2)**2+p2cm(3)**2)*
c     -             sqrt(p3cm(1)**2+p3cm(2)**2+p3cm(3)**2)
cc         fun(15) = -(dot(p3cm,p2cm)-p3cm(0)*p2cm(0))/xminnie
c         fun(15) = tridot(p3cm,p2cm)/xminnie
c         fun(15)=p3cm(3)/tridot(p3cm,p3cm)
c      end if

      fun(16) = fun(5)



      fun(17) = sqrt(qa(1)**2+qa(2)**2)
      fun(18) = qa(0)
c      print*,'masses=',mw,mz
      fun(19) = fun(4)/mw
      fun(20) = fun(4)/mz
      fun(21) = fun(1)/mw
      fun(22) = fun(1)/mz

      fun(23) = x1pdf   ! I do not use this fun , but fun(12) and fun(14)

      ccs = cstar(p3,p4)
      if (ccs.ge.0.d0) then 
         fun(23) = fun(5)
         fun(25) = fun(5)
      else
         fun(24) = fun(5)
         fun(26) = fun(5)
      endif


c      print*,'finisce'
      ifirst = 0
      return
      end
*********************************************************************      
      subroutine resetname(name)
      integer i
      character*(*) name
      do i = 1,len(name)
         name(i:i) = ' '
      enddo
      return
      end
**
      function ifirstemptychar(name)
      integer i,k
      character*(*) name
      i = 1
      do k = 1,len(name)
         if (name(k:k).ne.' ') i=k
      enddo
      ifirstemptychar = i ! + 1
      return
      end
**
      subroutine merge_strings(str1,str2,merged)
      character*(*) merged,str1,str2
      integer n1
      n1 = ifirstemptychar(str1)
      if (str1(n1:n1).eq.' ') then
         merged=str2
      else
         merged=str1(:n1)//str2
      endif
      return
      end
**
      subroutine itoa(int,a)
      implicit integer (h-n)
      implicit character (a-g)
      character*(*) a
      parameter (nbase=10)

      call resetname(a)
      ichar0 = ichar('0')
      idiv = nbase
      i    = int
      ncifre    = 1
      do while((i/idiv).gt.0)
         idiv = idiv * nbase
         ncifre = ncifre+1
      enddo
      idiv = idiv/nbase

      do k =1,ncifre
         icifra = i/idiv
         ia = icifra + ichar0
         a(k:k) = char(ia)
         i = i - icifra*idiv
         idiv = idiv/nbase
      enddo      
      return
      end
************************************************+
      subroutine fptoa(d,a)
      implicit integer (h-n)
      double precision d,dec,tmp
      character*(*) a
      character*30 adec
      parameter (ndecimali=1000000000)
      parameter (nbase=10)

      intpart = d
      call itoa(intpart,a)

      dec = d-intpart
      idecimali = dec*ndecimali

      tmp = intpart
      tmp = tmp*ndecimali
      tmp = tmp + idecimali

      if (anint(d*ndecimali-tmp).gt.0) idecimali = idecimali + 1

      call itoa(idecimali,adec)

      nleadzeri = 0
      i = 1
      icont = 1
      ictrl = 0
      do while(i.lt.ndecimali.and.idecimali.gt.0)
         i = i*nbase
         ictrl = dec * i - ictrl
         if (ictrl.ge.1) icont = 0
         if (ictrl.eq.0) nleadzeri = nleadzeri + 1*icont
      enddo

      ladec = len(adec)
      istop = 0
      do k=1,ladec
         i = ladec+1 - k
         if (adec(i:i).ne.' ') then
            if (adec(i:i).ne.'0') istop = 1
            if (adec(i:i).eq.'0'.and.istop.eq.0) adec(i:i)=''
         endif
      enddo
      
      call merge_strings(a,'.',a)

      if (nleadzeri.gt.0) then
         do k=1,nleadzeri
            call merge_strings(a,'0',a)
         enddo
      endif

      call merge_strings(a,adec,a)
      return
      end
*****************************************************
      subroutine resetalldistributions
      include 'shared.inc'
      include 'nbin.inc'
      common/distr/distr,s_distr
      common/bmibma/bmi,bma
      common/fordebugging/idebugging
      common/ionlyfirsttime/ifirst
c      data ifirst /0/ 
      dimension s_distr(ndistr,ibin)
      dimension distr(ndistr,ibin),bmi(ndistr),bma(ndistr),fun(ndistr)
c      data ((distr(i,j),i=1,ndistr),j=1,ibin) /iperj * 0.d0/	
c      data ((s_distr(i,j),i=1,ndistr),j=1,ibin) /iperj * 0.d0/	
      integer*8 nc,passnc
      common/passnormalization/passnc


      do i = 1,ndistr
         do j = 1,ibin
            distr(i,j) = 0.d0
            s_distr(i,j) = 0.d0
         enddo
      enddo

      return
      end

*****************************************************
      subroutine distributions(sd,nc,p3,p4,qph)
      include 'shared.inc'
      include 'nbin.inc'
      common/distr/distr,s_distr
      common/bmibma/bmi,bma
      common/fordebugging/idebugging
      common/ionlyfirsttime/ifirst
      data ifirst /0/ 
      dimension s_distr(ndistr,ibin),s_ddistr(1,ibin,ibin)
      dimension ddistr(1,ibin,ibin)
      dimension distr(ndistr,ibin),bmi(ndistr),bma(ndistr),fun(ndistr)
      data ((distr(i,j),i=1,ndistr),j=1,ibin) /iperj * 0.d0/	
      data ((s_distr(i,j),i=1,ndistr),j=1,ibin) /iperj * 0.d0/	
      integer*8 nc,passnc
      common/passnormalization/passnc
      passnc = nc

      call get_functions(p3,p4,qph,fun)

***********************************************************************
***********************************************************************
! filling distr. components when entering the first time
      if (ifirst.eq.0) then
         do k = 1,ndistr-1
            d = (bma(k) - bmi(k))/nbin
            y = bmi(k)
            do i = 1,nbin
               y = y + d
               distr(k,i*3-2)   = y - d
               s_distr(k,i*3-2) = y - d
            enddo
         enddo

            d1 = (bma(12) - bmi(12))/nbin
            d2 = (bma(12) - bmi(14))/nbin
            y1 = bmi(12)
            y2 = bmi(14)
            do i = 1,nbin
               y1 = y1 + d1
               do j = 1,nbin
                  y2 = y2 + d2

                  ddistr(1,i*3-2,j*3-2)   = y1 - d1
                  s_ddistr(1,i*3-2,j*3-2) = y1 - d1
                  ddistr(1,i*3-2,j*3-2)   = y2 - d2
                  s_ddistr(1,i*3-2,j*3-2) = y2 - d2

               enddo
               y2 = bmi(14)

            enddo

      endif
***************************************
      DO k = 1,ndistr-1
         bmax = bma(k)
         d = (bma(k) - bmi(k))/nbin
         x = bmi(k)
         if (fun(k).lt.bma(k).and.fun(k).gt.bmi(k)) then
            i = (fun(k) - bmi(k))/d
            i = (i+1)*3
            distr(k,i-1) = distr(k,i-1) + sd
            distr(k,i)   = distr(k,i)   + sd**2
         endif
      ENDDO


         d1 = (bma(12) - bmi(12))/nbin
         d2 = (bma(14) - bmi(14))/nbin

         if (fun(12).lt.bma(12).and.fun(12).gt.bmi(12).and.
     -       fun(14).lt.bma(14).and.fun(14).gt.bmi(14) ) then
            i1 = (fun(12) - bmi(12))/d1
            i1 = (i1+1)*3
            i2 = (fun(14) - bmi(14))/d2
            i2 = (i2+1)*3
            ddistr(1,i1-1,i2-1) = ddistr(1,i1-1,i2-1) + sd
            ddistr(1,i1,i2)   = ddistr(1,i1,i2)   + sd**2
         endif


      ifirst = 1          
      return
      end
*********************************************************************
      subroutine writedistributions
      include 'shared.inc'
      include 'nbin.inc'
      dimension s_distr(ndistr,ibin),bmi(ndistr),bma(ndistr)
      dimension distr(ndistr,ibin),s_ddistr(1,ibin,ibin)
      dimension ddistr(1,ibin,ibin)
      real*4 xplot(nbin),yplot(nbin)
      character*150 filesnames(ndistr),afbfile,afbtotfile
      common/aawrite/ifirst_write,filesnames
      common/write_distributions_extra/afbfile,afbtotfile
      common/distr/distr,s_distr
      common/bmibma/bmi,bma
      data ifirst_write /1/
      integer*8 nc,passnc
      common/passnormalization/passnc
      nc = passnc
*****
      if (ifirst_write.eq.1) then 
         call setfilesnames
      endif
*****      
      do k = 1,ndistr-1

         open(12,file=filesnames(k),status='unknown')
         do i=1,nbin
            su  = distr(k,3*i-1)
            su2 = distr(k,3*i)
            s_distr(k,3*i-1) = su / nc
            argument = abs((su2/nc-s_distr(k,3*i-1)**2)/nc)
            s_distr(k,3*i) = sqrt(argument)
            sezd  = s_distr(k,i*3-1) /(bma(k)-bmi(k))*nbin
            esezd = s_distr(k,i*3) /(bma(k)-bmi(k))*nbin
            write(12,*)s_distr(k,i*3-2),sezd,esezd
c            if (k.eq.1) then
c               xplot(i) = s_distr(k,i*3-2)
c               yplot(i) = sezd
c            endif
	 enddo
         close(12)
      enddo

         open(12,file=filesnames(ndistr),status='unknown')
         do i=1,nbin
            do j=1,nbin
               su  = ddistr(1,3*i-1,3*j-1)
               su2 = ddistr(1,3*i,3*j)
               s_ddistr(1,3*i-1,3*j-1) = su / nc
               argument = abs((su2/nc-s_ddistr(1,3*i-1,3*j-1)**2)/nc)
               s_ddistr(1,3*i,3*j) = sqrt(argument)
               sezd  =s_ddistr(1,i*3-1,j*3-1)/(bma(12)-bmi(12))/
     -                                        (bma(14)-bmi(14))*nbin**2
               esezd = s_ddistr(1,i*3,j*3) /(bma(12)-bmi(12))/
     -                                      (bma(14)-bmi(14))*nbin**2

               write(12,*)s_ddistr(1,i*3-2,j*3-2),sezd,esezd
            enddo
	 enddo
         close(12)


*********
* PLOTTING via fifo with gnuplot at runtime!! use gnuplot < fifoplot
c      iplot = 3
c      if (ifirst_write.eq.1) open(30,file='fifoplot',status='unknown')
c      write(30,*)'plot ''1.dat'', ''2.dat'''
c      call flush(30)
*********
      ifirst_write = 0
      return
      end
*************************************************
      subroutine fotonepiuenergetico(q,q1)
! written by CMCC, last modified 9/10/2005
      implicit real*8 (a-h,o-z)
      dimension q(40,0:3),q1(0:3)
*  LEADING ENERGETIC PHOTON IS EXTRACTED
      q1(0) = 0.d0
      q1(1) = 0.d0
      q1(2) = 0.d0
      q1(3) = 0.d0

      ENHARD = 0.d0
      J=1
*
      DO I = 1,40
         if (q(i,0).gt.0.d0) then
            ENPHOT = Q(I,0)
            IF (ENHARD.GE.ENPHOT) THEN
               ENHARD = ENHARD
            ELSE
               J = I
               ENHARD = ENPHOT
            ENDIF
         endif
      ENDDO
*
      if (enhard.gt.0.d0) then
         DO I = 0,3
            Q1(I) = Q(J,I)
         ENDDO
      endif
      return
      end
***************************
      subroutine duefotonipiuenergetici(q,q1,q2)
      implicit double precision (a-h,o-z)
      dimension q(40,0:3),q1(0:3),q2(0:3)
*  leading energetic photon is extracted
        do i = 0,3
           q1(i) = 0.d0
           q2(i) = 0.d0
        enddo
        sum = 0.d0
	do k=1,40
          sum = sum + q(k,0)
        enddo
        if (sum.lt.1.d-8) return

        enphot = q(1,0)
        enhard = enphot
        j=1
        do i = 1,40
           enphot = q(i,0)
           if (enhard.ge.enphot) then
              enhard = enhard
           else
              j = i
              enhard = enphot
           endif
        enddo 
        enhardl = enhard
        jl = j
        if (j.gt.0) then
           do i = 0,3
              q1(i) = q(j,i)
           enddo
        endif
*  next-to-leading energetic photon is extracted
      enhard = 0.d0
      j = 0
      do i = 1,40
         if (i.ne.jl) then
            enphot = q(i,0)
            if (enhard.ge.enphot) then
               enhard = enhard
            else
               j = i
               enhard = enphot
            endif
         endif
      enddo 
      enhard = enhard
      if (j.gt.0) then
         do i = 0,3
            q2(i) = q(j,i)
         enddo
      endif
      return
      end
*********************************************
