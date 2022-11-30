      subroutine pdfsampling(s,spdf,w)
      include 'shared.inc'
      call get_x1x2pdf(s,pdfjac)
      if (boson.eq.'W') then 
         spdf = mw * rescalepdfscale
      else
         spdf = mz * rescalepdfscale
      endif
      call get_is(x1pdf,x2pdf,spdf,wpdfpdf)
      w = wpdfpdf*pdfjac
      return
      end
*****************************************************
      subroutine bwstandard(q2max,q2min,s,w)
      implicit double precision (a-h,o-z)
      real*4 csi(1)
      common/mvgvbwaMC/am,ag
      common/bwsamplecmn/aa,am2,amag,ooamag,ooaa,ifirsta
      common/bwsamplecmn2/emax,emin,F1min,F2min,d1,d2,ds
      data ifirsta /0/
      if (ifirsta.eq.0) then
         am2    = am*am
         amag   = am*ag
         ooamag = 1.d0/amag
         aa     = am2/amag
         ooaa   = 1.d0/aa
         emin   = (q2min - am2)*ooamag
         emax   = (q2max - am2)*ooamag
         
         F1min = aa*atan(emin)
         F2min = 0.5d0 * log(emin*emin+1.d0)
         
         d1   = aa * atan(emax) - F1min
         d2   = 0.5d0 * log(emax*emax+1.d0) - F2min
         ds   = d1 + d2         
         ifirsta = 1
      endif

      call wraprng(csi,1)
      x = d1 * csi(1) + F1min
      e = tan(x*ooaa)
      s = amag * e + am2      
      w = d1 * ((s-am2)*(s-am2) + amag*amag )
      return
      end
      

      function bwnormstandard(q2max,q2min)
      implicit double precision (a-h,o-z)
      common/mvgvbwaMC/am,ag
      common/bwsamplecmn/aa,am2,amag,ooamag,ooaa,ifirst
      common/bwsamplecmn2/emax,emin,F1min,F2min,d1,d2,ds
      data ifirst /0/
      if (ifirst.eq.0) then
         am2    = am*am
         amag   = am*ag
         ooamag = 1.d0/amag
         aa     = am2/amag
         ooaa   = 1.d0/aa
         emin   = (q2min - am2)*ooamag
         emax   = (q2max - am2)*ooamag
         
         F1min = aa*atan(emin)
         F2min = 0.5d0 * log(emin*emin+1.d0)
         
         d1   = aa * atan(emax) - F1min
         d2   = 0.5d0 * log(emax*emax+1.d0) - F2min
         ds   = d1 + d2         
         
         ifirst = 1
      endif

      bwnormstandard = d1
      
      return
      end
      
      
      function bwnormaMC(q2max,q2min)
      implicit double precision (a-h,o-z)
      common/mvgvbwaMC/am,ag
      common/bwsamplecmn/aa,am2,amag,ooamag,ooaa,ifirst
      common/bwsamplecmn2/emax,emin,F1min,F2min,d1,d2,ds
      data ifirst /0/
      if (ifirst.eq.0) then
         am2    = am*am
         amag   = am*ag
         ooamag = 1.d0/amag
         aa     = am2/amag
         ooaa   = 1.d0/aa
         emin   = (q2min - am2)*ooamag
         emax   = (q2max - am2)*ooamag
         
         F1min = aa*atan(emin)
         F2min = 0.5d0 * log(emin*emin+1.d0)
         
         d1   = aa * atan(emax) - F1min
         d2   = 0.5d0 * log(emax*emax+1.d0) - F2min
         ds   = d1 + d2         
         
         ifirst = 1
      endif

      bwnormaMC = ds
      
      return
      end
      
      subroutine bwsampleaMC(q2max,q2min,s,w)
      implicit double precision (a-h,o-z)
      real*4 csi(1)
      common/mvgvbwaMC/am,ag
      common/bwsamplecmn/aa,am2,amag,ooamag,ooaa,ifirst
      common/bwsamplecmn2/emax,emin,F1min,F2min,d1,d2,ds
      data ifirst /0/
! from AcerMC !
      
      if (ifirst.eq.0) then
         am2    = am*am
         amag   = am*ag
         ooamag = 1.d0/amag
         aa     = am2/amag
         ooaa   = 1.d0/aa
         emin   = (q2min - am2)*ooamag
         emax   = (q2max - am2)*ooamag
         
         F1min = aa*atan(emin)
         F2min = 0.5d0 * log(emin*emin+1.d0)
         
         d1   = aa * atan(emax) - F1min
         d2   = 0.5d0 * log(emax*emax+1.d0) - F2min
         ds   = d1 + d2
         
         ifirst = 1
      endif
      
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
***************************************************************************
      subroutine get_x1x2pdf(s,w)
      include 'shared.inc'
      real*4 rnd(2),csi(1)
      common/yminimo/ymin

      common/mvgvbwaMC/amaMC,agaMC
      
      w = 1.d0
      if (i_pdf.ne.0) then
! I change the integration variables: from x1,x2 to x1,y=x1*x2 (--> x2=y/x1).
! The Jacobian of the transformation is 1/x1.
! Then I sample x1 (xmin < x1 < xmax) according to 1/x1 and y
! (x1*xmin <y < x1*xmax) according to a multi-channel BW or 1/y
         if (boson.eq.'Z') then
            am = mz
            ag = gz
            a  = 2.3d0
c            a  = 4.d0
         else
            am = mw
            ag = gw
            a  = 1.3d0
c            a  = 4.d0            
         endif

         amaMC = am/sqrt(s)         
         agaMC = ag/sqrt(s)

         ymax = xmax*xmax
         
c         anbw = 1.d0/(s*am*ag)*( atan( ymax*s/am/ag-am/ag )
c     >        - atan( ymin*s/am/ag-am/ag ) )

         anbw = bwnormaMC(ymax,ymin)
c         anbw = bwnormstandard(ymax,ymin)
         
         uma  = 1.d0 - a
         an13 = 1.d0/uma*(ymax**uma - ymin**uma)

         pbw  = 0.5d0
c         pbw  = 0.d0
         pomy = 1.d0 - pbw

*     for y, I use a multi-channel sampling, in order to flatten both BW and
*     1/y (implied by PDFs) peaks.....         
         call wraprng(rnd,2)
         if (rnd(1).lt.pbw) then             
c            y    = 1.d0*rnd(2)
c            tmp  = atan(ymin*s/am/ag-am/ag)
c            y    = ag*am/s*( tan(s*am*ag*anbw*y+tmp)+am/ag )
            call bwsampleaMC(ymax,ymin,y,wbw)
c          call bwstandard(ymax,ymin,y,wbw)
            wbw = 1.d0/wbw
         else
c            y = ymin*exp(aln*rnd(2))            
c            y = (ymax-ymin)*rnd(2) + ymin ! flat...
            y = uma * an13 * rnd(2) + ymin**uma ! as y**-a
            y = y **(1.d0/uma)
         wbw  = y/anbw/((y-amaMC*amaMC)**2+amaMC*amaMC*agaMC*agaMC)
c         wbw  = 1.d0/anbw/((y-amaMC*amaMC)**2+amaMC*amaMC*agaMC*agaMC)
         endif

c         womy = 1.d0/y/aln
c         womy = ymax-ymin  ! flat...
         womy  = 1.d0/an13/y**a

         regulator  = 1.d0/(pbw*wbw + pomy*womy)         
         w = w*regulator

c         x1min = y/xmax
         x1min = max(y/xmax,xmin)
c         x1max = xmax
         x1max = min(y/xmin,xmax)

!! conviene generare come 1/x*1/(1-x), almeno per i g indotti, mh... no
         aln = log(x1max/x1min)
         call wraprng(csi,1)
         x1  = x1min*exp(aln*csi(1))
         w   = w * aln * x1

c     x1 as 1/x/(1-x)
c         x1max = x1max - 1d-8         
c         aln = log(x1max/(1.d0 - x1max)) - log(x1min/(1.d0 - x1min))
c         call wraprng(csi,1)
c         ee = exp(aln*csi(1))*(x1min/(1.d0 - x1min))
c         x1 = ee/(1.d0 + ee)
c         w = w *aln * x1*(1.d0 - x1)
         
         yjac = 1/x1
         w    = w*yjac
         
         x1pdf = x1
         x2pdf = y/x1
*[[[[[[[[[[[[[
! simmetrizzo x1 x2, non e' cosi' necessario...
         rr = 1.d0/x1pdf
         call wraprng(csi,1)
         if (csi(1).lt.0.5) then
            x1    = x1pdf
            x1pdf = x2pdf
            x2pdf = x1
            rr = 1.d0/x2pdf
         endif
         w = w * 2.d0 /(1.d0/x1pdf+1.d0/x2pdf) * rr
*]]]]]]]]]]]]]
      else  ! if the partonic cross section is selected
         x1pdf  = 1.d0        
         x2pdf  = 1.d0
         w = 1.d0       
      endif
      return
      end
*******************      
      subroutine get_x1x2pdf_ORIG(s,w)
      include 'shared.inc'
      real*4 rnd(2),csi(1)
      common/yminimo/ymin
      w = 1.d0
      if (i_pdf.ne.0) then
! I change the integration variables: from x1,x2 to x1,y=x1*x2 (--> x2=y/x1).
! The Jacobian of the transformation is 1/x1.
! Then I sample x1 (xmin < x1 < xmax) according to 1/x1 and y
! (x1*xmin <y < x1*xmax) according to a multi-channel BW or 1/y
         if (boson.eq.'Z') then
            am = mz
            ag = gz
            a  = 2.3d0
c            a  = 4.d0
         else
            am = mw
            ag = gw
            a  = 1.3d0
         endif

         ymax = xmax*xmax

         anbw = 1.d0/(s*am*ag)*( atan( ymax*s/am/ag-am/ag )
     >        - atan( ymin*s/am/ag-am/ag ) )

         uma  = 1.d0 - a
         an13 = 1.d0/uma*(ymax**uma - ymin**uma)

         pbw  = 0.5d0
c         pbw  = 0.d0
         pomy = 1.d0 - pbw

*     for y, I use a multi-channel sampling, in order to flatten both BW and
*     1/y (implied by PDFs) peaks.....         
         call wraprng(rnd,2)
         if (rnd(1).lt.pbw) then             
            y    = 1.d0*rnd(2)
            tmp  = atan(ymin*s/am/ag-am/ag)
            y    = ag*am/s*( tan(s*am*ag*anbw*y+tmp)+am/ag )
         else
c            y = ymin*exp(aln*rnd(2))            
c            y = (ymax-ymin)*rnd(2) + ymin ! flat...
            y = uma * an13 * rnd(2) + ymin**uma ! as y**-a
            y = y **(1.d0/uma)
         endif

         wbw  = 1.d0/anbw/((y*s-am*am)**2+am*am*ag*ag)
c         womy = 1.d0/y/aln
c         womy = ymax-ymin  ! flat...
         womy  = 1.d0/an13/y**a

         regulator  = 1.d0/(pbw*wbw + pomy*womy)         
         w = w*regulator

c         x1min = y/xmax
         x1min = max(y/xmax,xmin)
c         x1max = xmax
         x1max = min(y/xmin,xmax)

         aln = log(x1max/x1min)
         call wraprng(csi,1)
         x1  = x1min*exp(aln*csi(1))
         w   = w * aln * x1

         yjac = 1/x1
         w    = w*yjac
         
         x1pdf = x1
         x2pdf = y/x1
*[[[[[[[[[[[[[
! simmetrizzo x1 x2, non e' cosi' necessario...
         rr = 1.d0/x1pdf
         call wraprng(csi,1)
         if (csi(1).lt.0.5) then
            x1    = x1pdf
            x1pdf = x2pdf
            x2pdf = x1
            rr = 1.d0/x2pdf
         endif
         w = w * 2.d0 /(1.d0/x1pdf+1.d0/x2pdf) * rr
*]]]]]]]]]]]]]
      else  ! if the partonic cross section is selected
         x1pdf  = 1.d0        
         x2pdf  = 1.d0
         w = 1.d0       
      endif
      return
      end
********************************************************************
      subroutine get_is(x1,x2,s,weight)
* This subroutine selects the initial state particles (storing the masses,
* the charges and type in the common string 'process')
      include 'shared.inc'
      common/getisfirst/igetis
      data igetis /1/
      if (igetis.eq.1) then
         do j = -ifl,ifl
            do i = -ifl,ifl
               pdfpdf(i,j) = 0.d0
               sezij(i,j)  = 0.d0
               sezij2(i,j) = 0.d0
            enddo
         enddo
         igetis = 0
      endif
      call get_pdf(x1,x2,s)

      if (abs(pdf1(0)).gt.0.d0.or.abs(pdf2(0)).gt.0.d0) then
! a bit (~10%) slower, but with photon induced proc.
        call get_is_WZ(weight)
      else
! faster, but without photon induced proc.
c         if (pdf1(0).lt.0.d0) pdf1(0) = 0.d0
c         if (pdf2(0).lt.0.d0) pdf2(0) = 0.d0
         if (boson.eq.'W') then
            call get_is_W(x1,x2,s,weight)
         else
            call get_is_Z(x1,x2,s,weight)
         endif
      endif
      return
      end
********************************************************************
      subroutine get_is_WZ(weight)
! with photon induced, both for W and Z
      include 'shared.inc'
      real*4 rnd(1)
      common/getiswzcommon/ifirst,ii(-ifl:ifl,-ifl:ifl)
      parameter (nii = (2*ifl+1)*(2*ifl+1))
      data ii,ifirst /nii*0,1/
      common/previouscharge/pch
      data pch /-10.d0/
      if (boson.eq.'W') then
         if (abs(pch-chfs).gt.0.1d0) then ! to avoid to reinitialize every time
! init for W+
            ii(0,2)  = 0
            ii(0,4)  = 0
            ii(0,-1) = 0
            ii(0,-3) = 0
            ii(0,-5) = 0
            ii(2,0)  = 0
            ii(4,0)  = 0
            ii(-1,0) = 0
            ii(-3,0) = 0
            ii(-5,0) = 0
            ii(2,-1) = 0
            ii(2,-3) = 0
            ii(2,-5) = 0
            ii(4,-1) = 0
            ii(4,-3) = 0
            ii(4,-5) = 0
            ii(-1,2) = 0
            ii(-3,2) = 0
            ii(-5,2) = 0
            ii(-1,4) = 0
            ii(-3,4) = 0
            ii(-5,4) = 0
! init for W-
            ii(0,-2) = 0
            ii(0,-4) = 0
            ii(0,1)  = 0
            ii(0,3)  = 0
            ii(0,5)  = 0
            ii(-2,0) = 0
            ii(-4,0) = 0
            ii(1,0)  = 0
            ii(3,0)  = 0
            ii(5,0)  = 0
            ii(-2,1) = 0
            ii(-2,3) = 0
            ii(-2,5) = 0
            ii(-4,1) = 0
            ii(-4,3) = 0
            ii(-4,5) = 0
            ii(1,-2) = 0
            ii(3,-2) = 0
            ii(5,-2) = 0
            ii(1,-4) = 0
            ii(3,-4) = 0
            ii(5,-4) = 0
            if (chfs.gt.0.d0) then ! W+
               ii(0,2)  = iphinduced
               ii(0,4)  = iphinduced
               ii(0,-1) = iphinduced
               ii(0,-3) = iphinduced
               ii(0,-5) = iphinduced
               ii(2,0)  = iphinduced
               ii(4,0)  = iphinduced
               ii(-1,0) = iphinduced
               ii(-3,0) = iphinduced
               ii(-5,0) = iphinduced
               ii(2,-1) = 1
               ii(2,-3) = 1
               ii(2,-5) = 1
               ii(4,-1) = 1
               ii(4,-3) = 1
               ii(4,-5) = 1
               ii(-1,2) = 1
               ii(-3,2) = 1
               ii(-5,2) = 1
               ii(-1,4) = 1
               ii(-3,4) = 1
               ii(-5,4) = 1
            else                ! W-
               ii(0,-2) = iphinduced
               ii(0,-4) = iphinduced
               ii(0,1)  = iphinduced
               ii(0,3)  = iphinduced
               ii(0,5)  = iphinduced
               ii(-2,0) = iphinduced
               ii(-4,0) = iphinduced
               ii(1,0)  = iphinduced
               ii(3,0)  = iphinduced
               ii(5,0)  = iphinduced
               ii(-2,1) = 1
               ii(-2,3) = 1
               ii(-2,5) = 1
               ii(-4,1) = 1
               ii(-4,3) = 1
               ii(-4,5) = 1
               ii(1,-2) = 1
               ii(3,-2) = 1
               ii(5,-2) = 1
               ii(1,-4) = 1
               ii(3,-4) = 1
               ii(5,-4) = 1
            endif
         endif
         pch = chfs
      else ! if boson = Z
         if (ifirst.eq.1) then
            ii(0,0)  = iphinduced
            ii(1,-1) = 1  !* 0
            ii(2,-2) = 1  !* 0
            ii(3,-3) = 1  !* 0
            ii(4,-4) = 1  !* 0
            ii(5,-5) = 1  !* 0
            ii(-1,1) = 1  !* 0
            ii(-2,2) = 1  !* 0
            ii(-3,3) = 1  !* 0
            ii(-4,4) = 1  !* 0
            ii(-5,5) = 1  !* 0
            ii(0,1)  = iphinduced  !* 0
            ii(0,2)  = iphinduced  !* 0
            ii(0,3)  = iphinduced  !* 0
            ii(0,4)  = iphinduced  !* 0
            ii(0,5)  = iphinduced  !* 0
            ii(0,-1) = iphinduced  !* 0
            ii(0,-2) = iphinduced  !* 0
            ii(0,-3) = iphinduced  !* 0
            ii(0,-4) = iphinduced  !* 0
            ii(0,-5) = iphinduced  !* 0
            ii(1,0)  = iphinduced  !* 0
            ii(2,0)  = iphinduced  !* 0
            ii(3,0)  = iphinduced  !* 0 
            ii(4,0)  = iphinduced  !* 0
            ii(5,0)  = iphinduced  !* 0
            ii(-1,0) = iphinduced  !* 0
            ii(-2,0) = iphinduced  !* 0
            ii(-3,0) = iphinduced  !* 0
            ii(-4,0) = iphinduced  !* 0
            ii(-5,0) = iphinduced  !* 0
         endif
      endif
      sumpdfpdf = 0.d0
      weight    = 0.d0
      Vud2      = 1.d0
      do j = -ifl,ifl
         do i = -ifl,ifl
            if (ii(i,j).gt.0) then
               if (boson.eq.'W') then
                  Vud2 = myckm(i,j)*myckm(i,j)
                  pdfpdf(-i,-j) = 0.d0
               endif
               if (i_pdf.eq.0) Vud2 = 1.d0
c               if (pdf1(i).lt.0.d0) pdf1(i) = 0.d0
c               if (pdf2(j).lt.0.d0) pdf2(j) = 0.d0
               pdfpdf(i,j) = pdf1(i)*pdf2(j)*Vud2
               sumpdfpdf   = sumpdfpdf + abs(pdfpdf(i,j))
            endif
         enddo
      enddo
* It can happen that the pdf product for a given x1, x2 is zero for all
* the particles. In this case, no process can occur: I keep ipart1 
* and ipart2 = ifl+1 and the routine exits
      if (sumpdfpdf.eq.0.d0) then        
         ipart1 = ifl + 1
         ipart2 = ifl + 1
         return
      endif
      weight = sumpdfpdf
* choose the process...
      i_stop = 0
      psum   = 0.d0
      call wraprng(rnd,1)
      j = -ifl-1
      do while(i_stop.eq.0.and.j.lt.ifl)
         j = j + 1
         i = -ifl - 1
         do while(i_stop.eq.0.and.i.lt.ifl)
            i = i + 1
            if (ii(i,j).gt.0) then
               psum   = psum + abs(pdfpdf(i,j))
               psum_n = psum/sumpdfpdf               
               if (psum_n.gt.rnd(1)) then
                  m1  = mq(i)
                  m2  = mq(j)
                  ch1 = chq(i)
                  ch2 = chq(j)
                  ipart1 = i
                  ipart2 = j               
                  i_stop = 1
               endif
            endif
         enddo
      enddo
      weight = weight * pdfpdf(ipart1,ipart2)/abs(pdfpdf(ipart1,ipart2))
      ifirst = 0
      return
      end
************************
      subroutine get_is_W(x1,x2,s,weight)
      include 'shared.inc'
      real*4 rnd(1)
      common/getiswfirstoriginal/ifirst
      data ifirst /1/
      if (ifirst.eq.1) then
         ifirst = 0
      endif
      sumpdfpdf = 0.d0
      psum      = 0.d0
      weight    = 0.d0
      do j = 1,ifl,2
         do i = 2, ifl,2
            if (chfs.lt.0.d0) then
               i1 =  j
               i2 = -i
            else
               i1 = -j
               i2 =  i
            endif
            Vud2 = myckm(i1,i2)**2
            if (i_pdf.eq.0) Vud2 = 1.d0

c            if (pdf1(i1).lt.0.d0) pdf1(i1) = 0.d0
c            if (pdf2(i2).lt.0.d0) pdf2(i2) = 0.d0
c            if (pdf1(i2).lt.0.d0) pdf1(i2) = 0.d0
c            if (pdf2(i1).lt.0.d0) pdf2(i1) = 0.d0

            pdfpdf(-i1,-i2) = 0.d0
            pdfpdf(-i2,-i1) = 0.d0

            pdfpdf(i1,i2) = pdf1(i1)*pdf2(i2) * Vud2
            pdfpdf(i2,i1) = pdf1(i2)*pdf2(i1) * Vud2
            sumpdfpdf = sumpdfpdf+abs(pdfpdf(i1,i2))+abs(pdfpdf(i2,i1))
         enddo
      enddo 
* It can happen that the pdf product for a given x1, x2 is zero for all
* the particles. In this case, no process can occur: I keep ipart1 
* and ipart2 = ifl+1 and the routine exits
      if (sumpdfpdf.eq.0.d0) then        
         ipart1 = ifl + 1
         ipart2 = ifl + 1
         return
      endif
      weight = sumpdfpdf
* choose the process...
      i_stop = 0     
      j= -ifl - 1
      call wraprng(rnd,1)
      do while(i_stop.eq.0.and.j.lt.ifl)
         j =  j + 1         
         i = -ifl-1
         do while(i_stop.eq.0.and.i.lt.ifl)
               i = i + 1
               if (i.ne.j.and.i.ne.0.and.j.ne.0) then
               psum   = psum + abs(pdfpdf(i,j))
               psum_n = psum/sumpdfpdf
               if (psum_n.gt.rnd(1)) then
                  m1  = mq(i)
                  m2  = mq(j)
                  ch1 = chq(i)
                  ch2 = chq(j)
                  ipart1 = i
                  ipart2 = j               
                  i_stop = 1
               endif
            endif
         enddo
      enddo
      weight = weight * pdfpdf(ipart1,ipart2)/abs(pdfpdf(ipart1,ipart2))
      return
      end
********************************************
      subroutine get_is_Z(x1,x2,s,weight)
      include 'shared.inc'
      real*4 rnd(1)
      sumpdfpdf = 0.d0
      psum      = 0.d0
      weight    = 0.d0
      do j = -ifl, ifl
         i = -j
*     no FCNC !!!!!!!!!!!!!!!	 
!         pdf1i = 
!         pdf2j = 
c         if (pdf1i.lt.0.d0) pdf1i = 0.d0
c         if (pdf2j.lt.0.d0) pdf2j = 0.d0
         pdfpdf(i,j) = pdf1(i)*pdf2(j)! pdf1i*pdf2j
         sumpdfpdf = sumpdfpdf + abs(pdfpdf(i,j))
      enddo 
* It can happen that the pdf product for a given x1, x2 is zero for all
* the particles. In this case, no process can occur: I keep ipart1 
* and ipart2 = ifl+1 and the routine exits
      if (sumpdfpdf.eq.0.d0) then        
         ipart1 = ifl+1
         ipart2 = ifl+1
         return
      endif
      weight = sumpdfpdf
* choose the process...
      i_stop = 0
      j=-ifl-1
      call wraprng(rnd,1)
      do while(i_stop.eq.0.and.j.lt.ifl)
         j =  j + 1
         i = -j
         psum   = psum + abs(pdfpdf(i,j))
         psum_n = psum/sumpdfpdf

         if (psum_n.gt.rnd(1)) then ! i,j correspond to 
                                ! the occurring
                                ! process,  when 
                                ! the do-loop stops
            m1  = mq(i)
            m2  = mq(j)
            ch1 = chq(i)
            ch2 = chq(j)               
            ipart1 = i
            ipart2 = j               
            i_stop = 1
         endif
      enddo
      weight = weight * pdfpdf(ipart1,ipart2)/abs(pdfpdf(ipart1,ipart2))
      return
      end
***********************************************************
      subroutine get_pdf(x1,x2,s)
      include 'shared.inc'
      parameter (iiiii=2*ifl + 1)
      data pdf1 /iiiii*0.d0/
      data pdf2 /iiiii*0.d0/
      real*4 csi(1)

      if (i_pdf.ne.0) then   ! if the PDF cross section is selected....
* PDF called for lepton 1

         call pdf_wrapper(x1,s,electron,positron,photon)
         if (hadr1.eq.2212) then	   ! if the lepton is an electron....
         
              
            pdf1(1)  = electron
c            pdf1(2)  = muon  
c            pdf1(3)  = tau      
            pdf1(-1) = positron
c            pdf1(-2) = amuon
c            pdf1(-3) = atau

            pdf1(0)  = photon
         else ! else, if it is a positron....
            pdf1(-1) = electron
c            pdf1(-2) = muon
c            pdf1(-3) = tau
            pdf1(1)  = positron
c            pdf1(2)  = amuon    
c            pdf1(3)  = atau       

            pdf1(0)  = photon
         endif
c         do i= -ifl, ifl ! checked after!!
c            if (pdf1(i).lt.0.d0) pdf1(i)=0.d0
c         enddo    
* PDF called for lepton 2

         call pdf_wrapper(x2,s,electron,positron,photon)    
         if (hadr2.eq.2212) then  ! if the lepton is an electron....
            pdf2(1)  = electron
c            pdf2(2)  = muon 
c            pdf2(3)  = tau       
            pdf2(-1) = positron
c            pdf2(-2) = amuon
c            pdf2(-3) = atau
            
            pdf2(0)  = photon
         else ! else, if it is a positron....
            pdf2(-1) = electron
c            pdf2(-2) = muon
c            pdf2(-3) = tau
            pdf2(1)  = positron
c            pdf2(2)  = amuon  
c            pdf2(3)  = atau         

            pdf2(0)  = photon
         endif     
         
c         do i = -ifl, ifl ! checked after!!
c            if (pdf2(i).lt.0.d0) pdf2(i)=0.d0
c         enddo    
      else  ! if the partonic cross section is selected....
         pdf1(iquark1)  = 1.d0
         pdf2(iquark2)  = 1.d0
      endif
      return
      end
