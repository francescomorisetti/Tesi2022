      subroutine initPDFhorace
      include 'shared.inc'
      character*50  pdfname
      integer pdfname_length
      data pdfname
     >     /'                                                  '/
      
      common/parametersforPDF/iset
      CHARACTER prefix*50
      common/mstwcmn/prefix
      common/yminimo/ymin       ! for x1,x2pdf sampling, changed in cuts!

      logical  has_photon,hphot
      external has_photon
      common/hasthephoton/hphot
c
      commentonpdfs = 'no comment (on PDFs)' ! max 50 characters

*     setting PDFs parameters.... 
      if (iwhich.ge.1) then
         qpdfmin = 1.25d0
         qpdfmax = 1.d10
         xpdfmax = 1.d0
         xpdfmin = 1.d-6
         
!!            LHAPDF interface
         if (iwhich.eq.1) pdfname="NLL_DELTA_MSBAR"
         if (iwhich.eq.2) pdfname="LL_ALGMU"
         
	print*,pdfname
	
        call initfromgrid_name(pdfname)

	print*,'ciao'
         commentonpdfs = pdfname(1:40)//' (eMELA)'
c         call initPDF(ipdfreplica)
c         hphot = has_photon()

      endif
***
      ymin = xpdfmin*xpdfmin
      return
      end
*****************************************
c      subroutine pdf_wrapper(x,s,up,down,aup,adown,str,ch,bot,glu,phot)
       subroutine pdf_wrapper(x,s,el,pos,phot)
      
* This subroutine returns the pdf's for the ---> ELECTRON <---
* inputs are:
* x   -  momentum fraction
* s   -  pdf scale (GeV)
* outputs are:

      include 'shared.inc'
      double precision f(-1:1)
      common/parametersforPDF/iset
      CHARACTER prefix*50
      common/mstwcmn/prefix
      common/pdfwrappercmn/ifirst

      data ifirst /0/
      logical hphot
      common/hasthephoton/hphot

      if (ifirst.eq.0) then
         ifirst = 1
      endif
      
c      if (iwhich.eq.1.or.iwhich.eq.2) then
c      if (hphot) then
c         call evolvePDFphoton(x,s,f,xphot)
      s=100.d0
      call elpdfq2(0,11,x,1.d0-x,s**2,1.d0,f(1))
      call elpdfq2(0,-11,x,1.d0-x,s**2,1.d0,f(-1))
      call elpdfq2(0,22,x,1.d0-x,s**2,1.d0,f(0))
      xphot = 0.d0
      
      el  =  f(1)/x
      
      pos =  f(-1)/x
      
      phot   =  f(0)/x
      
      
c      phot  =  xphot/x * iphinduced
      return
      end
   
*************************************************************      
      subroutine printdist
      
     
      double precision f(-1:1)
      double precision X 
      double precision Q

      
      open (unit=1, file='elpdf.dat')
      open (unit=2, file='pospdf.dat')
      open (unit=3, file='phpdf.dat')
      
      X=1.d-3
      Q=100.d0
      
      do i=1,999
         
          
         
         call elpdfq2(0,11,X,1.d0-X,Q**2,1.d0,f(1))
         write(1,*) X, f(1)
         
         call elpdfq2(0,-11,X,1.d0-X,Q**2,1.d0,f(-1))
         write(2,*) X, f(-1)
         
         call elpdfq2(0,22,X,1.d0-X,Q**2,1.d0,f(0))
         write(3,*) X, f(0)
         
         X=X+1.0d-3
         
         
      end do
      
      close(1)
      close(2)
      close(3)
      
      return
      end
      
c      endif

cc OLD STUFF
c$$$      if (iwhich.eq.-5) then
c$$$* standalone mstw2008
c$$$c         print*,has_photon()
c$$$         
c$$$c         call evolvePDF(x,s,f)
c$$$
c$$$         CALL GetAllPDFsAlt(prefix,iset,x,s,f,xphoton)
c$$$
c$$$         
c$$$         up    =  f(2)/x
c$$$         down  =  f(1)/x
c$$$         aup   =  f(-2)/x
c$$$         adown =  f(-1)/x
c$$$         str   =  f(3)/x
c$$$         ch    =  f(4)/x
c$$$         bot   =  f(5)/x
c$$$         astr  =  f(-3)/x
c$$$         ach   =  f(-4)/x
c$$$         abot  =  f(-5)/x
c$$$         glu   =  f(0)/x
c$$$         phot  =  xphot/x * iphinduced
c$$$         return
c$$$      endif
c$$$      if (iwhich.eq.1) then
c$$$! updated 29/08/2005 for mrst 2004 qed parameterization
c$$$*  iwhich = 1 --> MRST QED parametrization....
c$$$! MRST subroutines return x * PDF and NOT PDF
c$$$!         call mrst2001(x,s,imrst,uv,dv,us,ds,stra,cha,bottom,gluon)
c$$$         call mrstqed(x,s,iset,uv,dv,us,ds,stra,cha,bottom,gluon,
c$$$     .        photon)
c$$$         up    = (uv+us)/x
c$$$         down  = (dv+ds)/x
c$$$         aup   = us/x
c$$$         adown = ds/x
c$$$         str   = stra/x
c$$$         ch    = cha/x
c$$$         bot   = bottom/x
c$$$         astr  = stra/x
c$$$         ach   = cha/x
c$$$         abot  = bottom/x
c$$$         glu   = gluon/x
c$$$         phot  = photon/x * iphinduced
c$$$         return
c$$$      endif
c$$$      if (iwhich.eq.2) then 
c$$$*     iwhich = 2 --> CTEQ6 parametrization....
c$$$         up    =  Ctq6Pdf(1,x,s)
c$$$         down  =  Ctq6Pdf(2,x,s)
c$$$         aup   =  Ctq6Pdf(-1,x,s)
c$$$         adown =  Ctq6Pdf(-2,x,s)
c$$$         str   =  Ctq6Pdf(3,x,s)
c$$$         ch    =  Ctq6Pdf(4,x,s)
c$$$         bot   =  Ctq6Pdf(5,x,s)
c$$$         astr  =  Ctq6Pdf(-3,x,s)
c$$$         ach   =  Ctq6Pdf(-4,x,s)
c$$$         abot  =  Ctq6Pdf(-5,x,s)
c$$$         glu   =  Ctq6Pdf(0,x,s)
c$$$         phot  = 0.d0 * iphinduced
c$$$         return
c$$$      endif

     
