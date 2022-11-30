! begin HORACE
	double complex function B00res(p,mw,gw)
	implicit none
	double precision p,mw,gw,mw2
	double complex gwmw, B0, B0res,ff2,ff1
	integer ier
	
	external B0,B0res

	gwmw = dcmplx(0d0,gw*mw)
!	gwmw = dcmplx(0d0,1.d-8)
	mw2=mw*mw

	call ljffxb0(ff1,0.d0,mw2,mw2,ier)
	call ljffxb0(ff2,0.d0,0.d0,mw2,ier)

	B00res = 1d0/3d0*(
     -           mw2/2d0*ff1+
     -           (p-mw2)/4d0/p*((-p+mw2)*B0res(p,mw,gw)-
     -                           mw2*ff2)) +
     -           mw2/3d0-p/18d0
	end function B00res
!***********************************************************************
	double complex function B1res(p,mw,gw)
	implicit none
	double precision p,mw,gw,mw2
	double complex gwmw,B0,B0res,ff1
	integer ier
	
	external B0,B0res

	gwmw = dcmplx(0d0,gw*mw)
	mw2=mw*mw

	call ljffxb0(ff1,0.d0,mw2,mw2,ier)

	B1res = 0.5d0/p*((-p+mw2)*B0res(p,mw,gw)-
     -                     mw2*ff1-mw2)
	end function B1res
!***********************************************************************
	double complex function B1resWA(p,mw,gw)
	implicit none
	double precision p,mw,gw,mw2
	double complex gwmw,B0,B0res,ff1
	integer ier

	external B0,B0res

	gwmw = dcmplx(0d0,gw*mw)
	mw2=mw*mw

	call ljffxb0(ff1,0.d0,mw2,mw2,ier)
	
	B1resWA = 0.5d0/p*((-p-mw2)*B0res(p,mw,gw)+
     -                     mw2*ff1+mw2)
	end function B1resWA
!***********************************************************************
	double complex function B11res(p,mw,gw)
	implicit none
	double precision p,mw,gw,mw2
	double complex gwmw, B0, B0res,B11, ff1,ff2
	integer ier

	external B0,B0res,B11
	
	gwmw = dcmplx(0d0,gw*mw)
	mw2=mw*mw

cccc this is valid for a B11(p,0,mw2)
cccc  it is NOT valid for a B11(p,mw2,0)

	call ljffxb0(ff1,0.d0,mw2,mw2,ier)
	call ljffxb0(ff2,0.d0,0.d0,mw2,ier)

	
	B11res = 1d0/3d0/p*(
     -           mw2*ff1 +
     -           (p-mw2)*B0res(p,mw,gw) +
     -           p/6d0+mw2/2d0 -
     -           mw2/p*(p-mw2)*(B0res(p,mw,gw)-ff2))
	end function B11res
!***********************************************************************

!***********************************************************************
	double complex function B11resWA(p,mw,gw)
	implicit none
	double precision p,mw,gw,mw2
	double complex gwmw, B0, B0res,B11, ff1
	integer ier

	external B0,B0res,B11
	
	gwmw = dcmplx(0d0,gw*mw)
	mw2=mw*mw

cccc this is valid for a B11(p,mw2,0)
cccc  it is NOT valid for a B11(p,0,mw2)

	call ljffxb0(ff1,0.d0,0.d0,mw2,ier)

	B11resWA = 1d0/3d0/p/p*(
     -           p*p*B0res(p,mw,gw) +
     -           p*p/6d0-mw2*p/2d0 +
     -           mw2*(p+mw2)*(B0res(p,mw,gw)-ff1))
	

	end function B11resWA
!***********************************************************************
	double complex function B0res(p,mw,gw)
	implicit none
	double precision gw,mw,p,mw2
	double precision getmudim,getdelta
	double complex gwmw,b0res2

	external getmudim,getdelta
	gwmw = dcmplx(0d0,gw*mw)
!	gwmw = dcmplx(0d0,1.d-8)
	mw2=mw*mw
c	print*,'B0res---',getdelta(),getmudim(),p

	B0res = getdelta()+2+log(getmudim()/mw2)+
     -          (mw2/p-1d0)*log((mw2-p-gwmw)/mw2)
        end function B0res
!***********************************************************************
	double complex function C0resonant(p12,p1p2,p22,m02,m12,MW,GW)
	double precision p12,p1p2,p22,m02,m12,pi2
	double precision MW,GW
	double complex cli2,dummy

	integer ier,ipi12
c	parameter (pi2=3.1415926**2)
	parameter (pi2=16.d0*atan(1.d0)**2)
	
	call ljffzxdl(cli2,ipi12,dummy,(mw**2-p22)/mw**2,1,ier)
	
	C0resonant = 1d0/p22*
     -                    (log(p22/m12)*log((MW**2-dcmplx(0d0,GW*MW)-
     -                                       p22)/MW**2)+
     -                     dble(cli2)+ipi12*pi2/12d0-pi2/6d0)

	end function C0resonant
!***********************************************************************
	function D0_irHOR(cp1, cp2, cp3, cp4, cp1p2, cp2p3, m3c)
	implicit none
        double complex  D0_irHOR
        double complex  cp1, cp2, cp3, cp4, cp1p2, cp2p3,gzmz
        double precision p1, p2, p3, p4, p1p2, p2p3, m3
 	double precision d, delta
        double complex  m3c, cm3, m1_, m4_, ieps
        double complex  xs, x2, x3, y, c, f
        double complex  logxs, logx2, logx3, log1x2, log1x3, logy
        double precision mgamma2, pi, getlambda

        double complex  myln, spence, bdK, addeps, mycln
        external myln, spence, bdK, addeps, mycln

	character*1 boson
	common/vectorboson/boson
	integer ifirst,k
	common/d0irifirst/pi,ieps,ifirst
	data ifirst /1/
	if (ifirst.eq.1) then
	   do k = 1,50
c	      print*,'ATTENTION IN THE DO_ir !!'
	   enddo
           pi = 3.1415926535897932384626433832795029d0
           ieps = dcmplx(0d0,1d-20)
	   ifirst = 0
	endif

	mgamma2 = getlambda()
        delta = dble(mgamma2)
	m3 = dble(m3c)
        cm3 = m3**2

	gzmz = dcmplx(0.d0,-2.d0*m3*imag(m3c))
        
        p1 = dble(cp1)
        p2 = dble(cp2)
        p3 = dble(cp3)
        p4 = dble(cp4)
        p1p2 = dble(cp1p2)
        p2p3 = dble(cp2p3)
        
        m1_ = sqrt(p1)
        m4_ = sqrt(p4)

        d = p2p3 - (m1_ - m4_)**2

        f = .5D0/m1_/m4_/(p1p2 - cm3)
        
        if(abs(d).gt.0d0) then
          xs = bdK(p2p3, m1_, m4_)
          logxs = log(xs)
          f = f*2D0*xs/(1D0 - xs**2)
        endif
        
        ! massless case
        if (m3.eq.(0.D0,0.D0)) then
        !write(*,*) ' massless'
          if(p1.eq.p2 .and. p3.eq.p4) then
            D0_irHOR = 2D0*f*myln(-delta/p1p2, 1D0)
            if(abs(d).gt.0d0) D0_irHOR = -logxs*D0_irHOR
            return
          endif
        
          if(p1.eq.p2) then
        !write(*,*) 'p1=p2',1D0/(p2p3-p4)/p1p2
        !write(*,*) 'p1=p2',-m1_*m4_/(p2p3-p4)
        !write(*,*) 'p1=p2',mycln(-m1_*m4_/(p2p3-p4),1D0)
        !write(*,*) 'p1=p2',2D0*myln((p3-p4)/p1p2,1D0)
            D0_irHOR = 1D0/(p2p3-p4)/p1p2*
     +                 (mycln(-m1_*m4_/(p2p3-p4),1D0)*
     +                     (mycln(-m1_*m4_/(p2p3-p4),1D0)+log(delta/p4)+
     +                      2D0*myln((p3-p4)/p1p2,1D0))+
     +                                       pi**2/6D0)
	    print*,'p1eqp2'
            return
          endif
        
          if(p3.eq.p4) then
            print*, 'p3=p4'
            D0_irHOR = 1D0/(p2p3-p1)/p1p2*
     +                 (mycln(-m4_*m1_/(p2p3-p1),1D0)*
     +                     (mycln(-m4_*m1_/(p2p3-p1),1D0)+log(delta/p1)+
     +                      2D0*myln((p2-p1)/p1p2,1D0))+
     +                                       pi**2/6D0)
            return
          endif
        
          y = m1_/m4_*(p3 - p4 + IEPS)/ (p2 - p1 + IEPS)
          logy = log(y)
          c = mycln(delta/m1_/m4_, 0D0) +
     +       myln((p2 - p1)/p1p2, p1 - p2) + myln((p3 - p4)/p1p2, p4 - p3)
        
          if(abs(d).gt.0d0) then
            D0_irHOR = f*(pi**2/6D0 +logxs*(-.5D0*logxs+2D0*log(1D0-xs**2)-c)+
     +         spence(xs**2, 0D0) + .5D0*logy**2 -spence(xs/y, 0D0) - 
     +         (logxs + log(1D0/y))*log(1D0 - xs/y) -spence(xs*y, 0D0) - 
     +         (logxs + logy)*log(1D0 - xs*y))
          print*, ' massless',D0_irHOR
            return
          endif
          D0_irHOR = f*(c - 2D0 - (1D0 + y)/(1D0 - y)*logy)
          return
        endif
        
        x2 = bdK(p2, m1_, dcmplx(m3))
        x3 = bdK(p3, m4_, dcmplx(m3))
        logx2 = log(x2)
        logx3 = log(x3)
        log1x3 = log(1D0/x3)
c        c = mycln(m3*sqrt(delta)/(m3c**2 - p1p2), 1D0)
        c = log(m3*sqrt(delta)/(m3c**2 - p1p2))

        if(abs(d).gt.0d0) then
          log1x2 = log(1D0/x2)
          D0_irHOR = f*(.5D0*pi**2 +2D0*log(xs)*(log(1D0 - xs**2) - c) +
     +       spence(xs**2, 0D0)+logx2**2 + logx3**2 -spence(xs/x2/x3, 0D0) - 
     +       (logxs+log1x2+log1x3)*log(1D0-xs/x2/x3) -spence(xs*x2/x3, 0D0) - 
     +       (logxs + logx2 +log1x3)*log(1D0-xs*x2/x3)-spence(xs/x2*x3, 0D0) - 
     +       (logxs+log1x2+logx3)*log(1D0 - xs/x2*x3) -spence(xs*x2*x3, 0D0) - 
     +       (logxs + logx2 + logx3)*log(1D0 - xs*x2*x3))
          return
        endif

        D0_irHOR = f*(2D0*c - (1D0 + x2/x3)/(1D0 - x2/x3)*(logx2 + log1x3) - 
     +     (1D0 + x2*x3)/(1D0 - x2*x3)*(logx2 + logx3) - 2D0)

        end function D0_irHOR
!**********************************************************************
        function addeps(k)
        implicit none
        double complex addeps
        double complex k
        addeps = k*dcmplx(1D0, -sign(1d-8, dreal(k)))
        end function addeps        
!*********************************************************************
        subroutine roots(p, m1, m2, x1, x2, y1, y2, r)
        implicit none
        double precision p, m1, m2
        double complex x1, x2, y1, y2, r
        double precision q,eps

	eps = 1d-20

        r = sqrt(dcmplx(p*(p - 2D0*(m1 + m2)) + (m1 - m2)**2))
        q = p + m1 - m2
        x1 = (q + r)/2D0/p
        x2 = (q - r)/2D0/p
        if(abs(x2) > abs(x1)) then
          x1 = m1/p/x2
        else if(abs(x1) > abs(x2)) then
          x2 = m1/p/x1
        endif
        x1 = x1 + dcmplx(0D0,  abs(p*x1)/p*EPS)
        x2 = x2 + dcmplx(0D0, -abs(p*x2)/p*EPS)
        q = p - m1 + m2
        y2 = (q + r)/2D0/p
        y1 = (q - r)/2D0/p
        if(abs(y2) > abs(y1)) then
          y1 = m2/p/y2
        else if(abs(y1) > abs(y2)) then
          y2 = m2/p/y1
        endif
        y1 = y1 + dcmplx(0D0, -abs(p*y1)/p*EPS)
        y2 = y2 + dcmplx(0D0,  abs(p*y2)/p*EPS)
        end subroutine roots
!***********************************************************************
        function myln(x, isig)
        implicit none
        double complex myln
        double precision x, isig, pi

	pi = 3.1415926535897932384626433832795029d0

        if(x > 0D0) then
          myln = log(x)
        else
          myln = log(-x) +dcmplx(0D0, sign(pi, isig))
        endif
        end function myln
!*********************************************************************
        function mycln(z, isig)
        implicit none
        double complex mycln
        double complex z
        double precision isig, pi

	pi = 3.1415926535897932384626433832795029d0

        if(dimag(z).eq.0D0 .and. dble(z).lt.0D0) then
#ifdef WARNINGS
          if(isig.eq.0D0) print *, "mycln: argument on cut"
#endif
          mycln = log(-z) +dcmplx(0D0, sign(pi, isig))
        else
          mycln = log(z)
        endif
        end function mycln
!*********************************************************************
        function spence(z, isig)
        implicit none
        double complex spence
        double complex z
        double precision isig,pi
        double complex z1, ieps
        double precision az1
        double complex li2series, mycln
        external li2series, mycln
	parameter (pi = 3.1415926535897932384626433832795029d0)
	parameter (ieps = dcmplx(0d0,1d-20))

        z1 = 1D0 - z
        az1 = abs(z1)
#ifdef WARNINGS
        if(isig.eq.0D0 .and.dimag(z).eq.0D0 .and. abs(dble(z1)) < 
     +   CALACC)   print *, "spence: argument on cut"
#endif
        if(az1 < 1D-15) then
          spence = pi**2/6D0
        else if(dble(z) < .5D0) then
          if(abs(z) < 1D0) then
            spence = li2series(z, isig)
          else
            spence = -pi**2/6D0 -.5D0*mycln(-z, -isig)**2-li2series(1D0/z,-isig)
          endif
        else
          if(az1 < 1D0) then
            spence = pi**2/6D0 -mycln(z, isig)*mycln(z1, -isig)-li2series(z1,-isig)
          else
            spence = pi**2/3D0 +
     +         .5D0*mycln(-z1, -isig)**2 - mycln(z, isig)*mycln(z1, -isig) +
     +         li2series(1D0/z1, isig)
          endif
        endif
        end function spence
!*********************************************************************
        function li2series(z, isig)
        implicit none
        double complex li2series
        double complex z
        double precision isig
        double complex xm, x2, new
        integer j
        double complex mycln
        external mycln
        ! these are the even-n Bernoulli numbers, already divided by (n + 1)!
        ! as in Table[BernoulliB[n]/(n + 1)!, {n, 2, 50, 2}]
        double precision b(25)
        data b / 0.02777777777777777777777777777777777777777778774D0, 
     &    -0.000277777777777777777777777777777777777777777778D0, 
     &    4.72411186696900982615268329554043839758125472D-6, 
     &    -9.18577307466196355085243974132863021751910641D-8, 
     &    1.89788699889709990720091730192740293750394761D-9,  
     &    -4.06476164514422552680590938629196667454705711D-11,  
     &    8.92169102045645255521798731675274885151428361D-13,  
     &    -1.993929586072107568723644347793789705630694749D-14,  
     &    4.51898002961991819165047655285559322839681901D-16,  
     &    -1.035651761218124701448341154221865666596091238D-17,  
     &    2.39521862102618674574028374300098038167894899D-19,  
     &    -5.58178587432500933628307450562541990556705462D-21,  
     &    1.309150755418321285812307399186592301749849833D-22,  
     &    -3.087419802426740293242279764866462431595565203D-24,  
     &    7.31597565270220342035790560925214859103339899D-26,  
     &    -1.740845657234000740989055147759702545340841422D-27,  
     &    4.15763564461389971961789962077522667348825413D-29,  
     &    -9.96214848828462210319400670245583884985485196D-31,  
     &    2.394034424896165300521167987893749562934279156D-32,  
     &    -5.76834735536739008429179316187765424407233225D-34,  
     &    1.393179479647007977827886603911548331732410612D-35,  
     &    -3.372121965485089470468473635254930958979742891D-37,  
     &    8.17820877756210262176477721487283426787618937D-39,  
     &    -1.987010831152385925564820669234786567541858996D-40,  
     &    4.83577851804055089628705937311537820769430091D-42 /
        xm = -mycln(1D0 - z, -isig)
        x2 = xm**2
        li2series = xm - x2/4D0
        do j = 1, 25
          xm = xm*x2
          new = li2series + xm*b(j)
          if(new == li2series) return
          li2series = new
        enddo
#ifdef WARNINGS
        print *, "li2series: bad convergence"
#endif
        end function li2series
!*********************************************************************
        function contspence(z1, z2)
        implicit none
        double complex contspence, z1, z2, z12, ris
        double complex spence
        external  spence
        
        z12 = z1*z2
        !ris = spence(1D0-z12,0D0)+log(1D0-z12)*(log(z12)-Log(z1)-log(z2))
        !write(*,*) 'contspence=',ris
        
        contspence=spence(1D0-z12,0D0)+log(1D0-z12)*(log(z12)-Log(z1)-log(z2))
        
        end function contspence
!**********************************************************************
        function cspence(z1, z2, im1, im2)
        implicit none
        double complex cspence
        double complex z1, z2
        double precision im1, im2, pi
        double complex mycln, spence, ieps
        integer eta
        external mycln, spence, eta
        double complex z12
        double precision im12
        integer etas

	pi = 3.1415926535897932384626433832795029d0
	ieps = dcmplx(0d0,1d-20)

        
        if (cdabs(z1)<1D-15.or.cdabs(z2)<1D-15) then
          cspence = pi**2/6D0
          return
        end if
        z12 = z1*z2
        im12 = im2*sign(1D0, dble(z1))
        if(dble(z12) > .5D0) then
          cspence = spence(1D0 - z12, 0D0)
          etas = eta(z1, z2, im1, im2, im12)
          if(abs(etas).gt.0) cspence = cspence +etas*dcmplx(0D0, 2D0*pi)*
     +       mycln(1D0 - z12, -im12)
        else if(abs(z12) < 1D-4) then
          cspence = pi**2/6D0 -spence(z12, 0D0) + (mycln(z1, im1) + 
     +      mycln(z2, im2))*z12* 
     +       (1D0 + z12*(.5D0 + z12*(1D0/3D0 + z12/4D0)))
        else
          cspence = pi**2/6D0 -spence(z12, 0D0) - 
     +       (mycln(z1, im1) + mycln(z2, im2))*mycln(1D0 - z12, 0D0)
        endif
        end function cspence
!*********************************************************************
        function eta(c1, c2, im1x, im2x, im12x)
        implicit none
        integer eta
        double complex c1, c2
        double precision im1x, im2x, im12x
        double precision im1, im2, im12
        im1 = dimag(c1)
        if(im1 .eq. 0D0) im1 = im1x
        im2 = dimag(c2)
        if(im2 .eq. 0D0) im2 = im2x
        im12 = dimag(c1*c2)
        if(im12 .eq. 0D0) im12 = im12x
        if(im1 < 0D0 .and. im2 < 0D0 .and. im12 > 0D0) then
          eta = 1
        else if(im1 > 0D0 .and. im2 > 0D0 .and. im12 < 0D0) then
            eta = -1
          else
            eta = 0
#ifdef WARNINGS
            if(.not. (im2 .eq. 0D0 .and. dble(c2) > 0D0 .or.
     +          im1 .eq. 0D0 .and. dble(c1) > 0D0) .and. 
     +         (im1 .eq. 0D0 .and. dble(c1) < 0D0 .or.
     +         im2 .eq. 0D0 .and. dble(c2) < 0D0 .or.
     +         im12 .eq. 0D0 .and. dble(c1*c2) < 0D0)
     +         ) print *, 'eta not defined'
#endif
          endif
        end function eta
!********************************************************************
        function eta_n(c1, c2)
        implicit none
        integer eta_n
        double complex c1, c2
        integer eta
        external eta
        eta_n = eta(c1, c2, 0D0, 0D0, 0D0)
        end function eta_n
!*********************************************************************
        function eta_tilde(c1, c2, im1x, im2x)
       	implicit none
        integer eta_tilde
        double complex c1, c2
        double precision im1x, im2x
        double precision im1, im2
        integer eta
        external eta
        im1 = dimag(c1)
        if(im1 .eq. 0D0) im1 = im1x
        im2 = dimag(c2)
        if(abs(im2).gt.0D0) then
          eta_tilde = eta(c1, c2, im1x, 0D0, 0D0)
        else if(dble(c2) > 0D0) then
          eta_tilde = 0
        else if(im1 > 0D0 .and. im2x > 0D0) then
          eta_tilde = -1
        else if(im1 < 0D0 .and. im2x < 0D0) then
          eta_tilde = 1
        else
          eta_tilde = 0
#ifdef WARNINGS
          if(im1 .eq. 0D0 .and. dble(c1) < 0D0 .or.
     +       im2x .eq. 0D0 .and. dble(c1*c2) < 0D0
     +      ) print *, 'eta_tilde not defined'
#endif
        endif
        end function eta_tilde
!**********************************************************************
        function bdK(x, m1, m2)
! this is actually -K from the Beenakker/Denner paper for D0_ir
        implicit none
        double complex bdK
        double complex m1, m2, d
        double precision x
        double complex zz,ieps
	parameter (ieps = dcmplx(0d0,1d-20))
        d = x - (m1 - m2)**2
        !write(*,*) 'bdk=',x,m1,m2,d
        if(d .eq. dcmplx(0D0,0D0)) then
          bdK = dcmplx(1D0,0D0)
        else
          zz = sqrt(1D0 - 4D0*m1*m2/(d + IEPS))
          bdK = (zz - 1D0)/(zz + 1D0)
        endif        
        end function bdK
! end HORACE
