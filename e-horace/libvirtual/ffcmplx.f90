!	ffcheck.F
!		this is an adaptation of bcanew.f by Ansgar Denner
!		to the conventions of LoopTools
!		last modified 19 Jul 99 th
! functions visible from the outside are
!   A0_check,
!   B0_check, DB0_check,
!   B1_check, DB1_check,
!   B00_check, DB00_check,
!   B11_check, DB11_check,
!   C0_check,
!   D0_check
! and the subroutines for testing:
!   check_A, check_B, check_C, and check_D
!#define WARNINGS
#define ROOT_CALACC 1D-6
#define CALACC 1D-12
#define EPS 1D-20
#define ONEpEPS dcmplx(1D0, EPS)
#define ONEmEPS dcmplx(1D0, -EPS)
#define IEPS dcmplx(0D0, EPS)
!#define pi 3.14159265358979323846264338328D0
!***********************************************************************
subroutine roots(p, m1, m2, x1, x2, y1, y2, r)
use global
implicit none
real*16 :: p, m1, m2
double complex :: x1, x2, y1, y2, r
real*16 :: q
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
function fpv(n, x, y)
use global
implicit none
double complex :: fpv
integer :: n
double complex :: x, y
integer :: m
double complex :: xm
if(abs(x) < 10D0) then
  if(n == 0) then
    fpv = -log(-y/x)
  else if(x == dcmplx(0D0)) then
    fpv = -1D0/n
  else
    fpv = 0D0
    xm = 1D0
    do m = 0, n - 1
      fpv = fpv - xm/(n - m)
      xm = xm*x
    enddo
    fpv = fpv - xm*log(-y/x)
  endif
else
  fpv = 0D0
  xm = 1D0
  do m = 1, 30
    xm = xm/x
    fpv = fpv + xm/(m + n)
    if(abs(xm/fpv) < CALACC**2) return
  enddo
endif
end function fpv

!***********************************************************************
function yfpv(n, x, y)
use global
implicit none
double complex :: yfpv
integer :: n
double complex :: x, y
double complex :: fpv
external fpv
if(abs(y) == 0D0) then
  yfpv = 0D0
else
  yfpv = y*fpv(n, x, y)
endif
end function yfpv

!***********************************************************************
function xlogx(x)
use global
implicit none
double complex :: xlogx
double complex :: x
if(abs(x) == 0D0) then
  xlogx = 0D0
else
  xlogx = x*log(x)
endif
end function xlogx

!***********************************************************************
function ln(x, isig)
use global
implicit none
double complex :: ln
real*16 :: x, isig
if(x > 0D0) then
  ln = log(x)
else
  ln = log(-x) +dcmplx(0D0, sign(pi, isig))
endif
end function ln

!***********************************************************************
function cln(z, isig)
use global
implicit none
double complex :: cln
double complex :: z
real*16 :: isig
if(dimag(z) == 0D0 .and. dble(z) <= 0D0) then
#ifdef WARNINGS
  if(isig == 0D0) print *, "cln: argument on cut"
#endif
  cln = log(-z) +dcmplx(0D0, sign(pi, isig))
else
  cln = log(z)
endif
end function cln

!***********************************************************************
function spence(z, isig)
use global
implicit none
double complex :: spence
double complex :: z
real*16 :: isig
double complex :: z1
real*16 :: az1
double complex :: li2series, cln
external li2series, cln
z1 = 1D0 - z
az1 = abs(z1)
#ifdef WARNINGS
if(isig == 0D0 .and.dimag(z) == 0D0 .and. abs(dble(z1)) < CALACC) &
  print *, "spence: argument on cut"
#endif
if(az1 < 1D-15) then
  spence = pi**2/6D0
else if(dble(z) < .5D0) then
  if(abs(z) < 1D0) then
    spence = li2series(z, isig)
  else
    spence = -pi**2/6D0 -.5D0*cln(-z, -isig)**2 - li2series(1D0/z, -isig)
  endif
else
  if(az1 < 1D0) then
    spence = pi**2/6D0 -cln(z, isig)*cln(z1, -isig) - li2series(z1, -isig)
  else
    spence = pi**2/3D0 +&
      &.5D0*cln(-z1, -isig)**2 - cln(z, isig)*cln(z1, -isig) +&
      &li2series(1D0/z1, isig)
  endif
endif
end function spence

!***********************************************************************
function li2series(z, isig)
use global
implicit none
double complex :: li2series
double complex :: z
real*16 :: isig
double complex :: xm, x2, new
integer :: j
double complex :: cln
external cln
! these are the even-n Bernoulli numbers, already divided by (n + 1)!
! as in Table[BernoulliB[n]/(n + 1)!, {n, 2, 50, 2}]
real*16 :: b(25)
data b / 0.02777777777777777777777777777777777777777778774D0, &
  -0.000277777777777777777777777777777777777777777778D0, &
  4.72411186696900982615268329554043839758125472D-6, &
  -9.18577307466196355085243974132863021751910641D-8, &
  1.89788699889709990720091730192740293750394761D-9, &
  -4.06476164514422552680590938629196667454705711D-11, &
  8.92169102045645255521798731675274885151428361D-13, &
  -1.993929586072107568723644347793789705630694749D-14, &
  4.51898002961991819165047655285559322839681901D-16, &
  -1.035651761218124701448341154221865666596091238D-17, &
  2.39521862102618674574028374300098038167894899D-19, &
  -5.58178587432500933628307450562541990556705462D-21, &
  1.309150755418321285812307399186592301749849833D-22, &
  -3.087419802426740293242279764866462431595565203D-24, &
  7.31597565270220342035790560925214859103339899D-26, &
  -1.740845657234000740989055147759702545340841422D-27, &
  4.15763564461389971961789962077522667348825413D-29, &
  -9.96214848828462210319400670245583884985485196D-31, &
  2.394034424896165300521167987893749562934279156D-32, &
  -5.76834735536739008429179316187765424407233225D-34, &
  1.393179479647007977827886603911548331732410612D-35, &
  -3.372121965485089470468473635254930958979742891D-37, &
  8.17820877756210262176477721487283426787618937D-39, &
  -1.987010831152385925564820669234786567541858996D-40, &
  4.83577851804055089628705937311537820769430091D-42 /
xm = -cln(1D0 - z, -isig)
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

!***********************************************************************
function contspence(z1, z2)
use global
implicit none
double complex :: contspence, z1, z2, z12, ris
double complex :: spence
external  spence

z12 = z1*z2
!ris = spence(1D0-z12,0D0)+log(1D0-z12)*(log(z12)-Log(z1)-log(z2))
!write(*,*) 'contspence=',ris

contspence = spence(1D0-z12,0D0)+log(1D0-z12)*(log(z12)-Log(z1)-log(z2))

end function contspence

!***********************************************************************
function cspence(z1, z2, im1, im2)
use global
implicit none
double complex :: cspence
double complex :: z1, z2
real*16 :: im1, im2
double complex :: cln, spence
integer :: eta
external cln, spence, eta
double complex :: z12
real*16 :: im12
integer :: etas

if (cdabs(z1)<1D-15.or.cdabs(z2)<1D-15) then
  cspence = pi**2/6D0
  return
end if
z12 = z1*z2
im12 = im2*sign(1D0, dble(z1))
if(dble(z12) > .5D0) then
  cspence = spence(1D0 - z12, 0D0)
  etas = eta(z1, z2, im1, im2, im12)
  if(etas /= 0) cspence = cspence +etas*dcmplx(0D0, 2D0*pi)*&
    &cln(1D0 - z12, -im12)
else if(abs(z12) < 1D-4) then
  cspence = pi**2/6D0 -spence(z12, 0D0) + (cln(z1, im1) + cln(z2, im2))*z12* &
    (1D0 + z12*(.5D0 + z12*(1D0/3D0 + z12/4D0)))
else
  cspence = pi**2/6D0 -spence(z12, 0D0) - &
    (cln(z1, im1) + cln(z2, im2))*cln(1D0 - z12, 0D0)
endif
end function cspence

!***********************************************************************
function eta(c1, c2, im1x, im2x, im12x)
use global
implicit none
integer :: eta
double complex :: c1, c2
real*16 :: im1x, im2x, im12x
real*16 :: im1, im2, im12
im1 = dimag(c1)
if(im1 == 0D0) im1 = im1x
im2 = dimag(c2)
if(im2 == 0D0) im2 = im2x
im12 = dimag(c1*c2)
if(im12 == 0D0) im12 = im12x
if(im1 < 0D0 .and. im2 < 0D0 .and. im12 > 0D0) then
  eta = 1
else if(im1 > 0D0 .and. im2 > 0D0 .and. im12 < 0D0) then
    eta = -1
  else
    eta = 0
#ifdef WARNINGS
    if(.not. (im2 == 0D0 .and. dble(c2) > 0D0 .or.&
      &im1 == 0D0 .and. dble(c1) > 0D0) .and. &
      (im1 == 0D0 .and. dble(c1) < 0D0 .or.&
      &im2 == 0D0 .and. dble(c2) < 0D0 .or.&
      &im12 == 0D0 .and. dble(c1*c2) < 0D0)) print *, "eta not defined"
#endif
  endif
end function eta

!***********************************************************************
function eta_n(c1, c2)
use global
implicit none
integer :: eta_n
double complex :: c1, c2
integer :: eta
external eta
eta_n = eta(c1, c2, 0D0, 0D0, 0D0)
end function eta_n

!***********************************************************************
function eta_tilde(c1, c2, im1x, im2x)
use global
implicit none
integer :: eta_tilde
double complex :: c1, c2
real*16 :: im1x, im2x
real*16 :: im1, im2
integer :: eta
external eta
im1 = dimag(c1)
if(im1 == 0D0) im1 = im1x
im2 = dimag(c2)
if(im2 /= 0D0) then
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
  if(im1 == 0D0 .and. dble(c1) < 0D0 .or.&
    &im2x == 0D0 .and. dble(c1*c2) < 0D0) print *, "eta_tilde not defined"
#endif
endif
end function eta_tilde

!***********************************************************************
function bdK(x, m1, m2)
! this is actually -K from the Beenakker/Denner paper for D0_ir
use global
implicit none
double complex :: bdK
double complex :: m1, m2, d
real*16 :: x
double complex :: zz
d = x - (m1 - m2)**2
!write(*,*) 'bdk=',x,m1,m2,d
if(d == dcmplx(0D0,0D0)) then
  bdK = dcmplx(1D0,0D0)
!write(*,*) 'bdK1=',bdK
else
  zz = sqrt(1D0 - 4D0*m1*m2/(d + IEPS))
  bdK = (zz - 1D0)/(zz + 1D0)
!write(*,*) 'bdK2=',zz,sqrt(1D0 - 4D0*m1*m2/(d - IEPS))
endif


end function bdK

!***********************************************************************
function A0_check(m)
use global
implicit none
double complex :: A0_check
real*16 :: m
real*16 :: mudim, divergence
common /cutoff/ mudim, divergence
if(m == 0D0) then
  A0_check = 0D0
else
  A0_check = m*(1D0 - log(m/mudim) + divergence)
endif
end function A0_check

!***********************************************************************
function B0_check(p, m1, m2)
use global
implicit none
double complex :: B0_check
real*16 :: p, m1, m2
real*16 :: mudim, divergence
common /cutoff/ mudim, divergence
double complex :: fpv, xlogx
external fpv, xlogx
double complex :: x1, x2, y1, y2, r
real*16 :: minacc
minacc = CALACC*(m1 + m2)
! general case
if(abs(p) > minacc) then
  call roots(p, m1, m2, x1, x2, y1, y2, r)
  if(abs(y1) > .5D0 .and. abs(y2) > .5D0) then
    B0_check = -log(m2/mudim) -fpv(1, x1, y1) - fpv(1, x2, y2)
  else if(abs(x1) < 10D0 .and. abs(x2) < 10D0) then
    B0_check = 2D0 - log(p*ONEmEPS/mudim) +&
      &xlogx(-x1) + xlogx(-x2) - xlogx(y1) - xlogx(y2)
  else if(abs(x1) > .5D0 .and. abs(x2) > .5D0) then
    B0_check = -log(m1/mudim) -fpv(1, y1, x1) - fpv(1, y2, x2)
  else
    print *, "B0(", p, ",", m1, ",", m2, ") not defined"
    B0_check = dcmplx(999D300)
  endif
! zero momentum
else if(abs(m1 - m2) > minacc) then
  x2 = ONEmEPS*m1/(m1 - m2)
  y2 = ONEmEPS*m2/(m2 - m1)
  if(abs(y2) > .5D0) then
    B0_check = -log(m2/mudim) - fpv(1, x2, y2)
  else
    B0_check = -log(m1/mudim) - fpv(1, y2, x2)
  endif
else
  B0_check = -log(m2/mudim)
endif
B0_check = B0_check + divergence
end function B0_check

!***********************************************************************
function DB0_check(p, m1, m2)
use global
implicit none
double complex :: DB0_check
real*16 :: p, m1, m2
real*16 :: delta
common /ffcut/ delta
double complex :: fpv, yfpv
external fpv, yfpv
double complex :: x1, x2, y1, y2, r
real*16 :: minacc
minacc = CALACC*(m1 + m2)
if(abs(p) > minacc) then
! IR divergent case
  if(m1*m2 == 0D0 .and. abs(p - m1 - m2) < CALACC) then
    DB0_check = -(1D0 + .5D0*log(delta/p))/p
    return
  endif
! general case
  call roots(p, m1, m2, x1, x2, y1, y2, r)
  if(abs(x1 - x2) > ROOT_CALACC*abs(x1 + x2)) then
    DB0_check = (yfpv(1, x2, y2) - yfpv(1, x1, y1))/r
  else if(abs(x1) > 10D0) then
    DB0_check = dble(-(.5D0 + (1D0 - 2D0*x1)*fpv(2, x1, y1))/x1**2)/p
  else if(abs(y1) > CALACC .and. abs(x1) > CALACC) then
    DB0_check = -dble(2D0 + (1D0 - 2D0*x1)*fpv(0, x1, y1))/p
  else
    print *, "DB0(", p, ",", m1, ",", m2, ") not defined"
    DB0_check = dcmplx(999D300)
  endif
! zero momentum
else if(abs(m1 - m2) > minacc) then
  x2 = ONEmEPS*m1/(m1 - m2)
  y2 = ONEmEPS*m2/(m2 - m1)
  if(abs(x2) < 10D0) then
    DB0_check = (.5D0 + yfpv(1, x2, y2))/(m1 - m2)
  else
    DB0_check = (.5D0 + yfpv(2, x2, y2))/m1
  endif
else
  DB0_check = 1D0/6D0/m1
endif
end function DB0_check

!***********************************************************************
function B1_check(p, m1, m2)
use global
implicit none
double complex :: B1_check
real*16 :: p, m1, m2
real*16 :: mudim, divergence
common /cutoff/ mudim, divergence
double complex :: fpv, xlogx
external fpv, xlogx
double complex :: x1, x2, y1, y2, r
real*16 :: minacc
minacc = CALACC*(m1 + m2)
! general case
if(abs(p) > minacc) then
  call roots(p, m1, m2, x1, x2, y1, y2, r)
  if(abs(y1) > .5D0 .and. abs(y2) > .5D0) then
    B1_check = .5D0*(log(m2/mudim) +fpv(2, x1, y1) + fpv(2, x2, y2))
  else if(abs(x1) < 10D0 .and. abs(x2) < 10D0) then
    B1_check = -.5D0*(1D0 - log(p*ONEmEPS/mudim) +&
      &x1*xlogx(-x1) + x1 + x2*xlogx(-x2) + x2 - &
      (1D0 + x1)*xlogx(y1) - (1D0 + x2)*xlogx(y2))
  else if(abs(x1) > .5D0 .and. abs(x2) > .5D0) then
    B1_check = .5D0*(log(m1/mudim) + 1D0 + &
      (1D0 + x1)*fpv(1, y1, x1) + (1D0 + x2)*fpv(1, y2, x2))
  else
    print *, "B1(", p, ",", m1, ",", m2, ") not defined"
    B1_check = dcmplx(999D300)
  endif
! zero momentum
else if(abs(m1 - m2) > minacc) then
  x2 = ONEmEPS*m1/(m1 - m2)
  y2 = ONEmEPS*m2/(m2 - m1)
  if(abs(y2) > .5D0) then
    B1_check = .5D0*(log(m2/mudim) + fpv(2, x2, y2))
  else
    B1_check = .5D0*(log(m1/mudim) + (1D0 + x2)*fpv(1, y2, x2) + .5D0)
  endif
else
  B1_check = .5D0*log(m2/mudim)
endif
B1_check = B1_check - .5D0*divergence
end function B1_check

!***********************************************************************
function DB1_check(p, m1, m2)
use global
implicit none
double complex :: DB1_check
real*16 :: p, m1, m2
real*16 :: delta
common /ffcut/ delta
double complex :: fpv, yfpv
external fpv, yfpv
double complex :: x1, x2, y1, y2, r
real*16 :: minacc
minacc = CALACC*(m1 + m2)
if(abs(p) > minacc) then
! IR divergent case
  if(m2 == 0D0 .and. abs(p - m1) < CALACC) then
    DB1_check = .5D0*(3D0 + log(delta/p))/p
    return
  endif
! general case
  call roots(p, m1, m2, x1, x2, y1, y2, r)
  if(abs(x1 - x2) > ROOT_CALACC*abs(x1 + x2)) then
    DB1_check = (yfpv(2, x1, y1) - yfpv(2, x2, y2))/r
  else if(abs(x1) > 10D0) then
    DB1_check = dble((2D0/3D0 + (2D0 - 3D0*x1)*fpv(3, x1, y1))/x1**2)/p
  else if(abs(y1) > CALACC) then
    DB1_check = dble(3D0/2D0 + (2D0 - 3D0*x1)*fpv(1, x1, y1))/p
  else
    print *, "DB1(", p, ",", m1, ",", m2, ") not defined"
    DB1_check = dcmplx(999D300)
  endif
! zero momentum
else if(abs(m1 - m2) > minacc) then
  x2 = ONEmEPS*m1/(m1 - m2)
  y2 = ONEmEPS*m2/(m2 - m1)
  if(abs(x2) < 10D0) then
    DB1_check = -(1D0/3D0 + yfpv(2, x2, y2))/(m1 - m2)
  else
    DB1_check = -(1D0/3D0 + yfpv(3, x2, y2))/m1
  endif
else
  DB1_check = -1D0/12D0/m1
endif
end function DB1_check

!***********************************************************************
function B00_check(p, m1, m2)
use global
implicit none
double complex :: B00_check
real*16 :: p, m1, m2
double complex :: A0_check, B0_check, B1_check
external A0_check, B0_check, B1_check
if(abs(p) > CALACC*(abs(p) + m1 + m2)) then
  B00_check = ((p + m1 - m2)*B1_check(p, m1, m2) +&
    &2D0*m1*B0_check(p, m1, m2) + A0_check(m2) +m1 + m2 - p/3D0)/6D0
else
  B00_check = (2D0*(m2*B0_check(0D0, m1, m2) +A0_check(m1)) + m1 + m2)/8D0
endif
end function B00_check

!***********************************************************************
function DB00_check(p, m1, m2)
use global
implicit none
double complex :: DB00_check
real*16 :: p, m1, m2
double complex :: B1_check, DB0_check, DB1_check
external B1_check, DB0_check, DB1_check
DB00_check = 1D0/6D0*(2D0*m1*DB0_check(p, m1, m2) +B1_check(p, m1, m2) + &
  (p + m1 - m2)*DB1_check(p, m1, m2)) - 1D0/18D0
end function DB00_check

!***********************************************************************
function B11_check(p, m1, m2)
use global
implicit none
double complex :: B11_check
real*16 :: p, m1, m2
real*16 :: mudim, divergence
common /cutoff/ mudim, divergence
double complex :: fpv, yfpv, xlogx
external fpv, yfpv, xlogx
double complex :: x1, x2, y1, y2, r
! general case
if(abs(p) > ROOT_CALACC*(m1 + m2)) then
  call roots(p, m1, m2, x1, x2, y1, y2, r)
  if(abs(y1) > .5D0 .and. abs(y2) > .5D0) then
    B11_check = (-log(m2/mudim) -fpv(3, x1, y1) - fpv(3, x2, y2))/3D0
  else if(abs(x1) < 10D0 .and. abs(x2) < 10D0) then
    x1 = x1**2*(xlogx(-x1) + 1D0) + .5D0*x1 - (1D0 + x1*(1D0 + x1))*xlogx(y1)
    x2 = x2**2*(xlogx(-x2) + 1D0) + .5D0*x2 - (1D0 + x2*(1D0 + x2))*xlogx(y2)
    B11_check = (2D0/3D0 - log(p*ONEmEPS/mudim) +x1 + x2)/3D0
  else if(abs(x1) > .5D0 .and. abs(x2) > .5D0) then
    x1 = (1D0 + x1*(1D0 + x1))*fpv(1, y1, x1) + .5D0*x1
    x2 = (1D0 + x2*(1D0 + x2))*fpv(1, y2, x2) + .5D0*x2
    B11_check = (-4D0/3D0 - log(m1/mudim) - x1 - x2)/3D0
  else
    print *, "B11(", p, ",", m1, ",", m2, ") not defined"
    B11_check = dcmplx(999D300)
  endif
! zero momentum
else if(abs(m1 - m2) > CALACC*(m1 + m2)) then
  x2 = ONEmEPS*m1/(m1 - m2)
  y2 = ONEmEPS*m2/(m2 - m1)
  if(abs(y2) > .5D0) then
    B11_check = (-log(m2/mudim) - fpv(3, x2, y2))/3D0
  else
    B11_check = (-log(m1/mudim) + (1D0 + x2*(1D0 + x2))*yfpv(0, x2, y2) +&
      &x2*(x2 + .5D0) + 1D0/3D0)/3D0
  endif
else
  B11_check = -log(m2/mudim)/3D0
endif
B11_check = B11_check + divergence/3D0
end function B11_check

!***********************************************************************
function DB11_check(p, m1, m2)
use global
implicit none
double complex :: DB11_check
real*16 :: p, m1, m2
double complex :: fpv, yfpv
external fpv, yfpv
double complex :: x1, x2, y1, y2, r
real*16 :: minacc
minacc = CALACC*(m1 + m2)
! general case
if(abs(p) > minacc) then
  call roots(p, m1, m2, x1, x2, y1, y2, r)
  if(abs(x1 - x2) > ROOT_CALACC*abs(x1 + x2)) then
    DB11_check = (yfpv(3, x2, y2) - yfpv(3, x1, y1))/r
    return
  else if(abs(x1) > 10D0) then
    DB11_check = dble((-3D0/4D0 + (4D0*x1 - 3D0)*fpv(4, x1, y1))/x1**2)/p
    return
  else if(abs(y1) > CALACC) then
    DB11_check = dble(-4D0/3D0 + (4D0*x1 - 3D0)*fpv(2, x1, y1))/p
    return
  endif
endif
! zero momentum
print *, "DB11(", p, ",", m1, ",", m2, ") not defined"
DB11_check = dcmplx(999D300)
end function DB11_check

!***********************************************************************
function C0_check(p1, p2, p1p2, m1, m2, m3)
use global
implicit none
double complex :: C0_check
real*16 :: p1, p2, p1p2, m1, m2, m3
double complex :: C0_ir, C0_reg
external C0_ir, C0_reg
if(m1 == 0D0 .and. (abs(p1 - m2) + abs(p1p2 - m3)) < CALACC) then
  C0_check = C0_ir(p2, p1, p1p2)
  return
endif
if(m2 == 0D0 .and. (abs(p1 - m1) + abs(p2 - m3)) < CALACC) then
  C0_check = C0_ir(p1p2, p1, p2)
  return
endif
if(m3 == 0D0 .and. (abs(p2 - m2) + abs(p1p2 - m1)) < CALACC) then
  C0_check = C0_ir(p1, p2, p1p2)
  return
endif
C0_check = C0_reg(p1, p2, p1p2, m1, m2, m3)
end function C0_check

!***********************************************************************
function C0_reg(p1, p2, p1p2, m1, m2, m3)
use global
implicit none
double complex :: C0_reg
real*16 :: p1, p2, p1p2, m1, m2, m3
real*16 :: q(5), m(5), mki, mkj, mij, qijk, ar
double complex :: a, b, h, h0, h1, h2, h3, h4
double complex :: y1, y2, y3, y4, x1, x2, x3, x4
integer :: i, j, k
double complex :: ln, spence
integer :: eta_n
external ln, spence, eta_n
q(1) = p1
q(2) = p2
q(3) = p1p2
q(4) = q(1)
q(5) = q(2)
m(1) = m1
m(2) = m2
m(3) = m3
m(4) = m(1)
m(5) = m(2)
C0_reg = dcmplx(0D0)
! all mom-squares != 0
if(p1*p2*p1p2 /= 0D0) then
  a = sqrt(dcmplx((p2 - p1 - p1p2)**2 - 4D0*p1*p1p2))
  do i = 1, 3
    j = i + 1
    k = i + 2
    mki = m(k) - m(i)
    mkj = m(k) - m(j)
    mij = m(i) - m(j)
    qijk = q(i) - q(j) - q(k)
    h2 = .5D0/a/q(i)
    h = q(i)*(qijk + mki + mkj) - mij*(q(j) - q(k))
    y1 = h2*(h + a*(q(i) - mij))
    y2 = h2*(h + a*(-q(i) - mij))
    b = sqrt(dcmplx((q(i) - mij)**2 - 4D0*q(i)*m(j)))
    y3 = h2*(h + a*b)
    y4 = h2*(h - a*b)
    h0 = q(i)*(q(j)*q(k) + qijk*m(k) + mki*mkj) -mij*(q(j)*mki - q(k)*mkj)
    qijk = q(j) - q(k) - q(i)
    h3 = h0 + q(j)*qijk*m(i) + q(k)*(q(k) - q(i) - q(j))*m(j)
    if(abs(y3) < abs(y4)) then
      y3 = h3/a**2/q(i)/y4
    else
      y4 = h3/a**2/q(i)/y3
    endif
    if(a*b /= dcmplx(0D0)) then
      y3 = y3 + IEPS/a/b*abs(a*b*y3)
      y4 = y4 - IEPS/a/b*abs(a*b*y4)
    else
      y3 = y3*ONEpEPS
      y4 = y4*ONEmEPS
    endif
    h1 = h2*(h - a*(q(i) - mij))
    if(abs(y1) < abs(h1)) then
      h3 = h0 + q(j)*qijk*m(i) + (q(k)*(q(i) + q(j)) - (q(i) - q(j))**2)*m(j)
      y1 = h3/a**2/q(i)/h1
    endif
    h1 = h2*(h + a*(q(i) + mij))
    if(abs(y2) < abs(h1)) then
      h3 = h0 + q(k)*(q(k) - q(i) - q(j))*m(j) + &
        (q(j)*(q(k) + q(i)) - (q(k) - q(i))**2)*m(i)
      y2 = h3/a**2/q(i)/h1
    endif
    C0_reg = C0_reg +spence(y2/y3, 0D0) + spence(y2/y4, 0D0) -&
      &spence(y1/y3, 0D0) - spence(y1/y4, 0D0)
    if(dimag(a) /= 0D0) then
      h3 = IEPS*abs(b)/b
      x1 = -.5D0*(q(i) - mij + b)/q(i) - h3
      x2 = -.5D0*(q(i) - mij - b)/q(i) - h3
      x3 = -.5D0*(-q(i) - mij + b)/q(i) - h3
      x4 = -.5D0*(-q(i) - mij - b)/q(i) - h3
      h3 = 1D0/y3
      h4 = 1D0/y4
      h = log(y1)*(eta_n(x1, x2) + eta_n(h3, h4) -&
        &eta_n(x1, h3) - eta_n(x2, h4) ) -&
        &log(y2)*(eta_n(x3, x4) + eta_n(h3, h4) -&
        &eta_n(x3, h3) - eta_n(x4, h4) ) +&
        &log(y3)*(eta_n(x1, h3) - eta_n(x3, h3)) +&
        &log(y4)*(eta_n(x2, h4) - eta_n(x4, h4))
      if(dimag(a) > 0D0 .and. q(i) < 0D0) h = h - log(y1/y2)
      C0_reg = C0_reg +dcmplx(0D0, 2D0*pi)*h
    endif
  enddo
  C0_reg = C0_reg/a
  return
endif
! one mom-square zero
if((p2*p1 + p1p2*p2 + p1*p1p2) /= 0D0) then
  if(p1 /= 0D0) then
    if(p2 == 0D0) then
      m(1) = m2
      m(2) = m3
      m(3) = m1
      q(1) = p2
      q(2) = p1p2
      q(3) = p1
    else
      m(1) = m3
      m(2) = m1
      m(3) = m2
      q(1) = p1p2
      q(2) = p1
      q(3) = p2
    endif
    m(4) = m(1)
    m(5) = m(2)
    q(4) = q(1)
    q(5) = q(2)
  endif
  ar = q(2) - q(3)
  do i = 2, 3
    j = i + 1
    k = i + 2
    mki = m(k) - m(i)
    mkj = m(k) - m(j)
    mij = m(i) - m(j)
    qijk = q(i) - q(j) - q(k)
    if(i == 2) then
      y1 = 2D0*q(2)*(mki + ar)
      y2 = 2D0*q(2)*mki
    else
      y1 = 2D0*q(3)*mkj
      y2 = 2D0*q(3)*(mkj - ar)
    endif
    h = q(i)*(qijk + mki + mkj) - mij*(q(j) - q(k))
    b = sqrt(dcmplx((q(i) - mij)**2 - 4D0*q(i)*m(j)))
    y3 = h + ar*b
    y4 = h - ar*b
    h0 = q(i)*(q(j)*q(k) + qijk*m(k) + mki*mkj) -mij*(q(j)*mki - q(k)*mkj)
    h3 = h0 + q(j)*(q(j) - q(k) - q(i))*m(i) +&
      &q(k)*(q(k) - q(i) - q(j))*m(j)
    h3 = 4D0*h3*q(i)
    if(abs(y3) < abs(y4)) then
      y3 = h3/y4
    else
      y4 = h3/y3
    endif
    qijk = ar/q(i)
    if(qijk /= 0D0) then
      y3 = y3 + IEPS/qijk*abs(qijk*y3)
      y4 = y4 - IEPS/qijk*abs(qijk*y4)
    else
      y3 = y3*ONEpEPS
      y4 = y4*ONEmEPS
    endif
    C0_reg = C0_reg +spence(y2/y3, 0D0) + spence(y2/y4, 0D0) -&
      &spence(y1/y3, 0D0) - spence(y1/y4, 0D0)
  enddo
  C0_reg = C0_reg/ar
  return
endif
! two mom-squares zero
if(p1p2 == 0D0) then
  if(p2 /= 0D0) then
    m(1) = m3
    m(2) = m1
    m(3) = m2
    q(1) = p1p2
    q(2) = p1
    q(3) = p2
  else
    m(1) = m2
    m(2) = m3
    m(3) = m1
    q(1) = p2
    q(2) = p1p2
    q(3) = p1
  endif
  m(4) = m(1)
  m(5) = m(2)
  q(4) = q(1)
  q(5) = q(2)
endif
mki = m(2) - m(3)
mkj = m(2) - m(1)
mij = m(3) - m(1)
if(m(2) /= m(3)) then
  y1 = -q(3) - mkj
  y2 = -mkj
  qijk = -mkj - q(3)*m(2)/mki
  y3 = qijk - IEPS*sign(1D0, -q(3)/mki)*abs(qijk)
  C0_reg = C0_reg + spence(y2/y3, 0D0) - spence(y1/y3, 0D0)
endif
b = sqrt(dcmplx((q(3) - mij)**2 - 4D0*q(3)*m(1)))
h = q(3)*(q(3) + mki + mkj)
y1 = 2D0*q(3)*mkj
y2 = 2D0*q(3)*(q(3) + mkj)
y3 = h - q(3)*b
y4 = h + q(3)*b
h0 = 4D0*q(3)**2*(q(3)*m(2) + mki*mkj)
if(abs(y3) < abs(y4)) then
  y3 = h0/y4
else
  y4 = h0/y3
endif
y3 = y3 - IEPS*abs(y3)
y4 = y4 + IEPS*abs(y4)
C0_reg = -(C0_reg +spence(y2/y3, 0D0) + spence(y2/y4, 0D0) -&
  &spence(y1/y3, 0D0) - spence(y1/y4, 0D0))/q(3)
end function C0_reg

!***********************************************************************
function C0_ir(p2, p1, p1p2)
use global
implicit none
double complex :: C0_ir
real*16 :: p2, p1, p1p2
real*16 :: delta
common /ffcut/ delta
double complex :: spence, ln
external spence, ln
real*16 :: a, h1, h2, h3, ps
double complex :: c
if(p1p2 == 0D0 .or. p1 == 0D0) then
  print *, "C0_ir: mass singular case"
  C0_ir = dcmplx(999D300)
  return
endif
if(p2 == 0D0) then
  C0_ir = -log(p1p2*p1/delta**2)*log(p1/p1p2)/4D0/(p1 - p1p2)
  return
endif
ps = p2 - p1 - p1p2
a = ps**2 - 4D0*p1*p1p2
if(a < 0D0) print *, "C0_ir: complex square root not implemented"
a = sqrt(a)
if(ps <= 0D0) then
  h1 = .5D0*(a - ps)
else
  h1 = -2D0*p1*p1p2/(a + ps)
endif
ps = p2 - p1 + p1p2
if(ps <= 0D0) then
  h2 = .5D0*(a - ps)
else
  h2 = -2D0*p2*p1p2/(a + ps)
endif
ps = p2 + p1 - p1p2
if(ps <= 0D0) then
  h3 = .5D0*(a - ps)
else
  h3 = -2D0*p1*p2/(a + ps)
endif
c = ln(-a/p2, -1D0)
C0_ir = (-pi**2/6D0 +&
  &spence(dcmplx(h2/a), -1D0) + spence(dcmplx(h3/a), -1D0) -&
  &.5D0*(ln(-h2/p2, -1D0)**2 + ln(-h3/p2, -1D0)**2) +&
  &.25D0*(ln(-p1/p2, -1D0)**2 + ln(-p1p2/p2, -1D0)**2) -&
  &c*(ln(-h1/p2, -1D0) - c) +ln(-delta/p2, -1D0)*ln(h1/sqrt(p1*p1p2), 1D0))/a
end function C0_ir

!***********************************************************************
function D0_check(cp1, cp2, cp3, cp4, cp1p2, cp2p3, &
  inm1, inm2, inm3, inm4)
use global
implicit none
double complex :: D0_check, D0_check2
real*16 :: p1, p2, p3, p4, p1p2, p2p3
double complex  :: cp1, cp2, cp3, cp4, cp1p2, cp2p3
double complex  :: cm1, cm2, cm3, cm4,inm1,inm2,inm3,inm4
double complex :: D0_ir, D0_m0, D0_reg
external D0_ir, D0_m0, D0_reg
! note: it is important to check _all_ IR divergent cases first,
! and only then the other zero-mass cases, otherwise the
! doubly-IR-divergent box is not caught


if (cp1p2-inm3**2.eq.dcmplx(0D0,0D0)) then
!write(*,*) 'reverse'
 cm1 = inm1
 cm2 = inm3
 cm3 = inm2
 cm4 = inm4
 p1 = qext(cp1p2)
 p2 = qext(cp2)
 p3 = qext(cp2p3)
 p4 = qext(cp4)
 p1p2 = qext(cp1)
 p2p3 = qext(cp3)
else
!write(*,*) 'straight'
 cm1 = inm1
 cm2 = inm2
 cm3 = inm3
 cm4 = inm4
 p1 = qext(cp1)
 p2 = qext(cp2)
 p3 = qext(cp3)
 p4 = qext(cp4)
 p1p2 = qext(cp1p2)
 p2p3 = qext(cp2p3)
 end if

!write(*,*) 'ciao ciao'
!write(*,*) 'cp1=',p1
!write(*,*) 'cp2=',p2
!write(*,*) 'cp3=',p3
!write(*,*) 'cp4=',p4
!write(*,*) 'cp1p2=',p1p2
!write(*,*) 'cp2p3=',p2p3
!write(*,*) 'cm1=',cm1
!write(*,*) 'cm2=',cm2
!write(*,*) 'cm3=',cm3
!write(*,*) 'cm4=',cm4

!write(*,*) 'diff1=',dble(p1)-dreal(cm2)**2
!write(*,*) 'diff1=',dble(p4)-dreal(cm4)**2
!write(*,*) 'diff1=',dble(p1p2)-dreal(cm3)**2


!!$if(cm1.eq.(0D0,0D0) .and. (abs(cm2 - cp1) + abs(cm3 - cp1p2)) < CALACC) then
!  p1 = qext(cp1)
!  p2 = qext(cp2)
!  p3 = qext(cp3)
!  p4 = qext(cp4)
!  p1p2 = qext(cp1p2)
!  p2p3 = qext(cp2p3)
!!$
!!$  D0_check = D0_ir(p1, p2, p3, p1p2, p4, p2p3, cm4)
!!$  return
!!$endif
!!$if(m2 == 0D0 .and. (abs(m1 - p1) + abs(m3 - p2)) < CALACC) then
!!$  D0_check = D0_ir(p2, p3, p4, p1, p2p3, p1p2, m4)
!!$  return
!!$endif
!!$if(m3 == 0D0 .and. (abs(m2 - p2) + abs(m4 - p3)) < CALACC) then
!!$  D0_check = D0_ir(p3, p4, p1, p2, p1p2, p2p3, m1)
!!$  return
!!$endif
!!$if(m4 == 0D0 .and. (abs(m1 - p4) + abs(m3 - p3)) < CALACC) then
!!$  D0_check = D0_ir(p4, p1, p2, p3, p2p3, p1p2, m2)
!!$  return
!!$endif



if(cm1 == dcmplx(0D0,0D0)) then
write(*,*) 'fin qui'
  D0_check = D0_m0(p3, p4, p1, p2, p1p2, p2p3, cm3, cm4, cm2)
  return
endif
if(cm2 == dcmplx(0D0,0D0)) then
  D0_check = D0_m0(p4, p1, p2, p3, p2p3, p1p2, cm4, cm1, cm3)
  return
endif
if(cm3 == dcmplx(0D0,0D0)) then
  D0_check = D0_m0(p1, p2, p3, p4, p1p2, p2p3, cm1, cm2, cm4)
  return
endif
if(cm4 == dcmplx(0D0,0D0)) then
  D0_check = D0_m0(p2, p3, p4, p1, p2p3, p1p2, cm2, cm3, cm1)
  return
endif

 cm1 = inm1
 cm2 = inm2
 cm3 = inm3
 cm4 = inm4
 p1 = qext(cp1)
 p2 = qext(cp2)
 p3 = qext(cp3)
 p4 = qext(cp4)
 p1p2 = qext(cp1p2)
 p2p3 = qext(cp2p3)

D0_check = D0_reg(p1, p2, p3, p4, p1p2, p2p3, cm1, cm2, cm3, cm4)





end function D0_check



!***********************************************************************
function D0_ir(cp1, cp2, cp3, cp4, cp1p2, cp2p3, m3)
use global
implicit none
double complex :: D0_ir
double complex :: cp1, cp2, cp3, cp4, cp1p2, cp2p3
real*16 :: p1, p2, p3, p4, p1p2, p2p3
real*16 :: d, delta
double complex :: m3, cm3, m1_, m4_
double complex :: xs, x2, x3, y, c, f
double complex :: logxs, logx2, logx3, log1x2, log1x3, logy
real*8 :: mgamma
common /cutoffIR/ mgamma
double complex :: ln, spence, bdK, addeps, cln
external ln, spence, bdK, addeps, cln

delta = qext(mgamma**2)
cm3 = m3**2

p1 = qext(cp1)
p2 = qext(cp2)
p3 = qext(cp3)
p4 = qext(cp4)
p1p2 = qext(cp1p2)
p2p3 = qext(cp2p3)

m1_ = sqrt(p1)
m4_ = sqrt(p4)
d = p2p3 - (m1_ - m4_)**2
f = .5D0/m1_/m4_/(p1p2 - cm3)


if(d /= 0D0) then
  xs = bdK(p2p3, m1_, m4_)
  logxs = log(xs)
  f = f*2D0*xs/(1D0 - xs**2)
endif
!write(*,*) 'f=',f

! massless case
if (m3.eq.(0.D0,0.D0)) then
!write(*,*) ' massless'
  if(p1.eq.p2 .and. p3.eq.p4) then
    D0_ir = 2D0*f*ln(-delta/p1p2, 1D0)
    if(d /= 0D0) D0_ir = -logxs*D0_ir
    return
  endif

  if(p1.eq.p2) then
!write(*,*) 'p1=p2',1D0/(p2p3-p4)/p1p2
!write(*,*) 'p1=p2',-m1_*m4_/(p2p3-p4)
!write(*,*) 'p1=p2',cln(-m1_*m4_/(p2p3-p4),1D0)
!write(*,*) 'p1=p2',2D0*ln((p3-p4)/p1p2,1D0)
    D0_ir = 1D0/(p2p3-p4)/p1p2*&
              (cln(-m1_*m4_/(p2p3-p4),1D0)*&
                  (cln(-m1_*m4_/(p2p3-p4),1D0)+log(delta/p4)+&
                   2D0*ln((p3-p4)/p1p2,1D0))+&
                                    pi**2/6D0)
    return
  endif

  if(p3.eq.p4) then
write(*,*) 'p3=p4'
    D0_ir = 1D0/(p2p3-p1)/p1p2*&
              (cln(-m4_*m1_/(p2p3-p1),1D0)*&
                  (cln(-m4_*m1_/(p2p3-p1),1D0)+log(delta/p1)+&
                   2D0*ln((p2-p1)/p1p2,1D0))+&
                                    pi**2/6D0)
    return
  endif

  y = m1_/m4_*(p3 - p4 + IEPS)/ (p2 - p1 + IEPS)
  logy = log(y)
  c = cln(delta/m1_/m4_, 0D0) +&
    &ln((p2 - p1)/p1p2, p1 - p2) + ln((p3 - p4)/p1p2, p4 - p3)

  if(d /= 0D0) then
    D0_ir = f*(pi**2/6D0 +logxs*(-.5D0*logxs + 2D0*log(1D0 - xs**2) - c) +&
      &spence(xs**2, 0D0) + .5D0*logy**2 -spence(xs/y, 0D0) - &
      (logxs + log(1D0/y))*log(1D0 - xs/y) -spence(xs*y, 0D0) - &
      (logxs + logy)*log(1D0 - xs*y))
!write(*,*) ' massless',D0_ir
    return
  endif
  D0_ir = f*(c - 2D0 - (1D0 + y)/(1D0 - y)*logy)
  return
endif


! massive case
!write(*,*) ' massive'


x2 = bdK(p2, m1_, m3)
x3 = bdK(p3, m4_, m3)
logx2 = log(x2)
logx3 = log(x3)
log1x3 = log(1D0/x3)
c = cln(m3*sqrt(delta)/(cm3 - p1p2), 1D0)

if(d /= 0D0) then
  log1x2 = log(1D0/x2)
  D0_ir = f*(.5D0*pi**2 +2D0*log(xs)*(log(1D0 - xs**2) - c) +&
    &spence(xs**2, 0D0) + logx2**2 + logx3**2 -spence(xs/x2/x3, 0D0) - &
    (logxs + log1x2 + log1x3)*log(1D0 - xs/x2/x3) -spence(xs*x2/x3, 0D0) - &
    (logxs + logx2 + log1x3)*log(1D0 - xs*x2/x3) -spence(xs/x2*x3, 0D0) - &
    (logxs + log1x2 + logx3)*log(1D0 - xs/x2*x3) -spence(xs*x2*x3, 0D0) - &
    (logxs + logx2 + logx3)*log(1D0 - xs*x2*x3))
  return
endif
D0_ir = f*(2D0*c - (1D0 + x2/x3)/(1D0 - x2/x3)*(logx2 + log1x3) - &
  (1D0 + x2*x3)/(1D0 - x2*x3)*(logx2 + logx3) - 2D0)
end function D0_ir



!*************************************************************************
 function D0_dresIRfin(cp1, cp2, cp3, cp4, cp1p2, cp2p3, m1, m2, m3, m4)
 implicit none
 complex*16 :: D0_dresIRfin,cp1,cp2,cp3,cp4,cp1p2,cp2p3,m1,m2,m3,m4,mv,mf
 complex*16 :: s,d1,d2,zeta,beta,y0,xs
 complex*16 :: addeps,contspence
 external addeps,contspence

   if (dimag(m2) /= 0D0.and.m2.eq.m3) then
write(*,*) 'caso 23'
     s = cp2
     mv = m2
mv = dcmplx(dreal(mv),0D0)
m2 = dcmplx(dreal(m2),0D0)
m3 = dcmplx(dreal(m3),0D0)
     mf = m4
     if (abs(cp3)<1.d-1) then
     d1 = cp1-m2**2
     d2 = cp1p2-m3**2
     zeta = 1D0-cp2p3/mv**2
     else
     d2 = cp1-m2**2
     d1 = cp1p2-m3**2
     zeta = 1D0-cp3/mv**2
     end if

     y0 = d1/d2
     beta = cdsqrt(1D0-4D0*mv**2/s)
     xs = (beta-1D0)/(beta+1D0)+dcmplx(0D0,EPS)

   else if (dimag(m3) /= 0D0.and.m3.eq.m4) then
write(*,*) 'caso 34'
     s = cp3
     mv = m3
mv = dcmplx(dreal(mv),0D0)
m3 = dcmplx(dreal(m3),0D0)
m4 = dcmplx(dreal(m4),0D0)
     mf = m2
     if (abs(cp2)<1.d-1) then
     d1 = cp4-m3**2
     d2 = cp1p2-m4**2
     zeta = 1D0-cp2p3/mv**2
     else
     d2 = cp4-m3**2
     d1 = cp1p2-m4**2
     zeta = 1D0-cp2/mv**2
     end if
     y0 = d1/d2
     beta = cdsqrt(1D0-4D0*mv**2/s)
     xs = (beta-1D0)/(beta+1D0)+dcmplx(0D0,EPS)

   else if (dimag(m2) /= 0D0.and.m2.eq.m4) then
write(*,*) 'caso 24'
     s = cp2p3
     mv = m2
mv = dcmplx(dreal(mv),0D0)
m2 = dcmplx(dreal(m2),0D0)
m4 = dcmplx(dreal(m4),0D0)
     mf = m3
     if (abs(cp2)<1.d-1) then
     d1 = cp4-m2**2
     d2 = cp1-m4**2
     zeta = 1D0-cp3/mv**2
     else
     d2 = cp4-m2**2
     d1 = cp1-m4**2
     zeta = 1D0-cp2/mv**2
     end if
     y0 = d1/d2
     beta = cdsqrt(1D0-4D0*mv**2/s)
     xs = (beta-1D0)/(beta+1D0)+dcmplx(0D0,EPS)

  end if

write(*,*) ' s=',s
write(*,*) 'd1=',d1
write(*,*) 'd2=',d2
write(*,*) 'norm=',1/mv**2/(d2-zeta*d1)
write(*,*) 'spence1=',-2D0*contspence(y0,zeta)
write(*,*) 'spence2=',contspence(y0,xs)
write(*,*) 'spence3=',contspence(y0,1D0/xs)
write(*,*) 'spence4=',-contspence(1D0/zeta,xs)
write(*,*) 'spence5=',-contspence(1D0/zeta,1D0/xs)

write(*,*) 'spence1=',-2D0*contspence(1D0/y0,1D0/zeta)
write(*,*) 'spence2=',contspence(1D0/y0,xs)
write(*,*) 'spence3=',contspence(1D0/y0,1D0/xs)
write(*,*) 'spence4=',-contspence(zeta,xs)
write(*,*) 'spence5=',-contspence(zeta,1D0/xs)
write(*,*) 'log    =',log(mv**2/mf**2)*(log(zeta)+log(y0)) 
write(*,*) 'log    =',log(mv**2/mf**2)*(log(zeta)-log(y0)) 
!write(*,*) 'log 2=',log(zeta)
!write(*,*) 'log 3=',log(y0)

 D0_dresIRfin = 1/mv**2/(d2-zeta*d1)*&
                (-2D0*contspence(y0,zeta)+&
                      contspence(y0,xs)+&
                      contspence(y0,1D0/xs)-&
                      contspence(1D0/zeta,xs)-&
                      contspence(1D0/zeta,1D0/xs)+&
                  log(mv**2/mf**2)*(log(y0)+log(zeta)))


end function D0_dresIRfin



!*************************************************************************
 function D0_dres(cp1, cp2, cp3, cp4, cp1p2, cp2p3, m1, m2, m3, m4)
 use global
 implicit none
 complex*16 :: D0_dres,cp1,cp2,cp3,cp4,cp1p2,cp2p3,m1,m2,m3,m4,mv
 real*8 :: s12,s13,s14,s23,s24,s34
 real*8 :: q1q2,q1k1,q1k2,q2k1,q2k2,k1k2
 complex*16 :: a,b,c,d,s, x1,x2,y1,y2, norm, ris
 complex*16 :: addeps,cspence,contspence,ln
 complex*16 :: a1,a2,a3
 external addeps,cspence,contspence,ln


!write(*,*) 'cp1=',cp1
!write(*,*) 'cp2=',cp2
!write(*,*) 'cp3=',cp3
!write(*,*) 'cp4=',cp4
!write(*,*) 'c12=',cp1p2
!write(*,*) 'c23=',cp2p3
!write(*,*) 'm1=',m1
!write(*,*) 'm2=',m2
!write(*,*) 'm3=',m3
!write(*,*) 'm4=',m4


         s12=dreal((-cp1 + m1**2 + m2**2)/2D0)
         s13=dreal((-cp1p2 + m1**2 + m3**2)/2D0)
         s14=dreal((-cp4 + m1**2 + m4**2)/2D0)
         s23=dreal((-cp2 + m2**2 + m3**2)/2D0)
         s24=dreal((-cp2p3 + m2**2 + m4**2)/2D0)
         s34=dreal((-cp3 + m3**2 + m4**2)/2D0)



      if ((dimag(m1)/=0D0).and.(m1.eq.m2)) then

         q1q2 = s12
         q1k1 = s14
         q1k2 = s13
         q2k1 = s24
         q2k2 = s23
         k1k2 = s34
         mv=m1

 else if ((dimag(m1)/=0D0).and.(m1.eq.m3)) then
         q1q2 = s13
         q1k1 = s14
         q1k2 = s12
         q2k1 = s34
         q2k2 = s23
         k1k2 = s24
         mv=m1

 else if ((dimag(m1)/=0D0).and.(m1.eq.m4)) then
         q1q2 = s14
         q1k1 = s13
         q1k2 = s12
         q2k1 = s34
         q2k2 = s24
         k1k2 = s23
         mv=m1

 else if ((dimag(m2)/=0D0).and.(m2.eq.m3)) then
         q1q2 = s23
         q1k1 = s24
         q1k2 = s12
         q2k1 = s34
         q2k2 = s13
         k1k2 = s14
         mv=m2

 else if ((dimag(m2)/=0D0).and.(m2.eq.m4)) then
         q1q2 = s24
         q1k1 = s23
         q1k2 = s12
         q2k1 = s34
         q2k2 = s14
         k1k2 = s13
         mv=m2

 else if ((dimag(m3)/=0D0).and.(m3.eq.m4)) then
         q1q2 = s34
         q1k1 = s23
         q1k2 = s13
         q2k1 = s24
         q2k2 = s14
         k1k2 = s12
         mv=m3

      end if




     a=2D0*(q2k1-k1k2)
     b=mv**2+4D0/mv**2*(q1k2*q2k1-q1q2*k1k2)
     c=2D0*(q1k2-k1k2)
     d=2D0*k1k2
     x1 = (-b+cdsqrt(b**2-4D0*a*c-dcmplx(0D0,4D0)*a*d*EPS))/2D0/a
     x2 = (-b-cdsqrt(b**2-4D0*a*c-dcmplx(0D0,4D0)*a*d*EPS))/2D0/a
     s=-2D0*q1q2/mv**2
     y1 = (-s-cdsqrt(s**2-4D0+dcmplx(0D0,EPS)))/2D0
     y2 = (-s+cdsqrt(s**2-4D0+dcmplx(0D0,EPS)))/2D0

! write(*,*) 's12=',s12
! write(*,*) 's13=',s13
! write(*,*) 's14=',s14
! write(*,*) 's23=',s23
! write(*,*) 's24=',s24
! write(*,*) 's34=',s34
! write(*,*) 'q1k2=',q1k2
! write(*,*) 'k1k2=',k1k2
! write(*,*) ' b^2=',b**2
! write(*,*) ' 4ac=',4D0*a*c
! write(*,*) '   c=',c
! write(*,*) '   d=',d
!if (abs(b**2-4D0*a*c)<1.d3) then
! write(*,*) 'disc=',b**2-4D0*a*c
!end if
! write(*,*) '  x1=',x1
! write(*,*) '  x2=',x2
! write(*,*) '   s=',s
! write(*,*) '  y1=',y1
! write(*,*) '  y2=',y2
! write(*,*) '  mv=',mv



     norm = a*(x1-x2)*mv**2
     ris =-(contspence(y1,-x1)+&
            contspence(1/y1,-x1)-&
            contspence(2D0*q2k1/mv**2-dcmplx(0D0,EPS),-x1)-&
            contspence(mv**2/2D0/q1k2+dcmplx(0D0,EPS),-x1)+&
                       log(-x1)*ln(qext(q1k2/k1k2),1D0))+&
           (contspence(y1,-x2)+&
            contspence(1/y1,-x2)-&
            contspence(2D0*q2k1/mv**2-dcmplx(0D0,EPS),-x2)-&
            contspence(mv**2/2D0/q1k2+dcmplx(0D0,EPS),-x2)+&
                       log(-x2)*ln(qext(q1k2/k1k2),1D0))


!write(*,*) 'chap=',ris/norm

     D0_dres = ris/norm/(-2d0)


 end function D0_dres


!***********************************************************************
function k2r(k)
use global
implicit none
double complex :: k2r
double complex :: k
k2r = .5D0*k*(1D0 + cdsqrt(dcmplx(1D0 - 4D0/k**2)))
end function k2r

!***********************************************************************
function addeps(k)
use global
implicit none
double complex :: addeps
double complex :: k
addeps = k*dcmplx(1D0, -sign(EPS, dreal(k)))
end function addeps

!***********************************************************************
function D0_m0(p1, p2, p3, p4, p1p2, p2p3, &
  cm1, cm2, cm4)
use global
implicit none
double complex :: D0_m0
real*16 :: p1, p2, p3, p4, p1p2, p2p3
double complex :: cm1, cm2, cm4
double complex :: D0_m00, cspence, cln, k2r
integer :: eta_tilde
external D0_m00, cspence, eta_tilde, cln, k2r
double complex :: m12, m22, m42
double complex :: k12, k13, k14, k23, k24, k34
real*16 :: ir12, ir14, ir24, ix1(2), ix4(2)
double complex :: r12, r14, r24, x4(2), x1
double complex :: a, b, c, d, disc
integer :: i
if(cm1 == dcmplx(0D0,0D0)) then
write(*,*) 'fin qui 21'
  D0_m0 = D0_m00(p1, p1p2, p3, p2p3, p2, p4, cm2, cm4)
  return
endif
if(cm2 == dcmplx(0D0,0D0)) then
write(*,*) 'fin qui 22'
  D0_m0 = D0_m00(p1, p2, p3, p4, p1p2, p2p3, cm1, cm4)
  return
endif
if(cm4 == dcmplx(0D0,0D0)) then
write(*,*) 'fin qui 23'
  D0_m0 = D0_m00(p4, p3, p2, p1, p1p2, p2p3, cm1, cm2)
  return
endif
m12 = cm1**2
m22 = cm2**2
m42 = cm4**2
k12 = (m12 + m22 - p1)/cm1/cm2
k13 = (m12 - p1p2)/m12
k14 = (m12 + m42 - p4)/cm1/cm4
k23 = (m22 - p2)/cm2/cm1
k24 = (m22 + m42 - p2p3)/cm2/cm4
k34 = (m42 - p3)/cm1/cm4
write(*,*) 'definiti i kij',k34
r12 = k2r(k12)
r14 = k2r(k14)
r24 = k2r(k24)
write(*,*) 'definiti i kij 2'
a = k34/r24 - k23
b = k13*(1D0/r24 - r24) + k12*k34 - k14*k23
c = k13*(k12 - r24*k14) + r24*k34 - k23
d = -k34*r24 + k23
disc = cdsqrt((k12*k34 - k13*k24 - k14*k23)**2 -&
  &4D0*(k13*(k13 - k23*(k12 - k14*k24)) +&
  &k23*(k23 - k24*k34) + k34*(k34 - k13*k14)))


x4(1) = .5D0/a*(-b + disc)
x4(2) = .5D0/a*(-b - disc)
if(cdabs(x4(1)) > cdabs(x4(2))) then
  x4(2) = c/a/x4(1)
else
  x4(1) = c/a/x4(2)
endif
!write(*,*) 'definiti i kij 3'
if(dreal(k12) < -2D0) then
  ir12 = sign(10D0, 1D0 - cdabs(r12))
else
  ir12 = 0D0
endif
if(dreal(k14) < -2D0) then
  ir14 = sign(10D0, 1D0 - cdabs(r14))
else
  ir14 = 0D0
endif
if(dreal(k24) < -2D0) then
  ir24 = sign(10D0, 1D0 - cdabs(r24))
else
  ir24 = 0D0
endif
!write(*,*) 'definiti i kij 4',k13
ix4(2) = sign(1D0, dble(d))
ix4(1) = -ix4(2)
ix1(1) = sign(1D0, ix4(1)*dble(r24))
ix1(2) = -ix1(1)
!b = dcmplx(k34/k13)
!c = dcmplx(k23/k13)
b = k34/k13
c = k23/k13
D0_m0 = dcmplx(0D0)
!write(*,*) 'definiti i kij 5'
do i = 1, 2
  x1 = -x4(i)/r24


  D0_m0 = D0_m0 + (2*i - 3)*( cspence(-x4(i), r14, -ix4(i), ir14) +&
    &cspence(-x4(i), 1D0/r14, -ix4(i), -ir14) -&
    &cspence(x1, r12, -ix1(i), ir12) -cspence(x1, 1D0/r12, -ix1(i), -ir12) -&
    &cspence(-x4(i),b,-ix4(i), -qext(k13))+cspence(x1,c,-ix1(i),-qext(k13)) -&
    &dcmplx(0D0, 2D0*pi)*eta_tilde(-x4(i), 1D0/r24, -ix4(i), -ir24)*( &
    cln((k12 - r24*(k14 + x4(i)) - x1)/d, &
    dble(-(r24 - 1d0/r24)*ix4(i)/d)) +cln(k13, -1D0) ) )
enddo
!write(*,*) 'definiti i kij 6'
D0_m0 = D0_m0/m12/cm2/cm4/a/(x4(1) - x4(2))
end function D0_m0

!***********************************************************************
function D0_m00(p1, p2, p3, p4, p1p2, p2p3, &
  cm1, cm4)
use global
implicit none
double complex :: D0_m00
real*16 :: p1, p2, p3, p4, p1p2, p2p3
double complex :: cm1, cm4
double complex :: D0_m000, cspence, k2r, addeps
external D0_m000, cspence, k2r, addeps
double complex :: m12, m42
double complex :: k12, k13, k14, k23, k24, k34
double complex :: k12c, k13c, k23c, k24c, k34c
double complex :: r14, x4(2)
double complex :: a, b, c, disc
integer :: i
if(cm1 == dcmplx(0D0,0D0)) then
  D0_m00 = D0_m000(p4, p1, p2, p3, p2p3, p1p2, cm4)
  return
endif
if(cm4 == dcmplx(0D0,0D0)) then
  D0_m00 = D0_m000(p1, p2, p3, p4, p1p2, p2p3, cm1)
  return
endif
write(*,*) 'fin qui 31'
m12 = cm1**2
m42 = cm4**2
k12 = (m12 - p1)/m12
k13 = (m12 - p1p2)/m12
k14 = (m12 + m42 - p4)/cm1/cm4
k23 = -p2/m12
k24 = (m42 - p2p3)/cm1/cm4
k34 = (m42 - p3)/cm1/cm4
a = k34*k24 - k23
b = k13*k24 + k12*k34 - k14*k23
c = k13*k12 - ONEmEPS*k23
disc = cdsqrt(b*b - 4D0*a*c)
x4(1) = .5D0/a*(-b + disc)
x4(2) = .5D0/a*(-b - disc)
write(*,*) 'fin qui 32'
if(abs(x4(1)) > abs(x4(2))) then
  x4(2) = c/a/x4(1)
else
  x4(1) = c/a/x4(2)
endif
write(*,*) 'fin qui 33'
k12c = addeps(k12)
k13c = addeps(k13)
k23c = addeps(k23)
k24c = addeps(k24)/k12c
k34c = addeps(k34)/k13c
write(*,*) 'fin qui 34'

c = log(k12c) + log(k13c) - log(k23c)
r14 = k2r(k14)
r14 = r14*dcmplx(1D0, sign(EPS, dble(1D0/r14 - r14)))
D0_m00 = dcmplx(0D0)
do i = 1, 2
  D0_m00 = D0_m00 + (2*i - 3)*( cspence(-x4(i), r14, 0D0, 0D0) +&
    &cspence(-x4(i), 1D0/r14, 0D0, 0D0) -cspence(-x4(i), k34c, 0D0, 0D0) -&
    &cspence(-x4(i), k24c, 0D0, 0D0) +log(-x4(i))*c )
enddo
D0_m00 = D0_m00/m12/cm1/cm4/a/(x4(1) - x4(2))
end function D0_m00

!***********************************************************************
function D0_m000(p1, p2, p3, p4, p1p2, p2p3, cm1)
use global
implicit none
double complex :: D0_m000
real*16 :: p1, p2, p3, p4, p1p2, p2p3
double complex :: cm1,m12
double complex :: D0_m0000, cspence, addeps
external D0_m0000, cspence, addeps
double complex :: k12, k13, k14, k23, k24, k34
double complex :: k12c, k13c, k14c, k23c, k24c, k34c
double complex :: a, b
double complex :: c, disc, x4(2)
integer :: i
if(cm1 == dcmplx(0D0,0D0)) then
  D0_m000 = D0_m0000(p1, p2, p3, p4, p1p2, p2p3)
  return
endif
m12 = cm1**2
k12 = (m12 - p1)/m12
k13 = (m12 - p1p2)/m12
k14 = (m12 - p4)/m12
k23 = -p2/m12
k24 = -p2p3/m12
k34 = -p3/m12
a = k34*k24
b = k13*k24 + k12*k34 - k14*k23
c = k13*k12 - ONEmEPS*k23
disc = cdsqrt(b*b - 4D0*a*c)
x4(1) = .5D0/a*(-b + disc)
x4(2) = .5D0/a*(-b - disc)
if(abs(x4(1)) > abs(x4(2))) then
  x4(2) = c/a/x4(1)
else
  x4(1) = c/a/x4(2)
endif
k12c = addeps(k12)
k13c = addeps(k13)
k23c = addeps(k23)
k14c = addeps(k14)
k24c = addeps(k24)/k12c
k34c = addeps(k34)/k13c
c = log(k12c) + log(k13c) - log(k23c)
D0_m000 = dcmplx(0D0)
do i = 1, 2
  D0_m000 = D0_m000 + (2*i - 3)*( cspence(-x4(i), k14c, 0D0, 0D0) -&
    &cspence(-x4(i), k34c, 0D0, 0D0) -cspence(-x4(i), k24c, 0D0, 0D0) +&
    &log(-x4(i))*c )
enddo
D0_m000 = D0_m000/m12/a/(x4(1) - x4(2))
end function D0_m000

!***********************************************************************
function D0_m0000(p1, p2, p3, p4, p1p2, p2p3)
use global
implicit none
double complex :: D0_m0000
real*16 :: p1, p2, p3, p4, p1p2, p2p3
double complex :: cspence, addeps
external cspence, addeps
real*16 :: m2
real*16 :: k12, k13, k14, k23, k24, k34
double complex :: k12c, k13c, k14c, k23c, k24c, k34c
real*16 :: a, b
double complex :: c, disc, x4(2)
integer :: i
m2 = abs(p2p3)
k12 = -p1/m2
k13 = -p1p2/m2
k14 = -p4/m2
k23 = -p2/m2
k24 = -p2p3/m2
k34 = -p3/m2
a = k34*k24
b = k13*k24 + k12*k34 - k14*k23
c = k13*k12 + IEPS*k23
disc = cdsqrt(b*b - 4D0*a*c)
x4(1) = .5D0/a*(-b + disc)
x4(2) = .5D0/a*(-b - disc)
if(abs(x4(1)) > abs(x4(2))) then
  x4(2) = c/a/x4(1)
else
  x4(1) = c/a/x4(2)
endif
k12c = addeps(k12)
k13c = addeps(k13)
k23c = addeps(k23)
k14c = addeps(k14)
k24c = addeps(k24)/k12c
k34c = addeps(k34)/k13c
c = log(k12c) + log(k13c) - log(k23c) - log(k14c)
D0_m0000 = dcmplx(0D0)
do i = 1, 2
  disc = log(-x4(i))
  D0_m0000 = D0_m0000 + (2*i - 3)*( -cspence(-x4(i), k34c, 0D0, 0D0) -&
    &cspence(-x4(i), k24c, 0D0, 0D0) +disc*(c - .5D0*disc) )
enddo
D0_m0000 = D0_m0000/m2**2/a/(x4(1) - x4(2))
end function D0_m0000

!***********************************************************************
function D0_reg(p1, p2, p3, p4, p1p2, p2p3, &
  cm1, cm2, cm3, cm4)
use global
implicit none
double complex :: D0_reg
real*16 :: p1, p2, p3, p4, p1p2, p2p3
double complex :: cm1, cm2, cm3, cm4
double complex :: cspence, cln, k2r
integer :: eta
external cspence, cln, eta, k2r
double complex :: m12, m22, m32, m42
double complex :: k12, k13, k14, k23, k24, k34
double complex :: ir12, ir14, ir23, ir24, ir34
double complex :: r12, r14, r13, r23, r24, r34
double complex :: x(2, 4), s(4)
real*16 :: ix(2, 4), is(4)
double complex :: a, b, c, disc
integer :: j, k
m12 = cm1**2
m22 = cm2**2
m32 = cm3**2
m42 = cm4**2
k12 = (m12 + m22 - p1)/cm1/cm2
k13 = (m12 + m32 - p1p2)/cm1/cm3
k14 = (m12 + m42 - p4)/cm1/cm4
k23 = (m22 + m32 - p2)/cm2/cm3
k24 = (m22 + m42 - p2p3)/cm2/cm4
k34 = (m32 + m42 - p3)/cm3/cm4
write(*,*) 'D0_reg  definiti kij ',k12,k13,k14,k23,k24,k34
#ifdef WARNINGS
if(dreal(k13) < -2D0) print *, "D0_reg: case k13 < 0 not implemented."
#endif
r12 = k2r(k12)
r13 = 1D0/k2r(k13)
r14 = k2r(k14)
r23 = k2r(k23)
r24 = 1D0/k2r(k24)
r34 = k2r(k34)
a = k34/r24 - k23 + (k12 - k14/r24)*r13
b = (1D0/r13 - r13)*(1D0/r24 - r24) + k12*k34 - k14*k23
c = k34*r24 - k23 + (k12 - k14*r24)/r13
disc = cdsqrt(b*b - 4D0*a*c)
x(1, 4) = .5D0/a*(-b + disc)
x(2, 4) = .5D0/a*(-b - disc)
if(cdabs(x(1, 4)) > cdabs(x(2, 4))) then
  x(2, 4) = c/a/x(1, 4)
else
  x(1, 4) = c/a/x(2, 4)
endif
if(dreal(k12) < -2D0) then
  ir12 = sign(10D0, 1D0 - cdabs(r12))
else
  ir12 = 0D0
endif
if(dreal(k14) < -2D0) then
  ir14 = sign(10D0, 1D0 - cdabs(r14))
else
  ir14 = 0D0
endif
if(dreal(k23) < -2D0) then
  ir23 = sign(10D0, 1D0 - cdabs(r23))
else
  ir23 = 0D0
endif
if(dreal(k24) < -2D0) then
  ir24 = sign(10D0, 1D0 - cdabs(r24))
else if(dreal(k24) == -2D0) then
  ir24 = 10D0
else
  ir24 = 0D0
endif
if(dreal(k34) < -2D0) then
  ir34 = sign(10D0, 1D0 - cdabs(r34))
else
  ir34 = 0D0
endif
x(1, 1) = x(1, 4)/r24
x(2, 1) = x(2, 4)/r24
x(1, 2) = x(1, 4)/r24*r13
x(2, 2) = x(2, 4)/r24*r13
x(1, 3) = x(1, 4)*r13
x(2, 3) = x(2, 4)*r13
if(dble(x(1, 4)) > 0D0) then
  ix(1, 4) = 1D0
else
  ix(1, 4) = 0D0
endif
if(dble(x(2, 4)) > 0D0) then
  ix(2, 4) = -1D0
else
  ix(2, 4) = 0D0
endif
ix(1, 1) = ix(1, 4) + ir24
if(dble(x(1, 1)) <= 0D0) ix(1, 1) = -ix(1, 1)
ix(2, 1) = ix(2, 4) + ir24
if(dble(x(2, 1)) <= 0D0) ix(2, 1) = -ix(2, 1)
ix(1, 3) = ix(1, 4)
ix(2, 3) = ix(2, 4)
ix(1, 2) = ix(1, 1)
ix(2, 2) = ix(2, 1)

s(1) = r12
s(2) = r23
s(3) = r34
s(4) = r14
is(1) = ir12
is(2) = ir23
is(3) = ir34
is(4) = ir14
D0_reg = dcmplx(0D0)
do k = 1, 2
  do j = 1, 4
    D0_reg = D0_reg - (2*mod(j + k, 2) - 1)*( &
      cspence(-x(k, j), s(j), -ix(k, j), is(j)) +&
      &cspence(-x(k, j), 1D0/s(j), -ix(k, j), -is(j)) )
  enddo
  b = 1D0 + (k34 + x(k, 3))*x(k, 3)
  D0_reg = D0_reg + (2*mod(k, 2) - 1)*( &
    eta(-x(k, 4), 1D0/r24, -ix(k, 4), -ir24, -ix(k, 1))*dcmplx(0D0, 2D0*pi)*&
    &cln((1D0 + (k14 + x(k, 4))*x(k, 4))/b, -dble(b)) )
enddo
D0_reg = D0_reg/cm1/cm2/cm3/cm4/disc
end function D0_reg

!***********************************************************************
subroutine check_A(result, func, name, m)
use global
implicit none
double complex :: result, func
external func
character :: name*(*)
real*16 :: m
real*16 :: maxdev
common /ffcheck/ maxdev
double complex :: check
check = func(m)
if(abs(result - check)/abs(result) > maxdev) then
  print *, "deviation in ", name
  print *, "  m = ", m
  print *, "FF's result: ", result
  print *, "check:       ", check
endif
end subroutine check_A

!***********************************************************************
subroutine check_B(result, func, name, p, m1, m2)
use global
implicit none
double complex :: result, func
external func
character :: name*(*)
real*16 :: p, m1, m2
real*16 :: maxdev
common /ffcheck/ maxdev
double complex :: check
check = func(p, m1, m2)
if(abs(result - check)/abs(result) > maxdev) then
  print *, "deviation in ", name
  print *, "  p  = ", p
  print *, "  m1 = ", m1
  print *, "  m2 = ", m2
  print *, "FF's result: ", result
  print *, "check:       ", check
endif
end subroutine check_B

!***********************************************************************
subroutine check_C(result, func, name, p1, p2, p1p2, m1, m2, m3)
use global
implicit none
double complex :: result, func
external func
character :: name*(*)
real*16 :: p1, p2, p1p2, m1, m2, m3
real*16 :: maxdev
common /ffcheck/ maxdev
double complex :: check
check = func(p1, p2, p1p2, m1, m2, m3)
if(abs(result - check)/abs(result) > maxdev) then
  print *, "deviation in ", name
  print *, "  p1   = ", p1
  print *, "  p2   = ", p2
  print *, "  p1p2 = ", p1p2
  print *, "  m1   = ", m1
  print *, "  m2   = ", m2
  print *, "  m3   = ", m3
  print *, "FF's result: ", result
  print *, "check:       ", check
endif
end subroutine check_C

!***********************************************************************
subroutine check_D(result, func, name, &
  p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4)
use global
implicit none
double complex :: result, func
external func
character :: name*(*)
real*16 :: p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4
real*16 :: maxdev
common /ffcheck/ maxdev
double complex :: check
check = func(p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4)
if(abs(result - check)/abs(result) > maxdev) then
  print *, "deviation in ", name
  print *, "  p1   = ", p1
  print *, "  p2   = ", p2
  print *, "  p3   = ", p3
  print *, "  p4   = ", p4
  print *, "  p1p2 = ", p1p2
  print *, "  p2p3 = ", p2p3
  print *, "  m1   = ", m1
  print *, "  m2   = ", m2
  print *, "  m3   = ", m3
  print *, "  m4   = ", m4
  print *, "FF's result: ", result
  print *, "check:       ", check
endif
end subroutine check_D

