!===============================================================================
! OLAM version 3.3

! Copyright (C) 2002-2008; All Rights Reserved; 
! Duke University, Durham, North Carolina, USA 

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

! It is requested that any scientific publications based on application of OLAM
! include the following acknowledgment:  "OLAM was developed at the 
! Edmund T. Pratt Jr. School of Engineering, Duke University."

! For additional information, including published references, please contact
! the software authors, Robert L. Walko (robert.walko@duke.edu)
! or Roni Avissar (avissar@duke.edu).
!===============================================================================
real function gammp(a,x)

implicit none

real,intent(in) :: a,x

real :: gln,gammcf

if (x < a + 1.) then
   call gser(gammp,a,x,gln)
else
   call gcf(gammcf,a,x,gln)
   gammp = 1. - gammcf
endif

return
end function gammp

!===============================================================================

real function gammq(a,x)

implicit none

real, intent(in) :: a,x

real :: gamser,gln

if (x < a + 1.) then
   call gser(gamser,a,x,gln)
   gammq = 1. - gamser
else
   call gcf(gammq,a,x,gln)
endif

return
end function gammq

!===============================================================================

subroutine gcf(gammcf,a,x,gln)

implicit none

real, intent(out) :: gammcf,gln
real, intent(in) :: a,x

integer, parameter :: itmax=100
real, parameter :: eps=3.e-7
real,external :: gammln

real :: gold,a0,a1,b0,b1,fac,an,ana,anf,gaccel
integer :: n

gln = gammln(a)
gold = 0.
a0 = 1.
a1 = x
b0 = 0.
b1 = 1.
fac = 1.
do n = 1, itmax
    an = float(n)
    ana = an - a
    a0 = (a1 + a0 * ana) * fac
    b0 = (b1 + b0 * ana) * fac
    anf = an * fac
    a1 = x * a0 + anf * a1
    b1 = x * b0 + anf * b1
    if (a1 /= 0.) then
        fac = 1. / a1
        gaccel = b1 * fac
        if (abs((gaccel - gold) / gaccel) < eps) goto 20
        gold = gaccel
    endif
enddo

20 continue

gammcf = exp(-x + a * alog(x) - gln) * gaccel

if ((-x + a * log(x) - gln) > -38.) then
  gammcf = exp(-x + a * alog(x) - gln) * gaccel
else
  gammcf = 0.
endif

return
end subroutine gcf

!===============================================================================

subroutine gser(gamser,a,x,gln)

implicit none

real, intent(out) :: gamser,gln
real, intent(in) :: a,x

integer, parameter :: itmax=100
real, parameter :: eps=3.e-7
real,external ::gammln

real :: ap,sum,del
integer :: n

gln = gammln(a)

if (x <= 0.) then
    gamser = 0.
    return
endif

ap = a
sum = 1. / a
del = sum

do n = 1, itmax
    ap = ap + 1.
    del = del * x / ap
    sum =  sum  +  del
    if (abs(del) < abs(sum) * eps) goto 20
enddo

20 continue

if ((-x + a * log(x) - gln) > -38.) then
  gamser = sum * exp(-x + a * log(x) - gln)
else
  gamser = 0.
endif

return
end subroutine gser

!===============================================================================

real function gammln(xx)

implicit none

real, intent(in) :: xx

real(kind=8) :: cof(6),stp
data cof, stp/76.18009173d0, -86.50532033d0, 24.01409822d0,  &
     -1.231739516d0, .120858003d-2, -.536382d-5, 2.50662827465d0/
real(kind=8), parameter :: half=0.5d0, one=1.0d0, fpf=5.5d0

real :: x,tmp,ser
integer :: j

x = xx - one
tmp = x + fpf
tmp = (x + half) * log(tmp) - tmp
ser = one
do j = 1,6
    x = x + one
    ser = ser + cof(j) / x
enddo

gammln = tmp + log(stp * ser)

return
end function gammln

!===============================================================================

subroutine avint(x,y,n,xlo,xup,ans)

use misc_coms, only: io6

implicit none

integer, intent(in) :: n
real, intent(in) :: x(n),y(n),xlo,xup
real, intent(out) :: ans

real(kind=8) ::r3,rp5,sum,syl,syl2,syl3,syu,syu2,syu3,x1,x2,x3  &
   ,x12,x13,x23,term1,term2,term3,a,b,c,ca,cb,cc

integer :: i,inlft,inrt,istart,istop
real :: slope,fl,fr

ans = 0.0

if (xlo < xup) goto 3
if (xlo == xup) goto 100
if (xlo > xup) goto 200
3 if (n < 2) goto 215

do i = 2,n
   if (x(i) <= x(i-1)) goto 210
   if (x(i) > xup) goto 6
enddo

6 continue

if (n >=3) goto 9

! special n=2 case

slope = (y(2) - y(1)) / (x(2) - x(1))
fl = y(1) + slope * (xlo - x(1))
fr = y(2) + slope * (xup - x(2))
ans = .5 * (fl + fr) * (xup - xlo)

return

9 continue

if (x(n-2) < xlo)  goto 205
if (x(3) > xup)    goto 205
i = 1
10 if (x(i) >= xlo) goto 15
i = i + 1

goto 10

15 inlft = i
i = n
20 if (x(i) <= xup) goto 25
i = i - 1

goto 20

25 inrt = i
if ((inrt - inlft) < 2) goto 205
istart = inlft
if (inlft == 1) istart = 2
istop  = inrt
if (inrt == n)  istop  = n-1

r3 = 3.0d0
rp5 = 0.5d0
sum = 0.0
syl = xlo
syl2 = syl * syl
syl3 = syl2 * syl

do i = istart,istop
   x1 = x(i-1)
   x2 = x(i)
   x3 = x(i+1)
   x12 = x1 - x2
   x13 = x1 - x3
   x23 = x2 - x3
   term1 =  dble(y(i-1)) / (x12 * x13)
   term2 = -dble(y(i))   / (x12 * x23)
   term3 =  dble(y(i+1)) / (x13 * x23)
   a = term1 + term2 + term3
   b = -(x2 + x3) * term1 - (x1 + x3) * term2 - (x1 + x2) * term3
   c = x2 * x3 * term1 + x1 * x3 * term2 + x1 * x2 * term3

   if (i <= istart) goto 30
   if (i > istart) goto 35

30 ca = a
   cb = b
   cc = c

   goto 40

35 ca = .5 * (a + ca)
   cb = .5 * (b + cb)
   cc = .5 * (c + cc)
40 syu  = x2
   syu2 = syu * syu
   syu3 = syu2 * syu
   sum  = sum + ca * (syu3 - syl3) / r3 + cb * rp5 * (syu2 - syl2)  &
        + cc * (syu - syl)
   ca   = a
   cb   = b
   cc   = c
   syl  = syu
   syl2 = syu2
   syl3 = syu3
enddo

syu = xup
ans = sum + ca * (syu**3 - syl3) / r3 + cb * rp5 * (syu**2 - syl2)  &
          + cc * (syu - syl)
100 return
200 write(io6,*) 'Upper limit of integration not greater than lower limit.'
stop 'avint2'
205 write(io6,*) 'Less than 3 function values between integration limits.'
stop 'avint3'
210 write(io6,*) 'Abscissas not strictly increasing.'
stop 'avint4'
215 write(io6,*) 'Less than 2 function values were supplied.'
stop 'avint5'

end subroutine avint
