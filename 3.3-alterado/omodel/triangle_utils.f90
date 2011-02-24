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
subroutine tri_neighbors()

use mem_ijtabs, only: itab_m, itab_u, itab_w
use mem_grid,   only: mua, mwa
use misc_coms,  only: io6

implicit none

integer :: iu,iw
integer :: iu1,iu2,iu3,iu4
integer :: iw1,iw2,iw3,iw4,iw5,iw6
integer :: iw1_iu1,iw1_iu2,iw1_iu3
integer :: iw2_iu1,iw2_iu2,iw2_iu3
integer :: j,iunow,iu0,ntpn,im

! Loop over W points

do iw = 2,mwa
   iu1 = itab_w(iw)%iu1
   iu2 = itab_w(iw)%iu2
   iu3 = itab_w(iw)%iu3

! Fill M points for current W point

   if (iw == itab_u(iu2)%iw1) then
      itab_w(iw)%im1   = itab_u(iu2)%im1
      itab_w(iw)%im3   = itab_u(iu2)%im2
   else
      itab_w(iw)%im1   = itab_u(iu2)%im2
      itab_w(iw)%im3   = itab_u(iu2)%im1
   endif   	

   if (iw == itab_u(iu1)%iw1) then
      itab_w(iw)%im2   = itab_u(iu1)%im2
   else
      itab_w(iw)%im2   = itab_u(iu1)%im1
   endif

! Fill W points and inner U directions for current W point

   if (iw == itab_u(iu1)%iw1) then
      itab_w(iw)%iw1   = itab_u(iu1)%iw2
      itab_w(iw)%diru1 = -1.
   else
      itab_w(iw)%iw1   = itab_u(iu1)%iw1
      itab_w(iw)%diru1 = 1.
   endif

   if (iw == itab_u(iu2)%iw1) then
      itab_w(iw)%iw2   = itab_u(iu2)%iw2
      itab_w(iw)%diru2 = -1.
   else
      itab_w(iw)%iw2   = itab_u(iu2)%iw1
      itab_w(iw)%diru2 = 1.
   endif   	

   if (iw == itab_u(iu3)%iw1) then
      itab_w(iw)%iw3   = itab_u(iu3)%iw2
      itab_w(iw)%diru3 = -1.
   else
      itab_w(iw)%iw3   = itab_u(iu3)%iw1
      itab_w(iw)%diru3 = 1.
   endif   	

! Fill outer U points for current W point

   iw1 = itab_w(iw)%iw1
   iw2 = itab_w(iw)%iw2
   iw3 = itab_w(iw)%iw3

! This should work ok for iw1 = 1, but may bypass if desired

   if (iu1 == itab_w(iw1)%iu1) then
      itab_w(iw)%iu4 = itab_w(iw1)%iu2
      itab_w(iw)%iu5 = itab_w(iw1)%iu3
   elseif (iu1 == itab_w(iw1)%iu2) then
      itab_w(iw)%iu4 = itab_w(iw1)%iu3
      itab_w(iw)%iu5 = itab_w(iw1)%iu1
   else
      itab_w(iw)%iu4 = itab_w(iw1)%iu1
      itab_w(iw)%iu5 = itab_w(iw1)%iu2
   endif

! This should work ok for iw2 = 1, but may bypass if desired

   if (iu2 == itab_w(iw2)%iu1) then
      itab_w(iw)%iu6 = itab_w(iw2)%iu2
      itab_w(iw)%iu7 = itab_w(iw2)%iu3
   elseif (iu2 == itab_w(iw2)%iu2) then
      itab_w(iw)%iu6 = itab_w(iw2)%iu3
      itab_w(iw)%iu7 = itab_w(iw2)%iu1
   else
      itab_w(iw)%iu6 = itab_w(iw2)%iu1
      itab_w(iw)%iu7 = itab_w(iw2)%iu2
   endif

! This should work ok for iw3 = 1, but may bypass if desired

   if (iu3 == itab_w(iw3)%iu1) then
      itab_w(iw)%iu8 = itab_w(iw3)%iu2
      itab_w(iw)%iu9 = itab_w(iw3)%iu3
   elseif (iu3 == itab_w(iw3)%iu2) then
      itab_w(iw)%iu8 = itab_w(iw3)%iu3
      itab_w(iw)%iu9 = itab_w(iw3)%iu1
   else
      itab_w(iw)%iu8 = itab_w(iw3)%iu1
      itab_w(iw)%iu9 = itab_w(iw3)%iu2
   endif

enddo

! Loop over U points

do iu = 2,mua

   iw1 = itab_u(iu)%iw1
   iw2 = itab_u(iu)%iw2

   itab_u(iu)%mrlu = max(itab_w(iw1)%mrlw,itab_w(iw2)%mrlw)

! This should work ok for iw1 = 1, but may bypass if desired

   iw1_iu1 = itab_w(iw1)%iu1
   iw1_iu2 = itab_w(iw1)%iu2
   iw1_iu3 = itab_w(iw1)%iu3

! Fill IU1 and IU2 for current U point

   if (iw1_iu1 == iu) then
      itab_u(iu)%iu1 = iw1_iu2	
      itab_u(iu)%iu2 = iw1_iu3
   elseif (iw1_iu2 == iu) then
      itab_u(iu)%iu1 = iw1_iu3	
      itab_u(iu)%iu2 = iw1_iu1
   else
      itab_u(iu)%iu1 = iw1_iu1	
      itab_u(iu)%iu2 = iw1_iu2
   endif

! This should work ok for iw2 = 1, but may bypass if desired

   iw2_iu1 = itab_w(iw2)%iu1
   iw2_iu2 = itab_w(iw2)%iu2
   iw2_iu3 = itab_w(iw2)%iu3

! Fill IU3 and IU4 for current U point

   if (iw2_iu1 == iu) then
      itab_u(iu)%iu3 = iw2_iu3	
      itab_u(iu)%iu4 = iw2_iu2
   elseif (iw2_iu2 == iu) then
      itab_u(iu)%iu3 = iw2_iu1	
      itab_u(iu)%iu4 = iw2_iu3
   else
      itab_u(iu)%iu3 = iw2_iu2	
      itab_u(iu)%iu4 = iw2_iu1
   endif

   iu1 = itab_u(iu)%iu1
   iu2 = itab_u(iu)%iu2
   iu3 = itab_u(iu)%iu3
   iu4 = itab_u(iu)%iu4

! Fill IW3 and DIRU1 for current U point
! This should work ok for iu1 = 1, but may bypass if desired

   if (itab_u(iu1)%iw1 == iw1) then
      itab_u(iu)%iw3 = itab_u(iu1)%iw2
      itab_u(iu)%diru1 = -1.
   else
      itab_u(iu)%iw3 = itab_u(iu1)%iw1
      itab_u(iu)%diru1 = 1.
   endif

! Fill IW4 and DIRU2 for current U point
! This should work ok for iu2 = 1, but may bypass if desired

   if (itab_u(iu2)%iw1 == iw1) then
      itab_u(iu)%iw4 = itab_u(iu2)%iw2
      itab_u(iu)%diru2 = -1.
   else
      itab_u(iu)%iw4 = itab_u(iu2)%iw1
      itab_u(iu)%diru2 = 1.
   endif

! Fill IW5 and DIRU3 for current U point
! This should work ok for iu3 = 1, but may bypass if desired

   if (itab_u(iu3)%iw1 == iw2) then
      itab_u(iu)%iw5 = itab_u(iu3)%iw2
      itab_u(iu)%diru3 = -1.
   else
      itab_u(iu)%iw5 = itab_u(iu3)%iw1
      itab_u(iu)%diru3 = 1.
   endif

! Fill IW6 and DIRU4 for current U point
! This should work ok for iu4 = 1, but may bypass if desired

   if (itab_u(iu4)%iw1 == iw2) then
      itab_u(iu)%iw6 = itab_u(iu4)%iw2
      itab_u(iu)%diru4 = -1.
   else
      itab_u(iu)%iw6 = itab_u(iu4)%iw1
      itab_u(iu)%diru4 = 1.
   endif

   iw3 = itab_u(iu)%iw3
   iw4 = itab_u(iu)%iw4
   iw5 = itab_u(iu)%iw5
   iw6 = itab_u(iu)%iw6

! Fill IU5 and IU6 for current U point
! This should work ok for iw3 = 1, but may bypass if desired

   if (iu1 == itab_w(iw3)%iu1) then
      itab_u(iu)%iu5 = itab_w(iw3)%iu2
      itab_u(iu)%iu6 = itab_w(iw3)%iu3
   elseif (iu1 == itab_w(iw3)%iu2) then
      itab_u(iu)%iu5 = itab_w(iw3)%iu3
      itab_u(iu)%iu6 = itab_w(iw3)%iu1
   else
      itab_u(iu)%iu5 = itab_w(iw3)%iu1
      itab_u(iu)%iu6 = itab_w(iw3)%iu2
   endif

! Fill IU7 and IU8 for current U point
! This should work ok for iw4 = 1, but may bypass if desired

   if (iu2 == itab_w(iw4)%iu1) then
      itab_u(iu)%iu7 = itab_w(iw4)%iu2
      itab_u(iu)%iu8 = itab_w(iw4)%iu3
   elseif (iu2 == itab_w(iw4)%iu2) then
      itab_u(iu)%iu7 = itab_w(iw4)%iu3
      itab_u(iu)%iu8 = itab_w(iw4)%iu1
   else
      itab_u(iu)%iu7 = itab_w(iw4)%iu1
      itab_u(iu)%iu8 = itab_w(iw4)%iu2
   endif

! Fill IU9 and IU10 for current U point
! This should work ok for iw5 = 1, but may bypass if desired

   if (iu3 == itab_w(iw5)%iu1) then
      itab_u(iu)%iu9  = itab_w(iw5)%iu3
      itab_u(iu)%iu10 = itab_w(iw5)%iu2
   elseif (iu3 == itab_w(iw5)%iu2) then
      itab_u(iu)%iu9  = itab_w(iw5)%iu1
      itab_u(iu)%iu10 = itab_w(iw5)%iu3
   else
      itab_u(iu)%iu9  = itab_w(iw5)%iu2
      itab_u(iu)%iu10 = itab_w(iw5)%iu1
   endif

! Fill IU11 and IU12 for current U point
! This should work ok for iw6 = 1, but may bypass if desired

   if (iu4 == itab_w(iw6)%iu1) then
      itab_u(iu)%iu11 = itab_w(iw6)%iu3
      itab_u(iu)%iu12 = itab_w(iw6)%iu2
   elseif (iu4 == itab_w(iw6)%iu2) then
      itab_u(iu)%iu11 = itab_w(iw6)%iu1
      itab_u(iu)%iu12 = itab_w(iw6)%iu3
   else
      itab_u(iu)%iu11 = itab_w(iw6)%iu2
      itab_u(iu)%iu12 = itab_w(iw6)%iu1
   endif

enddo  ! end loop over U points

! Fill U and W points for M points (do this as loop over U points)

do iu = 2,mua
   do j = 1,2

      if (j == 1) im = itab_u(iu)%im1
      if (j == 2) im = itab_u(iu)%im2

      iw1 = itab_u(iu)%iw1
      iw2 = itab_u(iu)%iw2

      if (itab_m(im)%ntpn == 0   .or.  &
         (iw1 == 1 .and. j == 1) .or.  &  ! This and next line added for walls
         (iw2 == 1 .and. j == 2)) then

         iunow = iu
         iu0 = 0
         ntpn = 0

         do while (iunow /= iu0)

            iu0 = iu
!!                  ntpn = ntpn + 1

            if (itab_u(iunow)%im1 == im) then
               if (itab_u(iunow)%iw2 > 1) then
                  ntpn = ntpn + 1

                  itab_m(im)%iu(ntpn) = iunow  ! NEW

                  itab_m(im)%iw(ntpn) = itab_u(iunow)%iw2
                  iunow = itab_u(iunow)%iu3
               else
                  iunow = iu0  ! this section added for walls
               endif
            else
               if (itab_u(iunow)%iw1 > 1) then
                  ntpn = ntpn + 1

                  itab_m(im)%iu(ntpn) = iunow  ! NEW

                  itab_m(im)%iw(ntpn) = itab_u(iunow)%iw1
                  iunow = itab_u(iunow)%iu2
               else
                  iunow = iu0  ! this section added for walls
               endif
            endif

            itab_m(im)%ntpn = ntpn

         enddo

! Define extra (wrap-around) iu value for itab_m            

         if (ntpn > 0) itab_m(im)%iu(ntpn+1) = iunow

      endif

   enddo
enddo

return
end subroutine tri_neighbors

!===============================================================================

subroutine matrix_3x3(a11,a21,a31,a12,a22,a32,a13,a23,a33,b1,b2,b3,x1,x2,x3)

implicit none

real, intent(in) :: a11,a21,a31,a12,a22,a32,a13,a23,a33,b1,b2,b3
real, intent(out) :: x1,x2,x3

integer :: i
real, dimension(4,3) :: abr

abr(1,1) = a11
abr(2,1) = a21
abr(3,1) = a31
abr(4,1) = b1

abr(1,2) = a12
abr(2,2) = a22
abr(3,2) = a32
abr(4,2) = b2

abr(1,3) = a13
abr(2,3) = a23
abr(3,3) = a33
abr(4,3) = b3

! Interchange rows if necessary so that first row has 
! largest (magnitude) element of first column

if (abs(abr(1,2)) > abs(abr(1,1))) then
   do i = 1,4
      call rchange(abr(i,1),abr(i,2))
   enddo
endif

if (abs(abr(1,3)) > abs(abr(1,1))) then
   do i = 1,4
      call rchange(abr(i,1),abr(i,3))
   enddo
endif

! Add -abr(1,2)/abr(1,1) times first row to second row and
! add -abr(1,3)/abr(1,1) times first row to third row.

do i = 2,4
   abr(i,2) = abr(i,2) - abr(1,2)/abr(1,1)*abr(i,1)
   abr(i,3) = abr(i,3) - abr(1,3)/abr(1,1)*abr(i,1)
enddo

! Interchange rows 2 and 3 if necessary so that second row
! has larger (magnitude) element of second column

if (abs(abr(2,3)) > abs(abr(2,2))) then
   do i = 2,4
      call rchange(abr(i,2),abr(i,3))
   enddo
endif

! Add -abr(2,3)/abr(2,2) times second row to third row.

do i = 3,4
   abr(i,3) = abr(i,3) - abr(2,3)/abr(2,2)*abr(i,2)
enddo

! Back substitution

x3 = abr(4,3) / abr(3,3)
x2 = (abr(4,2) - abr(3,2) * x3) / abr(2,2)
x1 = (abr(4,1) - abr(2,1) * x2 - abr(3,1) * x3) / abr(1,1)

return
end subroutine matrix_3x3

!===============================================================================

subroutine matrix_2x2(a11,a21,a12,a22,b1,b2,x1,x2)

implicit none

real, intent(in) :: a11,a21,a12,a22,b1,b2
real, intent(out) :: x1,x2

integer :: i
real, dimension(3,2) :: abr

abr(1,1) = a11
abr(2,1) = a21
abr(3,1) = b1

abr(1,2) = a12
abr(2,2) = a22
abr(3,2) = b2

! Interchange rows if necessary so that first row has 
! largest (magnitude) element of first column

if (abs(abr(1,2)) > abs(abr(1,1))) then
   do i = 1,3
      call rchange(abr(i,1),abr(i,2))
   enddo
endif

! Add -abr(1,2)/abr(1,1) times first row to second row

do i = 2,3
   abr(i,2) = abr(i,2) - abr(1,2)/abr(1,1)*abr(i,1)
enddo

! Back substitution

x2 = abr(3,2) / abr(2,2)
x1 = (abr(3,1) - abr(2,1) * x2) / abr(1,1)

return
end subroutine matrix_2x2

!===============================================================================

subroutine rchange(r1,r2)

implicit none

real, intent(inout) :: r1,r2
real :: c

c = r1
r1 = r2
r2 = c

return
end subroutine rchange

!===============================================================================

subroutine matinv3x3(a11,a21,a31,a12,a22,a32,a13,a23,a33  &
                    ,b11,b21,b31,b12,b22,b32,b13,b23,b33)

implicit none
 
real, intent(in)  :: a11,a21,a31,a12,a22,a32,a13,a23,a33
real, intent(out) :: b11,b21,b31,b12,b22,b32,b13,b23,b33

real :: det

det = a11*(a22*a33-a23*a32)  &
    + a21*(a32*a13-a12*a33)  &
    + a31*(a12*a23-a22*a13)
    
if (abs(det) < 1.e-12) stop 'stop singular matrix'    
    
b11 = (a22*a33-a32*a23)/det
b21 = (a31*a23-a21*a33)/det
b31 = (a21*a32-a31*a22)/det
b12 = (a32*a13-a12*a33)/det
b22 = (a11*a33-a31*a13)/det
b32 = (a31*a12-a11*a32)/det
b13 = (a12*a23-a22*a13)/det
b23 = (a21*a13-a11*a23)/det
b33 = (a11*a22-a21*a12)/det

return
end subroutine matinv3x3

!===============================================================================

subroutine matinv2x2(a11,a21,a12,a22,b11,b21,b12,b22)

implicit none

real, intent(in)  :: a11,a21,a12,a22
real, intent(out) :: b11,b21,b12,b22

real :: det

det = a11*a22 - a21*a12

if (abs(det) < 1.e-12) stop 'stop singular matrix'    

b11 =  a22/det
b21 = -a21/det
b12 = -a12/det
b22 =  a11/det

return
end subroutine matinv2x2

!===============================================================================

subroutine unit_normal(px,py,pz,qx,qy,qz,rx,ry,rz,vx,vy,vz)

implicit none

real, intent(in)  :: px,py,pz,qx,qy,qz,rx,ry,rz
real, intent(out) :: vx,vy,vz

real :: v

! Find components (vx,vy,vz) of unit normal vector to plane
! that passes through points (p,q,r)

! Un-normalized V = (PQ) X (PR):

vx = (qy - py) * (rz - pz) - (qz - pz) * (ry - py)
vy = (qz - pz) * (rx - px) - (qx - px) * (rz - pz)
vz = (qx - px) * (ry - py) - (qy - py) * (rx - px)

! Magnitude of V:

v = sqrt(vx * vx + vy * vy + vz * vz)

! Normalized components:

vx = vx / v
vy = vy / v
vz = vz / v

return
end subroutine unit_normal




