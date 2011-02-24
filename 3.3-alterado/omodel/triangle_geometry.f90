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
subroutine triangle_geometry()

use mem_ijtabs,  only: itab_m, itab_u, itab_w
use mem_grid,    only: nza, mma, mua, mwa, xeu, yeu, zeu, xem, yem, zem,  &
                       xew, yew, zew, unx, uny, unz, wnx, wny, wnz,       &
                       utx, uty, utz, glonw, glatw, dtu, dniu, arw0, glonm, glatm
use misc_coms,   only: io6, mdomain, centlat, centlon
use consts_coms, only: erad, erad2, piu180
use mem_nudge,   only: nudflag, xenudp, yenudp, zenudp, itab_nudp

implicit none

integer :: im,iu,iw
integer :: im1,im2,im3
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9,iu10,iu11,iu12
integer :: iw1,iw2,iw3,iw4,iw5,iw6
integer :: itpn
integer :: inudp,inudp1,inudp2,inudp3,jnudp,jjnudp
integer :: j,jj,jmax,jnudpmax,jmaxneg,jminpos

real :: expansion
real :: raxis
real :: hper

real :: xiw,yiw,ziw
real :: xij(6),yij(6),zij(6)
real :: x1,x2,x3,y1,y2,y3
real :: scalprod,scalprod_max
real :: vecjx,vecjy,vecjz
real :: vecmx,vecmy,vecmz
real :: vecmjx,vecmjy,vecmjz
real :: xi,yi,xj(6),yj(6)
real :: vecprodz,vecprodz_maxneg,vecprodz_minpos

real :: ef,x12,y12,z12,x34,y34,z34
real :: b11,b21,b31,b12,b22,b32,b13,b23,b33  &
       ,c11,c21,c31,c12,c22,c32,c13,c23,c33  &
       ,d11,d21,d31,d12,d22,d32,d13,d23,d33
real :: dnupgf
real :: s1, s2, s3

ef = 1.01  ! radial expansion factor (from earth center) for defining 
           ! auxiliary point for computing unit normal to U face

! Loop over all U points

do iu = 2,mua

! Fill global index (replaced later if this run is parallel)

   itab_u(iu)%iuglobe = iu

! M-point indices of two end points of U segment

   im1 = itab_u(iu)%im1
   im2 = itab_u(iu)%im2

! Fill U point coordinates from M point coordinates

   xeu(iu) = .5 * (xem(im1) + xem(im2))
   yeu(iu) = .5 * (yem(im1) + yem(im2))
   zeu(iu) = .5 * (zem(im1) + zem(im2))

! If mdomain <= 1, push U point coordinates out to earth radius

   if (mdomain <= 1) then

      expansion = erad / sqrt(xeu(iu) ** 2  &
                            + yeu(iu) ** 2  &
                            + zeu(iu) ** 2  )

      xeu(iu) = xeu(iu) * expansion
      yeu(iu) = yeu(iu) * expansion
      zeu(iu) = zeu(iu) * expansion

   endif

! Length of U side

   dtu(iu) = sqrt( (xem(im1) - xem(im2))**2  &
                 + (yem(im1) - yem(im2))**2  &
                 + (zem(im1) - zem(im2))**2  )
                 
! Convert to geodesic arc length

!   dtu(iu) = erad2 * asin(dtu(iu) / erad2)

! Fill U face unit normal vector components

   if (mdomain <= 1) then   ! Spherical geometry case

      call unit_normal(xem(im2),yem(im2),zem(im2)           &
                      ,xem(im1),yem(im1),zem(im1)           &
                      ,ef*xem(im1),ef*yem(im1),ef*zem(im1)  &
                      ,unx(iu),uny(iu),unz(iu))   

   else                     ! Cartesian geometry case

      call unit_normal(xem(im2),yem(im2),zem(im2)         &
                      ,xem(im1),yem(im1),zem(im1)         &
                      ,xem(im1),yem(im1),zem(im1) + 1.e3  &
                      ,unx(iu),uny(iu),unz(iu))   

   endif                

! Fill U face unit horizontal tangential vector components

   utx(iu) = (xem(im2) - xem(im1)) / dtu(iu)
   uty(iu) = (yem(im2) - yem(im1)) / dtu(iu)
   utz(iu) = (zem(im2) - zem(im1)) / dtu(iu)

enddo

! Loop over all W points

do iw = 2,mwa

! Fill global index (replaced later if this run is parallel)

   itab_w(iw)%iwglobe = iw

! Indices of 3 M points surrounding W point

   im1 = itab_w(iw)%im1
   im2 = itab_w(iw)%im2
   im3 = itab_w(iw)%im3

! Indices of 3 U points surrounding W point

   iu1 = itab_w(iw)%iu1                   
   iu2 = itab_w(iw)%iu2                   
   iu3 = itab_w(iw)%iu3                   

! Fill W face unit normal vector components

   call unit_normal(xem(im1),yem(im1),zem(im1)  &
                   ,xem(im2),yem(im2),zem(im2)  &
                   ,xem(im3),yem(im3),zem(im3)  &
                   ,wnx(iw) ,wny(iw) ,wnz(iw)   )

! Fill W point coordinates from coordinates of 3 M points

! OPTION 1: barycentric coordinates for IW point

   xew(iw) = (xem(im1) + xem(im2) + xem(im3)) / 3.
   yew(iw) = (yem(im1) + yem(im2) + yem(im3)) / 3.
   zew(iw) = (zem(im1) + zem(im2) + zem(im3)) / 3.

! OPTION 2: Circumcentric coordinates for IW point

!   xew(iw) = wnx(iw) * erad
!   yew(iw) = wny(iw) * erad
!   zew(iw) = wnz(iw) * erad

! If mdomain <= 1, push W point coordinates out to earth radius

   if (mdomain <= 1) then

      expansion = erad / sqrt(xew(iw) ** 2  &
                            + yew(iw) ** 2  &
                            + zew(iw) ** 2  )

      xew(iw) = xew(iw) * expansion
      yew(iw) = yew(iw) * expansion
      zew(iw) = zew(iw) * expansion

   endif

! Fill latitude and longitude of W point

   if (mdomain <= 1) then
      raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)

      glatw(iw) = atan2(zew(iw),raxis)   * piu180
      glonw(iw) = atan2(yew(iw),xew(iw)) * piu180
   else
      glatw(iw) = centlat(1)  ! want it this way?
      glonw(iw) = centlon(1)  ! want it this way?
   endif

! IW triangle area at sea level

   hper = .5 * (dtu(iu1) + dtu(iu2) + dtu(iu3))  ! half perimeter
         
   arw0(iw) = sqrt(hper * (hper - dtu(iu1))  &
                        * (hper - dtu(iu2))  &
                        * (hper - dtu(iu3))  )

enddo

! Loop over all M points

do im = 2,mma

! Latitude and longitude at M points

   if (mdomain <= 1) then
      raxis = sqrt(xem(im) ** 2 + yem(im) ** 2)  ! dist from earth axis
      glatm(im) = atan2(zem(im),raxis)   * piu180
      glonm(im) = atan2(yem(im),xem(im)) * piu180
   else
      glatm(im) = centlat(1)  ! want it this way?
      glonm(im) = centlon(1)  ! want it this way?
   endif

! Fill global index (replaced later if this run is parallel)

   itab_m(im)%imglobe = im

! Set area of dual cell surrounding M to zero prior to summation over triangles

   itab_m(im)%arm = 0.

! Loop over all U neighbors of M

   do itpn = 1,itab_m(im)%ntpn

      iu = itab_m(im)%iu(itpn)

      iw1 = itab_u(iu)%iw1
      iw2 = itab_u(iu)%iw2

! Contribution to dual-cell area around M point of triangle M-IW1-IW2

      s1 = sqrt( (xew(iw1) - xew(iw2))**2  &
               + (yew(iw1) - yew(iw2))**2  &
               + (zew(iw1) - zew(iw2))**2  )

      s2 = sqrt( (xem(im) - xew(iw2))**2  &
               + (yem(im) - yew(iw2))**2  &
               + (zem(im) - zew(iw2))**2  )
            
      s3 = sqrt( (xew(iw1) - xem(im))**2  &
               + (yew(iw1) - yem(im))**2  &
               + (zew(iw1) - zem(im))**2  )
            
      hper = .5 * (s1 + s2 + s3)  ! half perimeter of triangle
 
      itab_m(im)%arm = itab_m(im)%arm  &
                     + sqrt(hper * (hper - s1) * (hper - s2) * (hper - s3))

   enddo

enddo

! Loop over all U points

do iu = 2,mua

   im1 = itab_u(iu)%im1
   im2 = itab_u(iu)%im2

   iu1  = itab_u(iu)%iu1 
   iu2  = itab_u(iu)%iu2 
   iu3  = itab_u(iu)%iu3 
   iu4  = itab_u(iu)%iu4
   iu5  = itab_u(iu)%iu5
   iu6  = itab_u(iu)%iu6
   iu7  = itab_u(iu)%iu7
   iu8  = itab_u(iu)%iu8
   iu9  = itab_u(iu)%iu9
   iu10 = itab_u(iu)%iu10
   iu11 = itab_u(iu)%iu11
   iu12 = itab_u(iu)%iu12

   iw1 = itab_u(iu)%iw1
   iw2 = itab_u(iu)%iw2
   iw3 = itab_u(iu)%iw3
   iw4 = itab_u(iu)%iw4
   iw5 = itab_u(iu)%iw5
   iw6 = itab_u(iu)%iw6

! Project neighbor U & W unit vectors onto U horizontal tangential vector

   if (iw1 > 1)                                       &
   call matrix_3x3(unx(iu1),unx(iu2),wnx(iw1)         &
                  ,uny(iu1),uny(iu2),wny(iw1)         &
                  ,unz(iu1),unz(iu2),wnz(iw1)         &
                  ,utx(iu),uty(iu),utz(iu)            &
                  ,b11,b12,b13                        )

   if (iw2 > 1)                                       &
   call matrix_3x3(unx(iu3),unx(iu4),wnx(iw2)         &
                  ,uny(iu3),uny(iu4),wny(iw2)         &
                  ,unz(iu3),unz(iu4),wnz(iw2)         &
                  ,utx(iu),uty(iu),utz(iu)            &
                  ,c11,c12,c13                        )

   if (iw1 > 1 .and. iw2 > 1) then
      itab_u(iu)%tuu1 = .5 * b11
      itab_u(iu)%tuu2 = .5 * b12
      itab_u(iu)%tuu3 = .5 * c11
      itab_u(iu)%tuu4 = .5 * c12
   elseif (iw1 > 1) then
      itab_u(iu)%tuu1 = b11
      itab_u(iu)%tuu2 = b12
   elseif (iw2 > 1) then
      itab_u(iu)%tuu3 = c11
      itab_u(iu)%tuu4 = c12
   endif

! Skip this U point if iw1 < 2 or iw2 < 2

   if (iw1 < 2 .or. iw2 < 2) cycle
   
   itab_u(iu)%mrlu = max(itab_w(iw1)%mrlw,itab_w(iw2)%mrlw)

! Project neighbor U & W unit vectors onto U unit vector

   if (iw3 > 1)                                &
   call matrix_3x3(unx(iu5),unx(iu6),wnx(iw3)  &
                  ,uny(iu5),uny(iu6),wny(iw3)  &
                  ,unz(iu5),unz(iu6),wnz(iw3)  &
                  ,unx(iu),uny(iu),unz(iu)     &
                  ,itab_u(iu)%fuu5             &
                  ,itab_u(iu)%fuu6             &
                  ,itab_u(iu)%fuw3             )

   if (iw4 > 1)                                &
   call matrix_3x3(unx(iu7),unx(iu8),wnx(iw4)  &
                  ,uny(iu7),uny(iu8),wny(iw4)  &
                  ,unz(iu7),unz(iu8),wnz(iw4)  &
                  ,unx(iu),uny(iu),unz(iu)     &
                  ,itab_u(iu)%fuu7             &
                  ,itab_u(iu)%fuu8             &
                  ,itab_u(iu)%fuw4             )

   if (iw5 > 1)                                 &
   call matrix_3x3(unx(iu9),unx(iu10),wnx(iw5)  &
                  ,uny(iu9),uny(iu10),wny(iw5)  &
                  ,unz(iu9),unz(iu10),wnz(iw5)  &
                  ,unx(iu),uny(iu),unz(iu)      &
                  ,itab_u(iu)%fuu9              &
                  ,itab_u(iu)%fuu10             &
                  ,itab_u(iu)%fuw5              )

   if (iw6 > 1)                                  &
   call matrix_3x3(unx(iu11),unx(iu12),wnx(iw6)  &
                  ,uny(iu11),uny(iu12),wny(iw6)  &
                  ,unz(iu11),unz(iu12),wnz(iw6)  &
                  ,unx(iu),uny(iu),unz(iu)       &
                  ,itab_u(iu)%fuu11              &
                  ,itab_u(iu)%fuu12              &
                  ,itab_u(iu)%fuw6               )

! Divide fuw3, fuw4, fuw5, and fuw6 by 2 for use with two W points

   if (iw3 > 1) itab_u(iu)%fuw3 = .5 * itab_u(iu)%fuw3
   if (iw4 > 1) itab_u(iu)%fuw4 = .5 * itab_u(iu)%fuw4
   if (iw5 > 1) itab_u(iu)%fuw5 = .5 * itab_u(iu)%fuw5
   if (iw6 > 1) itab_u(iu)%fuw6 = .5 * itab_u(iu)%fuw6

! Find earth velocity coefficients for Coriolis force

   call matinv3x3(unx(iu1),uny(iu1),unz(iu1)          &
                 ,unx(iu2),uny(iu2),unz(iu2)          &
                 ,wnx(iw1),wny(iw1),wnz(iw1)          &
                 ,b11,b21,b31,b12,b22,b32,b13,b23,b33 )

   call matinv3x3(unx(iu3),uny(iu3),unz(iu3)          &
                 ,unx(iu4),uny(iu4),unz(iu4)          &
                 ,wnx(iw2),wny(iw2),wnz(iw2)          &
                 ,c11,c21,c31,c12,c22,c32,c13,c23,c33 )

   itab_u(iu)%vxu1_u = b11
   itab_u(iu)%vxu2_u = b21
   itab_u(iu)%vxw1_u = b31 / 2.

   itab_u(iu)%vyu1_u = b12
   itab_u(iu)%vyu2_u = b22
   itab_u(iu)%vyw1_u = b32 / 2.

   itab_u(iu)%vxu3_u = c11
   itab_u(iu)%vxu4_u = c21
   itab_u(iu)%vxw2_u = c31 / 2.

   itab_u(iu)%vyu3_u = c12
   itab_u(iu)%vyu4_u = c22
   itab_u(iu)%vyw2_u = c32 / 2.

! U point PGF coefficients for full 6-point stencil

   x12 = xew(iw2) - xew(iw1)
   y12 = yew(iw2) - yew(iw1)
   z12 = zew(iw2) - zew(iw1)

   dniu(iu) = 1. / (unx(iu) * x12 + uny(iu) * y12 + unz(iu) * z12)

   dnupgf = (arw0(iw1) + arw0(iw2)) / dtu(iu)


!!!!!!!!! ALTERNATE FORM
   dnupgf = 1.
!!!!!!!!!!!!!!!!!!!!!!!!   

   if (iw3 > 1 .and. iw4 > 1 .and. iw5 > 1 .and. iw6 > 1) then

      x34 = xew(iw3) + xew(iw5) - xew(iw4) - xew(iw6)
      y34 = yew(iw3) + yew(iw5) - yew(iw4) - yew(iw6)
      z34 = zew(iw3) + zew(iw5) - zew(iw4) - zew(iw6)

      if (mdomain <= 1) then ! Spherical geometry case

         call matinv3x3(x12,y12,z12  &
                       ,x34,y34,z34  &
                       ,y12*z34-z12*y34,z12*x34-x12*z34,x12*y34-y12*x34   &
                       ,b11,b21,b31,b12,b22,b32,b13,b23,b33)

      else                   ! Cartesian geometry case

         call matinv2x2(x12,y12,x34,y34,b11,b21,b12,b22)

         b13 = 0.
         b23 = 0.
      endif
                 
      itab_u(iu)%pgc12 = dnupgf * (unx(iu) * b11 + uny(iu) * b12 + unz(iu) * b13)
      itab_u(iu)%pgc45 = dnupgf * (unx(iu) * b21 + uny(iu) * b22 + unz(iu) * b23)
      itab_u(iu)%pgc63 = itab_u(iu)%pgc45
      
   endif

! U point PGF coefficients for (1-2-4-5) 4-point stencil

   if (iw4 > 1 .and. iw5 > 1) then

      x34 = xew(iw5) - xew(iw4)
      y34 = yew(iw5) - yew(iw4)
      z34 = zew(iw5) - zew(iw4)

      if (mdomain <= 1) then ! Spherical geometry case

         call matinv3x3(x12,y12,z12  &
                       ,x34,y34,z34  &
                       ,y12*z34-z12*y34,z12*x34-x12*z34,x12*y34-y12*x34   &
                       ,b11,b21,b31,b12,b22,b32,b13,b23,b33)

      else                   ! Cartesian geometry case

         call matinv2x2(x12,y12,x34,y34,b11,b21,b12,b22)

         b13 = 0.
         b23 = 0.
      endif

      itab_u(iu)%pgc12b = dnupgf * (unx(iu) * b11 + uny(iu) * b12 + unz(iu) * b13)
      itab_u(iu)%pgc45b = dnupgf * (unx(iu) * b21 + uny(iu) * b22 + unz(iu) * b23)

   endif

! U point PGF coefficients for (1-2-6-3) 4-point stencil

   if (iw3 > 1 .and. iw6 > 1) then

      x34 = xew(iw3) - xew(iw6)
      y34 = yew(iw3) - yew(iw6)
      z34 = zew(iw3) - zew(iw6)

      if (mdomain <= 1) then ! Spherical geometry case

         call matinv3x3(x12,y12,z12  &
                       ,x34,y34,z34  &
                       ,y12*z34-z12*y34,z12*x34-x12*z34,x12*y34-y12*x34   &
                       ,b11,b21,b31,b12,b22,b32,b13,b23,b33)

      else                   ! Cartesian geometry case

         call matinv2x2(x12,y12,x34,y34,b11,b21,b12,b22)

         b13 = 0.
         b23 = 0.
      endif
                 
      itab_u(iu)%pgc12c = dnupgf * (unx(iu) * b11 + uny(iu) * b12 + unz(iu) * b13)
      itab_u(iu)%pgc63c = dnupgf * (unx(iu) * b21 + uny(iu) * b22 + unz(iu) * b23)
      
   endif

! U point PGF coefficients for (1-2) 2-point stencil

   if (mdomain <= 1) then ! Spherical geometry
      itab_u(iu)%pgc12d = dnupgf / (unx(iu) * x12 + uny(iu) * y12 + unz(iu) * z12)
   else                   ! Cartesian geometry
      itab_u(iu)%pgc12d = dnupgf / (unx(iu) * x12 + uny(iu) * y12)
   endif

enddo

! Loop over all W points

do iw = 2,mwa

! Indices of 3 M points, 9 U points, and 3 W points surrounding W point

   im1 = itab_w(iw)%im1                   
   im2 = itab_w(iw)%im2                   
   im3 = itab_w(iw)%im3                   

   iu1 = itab_w(iw)%iu1                   
   iu2 = itab_w(iw)%iu2                   
   iu3 = itab_w(iw)%iu3                   
   iu4 = itab_w(iw)%iu4                   
   iu5 = itab_w(iw)%iu5                   
   iu6 = itab_w(iw)%iu6                   
   iu7 = itab_w(iw)%iu7                   
   iu8 = itab_w(iw)%iu8                   
   iu9 = itab_w(iw)%iu9

   iw1 = itab_w(iw)%iw1                   
   iw2 = itab_w(iw)%iw2                   
   iw3 = itab_w(iw)%iw3                   

! Project neighbor U & W unit vectors onto W unit vector

   if (iw1 > 1)                                &
   call matrix_3x3(unx(iu4),unx(iu5),wnx(iw1)  &
                  ,uny(iu4),uny(iu5),wny(iw1)  &
                  ,unz(iu4),unz(iu5),wnz(iw1)  &
                  ,wnx(iw),wny(iw),wnz(iw)     &
                  ,itab_w(iw)%fwu4             &
                  ,itab_w(iw)%fwu5             &
                  ,itab_w(iw)%fww1             )

   if (iw2 > 1)                                &
   call matrix_3x3(unx(iu6),unx(iu7),wnx(iw2)  &
                  ,uny(iu6),uny(iu7),wny(iw2)  &
                  ,unz(iu6),unz(iu7),wnz(iw2)  &
                  ,wnx(iw),wny(iw),wnz(iw)     &
                  ,itab_w(iw)%fwu6             &
                  ,itab_w(iw)%fwu7             &
                  ,itab_w(iw)%fww2             )

   if (iw3 > 1)                                &
   call matrix_3x3(unx(iu8),unx(iu9),wnx(iw3)  &
                  ,uny(iu8),uny(iu9),wny(iw3)  &
                  ,unz(iu8),unz(iu9),wnz(iw3)  &
                  ,wnx(iw),wny(iw),wnz(iw)     &
                  ,itab_w(iw)%fwu8             &
                  ,itab_w(iw)%fwu9             &
                  ,itab_w(iw)%fww3             )

! Divide fwu4, fwu5, fwu6, fwu7, fwu8, and fwu9 by 2 for use with two W points
                  
   if (iw1 > 1) itab_w(iw)%fwu4 = .5 * itab_w(iw)%fwu4
   if (iw1 > 1) itab_w(iw)%fwu5 = .5 * itab_w(iw)%fwu5
   if (iw2 > 1) itab_w(iw)%fwu6 = .5 * itab_w(iw)%fwu6
   if (iw2 > 1) itab_w(iw)%fwu7 = .5 * itab_w(iw)%fwu7
   if (iw3 > 1) itab_w(iw)%fwu8 = .5 * itab_w(iw)%fwu8
   if (iw3 > 1) itab_w(iw)%fwu9 = .5 * itab_w(iw)%fwu9

! Find horizontal unit vector components for each pair of U's in this W

   call matinv3x3(unx(iu1),uny(iu1),unz(iu1)  &
                 ,unx(iu2),uny(iu2),unz(iu2)  &
                 ,wnx(iw) ,wny(iw) ,wnz(iw)   &
                 ,b11,b21,b31,b12,b22,b32,b13,b23,b33)
                 
   call matinv3x3(unx(iu2),uny(iu2),unz(iu2)  &
                 ,unx(iu3),uny(iu3),unz(iu3)  &
                 ,wnx(iw) ,wny(iw) ,wnz(iw)   &
                 ,c11,c21,c31,c12,c22,c32,c13,c23,c33)

   call matinv3x3(unx(iu3),uny(iu3),unz(iu3)  &
                 ,unx(iu1),uny(iu1),unz(iu1)  &
                 ,wnx(iw) ,wny(iw) ,wnz(iw)   &
                 ,d11,d21,d31,d12,d22,d32,d13,d23,d33)

! Coefficients for 3D wind (for HORIZONTAL wind, don't use tri*4)

   itab_w(iw)%vxu1 = (b11 + d21) / 3.
   itab_w(iw)%vxu2 = (b21 + c11) / 3.
   itab_w(iw)%vxu3 = (c21 + d11) / 3.
   itab_w(iw)%vxw  = (b31 + c31 + d31) / 6.  ! Includes div by 2 for 2 W's

   itab_w(iw)%vyu1 = (b12 + d22) / 3.
   itab_w(iw)%vyu2 = (b22 + c12) / 3.
   itab_w(iw)%vyu3 = (c22 + d12) / 3.
   itab_w(iw)%vyw  = (b32 + c32 + d32) / 6.

   itab_w(iw)%vzu1 = (b13 + d23) / 3.
   itab_w(iw)%vzu2 = (b23 + c13) / 3.
   itab_w(iw)%vzu3 = (c23 + d13) / 3.
   itab_w(iw)%vzw  = (b33 + c33 + d33) / 6.

!!!!!!!!!!!! special (new form that omits umc(k,*) from its own fcor)
   itab_w(iw)%vxu1_w = (b11 + d21) / 3. ! Includes div by 3 for 3 tri vertices
   itab_w(iw)%vxu2_w = (b21 + c11) / 3. ! Includes div by 3 for 3 tri vertices
   itab_w(iw)%vxu3_w = (c21 + d11) / 3. ! Includes div by 3 for 3 tri vertices

   itab_w(iw)%vyu1_w = (b12 + d22) / 3. ! Includes div by 3 for 3 tri vertices
   itab_w(iw)%vyu2_w = (b22 + c12) / 3. ! Includes div by 3 for 3 tri vertices
   itab_w(iw)%vyu3_w = (c22 + d12) / 3. ! Includes div by 3 for 3 tri vertices
!!!!!!!!!!!! end special


!!!!!!! SAMPLE ONLY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! (xe,ye,ze) components of 3D wind

!   um = vxu1 * umc(k,iu1) + vxu2 * umc(k,iu2) + vxu3 * umc(k,iu3)  &
!      + vxw * (wmc(k,iw) + wmc(k+1,iw))
!   vm = vyu1 * umc(k,iu1) + vyu2 * umc(k,iu2) + vyu3 * umc(k,iu3)  &
!      + vyw * (wmc(k,iw) + wmc(k+1,iw))
!   wm = vzu1 * umc(k,iu1) + vzu2 * umc(k,iu2) + vzu3 * umc(k,iu3)  &
!      + vzw * (wmc(k,iw) + wmc(k+1,iw))

! (xe,ye,ze) components of HORIZONTAL wind

!   um = vxu1 * umc(k,iu1) + vxu2 * umc(k,iu2) + vxu3 * umc(k,iu3)
!   vm = vyu1 * umc(k,iu1) + vyu2 * umc(k,iu2) + vyu3 * umc(k,iu3)
!   wm = vzu1 * umc(k,iu1) + vzu2 * umc(k,iu2) + vzu3 * umc(k,iu3)

!!!!!!!!!!! END SAMPLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

enddo

! If using data assimilation with polygon nudging points, assign nudging
! point indices and coefficients to all W points

if (nudflag > 0) then

! Loop over all W points

   do iw = 2,mwa

! Get primary nudging point index for current W point and compute its
! x,y components on a polar stereographic plane tangent at W point
! (W point is at 0,0)

      inudp = itab_w(iw)%inudp(1)

      call e_ps(xenudp(inudp),yenudp(inudp),zenudp(inudp)  &
               ,glatw(iw),glonw(iw),xi,yi)

! Initialize vecprodz_minpos and vecprodz_maxneg

      vecprodz_maxneg = -1.e15
      vecprodz_minpos =  1.e15

! Loop through nearest polygon neighbors (j, jnudp) of nudging point INUDP

      do j = 1,6

! Get nudging point index (jnudp) for current polygon neighbor of inudp.
! Skip this j if jnudp < 2

         jnudp = itab_nudp(inudp,j)

         if (jnudp < 2) cycle

! Compute x,y components of jnudp polygon center on a polar stereographic 
! plane tangent at W point

         call e_ps(xenudp(jnudp),yenudp(jnudp),zenudp(jnudp)  &
            ,glatw(iw),glonw(iw),xj(j),yj(j))

! Compute z component (in polar stereographic space) of vector product 
! of the vector from inudp to iw and the vector from inudp to jnudp.

         vecprodz = (0. - xi) * (yj(j) - yi) - (0. - yi) * (xj(j) - xi)
         
! Compute scalar product of the vector from inudp to iw and the vector from
! inudp to jnudp in polar stereographic space.

         scalprod = (0. - xi) * (xj(j) - xi) + (0. - yi) * (yj(j) - yi) 

! Check whether scalar product is positive for current J point.  If so, 
! J point is a candidate for the nudging triad for IW.

         if (scalprod > 0.) then

! Identify maximum negative vecprodz among all jnudp polygon neighbors of
! inudp.  This jnudp will be second point of nudging triad for IW         

            if (vecprodz < 0. .and. vecprodz > vecprodz_maxneg) then
               vecprodz_maxneg = vecprodz
               jmaxneg = j
               itab_w(iw)%inudp(2) = jnudp
            endif

! Identify minimum positive vecprodz among all jnudp polygon neighbors of
! inudp.  This jnudp will be third point of nudging triad for IW         

            if (vecprodz >= 0. .and. vecprodz < vecprodz_minpos) then
               vecprodz_minpos = vecprodz
               jminpos = j
               itab_w(iw)%inudp(3) = jnudp
            endif
         endif

      enddo

! Lastly, fill 3 nudging weight coefficients for this W point.
! Weights are computed in 2_d polar stereographic space.

! Invert matrix of coordinates

      call matinv3x3(1.,xi,yi  &
                    ,1.,xj(jmaxneg),yj(jmaxneg)  &
                    ,1.,xj(jminpos),yj(jminpos)  &
                 ,b11,b21,b31,b12,b22,b32,b13,b23,b33)

! Assign coefficients

      itab_w(iw)%fnudp(1) = b11
      itab_w(iw)%fnudp(2) = b21
      itab_w(iw)%fnudp(3) = b31

   enddo

endif

return
end subroutine triangle_geometry

