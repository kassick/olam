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
subroutine spring_dynamics()

use mem_ijtabs,  only: itab_m, itab_u, itab_w
use mem_grid,    only: nma, nua, xem, yem, zem, impent, mrows
use consts_coms, only: pi2, erad, erador5
use misc_coms,   only: io6, nxp

implicit none

integer, parameter :: niter = 2000
real, parameter :: relax = .04, beta = 1.

! Automatic arrays

real :: dxem(nma)
real :: dyem(nma)
real :: dzem(nma)

integer :: iu,im1,im2,im,iter,ipent,iw1,iw2,mrow1,mrow2

real :: dist00,dist0,dist,expansion,frac_change,dx2,dy2,dz2

! special
!RETURN
! end special



! Compute mean length of U segments for global grid

dist00 = beta * pi2 * erad / (5. * real(nxp))

! Main iteration loop 

do iter = 1,niter

! Loop over all M points

   do im = 2,nma

! Initialize displacement sums to zero

      dxem(im) = 0.
      dyem(im) = 0.
      dzem(im) = 0.

   enddo

! Loop over all U points

   do iu = 2,nua

      im1 = itab_u(iu)%im1
      im2 = itab_u(iu)%im2

! Compute current distance between IM1 and IM2

      dist = sqrt((xem(im1) - xem(im2)) ** 2  &
                + (yem(im1) - yem(im2)) ** 2  &
                + (zem(im1) - zem(im2)) ** 2)

! Compute fractional change to dist that gives dist0

!-------------------------------------------------------------------------
!  SPECIAL - VARIABLE EQUILIBRIUM SPRING DISTANCE
!-------------------------------------------------------------------------

! Distance for any MRL value

      dist0 = dist00 * 2. ** real(1 - itab_u(iu)%mrlu)

! Reduced distance just outside MRL boundary

      iw1 = itab_u(iu)%iw1
      iw2 = itab_u(iu)%iw2

      mrow1 = itab_w(iw1)%mrow
      mrow2 = itab_w(iw2)%mrow

      if (mrow1 > 0 .and. mrow1 <= mrows .and.  &
          mrow2 > 0 .and. mrow2 <= mrows) then

         dist0 = dist0 * (.5 + .25/real(mrows) * (mrow1 + mrow2 - 1))      

! TESTS...
!         if     (mrow1 + mrow2 == 2) then
!            dist0 = dist0 * .52      
!         elseif (mrow1 + mrow2 == 3) then
!            dist0 = dist0 * .55      
!         elseif (mrow1 + mrow2 == 4) then
!            dist0 = dist0 * .60      
!         elseif (mrow1 + mrow2 == 5) then
!            dist0 = dist0 * .66      
!         elseif (mrow1 + mrow2 == 6) then
!            dist0 = dist0 * .73      
!         elseif (mrow1 + mrow2 == 7) then
!            dist0 = dist0 * .81      
!         elseif (mrow1 + mrow2 == 8) then
!            dist0 = dist0 * .90      
!         endif
         
      endif
      
!-------------------------------------------------------------------------
!  END SPECIAL
!-------------------------------------------------------------------------



      frac_change = (dist0 - dist) / dist

! Compute components of displacement that gives dist0

      dx2 = (xem(im2) - xem(im1)) * frac_change
      dy2 = (yem(im2) - yem(im1)) * frac_change
      dz2 = (zem(im2) - zem(im1)) * frac_change

! Add components of displacement to displacement of both M points

      dxem(im1) = dxem(im1) - dx2
      dyem(im1) = dyem(im1) - dy2
      dzem(im1) = dzem(im1) - dz2

      dxem(im2) = dxem(im2) + dx2
      dyem(im2) = dyem(im2) + dy2
      dzem(im2) = dzem(im2) + dz2

   enddo

! Loop over all M points

   do im = 2,nma

! For now, prevent either polar M point from moving

      if (im == impent(1 )) cycle
      if (im == impent(12)) cycle

! For preventing all pentagonal points from moving:
!     if (any(im == impent(1:12)) cycle

! Apply fraction of displacement to coordinates of M points

      xem(im) = xem(im) + relax * dxem(im)
      yem(im) = yem(im) + relax * dyem(im)
      zem(im) = zem(im) + relax * dzem(im)

! Push M point coordinates out to earth radius

      expansion = erad / sqrt(xem(im) ** 2 + yem(im) ** 2 + zem(im) ** 2)

      xem(im) = xem(im) * expansion
      yem(im) = yem(im) * expansion
      zem(im) = zem(im) * expansion

   enddo
enddo

return
end subroutine spring_dynamics

!===============================================================================

subroutine pcvt()

! Determine the coordinates of the centroids of all Voronoi cells of the OLAM
! dual grid.
   
use mem_ijtabs,  only: itab_m, itab_w
use mem_grid,    only: nma, nwa, xem, yem, zem, xew, yew, zew, wnx, wny, wnz
use consts_coms, only: erad, piu180
use misc_coms,   only: io6, mdomain

implicit none

integer, parameter :: niter = 2000
real, parameter :: frac = .1

integer :: i,ip,iw,im,im1,im2,im3,iter,ntpn

real :: xw(7),yw(7)
real :: raxis,glatm,glonm,area,xc,yc,xec,yec,zec,expansion

! Main iteration loop

! special
RETURN
! end special

do iter = 1,niter

! Loop over all W points

   do iw = 2,nwa

! Indices of 3 M points surrounding W point

      im1 = itab_w(iw)%im1
      im2 = itab_w(iw)%im2
      im3 = itab_w(iw)%im3

! Fill W point coordinates from coordinates of 3 M points

! OPTION 1: barycentric coordinates for IW point

      xew(iw) = (xem(im1) + xem(im2) + xem(im3)) / 3.
      yew(iw) = (yem(im1) + yem(im2) + yem(im3)) / 3.
      zew(iw) = (zem(im1) + zem(im2) + zem(im3)) / 3.

! OPTION 2: Circumcentric coordinates for IW point
! Fill W face unit normal vector components

      call unit_normal(xem(im1),yem(im1),zem(im1)  &
                      ,xem(im2),yem(im2),zem(im2)  &
                      ,xem(im3),yem(im3),zem(im3)  &
                      ,wnx(iw) ,wny(iw) ,wnz(iw)   )


!      xew(iw) = wnx(iw) * erad
!      yew(iw) = wny(iw) * erad
!      zew(iw) = wnz(iw) * erad

! If mdomain <= 1, push W point coordinates out to earth radius

      if (mdomain <= 1) then

         expansion = erad / sqrt(xew(iw) ** 2  &
                               + yew(iw) ** 2  &
                               + zew(iw) ** 2  )

         xew(iw) = xew(iw) * expansion
         yew(iw) = yew(iw) * expansion
         zew(iw) = zew(iw) * expansion

      endif

   enddo

! Loop over all M points

   do im = 2,nma

! Determine latitude and longitude of M point

      raxis = sqrt(xem(im) ** 2 + yem(im) ** 2)  ! dist from earth axis

      glatm = atan2(zem(im),raxis)   * piu180
      glonm = atan2(yem(im),xem(im)) * piu180

! Loop over all W points adjacent to current M point

      ntpn = itab_m(im)%ntpn

      do i = 1,ntpn
         iw = itab_m(im)%iw(i)

! Determine local PS coordinates of current W point

         call e_ps(xew(iw),yew(iw),zew(iw),glatm,glonm,xw(i),yw(i))

      enddo
   
! Compute Voronoi cell area for current M point

      area = 0.

      do i = 1,ntpn
         ip = i + 1
         if (ip > ntpn) ip = 1
   
         area = area + .5 * (xw(i) * yw(ip) - xw(ip) * yw(i))
      enddo
   
! Compute Voronoi cell centroid (in local PS coordinates) for current M point

     xc = 0.
     yc = 0.

     do i = 1,ntpn
        ip = i + 1
        if (ip > ntpn) ip = 1

         xc = xc + (xw(i) + xw(ip)) * (xw(i) * yw(ip) - xw(ip) * yw(i))
         yc = yc + (yw(i) + yw(ip)) * (xw(i) * yw(ip) - xw(ip) * yw(i))
     enddo

     xc = xc / (6. * area)
     yc = yc / (6. * area)

! Transform local PS coordinates of centroid to Earth coordinates

      call ps_e(xec,yec,zec,glatm,glonm,xc,yc)

! Move current M point (fractionally) to centroid location

      xem(im) = (1. - frac) * xem(im) + frac * xec
      yem(im) = (1. - frac) * yem(im) + frac * yec
      zem(im) = (1. - frac) * zem(im) + frac * zec

! Adjust M point coordinates to earth radius in case required

      expansion = erad / sqrt(xem(im) ** 2 + yem(im) ** 2 + zem(im) ** 2)

      xem(im) = xem(im) * expansion
      yem(im) = yem(im) * expansion
      zem(im) = zem(im) * expansion

   enddo

enddo

return
end subroutine pcvt
