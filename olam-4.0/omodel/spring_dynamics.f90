!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================
subroutine spring_dynamics(ngr,nconcave)

use mem_ijtabs,  only: itab_md, itab_ud, itab_wd
use mem_grid,    only: nma, nua, xem, yem, zem, impent, mrows
use consts_coms, only: pi2, erad, erador5
use misc_coms,   only: io6, nxp, ngrids
use oname_coms,  only: nl

!$ use omp_lib

implicit none

integer, intent(in) :: ngr,nconcave

integer :: niter
real, parameter :: relax = .04, beta = 1.

! Automatic arrays

real :: dxem(nma)
real :: dyem(nma)
real :: dzem(nma)

integer :: iu,im1,im2,im,iter,ipent,iw1,iw2,mrow1,mrow2,j,iw,npoly
integer :: iu1,iu2,iu3,iu4
integer :: im3,im4

integer :: iskip
real :: xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2

real :: dist1,dist2,dist3,dist4
real :: cosphi3, cosphi4
real :: coslim3, coslim4, ratio

real :: dist00,dist0,dist,expansion,frac_change
real :: dx2,dy2,dz2

! Array MOVEM is used to flag selected M points that will be allowed to move
! in the spring dynamics procedure.  If MOVEM = 1 the point is allowed to move,
! and if MOVEM = 0 the point remains stationary.  This selective movement is
! optional for static mesh refinements, but it will be a requirement for 
! dynamically adaptive refinements.  

integer :: movem(nma)

! special
!RETURN
! end special

if (nl%moveall(ngr) == 0 .and. ngr > 1 .and. nconcave == 3) then

! Set MOVEM flag = 1 only for M points that have W neighbors with MROW values
! of -3, -2, -1, or 1.

   movem(:) = 0

   do im = 2,nma
      npoly = itab_md(im)%npoly
         
      do j = 1,npoly
         iw = itab_md(im)%iw(j)

         if (itab_wd(iw)%mrow == -3 .or. &
              itab_wd(iw)%mrow == -2 .or. &
              itab_wd(iw)%mrow == -1 .or. &
              itab_wd(iw)%mrow == 1) movem(im) = 1
      enddo
   enddo

else

! The default is to set MOVEM flag = 1 for all points

   movem(:) = 1

endif

! Set value for niter based on ngr and ngrids

if (ngr > 1 .and. ngr < ngrids) then
   niter = 2000
else
   niter = 2000
endif

! Compute mean length of U segments for global grid

dist00 = beta * pi2 * erad / (5. * real(nxp))

! Main iteration loop 

do iter = 1,niter

! Loop over all M points

!$omp parallel

!$omp do
   do im = 2,nma

! Initialize displacement sums to zero

      dxem(im) = 0.
      dyem(im) = 0.
      dzem(im) = 0.

   enddo
!$omp end do

! Loop over all U points

!$omp do private (im1,im2,dist,dist0,iw1,iw2,mrow1,mrow2,iu1,iu3,im3,im4, &
!$omp             dist1,dist2,dist3,dist4,cosphi3,cosphi4, &
!$omp             ratio,frac_change,dx2,dy2,dz2)
   do iu = 2,nua

      im1 = itab_ud(iu)%im(1)
      im2 = itab_ud(iu)%im(2)

      if (ngr > 1 .and. movem(im1) == 0 .and. movem(im2) == 0) cycle

! Compute current distance between IM1 and IM2

      dist = sqrt((xem(im1) - xem(im2)) ** 2  &
                + (yem(im1) - yem(im2)) ** 2  &
                + (zem(im1) - zem(im2)) ** 2)

! Compute target distance for any MRL value

      dist0 = dist00 * 2. ** real(1 - itab_ud(iu)%mrlu)

! Modified distance in MRL border zone

      iw1 = itab_ud(iu)%iw(1)
      iw2 = itab_ud(iu)%iw(2)

      mrow1 = itab_wd(iw1)%mrow
      mrow2 = itab_wd(iw2)%mrow

      if (nconcave == 3 .and. mrow1 /= 0 .and. mrow2 /= 0) then

         if     (mrow1 + mrow2 == -4) then
            dist0 = dist0 *  7. / 6.  !* .90
         elseif (mrow1 + mrow2 == -3) then
            dist0 = dist0 *  8. / 6.  !* .90
         elseif (mrow1 + mrow2 == -2) then
            dist0 = dist0 *  9. / 6.  !* .90
         elseif (mrow1 + mrow2 ==  0) then
            dist0 = dist0 * 10. / 6.  !* .90
         elseif (mrow1 + mrow2 ==  2) then
            dist0 = dist0 * 11. / 12. ! * .90
         endif

      elseif (mrow1 > 0 .and. mrow1 <= mrows .and.  &
              mrow2 > 0 .and. mrow2 <= mrows) then

         dist0 = dist0 * (.5 + .25/real(mrows) * (mrow1 + mrow2 - 1))      

      endif ! nconcave

! Adjustment of dist0 based on opposite angles of triangles

      iu1 = itab_ud(iu)%iu(1)
      iu3 = itab_ud(iu)%iu(3)

! Determine im3 and im4 points for current IU

      if (itab_ud(iu1)%im(1) == im1) then
         im3 = itab_ud(iu1)%im(2)
      else
         im3 = itab_ud(iu1)%im(1)
      endif

      if (itab_ud(iu3)%im(1) == im1) then
         im4 = itab_ud(iu3)%im(2)
      else
         im4 = itab_ud(iu3)%im(1)
      endif

! Compute current IU1 length (distance between IM1 and IM3)

      dist1 = sqrt((xem(im1) - xem(im3)) ** 2  &
                 + (yem(im1) - yem(im3)) ** 2  &
                 + (zem(im1) - zem(im3)) ** 2)

! Compute current IU2 length (distance between IM2 and IM3)

      dist2 = sqrt((xem(im2) - xem(im3)) ** 2  &
                 + (yem(im2) - yem(im3)) ** 2  &
                 + (zem(im2) - zem(im3)) ** 2)

! Compute current IU3 length (distance between IM1 and IM4)

      dist3 = sqrt((xem(im1) - xem(im4)) ** 2  &
                 + (yem(im1) - yem(im4)) ** 2  &
                 + (zem(im1) - zem(im4)) ** 2)

! Compute current IU4 length (distance between IM2 and IM4)

      dist4 = sqrt((xem(im2) - xem(im4)) ** 2  &
                 + (yem(im2) - yem(im4)) ** 2  &
                 + (zem(im2) - zem(im4)) ** 2)

! Compute cosine of angles at IM3 and IM4

      cosphi3 = (dist1**2 + dist2**2 - dist**2) / (2. * dist1 * dist2) 
      cosphi4 = (dist3**2 + dist4**2 - dist**2) / (2. * dist3 * dist4) 

! Ratio of smaller cosine to limiting value of cos(72 deg)

      ratio = min(cosphi3,cosphi4) / .309
      
! Decrease dist0 if ratio < 1
! (NOTE: could just use .309 limit regardless of npoly)

      if (ratio < 1.) then
         dist0 = dist0 * ratio
!alt     dist0 = dist0 * (1. - ratio * frac)
      endif

! Fractional change to dist that would make it equal dist0

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
!$omp end do

! Loop over all M points

!$omp do private(expansion)
   do im = 2,nma

      if (ngr > 1 .and. movem(im) == 0) cycle

! For now, prevent either polar M point from moving

!      if (im == impent(1 )) cycle
!      if (im == impent(12)) cycle

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
!$omp end do

!$omp end parallel

! Section for plotting grid at intermediate stages of spring dynamics adjustment

   CYCLE ! Bypass grid plotting

   if (ngr > 1 .and. mod(iter,20) == 1) then

! Plot grid lines

      call o_reopnwk()
      call plotback()

      call oplot_set(1)
 
      do iu = 2,nua
         im1 = itab_ud(iu)%im(1)
         im2 = itab_ud(iu)%im(2)

         call oplot_transform(1,xem(im1),yem(im1),zem(im1),xp1,yp1)
         call oplot_transform(1,xem(im2),yem(im2),zem(im2),xp2,yp2)

         call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

         if (iskip == 1) cycle

         call o_frstpt (xq1,yq1)
         call o_vector (xq2,yq2)
      enddo

      call o_frame()
      call o_clswk()

   endif ! mod(iter,*)

enddo ! iter

return
end subroutine spring_dynamics
