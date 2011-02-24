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
subroutine interp_hnu_ll(nlon,nlat,nlevin,nlevout,field,field_ll)

use mem_grid,   only: mma, mza, mwa, mua, zm, zt, lpw, lpu,  &
                      xem, yem, zem, xeu, yeu, zeu, xew, yew, zew, glatw, glonw
use mem_ijtabs, only: itab_m, itab_w
use misc_coms,  only: io6
use consts_coms, only: erad,piu180,pio180

implicit none

integer, intent(in) :: nlon,nlat,nlevin,nlevout
real, intent(in) :: field(nlevin,mua)
real, intent(out) :: field_ll(nlon,nlat,nlevout)

integer :: k,im,iw,ntpn,lpwmax,itpn,itpn_next,iw_next,ilat,ilon,kll
integer :: k1,k2,koff,lpumax,iu,iu_next,iu1,iu2,iu3

real :: x(3),y(3),z(3)
real :: xw(8),yw(8),xu(8),yu(8)

real :: a(nlevin),b(nlevin),c(nlevin)
real :: field_avg(nlevin)

real :: raxis,glatm,glonm,radmax,rad
real :: deglat,abslat,aminlat,amaxlat,aminlon,amaxlon,deglon
real :: v0x,v0y,v1x,v1y,v2x,v2y,dot00,dot01,dot02,dot11,dot12,denomi,u,v
real :: spacelat,spacelon,qlat,qlon,qx,qy

if (nlevin == 1) then
   k1 = 1
   k2 = 1
   koff = 0
else
   k2 = mza-1
   koff = 1
endif

x(1) = 0.
y(1) = 0.

!-----------------------------------------------------------
! First loop is over M points for interpolating U points
!-----------------------------------------------------------

do im = 2,mma
   ntpn = itab_m(im)%ntpn

! Find lat/lon of current M point

   raxis = sqrt(xem(im) ** 2 + yem(im) ** 2)  ! dist from earth axis

   glatm = atan2(zem(im),raxis)   * piu180
   glonm = atan2(yem(im),xem(im)) * piu180

! Initialize maximum radius, maximum lpu, and field average

   radmax = 0.
   lpumax = 2
   field_avg(1:nlevin) = 0.

! Loop over all U points that surround current M point

   do itpn = 1,ntpn

! Current W point index   

      iu = itab_m(im)%iu(itpn)
      
! Skip this IM point if iu < 2

      if (iu < 2) go to 9

! Transform current U point to PS coordinates tangent at M point

      call e_ps(xeu(iu),yeu(iu),zeu(iu),glatm,glonm,xu(itpn),yu(itpn))

      rad = sqrt(xu(itpn)**2 + yu(itpn)**2)
      if (radmax < rad)     radmax = rad
      if (lpumax < lpu(iu)) lpumax = lpu(iu)

      if (nlevin > 1) then
         k1 = lpumax
      endif

      do k = k1,k2
         field_avg(k) = field_avg(k) + field(k,iu) / real(ntpn)
      enddo

   enddo
   
! Find maximum range of latitude and longitude inside the polygon of U points

   deglat = piu180 * radmax / erad
   
   aminlat = max(-90.,glatm - deglat)
   amaxlat = min( 90.,glatm + deglat)
   
   abslat = max(abs(aminlat),abs(amaxlat))
   
   if (abslat > 89.9) then
   
      aminlon =   -0.01
      amaxlon =  360.01
      
   else

      deglon = min(180.,deglat / cos(abslat * pio180))
      aminlon = glonm - deglon
      amaxlon = glonm + deglon
      
! transform aminlon and amaxlon to range 0,360

      if (aminlon < 0.) aminlon = aminlon + 360.
      if (amaxlon < 0.) amaxlon = amaxlon + 360.      
      
   endif

! Loop over all U points that surround current M point and fill field values

   do itpn = 1,ntpn
      itpn_next = itpn + 1
      if (itpn == ntpn) itpn_next = 1

      iu      = itab_m(im)%iu(itpn)
      iu_next = itab_m(im)%iu(itpn_next)

      x(2) = xu(itpn)
      y(2) = yu(itpn)
      
      x(3) = xu(itpn_next)
      y(3) = yu(itpn_next)

! Loop over vertical levels

      do k = k1,k2
         kll = k - koff

         z(1) = field_avg(k)
         z(2) = field(k,iu)
         z(3) = field(k,iu_next)

! Evaluate interpolation coefficients for current trio of points

         call matrix_3x3(1.  , x(1), y(1),  &
                         1.  , x(2), y(2),  &
                         1.  , x(3), y(3),  &
                         z(1), z(2), z(3),  &
                         a(k), b(k), c(k)   )
      enddo             

! Set up some triangle-check coefficients

      v0x = x(2) - x(1)
      v0y = y(2) - y(1)
      
      v1x = x(3) - x(1)
      v1y = y(3) - y(1) 
      
      dot00 = v0x * v0x + v0y * v0y
      dot01 = v0x * v1x + v0y * v1y
      dot11 = v1x * v1x + v1y * v1y

      denomi = 1. / (dot00 * dot11 - dot01 * dot01)
            
! Loop over all possible lat-lon points in range

      spacelat = 180. / real(nlat-1)
      spacelon = 360. / real(nlon)
      
      do ilat = 1,nlat
         qlat = -90. + spacelat * real(ilat-1)

         if (qlat < aminlat .or. qlat > amaxlat) cycle

         do ilon = 1,nlon  ! loop over longitude points
                  
            qlon = 0. + spacelon * real(ilon-1)

            if (amaxlon > aminlon .and. &
               (qlon < aminlon .or. qlon > amaxlon)) cycle
               
            if (amaxlon < aminlon .and. &
               qlon < aminlon .and. qlon > amaxlon) cycle
               
! Transform current lat-lon point to PS space

            call ll_xy(qlat,qlon,glatm,glonm,qx,qy)         
         
! Set up remaining triangle_check coefficients

            v2x = qx - x(1)
            v2y = qy - y(1)
            
            dot02 = v0x * v2x + v0y * v2y
            dot12 = v1x * v2x + v1y * v2y

            u = (dot11 * dot02 - dot01 * dot12) * denomi
            v = (dot00 * dot12 - dot01 * dot02) * denomi

! Check if current qx,qy point is inside or very near current triangle

            if (u > -.01 .and. v > -.01 .and. u + v < 1.01) then

! Point is inside or very near triangle; loop over vertical levels

               do k = k1,k2
                  kll = k - koff

! Interpolate to current field point

                  field_ll(ilon,ilat,kll) = a(k) + b(k) * qx + c(k) * qy
        
               enddo  ! k
               
            endif  ! q point inside triangle

         enddo ! ilon

      enddo  ! ilat

   enddo   ! itpn

9  continue

enddo   ! im

!-----------------------------------------------------------
! Second loop is over W points for interpolating U points
!-----------------------------------------------------------

do iw = 2,mwa
   iu1 = itab_w(iw)%iu1
   iu2 = itab_w(iw)%iu2
   iu3 = itab_w(iw)%iu3

! Transform U points to PS coordinates tangent at W point

   call e_ps(xeu(iu1),yeu(iu1),zeu(iu1),glatw(iw),glonw(iw),x(1),y(1))
   call e_ps(xeu(iu2),yeu(iu2),zeu(iu2),glatw(iw),glonw(iw),x(2),y(2))
   call e_ps(xeu(iu3),yeu(iu3),zeu(iu3),glatw(iw),glonw(iw),x(3),y(3))

   radmax = max(sqrt(x(1)**2 + y(1)**2),  &
                sqrt(x(2)**2 + y(2)**2),  &
                sqrt(x(3)**2 + y(3)**2)   )

   lpumax = max(lpu(iu1),lpu(iu2),lpu(iu3))

   if (nlevin > 1) then
      k1 = lpumax
   endif

! Find maximum range of latitude and longitude inside the polygon of W points

   deglat = piu180 * radmax / erad
   
   aminlat = max(-90.,glatw(iw) - deglat)
   amaxlat = min( 90.,glatw(iw) + deglat)
   
   abslat = max(abs(aminlat),abs(amaxlat))
   
   if (abslat > 89.9) then
   
      aminlon =   -0.01
      amaxlon =  360.01
      
   else

      deglon = min(180.,deglat / cos(abslat * pio180))
      aminlon = glonw(iw) - deglon
      amaxlon = glonw(iw) + deglon
      
! transform aminlon and amaxlon to range 0,360

      if (aminlon < 0.) aminlon = aminlon + 360.
      if (amaxlon < 0.) amaxlon = amaxlon + 360.      
      
   endif

! Loop over vertical levels

   do k = k1,k2
      kll = k - koff

      z(1) = field(k,iu1)
      z(2) = field(k,iu2)
      z(3) = field(k,iu3)

! Evaluate interpolation coefficients for current trio of points

      call matrix_3x3(1.  , x(1), y(1),  &
                      1.  , x(2), y(2),  &
                      1.  , x(3), y(3),  &
                      z(1), z(2), z(3),  &
                      a(k), b(k), c(k)   )
   enddo             

! Set up some triangle-check coefficients

   v0x = x(2) - x(1)
   v0y = y(2) - y(1)
      
   v1x = x(3) - x(1)
   v1y = y(3) - y(1) 
      
   dot00 = v0x * v0x + v0y * v0y
   dot01 = v0x * v1x + v0y * v1y
   dot11 = v1x * v1x + v1y * v1y

   denomi = 1. / (dot00 * dot11 - dot01 * dot01)
            
! Loop over all possible lat-lon points in range

   spacelat = 180. / real(nlat-1)
   spacelon = 360. / real(nlon)
      
   do ilat = 1,nlat
      qlat = -90. + spacelat * real(ilat-1)

      if (qlat < aminlat .or. qlat > amaxlat) cycle

      do ilon = 1,nlon  ! loop over longitude points
                  
         qlon = 0. + spacelon * real(ilon-1)

         if (amaxlon > aminlon .and. &
            (qlon < aminlon .or. qlon > amaxlon)) cycle
               
         if (amaxlon < aminlon .and. &
            qlon < aminlon .and. qlon > amaxlon) cycle
               
! Transform current lat-lon point to PS space

         call ll_xy(qlat,qlon,glatw(iw),glonw(iw),qx,qy)         
         
! Set up remaining triangle_check coefficients

         v2x = qx - x(1)
         v2y = qy - y(1)
            
         dot02 = v0x * v2x + v0y * v2y
         dot12 = v1x * v2x + v1y * v2y

         u = (dot11 * dot02 - dot01 * dot12) * denomi
         v = (dot00 * dot12 - dot01 * dot02) * denomi

! Check if current qx,qy point is inside or very near current triangle

         if (u > -.01 .and. v > -.01 .and. u + v < 1.01) then

! Point is inside or very near triangle; loop over vertical levels

            do k = k1,k2
               kll = k - koff

! Interpolate to current field point

               field_ll(ilon,ilat,kll) = a(k) + b(k) * qx + c(k) * qy
        
            enddo  ! k
               
         endif  ! q point inside triangle

      enddo ! ilon

   enddo  ! ilat

enddo   ! iw

return
end subroutine interp_hnu_ll

!================================================================================

subroutine interp_htw_ll(nlon,nlat,nlevin,nlevout,field,field_ll)

use mem_grid,   only: mma, mza, mwa, zm, zt, lpw, xem, yem, zem, xew, yew, zew
use mem_ijtabs, only: itab_m
use misc_coms,  only: io6
use consts_coms, only: erad,piu180,pio180

implicit none

integer, intent(in) :: nlon,nlat,nlevin,nlevout
real, intent(in) :: field(nlevin,mwa)
real, intent(out) :: field_ll(nlon,nlat,nlevout)

integer :: k,im,iw,ntpn,lpwmax,itpn,itpn_next,iw_next,ilat,ilon,kll
integer :: k1,k2,koff

real :: x(3),y(3),z(3)
real :: xw(8),yw(8)

real :: a(nlevin),b(nlevin),c(nlevin)
real :: field_avg(nlevin)

real :: raxis,glatm,glonm,radmax,rad
real :: deglat,abslat,aminlat,amaxlat,aminlon,amaxlon,deglon
real :: v0x,v0y,v1x,v1y,v2x,v2y,dot00,dot01,dot02,dot11,dot12,denomi,u,v
real :: spacelat,spacelon,qlat,qlon,qx,qy

if (nlevin == 1) then
   k1 = 1
   k2 = 1
   koff = 0
else
   k2 = mza-1
   koff = 1
endif

x(1) = 0.
y(1) = 0.

! Loop over M points for interpolating W points

do im = 2,mma

   ntpn = itab_m(im)%ntpn

! Find lat/lon of current M point

   raxis = sqrt(xem(im) ** 2 + yem(im) ** 2)  ! dist from earth axis

   glatm = atan2(zem(im),raxis)   * piu180
   glonm = atan2(yem(im),xem(im)) * piu180

! Initialize maximum radius, maximum lpw, and field average

   radmax = 0.
   lpwmax = 2
   field_avg(1:nlevin) = 0.

! Loop over all W points that surround current M point

   do itpn = 1,ntpn

! Current W point index   

      iw = itab_m(im)%iw(itpn)
      
! Skip this IM point if iw < 2

      if (iw < 2) go to 9

! Transform current W point to PS coordinates tangent at M point

      call e_ps(xew(iw),yew(iw),zew(iw),glatm,glonm,xw(itpn),yw(itpn))

      rad = sqrt(xw(itpn)**2 + yw(itpn)**2)
      if (radmax < rad)     radmax = rad
      if (lpwmax < lpw(iw)) lpwmax = lpw(iw)

      if (nlevin > 1) then
         k1 = lpwmax
      endif

      do k = k1,k2
         field_avg(k) = field_avg(k) + field(k,iw) / real(ntpn)
      enddo

   enddo
   
! Find maximum range of latitude and longitude inside the polygon of W points

   deglat = piu180 * radmax / erad
   
   aminlat = max(-90.,glatm - deglat)
   amaxlat = min( 90.,glatm + deglat)
   
   abslat = max(abs(aminlat),abs(amaxlat))
   
   if (abslat > 89.9) then
   
      aminlon =   -0.01
      amaxlon =  360.01
      
   else

      deglon = min(180.,deglat / cos(abslat * pio180))
      aminlon = glonm - deglon
      amaxlon = glonm + deglon
      
! transform aminlon and amaxlon to range 0,360

      if (aminlon < 0.) aminlon = aminlon + 360.
      if (amaxlon < 0.) amaxlon = amaxlon + 360.      
      
   endif

! Loop over all W points that surround current M point and fill field values

   do itpn = 1,ntpn
      itpn_next = itpn + 1
      if (itpn == ntpn) itpn_next = 1

      iw      = itab_m(im)%iw(itpn)
      iw_next = itab_m(im)%iw(itpn_next)

      x(2) = xw(itpn)
      y(2) = yw(itpn)
      
      x(3) = xw(itpn_next)
      y(3) = yw(itpn_next)

! Loop over vertical levels

      do k = k1,k2
         kll = k - koff

         z(1) = field_avg(k)
         z(2) = field(k,iw)
         z(3) = field(k,iw_next)

! Evaluate interpolation coefficients for current trio of points

         call matrix_3x3(1.  , x(1), y(1),  &
                         1.  , x(2), y(2),  &
                         1.  , x(3), y(3),  &
                         z(1), z(2), z(3),  &
                         a(k), b(k), c(k)   )
      enddo             

! Set up some triangle-check coefficients

      v0x = x(2) - x(1)
      v0y = y(2) - y(1)
      
      v1x = x(3) - x(1)
      v1y = y(3) - y(1) 
      
      dot00 = v0x * v0x + v0y * v0y
      dot01 = v0x * v1x + v0y * v1y
      dot11 = v1x * v1x + v1y * v1y

      denomi = 1. / (dot00 * dot11 - dot01 * dot01)
            
! Loop over all possible lat-lon points in range

      spacelat = 180. / real(nlat-1)
      spacelon = 360. / real(nlon)
      
      do ilat = 1,nlat
         qlat = -90. + spacelat * real(ilat-1)

         if (qlat < aminlat .or. qlat > amaxlat) cycle

         do ilon = 1,nlon  ! loop over longitude points
                  
            qlon = 0. + spacelon * real(ilon-1)

            if (amaxlon > aminlon .and. &
               (qlon < aminlon .or. qlon > amaxlon)) cycle
               
            if (amaxlon < aminlon .and. &
               qlon < aminlon .and. qlon > amaxlon) cycle
               
! Transform current lat-lon point to PS space

            call ll_xy(qlat,qlon,glatm,glonm,qx,qy)         
         
! Set up remaining triangle_check coefficients

            v2x = qx - x(1)
            v2y = qy - y(1)
            
            dot02 = v0x * v2x + v0y * v2y
            dot12 = v1x * v2x + v1y * v2y

            u = (dot11 * dot02 - dot01 * dot12) * denomi
            v = (dot00 * dot12 - dot01 * dot02) * denomi

! Check if current qx,qy point is inside or very near current triangle

            if (u > -.01 .and. v > -.01 .and. u + v < 1.01) then

! Point is inside or very near triangle; loop over vertical levels

               do k = k1,k2
                  kll = k - koff

! Interpolate to current field point

                  field_ll(ilon,ilat,kll) = a(k) + b(k) * qx + c(k) * qy
        
               enddo  ! k
               
            endif  ! q point inside triangle

         enddo ! ilon

      enddo  ! ilat

   enddo   ! itpn

9  continue

enddo   ! im

return
end subroutine interp_htw_ll


