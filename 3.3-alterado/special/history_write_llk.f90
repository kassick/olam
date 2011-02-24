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
subroutine history_write_ll(nlon,nlat)

use mem_ijtabs,  only: itab_w, itab_u, rotation_angle
use mem_basic,   only: uc, wc, rho, press, theta

use mem_grid,    only: xew, yew, zew, &
                       topm, glatw, glonw, lpw, mza, mua, mwa, zm, zt, lpu,  &
                       xem, yem, zem, xeu, yeu, zeu,  &
                       unx, uny, unz, utx, uty, utz, dzt

use mem_addsc,   only: addsc

use misc_coms,   only: io6, time8, naddsc, timmax8, dtlong

use consts_coms, only: p00, rocp, piu180, grav, erad, pio180, cp

implicit none

integer, intent(in) :: nlon,nlat

real :: lat(nlat)
real :: lon(nlon)
real :: lev(mza-2)

integer :: k,ka
real :: pressw
real :: vx,vy,vz,raxis,u,v
real :: tempk
real :: tuu1,tuu2,tuu3,tuu4,vc

real :: alon, alat
real :: fldval, phis

integer :: iu,iu1,iu2,iu3,iu4,im1,im2,im3,iw1,iw2,iw
integer :: ilat,ilon,ilatlon
integer :: kll

real :: scr(mza,mwa), scr1(mwa), scr2(mwa)

! Interpolate some model fields to uniform lat-lon grid with spacing 
! comparable to global grid (at equator).

real :: phis_ll(nlon,nlat)

real :: u_ll (nlon,nlat,mza-2)
real :: v_ll (nlon,nlat,mza-2)
real :: w_ll (nlon,nlat,mza-2)
real :: t_ll (nlon,nlat,mza-2)
real :: p_ll (nlon,nlat,mza-2)
real :: r_ll (nlon,nlat,mza-2)

real :: z3_ll(nlon,nlat,mza-2)
real :: uzonal(mza,mua), umerid(mza,mua)

real :: alat1,alat2

do ilat = 1,nlat
   lat(ilat) = -90.0 + 180.0 * real(ilat-1) / real(nlat-1)
enddo

do ilon = 1,nlon
   lon(ilon) = 360.0 * real(ilon-1) / real(nlon)
enddo

do k=2,mza-1
   lev(k-1) = zt(k)
enddo

scr  = 0.0
scr1 = 0.0
scr2 = 0.0

!------------------------------------------------------------
! UZONAL and UMERID wind components at U point
!------------------------------------------------------------

do iu = 2,mua
   iu1 = itab_u(iu)%iu1
   iu2 = itab_u(iu)%iu2
   iu3 = itab_u(iu)%iu3
   iu4 = itab_u(iu)%iu4

   tuu1 = itab_u(iu)%tuu1
   tuu2 = itab_u(iu)%tuu2
   tuu3 = itab_u(iu)%tuu3
   tuu4 = itab_u(iu)%tuu4

   do k = 2,mza - 1
      kll = k - 1

      fldval = uc(k,iu1) * tuu1  &
             + uc(k,iu2) * tuu2  &
             + uc(k,iu3) * tuu3  &
             + uc(k,iu4) * tuu4

      vx = unx(iu) * uc(k,iu) + utx(iu) * fldval
      vy = uny(iu) * uc(k,iu) + uty(iu) * fldval
      vz = unz(iu) * uc(k,iu) + utz(iu) * fldval

      raxis = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)  ! dist from earth axis

      uzonal(k,iu) = (vy * xeu(iu) - vx * yeu(iu)) / raxis
      umerid(k,iu) = vz * raxis / erad  &
         - (vx * xeu(iu) + vy * yeu(iu)) * zeu(iu) / (raxis * erad) 

   enddo
enddo

print*, 'interpolating uzonal to lat-lon'

call interp_hnu_ll(nlon,nlat,mza,mza-2,uzonal,u_ll)

print*, 'interpolating umerid to lat-lon'

call interp_hnu_ll(nlon,nlat,mza,mza-2,umerid,v_ll)

!------------------------------------------------------------
! PHIS (topography height)
!------------------------------------------------------------

print*, 'interpolating uzonal and umerid to lat-lon'

do iw = 2,mwa
   im1 = itab_w(iw)%im1
   im2 = itab_w(iw)%im2
   im3 = itab_w(iw)%im3
   
   ka = lpw(iw)
   
   scr1(iw) = (topm(im1) + topm(im2) + topm(im3)) / 3.  ! PHIS
enddo

call interp_htw_ll(nlon,nlat,1,1,scr1,phis_ll)

!------------------------------------------------------------
! Z3 height of model levels
!------------------------------------------------------------

print*, 'filling z3_ll values'

do k = 2,mza-1
   kll = k - 1
   do ilat = 1,nlat
      do ilon = 1,nlon
         z3_ll(ilon,ilat,kll) = zt(k)
      enddo
   enddo
enddo

!------------------------------------------------------------
! W vertical velocity at each K level
!------------------------------------------------------------

print*, 'interpolating w to lat-lon'

do iw = 2,mwa
   do k = 2,mza
      scr(k,iw) = .5 * (wc(k-1,iw) + wc(k,iw))
   enddo
enddo

call interp_htw_ll(nlon,nlat,mza,mza-2,scr,w_ll)

!------------------------------------------------------------
! Pressure
!------------------------------------------------------------

print*, 'interpolating pressure to lat-lon'

do iw = 2,mwa
   do k = 1,mza
      scr(k,iw) = press(k,iw)
   enddo
enddo

call interp_htw_ll(nlon,nlat,mza,mza-2,scr,p_ll)

!------------------------------------------------------------
! Density
!------------------------------------------------------------

print*, 'interpolating density to lat-lon'

do iw = 2,mwa
   do k = 1,mza
      scr(k,iw) = rho(k,iw)
   enddo
enddo

call interp_htw_ll(nlon,nlat,mza,mza-2,scr,r_ll)

!------------------------------------------------------------
! Temperature
!------------------------------------------------------------

print*, 'interpolating temperature to lat-lon'

do iw = 2,mwa
   do k = 1,mza
      scr(k,iw) = theta(k,iw) * (press(k,iw) / p00) ** rocp
   enddo
enddo

call interp_htw_ll(nlon,nlat,mza,mza-2,scr,t_ll)

return
end subroutine history_write_ll

