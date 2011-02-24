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
subroutine hurricane_frances_init()

use mem_ijtabs
use mem_basic
use misc_coms
use mem_grid
use consts_coms

! Use stronger initial perturbation than in build_frances directory
! NOW, Use stronger initial perturbation than in build_frances9200b directory

implicit none

integer :: iw,i,j,k,ngr,iu,iv,iw1,iw2,kbc,irad,iter
real :: timescale,vertfac,latfac,dtemp,rads,pkhyd
real :: rad,radfac,vtan
real :: zef,ref,xef,yef,fcentlat,fcentlon
real :: unxf,unyf,unzf,unxrad,unyrad,unzrad,unxtan,unytan,unztan

real :: ulat,ulon,val95,deltmp,exner,deltheta,temp,rhovs,val

real :: rhovs_eye0 ,rhovs_eye15 ,rhovs_eye30
real :: rhovs_wall0,rhovs_wall15,rhovs_wall30
real :: rhovs_260k0,rhovs_260k15,rhovs_260k30,rhovs_260k60,rhovs_260k75  &
       ,rhovs_260k100,rhovs_260k110

real, external :: rhovsl,rhovsi

! lat/lon coords of eye center (correct observed location)

fcentlat = 22.3
fcentlon = -71.4

! Find "earth" coordinates of hurricane center

zef = erad * sin(fcentlat * pio180)
ref = erad * cos(fcentlat * pio180)  ! distance from earth axis
xef = ref  * cos(fcentlon * pio180)
yef = ref  * sin(fcentlon * pio180)

! Components of unit vector outward normal to earth surface at hurricane center

unxf = xef / erad
unyf = yef / erad
unzf = zef / erad

! Compute saturation vapor density for various observed dewpoint values

rhovs_eye0  = rhovsl(29.0)  ! eye at surface
rhovs_eye15 = rhovsl(24.0)  ! eye at 1500 m
rhovs_eye30 = rhovsl(21.0)  ! eye at 3000 m

rhovs_wall0  = rhovsl(29.0)  ! eyewall at surface
rhovs_wall15 = rhovsl(24.0)  ! eyewall at 1500 m
rhovs_wall30 = rhovsl(19.0)  ! eyewall at 3000 m

rhovs_260k0    = rhovsl(27.0)   ! 260 km radius at surface
rhovs_260k15   = rhovsl(19.0)   ! 260 km radius at 1500 m
rhovs_260k30   = rhovsl(12.0)   ! 260 km radius at 3000 m
rhovs_260k60   = rhovsl(-3.0)   ! 260 km radius at 6000 m
rhovs_260k75   = rhovsl(-14.0)  ! 260 km radius at 7500 m
rhovs_260k100  = rhovsl(-29.0)  ! 260 km radius at 10000 m
rhovs_260k110  = rhovsl(-42.0)  ! 260 km radius at 11000 m
   
print*, 'rhovs ',rhovs_eye0,rhovs_eye30,rhovs_260k60,rhovs_260k100

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)  ! jend(1) = hardwired for mrl 1
!----------------------------------------------------------------------
call qsub('W',iw)

! Distance of this IW point from eye center

   rad = sqrt((xew(iw)-xef)**2 + (yew(iw)-yef)**2 + (zew(iw)-zef)**2)
   
! Skip hurricane assimilation for all points outside specified radius

   if (rad >= 400.e3) cycle

! Define radial profile of temperature perturbation field based on 
! observed Hurricane Frances temperatures at 700 mb on morning of 2 September 2004.
! This is done by piecewise linear interpolation.  Outside a radius of 150 km,
! linearly trend perturbation to zero.  
! Form for each segment is: t = t1 + (t2-t1) * (r-r1) / r2-r1)

   if (rad < 30.e3) then
      deltmp = 16. + (9. - 16.) *  (rad - 0.e3)  / (30.e3 - 0.e3)
   elseif (rad < 200.e3) then
      deltmp = 9.  + (7. - 9.)  * (rad - 30.e3)  / (200.e3 - 30.e3)
   else
      deltmp = 7.  + (0. - 7.)  * (rad - 200.e3) / (400.e3 - 200.e3)
   endif   

   do k = lpw(iw),mza-1

! Define vertically-dependent modulating factor that varies pertubation with height,
! with maximum value at 700 mb (3000 m).  Use exponent to control vertical profile
! of factor. Skip current level if above height where weight factor is positive

      vertfac = 1. - (abs(zt(k) - 3000.) / 10000.) ** 2.0

      if (vertfac <= 0.) cycle
         
!      vertfac = max(0.,(15000. - zt(k)) / 15000.)

! Compute theta perturbation from temperature perturbation

      exner = (press(k,iw) / p00) ** rocp   ! Defined WITHOUT CP factor
      deltheta = deltmp / exner
      
      thil(k,iw)  = thil(k,iw)  + deltheta * vertfac
      theta(k,iw) = theta(k,iw) + deltheta * vertfac 

   enddo

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hydrostatically balance modified initial profile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Select model level as top boundary condition for hydrostatic integration:
! Profile remains unchanged at and above this level 

   kbc = 35

! Carry out iterative hydrostatic balance procedure

   do iter = 1,100
      do k = kbc-1,1,-1

! Set vapor specific humidity to 90% of saturation inside radius of 300 km 
! and below 13 km.  Between this radius and 400 km, linearly trend to 
! reanalysis value.

         if (zt(k) < 1500.) then
            if (rad <= 25.e3) then
               val = rhovs_eye0 + (rhovs_eye15 - rhovs_eye0)  &
                  * (zt(k) - 0.) / (1500. - 0.)
               sh_w(k,iw) = val / rho(k,iw)
               sh_v(k,iw) = sh_w(k,iw)
            elseif (rad <= 120.e3) then
               val = rhovs_wall0 + (rhovs_wall15 - rhovs_wall0)  &
                  * (zt(k) - 0.) / (1500. - 0.)
               sh_w(k,iw) = val / rho(k,iw)
               sh_v(k,iw) = sh_w(k,iw)
            else
               val = rhovs_260k0 + (rhovs_260k15 - rhovs_260k0)  &
                  * (zt(k) - 0.) / (1500. - 0.)
               sh_w(k,iw) = val / rho(k,iw)
               sh_v(k,iw) = sh_w(k,iw)
            endif
         elseif (zt(k) < 3000.) then
            if (rad <= 25.e3) then
               val = rhovs_eye15 + (rhovs_eye30 - rhovs_eye15)  &
                  * (zt(k) - 1500.) / (3000. - 1500.)
               sh_w(k,iw) = val / rho(k,iw)
               sh_v(k,iw) = sh_w(k,iw)
            elseif (rad <= 120.e3) then
               val = rhovs_wall15 + (rhovs_wall30 - rhovs_wall15)  &
                  * (zt(k) - 1500.) / (3000. - 1500.)
               sh_w(k,iw) = val / rho(k,iw)
               sh_v(k,iw) = sh_w(k,iw)
            else
               val = rhovs_260k15 + (rhovs_260k30 - rhovs_260k15)  &
                  * (zt(k) - 1500.) / (3000. - 1500.)
               sh_w(k,iw) = val / rho(k,iw)
               sh_v(k,iw) = sh_w(k,iw)
            endif
         elseif (zt(k) < 6000.) then
               val = rhovs_260k30 + (rhovs_260k60 - rhovs_260k30)  &
                  * (zt(k) - 3000.) / (6000. - 3000.)
               sh_w(k,iw) = val / rho(k,iw)
               sh_v(k,iw) = sh_w(k,iw)
         elseif (zt(k) < 7500.) then
               val = rhovs_260k60 + (rhovs_260k75 - rhovs_260k60)  &
                  * (zt(k) - 6000.) / (7500. - 6000.)
               sh_w(k,iw) = val / rho(k,iw)
               sh_v(k,iw) = sh_w(k,iw)
         elseif (zt(k) < 10000.) then
               val = rhovs_260k75 + (rhovs_260k100 - rhovs_260k75)  &
                  * (zt(k) - 75000.) / (10000. - 75000.)
               sh_w(k,iw) = val / rho(k,iw)
               sh_v(k,iw) = sh_w(k,iw)
         elseif (zt(k) < 11000.) then
               val = rhovs_260k100 + (rhovs_260k110 - rhovs_260k100)  &
                  * (zt(k) - 10000.) / (11000. - 10000.)
               sh_w(k,iw) = val / rho(k,iw)
               sh_v(k,iw) = sh_w(k,iw)
         endif
      

!if (rad < 10.e3) then
!print*, 'initfr ',k,iw,val,sh_w(k,iw)
!endif

!  Compute density - (ASSUME MICPHYS LEVEL = 1 FOR NOW)

         rho(k,iw) = press(k,iw) ** cvocp * p00k  &
            / (theta(k,iw) * (rdry * (1. - sh_v(k,iw)) + rvap * sh_v(k,iw)))

! Hydrostatically integrate downward using weighting to damp oscillations

         pkhyd = press(k+1,iw)  &
            + gravo2 * (rho(k+1,iw) * dzt(k+1) + rho(k,iw) * dzt(k))
         press(k,iw) = .05 * press(k,iw) + .95 * pkhyd

      enddo
   enddo

enddo
call rsub('W_frances_a',7)

! Enhance vortex winds since reanalysis does not resolve them

! A comparison of dropwindsonde wind data and reanalysis fields at 1200 UTC
! on 02 September 2004 indicates that hurricane winds are represented 
! in the reanalysis reasonably well outside a radius of 400 km.  Thus, add a 
! perturbation to winds inside this radius such that the strength of the  
! perturbation smoothly tends toward zero at the 400 km radius.

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_u(7)%jend(1); iu = jtab_u(7)%iu(j)  ! jend(1) = hardwired for mrl 1
   iw1 = itab_u(iu)%iw1; iw2 = itab_u(iu)%iw2
!----------------------------------------------------------------------
call qsub('U',iu)
 
! Distance of this point from eye center

   rad = sqrt((xeu(iu)-xef)**2 + (yeu(iu)-yef)**2 + (zeu(iu)-zef)**2)

! Skip current iu column if at or outside limiting radius

   if (rad >= 400.e3) cycle

! Define radial profile of tangential wind perturbation field based on 
! observed Hurricane Frances winds at 700 mb level on morning of 2 September 2004.
! This is done by piecewise linear interpolation.  Outside a radius of 200 km,
! Linearly trend perturbation to zero.  
! Form for each segment is: v = v1 + (v2-v1) * (r-r1) / r2-r1)

   if (rad < 30.e3) then
      vtan = 0.  + (70. - 0.) *  (rad - 0.e3)   / (30.e3 - 0.e3)
   elseif (rad < 60.e3) then
      vtan = 70. + (70. - 70.) * (rad - 30.e3)  / (60.e3 - 30.e3)
   elseif (rad < 90.e3) then
      vtan = 55. + (55. - 70.) * (rad - 60.e3)  / (90.e3 - 60.e3)
   elseif (rad < 150.e3) then
      vtan = 46. + (46. - 55.) * (rad - 90.e3)  / (150.e3 - 90.e3)
   elseif (rad < 200.e3) then
      vtan = 34. + (34. - 46.) * (rad - 150.e3) / (200.e3 - 150.e3)
   else
      vtan = 34. + (0. - 34.)  * (rad - 200.e3) / (400.e3 - 200.e3)
   endif   

! Unit normal vector components from hurricane center to current IU point

   unxrad = (xeu(iu) - xef) / rad
   unyrad = (yeu(iu) - yef) / rad
   unzrad = (zeu(iu) - zef) / rad

! Unit vector components in direction of tangential vortex wind

   unxtan = unyf * unzrad - unzf * unyrad
   unytan = unzf * unxrad - unxf * unzrad
   unztan = unxf * unyrad - unyf * unxrad

   do k = 1,mza
      
! Define vertically-dependent modulating factor that varies pertubation with height,
! with maximum value at 700 mb (3000 m).  Use exponent to control vertical profile
! of factor. Skip current level if above height where weight factor is positive

      vertfac = 1. - (abs(zt(k) - 3000.) / 10000.) ** 2.0

      if (vertfac <= 0.) cycle
         
! Add enhanced vortex to winds interpolated from NCEP reanalysis

      uc(k,iu) = uc(k,iu) + vtan * vertfac &
         * (unx(iu) * unxtan + uny(iu) * unytan + unz(iu) * unztan)

      umc(k,iu) = uc(k,iu) * .5 * (rho(k,iw1) + rho(k,iw2))
      ump(k,iu) = umc(k,iu)
   enddo

enddo
call rsub('U_frances_a',7)

return
end subroutine hurricane_frances_init


!**************************************************************************

subroutine bubble_test()

use mem_ijtabs
use mem_basic
use misc_coms
use mem_grid
use consts_coms

implicit none

integer :: iw,i,j,k,ngr,iu,iv,iw1,iw2
real :: timescale,vertfac,latfac,dtemp
real :: rad,radfac,vtan
real :: zef,ref,xef,yef,fcentlat,fcentlon
real :: unxf,unyf,unzf,unxrad,unyrad,unzrad,unxtan,unytan,unztan

iw = 90
k = 2

wc(k,90) = 10.
wmc(k,90) = wc(k,90) * .5 * (rho(k,90) + rho(k+1,90))

return
end subroutine bubble_test


