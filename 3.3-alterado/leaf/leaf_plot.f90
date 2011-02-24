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
subroutine leaf_plot(iland, nlev_sfcwater, ntext_soil, leaf_class, ktrans, &
                     soil_water, soil_rfactor, soil_energy, soil_tempk,    &
                     soil_fracliq, hxferg, wxfer, qwxfer,                  &
                     psiplusz, hydraul_cond,                               &
                     sfcwater_mass, sfcwater_energy, sfcwater_depth,       &
                     sfcwater_tempk, sfcwater_fracliq,                     &
                     rshort_s, rshort_v, rshort_g, rshort,                 &
                     rlong_v, rlong_s, rlong_g,                            &
                     veg_height, veg_rough, veg_tai, veg_lai,              &
                     hcapveg, can_depth,                                   &
                     rhos, vels, prss, pcpg, qpcpg, dpcpg,                 &
                     sxfer_t, sxfer_r, ustar, snowfac, vf,                 &
                     surface_ssh, ground_shv,                              &
                     veg_water, veg_temp, can_temp, can_shv,               &
                     wshed, qwshed, transp, stom_resist,                   &
                     hxfergc, wxfergc, hxfersc, wxfersc, veg_fracarea,     &
                     hxferca, hxfervc, wxfervc, rdi, rb, time8             )

use leaf_coms,   only: nzg, nzs, slz, dslz
use consts_coms, only: cp
use oplot_coms,  only: op
use misc_coms,   only: io6, runtype

implicit none

integer, intent(in) :: iland              ! current land cell index
integer, intent(in) :: nlev_sfcwater      ! # active levels of surface water
integer, intent(in) :: ntext_soil   (nzg) ! soil textural class
integer, intent(in) :: leaf_class         ! leaf ("vegetation") class
integer, intent(in) :: ktrans             ! k index of soil layer supplying transp

real, intent(in) :: soil_water      (nzg) ! soil water content [vol_water/vol_tot]
real, intent(in) :: soil_rfactor    (nzg) ! soil thermal resistivity [K m^2/W]
real, intent(in) :: soil_energy     (nzg) ! soil energy [J/m^3]
real, intent(in) :: soil_tempk      (nzg) ! soil temperature [K]
real, intent(in) :: soil_fracliq    (nzg) ! fraction of soil moisture in liquid phase
real, intent(in) :: hxferg        (nzg+1) ! heat xfer between soil layers [J/m^2]
real, intent(in) :: wxfer         (nzg+1) ! soil water xfer [m]
real, intent(in) :: qwxfer        (nzg+1) ! soil energy xfer from water xfer [J/m^2]
real, intent(in) :: psiplusz        (nzg) ! soil water potential (including grav) [m]
real, intent(in) :: hydraul_cond    (nzg) ! soil hydraulic conductivity [m/s]
real, intent(in) :: sfcwater_mass   (nzs) ! surface water mass [kg/m^2]
real, intent(in) :: sfcwater_energy (nzs) ! surface water energy [J/kg]
real, intent(in) :: sfcwater_depth  (nzs) ! surface water depth [m]
real, intent(in) :: sfcwater_tempk  (nzs) ! surface water temperature [K]
real, intent(in) :: sfcwater_fracliq(nzs) ! fraction of sfc water in liquid phase(nzs)
real, intent(in) :: rshort_s        (nzs) ! l/w net rad flux to sfc water [W/m^2]
real, intent(in) :: rshort_v     ! s/w net rad flux to veg [W/m^2]
real, intent(in) :: rshort_g     ! s/w net rad flux to soil [W/m^2]
real, intent(in) :: rshort       ! s/w incident sfc rad flux [W/m^2]
real, intent(in) :: rlong_v      ! l/w net rad flux to veg [W/m^2]
real, intent(in) :: rlong_s      ! l/w net rad flux to sfc water [W/m^2]
real, intent(in) :: rlong_g      ! l/w net rad flux to soil [W/m^2]
real, intent(in) :: veg_height   ! veg height [m]
real, intent(in) :: veg_rough    ! veg roughness height [m]
real, intent(in) :: veg_tai      ! veg total area index
real, intent(in) :: veg_lai      ! veg leaf area index
real, intent(in) :: hcapveg      ! veg heat capacity [J/(m^2 K)]
real, intent(in) :: can_depth    ! canopy depth for heat & vap capacity [m]
real, intent(in) :: rhos         ! atmosphere air density [kg/m^3]
real, intent(in) :: vels         ! surface wind speed [m/s]
real, intent(in) :: prss         ! pressure [Pa]
real, intent(in) :: pcpg         ! new precip amount this leaf timestep [kg/m^2]
real, intent(in) :: qpcpg        ! new precip energy this leaf timestep [J/m^2]
real, intent(in) :: dpcpg        ! new precip depth this leaf timestep [m]
real, intent(in) :: sxfer_t      ! can-to-atm heat xfer this step [kg_air K/m^2]
real, intent(in) :: sxfer_r      ! can-to-atm vapor xfer this step [kg_vap/m^2]
real, intent(in) :: ustar        ! friction velocity [m/s]
real, intent(in) :: snowfac      ! fractional veg burial by snowcover
real, intent(in) :: vf           ! fractional coverage of non-buried part of veg
real, intent(in) :: surface_ssh  ! surface saturation spec hum [kg_vap/kg_air]
real, intent(in) :: ground_shv   ! soil vapor spec hum [kg_vap/kg_air]
real, intent(in) :: veg_water    ! veg sfc water content [kg/m^2]
real, intent(in) :: veg_temp     ! veg temperature [K]
real, intent(in) :: can_temp     ! canopy air temperature [K]
real, intent(in) :: can_shv      ! canopy air vapor spec hum [kg/kg]
real, intent(in) :: wshed        ! water shed from veg this timestep [kg/m^2]
real, intent(in) :: qwshed       ! water energy shed from veg this timestep [J/m^2]
real, intent(in) :: transp       ! transpiration xfer this LEAF timestep [kg/m^2]
real, intent(in) :: stom_resist  ! veg stomatal resistance [s/m]
real, intent(in) :: hxfergc      ! heat xfer from soil to can_air this step [J/m^2]
real, intent(in) :: wxfergc      ! vap xfer from soil to can_air this step [kg_vap/m^2]
real, intent(in) :: hxfersc      ! heat xfer from sfcwater to can_air this step [J/m^2]
real, intent(in) :: wxfersc      ! vap xfer from sfcwater to can_air this step [kg_vap/m^2]
real, intent(in) :: veg_fracarea ! veg fractional area ! veg fractional area
real, intent(in) :: hxferca      ! can_air-to-atm heat xfer this step [J/m^2]
real, intent(in) :: hxfervc      ! veg-to-can_air heat xfer this step [J/m^2]
real, intent(in) :: wxfervc      ! veg-to-can_air vapor xfer this step [kg_vap/m^2]
real, intent(in) :: rdi          ! (soil or surface water)-to-can_air conductance [m/s]
real, intent(in) :: rb           ! veg-to-can_air resistance [s/m]
real(kind=8), intent(in) :: time8 ! model time [s]

! Local variables

integer :: k

real :: psiz
real :: ybot
real :: dy_ss

real :: x1,x2,xv1,xv2,y1,y2,yc2,ys2,yv1,yv2,yg2

real :: xw1
real :: xw2
real :: xw3
real :: xw4
real :: xw5
real :: xw6
real :: xw7
real :: xw8
real :: xw9
real :: xw10
real :: xw11

real :: yw1
real :: yw2
real :: yw3
real :: yw4
real :: yw5
real :: yw6
real :: yw7
real :: yw8

real :: yk1
real :: yk2
real :: yk3
real :: yk4
real :: yk5

real :: ywk

real :: xtpn(4)
real :: ytpn(4)

real :: water
real :: energy

character(len=20) :: number

! Reopen the current graphics output workstation if it is closed
call o_reopnwk()

! Initialize plot coordinates

psiz = .007

x1 = .05
x2 = .95

xv1 = .30
xv2 = .70

y1 = .12
y2 = .95
yc2 = .90
ys2 = .70
yv1 = .75
yv2 = .85

! key coordinates

yk1 = .02
yk5 = .08

yk3 = .5 * (yk1 + yk5)
yk2 = .5 * (yk1 + yk3)
yk4 = .5 * (yk3 + yk5)

dy_ss = (ys2 - y1) / real(nzg+nzs)
yg2 = y1 + real(nzg) * dy_ss

! Plot white background

call plotback()

! Scale plotter coordinates for unit square

call o_set(0.,1.,0.,1.,0.,1.,0.,1.,1)

! Fill soil background color

xtpn(1) = x1
xtpn(2) = x2
xtpn(3) = x2
xtpn(4) = x1

ytpn(1) = y1
ytpn(2) = y1
ytpn(3) = ys2
ytpn(4) = ys2

call fillpolyg (4,xtpn,ytpn,121)

! Fill soil key background color

xtpn(1) = x1
xtpn(2) = x2
xtpn(3) = x2
xtpn(4) = x1

ytpn(1) = yk1
ytpn(2) = yk1
ytpn(3) = yk3
ytpn(4) = yk3

call fillpolyg (4,xtpn,ytpn,121)

! Fill canopy air background color

xtpn(1) = x1
xtpn(2) = x2
xtpn(3) = x2
xtpn(4) = x1

ytpn(1) = yg2
ytpn(2) = yg2
ytpn(3) = y2
ytpn(4) = y2

call fillpolyg (4,xtpn,ytpn,15)

! Fill atm background color

xtpn(1) = x1
xtpn(2) = x2
xtpn(3) = x2
xtpn(4) = x1

ytpn(1) = yc2
ytpn(2) = yc2
ytpn(3) = y2
ytpn(4) = y2

call fillpolyg (4,xtpn,ytpn,100)

! Fill snowcover background color

if (nlev_sfcwater > 0) then

   xtpn(1) = x1
   xtpn(2) = x2
   xtpn(3) = x2
   xtpn(4) = x1

   ytpn(1) = yg2
   ytpn(2) = yg2
   ytpn(3) = yg2 + real(nlev_sfcwater) * dy_ss
   ytpn(4) = yg2 + real(nlev_sfcwater) * dy_ss

   call fillpolyg (4,xtpn,ytpn,141)

endif

! Fill snowcover key background color

if (nlev_sfcwater > 0) then

   xtpn(1) = x1
   xtpn(2) = x2
   xtpn(3) = x2
   xtpn(4) = x1

   ytpn(1) = yk3
   ytpn(2) = yk3
   ytpn(3) = yk5
   ytpn(4) = yk5

   call fillpolyg (4,xtpn,ytpn,141)

endif

! Fill veg background color

if (leaf_class > 3) then
   xtpn(1) = xv1
   xtpn(2) = xv2
   xtpn(3) = xv2
   xtpn(4) = xv1

   ytpn(1) = yv1
   ytpn(2) = yv1
   ytpn(3) = yv2
   ytpn(4) = yv2

   call fillpolyg (4,xtpn,ytpn,108)
endif

! Set color for cell boundaries

call o_sflush()
call o_gsplci(143)
call o_gstxci(143)
call o_sflush()

! Vertical sides

call o_frstpt(x1,y1)
call o_vector(x1,y2)

call o_frstpt(x2,y1)
call o_vector(x2,y2)

! Vertical key sides

call o_frstpt(x1,yk1)
call o_vector(x1,yk5)

call o_frstpt(x2,yk1)
call o_vector(x2,yk5)

! Soil and snowcover layers

do k = 1,nzg+nlev_sfcwater+1
   ybot = y1 + real(k-1) * dy_ss

   call o_frstpt(x1,ybot)
   call o_vector(x2,ybot)
enddo

! Soil and snowcover key layers

call o_frstpt(x1,yk1)
call o_vector(x2,yk1)

call o_frstpt(x1,yk3)
call o_vector(x2,yk3)

call o_frstpt(x1,yk5)
call o_vector(x2,yk5)

! Top of canopy
   
call o_frstpt(x1,yc2)
call o_vector(x2,yc2)

! Top of atm layers
   
call o_frstpt(x1,y2)
call o_vector(x2,y2)

! vegetation border

if (leaf_class > 3) then
   call o_frstpt(xv1,yv1)
   call o_vector(xv2,yv1)
   call o_vector(xv2,yv2)
   call o_vector(xv1,yv2)
   call o_vector(xv1,yv1)
endif

! Sum water [kg/m^2] and energy [J/m^2] in land cell

water = 0.
energy = 0.

! Soil contribution

do k = 1,nzg
   water = water + 1000. * soil_water(k) * dslz(k)
   energy = energy + soil_energy(k) * dslz(k)
enddo

! Surface_water contribution

do k = 1,nlev_sfcwater
   water = water + sfcwater_mass(k)
   energy = energy + sfcwater_energy(k) * sfcwater_mass(k)
enddo

! Canopy air contribution

water = water + can_depth * rhos * can_shv
energy = energy + can_depth * rhos * cp * (can_temp - 273.15)

! Vegetation contribution

water = water + veg_water
energy = energy + hcapveg * (veg_temp - 273.15)

! Flush plot buffer and set font number to 4 (font 2 is similar but wider spacing)

call o_sflush()
call o_pcseti ('FN',4)

! Set color for printed numbers

call o_gsplci(10)
call o_gstxci(10)
call o_sflush()

! Print values in soil layers

xw1  = x1 +  1. * (x2 - x1) / 22.
xw2  = x1 +  3. * (x2 - x1) / 22.
xw3  = x1 +  5. * (x2 - x1) / 22.
xw4  = x1 +  7. * (x2 - x1) / 22.
xw5  = x1 +  9. * (x2 - x1) / 22.
xw6  = x1 + 11. * (x2 - x1) / 22.
xw7  = x1 + 13. * (x2 - x1) / 22.
xw8  = x1 + 15. * (x2 - x1) / 22.
xw9  = x1 + 17. * (x2 - x1) / 22.
xw10 = x1 + 19. * (x2 - x1) / 22.
xw11 = x1 + 21. * (x2 - x1) / 22.

do k = 1,nzg
   ybot = y1 + real(k-1) * dy_ss
   ywk = ybot + .50 * dy_ss

   write (number,'(i5)') k
   call o_plchhq(.03,ywk,trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f10.3)') dslz(k)
   call o_plchhq(xw1,ywk,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw1,yk2,'dslz',psiz,0.,0.)

   write (number,'(i5)') ntext_soil(k)
   call o_plchhq(xw2,ywk,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw2,yk2,'ntext',psiz,0.,0.)

   write (number,'(f10.6)') soil_water(k)
   call o_plchhq(xw3,ywk,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw3,yk2,'water',psiz,0.,0.)

   write (number,'(f12.3)') soil_energy(k) * 1.e-6  ! convert from J/m^3 to J/cm^3
   call o_plchhq(xw4,ywk,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw4,yk2,'energy',psiz,0.,0.)

   write (number,'(f12.3)') soil_tempk(k)
   call o_plchhq(xw5,ywk,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw5,yk2,'tempk',psiz,0.,0.)

   write (number,'(f12.3)') soil_fracliq(k)
   call o_plchhq(xw6,ywk,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw6,yk2,'fracliq',psiz,0.,0.)

   write (number,'(f12.3)') psiplusz(k)
   call o_plchhq(xw7,ywk,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw7,yk2,'psi+z',psiz,0.,0.)

   write (number,'(e12.3)') soil_rfactor(k)
   call o_plchhq(xw8,ywk,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw8,yk2,'rfactor',psiz,0.,0.)

   if (k == nzg) then
      write (number,'(f12.1)') rshort_g
      call o_plchhq(xw9,ywk,trim(adjustl(number)),psiz,0.,0.)
      call o_plchhq(xw9,yk2,'rshort_g',psiz,0.,0.)

      write (number,'(f12.1)') rlong_g
      call o_plchhq(xw10,ywk,trim(adjustl(number)),psiz,0.,0.)
      call o_plchhq(xw10,yk2,'rlong_g',psiz,0.,0.)

   endif

   write (number,'(f10.6)') soil_water(k) * 1000. * dslz(k)
   call o_plchhq(xw11,ywk,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw11,yk2,'water_a',psiz,0.,0.)

enddo

! Print values between soil layers

do k = 1,nzg + 1
   ybot = y1 + real(k-1) * dy_ss

   write (number,'(f12.3)') slz(k)
   call o_plchhq(xw1,ybot,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw1,yk1,'slz',psiz,0.,0.)
 
   write (number,'(e12.3)') hxferg(k)
   call o_plchhq(xw2,ybot,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw2,yk1,'hxferg',psiz,0.,0.)
 
   write (number,'(e12.3)') wxfer(k)
   call o_plchhq(xw3,ybot,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw3,yk1,'wxfer',psiz,0.,0.)
 
   write (number,'(e12.3)') qwxfer(k)
   call o_plchhq(xw4,ybot,trim(adjustl(number)),psiz,0.,0.)
   if (k == 1) call o_plchhq(xw4,yk1,'qwxfer',psiz,0.,0.)

   if (k > 1 .and. k <= nzg) then
      write (number,'(e12.3)') hydraul_cond(k)
      call o_plchhq(xw5,ybot,trim(adjustl(number)),psiz,0.,0.)
      call o_plchhq(xw5,yk1,'hydcond',psiz,0.,0.)
   endif
enddo

! Print values in snowcover layers

if (nlev_sfcwater > 0) then

   do k = 1,nlev_sfcwater
      ybot = yg2 + real(k-1) * dy_ss
      ywk = ybot + .50 * dy_ss

      write (number,'(i5)') k
      call o_plchhq(.03,ywk,trim(adjustl(number)),psiz,0.,0.)

      write (number,'(f12.3)') sfcwater_depth(k)
      call o_plchhq(xw1,ywk,trim(adjustl(number)),psiz,0.,0.)
      if (k == 1) call o_plchhq(xw1,yk4,'depth',psiz,0.,0.)

      write (number,'(f12.3)') sfcwater_mass(k)
      call o_plchhq(xw3,ywk,trim(adjustl(number)),psiz,0.,0.)
      if (k == 1) call o_plchhq(xw3,yk4,'mass',psiz,0.,0.)

      write (number,'(f12.3)') sfcwater_energy(k)
      call o_plchhq(xw4,ywk,trim(adjustl(number)),psiz,0.,0.)
      if (k == 1) call o_plchhq(xw4,yk4,'energy',psiz,0.,0.)

      write (number,'(f12.3)') sfcwater_tempk(k)
      call o_plchhq(xw5,ywk,trim(adjustl(number)),psiz,0.,0.)
      if (k == 1) call o_plchhq(xw5,yk4,'tempk',psiz,0.,0.)

      write (number,'(f12.3)') sfcwater_fracliq(k)
      call o_plchhq(xw6,ywk,trim(adjustl(number)),psiz,0.,0.)
      if (k == 1) call o_plchhq(xw6,yk4,'fracliq',psiz,0.,0.)
 
      write (number,'(f12.1)') rshort_s(k)
      call o_plchhq(xw9,ywk,trim(adjustl(number)),psiz,0.,0.)
      if (k == 1) call o_plchhq(xw9,yk4,'rshort_s',psiz,0.,0.)

      if (k == nlev_sfcwater) then
         write (number,'(f12.1)') rlong_s
         call o_plchhq(xw10,ywk,trim(adjustl(number)),psiz,0.,0.)
         call o_plchhq(xw10,yk4,'rlong_s',psiz,0.,0.)
      endif

   enddo

endif

! Print values at bottom of canopy air

yw1 = y1 + (real(nzg+nlev_sfcwater) + .30) * dy_ss
yw2 = y1 + (real(nzg+nlev_sfcwater) + .60) * dy_ss

write (number,'(f12.2)') hxfersc
call o_plchhq(xw1,yw1,'hxfersc = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.6)') wxfersc
call o_plchhq(xw3,yw1,'wxfersc = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.2)') hxfergc
call o_plchhq(xw5,yw1,'hxfergc = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.6)') wxfergc
call o_plchhq(xw7,yw1,'wxfergc = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.6)') rdi
call o_plchhq(xw9,yw1,'rdi = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.3)') surface_ssh * 1000.
call o_plchhq(xw11,yw1,'surface_ssh = '//trim(adjustl(number)),psiz,0.,0.)

if (nlev_sfcwater == 0) then
   write (number,'(f12.3)') ground_shv * 1000.
   call o_plchhq(xw11,yw2,'ground_shv = '//trim(adjustl(number)),psiz,0.,0.)
endif
       
! Print sums at canopy top

yw1 = yc2 - .012

write (number,'(f12.1)') energy
call o_plchhq(xw1,yw1,'energy = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.5)') water
call o_plchhq(xw3,yw1,'water = '//trim(adjustl(number)),psiz,0.,0.)


! Print values in canopy air

write (number,'(i5)') iland
call o_plchhq(.18,yc2 - .034,'iland = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(i5)') leaf_class
call o_plchhq(.18,yc2 - .046,'leaf_class = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.2)') can_depth
call o_plchhq(.18,yc2 - .058,'can_depth = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.2)') can_temp
call o_plchhq(.18,yc2 - .070,'can_temp = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.6)') can_shv
call o_plchhq(.18,yc2 - .082,'can_shv = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.6)') can_shv * can_depth * rhos
call o_plchhq(.18,yc2 - .094,'can_water = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.2)') snowfac
call o_plchhq(.18,yc2 - .106,'snowfac = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.2)') vf
call o_plchhq(.18,yc2 - .118,'vf = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.2)') real(time8)
call o_plchhq(.18,yc2 - .130,'time8 = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.3)') real(time8)/86400.
call o_plchhq(.18,yc2 - .142,'days = '//trim(adjustl(number)),psiz,0.,0.)

! Print values in atmosphere

yw1 = yc2 + .012
yw2 = yc2 + .024
yw3 = yc2 + .036

write (number,'(f12.3)') rhos
call o_plchhq(xw1,yw1,'rhos = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.3)') vels
call o_plchhq(xw1,yw2,'vels = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.3)') prss
call o_plchhq(xw1,yw3,'prss = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.3)') rshort
call o_plchhq(xw3,yw1,'rshort = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.6)') pcpg
call o_plchhq(xw3,yw2,'pcpg = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.3)') qpcpg
call o_plchhq(xw3,yw3,'qpcpg= '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.6)') dpcpg
call o_plchhq(xw5,yw1,'dpcpg = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.3)') sxfer_t
call o_plchhq(xw5,yw2,'sxfer_t = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.6)') sxfer_r
call o_plchhq(xw5,yw3,'sxfer_r = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.3)') hxferca
call o_plchhq(xw7,yw1,'hxferca = '//trim(adjustl(number)),psiz,0.,0.)

write (number,'(f12.3)') ustar
call o_plchhq(xw7,yw2,'ustar = '//trim(adjustl(number)),psiz,0.,0.)

! Check if vegetation exists in this cell

if (leaf_class > 3) then

   yw1 = yv2 - .012
   yw2 = yv2 - .024
   yw3 = yv2 - .036
   yw4 = yv2 - .048
   yw5 = yv2 - .060
   yw6 = yv2 - .072
   yw7 = yv2 - .084
   yw8 = yv2 - .096
   
   xw1 = xv1 + .10
   xw2 = xv1 + .27
   xw3 = xv2 + .07

! Print values in vegetation

   write (number,'(f12.1)') veg_height
   call o_plchhq(xw1,yw1,'height = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.2)') veg_tai
   call o_plchhq(xw1,yw2,'tai = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.2)') veg_lai
   call o_plchhq(xw1,yw3,'lai = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.2)') veg_temp
   call o_plchhq(xw1,yw4,'temp = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.6)') veg_water
   call o_plchhq(xw1,yw5,'water = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.1)') hcapveg
   call o_plchhq(xw2,yw1,'hcap = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.3)') veg_rough
   call o_plchhq(xw2,yw2,'rough = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.3)') veg_fracarea
   call o_plchhq(xw2,yw3,'fracarea = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.1)') rshort_v
   call o_plchhq(xw2,yw4,'rshort_v = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.1)') rlong_v
   call o_plchhq(xw2,yw5,'rlong_v = '//trim(adjustl(number)),psiz,0.,0.)

! Print fluxes to/from vegetation

   write (number,'(f12.1)') stom_resist
   call o_plchhq(xw3,yw1,'rc = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.1)') rb
   call o_plchhq(xw3,yw2,'rb = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.6)') transp
   call o_plchhq(xw3,yw3,'transp = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(i5)') ktrans
   call o_plchhq(xw3,yw4,'ktrans = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.1)') hxfervc
   call o_plchhq(xw3,yw5,'hxfervc = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.6)') wxfervc
   call o_plchhq(xw3,yw6,'wxfervc = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.6)') wshed
   call o_plchhq(xw3,yw7,'wshed = '//trim(adjustl(number)),psiz,0.,0.)

   write (number,'(f12.1)') qwshed
   call o_plchhq(xw3,yw8,'qwshed = '//trim(adjustl(number)),psiz,0.,0.)

endif

call o_frame()

! Close the current workstation if output
! is to a NCAR graphics meta file. This allows viewing the complete
! meta file (including the last frame) during a run and in case the
! simulation crashes.

if ((trim(runtype) .ne. 'PLOTONLY') .and. (op%plttype .eq. 0)) then
   call o_clswk()
endif

return
end subroutine leaf_plot

