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
Module leaf3_interface

Interface
   subroutine leaf3_land(iland, nlev_sfcwater, leaf_class, ntext_soil,        &
                      soil_water, soil_energy,                                &
                      sfcwater_mass, sfcwater_energy, sfcwater_depth,         &
                      rshort_s,  rshort_v, rshort_g, rshort,                  &
                      rlong_v, rlong_s, rlong_g,                              &
                      veg_height, veg_rough, veg_tai, veg_lai, veg_fracarea,  &
                      hcapveg, can_depth,                                     &
                      rhos, vels, prss, pcpg, qpcpg, dpcpg,                   &
                      sxfer_t, sxfer_r, ustar, snowfac, vf,                   &
                      surface_ssh, ground_shv, veg_water, veg_temp,           &
                      can_temp, can_shv, stom_resist, veg_ndvip, veg_ndvif,   &
                      veg_ndvic, veg_albedo, rough, lsl,                      &
		      timefac_ndvi, time8,                                    &
                      ed_patch)

     use ed_structure_defs, only: patch
     use leaf_coms, only: nzg,nzs
    
     integer, intent(in) :: lsl ! lowest soil layer
     integer, intent(in   ) :: iland             ! current land cell index
     integer, intent(inout) :: nlev_sfcwater     ! # active levels of surface water
     integer, intent(in   ) :: leaf_class        ! leaf ("vegetation") class
     integer, intent(in   ) :: ntext_soil  (nzg) ! soil textural class

     real, intent(inout) :: soil_water     (nzg) ! soil water content [vol_water/vol_tot]
     real, intent(inout) :: soil_energy    (nzg) ! soil energy [J/m^3]
     real, intent(inout) :: sfcwater_mass  (nzs) ! surface water mass [kg/m^2]
     real, intent(inout) :: sfcwater_energy(nzs) ! surface water energy [J/kg]
     real, intent(inout) :: sfcwater_depth (nzs) ! surface water depth [m]
     real, intent(in   ) :: rshort_s       (nzs) ! s/w net rad flux to sfc water [W/m^2]
     
     real, intent(in   ) :: rshort_v     ! s/w net rad flux to veg [W/m^2]
     real, intent(in   ) :: rshort_g     ! s/w net rad flux to soil [W/m^2]
     real, intent(in   ) :: rshort       ! s/w incident sfc rad flux [W/m^2]
     real, intent(in   ) :: rlong_v      ! l/w net rad flux to veg [W/m^2]
     real, intent(in   ) :: rlong_s      ! l/w net rad flux to sfc water [W/m^2]
     real, intent(in   ) :: rlong_g      ! l/w net rad flux to soil [W/m^2]
     real, intent(in   ) :: veg_height   ! veg height [m]
     real, intent(inout) :: veg_rough    ! veg roughness height [m]
     real, intent(inout) :: veg_tai      ! veg total area index
     real, intent(inout) :: veg_lai      ! veg leaf area index
     real, intent(inout) :: veg_fracarea ! veg fractional area
     real, intent(in   ) :: hcapveg      ! veg heat capacity [J/(m^2 K)]
     real, intent(in   ) :: can_depth    ! depth of canopy air [m]
     real, intent(in   ) :: timefac_ndvi ! frac of time from past to future NDVI obs
     real, intent(in   ) :: rhos         ! atmosphere air density [kg/m^3]
     real, intent(in   ) :: vels         ! surface wind speed [m/s]
     real, intent(in   ) :: prss         ! pressure [Pa]
     real, intent(inout) :: pcpg         ! new precip amount this leaf timestep [kg/m^2]
     real, intent(inout) :: qpcpg        ! new precip energy this leaf timestep [J/m^2]
     real, intent(inout) :: dpcpg        ! new precip depth this leaf timestep [m]
     real, intent(in   ) :: sxfer_t      ! can-to-atm heat xfer this step [kg_air K/m^2]
     real, intent(in   ) :: sxfer_r      ! can-to-atm vapor xfer this step [kg_vap/m^2]
     real, intent(in   ) :: ustar        ! friction velocity [m/s]
     real, intent(in   ) :: snowfac      ! fractional veg burial by snowcover
     real, intent(in   ) :: vf           ! fractional coverage of non-buried part of veg
     real, intent(inout) :: surface_ssh  ! surface saturation spec hum [kg_vap/kg_air]
     real, intent(inout) :: ground_shv   ! soil vapor spec hum [kg_vap/kg_air]
     real, intent(inout) :: veg_water    ! veg sfc water content [kg/m^2]
     real, intent(inout) :: veg_temp     ! veg temperature [K]
     real, intent(inout) :: can_temp     ! canopy air temperature [K]
     real, intent(inout) :: can_shv      ! canopy air vapor spec hum [kg/kg]
     real, intent(inout) :: stom_resist  ! veg stomatal resistance [s/m]
     real, intent(inout) :: veg_ndvip    ! past veg ndvi (obs time)
     real, intent(inout) :: veg_ndvif    ! past veg ndvi (obs time)
     real, intent(inout) :: veg_ndvic    ! current veg ndvi
     real, intent(inout) :: veg_albedo   ! veg albedo
     real, intent(inout) :: rough        ! net land cell roughness height [m]
     real(kind=8), intent(in) :: time8   ! model time
     
     type(patch), target, optional :: ed_patch

   end subroutine leaf3_land


   subroutine canopy(iland,                                               &
		  nlev_sfcwater,           ntext_soil,                    &
		  leaf_class,              ktrans,                        &
                  soil_water,              soil_fracliq,                  &
		  soil_tempk,              sfcwater_mass,                 &
		  sfcwater_tempk,          veg_height,                    &
		  veg_rough,               veg_tai,                       &
		  veg_lai,                 hcapveg,                       &
		  can_depth,               rhos,                          &
		  vels,                    prss,                          &
		  pcpg,                    qpcpg,                         &
                  rshort,                  rshort_v,                      &
		  rlong_v,                 sxfer_t,                       &
		  sxfer_r,                 ustar,                         &
                  snowfac,                 vf,                            &
		  surface_ssh,             ground_shv,                    &
                  veg_water,               veg_temp,                      &
		  can_temp,                can_shv,                       &
                  wshed, qwshed, transp, stom_resist,                     &
                  hxfergc, wxfergc,                                       &
		  hxfersc, wxfersc,                                       &
		  hxferca, hxfervc,                                       &
		  wxfervc, rdi,                                           &
		  rb, time8,                                              &
                  lsl, ed_transp, ed_patch                                )

     use leaf_coms,         only: nzg
     use ed_structure_defs, only: patch

     integer, intent(in)  :: iland           ! index of current land cell
     integer, intent(in)  :: nlev_sfcwater   ! # active levels of surface water
     integer, intent(in)  :: ntext_soil(nzg) ! soil textural class
     integer, intent(in)  :: leaf_class      ! leaf class (vegetation class)
     integer, intent(out) :: ktrans          ! k index of soil layer supplying transp

     integer, intent(in) :: lsl
     real, intent(in) :: soil_water(nzg)   ! soil water content [vol_water/vol_tot]
     real, intent(in) :: soil_fracliq(nzg) ! fraction of soil moisture in liquid phase
     real, intent(in) :: soil_tempk        ! soil temp [K]
     real, intent(in) :: sfcwater_mass     ! surface water mass [kg/m^2]
     real, intent(in) :: sfcwater_tempk    ! surface water temp [K]
     real, intent(in) :: veg_height  ! veg height [m]
     real, intent(in) :: veg_rough   ! veg roughess height [m]
     real, intent(in) :: veg_tai     ! veg total area index
     real, intent(in) :: veg_lai     ! veg leaf area index
     real, intent(in) :: hcapveg     ! veg heat capacity [J/(m^2 K)]
     real, intent(in) :: can_depth   ! depth of canopy air [m]
     real, intent(in) :: rhos        ! atmospheric air density [kg/m^3]
     real, intent(in) :: vels        ! surface wind speed [m/s]
     real, intent(in) :: prss        ! air pressure [Pa]
     real, intent(in) :: pcpg        ! new precip amount this leaf timestep [kg/m^2]
     real, intent(in) :: qpcpg       ! new precip energy this leaf timestep [J/m^2]
     real, intent(in) :: rshort      ! downward sfc s/w rad flux [W/m^2]
     real, intent(in) :: rshort_v    ! s/w rad flux absorbed by veg [W/m^2]
     real, intent(in) :: rlong_v     ! l/w rad flux absorbed by veg [W/m^2]
     real, intent(in) :: sxfer_t     ! surface heat xfer this step [kg_air K/m^2]
     real, intent(in) :: sxfer_r     ! surface vapor xfer this step [kg_vap/m^2]
     real, intent(in) :: ustar       ! friction velocity [m/s]
     real, intent(in) :: snowfac     ! fractional veg burial by snowcover
     real, intent(in) :: vf          ! fractional coverage of non-buried part of veg
     real, intent(in) :: surface_ssh ! surface sat spec hum [kg_vap/kg_air]
     real, intent(in) :: ground_shv  ! soil vapor spec hum [kg_vap/kg_air]
     
     real, intent(inout) :: veg_water   ! veg sfc water content [kg/m^2]
     real, intent(inout) :: veg_temp    ! veg temp [K]
     real, intent(inout) :: can_temp    ! canopy air temp [K]
     real, intent(inout) :: can_shv     ! canopy air vapor spec hum [kg_vap/kg_air]
     real, intent(inout) :: stom_resist ! veg stomatal resistance [s/m]
     
     real, intent(out) :: hxfergc ! (soil or sfc water)-to-can_air heat xfer this step [J/m^2]
     real, intent(out) :: wxfergc ! (soil or sfc water)-to-can_air vap xfer this step [kg_vap/m^2]
     real, intent(out) :: hxfersc ! (soil or sfc water)-to-can_air heat xfer this step [J/m^2]
     real, intent(out) :: wxfersc ! (soil or sfc water)-to-can_air vap xfer this step [kg_vap/m^2]
     real, intent(out) :: wshed   ! water shed from veg this LEAF timestep [kg/m^2]
     real, intent(out) :: qwshed  ! water energy shed from veg this LEAF timestep [J/m^2]
     real, intent(out) :: transp  ! transpiration xfer this LEAF timestep [kg_vap/m^2]
     real, intent(out) :: hxferca ! can_air-to-atm heat xfer this step [J/m^2]
     real, intent(out) :: hxfervc ! veg-to-can_air heat xfer this step [J/m^2]
     real, intent(out) :: wxfervc ! veg-to-can_air vapor xfer this step [kg_vap/m^2]
     real, intent(out) :: rdi     ! (soil or surface water)-to-can_air conductance [m/s]
     real, intent(out) :: rb      ! veg-to-can_air resistance [s/m]
     real(kind=8), intent(in) :: time8   ! model time [s]
     real :: ed_transp(nzg) ! Transpiration from each soil layer as computed by ED2 [kg/m2]
     type(patch), optional, target :: ed_patch

   end subroutine canopy


   subroutine soil(iland,                                 &
		leaf_class,      nlev_sfcwater,           &
		ntext_soil,      ktrans,                  &
                soil_tempk,      soil_fracliq,            &
		soil_rfactor,    hxfergc,                 &
		wxfergc,         rshort_g,                &
		rlong_g,         transp,                  &
                soil_energy,     soil_water,              &
		hxferg,          wxfer,                   &
		qwxfer,          psiplusz,                &
		hydraul_cond,    lsl,                     &
		ed_transp,       ed_patch                 )

     use leaf_coms,         only: nzg
     use ed_structure_defs, only: patch

     integer, intent(in) :: iland            ! current land cell number 
     integer, intent(in) :: leaf_class       ! leaf class
     integer, intent(in) :: nlev_sfcwater    ! # active levels of surface water
     integer, intent(in) :: ntext_soil (nzg) ! soil textural class
     
     integer, intent(in) :: ktrans           ! k index of soil layer supplying transpiration
     integer, intent(in) :: lsl
     real, intent(in) :: soil_tempk    (nzg) ! soil temperature (K)
     real, intent(in) :: soil_fracliq  (nzg) ! fraction of soil water that is liquid
     real, intent(in) :: soil_rfactor  (nzg) ! soil thermal resistance
     real, intent(in) :: hxfergc
     real, intent(in) :: wxfergc             ! water xfer from ground to canopy [kg/m^2]
     real, intent(in) :: rshort_g            ! s/w radiative flux abs by ground [W/m^2]
     real, intent(in) :: rlong_g             ! l/w radiative flux abs by ground [W/m^2]
     real, intent(in) :: transp              ! transpiration loss [kg/m^2]
     
     real, intent(inout) :: soil_energy(nzg) ! [J/m^3]
     real, intent(inout) :: soil_water (nzg) ! soil water content [vol_water/vol_tot]
     real, intent(inout) :: hxferg   (nzg+1) ! soil internal heat xfer (J/m^2]
     
     real, intent(out) :: wxfer       (nzg+1) ! soil water xfer [m]
     real, intent(out) :: qwxfer      (nzg+1) ! soil energy xfer from water xfer [J/m^2] 
     real, intent(out) :: psiplusz    (nzg)   ! soil water potential (including grav) [m]
     real, intent(out) :: hydraul_cond(nzg)   ! soil hydraulic conductivity [m/s]
     real :: ed_transp(nzg) ! Transpiration from each soil layer as computed by ED2 [kg/m2]
     type(patch), target, optional :: ed_patch
     
   end subroutine soil

End Interface
End Module leaf3_interface
