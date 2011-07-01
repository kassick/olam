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
subroutine leaf3()

use leaf_coms, only: nzg, nzs, mwl, iupdndvi, s1900_ndvi, indvifile, nndvifiles

use mem_leaf,  only: land, itab_wl
use misc_coms, only: io6, time8, s1900_sim, iparallel

use ed_structure_defs, only: patch,site
use ed_options,        only: ied_offline
use leaf3_interface,   only: leaf3_land
use mem_para,          only: myrank

implicit none

! Local variables

integer :: iwl
real :: timefac_ndvi  ! fraction of elapsed time from past to future NDVI obs

type(patch), pointer :: ed_patch
type(site),  pointer :: ed_site

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAM_LEAF3_IN,rst_buf)
#endif


! Time interpolation factor for updating NDVI

timefac_ndvi = 0.

if (iupdndvi == 1 .and. nndvifiles > 1) then
   timefac_ndvi = (s1900_sim             - s1900_ndvi(indvifile-1))  &
                / (s1900_ndvi(indvifile) - s1900_ndvi(indvifile-1))
endif

! Loop over ALL LAND CELLS

ed_site => land%first_site

do iwl = 2,mwl

! Skip IWL cell if running in parallel and primary rank of IWL /= MYRANK

   if (iparallel == 1 .and. itab_wl(iwl)%irank /= myrank) cycle

   if (land%ed_flag(iwl) == 0 .and. ied_offline == 0) then 
! Update LAND fields

   call leaf3_land(iwl               ,                                    &
      land%nlev_sfcwater        (iwl), land%leaf_class           (iwl), &
      land%ntext_soil     (1:nzg,iwl), land%soil_water     (1:nzg,iwl), &
      land%soil_energy    (1:nzg,iwl), land%sfcwater_mass  (1:nzs,iwl), &
      land%sfcwater_energy(1:nzs,iwl), land%sfcwater_depth (1:nzs,iwl), &
      land%rshort_s       (1:nzs,iwl), land%rshort_v             (iwl), &
      land%rshort_g             (iwl), land%rshort               (iwl), &
      land%rlong_v              (iwl), land%rlong_s              (iwl), &
      land%rlong_g              (iwl), land%veg_height           (iwl), &
      land%veg_rough            (iwl), land%veg_tai              (iwl), &
      land%veg_lai              (iwl), land%veg_fracarea         (iwl), &
      land%hcapveg              (iwl), land%can_depth            (iwl), &
      land%rhos                 (iwl), land%vels                 (iwl), &
      land%prss                 (iwl), land%pcpg                 (iwl), &
      land%qpcpg                (iwl), land%dpcpg                (iwl), &
      land%sxfer_t              (iwl), land%sxfer_r              (iwl), &
      land%ustar                (iwl), land%snowfac              (iwl), &
      land%vf                   (iwl), land%surface_ssh          (iwl), &
      land%ground_shv           (iwl), land%veg_water            (iwl), &
      land%veg_temp             (iwl), land%can_temp             (iwl), &
      land%can_shv              (iwl), land%stom_resist          (iwl), &
      land%veg_ndvip            (iwl), land%veg_ndvif            (iwl), &
      land%veg_ndvic            (iwl), land%veg_albedo           (iwl), &
      land%rough                (iwl), land%lsl                  (iwl), &
      timefac_ndvi                     , time8                              )

   elseif(land%ed_flag(iwl) == 1)then

! in this case, ED is being run.  Therefore, we want to loop over patches 
! and send the patch variables.

      ed_site%omean_precip = ed_site%omean_precip + land%pcpg(ed_site%iland)
      ed_site%omean_qprecip = ed_site%omean_qprecip + land%qpcpg(ed_site%iland)

      ed_patch => ed_site%oldest_patch
      
      do while(associated(ed_patch))

         ed_site%omean_netrad = ed_site%omean_netrad +   &
              ed_patch%area * (  &
              (1.0 - ed_patch%albedo_beam) *   &
              (land%rshort(ed_site%iland) -   &
              land%rshort_diffuse(ed_site%iland)) +  &
              (1.0 - ed_patch%albedo_diffuse) *   &
              land%rshort_diffuse(ed_site%iland) ) +   &
              land%rlong(ed_site%iland) *  &
              (1.0 - ed_patch%rlong_albedo) -   &
              ed_patch%rlongup

         call leaf3_land(iwl                                              , &
              ed_patch%nlev_sfcwater       , land%leaf_class         (iwl), &
              ed_patch%ntext_soil   (1:nzg), ed_patch%soil_water     (1:nzg), &
              ed_patch%soil_energy  (1:nzg), ed_patch%sfcwater_mass  (1:nzs), &
              ed_patch%sfcwater_energy(1:nzs), ed_patch%sfcwater_depth   (1:nzs), &
              ed_patch%rshort_s     (1:nzs), land%rshort_v           (iwl), &
              ed_patch%rshort_g            , land%rshort             (iwl), &
              land%rlong_v          (iwl), ed_patch%rlong_s               , &
              ed_patch%rlong_g             , ed_patch%veg_height            , &
              ed_patch%veg_rough           , ed_patch%lai                   , &
              land%veg_lai          (iwl), land%veg_fracarea       (iwl), &
              land%hcapveg          (iwl), ed_patch%can_depth             , &
              land%rhos             (iwl), land%vels               (iwl), &
              land%prss             (iwl), land%pcpg               (iwl), &
              land%qpcpg            (iwl), land%dpcpg              (iwl), &
              ed_patch%sxfer_t             , ed_patch%sxfer_r               , &
              ed_patch%ustar               , ed_patch%snowfac               , &
              land%vf               (iwl), ed_patch%surface_ssh           , &
              ed_patch%ground_shv          , land%veg_water          (iwl), &
              land%veg_temp         (iwl), ed_patch%can_temp              , &
              ed_patch%can_shv             , land%stom_resist        (iwl), &
              land%veg_ndvip        (iwl), land%veg_ndvif          (iwl), &
              land%veg_ndvic        (iwl), land%veg_albedo         (iwl), &
              ed_patch%rough               , land%lsl                (iwl), &
              timefac_ndvi                   , &
              time8                        , ed_patch                         )

         ed_patch => ed_patch%younger
         
      enddo
      
   endif

! Zero out LAND%SXFER_T(iwl) and LAND%SXFER_R(iwl) now that they have 
! been applied to the canopy

   land%sxfer_t(iwl) = 0.
   land%sxfer_r(iwl) = 0.

   if (land%ed_flag(iwl) == 1) then
      ed_patch => ed_site%oldest_patch
      do while(associated(ed_patch))
         ed_patch%sxfer_t = 0.0
         ed_patch%sxfer_r = 0.0
         ed_patch => ed_patch%younger
      enddo
      ed_site => ed_site%next_site
   endif

enddo

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_LEAF3_OUT,rst_buf)
#endif

return
end subroutine leaf3

!===============================================================================

subroutine leaf3_land(iwl, nlev_sfcwater, leaf_class, ntext_soil,           &
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
                      timefac_ndvi, time8,      &
                      ed_patch)
                              
use leaf_coms,         only: nzg, nzs, soil_rough, snow_rough, slcpd, dt_leaf
use misc_coms,         only: io6
use ed_structure_defs, only: patch
use leaf3_interface,   only: canopy, soil

implicit none

integer, intent(in   ) :: iwl             ! current land cell index
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
real, intent(in   ) :: can_depth    ! canopy depth for heat and vap capacity [m]
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
integer, intent(in) :: lsl          ! Lowest soil layer
real(kind=8), intent(in) :: time8   ! model time [s]

type(patch), target, optional :: ed_patch

! Local arrays

real :: soil_tempk      (nzg) ! soil temperature [K]
real :: soil_fracliq    (nzg) ! fraction of soil moisture in liquid phase
real :: soil_rfactor    (nzg) ! soil thermal resistivity [K m^2/W]
real :: sfcwater_tempk  (nzs) ! surface water temperature [K]
real :: sfcwater_fracliq(nzs) ! fraction of sfc water in liquid phase
real :: hxferg        (nzg+1) ! heat xfer between soil layers [J/m^2]

real :: wxfer       (nzg+1) ! soil water xfer [m]
real :: qwxfer      (nzg+1) ! soil energy xfer from water xfer [J/m^2] 
real :: psiplusz    (nzg)   ! soil water potential (including grav) [m]
real :: hydraul_cond(nzg)   ! soil hydraulic conductivity [m/s]

! Local variables

integer :: k     ! vertical index over soil layers
integer :: nlsw1 ! maximum of (1,nlev_sfcwater)

integer :: ktrans ! vertical index of soil layer supplying transpiration

real :: transp  ! transpiration xfer this LEAF timestep [kg/m^2]
real, dimension(nzg) :: ed_transp ! transpired water from each soil level; ED2 cells only [kg/m^2]
real :: hxfergc ! heat xfer from ground (soil) to can_air this step [J/m^2]
real :: wxfergc ! vapor xfer from ground (soil) to can_air this step [kg_vap/m^2]
real :: hxfersc ! heat xfer from sfcwater to can_air this step [J/m^2]
real :: wxfersc ! vapor xfer from sfcwater to can_air this step [kg_vap/m^2]
real :: wshed   ! water shed from veg this timestep [kg/m^2]
real :: qwshed  ! water energy shed from veg this timestep [J/m^2]
real :: hxferca ! can_air-to-atm heat xfer this step [J/m^2]
real :: hxfervc ! veg-to-can_air heat xfer this step [J/m^2]
real :: wxfervc ! veg-to-can_air vapor xfer this step [kg_vap/m^2]
real :: rdi     ! (soil or surface water)-to-can_air conductance [m/s]
real :: rb      ! veg-to-can_air resistance [s/m]

real, external :: rhovsil ! function to compute sat spec hum over liquid or ice

!=================================================================
! TEMPORARY RUNOFF VARIABLES
real :: runoff
real :: qrunoff
!=================================================================

if(.not.present(ed_patch))then

! Diagnose soil temperature and liquid fraction

do k = 1,nzg
   call qwtk(soil_energy(k),soil_water(k)*1.e3,  &
      slcpd(ntext_soil(k)),soil_tempk(k),soil_fracliq(k))
enddo

! Diagnose surface water temperature and liquid fraction

do k = 1,nlev_sfcwater
   call qtk(sfcwater_energy(k),sfcwater_tempk(k),sfcwater_fracliq(k))
enddo

! Update vegetation TAI, LAI, fractional area, albedo, and roughness

call vegndvi(iwl,                      &
             leaf_class,   timefac_ndvi, &
             veg_height,   veg_ndvip,    &
             veg_ndvif,    veg_ndvic,    &
             veg_tai,      veg_lai,      &
             veg_fracarea, veg_albedo,   &
             veg_rough                   ) 

else
   soil_tempk(lsl:nzg) = ed_patch%soil_tempk(lsl:nzg)
   soil_fracliq(lsl:nzg) = ed_patch%soil_fracliq(lsl:nzg)
   sfcwater_tempk(lsl:nlev_sfcwater) =   &
        ed_patch%sfcwater_tempk(lsl:nlev_sfcwater)
   sfcwater_fracliq(lsl:nlev_sfcwater) =   &
        ed_patch%sfcwater_fracliq(lsl:nlev_sfcwater)
   veg_fracarea = 1.0
endif

! Compute roughness length based on vegetation and snow.

rough = max(soil_rough,veg_rough) * (1. - snowfac) + snow_rough

! Evaluate turbulent exchanges of heat and moisture between vegetation and canopy air
! and also between soil or snow surface and canopy air.  Evaluate transfer of
! precipitation moisture and heat to vegetation and shed from vegetation to surface.
! Update vegetation and canopy air temperatures resulting from these 
! plus radiative fluxes.

nlsw1 = max(nlev_sfcwater,1)

call canopy(iwl,                                        &
            nlev_sfcwater,         ntext_soil,            &
            leaf_class,            ktrans,                &
            soil_water,            soil_fracliq,          &
            soil_tempk    (nzg),   sfcwater_mass (nlsw1), &
            sfcwater_tempk(nlsw1), veg_height,            &
            veg_rough,             veg_tai,               &
            veg_lai,               hcapveg,               &
            can_depth,             rhos,                  &
            vels,                  prss,                  &
            pcpg,                  qpcpg,                 &
            rshort,                rshort_v,              &
            rlong_v,               sxfer_t,               &
            sxfer_r,               ustar,                 &
            snowfac,               vf,                    &
            surface_ssh,           ground_shv,            &
            veg_water,             veg_temp,              &
            can_temp,              can_shv,               &
            wshed,                 qwshed,                &
            transp,                stom_resist,           &
            hxfergc,               wxfergc,               &
            hxfersc,               wxfersc,               &
            hxferca,               hxfervc,               &
            wxfervc,               rdi,                   &
            rb,                    time8,                 &
            lsl,                   ed_transp,             &
            ed_patch               )

! CALL SFCWATER:
!  1. Compute soil and sfcwater heat conductivities
!  2. Compute surface heat flux for top layer of soil or sfcwater
!  3. Evaluate internal and bottom fluxes for sfcwater
!  4. Update sfcwater layer energies due to heat flux and solar radiation
!  5. Evaluate melting and percolation of liquid through sfcwater layers

call sfcwater(iwl,                              &
              nlev_sfcwater,    ntext_soil,       &
              soil_rfactor,     soil_water,       &
              soil_energy,      sfcwater_mass,    &
              sfcwater_energy,  sfcwater_depth,   &
              soil_tempk,       soil_fracliq ,    &
              sfcwater_tempk,   sfcwater_fracliq, &
              rshort_s,         hxfersc,          &
              wxfersc,          rlong_g,          &
              rlong_s,          pcpg,             &
              qpcpg,            dpcpg,            &
              wshed,            qwshed,           &
              veg_fracarea,     lsl,time8             )

call soil(iwl,                       &
          leaf_class,   nlev_sfcwater, &
          ntext_soil,   ktrans,        &
          soil_tempk,   soil_fracliq,  &
          soil_rfactor, hxfergc,       &
          wxfergc,      rshort_g,      &
          rlong_g,      transp,        &
          soil_energy,  soil_water,    &
          hxferg,       wxfer,         &
          qwxfer,       psiplusz,      &
          hydraul_cond, lsl,           &
          ed_transp,    ed_patch       )

! Compute surface and ground vap mxrat for next timestep; put into surface_ssh 
! and ground_ssv.

nlsw1 = max(nlev_sfcwater,1)

call grndvap(iwl,                                     &
             nlev_sfcwater,          ntext_soil  (nzg), &
             soil_water     (nzg),   soil_energy (nzg), &
             sfcwater_energy(nlsw1), rhos,              &
             can_shv,                ground_shv,        &
             surface_ssh                                )
             
         
!if (iwl == 102) then
!   write(io6,*) 'sfcw0a ',ground_shv,surface_ssh,soil_water(nzg),ntext_soil(nzg)
!   write(io6,*) 'sfcw1a ',nlev_sfcwater,soil_energy(nzg),  &
!      sfcwater_energy(nlsw1),rhos,can_shv
!endif

if (present(ed_patch)) then

   do k = lsl,nzg
      call qwtk(soil_energy(k),soil_water(k)*1.e3,  &
           slcpd(ntext_soil(k)),soil_tempk(k),soil_fracliq(k))
   enddo
   
   ! Diagnose surface water temperature and liquid fraction
   
   do k = 1,nlev_sfcwater
      call qtk(sfcwater_energy(k),sfcwater_tempk(k),sfcwater_fracliq(k))
   enddo
   
   do k = lsl,nzg
      ed_patch%soil_tempk(k) = soil_tempk(k)
      ed_patch%soil_fracliq(k) = soil_fracliq(k)
   enddo
   
   do k = 1,nlev_sfcwater
      ed_patch%sfcwater_tempk(k) = sfcwater_tempk(k)
      ed_patch%sfcwater_fracliq(k) = sfcwater_fracliq(k)
   enddo

endif

!-----------------------------------------------------------------------------
! TEMPORARY UNTIL FULL LEAF-HYDRO MODEL IS COMPLETED WITH STREAM/RIVER RUNOFF:
! Simple representation of runoff

call remove_runoff(iwl,            nlev_sfcwater,    &
                   sfcwater_fracliq, sfcwater_mass,    &
                   sfcwater_tempk,   sfcwater_energy,  &
                   sfcwater_depth,   runoff,           &
                   qrunoff                             )
!-----------------------------------------------------------------------------

if(present(ed_patch))then
   ed_patch%omean_runoff = ed_patch%omean_runoff + runoff / dt_leaf
   ed_patch%omean_qrunoff = ed_patch%omean_qrunoff + qrunoff / dt_leaf
endif

! Call land patch plot routine for selected iwl values.
! Use '(iwl == 0)' to turn off this plot call.

if (iwl == 0)            &
   call leaf_plot(iwl,                             &
                  nlev_sfcwater,    ntext_soil,      &
                  leaf_class,       ktrans,          &
                  soil_water,       soil_rfactor,    &
                  soil_energy,      soil_tempk,      &
                  soil_fracliq,     hxferg,          &
                  wxfer,            qwxfer,          &
                  psiplusz,         hydraul_cond,    &
                  sfcwater_mass,    sfcwater_energy, &
                  sfcwater_depth,   sfcwater_tempk,  &
                  sfcwater_fracliq, rshort_s,        &
                  rshort_v,         rshort_g,        &
                  rshort,           rlong_v,         &
                  rlong_s,          rlong_g,         &
                  veg_height,       veg_rough,       &
                  veg_tai,          veg_lai,         &
                  hcapveg,          can_depth,       &
                  rhos,             vels,            &
                  prss,             pcpg,            &
                  qpcpg,            dpcpg,           &
                  sxfer_t,          sxfer_r,         &
                  ustar,            snowfac,         &
                  vf,               surface_ssh,     &
                  ground_shv,       veg_water,       &
                  veg_temp,         can_temp,        &
                  can_shv,          wshed,           &
                  qwshed,           transp,          &
                  stom_resist,      hxfergc,         &
                  wxfergc,          hxfersc,         &
                  wxfersc,          veg_fracarea,    &
                  hxferca,          hxfervc,         &
                  wxfervc,          rdi,             &
                  rb,               time8            )

! At this point, any precipitation has been added to vegetation and surface
! water, so zero out precipitation variables prior to new contribution by
! surface flux subroutines

if(present(ed_patch))then
   if(associated(ed_patch%younger))return
endif

pcpg  = 0.
qpcpg = 0.
dpcpg = 0.

return
end subroutine leaf3_land

!===============================================================================

subroutine canopy(iwl, nlev_sfcwater, ntext_soil, leaf_class, ktrans,   &
                  soil_water, soil_fracliq, soil_tempk,                   &
                  sfcwater_mass, sfcwater_tempk,                          &
                  veg_height, veg_rough, veg_tai, veg_lai,                &
                  hcapveg, can_depth,                                     &
                  rhos, vels, prss, pcpg, qpcpg,                          &
                  rshort, rshort_v, rlong_v, sxfer_t, sxfer_r, ustar,     &
                  snowfac, vf, surface_ssh, ground_shv,                   &
                  veg_water, veg_temp, can_temp, can_shv,                 &
                  wshed, qwshed, transp, stom_resist,                     &
                  hxfergc, wxfergc, hxfersc, wxfersc, hxferca,            &
                  hxfervc, wxfervc, rdi, rb, time8,                       &
                  lsl, ed_transp, ed_patch)

use leaf_coms,   only: nzg, soil_rough, dt_leaf,  &
                       slpots, slmsts, slbs, kroot, rcmin, soilcp, dslz

use consts_coms, only: cp, vonk, eps_vap, alvl, cliq, cice, alli, rvap

use ed_structure_defs
use misc_coms, only: io6
use mem_leaf,  only: itab_wl

implicit none

integer, intent(in)  :: iwl           ! index of current land cell
integer, intent(in)  :: nlev_sfcwater   ! # active levels of surface water
integer, intent(in)  :: ntext_soil(nzg) ! soil textural class
integer, intent(in)  :: leaf_class      ! leaf class (vegetation class)
integer, intent(in)  :: lsl
integer, intent(out) :: ktrans          ! k index of soil layer supplying transp

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
real, intent(in) :: can_depth   ! canopy depth for heat and vap capacity [m]
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

real, intent(out) :: hxfergc ! soil-to-can_air heat xfer this step [J/m^2]
real, intent(out) :: wxfergc ! soil-to-can_air vap xfer this step [kg_vap/m^2]
real, intent(out) :: hxfersc ! sfc_water-to-can_air heat xfer this step [J/m^2]
real, intent(out) :: wxfersc ! sfc_water-to-can_air vap xfer this step [kg_vap/m^2]
real, intent(out) :: wshed   ! water shed from veg this LEAF timestep [kg/m^2]
real, intent(out) :: qwshed  ! water energy shed from veg this LEAF timestep [J/m^2]
real, intent(out) :: transp  ! transpiration xfer this LEAF timestep [kg_vap/m^2]
real, intent(out) :: hxferca ! can_air-to-atm heat xfer this step [J/m^2]
real, intent(out) :: hxfervc ! veg-to-can_air heat xfer this step [J/m^2]
real, intent(out) :: wxfervc ! veg-to-can_air vapor xfer this step [kg_vap/m^2]
real, intent(out) :: rdi     ! (soil or surface water)-to-can_air conductance [m/s]
real, intent(out) :: rb      ! veg-to-can_air resistance [s/m]
real(kind=8), intent(in) :: time8   ! model time [s]

real, dimension(nzg) :: ed_transp ! transpired water from each soil level; ED2 cells only [kg/m^2]
type(patch), target, optional :: ed_patch

! Local parameters

real, parameter :: exar = 2.5  ! for computing rasveg
real, parameter :: covr = 2.16 ! scaling tai value for computing wtveg
real, parameter :: c1 = 116.6  ! for computing rb

!     Note: c1=261.5*sqrt((1.-exp(-2.*exar))/(2.*exar)) 
!     from Lee's dissertation, Eq. 3.36.  The factor of 261.5 is
!     100 * ln((h-d)/zo) / vonk   where d = .63 * h and zo = .13 * h.
!     The factor of 100 is 1/L in Eq. 3.37.  Thus, c1 * ustar is the
!     total expression inside the radical in Eq. 3.37.
!bob      parameter(exar=3.5,covr=2.16,c1=98.8)

! Intercept, slope parameters for stomatal conductance factors

real, parameter :: brad = 196.   , srad = .047    ! for s/w radiative flux
real, parameter :: btlo = 281.5  , stlo = .26     ! for low canopy temperature
real, parameter :: bthi = 310.1  , sthi = -.124   ! for high canopy temperature
real, parameter :: bvpd = 4850.  , svpd = -.0051  ! for vapor pressure deficit
real, parameter :: bswp = -1.07e6, sswp = 7.42e-6 ! for soil water potential

! Local variables

integer :: k        ! loop index over soil layers
integer :: nts      ! soil textural class for current soil layer

real :: factv       ! for computing rasveg
real :: aux         ! for computing rasveg
real :: rasveg      ! full-veg value of rd [s/m]
real :: wtveg       ! weighting of rasveg in computing rd
real :: rasgnd      ! not used
real :: c3          ! veg_sfc-to-can_air vapor density difference [kg_vap/m^3]
real :: fracliqv    ! fraction of veg surface water in liquid phase
real :: fthi        ! high-temp environ factor for stomatal resist
real :: ftlo        ! low-temp environ factor for stomatal resist
real :: frad        ! s/w radiative environ factor for stomatal resist
real :: fswp        ! soil water potential environ factor for stomatal resist
real :: fvpd        ! vap press deficit environ factor for stomatal resist
real :: qwtot       ! total internal energy of veg plus new precip on veg [J/m^2]
real :: esat_veg    ! saturation vapor pressure at veg temp [Pa]
real :: veg_rhovs   ! saturation vapor density at veg temp [kg_vap/m^3]
real :: veg_rhovsp  ! saturation vapor density gradient at veg temp [kg_vap/(K m^3)]
real :: e_can       ! vapor pressure of canopy air [Pa]
real :: e_leaf      ! vapor pressure at leaf surface [Pa]
real :: vpd         ! vapor pressure deficit across stomata [Pa]
real :: rc_inf      ! asymptotic stomatal resistance at current veg environment [s/m]
real :: sigmaw      ! fractional coverage of leaf surface by veg_water
real :: slai        ! effective veg lai uncovered by surface water (snowcover)
real :: stai        ! effective veg tai uncovered by surface water (snowcover)
real :: slpotv      ! soil water potential [m]
real :: swp         ! soil water potential factor for stomatal control
real :: tvegc       ! vegetation temperature (deg C)
real :: tvegk       ! vegetation temperature (K)
real :: wtroot      ! not used
real :: zognd       ! soil roughness height [m]
real :: zoveg       ! vegetation roughness height [m]
real :: zdisp       ! vegetation displacement height remaining after snowcover [m]
real :: zveg        ! vegetation height remaining after snowcover [m]
real :: wxfer       ! (saturated soil or sfc water)-to-can_air vap xfer this step [kg_vap/m^2]
real :: transp_test ! test value of transpiration flux [kg_vap/(m^2 s)]
real :: rc          ! stomatal resistance [s/m]

real :: canair
real :: canhcap
real :: f1
real :: f2
real :: f3
real :: a1
real :: a2
real :: a3
real :: a4
real :: b1
real :: b2
real :: b3
real :: b4
real :: b5
real :: b6
real :: b7
real :: b8
real :: b9
real :: evapotransp

real, external :: rhovsil ! function to compute sat vapor density
real, external :: eslf    ! function to compute sat vapor pressure

! Can_air-to-atm heat xfer in [J/m^2]

hxferca = cp * sxfer_t

! Canopy air (vapor) and heat capacity

canair = rhos * can_depth
canhcap = cp * canair

! Initialize wshed = qwshed = 0

wshed  = 0.
qwshed = 0.

! Check whether this land cell has exposed vegetation

if ((.not.present(ed_patch)).and.(veg_tai < .1 .or. snowfac > .9)) then

! If the TAI is very small or if snow mostly covers the vegetation, 
! then BYPASS THE VEGETATION COMPUTATIONS.

! Aerodynamic conductance for bare soil or snow based on Garratt.

   rdi = .2 * ustar

! Set transpiration to zero

   transp = 0.
   ktrans = 0

! Check if any surface water layers currently exist

   if (nlev_sfcwater >= 1) then

! If surface water is present, compute heat and vapor xfer between 
! surface water and canopy air.

      hxfergc = 0.
      hxfersc = cp * rhos * dt_leaf * rdi * (sfcwater_tempk - can_temp)

      wxfergc = 0.
      wxfersc = min(sfcwater_mass,                                 &
                    rhos * dt_leaf * rdi * (surface_ssh - can_shv) )

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw01 ',wxfersc,sfcwater_mass,rhos,dt_leaf,rdi
!   write(io6,*) 'sfcw02 ',surface_ssh, can_shv
!endif

! Update canopy temperature and vapor specific humidity, and set vegetation
! temperature to canopy air value

!if (itab_wl(iwl)%iwglobe == 821 .or. itab_wl(iwl)%iwglobe == 822 .or.  &
!     iwl == 821 .or. iwl == 822) then
!   write(io6,*) ' '
!   write(io6,*) 'can0 ',iwl,can_temp,dble(can_temp)
!endif

      can_temp = can_temp + (hxfersc - hxferca) / canhcap
      can_shv  = can_shv  + (wxfersc - sxfer_r) / canair
      veg_temp = can_temp

!if (itab_wl(iwl)%iwglobe == 821 .or. itab_wl(iwl)%iwglobe == 822 .or.  &
!     iwl == 821 .or. iwl == 822) then
!   write(io6,*) ' '
!   write(io6,*) 'can1 ',can_temp,dble(can_temp)
!   write(io6,*) 'can2 ',hxfersc,dble(hxfersc)
!   write(io6,*) 'can3 ',hxferca,dble(hxferca)
!   write(io6,*) 'can4 ',canhcap,dble(canhcap)
!   write(io6,*) 'can5 ',rdi,dble(rdi)
!   write(io6,*) 'can6 ',dt_leaf,dble(dt_leaf)
!   write(io6,*) 'can7 ',rhos,dble(rhos)
!   write(io6,*) 'can8 ',cp,dble(cp)
!   write(io6,*) 'can9 ',sfcwater_tempk,dble(sfcwater_tempk)
!   write(io6,*) ' '
!endif

      return

   else

! If surface water is absent, compute heat xfer between soil and canopy air

      hxfergc = cp * rhos * dt_leaf * rdi * (soil_tempk - can_temp)
      hxfersc = 0.

! Compare saturation vapor specific humidity at soil temperature against 
! canopy air vapor specific humidity

      if (surface_ssh < can_shv) then

! If saturation vapor specific humidity at soil temperature is less than 
! canopy air vapor specific humidity, compute dew formation contribution 
! to surface water

         wxfergc = 0.
         wxfersc = rhos * dt_leaf * rdi * (surface_ssh - can_shv)

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw03 ',wxfersc,rhos,dt_leaf,rdi
!   write(io6,*) 'sfcw04 ',surface_ssh, can_shv
!endif

      elseif (ground_shv > can_shv) then

! Else, if equilibrium soil vapor specific humidity is greater than 
! canopy air vapor specific humidity, compute evaporation from soil

         wxfergc = max(0.,rhos * dt_leaf * rdi * (ground_shv - can_shv))
         wxfersc = 0.

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw05 ',wxfergc,rhos,dt_leaf,rdi
!   write(io6,*) 'sfcw06 ',ground_shv, can_shv
!endif

      else

! If neither of the above is true, both vapor fluxes are zero.

         wxfergc = 0.
         wxfersc = 0.

      endif

! Update canopy temperature and vapor specific humidity, and set vegetation
! temperature to canopy air value


! if (itab_wl(iwl)%iwglobe == 821 .or. itab_wl(iwl)%iwglobe == 822 .or.  &
!      iwl == 821 .or. iwl == 822) then
!    write(io6,*) ' '
!    write(io6,*) 'can10 ',iwl,can_temp,dble(can_temp)
! endif

      can_temp = can_temp + (hxfergc + hxfersc - hxferca) / canhcap
      can_shv  = can_shv  + (wxfergc + wxfersc - sxfer_r) / canair
      veg_temp = can_temp


! if (itab_wl(iwl)%iwglobe == 821 .or. itab_wl(iwl)%iwglobe == 822 .or.  &
!      iwl == 821 .or. iwl == 822) then
!    write(io6,*) ' '
!    write(io6,*) 'can11 ',can_temp,dble(can_temp)
!    write(io6,*) 'can12 ',hxfergc,dble(hxfergc)
!    write(io6,*) 'can13 ',hxfersc,dble(hxfersc)
!    write(io6,*) 'can14 ',hxferca,dble(hxferca)
!    write(io6,*) 'can15 ',canhcap,dble(canhcap)
!    write(io6,*) 'can16 ',rdi,dble(rdi)
!    write(io6,*) 'can17 ',dt_leaf,dble(dt_leaf)
!    write(io6,*) 'can18 ',rhos,dble(rhos)
!    write(io6,*) 'can19 ',cp,dble(cp)
!    write(io6,*) ' '
! endif

      return

   endif

else

! If vegetation is sufficiently abundant and not covered by snow,
! COMPUTE CANOPY XFER WITH VEGETATION INFLUENCE

! Compute ground-canopy resistance rd.  Assume zognd not affected by snow.
! Assume (zoveg,zdisp) decrease linearly with snow depth, attaining
! the values (zognd,0) when veg covered.

   zognd = soil_rough
   zoveg = veg_rough * (1. - snowfac) + zognd * snowfac
   zveg  = veg_height * (1. - snowfac)
   zdisp = zveg * .63
!bob      rasgnd = log(zts / zognd) * log((zdisp + zoveg) / zognd)
!bob     +      / (vonk * vonk * vels)

! Aerodynamic resistance (rd) between surface and canopy air are weighted
! between areas without and with vegetation.

!bob   factv  = log(zts / zoveg) / (vonk * vonk * vels)
   factv  = 1. / (vonk * ustar)
   aux    = exp(exar * (1. - (zdisp + zoveg) / zveg))
   rasveg = factv * zveg / (exar * (zveg - zdisp)) * (exp(exar) - aux)
   wtveg  = max(0.,min(1., 1.1 * veg_tai / covr))
   rdi    = ustar / (5. * (1. - wtveg) + ustar * rasveg * wtveg)
   
! Check if any surface water layers currently exist

   if (nlev_sfcwater >= 1) then

! If surface water is present, compute heat and vapor xfer between 
! surface water and canopy air.

      hxfergc = 0.
      hxfersc = cp * rhos * dt_leaf * rdi * (sfcwater_tempk - can_temp)

      wxfergc = 0.
      wxfersc = min(sfcwater_mass,                                 &
                    rhos * dt_leaf * rdi * (surface_ssh - can_shv) )

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw07 ',wxfersc,wxfergc,rhos,dt_leaf,rdi
!   write(io6,*) 'sfcw08 ',surface_ssh, can_shv
!endif

   else

! If surface water is absent, compute heat xfer between soil and canopy air

      hxfergc = cp * rhos * dt_leaf * rdi * (soil_tempk - can_temp)
      hxfersc = 0.

! Compare saturation vapor specific humidity at soil temperature against 
! canopy air vapor specific humidity

      if (surface_ssh < can_shv) then

! If saturation vapor specific humidity at soil temperature is less than 
! canopy air vapor specific humidity, compute dew formation contribution 
! to surface water

         wxfergc = 0.
         wxfersc = rhos * dt_leaf * rdi * (surface_ssh - can_shv)

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw09 ',wxfersc,wxfergc,rhos,dt_leaf,rdi
!   write(io6,*) 'sfcw010 ',surface_ssh, can_shv
!endif

      elseif (ground_shv > can_shv) then

! Else, if equilibrium soil vapor specific humidity is greater than 
! canopy air vapor specific humidity, compute evaporation from soil

         wxfergc = max(0.,rhos * dt_leaf * rdi * (ground_shv - can_shv))
         wxfersc = 0.

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw011 ',wxfersc,wxfergc,rhos,dt_leaf,rdi
!   write(io6,*) 'sfcw012 ',ground_shv, can_shv
!endif

      else

! If neither of the above is true, both vapor fluxes are zero.

         wxfergc = 0.
         wxfersc = 0.

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw013 ',wxfersc,wxfergc
!endif

      endif

   endif

! Here, LEAF3 and ED2 branch off.

   if(.not.present(ed_patch))then
! TAI and LAI reduced by ground snowcover (for low vegetation)

   stai = veg_tai * (1. - snowfac)
   slai = veg_lai * (1. - snowfac)

! Evaluate vegetation Celsius temperature

   tvegc = veg_temp - 273.15

! Check if precipitation is occurring

   if (pcpg > 1.e-12) then
   
! If precipitation, add intercepted mass and energy to vegetation surface
!   (Ignore energy of water already on vegetation surface)

      veg_water = veg_water + pcpg * vf
      qwtot = hcapveg * tvegc + qpcpg * vf

! Compute equilbrium temperature of veg + precipitation

      call qwtk(qwtot,pcpg * vf,hcapveg,veg_temp,fracliqv)
      tvegc = veg_temp - 273.15
      
! Shed any excess intercepted precipitation and its energy

      if (veg_water > .2 * stai) then
         wshed = veg_water - .2 * stai

         if (fracliqv <= .0001) then
            qwshed = cice * tvegc * wshed
         else
            qwshed = (cliq * tvegc + fracliqv * alli) * wshed
         endif

         veg_water = veg_water - wshed
      endif
      
   endif

! Vegetation saturation vapor density and gradient

   veg_rhovs = rhovsil(tvegc)
   veg_rhovsp = rhovsil(tvegc+1.) - veg_rhovs

! Compute veg-canopy resistance rb.  Note that rb and rc are both defined
! WITHOUT the LAI factor; the LAI factor is included later in computing
! xfers that involve rb and/or rc.

   rb = (1. + .5 * veg_tai) / (.01 * sqrt(ustar * c1))

! Initialize soil water potential to -600 m (later to be converted to Pascals)
! prior to search for wettest soil layer in root zone.  This value is low 
! enough to effectively prevent transpiration.

   swp = -600.

! Initialize ktrans to zero prior to search for wettest soil layer in root zone
   
   ktrans = 0

! Loop over soil layers in the root zone

   do k = kroot(leaf_class),nzg
      nts = ntext_soil(k)

! Soil water potential

      slpotv = slpots(nts) * (slmsts(nts) / soil_water(k)) ** slbs(nts)

! Multiply by liquid fraction (ice is unavailable for transpiration)

      slpotv = slpotv * soil_fracliq(k)

! Find layer in root zone with highest slpotv AND soil_water above minimum soilcp
! Set ktrans to this layer

      if (slpotv > swp .and. soil_water(k) > soilcp(nts)) then
         swp = slpotv
         ktrans = k
      endif
   enddo

  swp = swp * 9810. ! convert from meters to Pascals (hydrostatic eqn for water)

! Bypass stomatal resistance computation if root zone is too dry
   
   if (ktrans < 1) then

      rc = 1.e18

   else

! Compute saturation vapor pressure at veg temp

      esat_veg  = eslf(tvegc)

! Compute vapor pressure in canopy air from equation of state

      e_can = can_shv * rhos * rvap * can_temp

! Compute vapor pressure at leaf surface using rc from previous timestep

      rc = stom_resist
      e_leaf = (rb * esat_veg + rc * e_can) / (rb + rc)
      vpd = max(0.,esat_veg - e_leaf)

! Evaluate 5 environmental factors and new rc

      ftlo = 1. + exp(-stlo * (veg_temp - btlo))
      fthi = 1. + exp(-sthi * (veg_temp - bthi))
      frad = 1. + exp(-srad * (rshort   - brad))
      fswp = 1. + exp(-sswp * (swp      - bswp))
      fvpd = 1. + exp(-svpd * (vpd      - bvpd))

! Compute asymptotoc value of stomatal resistance based on environmental factors

      rc_inf = ftlo * fthi * frad * fvpd * fswp * rcmin(leaf_class)

! Update stomatal conductance assuming 15-minute response time 
! (requires dt_leaf <= 900.) 

      rc = 1. / (1. / rc + .0011 * dt_leaf * (1. / rc_inf - 1. / rc))

! Limit maximum transpiration to be <= 500 W/m^2 by increasing rc if necessary.

      transp_test = alvl * veg_lai * (veg_rhovs - rhos * can_shv) / (rb + rc)
      if (transp_test > 500.) then
         rc = (rb + rc) * transp_test * .002 - rb
      endif      
      
   endif
      
   stom_resist = rc

!-------------------------------------------------------------------------------
! Begin implicit exchange of heat and moisture between vegetation and canopy air
!-------------------------------------------------------------------------------

! veg_water fractional coverage
   
   sigmaw = min(1.,(veg_water / (.2 * stai)) ** .66667)

! auxiliary quantities

   f1 = (hxfergc + hxfersc - hxferca) / canhcap
   f2 = (wxfergc + wxfersc - sxfer_r) / canair
   f3 = dt_leaf * (rshort_v + rlong_v) / hcapveg

   a1 = dt_leaf * 2. * stai / rb                  ! vap xfer coef (dew formation)
   a2 = cp * rhos * a1                            ! heat xfer coef
   a3 = a1 * sigmaw                               ! vap xfer coef (veg_water evap)
   a4 = dt_leaf * slai * (1.-sigmaw) / (rb + rc)  ! vap xfer coef (transp)

   b1 = canhcap + a2
   b2 = canhcap * a2
   b3 = hcapveg * b1 + b2
   b4 = b2 * (can_temp + f1)
   b5 = b1 * hcapveg * (veg_temp + f3)
   b6 = veg_rhovs + veg_rhovsp * ((b5 + b4) / b3 - veg_temp)
   b7 = veg_rhovsp * b1 * alvl / b3
   b8 = b7 + 1. / can_depth
   b9 = b6 - rhos * (can_shv + f2)

! First, assume case of condensation/deposition of vapor onto veg

   evapotransp = a1 * b9 / (1. + a1 * b8)
   
! Now test this assumption

   if (evapotransp <= 0.) then

! Assumption was correct.  Compute individual transpiration and evaporation

      wxfervc = evapotransp
      transp = 0.

   elseif (veg_water > 1.e-12) then

! Assumption was incorrect.  Now, consider the case where veg_water is present.

! Assume that veg_water does not all evaporate

      evapotransp = (a3 + a4) * b9 / (1. + (a3 + a4) * b8)
      wxfervc = evapotransp * a3 / (a3 + a4)

! Now test this assumption

      if (wxfervc <= veg_water) then

! Assumption was correct.  Compute transpiration

         transp = evapotransp - wxfervc

      else

! Assumption was incorrect.  All of veg_water evaporates this leaf step
!  Compute combined and individual transpiration and evaporation

         evapotransp = (veg_water + a4 * b9) / (1. + a4 * b8)
         wxfervc = veg_water
         transp = evapotransp - wxfervc
      
      endif

      if(a4 == 0.0)transp = 0.0

   else

! Original condensation/deposition assumption was incorrect, and there is 
!   no veg_water.

      veg_water = 0.
      evapotransp = a4 * b9 / (1. + a4 * b8)
      wxfervc = 0.
      transp = evapotransp

   endif

! Compute remaining quantities

   can_shv = can_shv + f2 + evapotransp / canair
   veg_temp = (b5 + b4 - b1 * alvl * evapotransp) / b3
   hxfervc = (b2 * veg_temp - b4) / b1
   can_temp = can_temp + f1 + hxfervc / canhcap
   veg_water = veg_water - wxfervc

   else

      ktrans = 0 
      transp = 0.
      ed_transp(:) = 0.

      call ed_canopy_update(ed_patch, vels, rhos, prss, pcpg, qpcpg,   &
           wshed, qwshed, canair, canhcap, dt_leaf, hxfergc, hxferca,  &
           wxfergc, hxfersc, wxfersc, sxfer_r, ed_transp, ntext_soil,  &
           soil_water, soil_fracliq, lsl)

   endif

endif

return
end subroutine canopy

!===============================================================================

subroutine sfcwater(iwl, nlev_sfcwater,ntext_soil,                       &
                    soil_rfactor, soil_water, soil_energy,                 &
                    sfcwater_mass, sfcwater_energy, sfcwater_depth,        &
                    soil_tempk, soil_fracliq, sfcwater_tempk,              &
                    sfcwater_fracliq, rshort_s,                            &
                    hxfersc, wxfersc, rlong_g, rlong_s,                    &
                    pcpg, qpcpg, dpcpg, wshed, qwshed,                     &
                    veg_fracarea, lsl, time8                               )

use leaf_coms, only: nzg, nzs, dt_leaf, slz, dslz, dslzi, dslzo2,  &
                     slmsts, soilcond0, soilcond1, soilcond2, slcpd
                       
use consts_coms, only: alvi, cice, cliq, alli
use misc_coms, only: io6

implicit none

integer, intent(in)    :: iwl              ! current land cell index
integer, intent(inout) :: nlev_sfcwater      ! # of active sfc water levels
integer, intent(in)    :: ntext_soil   (nzg) ! soil textural class

real, intent(out)   :: soil_rfactor    (nzg) ! soil thermal resistance [K m^2/W]
real, intent(inout) :: soil_water      (nzg) ! soil water content [water_vol/total_vol]
real, intent(inout) :: soil_energy     (nzg) ! soil internal energy [J/m^3]
real, intent(inout) :: sfcwater_mass   (nzs) ! surface water mass [kg/m^2]
real, intent(inout) :: sfcwater_energy (nzs) ! surface water energy [J/kg]
real, intent(inout) :: sfcwater_depth  (nzs) ! surface water depth [m]
real, intent(inout) :: soil_tempk      (nzg) ! soil temperature [K]
real, intent(inout) :: soil_fracliq    (nzg) ! fraction of soil water in liq phase
real, intent(inout) :: sfcwater_tempk  (nzs) ! surface water temp [K]
real, intent(inout) :: sfcwater_fracliq(nzs) ! fraction of sfc water in liq phase
real, intent(inout) :: rshort_s        (nzs) ! s/w net rad flux to sfc water [W/m^2]
real, intent(in)    :: hxfersc ! sfc_water-to-can_air heat xfer this step [J/m^2]
real, intent(in)    :: wxfersc ! sfc_water-to-can_air vap xfer this step [kg_vap/m^2]
real, intent(in)    :: rlong_g ! l/w net rad flux to soil [W/m^2]
real, intent(in)    :: rlong_s ! l/w net rad flux to sfc water [W/m^2]
real, intent(in)    :: pcpg    ! new pcp amount this leaf timestep [kg/m^2]
real, intent(in)    :: qpcpg   ! new pcp energy this leaf timestep [J/m^2]
real, intent(in)    :: dpcpg   ! new pcp depth this leaf timestep [m]
real, intent(in)    :: wshed   ! water shed from veg this timestep [kg/m^2]
real, intent(in)    :: qwshed  ! water energy shed from veg this timestep [J/m^2]
real, intent(in)    :: veg_fracarea ! veg fractional area
integer, intent(in) :: lsl

! Local variables

real :: hxfers        (nzs+1) ! sfcwater heat xfer [J/m2] 
real :: sfcwater_rfactor(nzs) ! sfcwater thermal resistivity [K m^2/W]
real :: mass_new        (nzs) ! mass of new set of sfcwater layers [kg/m^2]
real :: energy_new      (nzs) ! energy of new set of sfcwater layers [J/kg]
real :: depth_new       (nzs) ! depth of new set of sfcwater layers [m]
real :: rshort_snew     (nzs) ! s/w rad flux to new set of sfcwater layers [W/m^2]
real :: energy_per_m2   (nzs) ! sfcwater energy per square meter [J/m^2]

integer :: k         ! vertical index over sfcwater layers
integer :: kold      ! vertical index of adjacent lower sfcwater layer
integer :: icelayers ! # of sfcwater layers that contain some ice
integer :: maxlayers ! max allowed # of sfcwater layers (1 if none contain ice)
integer :: nlev_new  ! new number of sfcwater layers after adjustment

real :: hxfergs   ! energy transfer from soil to sfcwater this step [J/m^2]
real :: rfac      ! bounded sfcwater thermal resistivity at k=1 [K m^2/W]
real :: snden     ! sfcwater density [kg/m^3]
real :: vegfracc  ! 1 minus veg fractional area
real :: wfree     ! free liquid in sfcwater layer that can percolate out [kg/m^2]
real :: qwfree    ! energy carried by wfree [J/m^2]
real :: dwfree    ! depth carried by wfree [m]
real :: fracstep  ! ratio of leaf timestep to snow density exponential decay time
real :: totsnow   ! sum of mass over sfcwater layers [kg/m^2]
real :: wt        ! sfcwater(1) + soil(nzg) water masses (impl balance) [kg/m^2]
real :: qwt       ! sfcwater(1) + soil(nzg) energies (impl balance) [J/m^2]
real :: soilhcap  ! soil(nzg) heat capacity [J/(m^2 K)]
real :: soilcap   ! capacity of top soil layer to accept surface water [kg/m^2]
real :: sndenmin  ! minimum sfcwater density [kg/m^3]
real :: wtnew     ! weight for new sfcwater layer when adjusting layer thickness
real :: wtold     ! weight for old sfcwater layer when adjusting layer thickness
real :: dwtold    ! change in wtold for partial mass transfer from old layer
real :: wdiff     ! change in sfcwater mass when adjusting layer thickness [kg/m^2]
real :: soilcond  ! soil thermal conductivity [W/(K m)]
real :: waterfrac ! soil water fraction in soil layer [vol_water/vol_total]
real :: tempk     ! Kelvin temperature [K]
real :: fracliq   ! fraction of water in liquid phase returned from qwtk
real :: flmin     ! lower bound on sfcwater_fracliq(1) in balance with soil
real :: flmax     ! upper bound on sfcwater_fracliq(1) in balance with soil
real :: specvol   ! specific volume of sfcwater involved in vapor xfer [m^3/kg]
real(kind=8), intent(in) :: time8   ! model time [s]

! Local parameters

real, parameter :: sndenmax = 1000.   ! max sfcwater density [kg/m^3]
real, parameter :: snowmin = 11.       ! min sfcwater layer mass with multiple layers [kg/m^2] 
real, parameter :: snowmin_expl = 10.  ! min sfcwater mass for explicit heat xfer [kg/m^2]
real, parameter :: rfac_snowmin = .01 ! min sfcwater rfactor [K m^2/W]

real, save, dimension(10,10) :: thick  ! snowlayer thickness scaling factor

data thick(1:10, 1)/  1., .00, .00, .00, .00, .00, .00, .00, .00, .00/
data thick(1:10, 2)/ .50, .50, .00, .00, .00, .00, .00, .00, .00, .00/
data thick(1:10, 3)/ .25, .50, .25, .00, .00, .00, .00, .00, .00, .00/
data thick(1:10, 4)/ .17, .33, .33, .17, .00, .00, .00, .00, .00, .00/
data thick(1:10, 5)/ .10, .20, .40, .20, .10, .00, .00, .00, .00, .00/
data thick(1:10, 6)/ .07, .14, .29, .29, .14, .07, .00, .00, .00, .00/
data thick(1:10, 7)/ .05, .09, .18, .36, .18, .09, .05, .00, .00, .00/
data thick(1:10, 8)/ .03, .07, .13, .27, .27, .13, .07, .03, .00, .00/
data thick(1:10, 9)/ .02, .04, .09, .18, .34, .18, .09, .04, .02, .00/
data thick(1:10,10)/ .02, .03, .06, .13, .26, .26, .13, .06, .03, .02/

!if (iwl == 102) then
!   write(io6,*) 'sfcw0 ',nlev_sfcwater
!endif
         
! Compute soil heat resistance times HALF layer depth (soil_rfactor).

do k = 1,nzg
   waterfrac = soil_water(k) / slmsts(ntext_soil(k))
   soilcond =        soilcond0(ntext_soil(k))  &
      + waterfrac * (soilcond1(ntext_soil(k))  &
      + waterfrac *  soilcond2(ntext_soil(k))  )
   soil_rfactor(k) = dslzo2(k) / soilcond
enddo

! Check whether surface water was present at the beginning of this leaf step

if (nlev_sfcwater > 0) then

! Surface water was present.

! Compute snow heat resistance times HALF layer depth (sfcwater_rfactor).
! Sfcwater_tempk(k) should be correctly balanced value at this point, so 
! sfcwater_rfactor(k) should have correct value.  Formula applies to snow,
! so limit temperature to no greater than 273.15.

   do k = 1,nlev_sfcwater
      snden = sfcwater_mass(k) / sfcwater_depth(k)
      tempk = min(273.15,sfcwater_tempk(k))

      sfcwater_rfactor(k) = .5 * sfcwater_depth(k)  &
         / (1.093e-3 * exp(.028 * tempk) *          &
         (.030 + snden * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))
   enddo

! Zero out sfcwater internal heat transfer array at bottom and top surfaces.
! Energy transfer at bottom and top are applied separately.

   hxfers(1) = 0.
   hxfers(nlev_sfcwater+1) = 0. 

! Compute internal sfcwater energy xfer if at least two layers exist [J/m2]

   if (nlev_sfcwater >= 2) then	
      do k = 2,nlev_sfcwater
         hxfers(k) = dt_leaf * (sfcwater_tempk(k-1) - sfcwater_tempk(k))  &
                   / (sfcwater_rfactor(k-1) + sfcwater_rfactor(k))      
      enddo
   endif

! Diagnose sfcwater energy per square meter, and add contributions from 
! internal transfers of sensible heat and vapor (latent heat) and from
! shortwave radiative transfer.  This excludes energy transfer from internal
! gravitational draining of free water mass.

   do k = 1,nlev_sfcwater
      energy_per_m2(k) = sfcwater_energy(k) * sfcwater_mass(k)  &
         + hxfers(k) - hxfers(k+1) + dt_leaf * rshort_s(k)
         
!if (iwl == 102) then
!   write(io6,*) 'sfcw1 ',k,nlev_sfcwater,energy_per_m2(k),sfcwater_energy(k)
!   write(io6,*) 'sfcw2 ',sfcwater_mass(k),  &
!                 hxfers(k),hxfers(k+1),dt_leaf,rshort_s(k)
!endif
         
   enddo

! Evaluate conductive heat transfer between top soil layer and bottom sfcwater
! layer.  Impose minimum value on sfcwater_rfactor(1).  If minimum applies, 
! energy will be implicitly re-balanced later.  Apply heat transfer to bottom
! sfcwater layer and top soil layer.

   rfac = max(rfac_snowmin,sfcwater_rfactor(1))

   hxfergs = dt_leaf * (soil_tempk(nzg) - sfcwater_tempk(1))   &
           / (soil_rfactor(nzg) + rfac)

   energy_per_m2(1) = energy_per_m2(1) + hxfergs

   soil_energy(nzg) = soil_energy(nzg) - hxfergs * dslzi(nzg)

! Apply longwave radiative transfer and sfcwater-to-can_air sensible heat
! transfer to top sfcwater layer

   energy_per_m2(nlev_sfcwater) = energy_per_m2(nlev_sfcwater)  &
      + dt_leaf * rlong_s - hxfersc

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw3 ',energy_per_m2(nlev_sfcwater),rlong_s,hxfersc
!endif

endif

! If no sfcwater layers exist and there are net positive contributions to 
! sfcwater from precipitation, shedding of water from vegetation, and 
! vapor flux with canopy air, create a new sfcwater layer and initialize 
! prognostic sfcwater fields to zero.

if (nlev_sfcwater == 0 .and.  &
    pcpg * (1. - veg_fracarea) + wshed * veg_fracarea - wxfersc > 1.e-9) then

   nlev_sfcwater = 1

   sfcwater_mass  (1) = 0.
   sfcwater_energy(1) = 0.
   energy_per_m2  (1) = 0.
   sfcwater_depth (1) = 0.
endif

! Return if no sfcwater layers now exist

if (nlev_sfcwater < 1) return

! Sfcwater layers do exist

! Apply mass, energy, and depth contributions to top sfcwater layer from
! precipitation, shedding of water from vegetation, and vapor flux with 
! canopy air.  Get value for specific volume of sfcwater involved in vapor xfer.

specvol = .001
if (wxfersc > 0.) specvol =  &
   sfcwater_depth(nlev_sfcwater) / sfcwater_mass(nlev_sfcwater)

sfcwater_mass(nlev_sfcwater) = sfcwater_mass(nlev_sfcwater)  &
   + pcpg  * (1. - veg_fracarea)                             &
   + wshed * veg_fracarea                                    &
   - wxfersc

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw4 ',energy_per_m2(nlev_sfcwater)
!endif

energy_per_m2(nlev_sfcwater) = energy_per_m2(nlev_sfcwater)  &
   + qpcpg   * (1. - veg_fracarea)                           &
   + qwshed  * veg_fracarea                                  &
   - wxfersc * alvi 

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw4.1 ',energy_per_m2(nlev_sfcwater),qpcpg,qwshed,wxfersc
!   write(io6,*) 'sfcw5 ',veg_fracarea,alvi
!endif

sfcwater_depth(nlev_sfcwater) = sfcwater_depth(nlev_sfcwater)  &
   + dpcpg * (1. - veg_fracarea)                               &
   + wshed * veg_fracarea * .001                               &
   - wxfersc * specvol

! If nlev_sfcwater = 1, check whether any sfcwater mass still remains 
! (after possibly having all evaporated into canopy).  If mass is below
! threshold value, set sfcwater quantities to zero and return.

if (nlev_sfcwater == 1 .and. sfcwater_mass(nlev_sfcwater) < 1.e-9) then

   nlev_sfcwater = 0

   sfcwater_mass(1:nzs)   = 0.
   sfcwater_energy(1:nzs) = 0.
   sfcwater_depth(1:nzs)  = 0.

   return
   
endif

! Prepare to transfer water downward through snow layers by percolation.
! Fracliq is the fraction of liquid in the snowcover or surface water.
! wfree [kg/m^2] is the quantity of that liquid that is free (not attached to
! snowcover) and therefore available to drain into the layer below.
! soilcap [kg/m^2] is the amount of water that can fit in the unfilled pore
! space of the top soil layer.  wfree in the lowest sfcwater layer is limited
! by this value.

! First, prepare to sum sfcwater mass over all existing layers

totsnow = 0

! Loop downward through all existing sfcwater layers beginning with top layer

do k = nlev_sfcwater,1,-1

! Diagnose sfcwater density.  Make sure sfcwater_depth is not too near zero.

   sfcwater_depth(k) = max(sfcwater_depth(k), .001 * sfcwater_mass(k))

   snden = sfcwater_mass(k) / sfcwater_depth(k)

! Assume that as snow ages on ground, its density difference with a limiting
! maximum value (currently set to 400 kg/m^3) decays exponentially (with a
! decay time currently set to about 3 weeks).  If sfcwater density is less
! than this limiting value, apply the density increase for this timestep.

! This formulation and decay constants are very crude approximations to a few
! widely variable observations of snowpack density and are intended only as a 
! rough representation of the tendency for snowcover to compress with time.  
! A better formulation that accounts for several environmental factors would
! be desirable here.

   if (snden < 400.) then
      fracstep = .5e-6 * dt_leaf  ! .5e-6 is inverse decay time scale
      snden = snden * (1. - fracstep) + 400. * fracstep
      sfcwater_depth(k) = sfcwater_mass(k) / snden
   endif

! Diagnose sfcwater temperature and liquid fraction now that new mass and
! energy contributions have been applied.  Use qwtk instead of qtk in case
! sfcwater_mass(k) is too small; "dryhcap" = 100 is very small value    

   call qwtk(energy_per_m2(k),sfcwater_mass(k),100.,  &
             sfcwater_tempk(k),sfcwater_fracliq(k))

! If this is bottom layer, diagnose sfcwater_rfactor.  Since sfcwater_tempk(k)
! may not be a stable computation at this point, assume that tempk = 0, which
! gives the minimum rfactor for a given density.

   if (k == 1) then
      tempk = 273.15

      sfcwater_rfactor(k) = .5 * sfcwater_depth(k)  &
         / (1.093e-3 * exp(.028 * tempk) *          &
         (.030 + snden * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))

   endif

! If this is bottom layer and sfcwater rfactor is less than minimum stable 
! value, bring bottom sfcwater and top soil layer into thermal equilibrium 
! by exchanging heat between them.

   if (k == 1 .and.   &
       (sfcwater_mass(1) < snowmin_expl .or.  &
        sfcwater_rfactor(1) < rfac_snowmin)) then

! Combined sfcwater and soil water mass per square meter

      wt = sfcwater_mass(1) + soil_water(nzg) * 1.e3 * dslz(nzg)

! Combined sfcwater and soil energy per square meter.

      qwt = energy_per_m2(1) + soil_energy(nzg) * dslz(nzg)

! Soil heat capacity per square meter

      soilhcap = slcpd(ntext_soil(nzg)) * dslz(nzg)

! Diagnose equilibrium temperature and fractional liquid/ice water phases

      call qwtk(qwt,wt,soilhcap,tempk,fracliq)

! Diagnose new energy value for sfcwater based on qwt value.

      if (qwt < 0.) then
      
! Case of equilibrium temperature below 0 deg C.  Sfcwater fracliq = 0.

         sfcwater_fracliq(1) = 0.
         sfcwater_tempk(1) = tempk
         energy_per_m2(1) = sfcwater_mass(1) * cice * (tempk - 273.15)

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw6 ',energy_per_m2(1),sfcwater_mass(1),cice,tempk
!endif

      elseif (qwt > wt * alli) then
      
! Case of equilibrium temperature above 0 deg C.  Sfcwater fracliq = 1.

         sfcwater_fracliq(1) = 1.
         sfcwater_tempk(1) = tempk
         energy_per_m2(1) = sfcwater_mass(1) * (cliq * (tempk - 273.15) + alli)

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw7 ',energy_per_m2(1),sfcwater_mass(1),cliq,tempk,alli
!endif

      else
      
! Equilibrium temperature is 0 deg C.  Determine separate values for
! sfcwater_fracliq(1) and soil_fracliq(nzg) using constraint that the sum
! of (mass * fracliq) over both components is (wt * fracliq).

! Lower bound on sfcwater_fracliq(1): case with soil_water(nzg) all liquid

         flmin = (fracliq * wt - soil_water(nzg) * 1.e3 * dslz(nzg))  &
               / sfcwater_mass(1)         

! Upper bound on sfcwater_fracliq(1): case with soil_water(nzg) all ice

         flmax = fracliq * wt / sfcwater_mass(1)         

! New sfcwater_fracliq(1) value becomes closest value within bounds to old value.

         sfcwater_fracliq(1) = max(0.,flmin,min(1.0,flmax,sfcwater_fracliq(1)))
         sfcwater_tempk(1) = 273.15
         energy_per_m2(1) = sfcwater_mass(1) * sfcwater_fracliq(1) * alli

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw8 ',energy_per_m2(1),sfcwater_mass(1),sfcwater_fracliq(1),alli
!endif

      endif

! New energy value for soil is combined energy minus new sfcwater energy

      soil_energy(nzg) = (qwt - energy_per_m2(1)) * dslzi(nzg)

   else

! Current sfcwater layer is either not the bottom one or is thick enough
! to not require implicit thermal balance with top soil layer.

! Diagnose sfcwater temperature and liquid fraction.  Use qwtk instead of qtk
! in case sfcwater_mass(k) is too small; "dryhcap" = 100 is neglible value    
  
      call qwtk(energy_per_m2(k),sfcwater_mass(k),100.,  &
                sfcwater_tempk(k),sfcwater_fracliq(k))

   endif

! If liquid exists in current sfcwater layer, any low-density ice structure
! tends to collapse.  Increase density accordingly using simple linear relation.

   if (snden < 1.e3 * sfcwater_fracliq(k)) then
      snden = 1.e3 * sfcwater_fracliq(k)
      sfcwater_depth(k) = sfcwater_mass(k) / snden
   endif

! Assume that excess of sfcwater_fracliq over 10% is free to drain out of layer

   if (sfcwater_fracliq(k) > .10) then
      wfree = sfcwater_mass(k) * (sfcwater_fracliq(k) - .10) / .90

! Check whether this is lowest sfcwater layer

      if (k == 1) then

! This is lowest sfcwater layer.  Reduce wfree if necessary to maximum 
! amount that can percolate into top soil layer

         soilcap = 1.e3 * max (0.,dslz(nzg)  &
                 * (slmsts(ntext_soil(nzg)) - soil_water(nzg)))
         if (wfree > soilcap) wfree = soilcap

      endif

! Evaluate energy and depth transferred in wfree (which is in liquid phase)

      qwfree = wfree * (cliq * (sfcwater_tempk(k) - 273.15) + alli)
      dwfree = wfree * .001

! Check if essentially all of sfcwater_mass(k) will drain from layer

      if (wfree > .999 * sfcwater_mass(k)) then
      
! All sfcwater_mass(k) drains from layer.  Set layer quantities to zero to
! avoid truncation error.

         sfcwater_mass(k)  = 0.
         energy_per_m2(k)  = 0.
         sfcwater_depth(k) = 0.
          
      else

! Not all sfcwater_mass(k) drains from layer.  Drain mass, energy, and depth 
! of free water out of current layer

         sfcwater_mass(k)  = sfcwater_mass(k)  - wfree
         energy_per_m2(k)  = energy_per_m2(k)  - qwfree
         sfcwater_depth(k) = sfcwater_depth(k) - dwfree

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw9 ',k,energy_per_m2(k),qwfree
!endif

      endif

! Check whether this is lowest sfcwater layer

      if (k == 1) then

! This is lowest sfcwater layer.  Drained water goes to top soil layer

         soil_water(nzg) = soil_water(nzg) + 1.e-3 * wfree * dslzi(nzg)
         soil_energy(nzg) = soil_energy(nzg) + qwfree * dslzi(nzg)

      else

! This is not lowest sfcwater layer.  Drained water goes to next lower
! sfcwater layer

         sfcwater_mass(k-1)  = sfcwater_mass(k-1)  + wfree
         energy_per_m2(k-1)  = energy_per_m2(k-1)  + qwfree
         sfcwater_depth(k-1) = sfcwater_depth(k-1) + dwfree

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw10 ',k-1,energy_per_m2(k-1),qwfree
!endif

      endif

   endif

! Add remaining sfcwater mass in current layer to mass sum

   totsnow = totsnow + sfcwater_mass(k)

enddo

! Check whether any sfcwater mass remains (after possibly having all drained into soil)

if (totsnow < 1.e-9) then

! Total sfcwater mass is very close to zero.  Set sfcwater layer count to zero,
! set sfcwater arrays to exactly zero, and return

   nlev_sfcwater = 0

   sfcwater_mass(1:nzs)   = 0.
   sfcwater_energy(1:nzs) = 0.
   sfcwater_depth(1:nzs)  = 0.

   return
   
endif

! Sfcwater mass remains.  Re-distribute mass, energy, and depth among layers to
! maintain prescribed distribution of mass

! Find maximum number of layers for which thinnest layer (top and bottom) would 
! be thicker than snowmin.

maxlayers = 1

do while (maxlayers < nzs .and. maxlayers < 10 .and.  &
   thick(1,maxlayers+1) * totsnow > snowmin)

   maxlayers = maxlayers + 1
enddo

! Count up all existing snow layers that are not totally liquid

icelayers = 0

do k = 1,nlev_sfcwater
   if (sfcwater_mass(k) > 1.e-9 .and.  &
       energy_per_m2(k) < sfcwater_mass(k) * alli) then

      icelayers = icelayers + 1
   endif
enddo

! Determine new number of layers.  This number may be at most one greater than
! nlev_sfcwater and one greater than icelayers, but no greater than maxlayers.

nlev_new = min (nlev_sfcwater+1, icelayers+1, maxlayers)

! Set index for first old layer; set transfer weights for first old layer

kold = 1
wtold = 1. ! fraction of "remaining" mass in old layer

! Loop over new set of layers

do k = 1,nlev_new

! To begin current new layer, set fraction of "unfilled" mass in new layer to 1.

   wtnew = 1.

! Set mass of current new layer (already determined)

   mass_new(k) = totsnow * thick(k,nlev_new)

! Initialize energy, depth, and s/w flux of current new layer to zero

   energy_new (k) = 0.
   depth_new  (k) = 0.
   rshort_snew(k) = 0.

10 continue

! Compute "unfilled" mass in current new layer minus remaining mass in
! current old layer

   wdiff = wtnew * mass_new(k) - wtold * sfcwater_mass(kold)

! Check sign of wdiff

   if (wdiff > 0.) then

! If "unfilled" mass in current new layer exceeds remaining mass in current old
! layer, transfer all of remaining energy, depth, and s/w flux from old layer

      energy_new (k) = energy_new (k) + wtold * energy_per_m2 (kold)
      depth_new  (k) = depth_new  (k) + wtold * sfcwater_depth(kold)
      rshort_snew(k) = rshort_snew(k) + wtold * rshort_s      (kold)

! Reduce fraction of "unfilled" mass in current new layer, jump to next old
! layer, and return wtold to 1.0 to indicate no mass in old layer yet removed.

      wtnew = wtnew - wtold * sfcwater_mass(kold) / mass_new(k)
      kold = kold + 1
      wtold = 1.

! If old-layer counter does not exceed top old layer, repeat transfer operation 
! for current old layer

      if (kold <= nlev_sfcwater) go to 10

   else

! If "unfilled" mass in current new layer is less than remaining mass in current
! old layer, transfer only the portion of remaining energy, depth, and s/w flux
! from old layer that fit into new layer.

! Change in wtold

      dwtold = wtnew * mass_new(k) / sfcwater_mass(kold)

      wtold = wtold - dwtold

! Energy, depth, and s/w flux transfer to current new layer

      energy_new (k) = energy_new (k) + dwtold * energy_per_m2 (kold)
      depth_new  (k) = depth_new  (k) + dwtold * sfcwater_depth(kold)
      rshort_snew(k) = rshort_snew(k) + dwtold * rshort_s      (kold)

   endif

enddo

! Now that mass, energy, depth, and s/w flux have been transferred to new layers,
! copy back to original arrays

do k = 1,nlev_new
   sfcwater_mass(k)   = mass_new(k)
   sfcwater_energy(k) = energy_new(k) / sfcwater_mass(k)
   sfcwater_depth(k)  = depth_new(k)
   rshort_s(k)        = rshort_snew(k)

! Replace sfcwater energy limiter from earlier code versions with energy
! check and message.  If message ever gets printed, investigate reasons.

   if (sfcwater_energy(k) > 5.e5 .or. sfcwater_energy(k) < -2.e5) then
      write(io6,*) ' '
      write(io6,*) 'Sfcwater energy is outside allowable range. '
      write(io6,*) 'iwl,k,sfcwater_energy = ',iwl,k,sfcwater_energy(k)
      write(io6,*) 'p1',energy_new(k),sfcwater_mass(k),nlev_sfcwater,nlev_new
      write(io6,*) 'p2',kold, energy_per_m2(kold)
      stop 'stop sfcwater energy'
   endif

enddo

nlev_sfcwater = nlev_new

return
end subroutine sfcwater

!===============================================================================

subroutine soil(iwl, leaf_class, nlev_sfcwater, ntext_soil, ktrans,    &
                soil_tempk, soil_fracliq, soil_rfactor,                  &
                hxfergc, wxfergc, rshort_g, rlong_g, transp,             &
                soil_energy, soil_water, hxferg, wxfer, qwxfer,          &
                psiplusz, hydraul_cond, lsl, ed_transp, ed_patch         )

use leaf_coms, only: nzg, dslz, dslzi, slzt, dslzidt, dslztidt, dt_leaf,  &
                     slcons1, soilcp, slbs, slpots, slmsts, kroot

use consts_coms, only: cliq1000, alli1000, alvi
use misc_coms, only: io6

use ed_structure_defs

implicit none

integer, intent(in) :: iwl            ! current land cell number 
integer, intent(in) :: leaf_class       ! leaf class
integer, intent(in) :: nlev_sfcwater    ! # active levels of surface water
integer, intent(in) :: ntext_soil (nzg) ! soil textural class

integer, intent(in) :: ktrans           ! k index of soil layer supplying transpiration

real, intent(in) :: soil_tempk    (nzg) ! soil temperature (K)
real, intent(in) :: soil_fracliq  (nzg) ! fraction of soil water that is liquid
real, intent(in) :: soil_rfactor  (nzg) ! soil thermal resistance
real, intent(in) :: hxfergc             ! heat xfer from ground to canopy [kg/m^2]
real, intent(in) :: wxfergc             ! water xfer from ground to canopy [kg/m^2]
real, intent(in) :: rshort_g            ! s/w radiative flux abs by ground [W/m^2]
real, intent(in) :: rlong_g             ! l/w radiative flux abs by ground [W/m^2]
real, intent(in) :: transp              ! transpiration loss [kg/m^2]

real, intent(inout) :: soil_energy(nzg) ! [J/m^3]
real, intent(inout) :: soil_water (nzg) ! soil water content [vol_water/vol_tot]

real, intent(out) :: hxferg      (nzg+1) ! soil internal heat xfer (J/m^2]
real, intent(out) :: wxfer       (nzg+1) ! soil water xfer [m]
real, intent(out) :: qwxfer      (nzg+1) ! soil energy xfer from water xfer [J/m^2] 
real, intent(out) :: psiplusz    (nzg)   ! soil water potential (including grav) [m]
real, intent(out) :: hydraul_cond(nzg)   ! soil hydraulic conductivity [m/s]
type(patch), target, optional :: ed_patch
integer, intent(in) :: lsl
! Local variables

real :: half_soilair(nzg) ! half of available airspace in soil layer [m]
real :: soil_liq    (nzg) ! liq water in soil layer available for transport [m]

integer :: k     ! vertical index over soil layers
integer :: nts   ! soil textural class

real :: watermid ! soil water content midway between layers [vol_water/vol_tot]
real :: availwat ! soil water available for transpiration
real :: wg       ! for finding soil layer with max availwat
real :: wloss    ! soil water loss from transpiration [vol_water/vol_tot]
real :: qwloss   ! soil energy loss from transpiration [J/vol_tot]
real :: runoff   ! runoff loss [kg/m^2]
real :: psiplusz0 ! assumed saturated (water + grav) potential at imaginary 
                  !   layer below lowest soil layer; used for drainage flux

real, dimension(nzg) :: ed_transp

! Remove transpiration water from ktrans soil layer
! Units of wloss are [vol_water/vol_tot], of transp are [kg/m^2].

if ((.not. present(ed_patch)) .and. ktrans > 0) then

   wloss = transp * dslzi(ktrans) * 1.e-3

   soil_water(ktrans) = soil_water(ktrans) - wloss

   qwloss = wloss * (cliq1000 * (soil_tempk(ktrans) - 273.15) + alli1000)

   soil_energy(ktrans) = soil_energy(ktrans) - qwloss

elseif (present(ed_patch)) then

   do k = lsl, nzg
      wloss = ed_transp(k) * dslzi(k) * 1.e-3
      qwloss = wloss * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000)
      soil_water(k) = soil_water(k) - wloss
      soil_energy(k) = soil_energy(k) - qwloss
      ed_patch%omean_latflux =   &
           ed_patch%omean_latflux + qwloss * dslz(k) / dt_leaf
   enddo

endif

! Soil bottom, top, and internal sensible heat xfers [J/m2]

hxferg(1) = 0.
hxferg(nzg+1) = 0.

do k = 2,nzg
   hxferg(k) = dt_leaf * (soil_tempk(k-1) - soil_tempk(k))   &
             / (soil_rfactor(k-1) + soil_rfactor(k))      
enddo

! Update soil Q values [J/m3] from internal sensible heat and upward water 
! vapor (latent heat) xfers, and from longwave and shortwave fluxes at the
! top of the soil.  This excludes effects of dew/frost formation, 
! precipitation, shedding, and percolation, which were already applied 
! to the top soil layer in subroutine sfcwater.  Update top soil moisture 
! from evaporation only if sfcwater was absent.

do k = 1,nzg
   soil_energy(k) = soil_energy(k) + dslzi(k) * (hxferg(k) - hxferg(k+1))
enddo

soil_energy(nzg) = soil_energy(nzg)  &
   + dslzi(nzg) * (dt_leaf * (rshort_g + rlong_g) - hxfergc - wxfergc * alvi)

if (nlev_sfcwater == 0) then
   soil_water(nzg) = soil_water(nzg) - 1.e-3 * wxfergc * dslzi(nzg)
endif

! [12/07/04] Revisit the computation of water xfer between soil layers in
! relation to the extreme nonlinearity of hydraulic conductivity with respect
! to soil moisture.  What is the best value for hydraulic conductivity (or
! resistivity) at the interface between the two layers?  The answer is
! definitely not the average of the resistivity values of the layers
! because a very dry layer would shut down xfer of water into it.  The average 
! conductivity would, on the other hand, over-estimate water xfer between layers 
! when one is wet and the other dry.  A good compromise seems to be to average
! the fractional moisture content between the two layers and to apply this
! average value in computing hydraulic resistance for the bottom half of the
! upper layer and the top half of the lower layer.  Then, add these resistances
! to obtain the total hydraulic resistance between the two layers.

! Compute gravitational potential plus moisture potential z + psi [m],
! liquid water content [m], and half the remaining water capacity [m].

do k = 1,nzg
   nts = ntext_soil(k)

   psiplusz(k) = slpots(nts) * (slmsts(nts) / soil_water(k)) ** slbs(nts)  &
               + slzt(k) 

   soil_liq(k) = dslz(k)  &
      * min(soil_water(k) - soilcp(nts) , soil_water(k) * soil_fracliq(k))

   half_soilair(k) = (slmsts(nts) - soil_water(k)) * dslz(k) * .5
enddo

! Find amount of water transferred between soil layers (wxfer) [m]
! modulated by the liquid water fraction

do k = 2,nzg
   nts = ntext_soil(k)

   watermid = .5 * (soil_water(k) + soil_water(k-1))
   
   hydraul_cond(k) = slcons1(k,nts)  &
      * (watermid / slmsts(nts)) ** (2. * slbs(nts) + 3.)
      
   wxfer(k) = dslztidt(k) * hydraul_cond(k) * (psiplusz(k-1) - psiplusz(k)) &
      * .5 * (soil_fracliq(k) + soil_fracliq(k-1))

! Limit water transfers to prevent over-saturation and over-depletion
! Compute q transfers between soil layers (qwxfer) [J/m2]

   if (wxfer(k) > 0.) then
      wxfer(k) = min(wxfer(k),soil_liq(k-1),half_soilair(k))
   else
      wxfer(k) = - min(-wxfer(k),soil_liq(k),half_soilair(k-1))
   endif

   qwxfer(k) = wxfer(k) * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000)
enddo

wxfer(1)    = 0.
wxfer(nzg+1)  = 0.
qwxfer(1)   = 0.
qwxfer(nzg+1) = 0.

!--------------------------------------------------------------------------------
! TEMPORARY; OPTIONAL - Remove soil water from lowest layer to represent 
! runoff + infiltration.
 
! RUNOFF has units of [kg/m^2]

! Use psiplusz value "at k = 0" that represents saturated soil

! psiplusz0 = slpots(nts) + 2. * slzm(1) - slzt(1) 

! wxfer(1) = dslztidt(1) * hydraul_cond(1) * (psiplusz0 - psiplusz(1)) &
!       * soil_fracliq(1)
! qwxfer(1) = wxfer(1) * (cliq1000 * (soil_tempk(1) - 273.15) + alli1000)

!--------------------------------------------------------------------------------

! Update soil moisture (impose minimum value of soilcp) and q value.

do k = 1,nzg
   nts = ntext_soil(k)
   soil_water(k) = max(soilcp(nts),soil_water(k)  &
      - dslzi(k) * (wxfer(k+1) - wxfer(k)))
   soil_energy(k) = soil_energy(k) - dslzi(k) * (qwxfer(k+1) - qwxfer(k))
enddo

! Compute soil respiration if ED is being run

if (present(ed_patch)) then
   call soil_respiration(ed_patch)
endif

return
end subroutine soil

!===============================================================================

subroutine grndvap(iwl, nlev_sfcwater, nts, soil_water, soil_energy,    &
                   sfcwater_energy, rhos, can_shv, ground_shv, surface_ssh)

use leaf_coms,   only: nstyp, slbs, slmsts, slpots, slcpd, nzg
use consts_coms, only: grav, rvap
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iwl         ! current land cell number 
integer, intent(in) :: nlev_sfcwater ! # active levels of surface water
integer, intent(in) :: nts           ! soil textural class (local name)

real, intent(in)  :: soil_water      ! soil water content [vol_water/vol_tot]
real, intent(in)  :: soil_energy     ! [J/m^3]
real, intent(in)  :: sfcwater_energy ! [J/kg]
real, intent(in)  :: rhos            ! air density [kg/m^3]
real, intent(in)  :: can_shv         ! canopy vapor spec hum [kg_vap/kg_air]
real, intent(out) :: ground_shv      ! ground equilibrium spec hum [kg_vap/kg_air]
real, intent(out) :: surface_ssh     ! surface (saturation) spec hum [kg_vap/kg_air]

! Local parameter

real, parameter :: gorvap = grav / rvap  ! gravity divided by vapor gas constant

! Local variables

real :: slpotvn ! soil water potential [m]
real :: alpha   ! "alpha" term in Lee and Pielke (1993)
real :: beta    ! "beta" term in Lee and Pielke (1993)
real :: tempk   ! surface water temp [K]
real :: fracliq ! fraction of surface water in liquid phase
real, save, dimension(nstyp) :: sfldcap  ! soil water field capacity [vol_water/vol_tot]

data sfldcap/.135,.150,.195,.255,.240,.255,.322,.325,.310,.370,.367,.535/

real, external :: rhovsil  ! function to compute sat vapor density (over ice or liq)

! surface_ssh is the saturation mixing ratio of the top soil or snow surface
! and is used for dew formation and snow evaporation.

if (nlev_sfcwater > 0) then
   call qtk(sfcwater_energy,tempk,fracliq)
   surface_ssh = rhovsil(tempk-273.15) / rhos
else

! Without snowcover, ground_shv is the effective saturation mixing
! ratio of soil and is used for soil evaporation.  First, compute the
! "alpha" term or soil "relative humidity" and the "beta" term.

   call qwtk(soil_energy,soil_water*1.e3,slcpd(nts),tempk,fracliq)
   surface_ssh = rhovsil(tempk-273.15) / rhos

   slpotvn = slpots(nts) * (slmsts(nts) / soil_water) ** slbs(nts)
   alpha = exp(gorvap * slpotvn / tempk)
   beta = .25 * (1. - cos (min(1.,soil_water / sfldcap(nts)) * 3.14159)) ** 2
   ground_shv = surface_ssh * alpha * beta + (1. - beta) * can_shv

endif

return
end subroutine grndvap

!===============================================================================

subroutine vegndvi(iwl, leaf_class, timefac_ndvi, veg_height,         &
                   veg_ndvip, veg_ndvif, veg_ndvic,                     &
                   veg_tai, veg_lai, veg_fracarea, veg_albedo, veg_rough)

use leaf_coms, only: veg_frac, albv_brown, albv_green, sai, dead_frac,   &
                     fpar_max, veg_clump, glai_max, dfpardsr, fpar_min,  &
                     sr_max, sr_min, tai_max
use misc_coms, only: io6

implicit none

integer, intent(in) :: iwl      ! index of current land cell
integer, intent(in) :: leaf_class ! leaf class

real, intent(in)  :: timefac_ndvi ! frac of time from past to future NDVI obs
real, intent(in)  :: veg_height   ! veg height [m]
real, intent(in)  :: veg_ndvip    ! veg past ndvi (obs time)
real, intent(in)  :: veg_ndvif    ! veg future ndvi (obs time)
real, intent(out) :: veg_ndvic    ! veg current ndvi
real, intent(out) :: veg_tai      ! veg total area index
real, intent(out) :: veg_lai      ! veg leaf area index
real, intent(out) :: veg_fracarea ! veg fractional area
real, intent(out) :: veg_albedo   ! veg albedo
real, intent(out) :: veg_rough    ! veg roughness height [m]

! Local parameters

real, parameter :: bz         = .91       ! for computing veg roughness height
real, parameter :: hz         = .0075     ! for computing veg roughness height
real, parameter :: extinc_veg = .5        ! veg extinction coefficient for rad flux
real, parameter :: fpcon      = -.3338082 ! for computing veg leaf area index

! Local variables

real :: sr         ! simple ratio
real :: fpar       ! fraction of photosynthetically-active radiation
real :: dead_lai   ! dead-matter leaf area index
real :: green_frac ! lai fraction of tai

! Update vegetation TAI, LAI, fractional area, albedo, and roughness

! Compute LAI, vegetation roughness, albedo, vegfrac from time-dependent NDVI

if (tai_max(leaf_class) < .1) then

   veg_lai = 0.
   veg_tai = 0.
   veg_rough = 0.
   veg_albedo = 0.
   veg_fracarea = 0.

else

! Time-interpolate ndvi to get current value veg_ndvic5 for this land cell

   veg_ndvic = veg_ndvip + (veg_ndvif - veg_ndvip) * timefac_ndvi
      
! Limit ndvi to prevent values > .99 to prevent division by zero.

   if (veg_ndvic > .99) veg_ndvic = .99

! Compute "simple ratio" and limit between sr_min and sr_max(leaf_class).

   sr = (1. + veg_ndvic) / (1. - veg_ndvic)

   if (sr < sr_min) sr = sr_min
   if (sr > sr_max(leaf_class)) sr = sr_max(leaf_class)

! Compute fpar

   fpar = fpar_min + (sr - sr_min) * dfpardsr(leaf_class)

! Compute green leaf area index (veg_lai), dead leaf area index (dead_lai),
! total area index (tai), and green fraction   

   veg_lai = glai_max(leaf_class) * (veg_clump(leaf_class) * fpar / fpar_max  &
      + (1. - veg_clump(leaf_class)) * alog(1. - fpar) * fpcon)

   dead_lai = (glai_max(leaf_class) - veg_lai) * dead_frac(leaf_class)

   veg_tai = sai(leaf_class) + veg_lai + dead_lai
   green_frac = veg_lai / veg_tai

! Compute vegetation roughness height, albedo, and fractional area

   veg_rough = veg_height * (1. - bz * exp(-hz * veg_tai))
   veg_albedo = albv_green(leaf_class) * green_frac  &
              + albv_brown(leaf_class) * (1. - green_frac)
   veg_fracarea = veg_frac(leaf_class) * (1. - exp(-extinc_veg * veg_tai))

endif         
return
end subroutine vegndvi

! Remaining issues:
! 
! 1. Relationship between clumping, V, vegfrac
! 2. Impact of V on radiation
! 3. Build lookup tables, especially for things with exponentials?

!===============================================================================

subroutine sfcrad_land(iwl, leaf_class, ntext_soil, nlev_sfcwater,      &
                       sfcwater_energy, sfcwater_depth,                 &
                       soil_energy, soil_water,                         &
                       veg_temp, veg_fracarea, veg_height,              &
                       veg_albedo, rshort, rlong,                       &
                       rshort_s, rshort_g, rshort_v,                    &
                       rlong_g, rlong_s, rlong_v,                       &
                       rlongup, rlong_albedo, albedo_beam, snowfac, vf, &
                       cosz, xewl, yewl, zewl                           )

use leaf_coms,   only: nzs, slcpd, slmsts, emisv, emisg
use consts_coms, only: stefan, eradi
use misc_coms,   only: io6
use mem_radiate, only: sunx, suny, sunz

implicit none

integer, intent(in) :: iwl         ! index of current land cell
integer, intent(in) :: leaf_class    ! leaf class
integer, intent(in) :: ntext_soil    ! soil textural class
integer, intent(in) :: nlev_sfcwater ! # of active surface water layers

real, intent(in) :: sfcwater_energy(nzs) ! surface water internal energy [J/kg]
real, intent(in) :: sfcwater_depth(nzs)  ! surface water depth [m]
real, intent(in) :: soil_energy  ! soil internal energy [J/m^3]
real, intent(in) :: soil_water   ! soil water content [vol_water/vol_tot]
real, intent(in) :: veg_temp     ! veg temp [K]
real, intent(in) :: veg_fracarea ! veg fractional area coverage
real, intent(in) :: veg_height   ! veg height [m]
real, intent(in) :: veg_albedo   ! veg albedo
real, intent(in) :: rshort       ! downward surface incident s/w rad flux [W/m^2]
real, intent(in) :: rlong        ! downward surface incident l/w rad flux [W/m^2]
real, intent(in) :: xewl         ! land cell earth x coordinate
real, intent(in) :: yewl         ! land cell earth y coordinate
real, intent(in) :: zewl         ! land cell earth z coordinate

real, intent(out) :: rshort_s(nzs) ! s/w net rad flux to sfc water [W/m^2]
real, intent(out) :: rshort_g ! s/w net rad flux to soil [W/m^2]
real, intent(out) :: rshort_v ! s/w net rad flux to veg [W/m^2]
real, intent(out) :: rlong_g  ! l/w net rad flux to soil [W/m^2]
real, intent(out) :: rlong_s  ! l/w net rad flux to sfc water [W/m^2]
real, intent(out) :: rlong_v  ! l/w net rad flux to veg [W/m^2]
real, intent(out) :: rlongup  ! upward sfc l/w rad flux [W/m^2]
real, intent(out) :: rlong_albedo  ! albedo for upward sfc l/w
real, intent(out) :: albedo_beam   ! land cell albedo (beam)
real, intent(out) :: snowfac  ! fractional veg burial by snowcover
real, intent(out) :: vf       ! fractional coverage of non-buried part of veg
real, intent(out) :: cosz     ! solar zenith angle for land cell

! Local variables

integer :: k             ! vertical index over surface water layers

real :: soil_tempk       ! soil temp [K]
real :: soil_fracliq     ! fraction of soil water in liquid phase
real :: sfcwater_tempk   ! surface water temp [K]
real :: sfcwater_fracliq ! fraction of surface water in liquid phase
real :: rlonga_v         ! longwave radiative flux from atm  to veg  [W/m^2]
real :: rlonga_s         ! longwave radiative flux from atm  to snow [W/m^2]
real :: rlonga_g         ! longwave radiative flux from atm  to soil [W/m^2]
real :: rlongv_a         ! longwave radiative flux from veg  to atm  [W/m^2]
real :: rlongv_s         ! longwave radiative flux from veg  to snow [W/m^2]
real :: rlongv_g         ! longwave radiative flux from veg  to soil [W/m^2]
real :: rlongs_a         ! longwave radiative flux from snow to atm  [W/m^2]
real :: rlongs_v         ! longwave radiative flux from snow to veg  [W/m^2]
real :: rlongg_a         ! longwave radiative flux from soil to atm  [W/m^2]
real :: rlongg_v         ! longwave radiative flux from soil to veg  [W/m^2]
real :: fracabs (nzs)    ! fraction of rshort that is absorbed by snowcover
real :: vfc       ! 1 - vf
real :: fcpct     ! soil water fraction
real :: alg       ! soil albedo
real :: als       ! snowcover albedo (needs better formula based on age of snow)
real :: alv       ! veg albedo
real :: rad       ! fraction s/w rad absorbed into sfcwater + soil
real :: fractrans ! fraction of s/w rad flux transmitted through snowcover layer
real :: absg      ! fraction of rshort that is absorbed by ground
real :: algs      ! albedo from snow plus ground
real :: emv       ! veg emissivity
real :: emg       ! soil emissivity
real :: ems       ! surface water emissivity
real :: glong     ! soil l/w rad emission [W/m^2]
real :: slong     ! sfc water l/w rad emission [W/m^2]
real :: vlong     ! veg l/w rad emission [W/m^2]

! This routine is called twice by the radiation parameterization, performing
! exactly the same operations on both calls.  All that is used from the first
! call are net surface albedo and upward longwave radiative flux for each land
! cell.  The radiation parameterization carries out atmospheric radiative 
! transfer computations following this first call using the albedo and upward
! longwave flux.  The second call to this subroutine is made after the 
! atmospheric radiative fluxes are computed.  This call provides net radiative 
! fluxes to vegetation, snowcover, and soil in each land cell, plus functions 
! of snowcover, all of which are required in leaf3.

! Compute solar zenith angle for land cells

cosz = (xewl * sunx + yewl * suny + zewl * sunz) * eradi

if (nlev_sfcwater == 0) then

! Case with no surface water

! Shortwave radiation calculations

   snowfac = 0.
   alv = veg_albedo

   fcpct = soil_water / slmsts(ntext_soil)  ! soil water fraction

   if (fcpct > .5) then
      alg = .14                ! ground albedo
   else
      alg = .31 - .34 * fcpct  ! ground albedo
   endif

   vf = veg_fracarea
   vfc = 1. - vf

   albedo_beam = vf * alv + vfc * vfc * alg

   rshort_g = rshort * vfc * (1. - alg)
   rshort_s(1:nzs) = 0.
   rshort_v = rshort * vf * (1. - alv + vfc * alg)
!  rshort_a = rshort * albedo_beam

! Longwave radiation calculations

! Diagnose soil temperature and liquid fraction

   call qwtk(soil_energy,soil_water*1.e3,  &
             slcpd(ntext_soil),soil_tempk,soil_fracliq)

   emv    = emisv(leaf_class)
   emg    = emisg(ntext_soil)
   glong  = emg * stefan * soil_tempk ** 4
   vlong  = emv * stefan * veg_temp ** 4

   rlonga_v = rlong * vf * (emv + vfc * (1. - emg))
   rlonga_g = rlong * vfc * emg
   rlongv_g = vlong * vf * emg
   rlongv_a = vlong * vf * (2. - emg - vf + emg * vf)
   rlongg_v = glong * vf * emv
   rlongg_a = glong * vfc

   rlong_albedo = (vf * (1. - emv) + vfc * vfc * (1. - emg))

   rlongup = rlongv_a + rlongg_a

   rlong_g = rlonga_g - rlongg_a + rlongv_g - rlongg_v
   rlong_s = 0.
   rlong_v = rlonga_v - rlongv_a + rlongg_v - rlongv_g

else

! Case with surface water

! Diagnose surface water temperature and liquid fraction

   call qtk(sfcwater_energy(nlev_sfcwater),sfcwater_tempk,sfcwater_fracliq)

! Shortwave radiation calculations

   alv = veg_albedo

! Sfcwater albedo ALS ranges from wet-soil value .14 for all-liquid
! to .5 for all-ice

   als = .5 - .36 * sfcwater_fracliq
   rad = 1. - als    ! fraction shortwave absorbed into sfcwater + soil

   snowfac = 0.
   do k = nlev_sfcwater,1,-1

      snowfac = snowfac + sfcwater_depth(k)

! fractrans is fraction of shortwave entering each sfcwater layer that
! gets transmitted through that layer

      fractrans = exp(-20. * sfcwater_depth(k))

! fracabs(k) is fraction of total incident shortwave (at top of top sfcwater
! layer) that is absorbed in each sfcwater layer

      fracabs(k) = rad * (1. - fractrans)

! rad is fraction of total incident shortwave (at top of top sfcwater layer)
! that remains at bottom of current sfcwater layer

      rad = rad * fractrans
   enddo

   snowfac = min(.99, snowfac / max(.001,veg_height))

   vf = veg_fracarea * (1. - snowfac)
   vfc = 1. - vf

   fcpct = soil_water / slmsts(ntext_soil)
   if (fcpct > .5) then
      alg = .14
   else
      alg = .31 - .34 * fcpct
   endif

   absg = (1. - alg) * rad
   algs = 1. - absg
   do k = nlev_sfcwater,1,-1
      algs = algs - fracabs(k)
      rshort_s(k) = rshort * vfc * fracabs(k)
   enddo

   albedo_beam = vf * alv + vfc * vfc * algs

   rshort_g = rshort * vfc * absg
   rshort_s(nlev_sfcwater+1:nzs) = 0.
   rshort_v = rshort * vf * (1. - alv + vfc * algs)
!  rshort_a = rshort * albedo_beam

! Longwave radiation calculations

   emv   = emisv(leaf_class)
   ems   = 1.0
   slong = ems * stefan * sfcwater_tempk ** 4
   vlong = emv * stefan * veg_temp ** 4

   rlonga_v = rlong * vf * (emv + vfc * (1. - ems))
   rlonga_s = rlong * vfc * ems
   rlongv_s = vlong * vf * ems
   rlongv_a = vlong * vf * (2. - ems - vf + ems * vf)
   rlongs_v = slong * vf * emv
   rlongs_a = slong * vfc

   rlong_albedo = (vf * (1. - emv) + vfc * vfc * (1. - ems))

   rlongup = rlongv_a + rlongs_a

   rlong_g = 0.
   rlong_s = rlonga_s - rlongs_a + rlongv_s - rlongs_v
   rlong_v = rlonga_v - rlongv_a + rlongs_v - rlongv_s

endif

return
end subroutine sfcrad_land

!===============================================================================

subroutine remove_runoff(iwl, ksn, sfcwater_fracliq, sfcwater_mass,   &
     sfcwater_tempk, sfcwater_energy, sfcwater_depth, runoff, qrunoff)

  use leaf_coms, only: nzs, dt_leaf
  use ed_options, only: runoff_time
  use consts_coms, only: alli, cliq
  use misc_coms, only: io6

  implicit none

  integer, intent(in) :: iwl
  integer, intent(in) :: ksn

  real, intent(in)    :: sfcwater_fracliq(nzs)
  real, intent(in)    :: sfcwater_tempk  (nzs)
  real, intent(inout) :: sfcwater_mass   (nzs)
  real, intent(inout) :: sfcwater_energy (nzs)
  real, intent(inout) :: sfcwater_depth  (nzs)
  real, intent(out) :: runoff
  real, intent(out) :: qrunoff

  ! Get rid of runoff

  runoff = 0.0
  qrunoff = 0.0

  if(ksn >= 1)then
     if(sfcwater_fracliq(ksn) > 0.1)then
        call qtk(sfcwater_energy(ksn), sfcwater_tempk(ksn),   &
             sfcwater_fracliq(ksn))
        runoff = sfcwater_mass(ksn) * (sfcwater_fracliq(ksn) - 0.1) /   &
             (0.9 * runoff_time) * dt_leaf
        qrunoff = runoff * (cliq * (sfcwater_tempk(ksn) - 273.15) +  alli)
        
        sfcwater_energy(ksn) = (sfcwater_energy(ksn) *   &
             sfcwater_mass(ksn) - qrunoff ) / (sfcwater_mass(ksn) - runoff)
        sfcwater_mass(ksn) = sfcwater_mass(ksn) - runoff
        sfcwater_depth(ksn) = sfcwater_depth(ksn) - 0.001 * runoff
     endif
  endif
  
  return
end subroutine remove_runoff

