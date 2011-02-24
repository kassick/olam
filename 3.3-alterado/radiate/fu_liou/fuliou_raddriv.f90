subroutine fuliou_raddriv(iw,ka,nrad,koff)

  ! Modify CO2 concentration in subroutine set_default_options_fu.
  ! Not clear how to handle hydrometeors.

  use mem_grid, only: mza, zm, zt, glatw, glonw
  use mem_basic, only: theta, press, rho, sh_v
  use consts_coms, only: p00i, rocp, stefan, solar, cp, t00
  use mem_micro, only: sh_c, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, q6, q7
  use mem_radiate, only: rlongup, albedt, cosz, solfac, rshort, rshort_top,  &
       rshortup_top, fthrd, aero_opt_depth, rlong_albedo, rlong,   &
       rlongup_top, fthrd_lw
  use micro_coms, only: jhabtab, ncat, jnmb, emb2, parm, emb0, emb1, pwmasi, &
       dnfac, gnu
  use misc_coms, only: ilwrtyp, iswrtyp
  
  use aero_coms, only: nac, n_atau, itps, a_wli, a_taus

  implicit none

  interface
     subroutine fuliou(nrad, opress, otempk, osh_v, ozone, skin_temp, &
          alb, zm, solar1, cosz, emis, cloud_water_content, &
          cloud_effective_radius, cloud_water_phase, ztl, &
          flxus, flxds, flxul, flxdl, iw,  &
          nac, n_atau, itps, a_wli, a_taus, aprof)

       integer, intent(in) :: nrad
       real, dimension(nrad), intent(in) :: opress
       real, dimension(nrad), intent(in) :: otempk
       real, dimension(nrad), intent(in) :: osh_v
       real, dimension(nrad), intent(in) :: ozone
       real, dimension(nrad), intent(in) :: cloud_water_content
       real, dimension(nrad), intent(in) :: cloud_effective_radius
       integer, dimension(nrad), intent(in) :: cloud_water_phase
       real, dimension(nrad), intent(in) :: ztl
       real, dimension(nrad), intent(inout) :: flxus
       real, dimension(nrad), intent(inout) :: flxds
       real, dimension(nrad), intent(inout) :: flxul
       real, dimension(nrad), intent(inout) :: flxdl
       real, intent(in) :: skin_temp
       real, intent(in) :: alb
       real, intent(in) :: zm
       real, intent(in) :: solar1
       real, intent(in) :: cosz
       real, intent(in) :: emis
       integer, intent(in) :: nac
       integer, intent(in) :: n_atau
       integer, intent(in) :: iw
       integer, dimension(nac), optional, intent(in) :: itps
       real, dimension(n_atau), optional, intent(in) :: a_wli
       real, dimension(n_atau,nac), optional, intent(in) :: a_taus
       real, dimension(nrad,nac), optional, intent(in) :: aprof

     end subroutine fuliou
  end interface

  integer, intent(in) :: iw
  integer, intent(in) :: ka
  integer, intent(in) :: nrad
  integer, intent(in) :: koff

  integer :: k
  integer :: kindex
  real, dimension(nrad) :: tl ! K
  real, dimension(nrad) :: rl ! kg/m3
  real, dimension(nrad) :: pl ! mb
  real, dimension(nrad) :: shvl ! kg/kg
  real, dimension(nrad) :: dl ! kg/m3
  real, dimension(nrad) :: zml ! m
  real, dimension(nrad) :: ztl ! m
  real, dimension(nrad) :: o3l ! kg/kg
  real, dimension(nrad) :: dzl ! m
  real, dimension(nrad) :: flxds
  real, dimension(nrad) :: flxus
  real, dimension(nrad) :: flxdl
  real, dimension(nrad) :: flxul
  real, dimension(nrad) :: flxds_inter
  real, dimension(nrad) :: flxus_inter
  real, dimension(nrad) :: flxdl_inter
  real, dimension(nrad) :: flxul_inter
  real, dimension(nrad) :: net_flux
  real, dimension(nrad) :: net_flux_lw
  real :: emis
  real :: skin_tempk ! K
  real, dimension(nrad) :: cloud_effective_radius ! um
  real, dimension(nrad) :: cloud_water_content ! kg/m3
  integer, dimension(nrad) :: cloud_water_phase ! 1 - liquid, 2 - frozen
  integer :: mcat
  real, dimension(mza,ncat) :: rx
  real, dimension(mza,ncat) :: cx
  real, dimension(mza,ncat) :: emb
  integer, dimension(mza,ncat) :: jhcat
  real :: tairc
  real :: rhovslair, rho_vap
  real, external :: rhovsl
  real, dimension(nrad) :: relhum
  integer :: ns
  integer :: nt
  integer :: icat
  integer :: ihcat
  real :: parmi
  real :: tcoal
  real :: fracliq6
  real :: fracliq7
  integer, save :: kradcat(15) = (/1,2,3,6,10,12,13,5,5,3,4,8,8,6,7/)
  real :: cloud_liquid_water
  real :: cloud_frozen_water
  real :: dn
! Use the following dn limiters (rather than those in microphysics) for 
! consistency with fit coefficients
  real, save :: dnmin(7) = (/   1.,   10.,   1.,  125.,   10.,   10.,   10./)
  real, save :: dnmax(7) = (/1000.,10000., 125.,10000.,10000.,10000.,10000./)
  integer :: krc
  real, dimension(nrad,6) :: waso
  real, dimension(nrad,6) :: soot
  real, dimension(nrad,6,8) :: dust
  real, dimension(nrad,6) :: sea_salt
  real, dimension(nrad, max(nac,1)) :: aprof
  integer :: iac
  real, dimension(nrad) :: dl_mclat
  real, dimension(nrad) :: pl_mclat
  real, dimension(nrad) :: rl_mclat
  real, dimension(nrad) :: tl_mclat
  real, dimension(nrad) :: o3l_mclat
  real, dimension(nrad) :: zml_mclat
  real, dimension(nrad) :: ztl_mclat
  real, dimension(nrad) :: dzl_mclat
  real :: tairk(mza) ! air temperature [K]
  real :: rhov (mza) ! vapor density [kg_vap/m^3]

  !============================================================

  ! Return if cosz is small and not doing long wave.
  if(cosz(iw) < 0.03 .and. ilwrtyp /= 4)then
     rshort(iw) = 0.0
     rshort_top(iw) = 0.0
     rshortup_top(iw) = 0.0
     return
  endif

  ! Copy surface and vertical-column values from model to radiation memory 
  ! space.  In this loop, kindex ranges from nrad:2.

  do k = ka,mza-1
     kindex = nrad - k + koff + 1
     tairk(k) = theta(k,iw) * (press(k,iw) * p00i) ** rocp
     tl(kindex) = tairk(k)
     rhov(k) = max(0.,sh_v(k,iw)) * rho(k,iw)
     pl(kindex) = 0.01 * press(k,iw)
     shvl(kindex) = sh_v(k,iw)
     dl(kindex) = rho(k,iw)
     rl(kindex) = rhov(k)
     zml(kindex) = zm(k)
     ztl(kindex) = zt(k)
  enddo

  ! Set surface values
  zml(nrad) = zm(ka-1)
  ztl(nrad) = zt(ka-1)
  pl(nrad) = pl(nrad-1) + (zml(nrad) - ztl(nrad-2)) /   &
       (ztl(nrad-1) - ztl(nrad-2)) * (pl(nrad-1) - pl(nrad-2))
  emis = 1.0 - rlong_albedo(iw)
  skin_tempk = sqrt(sqrt(rlongup(iw) / (emis * stefan)))

  tl(nrad) = skin_tempk
  dl(nrad) = dl(nrad-1)
  shvl(nrad) = shvl(nrad-1)
  rl(nrad) = rl(nrad-1)

  dl_mclat(1:nrad) = dl(nrad:1:-1)
  pl_mclat(1:nrad) = 100.0 * pl(nrad:1:-1)
  rl_mclat(1:nrad) = rl(nrad:1:-1)
  tl_mclat(1:nrad) = tl(nrad:1:-1)
  ztl_mclat(1:nrad) = ztl(nrad:1:-1)
  zml_mclat(1:nrad) = zml(nrad:1:-1)
  o3l_mclat(1:nrad) = o3l(nrad:1:-1)
  dzl_mclat(1:nrad) = dzl(nrad:1:-1)

  ! Get the ozone
  call rad_mclat(iw,nrad,koff,glatw(iw),dl_mclat,pl_mclat,  &
       rl_mclat,tl_mclat,o3l_mclat,zml_mclat,  &
       ztl_mclat,dzl_mclat)

  dl(1:nrad) = dl_mclat(nrad:1:-1)
  pl(1:nrad) = 0.01 * pl_mclat(nrad:1:-1)
  rl(1:nrad) = rl_mclat(nrad:1:-1)
  tl(1:nrad) = tl_mclat(nrad:1:-1)
  ztl(1:nrad) = ztl_mclat(nrad:1:-1)
  zml(1:nrad) = zml_mclat(nrad:1:-1)
  o3l(1:nrad) = o3l_mclat(nrad:1:-1)
  dzl(1:nrad) = dzl_mclat(nrad:1:-1)
  o3l = o3l / dl

  cloud_effective_radius(1:nrad) = 0.0
  cloud_water_phase(1:nrad) = 0
  cloud_water_content(1:nrad) = 0.0
  
  ! Prep the hydrometeors -- this is basically from subroutine cloud_opt()
  call cloudprep_rad(iw,ka,mcat,jhcat,tairk,rhov,rx,cx,emb)

  do k = ka, mza-1
     kindex = nrad - k + koff + 1
     tairc = tl(kindex) - t00
     rhovslair = rhovsl(tairc)
     rho_vap = sh_v(k,iw)*rho(k,iw)
     relhum(kindex) = min(1.0, max(0.0, rho_vap/rhovslair))
     ns = max(1,nint(100. * relhum(kindex)))
     nt = max(1,min(31,-nint(tairc)))
  enddo

  tairc = tl(nrad) - 273.15
  rhovslair = rhovsl(tairc)
  rho_vap = sh_v(ka,iw)*rho(ka,iw)
  relhum(nrad) = min(1.0, max(0.0, rho_vap/rhovslair))

  ! Loop over all hydrometeor categories
  do icat = 1,ncat
     ! Evaluate hydrometeor mean mass emb   
  
     ! This section of code was developed from subroutine enemb in 
     ! omic_misc.f90 by removing parts that are not needed for radiation 
     ! calculations.  Arrays rx and cx are the bulk mass and number of 
     ! hydrometeors PER KG OF AIR, NOT PER M^3 AS IN THE ORIGINAL SUBROUTINE 
     ! MIC_COPY.
  
     if (jnmb(icat) == 2) then
  
        do k = ka,mza-1
           ihcat = jhcat(k,icat)
           emb(k,icat) = emb2(ihcat)
           cx(k,icat) = rx(k,icat) / emb(k,icat)
        enddo
  
     elseif (jnmb(icat) == 4) then
  
        parmi = 1. / parm(icat)
        do k = ka,mza-1
           emb(k,icat) = max(emb0(icat),min(emb1(icat),rx(k,icat) * parmi))
           cx(k,icat) = rx(k,icat) / emb(k,icat)
        enddo
  
     elseif (jnmb(icat) >= 5) then
  
        do k = ka,mza-1
           emb(k,icat) = max(emb0(icat),min(emb1(icat),rx(k,icat)  &
                / max(1.e-9,cx(k,icat))))
           cx(k,icat) = rx(k,icat) / emb(k,icat)
        enddo
  
     endif
  
  enddo

  ! Now, compute cloud info
  do k = ka,mza-1

     kindex = nrad - k + koff + 1

     ! Get liquid water fractions for graupel and hail
     call qtc(q6(k,iw), tcoal, fracliq6)
     call qtc(q7(k,iw), tcoal, fracliq7)


!     cloud_liquid_water = sh_c(k,iw)
!     cloud_frozen_water = 0.0
     cloud_liquid_water = sh_c(k,iw) + sh_r(k,iw) + sh_g(k,iw) *   &
          fracliq6 + sh_h(k,iw) * fracliq7

     cloud_frozen_water = sh_p(k,iw) + sh_s(k,iw) + sh_a(k,iw) +   &
          sh_g(k,iw) * (1.0 - fracliq6) + sh_h(k,iw) *   &
          (1.0 - fracliq7)

     if(cloud_liquid_water >= cloud_frozen_water)then
        cloud_water_phase(kindex) = 1
!        cloud_effective_radius(kindex) = 17.5
     else
        cloud_water_phase(kindex) = 2
!        cloud_effective_radius(kindex) = 100.0
     endif
     
     cloud_effective_radius(kindex) = 0.0
     do icat = 1,mcat
        if (jnmb(icat) > 0) then
           if (cx(k,icat) > 1.e-9) then
              ihcat = jhcat(k,icat)
              krc = kradcat(ihcat)
              dn = gnu(icat) * 1.e6 * dnfac(ihcat) * emb(k,icat) ** pwmasi(ihcat)
              dn = max(dnmin(icat),min(dnmax(icat),dn))  ! dn units are microns
              cloud_effective_radius(kindex) =   &
                   cloud_effective_radius(kindex) + dn * cx(k,icat)
           endif
        endif
     enddo

     cloud_effective_radius(kindex) = cloud_effective_radius(kindex) /   &
          max(1.0e-30,sum(cx(k,1:mcat)))

     if(cloud_water_phase(kindex) == 1)then
        cloud_effective_radius(kindex) = max(5.0,min(30.0,  &
             cloud_effective_radius(kindex)))
     else
        cloud_effective_radius(kindex) = max(20.0,min(180.0,  &
             cloud_effective_radius(kindex)))
     endif

     cloud_water_content(kindex) = (cloud_liquid_water +   &
          cloud_frozen_water) * rho(k,iw) * 1000.0

  enddo

  ! Take care of aerosols
  if(nac > 0)then
     ! Get the aerosols
     call getaer(nrad, iw, relhum, tl, pl, waso, soot, sea_salt)
     ! Get the dust
     call getdst(iw, nrad, pl, dust)

     ! Pack the aerosols
     aprof(:,:) = 0.0
     do iac = 1, nac
        if(itps(iac) == 10)then
           a_taus(1,iac) = sum(waso(1:nrad,6))
           if(a_taus(1,iac) > 0.0)aprof(1:nrad,iac) = waso(1:nrad,6) /   &
                a_taus(1,iac)
        elseif(itps(iac) == 11)then
           a_taus(1,iac) = sum(soot(1:nrad,6))
           if(a_taus(1,iac) > 0.0)aprof(1:nrad,iac) = soot(1:nrad,6) /   &
                a_taus(1,iac)
        elseif(itps(iac) == 12 .or. itps(iac) == 13)then
           a_taus(1,iac) = sum(sea_salt(1:nrad,6))
           if(a_taus(1,iac) > 0.0)aprof(1:nrad,iac) = sea_salt(1:nrad,6) /  &
                a_taus(1,iac)
        elseif(itps(iac) == 4)then
           a_taus(1,iac) = sum(dust(1:nrad,6,1:3))
           if(a_taus(1,iac) > 0.0)aprof(1:nrad,iac) =   &
                sum(dust(1:nrad,6,1:3)) / a_taus(1,iac)
        elseif(itps(iac) == 5)then
           a_taus(1,iac) = sum(dust(1:nrad,6,4:5))
           if(a_taus(1,iac) > 0.0)aprof(1:nrad,iac) =   &
                sum(dust(1:nrad,6,4:5)) / a_taus(1,iac)
        elseif(itps(iac) == 6)then
           a_taus(1,iac) = sum(dust(1:nrad,6,6))
           if(a_taus(1,iac) > 0.0)aprof(1:nrad,iac) = dust(1:nrad,6,6) /  &
                a_taus(1,iac)
        elseif(itps(iac) == 7)then
           a_taus(1,iac) = sum(dust(1:nrad,6,7))
           if(a_taus(1,iac) > 0.0)aprof(1:nrad,iac) = dust(1:nrad,6,7) /  &
                a_taus(1,iac)
        elseif(itps(iac) == 8)then
           a_taus(1,iac) = sum(dust(1:nrad,6,8))
           if(a_taus(1,iac) > 0.0)aprof(1:nrad,iac) = dust(1:nrad,6,8) /  &
                a_taus(1,iac)
        else
           print*,'bad aerosol constituent; fuliou_raddriv.f90'
           stop
        endif
     enddo

     aero_opt_depth(iw) = sum(a_taus(1,1:nac))
  endif

  flxus(:) = 0.0
  flxds(:) = 0.0
  flxul(:) = 0.0
  flxdl(:) = 0.0

  if(nac > 0)then
     call fuliou(nrad, pl, tl, shvl, o3l, skin_tempk, albedt(iw),   &
          zm(ka-1), solar * solfac, cosz(iw), emis,   &
          cloud_water_content,  &
          cloud_effective_radius, cloud_water_phase, ztl,   &
          flxus, flxds, flxul, flxdl, iw,  &
          nac, n_atau, itps, a_wli, a_taus, aprof)
  else
     call fuliou(nrad, pl, tl, shvl, o3l, skin_tempk, albedt(iw),   &
          zm(ka-1), solar * solfac, cosz(iw), emis,   &
          cloud_water_content,  &
          cloud_effective_radius, cloud_water_phase, ztl,   &
          flxus, flxds, flxul, flxdl, iw,  &
          nac, n_atau)
  endif

  ! Obtain fluxes to cell interfaces
  ! First, the surface grid cell is all set
  if(iswrtyp == 4)then
     flxus_inter(nrad) = flxus(nrad)
     flxds_inter(nrad) = flxds(nrad)
     rshort(iw) = max(0.0,flxds_inter(nrad))

     ! Then, extrapolate to the top of the atmosphere
     flxus_inter(1) = (flxus(1) * (2.0 * zm(mza) - zm(mza-1) - zm(mza-2)) +   &
          flxus(2) * (zm(mza-1)-zm(mza))) / (zm(mza) - zm(mza-2))
     flxds_inter(1) = (flxds(1) * (2.0 * zm(mza) - zm(mza-1) - zm(mza-2)) +   &
          flxds(2) * (zm(mza-1)-zm(mza))) / (zm(mza) - zm(mza-2))
     rshortup_top(iw) = flxus_inter(1)
     rshort_top(iw) = flxds_inter(1)

     ! Finally, interpolate between interior cells
     do k = 2, nrad-1
        kindex = mza - k + 1
        flxus_inter(k) = (flxus(k-1) * (zm(k-1) - zm(k)) +   &
             flxus(k) * (zm(k)-zm(k+1))) / (zm(k-1) - zm(k+1))
        flxds_inter(k) = (flxds(k-1) * (zm(k-1) - zm(k)) +   &
             flxds(k) * (zm(k)-zm(k+1))) / (zm(k-1) - zm(k+1))
     enddo
  endif

  if(ilwrtyp == 4)then

     flxul_inter(nrad) = flxul(nrad)
     flxdl_inter(nrad) = flxdl(nrad)
     rlong(iw) = flxdl_inter(nrad)

     flxul_inter(1) = (flxul(1) * (2.0 * zm(mza) - zm(mza-1) - zm(mza-2)) +   &
          flxul(2) * (zm(mza-1)-zm(mza))) / (zm(mza) - zm(mza-2))
     flxdl_inter(1) = (flxdl(1) * (2.0 * zm(mza) - zm(mza-1) - zm(mza-2)) +   &
          flxdl(2) * (zm(mza-1)-zm(mza))) / (zm(mza) - zm(mza-2))
     rlongup_top(iw) = flxul_inter(1)
     
     do k = 2, nrad-1
        kindex = mza - k + 1
        flxul_inter(k) = (flxul(k-1) * (zm(k-1) - zm(k)) +   &
             flxul(k) * (zm(k)-zm(k+1))) / (zm(k-1) - zm(k+1))
        flxdl_inter(k) = (flxdl(k-1) * (zm(k-1) - zm(k)) +   &
             flxdl(k) * (zm(k)-zm(k+1))) / (zm(k-1) - zm(k+1))
     enddo

  endif


  net_flux(:) = 0.0
  net_flux_lw(:) = 0.0
  if(iswrtyp == 4)then
     do k = 1, nrad-1
        net_flux(k) = flxus_inter(k+1) - flxus_inter(k) + flxds_inter(k) -   &
             flxds_inter(k+1)
     enddo
  endif

  if(ilwrtyp == 4)then
     do k = 1, nrad-1
        net_flux(k) = net_flux(k) + flxul_inter(k+1) - flxul_inter(k) +   &
             flxdl_inter(k) - flxdl_inter(k+1)
        net_flux_lw(k) = flxul_inter(k+1) - flxul_inter(k) +   &
             flxdl_inter(k) - flxdl_inter(k+1)
     enddo
  endif

  do k = ka, mza-1
     kindex = nrad - 1 - (k-ka)
     fthrd(k,iw) = fthrd(k,iw) + net_flux(kindex) /   &
          (dl(kindex) * dzl(kindex) * cp)
     fthrd_lw(k,iw) = fthrd_lw(k,iw) + net_flux_lw(kindex) /   &
          (dl(kindex) * dzl(kindex) * cp)
  enddo

  return
end subroutine fuliou_raddriv

!==================================================================

subroutine fuliou(nrad, opress, otempk, osh_v, ozone, skin_temp, alb, zm,   &
     solar1, cosz, emis, cloud_water_content, cloud_effective_radius,  &
     cloud_water_phase, ztl, &
     flxus, flxds, flxul, flxdl, iw,  &
     nac, n_atau,   &
     itps, a_wli, a_taus, aprof)

  use extras, only: getatmosphere
  use fulioumulti, only: fi, rad_multi_fu
  use fuinput, only: set_default_options_fu
  use fuoutput, only: fsfc
  use vla_fu, only: prepare_model_profile_fu, vla_interface_fu
  USE GENERATE_FULIOU_LEVELS ,only : gflq, generate_level_scheme
  USE EXTRAS       ,only : aer_scale_hgt
  USE CALIPSO_OUTPUT, only : pack_sky,print_pack_sky,skyp,SKYP_TYPE
  
  implicit none 

  integer, intent(in) :: nrad
  real, dimension(nrad), intent(in) :: opress
  real, dimension(nrad), intent(in) :: otempk
  real, dimension(nrad), intent(in) :: osh_v
  real, dimension(nrad), intent(in) :: ozone
  real, dimension(nrad), intent(in) :: cloud_water_content
  real, dimension(nrad), intent(in) :: cloud_effective_radius
  integer, dimension(nrad), intent(in) :: cloud_water_phase
  real, dimension(nrad), intent(in) :: ztl
  real, dimension(nrad), intent(inout) :: flxus
  real, dimension(nrad), intent(inout) :: flxds
  real, dimension(nrad), intent(inout) :: flxul
  real, dimension(nrad), intent(inout) :: flxdl
  real, intent(in) :: skin_temp
  real, intent(in) :: alb
  real, intent(in) :: zm
  real, intent(in) :: solar1
  real, intent(in) :: cosz
  real, intent(in) :: emis
  integer, intent(in) :: nac
  integer, dimension(nac), optional, intent(in) :: itps
  integer, intent(in) :: n_atau
  real, dimension(n_atau), optional, intent(in) :: a_wli
  real, dimension(n_atau,nac), optional, intent(in) :: a_taus
  real, dimension(nrad,nac), optional, intent(in) :: aprof
  integer, intent(in) :: iw
  TYPE (SKYP_TYPE) ut,tu
  integer k
  real psfc
  integer :: iac
  integer :: itau
  !=================================================================
  
  ! Sets some of the more obsure inputs to reasonable values.
  ! This is where you can modify CO2 concentrations!!!
  call set_default_options_fu

  !Input profile assignment
  call load_olam_data(nrad, opress, otempk, osh_v, ozone, skin_temp, alb,  &
       zm, solar1, cosz, emis, cloud_water_content,   &
       cloud_effective_radius, cloud_water_phase,   &
       fi%vi%nlev, fi%vi%pp, fi%vi%pt, fi%vi%ph, fi%vi%po, fi%pts, fi%sfcalb, &
       gflq%hsfc, fi%vi%hsfc, fi%ss, fi%u0, fi%ee, fi%fc(1)%dpi%plwc_iwc,  &
       fi%fc(1)%dpi%pre_de, fi%fc(1)%dpi%phase)

  !-------Cloud definition
  fi%fc(1)%dpi%ldpi = .true.
  fi%fc(1)%cldfrac = 1.0 ! Cloud Fraction (0-1) 
  fi%fc(1)%novl = 0
  
  !Aerosols ------------------------------------------------------------
!  fi%nac = 1 ! 2 aerosol types 
!  fi%itps(1) = 2 ! Continental see types (1-18)
!  fi%n_atau = 1	! 1 wavelength input for aerosols
!  fi%a_wli(1) = 0.535 ! AOT wavelength(microns) of a_taus
!  fi%a_taus(1,1) = 0.40 ! AOT for constituent 1
  fi%nac = nac
  do iac = 1, fi%nac
     fi%itps(iac) = itps(iac)
  enddo

  fi%n_atau = n_atau
  if(nac > 0)then
     do itau = 1, fi%n_atau
        fi%a_wli(1) = a_wli(itau)
        do iac = 1, fi%nac
           fi%a_taus(itau,iac) = a_taus(itau,iac)
        enddo
     enddo
  endif
  do iac = 1, fi%nac
     fi%aprofs(1:nrad,iac) = aprof(1:nrad,iac)
  enddo

  !----------------------------------------------------------------------
  
  ! Define model Fixed layer structure pre-cloud by fixed DZ intervals...
!  gflq%mode = 'CALIP'
!  gflq%mode = 'CERES'
!  call generate_level_scheme 
  call olam_level_scheme(ztl)

  ! CALL After all FI%VD and FI%VI structures are defined.
  call prepare_model_profile_fu

  ! uses FI%VO !! Assign Model ATM Profile and CLD Levels
  call vla_interface_fu     

  !Aerosol Profile (after fi%pp is created )-----------------------------
!  call aer_scale_hgt(fi%nv,fi%pp,3.0,fi%aprofs(1:fi%nv,1) )

  ! RADIATVE TRANSFER --------------------------------------------------
  call rad_multi_fu(iw)  ! CALL THE CODE !!!

  !--------------------------------------------------------------------
!  if(iw == 279)then
!     call pack_sky
!     call print_pack_sky
!     stop
!  endif

  call pack_sky_olam(nrad, flxus, flxds, flxul, flxdl)

!  ut= skyp
!  write(11) ut
!  print'(4f10.2)',fsfc%swdir
!  print'(4f10.2)',fsfc%dirsfc
  
!  print'(4f10.2)',fsfc%swdif
!  print'(4f10.2)',fsfc%difsfc
  
!  print'(4f10.2)',fsfc%swdir + fsfc%swdif
  ! Aerosols Only; Clouds and Aerosols; None; Clouds Only
!  print'(5f10.2)',cosz,fsfc%dirsfc + fsfc%difsfc

  return
end subroutine fuliou

!====================================================================

subroutine load_olam_data(nrad, opress, otempk, osh_v, ozone, skin_temp,   &
     alb, zm, solar1, cosz, emis, cloud_water_content,   &
     cloud_effective_radius, cloud_water_phase,  &
     nlev, pp, pt, ph, po, pts, sfcalb, hsfc_surf, hsfc_prof, ss, u0, ee, &
     plwc_iwc, pre_de, phase)
  
  use fuinput, only: mccx

  implicit none

  integer, intent(in) :: nrad
  real, dimension(nrad), intent(in) :: opress
  real, dimension(nrad), intent(in) :: otempk
  real, dimension(nrad), intent(in) :: osh_v
  real, dimension(nrad), intent(in) :: ozone
  real, dimension(nrad), intent(in) :: cloud_water_content
  real, dimension(nrad), intent(in) :: cloud_effective_radius
  real, dimension(nrad), intent(in) :: cloud_water_phase
  real, intent(in) :: skin_temp
  real, intent(in) :: alb
  real, intent(in) :: zm
  real, intent(in) :: solar1
  real, intent(in) :: cosz
  real, intent(in) :: emis

  integer, intent(out) :: nlev
  real, dimension(nrad), intent(out) :: pp
  real, dimension(nrad), intent(out) :: pt
  real, dimension(nrad), intent(out) :: ph
  real, dimension(nrad), intent(out) :: po
  real, dimension(nrad), intent(out) :: plwc_iwc
  real, dimension(nrad), intent(out) :: pre_de
  real, dimension(nrad), intent(out) :: phase
  real, intent(out) :: pts
  real, dimension(15,2,0:mccx), intent(out) :: sfcalb
  real, intent(out) :: hsfc_surf
  real, intent(out) :: hsfc_prof
  real, intent(out) :: ss
  real, intent(out) :: u0
  real, dimension(12), intent(out) :: ee

  integer :: olam_k
  integer :: vla_k
  integer :: k

  ! Number of OLAM levels
  nlev = nrad

  ! Atmosphere
  do k = 1, nrad
     pp(k) = opress(k) ! mb
     pt(k) = otempk(k) ! K
     ph(k) = max(1.0e-10,osh_v(k)) ! g/g
     po(k) = ozone(k) ! g/g
  enddo
  pts = skin_temp ! K

  ! Albedo
  !   index 1: spectral band [1-15]
  !   index 2: clear sky [1]; pristine sky [2]
  !   index 3: 0 -> no clouds, 1:mccx -> different cloud conditions
  sfcalb = alb

  hsfc_surf = zm
  hsfc_prof = zm

  ss = solar1
  u0 = cosz

  ! Spectral surface emissivity for 12 long wave bands
  ee(1:12) = emis

  ! Clouds
  do k = 1, nrad
     plwc_iwc(k) = cloud_water_content(k) ! [kg/m3]
     pre_de(k) = cloud_effective_radius(k) ! liquid: 5-30; ice: 20-180 [um]
     phase(k) = cloud_water_phase(k) ! 1=liquid; 2=ice
  enddo

  return
end subroutine load_olam_data

!======================================================================
subroutine olam_level_scheme(ztl)

  use fulioumulti, only: fi
  use generate_fuliou_levels

  implicit none

  real, dimension(fi%vi%nlev) :: ztl

  fi%hsfc = gflq%hsfc
  gflq%nlev = fi%vi%nlev
  pp(1:gflq%nlev) = fi%vi%pp(1:gflq%nlev)
  hh(gflq%nlev) = fi%hsfc
  hh(1:(gflq%nlev-1)) = ztl(1:(gflq%nlev-1))

  fi%vd%nfix = gflq%nlev
  fi%vd%pfix(1:fi%vd%nfix) = pp(1:gflq%nlev)
  fi%vd%nflo = 0

  fi%vd%nrep = fi%vd%nfix
  fi%vd%report(1:fi%vd%nfix) = fi%vd%pfix(1:fi%vd%nfix)
  fi%vi%hsfc = gflq%hsfc

  return
end subroutine olam_level_scheme

!====================================================================
subroutine pack_sky_olam(nrad, flxus, flxds, flxul, flxdl)
  
  use fuinput, only: fi
  use fuoutput, only: fo

  implicit none

  integer, intent(in) :: nrad
  real, dimension(nrad), intent(inout) :: flxus
  real, dimension(nrad), intent(inout) :: flxds
  real, dimension(nrad), intent(inout) :: flxul
  real, dimension(nrad), intent(inout) :: flxdl
  integer :: k
  integer :: kk

  do k = 1, fi%vd%nrep
     kk = fi%vo%ireport(k)

     flxus(k) = fo(2)%fus(kk)
     flxds(k) = fo(2)%fds(kk)

     flxul(k) = fo(2)%fuir(kk)
     flxdl(k) = fo(2)%fdir(kk)

  enddo

  return
end subroutine pack_sky_olam

