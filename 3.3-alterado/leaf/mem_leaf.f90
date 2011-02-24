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
Module mem_leaf

   use ed_structure_defs
   use leaf_coms, only: maxjml
   use max_dims,  only: maxremote

!----------------------------------------------------------------------------

   Type itab_ml_vars
      integer :: imglobe = 1
   End type

   Type itab_ul_vars
      integer :: im1 = 1, im2 = 1
      integer :: iw1 = 1, iw2 = 1
      integer :: irank = -1
      integer :: iuglobe = 1
   End type

   Type itab_wl_vars
      logical, allocatable :: send_l(:)

      integer :: im(maxjml) = 1
      integer :: iu(maxjml) = 1
      integer :: irank = -1
      integer :: iwglobe = 1
      integer :: jm = 0
   End type

   type (itab_ml_vars), allocatable :: itab_ml(:)
   type (itab_ul_vars), allocatable :: itab_ul(:)
   type (itab_wl_vars), allocatable :: itab_wl(:)

   type (itab_ml_vars), allocatable :: ltab_ml(:)
   type (itab_ul_vars), allocatable :: ltab_ul(:)
   type (itab_wl_vars), allocatable :: ltab_wl(:)

!----------------------------------------------------------------------------

   Type itabg_ml_vars
      integer :: iml_myrank = -1
      integer :: irank = -1
   End Type itabg_ml_vars

   Type itabg_ul_vars
      integer :: iul_myrank = -1
      integer :: irank = -1
   End Type itabg_ul_vars

   Type itabg_wl_vars
      integer :: iwl_myrank = -1
      integer :: irank = -1
   End Type itabg_wl_vars

   type (itabg_ml_vars), allocatable         :: itabg_ml(:)
   type (itabg_ul_vars), allocatable         :: itabg_ul(:)
   type (itabg_wl_vars), allocatable, target :: itabg_wl(:)

!----------------------------------------------------------------------------

   Type jtab_wl_mpi_vars
      integer, allocatable :: iwl(:)
      integer, allocatable :: jend(:)
   End Type

   type (jtab_wl_mpi_vars) :: jtab_wl_mpi(maxremote) 

!----------------------------------------------------------------------------

   Type land_vars

      integer, allocatable :: nlev_sfcwater(:)  ! # of active surface water levels
      integer, allocatable :: leaf_class   (:)  ! leaf ("vegetation") class
      integer, allocatable :: ntext_soil (:,:)  ! soil textural class

      real, allocatable :: soil_water     (:,:) ! soil water content [vol_water/vol_tot]
      real, allocatable :: soil_energy    (:,:) ! soil energy [J/m^3]
      real, allocatable :: sfcwater_mass  (:,:) ! surface water mass [kg/m^2]
      real, allocatable :: sfcwater_energy(:,:) ! surface water energy [J/kg]
      real, allocatable :: sfcwater_depth (:,:) ! surface water depth [m]
      real, allocatable :: rshort_s       (:,:) ! s/w net rad flux to sfc water [W/m^2]

      real, allocatable :: rshort        (:) ! downward can-top s/w flux [W/m^2]
      real, allocatable :: rshort_diffuse(:) ! downward diffuse can-top s/w flux [W/m2]
      real, allocatable :: rlong         (:) ! downward can-top l/w flux [W/m^2]
      real, allocatable :: rlongup       (:) ! upward can-top l/w flux [W/m^2]
      real, allocatable :: rlong_albedo  (:) ! net canopy/soil l/w albedo [0-1]
      real, allocatable :: albedo_beam   (:) ! net canopy/soil s/w beam albedo [0-1]
      real, allocatable :: albedo_diffuse(:) ! net canopy/soil s/w diffuse albedo [0-1]

      real, allocatable :: rshort_g    (:) ! s/w net rad flux to soil [W/m^2]
      real, allocatable :: rshort_v    (:) ! s/w net rad flux to veg [W/m^2]
      real, allocatable :: rlong_g     (:) ! l/w net rad flux to soil [W/m^2]
      real, allocatable :: rlong_s     (:) ! l/w net rad flux to sfc water [W/m^2]
      real, allocatable :: rlong_v     (:) ! l/w net rad flux to veg [W/m^2]
      real, allocatable :: area        (:) ! land cell area [m^2]
      real, allocatable :: rhos        (:) ! atmosphere air density [kg_air/m^3]
      real, allocatable :: vels        (:) ! surface wind speed [m/s]
      real, allocatable :: prss        (:) ! air pressure [Pa]
      real, allocatable :: ustar       (:) ! friction velocity [m/s]
      real, allocatable :: sxfer_t     (:) ! canair-to-atm heat xfer this step [kg_air K/m^2]
      real, allocatable :: sxfer_r     (:) ! canair-to-atm vapor xfer this step [kg_vap/m^2]
      real, allocatable :: can_depth   (:) ! canopy depth for heat & vap capacity [m]
      real, allocatable :: hcapveg     (:) ! veg heat capacity [J/(m^2 K)]
      real, allocatable :: can_temp    (:) ! canopy air temp [K]
      real, allocatable :: can_shv     (:) ! canopy air vapor spec hum [kg_vap/kg_air]
      real, allocatable :: surface_ssh (:) ! surface saturation spec hum [kg_vap/kg_air]
      real, allocatable :: ground_shv  (:) ! soil vapor spec hum [kg_vap/kg_air]
      real, allocatable :: rough       (:) ! land cell roughness height [m]
      real, allocatable :: veg_fracarea(:) ! veg fractional area
      real, allocatable :: veg_lai     (:) ! veg leaf area index
      real, allocatable :: veg_rough   (:) ! veg roughness height [m]
      real, allocatable :: veg_height  (:) ! veg height [m]
      real, allocatable :: veg_albedo  (:) ! veg albedo
      real, allocatable :: veg_tai     (:) ! veg total area index
      real, allocatable :: veg_water   (:) ! veg sfc water content [kg/m^2]
      real, allocatable :: veg_temp    (:) ! veg temp [K]
      real, allocatable :: veg_ndvip   (:) ! veg past ndvi (obs time)
      real, allocatable :: veg_ndvif   (:) ! veg future ndvi (obs time)
      real, allocatable :: veg_ndvic   (:) ! veg current ndvi
      real, allocatable :: stom_resist (:) ! veg stomatal resistance [s/m]
      real, allocatable :: snowfac     (:) ! frac veg burial by snowcover
      real, allocatable :: vf          (:) ! frac coverage of non-buried part of veg
      real, allocatable :: pcpg        (:) ! new pcp amount this leaf timestep [kg/m^2]
      real, allocatable :: qpcpg       (:) ! new pcp energy this leaf timestep [J/m^2]
      real, allocatable :: dpcpg       (:) ! new pcp depth this leaf timestep [m]

      real, allocatable :: xeml     (:) ! earth x coord of land M points
      real, allocatable :: yeml     (:) ! earth y coord of land M points
      real, allocatable :: zeml     (:) ! earth z coord of land M points

      real, allocatable :: xewl     (:) ! earth x coord of land W points
      real, allocatable :: yewl     (:) ! earth y coord of land W points
      real, allocatable :: zewl     (:) ! earth z coord of land W points

      real, allocatable :: cosz        (:) ! cosine of the solar zenith angle of the land cell

      integer, allocatable :: ed_flag  (:) ! 0 or 1 flag specifying if ED is being run in this cell.

      real, allocatable :: gpp(:) ! Gross primary productivity (umol/m2/s)

      real, allocatable :: rh(:) ! Heterotrophic respiration (umol/m2/s)

      real, allocatable :: nep(:) ! Net ecosystem productivity (umol/m2/s)

      real, allocatable :: agb(:) ! Site above ground biomass (kgC/m2)
      real, allocatable :: basal_area(:) ! Basal area (m2/ha)
      real, allocatable :: agb_growth(:) ! Net plant growth rate [kgC/m2/y]
      real, allocatable :: agb_mort(:) ! Net plant mortality rate [kgC/m2/y]
      real, allocatable :: agb_cut(:) ! Net plant cut rate [kgC/m2/y]
      real, allocatable :: agb_recruit(:) ! Net plant recruitment rate [kgC/m2/y]

      integer, allocatable :: lsl(:)       ! Lowest soil layer for a given iland cell.

      type(site), pointer :: first_site    ! The basic unit of the ED model is the site, which loosely corresponds to a grid cell.  Now, the entire ED memory structure resides on linked lists.  This is a pointer to the first site of the linked list. 

   End Type

   type (land_vars) :: land

Contains

!=========================================================================

   subroutine alloc_land_grid(mml,mul,mwl,nzg,igroupsize)

   implicit none

   integer, intent(in) :: mml,mul,mwl,nzg,igroupsize
   integer :: iwl

! Allocate and initialize land grid arrays

   allocate (itab_ml(mml))
   allocate (itab_ul(mul))
   allocate (itab_wl(mwl))

   do iwl = 1,mwl
      allocate (itab_wl(iwl)%send_l (igroupsize))

      itab_wl(iwl)%send_l (1:igroupsize) = .false.
   enddo

   allocate (land%leaf_class(mwl))
   allocate (land%area      (mwl))

   allocate (land%ntext_soil(nzg,mwl))

   allocate (land%xeml(mml))
   allocate (land%yeml(mml))
   allocate (land%zeml(mml))

   allocate (land%xewl(mwl))
   allocate (land%yewl(mwl))
   allocate (land%zewl(mwl))

   land%leaf_class(1:mwl) = 0.
   land%area      (1:mwl) = 0.

   land%ntext_soil(1:nzg,1:mwl) = 0.

   land%xeml(1:mml) = 0.
   land%yeml(1:mml) = 0.
   land%zeml(1:mml) = 0.

   land%xewl(1:mwl) = 0.
   land%yewl(1:mwl) = 0.
   land%zewl(1:mwl) = 0.

   return
   end subroutine alloc_land_grid

!=========================================================================

   subroutine alloc_leaf(mwl,nzg,nzs)

   implicit none

   integer, intent(in) :: mwl,nzg,nzs

! Allocate land arrays

   allocate (land%nlev_sfcwater (mwl))
   allocate (land%rhos          (mwl))
   allocate (land%vels          (mwl))
   allocate (land%prss          (mwl))
   allocate (land%ustar         (mwl))
   allocate (land%sxfer_t       (mwl))
   allocate (land%sxfer_r       (mwl))
   allocate (land%can_depth     (mwl))
   allocate (land%hcapveg       (mwl))
   allocate (land%can_temp      (mwl))
   allocate (land%can_shv       (mwl))
   allocate (land%surface_ssh   (mwl))
   allocate (land%ground_shv    (mwl))
   allocate (land%rough         (mwl))
   allocate (land%veg_fracarea  (mwl))
   allocate (land%veg_lai       (mwl))
   allocate (land%veg_rough     (mwl))
   allocate (land%veg_height    (mwl))
   allocate (land%veg_albedo    (mwl))
   allocate (land%veg_tai       (mwl))
   allocate (land%veg_water     (mwl))
   allocate (land%veg_temp      (mwl))
   allocate (land%veg_ndvip     (mwl))
   allocate (land%veg_ndvic     (mwl))
   allocate (land%veg_ndvif     (mwl))
   allocate (land%stom_resist   (mwl))

   allocate (land%rshort        (mwl))
   allocate (land%rshort_diffuse(mwl))
   allocate (land%rlong         (mwl))
   allocate (land%rlongup       (mwl))
   allocate (land%rlong_albedo  (mwl))
   allocate (land%albedo_beam   (mwl))
   allocate (land%albedo_diffuse(mwl))

   allocate (land%rshort_g      (mwl))
   allocate (land%rshort_v      (mwl))
   allocate (land%rlong_g       (mwl))
   allocate (land%rlong_s       (mwl))
   allocate (land%rlong_v       (mwl))
   allocate (land%snowfac       (mwl))
   allocate (land%vf            (mwl))
   allocate (land%pcpg          (mwl))
   allocate (land%qpcpg         (mwl))
   allocate (land%dpcpg         (mwl))

   allocate (land%ed_flag       (mwl))

   allocate (land%cosz          (mwl))
   allocate (land%gpp           (mwl))

   allocate (land%agb           (mwl))
   allocate (land%basal_area    (mwl))
   allocate (land%agb_growth    (mwl))
   allocate (land%agb_mort      (mwl))
   allocate (land%agb_cut       (mwl))
   allocate (land%agb_recruit   (mwl))

   allocate (land%rh            (mwl))
   allocate (land%nep           (mwl))

   allocate (land%soil_water (nzg,mwl))
   allocate (land%soil_energy(nzg,mwl))

   allocate (land%lsl(mwl))

   allocate (land%sfcwater_mass  (nzs,mwl))
   allocate (land%sfcwater_energy(nzs,mwl))
   allocate (land%sfcwater_depth (nzs,mwl))
   allocate (land%rshort_s       (nzs,mwl))

   nullify(land%first_site)

! Initialize land arrays

   land%nlev_sfcwater (1:mwl) = 0.
   land%rhos          (1:mwl) = 0.
   land%vels          (1:mwl) = 0.
   land%prss          (1:mwl) = 0.
   land%ustar         (1:mwl) = 0.
   land%sxfer_t       (1:mwl) = 0.
   land%sxfer_r       (1:mwl) = 0.
   land%can_depth     (1:mwl) = 0.
   land%hcapveg       (1:mwl) = 0.
   land%can_temp      (1:mwl) = 0.
   land%can_shv       (1:mwl) = 0.
   land%surface_ssh   (1:mwl) = 0.
   land%ground_shv    (1:mwl) = 0.
   land%rough         (1:mwl) = 0.
   land%veg_fracarea  (1:mwl) = 0.
   land%veg_lai       (1:mwl) = 0.
   land%veg_rough     (1:mwl) = 0.
   land%veg_height    (1:mwl) = 0.
   land%veg_albedo    (1:mwl) = 0.
   land%veg_tai       (1:mwl) = 0.
   land%veg_water     (1:mwl) = 0.
   land%veg_temp      (1:mwl) = 0.
   land%veg_ndvip     (1:mwl) = 0.
   land%veg_ndvic     (1:mwl) = 0.
   land%veg_ndvif     (1:mwl) = 0.
   land%stom_resist   (1:mwl) = 0.

   land%rshort        (1:mwl) = 0.
   land%rshort_diffuse(1:mwl) = 0.
   land%rlong         (1:mwl) = 0.
   land%rlongup       (1:mwl) = 0.
   land%rlong_albedo  (1:mwl) = 0.
   land%albedo_beam   (1:mwl) = 0.
   land%albedo_diffuse(1:mwl) = 0.

   land%rshort_g      (1:mwl) = 0.
   land%rshort_v      (1:mwl) = 0.
   land%rlong_g       (1:mwl) = 0.
   land%rlong_s       (1:mwl) = 0.
   land%rlong_v       (1:mwl) = 0.
   land%snowfac       (1:mwl) = 0.
   land%vf            (1:mwl) = 0.
   land%pcpg          (1:mwl) = 0.
   land%qpcpg         (1:mwl) = 0.
   land%dpcpg         (1:mwl) = 0.

   land%ed_flag       (1:mwl) = 0
   land%cosz          (1:mwl) = 0.

   land%agb           (1:mwl) = -99.9
   land%basal_area    (1:mwl) = -99.9
   land%agb_growth    (1:mwl) = -99.9
   land%agb_mort      (1:mwl) = -99.9
   land%agb_cut       (1:mwl) = -99.9
   land%agb_recruit   (1:mwl) = -99.9

   land%gpp           (1:mwl) = -99.9
   land%rh            (1:mwl) = -99.9
   land%nep           (1:mwl) = -99.9

   land%soil_water (1:nzg,1:mwl) = 0.
   land%soil_energy(1:nzg,1:mwl) = 0.

   land%sfcwater_mass  (1:nzs,1:mwl) = 0.
   land%sfcwater_energy(1:nzs,1:mwl) = 0.
   land%sfcwater_depth (1:nzs,1:mwl) = 0.
   land%rshort_s       (1:nzs,1:mwl) = 0.

   land%lsl(1:mwl) = 1

   return
   end subroutine alloc_leaf

!=========================================================================

   subroutine filltab_leaf(mml,mul,mwl,nzg,nzs)

   use var_tables, only: vtables
   use leaf_coms,  only: maxjml  
   use misc_coms,  only: iparallel, runtype

   implicit none
   
   integer, intent(in) :: mml,mul,mwl,nzg,nzs
   
   integer :: ndims
   integer :: idims(2)
   character(20) :: action

! LAND VARIABLES

   ndims = 1
   idims(1) = mml
   idims(2) = 1

   if (iparallel==1 .or. runtype=='PLOTONLY' .or. runtype=='PARCOMBINE') then

      if (iparallel == 1) then
         action = ':hist:nohist'
      else
         action = ':hist'
      endif

      ! THESE ONLY NEED TO BE WRITTEN TO HISTORY FILE FOR PARALLEL RUNS,
      ! AND READ FOR PLOTONLY OR PARCOMBINE RUNS.

!------------------------------------------------------------------------
      if (allocated(itab_ml)) call vtables(ndims,idims  &
           ,'IMGLOBE_L '//trim(action)    &
           ,ivara1=itab_ml(:)%imglobe)
!------------------------------------------------------------------------

      idims(1) = mul

!------------------------------------------------------------------------
      if (allocated(itab_ul)) call vtables(ndims,idims  &
           ,'IUGLOBE_L '//trim(action)    &
           ,ivara1=itab_ul(:)%iuglobe)
!------------------------------------------------------------------------
      if (allocated(itab_ul)) call vtables(ndims,idims  &
           ,'IRANKU_L '//trim(action)     &
           ,ivara1=itab_ul(:)%irank)
!------------------------------------------------------------------------

      idims(1) = mwl

!------------------------------------------------------------------------
      if (allocated(itab_wl)) call vtables(ndims,idims  &
           ,'IWGLOBE_L '//trim(action)    &
           ,ivara1=itab_wl(:)%iwglobe)
!------------------------------------------------------------------------
      if (allocated(itab_wl)) call vtables(ndims,idims  &
           ,'IRANKW_L '//trim(action)     &
           ,ivara1=itab_wl(:)%irank)
!------------------------------------------------------------------------
   endif

   ndims = 1
   idims(1) = mwl
   idims(2) = 1

!------------------------------------------------------------------------
   if (allocated(land%nlev_sfcwater)) call vtables(ndims,idims  &
               ,'LAND%NLEV_SFCWATER   :hist'    &
         ,ivara1=land%nlev_sfcwater)
!------------------------------------------------------------------------
   if (allocated(land%sxfer_t))       call vtables(ndims,idims  &
               ,'LAND%SXFER_T         :hist'    &
         ,rvara1=land%sxfer_t)
!------------------------------------------------------------------------
   if (allocated(land%sxfer_r))       call vtables(ndims,idims  &
               ,'LAND%SXFER_R         :hist'    &
         ,rvara1=land%sxfer_r)
!------------------------------------------------------------------------
   if (allocated(land%can_depth))    call vtables(ndims,idims  &
               ,'LAND%CAN_DEPTH       :hist'    &
         ,rvara1=land%can_depth)
!------------------------------------------------------------------------
   if (allocated(land%hcapveg))       call vtables(ndims,idims  &
               ,'LAND%HCAPVEG         :hist'    &
         ,rvara1=land%hcapveg)
!------------------------------------------------------------------------
   if (allocated(land%can_temp))      call vtables(ndims,idims  &
               ,'LAND%CAN_TEMP        :hist'    &
         ,rvara1=land%can_temp)
!------------------------------------------------------------------------
   if (allocated(land%can_shv))       call vtables(ndims,idims  &
               ,'LAND%CAN_SHV         :hist'    &
         ,rvara1=land%can_shv)
!------------------------------------------------------------------------
   if (allocated(land%surface_ssh))   call vtables(ndims,idims  &
               ,'LAND%SURFACE_SSH     :hist'    &
         ,rvara1=land%surface_ssh)
!------------------------------------------------------------------------
   if (allocated(land%ground_shv))    call vtables(ndims,idims  &
               ,'LAND%GROUND_SHV      :hist'    &
         ,rvara1=land%ground_shv)
!------------------------------------------------------------------------
   if (allocated(land%rough))         call vtables(ndims,idims  &
               ,'LAND%ROUGH           :hist'    &
         ,rvara1=land%rough)
!------------------------------------------------------------------------
   if (allocated(land%veg_fracarea))  call vtables(ndims,idims  &
               ,'LAND%VEG_FRACAREA    :hist'    &
         ,rvara1=land%veg_fracarea)
!------------------------------------------------------------------------
   if (allocated(land%veg_lai))       call vtables(ndims,idims  &
               ,'LAND%VEG_LAI         :hist'    &
         ,rvara1=land%veg_lai)
!------------------------------------------------------------------------
   if (allocated(land%veg_rough))     call vtables(ndims,idims  &
               ,'LAND%VEG_ROUGH       :hist'    &
         ,rvara1=land%veg_rough)
!------------------------------------------------------------------------
   if (allocated(land%veg_height))    call vtables(ndims,idims  &
               ,'LAND%VEG_HEIGHT      :hist'    &
         ,rvara1=land%veg_height)
!------------------------------------------------------------------------
   if (allocated(land%veg_albedo))    call vtables(ndims,idims  &
               ,'LAND%VEG_ALBEDO      :hist'    &
         ,rvara1=land%veg_albedo)
!------------------------------------------------------------------------
   if (allocated(land%veg_tai))       call vtables(ndims,idims  &
               ,'LAND%VEG_TAI         :hist'    &
         ,rvara1=land%veg_tai)
!------------------------------------------------------------------------
   if (allocated(land%veg_water))     call vtables(ndims,idims  &
               ,'LAND%VEG_WATER       :hist'    &
         ,rvara1=land%veg_water)
!------------------------------------------------------------------------
   if (allocated(land%veg_temp))      call vtables(ndims,idims  &
               ,'LAND%VEG_TEMP        :hist'    &
         ,rvara1=land%veg_temp)
!------------------------------------------------------------------------
   if (allocated(land%veg_ndvic))     call vtables(ndims,idims  &
               ,'LAND%VEG_NDVIC       :hist'    &
         ,rvara1=land%veg_ndvic)
!------------------------------------------------------------------------
   if (allocated(land%stom_resist))   call vtables(ndims,idims  &
               ,'LAND%STOM_RESIST     :hist'    &
         ,rvara1=land%stom_resist)
!------------------------------------------------------------------------
   if (allocated(land%rshort))        call vtables(ndims,idims  &
               ,'LAND%RSHORT          :hist'    &
         ,rvara1=land%rshort)
!------------------------------------------------------------------------
   if (allocated(land%rshort_diffuse)) call vtables(ndims,idims  &
               ,'LAND%RSHORT_DIFFUSE   :hist'    &
         ,rvara1=land%rshort_diffuse)
!------------------------------------------------------------------------
   if (allocated(land%rlong))         call vtables(ndims,idims  &
               ,'LAND%RLONG           :hist'    &
         ,rvara1=land%rlong)
!------------------------------------------------------------------------
   if (allocated(land%rlongup))       call vtables(ndims,idims  &
               ,'LAND%RLONGUP         :hist'    &
         ,rvara1=land%rlongup)
!------------------------------------------------------------------------
   if (allocated(land%rlong_albedo))  call vtables(ndims,idims  &
               ,'LAND%RLONG_ALBEDO    :hist'    &
         ,rvara1=land%rlong_albedo)
!------------------------------------------------------------------------
   if (allocated(land%albedo_beam))   call vtables(ndims,idims  &
               ,'LAND%ALBEDO_BEAM     :hist'    &
         ,rvara1=land%albedo_beam)
!------------------------------------------------------------------------
   if (allocated(land%albedo_diffuse))  call vtables(ndims,idims  &
               ,'LAND%ALBEDO_DIFFUSE    :hist'    &
         ,rvara1=land%albedo_diffuse)
!------------------------------------------------------------------------
   if (allocated(land%rshort_g))      call vtables(ndims,idims  &
               ,'LAND%RSHORT_G        :hist'    &
         ,rvara1=land%rshort_g)
!------------------------------------------------------------------------
   if (allocated(land%rshort_v))      call vtables(ndims,idims  &
               ,'LAND%RSHORT_V        :hist'    &
         ,rvara1=land%rshort_v)
!------------------------------------------------------------------------
   if (allocated(land%rlong_g))       call vtables(ndims,idims  &
               ,'LAND%RLONG_G         :hist'    &
         ,rvara1=land%rlong_g)
!------------------------------------------------------------------------
   if (allocated(land%rlong_s))       call vtables(ndims,idims  &
               ,'LAND%RLONG_S         :hist'    &
         ,rvara1=land%rlong_s)
!------------------------------------------------------------------------
   if (allocated(land%rlong_v))       call vtables(ndims,idims  &
               ,'LAND%RLONG_V         :hist'    &
         ,rvara1=land%rlong_v)
!------------------------------------------------------------------------
   if (allocated(land%snowfac))       call vtables(ndims,idims  &
               ,'LAND%SNOWFAC         :hist'    &
         ,rvara1=land%snowfac)
!------------------------------------------------------------------------
   if (allocated(land%vf))            call vtables(ndims,idims  &
               ,'LAND%VF              :hist'    &
         ,rvara1=land%vf)
!------------------------------------------------------------------------
   if (allocated(land%pcpg))          call vtables(ndims,idims  &
               ,'LAND%PCPG            :hist'    &
         ,rvara1=land%pcpg)
!------------------------------------------------------------------------
   if (allocated(land%qpcpg))         call vtables(ndims,idims  &
               ,'LAND%QPCPG           :hist'    &
         ,rvara1=land%qpcpg)
!------------------------------------------------------------------------
   if (allocated(land%dpcpg))         call vtables(ndims,idims  &
               ,'LAND%DPCPG           :hist'    &
         ,rvara1=land%dpcpg)
!------------------------------------------------------------------------
   if (allocated(land%lsl))           call vtables(ndims,idims  &
               ,'LAND%LSL             :hist'    &
         ,ivara1=land%lsl)
!--------------------------------------------------------------------------

   ndims = 2
   idims(1) = nzg
   idims(2) = mwl

!--------------------------------------------------------------------------
   if (allocated(land%ntext_soil))      call vtables(ndims,idims  &
               ,'LAND%NTEXT_SOIL        :hist'                   &
         ,ivara2=land%ntext_soil)
!--------------------------------------------------------------------------
   if (allocated(land%soil_water))      call vtables(ndims,idims  &
               ,'LAND%SOIL_WATER        :hist'    &
         ,rvara2=land%soil_water)
!--------------------------------------------------------------------------
   if (allocated(land%soil_energy))     call vtables(ndims,idims  &
               ,'LAND%SOIL_ENERGY       :hist'    &
         ,rvara2=land%soil_energy)
!--------------------------------------------------------------------------

   idims(1) = nzs
   idims(2) = mwl

!--------------------------------------------------------------------------
   if (allocated(land%sfcwater_mass))   call vtables(ndims,idims  &
               ,'LAND%SFCWATER_MASS     :hist'    &
         ,rvara2=land%sfcwater_mass)
!--------------------------------------------------------------------------
   if (allocated(land%sfcwater_energy)) call vtables(ndims,idims  &
               ,'LAND%SFCWATER_ENERGY   :hist'    &
         ,rvara2=land%sfcwater_energy)
!--------------------------------------------------------------------------
   if (allocated(land%sfcwater_depth))  call vtables(ndims,idims  &
               ,'LAND%SFCWATER_DEPTH    :hist'    &
         ,rvara2=land%sfcwater_depth)
!--------------------------------------------------------------------------
   if (allocated(land%rshort_s))        call vtables(ndims,idims  &
               ,'LAND%RSHORT_S          :hist'    &
         ,rvara2=land%rshort_s)
!--------------------------------------------------------------------------

   return
   end subroutine filltab_leaf

!=========================================================================

   subroutine filltab_ED()
   use var_tables, only: vtables_ED
   use leaf_coms, only: mwl
   implicit none
   
   integer :: ndims
   integer :: idims(2)

! LAND VARIABLES

   ndims = 1
   idims(1) = mwl
   idims(2) = 1

!------------------------------------------------------------------------
   if (allocated(land%gpp))         call vtables_ED(ndims,idims  &
               ,'LAND%GPP           :hist:mavg:yavg'    &
         ,rvara1=land%gpp)
!------------------------------------------------------------------------
   if (allocated(land%rh))          call vtables_ED(ndims,idims  &
               ,'LAND%RH            :hist:mavg:yavg'    &
         ,rvara1=land%rh)
!------------------------------------------------------------------------
   if (allocated(land%nep))         call vtables_ED(ndims,idims  &
               ,'LAND%NEP           :hist:mavg:yavg'    &
         ,rvara1=land%nep)
!------------------------------------------------------------------------
   if (allocated(land%agb))         call vtables_ED(ndims,idims  &
               ,'LAND%AGB           :yavg'    &
         ,rvara1=land%agb)
!------------------------------------------------------------------------
   if (allocated(land%basal_area))  call vtables_ED(ndims,idims  &
               ,'LAND%BASAL_AREA    :yavg'    &
         ,rvara1=land%agb)
!------------------------------------------------------------------------
   if (allocated(land%agb_growth))  call vtables_ED(ndims,idims  &
               ,'LAND%AGB_GROWTH    :yavg'    &
         ,rvara1=land%agb_growth)
!------------------------------------------------------------------------
   if (allocated(land%agb_mort))    call vtables_ED(ndims,idims  &
               ,'LAND%AGB_MORT      :yavg'    &
         ,rvara1=land%agb_mort)
!------------------------------------------------------------------------
   if (allocated(land%agb_cut))     call vtables_ED(ndims,idims  &
               ,'LAND%AGB_CUT       :yavg'    &
         ,rvara1=land%agb_cut)
!------------------------------------------------------------------------
   if (allocated(land%agb_recruit)) call vtables_ED(ndims,idims  &
               ,'LAND%AGB_RECRUIT   :yavg'    &
         ,rvara1=land%agb_recruit)
!------------------------------------------------------------------------

   return
   end subroutine filltab_ED

!===============================================================================

   subroutine fill_jland()

   use mem_ijtabs, only: mrls

   use misc_coms,  only: io6, iparallel
   
   use mem_para,   only: mgroupsize, myrank,  &
!                        send_ul, recv_ul,  &
                         send_wl, recv_wl,  &
                         send_wlf, recv_wlf,  &
!                        nsends_ul, nrecvs_ul,  &
                         nsends_wl, nrecvs_wl,  &
                         nsends_wlf, nrecvs_wlf

   use leaf_coms,   only: mml, mul, mwl, nml, nul, nwl
   
   implicit none
   
   integer :: jsend,iwl,jend,mrl

! Allocate and zero-fill JTAB_WL_MPI%JEND

   do jsend = 1,maxremote
      allocate (jtab_wl_mpi(jsend)%jend(mrls))
                jtab_wl_mpi(jsend)%jend(1:mrls) = 0
   enddo

! Return if run is not parallel (jtab not needed)

   if (iparallel == 0) return

! Compute and store JTAB_WL_MPI%JEND(1)

   do jsend = 1,nsends_wl(1)
      jtab_wl_mpi(jsend)%jend(1) = 0
      do iwl = 2,mwl
         if (itab_wl(iwl)%send_l(jsend)) then
            jtab_wl_mpi(jsend)%jend(1) = jtab_wl_mpi(jsend)%jend(1) + 1
         endif
      enddo
      jtab_wl_mpi(jsend)%jend(1) = max(1,jtab_wl_mpi(jsend)%jend(1))
   enddo

! Allocate and zero-fill JTAB_WL_MPI%IWL

   do jsend = 1,nsends_wl(1)
      jend = jtab_wl_mpi(jsend)%jend(1)
      allocate (jtab_wl_mpi(jsend)%iwl(jend))
                jtab_wl_mpi(jsend)%iwl(1:jend) = 0
   enddo

! Initialize JTAB_WL_MPI%JEND counters to zero

   do jsend = 1,nsends_wl(1)
      jtab_wl_mpi(jsend)%jend(1:mrls) = 0
   enddo

! Compute JTAB_WL_MPI%IWL

   do mrl = mrls,1,-1
      do iwl = 2,mwl
         do jsend = 1,nsends_wl(1)

            if (itab_wl(iwl)%send_l(jsend)) then
               jtab_wl_mpi(jsend)%jend(1:mrl) = jtab_wl_mpi(jsend)%jend(1:mrl) + 1
               jtab_wl_mpi(jsend)%iwl(jtab_wl_mpi(jsend)%jend(1)) = iwl
            endif

         enddo
      enddo
   enddo

   return
   end subroutine fill_jland

End Module mem_leaf
