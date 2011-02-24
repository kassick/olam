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
Module mem_sea

   use sea_coms, only: maxjms
   use max_dims, only: maxremote

!----------------------------------------------------------------------------

   Type itab_ms_vars
      integer :: imglobe = 1
   End type

   Type itab_us_vars
      integer :: im1 = 1, im2 = 1
      integer :: iw1 = 1, iw2 = 1
      integer :: irank = -1
      integer :: iuglobe = 1
   End type

   Type itab_ws_vars
      logical, allocatable :: send_s(:)

      integer :: im(maxjms) = 1
      integer :: iu(maxjms) = 1
      integer :: irank = -1
      integer :: iwglobe = 1
      integer :: jm = 0
   End type

   type (itab_ms_vars), allocatable :: itab_ms(:)
   type (itab_us_vars), allocatable :: itab_us(:)
   type (itab_ws_vars), allocatable :: itab_ws(:)

   type (itab_ms_vars), allocatable :: ltab_ms(:)
   type (itab_us_vars), allocatable :: ltab_us(:)
   type (itab_ws_vars), allocatable :: ltab_ws(:)

!----------------------------------------------------------------------------

   Type itabg_ms_vars
      integer :: ims_myrank = -1
      integer :: irank = -1
   End Type itabg_ms_vars

   Type itabg_us_vars
      integer :: ius_myrank = -1
      integer :: irank = -1
   End Type itabg_us_vars

   Type itabg_ws_vars
      integer :: iws_myrank = -1
      integer :: irank = -1
   End Type itabg_ws_vars

   type (itabg_ms_vars), allocatable         :: itabg_ms(:)
   type (itabg_us_vars), allocatable         :: itabg_us(:)
   type (itabg_ws_vars), allocatable, target :: itabg_ws(:)

!----------------------------------------------------------------------------

   Type jtab_ws_mpi_vars
      integer, allocatable :: iws(:)
      integer, allocatable :: jend(:)
   End Type

   type (jtab_ws_mpi_vars) :: jtab_ws_mpi(maxremote)

!----------------------------------------------------------------------------

   Type sea_vars

      real, allocatable :: area       (:) ! sea cell surface area [m^2]
      real, allocatable :: rhos       (:) ! air density [kg/m^3]
      real, allocatable :: ustar      (:) ! friction velocity [m/s]
      real, allocatable :: sxfer_t    (:) ! can_air-to-atm heat xfer this step [kg_air K/m^2]
      real, allocatable :: sxfer_r    (:) ! can_air-to-atm vapor xfer this step [kg_vap/m^2]
      real, allocatable :: can_depth  (:) ! "canopy" depth for heat & vap capacity [m]
      real, allocatable :: seatp      (:) ! past sea temperature (obs time) [K]
      real, allocatable :: seatf      (:) ! future sea temperature (obs time) [K]
      real, allocatable :: seatc      (:) ! current sea temperature [K]
      real, allocatable :: seaicep    (:) ! past seaice fraction (obs time) [0-1]
      real, allocatable :: seaicef    (:) ! future seaice fraction (obs time) [0-1]
      real, allocatable :: seaicec    (:) ! seaice fraction [0-1]
      real, allocatable :: can_temp   (:) ! "canopy" air temperature [K]
      real, allocatable :: can_shv    (:) ! "canopy" vapor spec hum [kg_vap/kg_air]
      real, allocatable :: surface_ssh(:) ! sea surface sat spec hum [kg_vap/kg_air]
      real, allocatable :: rough      (:) ! water surface roughness height [m]

      real, allocatable :: rshort        (:) ! downward can-top s/w flux [W/m^2]
      real, allocatable :: rshort_diffuse(:) ! downward diffuse can-top s/w flux [W/m2]
      real, allocatable :: rlong         (:) ! downward can-top l/w flux [W/m^2]
      real, allocatable :: rlongup       (:) ! upward can-top l/w flux [W/m^2]
      real, allocatable :: rlong_albedo  (:) ! water l/w albedo [0-1]
      real, allocatable :: albedo_beam   (:) ! water s/w beam albedo [0-1]
      real, allocatable :: albedo_diffuse(:) ! water s/w diffuse albedo [0-1]

      real, allocatable :: pcpg   (:) ! new pcp amount this timestep [kg/m^2]
      real, allocatable :: qpcpg  (:) ! new pcp energy this timestep [J/m^2]
      real, allocatable :: dpcpg  (:) ! new pcp depth this timestep [m]

      real, allocatable :: xems(:) ! earth x coord of sea cell M points
      real, allocatable :: yems(:) ! earth y coord of sea cell M points
      real, allocatable :: zems(:) ! earth z coord of sea cell M points


      real, allocatable :: xews(:) ! earth x coord of sea cell W points
      real, allocatable :: yews(:) ! earth y coord of sea cell W points
      real, allocatable :: zews(:) ! earth z coord of sea cell W points

   End Type
  
   type (sea_vars) :: sea

Contains

!=========================================================================

   subroutine alloc_sea_grid(mms,mus,mws,igroupsize)

   implicit none

   integer, intent(in) :: mms,mus,mws,igroupsize
   integer :: iws

! Allocate and initialize sea arrays

   allocate (itab_ms(mms))
   allocate (itab_us(mus))
   allocate (itab_ws(mws))

   do iws = 1,mws
      allocate (itab_ws(iws)%send_s (igroupsize))

      itab_ws(iws)%send_s (1:igroupsize) = .false.
   enddo

   allocate (sea%area(mws))

   allocate (sea%xems(mms))
   allocate (sea%yems(mms))
   allocate (sea%zems(mms))

   allocate (sea%xews(mws))
   allocate (sea%yews(mws))
   allocate (sea%zews(mws))

   sea%area(1:mws) = 0.

   sea%xems(1:mms) = 0.
   sea%yems(1:mms) = 0.
   sea%zems(1:mms) = 0.

   sea%xews(1:mws) = 0.
   sea%yews(1:mws) = 0.
   sea%zews(1:mws) = 0.

   return
   end subroutine alloc_sea_grid

!=========================================================================

   subroutine alloc_sea(mws)

   implicit none

   integer, intent(in) :: mws

! Allocate sea arrays

   allocate (sea%rhos          (mws))
   allocate (sea%ustar         (mws))
   allocate (sea%sxfer_t       (mws))
   allocate (sea%sxfer_r       (mws))
   allocate (sea%can_depth     (mws))
   allocate (sea%seatp         (mws))
   allocate (sea%seatf         (mws))
   allocate (sea%seatc         (mws))
   allocate (sea%seaicep       (mws))
   allocate (sea%seaicef       (mws))
   allocate (sea%seaicec       (mws))
   allocate (sea%can_temp      (mws))
   allocate (sea%can_shv       (mws))
   allocate (sea%surface_ssh   (mws))
   allocate (sea%rough         (mws))

   allocate (sea%rshort        (mws))
   allocate (sea%rshort_diffuse(mws))
   allocate (sea%rlong         (mws))
   allocate (sea%rlongup       (mws))
   allocate (sea%rlong_albedo  (mws))
   allocate (sea%albedo_beam   (mws))
   allocate (sea%albedo_diffuse(mws))

   allocate (sea%pcpg          (mws))
   allocate (sea%qpcpg         (mws))
   allocate (sea%dpcpg         (mws))

! Initialize sea arrays
   
   sea%rhos          (1:mws) = 0.
   sea%ustar         (1:mws) = 0.
   sea%sxfer_t       (1:mws) = 0.
   sea%sxfer_r       (1:mws) = 0.
   sea%can_depth     (1:mws) = 0.
   sea%seatp         (1:mws) = 0.
   sea%seatf         (1:mws) = 0.
   sea%seatc         (1:mws) = 0.
   sea%seaicep       (1:mws) = 0.
   sea%seaicef       (1:mws) = 0.
   sea%seaicec       (1:mws) = 0.
   sea%can_temp      (1:mws) = 0.
   sea%can_shv       (1:mws) = 0.
   sea%surface_ssh   (1:mws) = 0.
   sea%rough         (1:mws) = 0.

   sea%rshort        (1:mws) = 0.
   sea%rshort_diffuse(1:mws) = 0.
   sea%rlong         (1:mws) = 0.
   sea%rlongup       (1:mws) = 0.
   sea%rlong_albedo  (1:mws) = 0.
   sea%albedo_beam   (1:mws) = 0.
   sea%albedo_diffuse(1:mws) = 0.

   sea%pcpg          (1:mws) = 0.
   sea%qpcpg         (1:mws) = 0.
   sea%dpcpg         (1:mws) = 0.

   return
   end subroutine alloc_sea

!=========================================================================

   subroutine filltab_sea(mms,mus,mws)

   use var_tables, only: vtables
   use sea_coms,   only: maxjms  
   use misc_coms,  only: iparallel, runtype

   implicit none
   
   integer, intent(in) :: mms,mus,mws
   
   integer :: ndims
   integer :: idims(2)
   character(20) :: action

   ndims = 1
   idims(1) = mms
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
      if (allocated(itab_ms)) call vtables(ndims,idims  &
           ,'IMGLOBE_S '//trim(action)    &
           ,ivara1=itab_ms(:)%imglobe)
!------------------------------------------------------------------------

      idims(1) = mus

!------------------------------------------------------------------------
      if (allocated(itab_us)) call vtables(ndims,idims  &
           ,'IUGLOBE_S '//trim(action)    &
           ,ivara1=itab_us(:)%iuglobe)
!------------------------------------------------------------------------
      if (allocated(itab_us)) call vtables(ndims,idims  &
           ,'IRANKU_S '//trim(action)     &
           ,ivara1=itab_us(:)%irank)
!------------------------------------------------------------------------

      idims(1) = mws

!------------------------------------------------------------------------
      if (allocated(itab_ws)) call vtables(ndims,idims  &
           ,'IWGLOBE_S '//trim(action)    &
           ,ivara1=itab_ws(:)%iwglobe)
!------------------------------------------------------------------------
      if (allocated(itab_ws)) call vtables(ndims,idims  &
           ,'IRANKW_S '//trim(action)     &
           ,ivara1=itab_ws(:)%irank)
!------------------------------------------------------------------------
   endif

   ndims = 1
   idims(1) = mws
   idims(2) = 1

!------------------------------------------------------------------------
   if (allocated(sea%sxfer_t))   call vtables(ndims,idims,'SEA%SXFER_T  :hist' &
         ,rvara1=sea%sxfer_t)
!------------------------------------------------------------------------
   if (allocated(sea%sxfer_r))   call vtables(ndims,idims,'SEA%SXFER_R  :hist' &
         ,rvara1=sea%sxfer_r)
!------------------------------------------------------------------------
   if (allocated(sea%can_depth)) call vtables(ndims,idims,'SEA%CAN_DEPTH:hist' &
         ,rvara1=sea%can_depth)
!------------------------------------------------------------------------
   if (allocated(sea%seatc))     call vtables(ndims,idims,'SEA%SEATC    :hist' &
         ,rvara1=sea%seatc)
!------------------------------------------------------------------------
   if (allocated(sea%seaicec))   call vtables(ndims,idims,'SEA%SEAICEC  :hist' &
         ,rvara1=sea%seaicec)
!------------------------------------------------------------------------
   if (allocated(sea%can_temp))  call vtables(ndims,idims,'SEA%CAN_TEMP :hist' &
         ,rvara1=sea%can_temp)
!------------------------------------------------------------------------
   if (allocated(sea%can_shv))   call vtables(ndims,idims,'SEA%CAN_SHV  :hist' &
         ,rvara1=sea%can_shv)
!------------------------------------------------------------------------
   if (allocated(sea%surface_ssh)) call vtables(ndims,idims  &
               ,'SEA%SURFACE_SSH :hist' &
         ,rvara1=sea%surface_ssh)
!------------------------------------------------------------------------
   if (allocated(sea%rough))     call vtables(ndims,idims,'SEA%ROUGH    :hist' &
         ,rvara1=sea%rough)
!------------------------------------------------------------------------
   if (allocated(sea%rshort))     call vtables(ndims,idims &
               ,'SEA%RSHORT  :hist' &
         ,rvara1=sea%rshort)
!------------------------------------------------------------------------
   if (allocated(sea%rshort_diffuse))     call vtables(ndims,idims &
               ,'SEA%RSHORT_DIFFUSE  :hist' &
         ,rvara1=sea%rshort_diffuse)
!------------------------------------------------------------------------
   if (allocated(sea%rlong))     call vtables(ndims,idims &
               ,'SEA%RLONG  :hist' &
         ,rvara1=sea%rlong)
!------------------------------------------------------------------------
   if (allocated(sea%rlongup))     call vtables(ndims,idims &
               ,'SEA%RLONGUP  :hist' &
         ,rvara1=sea%rlongup)
!------------------------------------------------------------------------
   if (allocated(sea%rlong_albedo))     call vtables(ndims,idims  &
               ,'SEA%RLONG_ALBEDO    :hist' &
         ,rvara1=sea%rlong_albedo)
!------------------------------------------------------------------------
   if (allocated(sea%albedo_beam))     call vtables(ndims,idims  &
               ,'SEA%ALBEDO_BEAM    :hist' &
         ,rvara1=sea%albedo_beam)
!------------------------------------------------------------------------
   if (allocated(sea%albedo_diffuse))     call vtables(ndims,idims  &
               ,'SEA%ALBEDO_DIFFUSE:hist' &
         ,rvara1=sea%albedo_diffuse)
!------------------------------------------------------------------------
   if (allocated(sea%pcpg))     call vtables(ndims,idims  &
               ,'SEA%PCPG:hist' &
         ,rvara1=sea%pcpg)
!------------------------------------------------------------------------
   if (allocated(sea%qpcpg))     call vtables(ndims,idims  &
               ,'SEA%QPCPG:hist' &
         ,rvara1=sea%qpcpg)
!------------------------------------------------------------------------
   if (allocated(sea%dpcpg))     call vtables(ndims,idims  &
               ,'SEA%DPCPG:hist' &
         ,rvara1=sea%dpcpg)
!------------------------------------------------------------------------

   return
   end subroutine filltab_sea

!===============================================================================

   subroutine fill_jsea()

   use mem_ijtabs, only: mrls

   use misc_coms,  only: io6, iparallel

   use mem_para,   only: mgroupsize, myrank,  &
!                        send_us, recv_us,  &
                         send_ws, recv_ws,  &
                         send_wsf, recv_wsf,  &
!                        nsends_us, nrecvs_us,  &
                         nsends_ws, nrecvs_ws,  &
                         nsends_wsf, nrecvs_wsf

   use sea_coms,   only: mms, mus, mws, nms, nus, nws

   implicit none
   
   integer :: jsend,iws,jend,mrl

! Allocate and zero-fill JTAB_WS_MPI%JEND

   do jsend = 1,maxremote
      allocate (jtab_ws_mpi(jsend)%jend(mrls))
                jtab_ws_mpi(jsend)%jend(1:mrls) = 0
   enddo

! Return if run is not parallel (jtab not needed)

   if (iparallel == 0) return
   
! Compute and store JTAB_WS_MPI%JEND(1)

   do jsend = 1,nsends_ws(1)
      jtab_ws_mpi(jsend)%jend(1) = 0
      do iws = 2,mws
         if (itab_ws(iws)%send_s(jsend)) then
            jtab_ws_mpi(jsend)%jend(1) = jtab_ws_mpi(jsend)%jend(1) + 1
         endif
      enddo
      jtab_ws_mpi(jsend)%jend(1) = max(1,jtab_ws_mpi(jsend)%jend(1))
   enddo

! Allocate and zero-fill JTAB_WS_MPI%IWS

   do jsend = 1,nsends_ws(1)
      jend = jtab_ws_mpi(jsend)%jend(1)
      allocate (jtab_ws_mpi(jsend)%iws(jend))
                jtab_ws_mpi(jsend)%iws(1:jend) = 0
   enddo

! Initialize JTAB_WS_MPI%JEND counters to zero

   do jsend = 1,nsends_ws(1)
      jtab_ws_mpi(jsend)%jend(1:mrls) = 0
   enddo

! Compute JTAB_WS_MPI%IWS

   do mrl = mrls,1,-1
      do iws = 2,mws
         do jsend = 1,nsends_ws(1)

            if (itab_ws(iws)%send_s(jsend)) then
               jtab_ws_mpi(jsend)%jend(1:mrl) = jtab_ws_mpi(jsend)%jend(1:mrl) + 1
               jtab_ws_mpi(jsend)%iws(jtab_ws_mpi(jsend)%jend(1)) = iws
            endif

         enddo
      enddo
   enddo

   return
   end subroutine fill_jsea

End Module mem_sea
