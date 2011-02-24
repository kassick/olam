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
subroutine write_ed_output()

  use misc_coms, only: current_time, frqstate
  use mem_leaf, only: land

  implicit none
  
  integer, save :: ifirst = 1

  if(ifirst == 1)then
     call write_land_grid_info()
     ifirst = 0
  endif

  if(.not.associated(land%first_site))return

  call normalize_ed_output_vars()
  call write_ed_hdf5('hist')
!  call print_ed_soi()
  call zero_ed_output_vars()

  if(current_time%date == 1 .and. current_time%time < frqstate)then
     call normalize_ed_monthly_vars()
     call write_ed_hdf5('mavg')
     call zero_ed_monthly_vars()

     if(current_time%month == 6)then
        ! Do the annual write-out
        call update_ed_yearly_vars()
        call write_ed_hdf5('yavg')
        call zero_ed_yearly_vars()
     endif

!     if(current_time%date == 1)call ed_history_write()

  endif

  call ed_history_write()

  return
end subroutine write_ed_output

!==================================================================

subroutine write_land_grid_info()

  use mem_ijtabs, only: itab_w, itabg_w
  use mem_grid, only: mwa, nza, zem, xem, yem,xew,yew,zew, arw0,lpw, topm
  use consts_coms, only: erad,piu180
  use misc_coms, only: hfilepref, iparallel
  use mem_sea, only: sea, itab_ws
  use sea_coms, only: mws
  use mem_leaf, only: land, itab_wl, itabg_wl
  use leaf_coms, only: mwl
  use mem_sflux, only: mlandflux, landflux
  use mem_para,  only: myrank

  implicit none

  integer :: iw
  integer :: im
  real :: glatm1,glatm2,glatm3
  real :: glonm1,glonm2,glonm3
  real :: glatw,glonw
  character(len=256) :: fname
  real :: max_lon
  real :: min_lon
  real :: lon_diff1
  real :: lon_diff2
  real :: xe_mean
  real :: ye_mean
  real :: ze_mean
  real :: glat_mean
  real :: glon_mean
  integer :: iwl
  integer :: iws
  real, dimension(mwl) :: land_area
  integer :: ilf
  real :: arf_land
  integer :: im1,im2,im3
  character(len=10) :: post

!======================
! W-point areas
!======================

  write(post,'(i10)') myrank
  post = 'r'//trim(adjustl(post))

  if(iparallel == 0)then
     fname = trim(hfilepref)//'-area-grid-info.txt'
  else
     fname = trim(hfilepref)//'-area-grid-info_'//trim(post)//'.txt'
  endif

  open(12,file=trim(fname),status='replace',form='formatted')

  do iw = 1,mwa
     im1 = itab_w(iw)%im1
     im2 = itab_w(iw)%im2
     im3 = itab_w(iw)%im3
     write(12,*)arw0(iw),(topm(im1)+topm(im2)+topm(im3))/3.0
  enddo

  close(12)

!======================
! Land-point areas
!======================

  if(iparallel == 0)then
     fname = trim(hfilepref)//'-land-area-grid-info.txt'
  else
     fname = trim(hfilepref)//'-land-area-grid-info_'//trim(post)//'.txt'
  endif

  open(12,file=trim(fname),status='replace',form='formatted')


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! NOTE:  LAND_AREA SHOULD ALREADY BE AVAILABLE IN LAND%AREA(:), SO THERE SHOULD
! BE NO NEED TO EVALUATE IT IN THE FOLLOWING LOOP

  land_area(:) = 0.0
  do ilf = 2, mlandflux
     iw       = landflux(ilf)%iw
     iwl      = landflux(ilf)%iwl

! If run is parallel:  Skip flux cell if IW is on remote rank or convert 
! iw and iwl to local domain if IW is on myrank.

      if (iparallel == 1) then
         if (itabg_w(iw)%irank /= myrank) cycle

         iw  = itabg_w(iw)%iw_myrank
         iwl = itabg_wl(iwl)%iwl_myrank
      endif

     arf_land = landflux(ilf)%arf_land
     land_area(iwl) = land_area(iwl) + arf_land * arw0(iw)
  enddo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  do iwl = 1,mwl
     write(12,*)land_area(iwl)
  enddo

  close(12)

!======================
! W-point lats and lons
!======================

  if(iparallel == 0)then
     fname = trim(hfilepref)//'-atm-grid-info.txt'
  else
     fname = trim(hfilepref)//'-atm-grid-info_'//trim(post)//'.txt'
  endif

  open(12,file=trim(fname),status='replace',form='formatted')

  do iw = 2,mwa

     im = itab_w(iw)%im1
     glatm1 = asin(zem(im) / erad) * piu180
     glonm1 = atan2(yem(im),xem(im)) * piu180
     im = itab_w(iw)%im2
     glatm2 = asin(zem(im) / erad) * piu180
     glonm2 = atan2(yem(im),xem(im)) * piu180
     im = itab_w(iw)%im3
     glatm3 = asin(zem(im) / erad) * piu180
     glonm3 = atan2(yem(im),xem(im)) * piu180
     glon_mean = atan2(yem(im)/3.0,xem(im)/3.0) * piu180

     lon_diff1 = abs(glon_mean)
     lon_diff2 = abs(abs(glon_mean) - 180.0)

     ! Convert lons from -180:180 to 0:360
     if(glonm1 < 0.0)glonm1 = 360.0 + glonm1
     if(glonm2 < 0.0)glonm2 = 360.0 + glonm2
     if(glonm3 < 0.0)glonm3 = 360.0 + glonm3

     min_lon = min(glonm1, min(glonm2, glonm3))
     max_lon = max(glonm1, max(glonm2, glonm3))

     if(max_lon <= 180.0)then
        write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,0
     elseif(min_lon >= 180.0)then
        glonm1 = glonm1 - 360.0
        glonm2 = glonm2 - 360.0
        glonm3 = glonm3 - 360.0
        write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,0
     else
        if(lon_diff1 < lon_diff2)then
           if(glonm1 > 180.0)glonm1 = glonm1 - 360.0
           if(glonm2 > 180.0)glonm2 = glonm2 - 360.0
           if(glonm3 > 180.0)glonm3 = glonm3 - 360.0
           write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,0
        else
           write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,0
           glonm1 = glonm1 - 360.0
           glonm2 = glonm2 - 360.0
           glonm3 = glonm3 - 360.0
           write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,1
        endif
     endif
  enddo

  close(12)

!==============
! Wind info
!==============

  if(iparallel == 0)then
     fname = trim(hfilepref)//'-wind-grid-info.txt'
  else
     fname = trim(hfilepref)//'-wind-grid-info_'//trim(post)//'.txt'
  endif

  open(12,file=trim(fname),status='replace',form='formatted')

  do iw = 2,mwa

     write(12,'(3(i0,x),12e16.6)')  &
          itab_w(iw)%iu1,itab_w(iw)%iu2,itab_w(iw)%iu3,  &
          itab_w(iw)%vxu1,itab_w(iw)%vyu1,itab_w(iw)%vzu1,  &
          itab_w(iw)%vxu2,itab_w(iw)%vyu2,itab_w(iw)%vzu2,  &
          itab_w(iw)%vxu3,itab_w(iw)%vyu3,itab_w(iw)%vzu3,  &
          xew(iw),yew(iw),zew(iw)

  enddo

  close(12)

!===========================
! LAND info
!===========================

  if(iparallel == 0)then
     fname = trim(hfilepref)//'-land-grid-info.txt'
  else
     fname = trim(hfilepref)//'-land-grid-info_'//trim(post)//'.txt'
  endif

  open(12,file=trim(fname),status='replace',form='formatted')

  do iwl = 2,mwl

     glatm1 =     asin(land%zeml(itab_wl(iwl)%im(1)) / erad) * piu180
     glonm1 =    atan2(land%yeml(itab_wl(iwl)%im(1))  &
                      ,land%xeml(itab_wl(iwl)%im(1))) * piu180
     glatm2 =     asin(land%zeml(itab_wl(iwl)%im(2)) / erad) * piu180
     glonm2 =    atan2(land%yeml(itab_wl(iwl)%im(2))  &
                      ,land%xeml(itab_wl(iwl)%im(2))) * piu180
     glatm3 =     asin(land%zeml(itab_wl(iwl)%im(3)) / erad) * piu180
     glonm3 =    atan2(land%yeml(itab_wl(iwl)%im(3))  &
                      ,land%xeml(itab_wl(iwl)%im(3))) * piu180
     glon_mean = atan2(land%yeml(itab_wl(iwl)%im(3))/3.0  &
                      ,land%xeml(itab_wl(iwl)%im(3))/3.0) * piu180

     lon_diff1 = abs(glon_mean)
     lon_diff2 = abs(abs(glon_mean) - 180.0)

     ! Convert lons from -180:180 to 0:360
     if(glonm1 < 0.0)glonm1 = 360.0 + glonm1
     if(glonm2 < 0.0)glonm2 = 360.0 + glonm2
     if(glonm3 < 0.0)glonm3 = 360.0 + glonm3

     min_lon = min(glonm1, min(glonm2, glonm3))
     max_lon = max(glonm1, max(glonm2, glonm3))

     if(max_lon <= 180.0)then
        write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,0
     elseif(min_lon >= 180.0)then
        glonm1 = glonm1 - 360.0
        glonm2 = glonm2 - 360.0
        glonm3 = glonm3 - 360.0
        write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,0
     else
        if(lon_diff1 < lon_diff2)then
           if(glonm1 > 180.0)glonm1 = glonm1 - 360.0
           if(glonm2 > 180.0)glonm2 = glonm2 - 360.0
           if(glonm3 > 180.0)glonm3 = glonm3 - 360.0
           write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,0
        else
           write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,0
           glonm1 = glonm1 - 360.0
           glonm2 = glonm2 - 360.0
           glonm3 = glonm3 - 360.0
           write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,1
        endif
     endif
  enddo

  close(12)

  if(iparallel == 0)then
     fname = trim(hfilepref)//'-sea-grid-info.txt'
  else
     fname = trim(hfilepref)//'-sea-grid-info_'//trim(post)//'.txt'
  endif

  open(12,file=trim(fname),status='replace',form='formatted')

  do iws = 2,mws

     glatm1 =     asin(sea%zems(itab_ws(iws)%im(1)) / erad) * piu180
     glonm1 =    atan2(sea%yems(itab_ws(iws)%im(1))  &
                      ,sea%xems(itab_ws(iws)%im(1))) * piu180
     glatm2 =     asin(sea%zems(itab_ws(iws)%im(2)) / erad) * piu180
     glonm2 =    atan2(sea%yems(itab_ws(iws)%im(2))  &
                      ,sea%xems(itab_ws(iws)%im(2))) * piu180
     glatm3 =     asin(sea%zems(itab_ws(iws)%im(3)) / erad) * piu180
     glonm3 =    atan2(sea%yems(itab_ws(iws)%im(3))  &
                      ,sea%xems(itab_ws(iws)%im(3))) * piu180
     glon_mean = atan2(sea%yems(itab_ws(iws)%im(3))/3.0  &
                      ,sea%xems(itab_ws(iws)%im(3))/3.0) * piu180

     lon_diff1 = abs(glon_mean)
     lon_diff2 = abs(abs(glon_mean) - 180.0)

     ! Convert lons from -180:180 to 0:360
     if(glonm1 < 0.0)glonm1 = 360.0 + glonm1
     if(glonm2 < 0.0)glonm2 = 360.0 + glonm2
     if(glonm3 < 0.0)glonm3 = 360.0 + glonm3

     min_lon = min(glonm1, min(glonm2, glonm3))
     max_lon = max(glonm1, max(glonm2, glonm3))

     if(max_lon <= 180.0)then
        write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,0
     elseif(min_lon >= 180.0)then
        glonm1 = glonm1 - 360.0
        glonm2 = glonm2 - 360.0
        glonm3 = glonm3 - 360.0
        write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,0
     else
        if(lon_diff1 < lon_diff2)then
           if(glonm1 > 180.0)glonm1 = glonm1 - 360.0
           if(glonm2 > 180.0)glonm2 = glonm2 - 360.0
           if(glonm3 > 180.0)glonm3 = glonm3 - 360.0
           write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,0
        else
           write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,0
           glonm1 = glonm1 - 360.0
           glonm2 = glonm2 - 360.0
           glonm3 = glonm3 - 360.0
           write(12,'(6f8.2,i2)')glatm1,glatm2,glatm3,glonm1,glonm2,glonm3,1
        endif
     endif
  enddo

  close(12)

  return
end subroutine write_land_grid_info

!==================================================================

subroutine normalize_ed_output_vars()

  use misc_coms, only: time8, dtlm, frqstate
  use ed_options, only: frq_phenology
  use mem_leaf, only: land
  use ed_structure_defs
  use leaf_coms, only: dt_leaf

  implicit none

  real :: tfact1
  real :: tfact2
  type(site), pointer :: cs
  type(patch), pointer :: cp
  type(cohort), pointer :: cc
  real :: storage_decay

  tfact1 = dt_leaf / frqstate

  ! Averaged over the state printout time step

  cs => land%first_site
  do while(associated(cs))
     
     ! Initialize the land output arrays.
     land%gpp(cs%iland) = 0.0
     land%rh(cs%iland)  = 0.0
     land%nep(cs%iland) = 0.0
     land%veg_lai(cs%iland) = 0.0

     cs%omean_precip = cs%omean_precip / frqstate
     cs%omean_qprecip = cs%omean_qprecip / frqstate
     cs%omean_netrad = cs%omean_netrad * tfact1

     cp => cs%oldest_patch
     do while(associated(cp))
        
        ! Normalize to obtain [umol/m2/s]
        cp%omean_rh = cp%omean_rh * tfact1
        
        cp%omean_runoff = cp%omean_runoff * tfact1
        cp%omean_wflux = cp%omean_wflux * tfact1
        cp%omean_latflux = cp%omean_latflux * tfact1
        cp%omean_qrunoff = cp%omean_qrunoff * tfact1
        cp%omean_hflux = cp%omean_hflux * tfact1

        storage_decay = 0.0
        cp%omean_nep = 0.0
        
        cs%mmean_rh = cs%mmean_rh + cp%omean_rh * cp%area
        
        cc => cp%tallest
        do while(associated(cc))
           
           ! Normalize to obtain [umol/m2/s]
           cc%omean_gpp = cc%omean_gpp * tfact1
           cc%omean_leaf_resp = cc%omean_leaf_resp * tfact1
           cc%omean_root_resp = cc%omean_root_resp * tfact1
           
           ! Calculate growth, storage and vleaf respiration
           !  need to convert units from kgC/plant/day to umolC/m2/s.
           storage_decay = (cc%growth_respiration +   &
                cc%storage_respiration   &
                + cc%vleaf_respiration) & 
                * cc%nplant / (frq_phenology * 1.2e-8)
           
           ! Update plant contribution to NEP
           cp%omean_nep = cp%omean_nep + cc%omean_gpp -   &
                cc%omean_leaf_resp - cc%omean_root_resp - storage_decay
           
           cs%mmean_gpp = cs%mmean_gpp + cc%omean_gpp * cp%area
           cs%mmean_plresp = cs%mmean_plresp + (cc%omean_leaf_resp +   &
                cc%omean_root_resp + storage_decay) * cp%area
           
           ! Fill the land%.
           land%gpp(cs%iland) = land%gpp(cs%iland) + cp%area * cc%omean_gpp
           land%veg_lai(cs%iland) = land%veg_lai(cs%iland) + cp%area * cc%lai

           cc => cc%shorter
        enddo
        
        ! Update heterotrophic component to NEP
        cp%omean_nep = cp%omean_nep - cp%omean_rh

        ! Fill the land%.
        
        land%rh(cs%iland) = land%rh(cs%iland) + cp%area * cp%omean_rh
        land%nep(cs%iland) = land%nep(cs%iland) + cp%area * cp%omean_nep
        
        cp => cp%younger
     enddo
     cs => cs%next_site
  enddo

  return
end subroutine normalize_ed_output_vars

!=======================================================================

subroutine write_ed_hdf5(avg_flag)

  use var_tables, only: num_ED, vtab_ED
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
  use misc_coms, only: io6,hfilepref, time8, iyear1, imonth1, idate1, itime1,   &
       iclobber
  use leaf_coms, only: mwl

  implicit none

  integer :: nvtot
  integer :: nv
  character(len=128) :: hnamel
  character(len=128) :: hnamelh
  integer, save :: ioaunt = 10
  integer :: lenl
  logical :: exans
  integer :: nvcnt
  character(len=32) :: varn
  integer :: ndims
  integer, dimension(2) :: idims
  character(len=4), intent(in) :: avg_flag
  character(len=256) :: type
  
  if(avg_flag == 'hist')then
     type = 'ED'
  elseif(avg_flag == 'mavg')then
     type = 'ED-mm'
  elseif(avg_flag == 'yavg')then
     type = 'ED-ym'
  endif

  ! Construct h5 file name and open the file
  
  call makefnam8(hnamel,hfilepref,time8,iyear1,imonth1,idate1,  &
       itime1*100,trim(type),'$','h5')  ! '$' argument suppresses grid# encoding
  inquire(file=hnamel,exist=exans)
  if (exans .and. iclobber == 0) then
     write(io6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*)'!!!   trying to open file name :'
     write(io6,*)'!!!       ',hnamel
     write(io6,*)'!!!   but it already exists. run is ended.'
     write(io6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'write_ed'
  endif

!  print*,'+++++++++++++++++++++++++++++++++++++'
!  print*,'write_ed:open file:',trim(hnamel)
!  print*,'+++++++++++++++++++++++++++++++++++++'
  call shdf5_open(hnamel,'W',iclobber)
  ndims = 1
  idims(1) = 1
  call shdf5_orec(ndims,idims,'mwl',ivars=mwl)
!  Loop through the main variable table and write those variables
!     with the correct flag set

  do nv = 1,num_ED
     
     if (   &
          (vtab_ED(nv)%ihist == 1 .and. avg_flag == 'hist') .or.  &
          (vtab_ED(nv)%mavg == 1  .and. avg_flag == 'mavg') .or.  &
          (vtab_ED(nv)%yavg == 1  .and. avg_flag == 'yavg')       &
          )then
        
        varn = vtab_ED(nv)%name

        ndims      = vtab_ED(nv)%ndims
        idims(1:2) = vtab_ED(nv)%idims(1:2)
        
        ! Write the header information out to the file.
        
        if     (associated(vtab_ED(nv)%ivar1_p)) then
           call shdf5_orec(ndims,idims,trim(varn),ivara=vtab_ED(nv)%ivar1_p)
        elseif (associated(vtab_ED(nv)%rvar1_p)) then
           call shdf5_orec(ndims,idims,trim(varn),rvara=vtab_ED(nv)%rvar1_p)
        elseif (associated(vtab_ED(nv)%ivar2_p)) then
           call shdf5_orec(ndims,idims,trim(varn),ivara=vtab_ED(nv)%ivar2_p)
        elseif (associated(vtab_ED(nv)%rvar2_p)) then
           call shdf5_orec(ndims,idims,trim(varn),rvara=vtab_ED(nv)%rvar2_p)
        elseif (associated(vtab_ED(nv)%dvar2_p)) then
           call shdf5_orec(ndims,idims,trim(varn),dvara=vtab_ED(nv)%dvar2_p)
        endif
        
!        write(io6,*)'wrote  :',nvcnt,varn
        
     endif
     
  enddo

  call shdf5_close()

  return
end subroutine write_ed_hdf5

!=======================================================================

subroutine zero_ed_output_vars()

  use mem_leaf, only: land
  use ed_structure_defs

  implicit none

  type(site), pointer :: cs
  type(patch), pointer :: cp
  type(cohort), pointer :: cc

  cs => land%first_site
  do while(associated(cs))

     cs%omean_precip = 0.0
     cs%omean_qprecip = 0.0
     cs%omean_netrad = 0.0

     cp => cs%oldest_patch
     do while(associated(cp))

        cp%omean_rh = 0.0
        cp%omean_runoff = 0.0
        cp%omean_wflux = 0.0
        cp%omean_latflux = 0.0
        cp%omean_qrunoff = 0.0
        cp%omean_hflux = 0.0

        cc => cp%tallest
        do while(associated(cc))
           
           cc%omean_gpp = 0.0
           cc%omean_leaf_resp = 0.0
           cc%omean_root_resp = 0.0

           cc => cc%shorter
        enddo

        cp => cp%younger
     enddo
     cs => cs%next_site
  enddo

  return
end subroutine zero_ed_output_vars

!==============================================================================
subroutine normalize_ed_daily_vars(cs)
  
  use ed_options, only: frq_phenology
  use ed_structure_defs
  use leaf_coms, only: dt_leaf

  implicit none
  
  type(site), target :: cs
  type(patch), pointer :: cp
  type(cohort), pointer :: cc
  real :: tfact2


  tfact2 = dt_leaf / frq_phenology

  cp => cs%oldest_patch
  do while(associated(cp))

     cp%dmean_A_decomp = cp%dmean_A_decomp * tfact2
     cp%dmean_Af_decomp = cp%dmean_Af_decomp * tfact2
     
     cc => cp%tallest
     do while(associated(cc))
        
        ! Cohort variables
        cc%dmean_gpp = cc%dmean_gpp * tfact2
        cc%dmean_gpp_pot = cc%dmean_gpp_pot * tfact2
        cc%dmean_gpp_max = cc%dmean_gpp_max * tfact2
        cc%dmean_leaf_resp = cc%dmean_leaf_resp * tfact2
        cc%dmean_root_resp = cc%dmean_root_resp * tfact2

        cc => cc%shorter
     enddo
     
     cp => cp%younger
  enddo

  return
end subroutine normalize_ed_daily_vars

!==============================================================================

subroutine reinitialize_ed_daily_vars(cs)
  
  use ed_structure_defs

  implicit none
  
  type(site), target :: cs
  type(patch), pointer :: cp
  type(cohort), pointer :: cc

  cp => cs%oldest_patch
  do while(associated(cp))

     cp%dmean_A_decomp = 0.0
     cp%dmean_Af_decomp = 0.0
     
     cc => cp%tallest
     do while(associated(cc))
        
        ! Cohort variables
        cc%dmean_gpp = 0.0
        cc%dmean_gpp_pot = 0.0
        cc%dmean_gpp_max = 0.0
        cc%dmean_leaf_resp = 0.0
        cc%dmean_root_resp = 0.0
        
        cc => cc%shorter
     enddo
     
     cp => cp%younger
  enddo

  return
end subroutine reinitialize_ed_daily_vars

!=================================================================

subroutine normalize_ed_monthly_vars()

  use misc_coms, only: frqstate
  use mem_leaf, only: land
  use ed_structure_defs

  implicit none
  
  type(site), pointer :: cs

  ! Change units on carbon fluxes to umol/m2/month/frqstate to tC/ha/month.

  cs => land%first_site
  do while(associated(cs))
     cs%mmean_gpp = cs%mmean_gpp * frqstate * 1.2e-7
     cs%mmean_plresp = cs%mmean_plresp * frqstate * 1.2e-7
     cs%mmean_rh = cs%mmean_rh * frqstate * 1.2e-7
     cs%mmean_nep = cs%mmean_gpp - cs%mmean_plresp - cs%mmean_rh
     land%gpp(cs%iland) = cs%mmean_gpp
     land%rh(cs%iland) = cs%mmean_rh
     land%nep(cs%iland) = cs%mmean_nep
     cs => cs%next_site
  enddo

  return
end subroutine normalize_ed_monthly_vars

!=================================================================

subroutine update_ed_yearly_vars()

  use ed_structure_defs
  use mem_leaf, only: land
  
  implicit none

  type(site), pointer :: cs
  type(patch), pointer :: cp
  type(cohort), pointer :: cc
  real, external :: ed_agb

  ! For all agb's, change units from kgC/m2/y to tC/ha/y.
  cs => land%first_site
  do while(associated(cs))
     land%agb(cs%iland) = sum(cs%agb(1:n_pft, 1:n_dbh)) * 10.0
     land%basal_area(cs%iland) = sum(cs%basal_area(1:n_pft, 1:n_dbh))
     land%agb_growth(cs%iland) = sum(cs%agb_growth(1:n_pft, 1:n_dbh)) * 10.0
     land%agb_mort(cs%iland) = sum(cs%agb_mort(1:n_pft, 1:n_dbh)) * 10.0
     land%agb_cut(cs%iland) = sum(cs%agb_cut(1:n_pft, 1:n_dbh)) * 10.0
     land%agb_recruit(cs%iland) = 0.0
     cp => cs%oldest_patch
     do while(associated(cp))
        cc => cp%tallest
        do while(associated(cc))
           if(cc%new_recruit_flag == 1)then
              land%agb_recruit(cs%iland) = land%agb_recruit(cs%iland) +   &
                   ed_agb(cc%bdead, cc%balive, cc%bleaf, cc%pft, cc%hite,   &
                   cc%bstorage) * cp%area * 10.0 * cc%nplant
              cc%new_recruit_flag = 0
           endif
           cc%first_census = 1
           cc => cc%shorter
        enddo
        cp => cp%younger
     enddo
     cs => cs%next_site
  enddo

  return
end subroutine update_ed_yearly_vars


!=================================================================

subroutine zero_ed_monthly_vars()

  use mem_leaf, only: land
  use ed_structure_defs

  implicit none
  
  type(site), pointer :: cs

  cs => land%first_site
  do while(associated(cs))
     cs%mmean_gpp = 0.0
     cs%mmean_plresp = 0.0
     cs%mmean_rh = 0.0
     cs%mmean_nep = 0.0
     cs => cs%next_site
  enddo

  return
end subroutine zero_ed_monthly_vars

!=================================================================

subroutine zero_ed_yearly_vars()

  use mem_leaf, only: land
  use ed_structure_defs

  implicit none
  
  type(site), pointer :: cs

  cs => land%first_site
  do while(associated(cs))
     cs%agb_growth(1:n_pft, 1:n_dbh) = 0.0
     cs%agb_mort(1:n_pft, 1:n_dbh) = 0.0
     cs%agb_cut(1:n_pft, 1:n_dbh) = 0.0
     cs%agb_recruit(1:n_pft, 1:n_dbh) = 0.0
     cs%basal_area_growth(1:n_pft, 1:n_dbh) = 0.0
     cs%basal_area_mort(1:n_pft, 1:n_dbh) = 0.0
     cs%basal_area_cut(1:n_pft, 1:n_dbh) = 0.0
     cs%basal_area_recruit(1:n_pft, 1:n_dbh) = 0.0
     cs => cs%next_site
  enddo

  return
end subroutine zero_ed_yearly_vars

!====================================================================
subroutine print_ed_budgets()

  use mem_leaf, only: land
  use ed_structure_defs
  use misc_coms, only: frqstate, current_time, dtlm
  use pft_coms, only: c2n_slow, c2n_structural, c2n_leaf, c2n_storage,   &
       c2n_stem, c2n_recruit, n_pft
  use decomposition_coms, only: f_labile
  use leaf_coms, only: nzg, slmsts
  
  implicit none

  type(site), pointer :: cs
  type(patch), pointer :: cp
  type(cohort), pointer :: cc

  integer, save :: ifirst=1
  real, save :: first_soil_carbon
  real, save :: first_plant_biomass
  real, save :: first_soil_nitrogen
  real, save :: first_plant_nitrogen

  real :: soil_carbon
  real :: plant_biomass
  real :: soil_nitrogen
  real :: plant_nitrogen
  real :: stored_seed_carbon
  real :: stored_seed_nitrogen

  real, save :: total_gpp=0.0
  real, save :: total_resp= 0.0
  real, save :: total_rh=0.0
  real, save :: total_maintenance=0.0

  integer :: ipft
  real, dimension(n_pft) :: lai
  real, dimension(n_pft) :: agb
  real, external :: ed_agb
  real, dimension(nzg) :: soil_water
  
  ! Compute the site average
  cs => land%first_site
  soil_carbon = 0.0
  soil_nitrogen = 0.0
  plant_biomass = 0.0
  plant_nitrogen = 0.0
  stored_seed_carbon = 0.0
  stored_seed_nitrogen = 0.0
  lai(1:n_pft) = 0.0
  agb(1:n_pft) = 0.0

  cp => cs%oldest_patch
  do while(associated(cp))
     soil_carbon = soil_carbon + cp%area * (cp%fast_soil_C +   &
          cp%slow_soil_C + cp%structural_soil_C)
     soil_nitrogen = soil_nitrogen + cp%area * (cp%fast_soil_N +   &
          cp%slow_soil_C / c2n_slow + cp%structural_soil_C /   &
          c2n_structural + cp%mineralized_soil_N)
     soil_water(nzg) = soil_water(nzg) + cp%area * cp%soil_water(nzg) /   &
          slmsts(cp%ntext_soil(nzg))

     cc => cp%tallest
     do while(associated(cc))
        plant_biomass = plant_biomass + cp%area * cc%nplant * (cc%balive +  &
             cc%bdead + cc%bstorage)
        total_gpp = total_gpp + cp%area * cc%omean_gpp * frqstate * 1.2e-8
        total_resp = total_resp + cp%area * cc%omean_gpp * frqstate * 1.2e-8
        if(current_time%time < dtlm(1))then
           total_maintenance = total_maintenance + cp%area *   &
                cc%maintenance_costs * cc%nplant
        endif
        plant_nitrogen = plant_nitrogen + cp%area * cc%nplant * (  &
             f_labile(cc%pft) * cc%balive / c2n_leaf(cc%pft) +   &
             cc%bstorage / c2n_storage +   &
             ((1.0 - f_labile(cc%pft)) * cc%balive + cc%bdead) / c2n_stem)
        lai(cc%pft) = lai(cc%pft) + cp%area * cc%lai
        agb(cc%pft) = agb(cc%pft) + cp%area * ed_agb(cc%bdead, cc%balive,   &
             cc%bleaf, cc%pft, cc%hite, cc%bstorage) * cc%nplant

        cc => cc%shorter
     enddo
     total_resp = total_resp - cp%area * cp%omean_nep * frqstate * 1.2e-8
     total_rh = total_rh + cp%area * cp%omean_rh * frqstate * 1.2e-8
     do ipft = 1,n_pft
        stored_seed_carbon = stored_seed_carbon + cp%area * cp%repro(ipft)
        stored_seed_nitrogen = stored_seed_nitrogen + cp%area *   &
             cp%repro(ipft) / c2n_recruit(ipft)
     enddo
     cp => cp%younger
  enddo

  if(ifirst == 1)then
     ifirst = 0
     open(12,file='ed_budgets.pft1-4.txt',form='formatted',status='replace')
     first_plant_biomass = plant_biomass
     first_plant_nitrogen = plant_nitrogen
     first_soil_carbon = soil_carbon
     first_soil_nitrogen = soil_nitrogen
  else
     open(12,file='ed_budgets.pft1-4.txt',form='formatted',status='old',  &
          position='append')
  endif

  if(current_time%time < dtlm(1))write(12,'(19e14.5)')  &
       soil_carbon,   &
       plant_biomass,  &
       total_gpp,  &
       total_resp,  &
       total_rh,  &
       stored_seed_carbon,  &
       soil_nitrogen, &
       plant_nitrogen,  &
       stored_seed_nitrogen, &
       total_maintenance, &
       soil_water(nzg),  &
       lai(1),  &
       agb(1), &
       lai(2),  &
       agb(2), &
       lai(3),  &
       agb(3), &
       lai(4),  &
       agb(4)


  close(12)

  return
end subroutine print_ed_budgets
!====================================================================
subroutine print_ed_soi()

  use mem_leaf, only: land
  use ed_structure_defs
  use misc_coms, only: frqstate, current_time, dtlm,io6
  use pft_coms, only: c2n_slow, c2n_structural, c2n_leaf, c2n_storage,   &
       c2n_stem, c2n_recruit, n_pft
  use decomposition_coms, only: f_labile
  use leaf_coms, only: nzg, slmsts, dslz
  use misc_coms, only: hfilepref
  use consts_coms, only: cliq, cice, alli
  
  implicit none

  type(site), pointer :: cs
  type(patch), pointer :: cp
  type(cohort), pointer :: cc

  integer, save :: ifirst=1

  real :: soil_carbon
  real :: soil_nitrogen
  real :: soil_water_fraction
  real :: patch_soil_water
  real :: patch_soil_water_max
  integer :: k

  real, dimension(n_pft) :: plant_biomass
  real, dimension(n_pft) :: lai
  real, dimension(n_pft) :: agb
  real, external :: ed_agb
  real :: surface_water
  real :: sfcwater_energy
  real :: site_soil_water
  real :: site_soil_energy
  real :: veg_water
  real :: canopy_water
  real :: runoff
  real :: wflux
  real :: latflux
  real :: qrunoff
  real :: hflux
  real :: canair_heat
  real :: veg_heat
  real :: site_gpp
  real :: site_nep
  real :: site_rh
  real :: plant_nitrogen

  ! Compute the site average
  cs => land%first_site
  soil_carbon = 0.0
  soil_nitrogen = 0.0
  soil_water_fraction = 0.0
  lai(1:n_pft) = 0.0
  agb(1:n_pft) = 0.0
  plant_biomass(1:n_pft) = 0.0
  surface_water = 0.0
  sfcwater_energy = 0.0
  site_soil_water = 0.0
  site_soil_energy = 0.0
  veg_water = 0.0
  canopy_water = 0.0
  wflux = 0.0
  latflux = 0.0
  runoff = 0.0
  hflux = 0.0
  qrunoff = 0.0
  canair_heat = 0.0
  veg_heat = 0.0
  site_gpp = 0.0
  site_rh = 0.0
  site_nep = 0.0
  plant_nitrogen = 0.0

  cp => cs%oldest_patch
  do while(associated(cp))
     soil_carbon = soil_carbon + cp%area * (cp%fast_soil_C +   &
          cp%slow_soil_C + cp%structural_soil_C)
     soil_nitrogen = soil_nitrogen + cp%area * (cp%fast_soil_N +   &
          cp%slow_soil_C / c2n_slow + cp%structural_soil_C /   &
          c2n_structural + cp%mineralized_soil_N)
     patch_soil_water = 0.0
     patch_soil_water_max = 0.0
     site_nep = site_nep + cp%area * cp%omean_nep
     site_rh = site_rh + cp%area * cp%omean_rh

     do k=1,nzg
        patch_soil_water = patch_soil_water + cp%soil_water(k) *   &
             dslz(k)
        patch_soil_water_max = patch_soil_water_max +   &
             slmsts(cp%ntext_soil(k)) * dslz(k)
        site_soil_water = site_soil_water + cp%area * cp%soil_water(k) *   &
             dslz(k)
        site_soil_energy = site_soil_energy + cp%area * cp%soil_energy(k) *   &
             dslz(k)
     enddo
     soil_water_fraction = soil_water_fraction + cp%area *   &
          patch_soil_water / patch_soil_water_max

     canair_heat = canair_heat + (cp%can_temp-273.15) * cp%area * 1004.0 *land%rhos(cs%iland) * cp%can_depth

     do k = 1,cp%nlev_sfcwater
        surface_water = surface_water + cp%sfcwater_mass(k) * cp%area
        sfcwater_energy = sfcwater_energy + cp%sfcwater_energy(k) * cp%sfcwater_mass(k) * cp%area
     enddo

     canopy_water = canopy_water + cp%area * cp%can_shv *   &
          land%rhos(cs%iland) * cp%can_depth

     wflux = wflux + cp%area * cp%omean_wflux
     latflux = latflux + cp%area * cp%omean_latflux
     runoff = runoff + cp%area * cp%omean_runoff
     hflux = hflux + cp%area * cp%omean_hflux
     qrunoff = qrunoff + cp%area * cp%omean_qrunoff

     cc => cp%tallest
     do while(associated(cc))
        plant_biomass(cc%pft) = plant_biomass(cc%pft) + cp%area *   &
             cc%nplant * (cc%balive + cc%bdead + cc%bstorage)
        lai(cc%pft) = lai(cc%pft) + cp%area * cc%lai
        agb(cc%pft) = agb(cc%pft) + cp%area * ed_agb(cc%bdead, cc%balive,   &
             cc%bleaf, cc%pft, cc%hite, cc%bstorage) * cc%nplant
        veg_water = veg_water + cp%area * cc%veg_water
        if(cc%veg_temp > 273.15)then
!           veg_heat = veg_heat + (cc%hcapveg * (cc%veg_temp - 273.15))* cp%area
           veg_heat = veg_heat + (cc%hcapveg * (cc%veg_temp - 273.15) +  &
                cc%veg_water * (cliq * (cc%veg_temp - 273.15) + alli) ) *   &
                cp%area
        else
!           veg_heat = veg_heat + (cc%hcapveg * (cc%veg_temp - 273.15) ) * cp%area
           veg_heat = veg_heat + (cc%hcapveg * (cc%veg_temp - 273.15) +  &
                cc%veg_water * cice * (cc%veg_temp - 273.15) ) * cp%area
        endif
        site_gpp = site_gpp + cp%area * cc%omean_gpp
        plant_nitrogen = plant_nitrogen + cp%area * cc%nplant *   &
             (cc%balive * f_labile(cc%pft) / c2n_leaf(cc%pft) +   &
             (cc%balive * (1.0 - f_labile(cc%pft)) + cc%bdead) / c2n_stem +  &
             cc%bstorage / c2n_storage)
        cc => cc%shorter
     enddo

     cp => cp%younger
  enddo

  if(ifirst == 1)then
     ifirst = 0
     open(12, file=trim(hfilepref)//'-ed_budgets.txt', form='formatted',  &
          status='replace')
     write(12,'(35a)')'soil.c ','soil.n ','soil.w ','agb1 ','agb2 ',  &
          'agb3 ','agb4 ','b1 ','b2 ','b3 ','b4 ','lai1 ','lai2 ','lai3 ',  &
          'lai4 ','surface.water ','precip ','soil.water ','can.water ',  &
          'veg.water ','wflux ','runoff ','soil.energy ','sfcwater.energy ', &
          'canair.heat ','veg.heat ','netrad ','q.precip ','q.runoff ',  &
          'hflux ','latflux ','gpp ','rh ','nep ','plant.n '
  else
     open(12, file=trim(hfilepref)//'-ed_budgets.txt', form='formatted',  &
          status='old', position='append')
  endif

  if(current_time%time < dtlm(1))  &
       write(12,'(35e14.5)')  &
       soil_carbon,   &
       soil_nitrogen,  &
       soil_water_fraction,  &
       agb(1:4), &
       plant_biomass(1:4),  &
       lai(1:4),  &
       surface_water,  &
       cs%omean_precip,  &
       site_soil_water, &
       canopy_water,  &
       veg_water,  &
       wflux,  &
       runoff, &
       site_soil_energy, &
       sfcwater_energy, &
       canair_heat,  &
       veg_heat,  &
       cs%omean_netrad, &
       cs%omean_qprecip,  &
       qrunoff,  &
       hflux,  &
       latflux,  &
       site_gpp,  & !umol/m2/s
       site_rh,  & ! umol/m2/s
       site_nep, &  ! umol/m2/s
       plant_nitrogen ! kgN/m2

  if( (agb(1)+1.0 == 0.0 .and. agb(1) - 1.0 == 0.0)  .or. &
       (agb(2)+1.0 == 0.0 .and. agb(2) - 1.0 == 0.0)  .or. &
       (agb(3)+1.0 == 0.0 .and. agb(3) - 1.0 == 0.0)  .or. &
       (agb(4)+1.0 == 0.0 .and. agb(4) - 1.0 == 0.0) )then
     write(io6,*)'bad agb'
     cp => cs%oldest_patch
     do while(associated(cp))
        write(io6,*)'Patch:',cp%dist_type,cp%age,cp%area
        cc => cp%tallest
        do while(associated(cc))
           write(io6,*)'Cohort:',cc%pft,cc%hite,cc%nplant,cc%balive,cc%bdead,cc%bstorage
           cc => cc%shorter
        enddo
        cp => cp%younger
     enddo
     stop
  endif

  close(12)

  return
end subroutine print_ed_soi
!===============================================================
subroutine write_lite_files_2d()
  
  use var_tables, only: num_var, vtab_r
  use misc_coms,  only: io6, ioutput, hfilepref, time8, iyear1, imonth1, idate1,  &
       itime1, iclobber, iparallel
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
  use mem_basic,  only: uc,wc,sh_w,sh_v,theta,press
  use mem_radiate, only: rshort,rlong,albedt,rlongup,rshort_top,rshortup_top,rlongup_top, aero_opt_depth
  use mem_turb, only: sflux_t,sflux_r,ustar
  use mem_grid, only: mwa, lpw, mua, lpu
  use mem_para, only: myrank

  implicit none

  ! This routine writes the chosen variables on the history file.

  character(len=10) :: post
  character(len=128) :: hnamel,hnamelh
  character(len=32) :: varn
  logical exans
  integer, save :: ioaunt=10,nvtot
  integer :: nv,nvcnt,lenl
  
  integer :: ndims,idims(2)
  integer :: i,j
  real, dimension(mua) :: surf_u
  real, dimension(mwa) :: surf_var
  integer :: iu,iw
  
  if (ioutput == 0) return

  ! Construct h5 file name and open the file

  if(iparallel == 0) then
     post = '$'
  else
     write(post,'(i10)')myrank
     post = 'r'//trim(adjustl(post))
  endif

  call makefnam8(hnamel,hfilepref,time8,iyear1,imonth1,idate1,  &
       itime1*100,'S',post,'h5')  ! '$' argument suppresses grid# encoding

  lenl = len_trim(hnamel)
  
  inquire(file=hnamel,exist=exans)
  if (exans .and. iclobber == 0) then
     write(io6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*)'!!!   trying to open file name :'
     write(io6,*)'!!!       ',hnamel
     write(io6,*)'!!!   but it already exists. run is ended.'
     write(io6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'history_write'
  endif
  
  write(io6,*)'+++++++++++++++++++++++++++++++++++++'
  write(io6,*)'history_write:open file:',trim(hnamel)
  write(io6,*)'+++++++++++++++++++++++++++++++++++++'
  call shdf5_open(hnamel,'W',iclobber)
  
  !  Loop through the main variable table and write those variables
  !     with the correct flag set
  
  varn = 'UC_LP'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mua
  idims(2) = mua
  do iu=2,mua
     surf_u(iu) = uc(lpu(iu),iu)
  enddo
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=surf_u)

  varn = 'WC_LP'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  do iw=2,mwa
     surf_var(iw) = wc(lpw(iw),iw)
  enddo
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=surf_var)

  varn = 'SH_W_LP'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  do iw=2,mwa
     surf_var(iw) = sh_w(lpw(iw),iw)
  enddo
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=surf_var)

  varn = 'SH_V_LP'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  do iw=2,mwa
     surf_var(iw) = sh_v(lpw(iw),iw)
  enddo
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=surf_var)

  varn = 'THETA_LP'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  do iw=2,mwa
     surf_var(iw) = theta(lpw(iw),iw)
  enddo
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=surf_var)

  varn = 'PRESS_LP'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  do iw=2,mwa
     surf_var(iw) = sngl(press(lpw(iw),iw))
  enddo
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=surf_var)

  varn = 'RSHORT'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=rshort)

  varn = 'RSHORT_TOP'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=rshort_top)

  varn = 'RSHORTUP_TOP'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=rshortup_top)

  varn = 'RLONGUP_TOP'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=rlongup_top)

  varn = 'RLONG'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=rlong)

  varn = 'ALBEDT'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=albedt)

  varn = 'RLONGUP'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=rlongup)

  varn = 'SFLUX_T'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=sflux_t)

  varn = 'SFLUX_R'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=sflux_r)

  varn = 'USTAR'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=ustar)

  if(allocated(aero_opt_depth))then
     varn = 'AERO_T'
     write(io6,*)'writing:',varn
     ndims      = 1
     idims(1) = mwa
     idims(2) = mwa
     write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
     call shdf5_orec(ndims,idims,trim(varn),rvara=aero_opt_depth)
  endif

  call shdf5_close()
  
  print 12,time8,hnamelh
12 format(/,1X,79('*'),/, '  History write         ','  Time8 = ',F11.0  &
        ,/,'      Header file name - ',A60,/,1X,79('*'))
  
  return
end subroutine write_lite_files_2d

!===============================================================
subroutine write_lite_files_3d()
  
  use var_tables, only: num_var, vtab_r
  use misc_coms,  only: io6,ioutput, hfilepref, time8, iyear1, imonth1, idate1,  &
       itime1, iclobber,iparallel
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
  use mem_basic,  only: uc,wc,sh_w,sh_v,theta,press,rho
  use mem_radiate, only: rshort,rlong,albedt,rlongup
  use mem_turb, only: sflux_t,sflux_r,ustar
  use mem_grid, only: nza, mwa, mua
  use mem_micro, only: sh_c,accpa,accpg,accph,accpp,accpr,accps
  use mem_cuparm, only: aconpr
  use leaf_coms, only: mwl, nzg
  use mem_leaf, only: land
  use mem_para, only: myrank

  implicit none

  ! This routine writes the chosen variables on the history file.


  character(len=10) :: post
  character(len=128) :: hnamel,hnamelh
  character(len=32) :: varn
  logical exans
  integer, save :: ioaunt=10,nvtot
  integer :: nv,nvcnt,lenl
  
  integer :: ndims,idims(2)
  integer :: i,j
  real, dimension(nza,mwa) :: press_s
  real, dimension(nza,mwa) :: rho_s
  
  if (ioutput == 0) return

  ! Construct h5 file name and open the file

  if(mod(time8,86400.d0) == 0.0d0)return

  if(iparallel == 0) then
     post = '$'
  else
     write(post,'(i10)')myrank
     post = 'r'//trim(adjustl(post))
  endif

  call makefnam8(hnamel,hfilepref,time8,iyear1,imonth1,idate1,  &
       itime1*100,'R',post,'h5')  ! '$' argument suppresses grid# encoding

  lenl = len_trim(hnamel)
  
  inquire(file=hnamel,exist=exans)
  if (exans .and. iclobber == 0) then
     write(io6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*)'!!!   trying to open file name :'
     write(io6,*)'!!!       ',hnamel
     write(io6,*)'!!!   but it already exists. run is ended.'
     write(io6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'history_write'
  endif
  
  write(io6,*)'+++++++++++++++++++++++++++++++++++++'
  write(io6,*)'history_write:open file:',trim(hnamel)
  write(io6,*)'+++++++++++++++++++++++++++++++++++++'
  call shdf5_open(hnamel,'W',iclobber)
  
  !  Loop through the main variable table and write those variables
  !     with the correct flag set
  
  varn = 'ACCPA'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1)
  call shdf5_orec(ndims,idims,trim(varn),rvara=accpa)

  varn = 'ACCPG'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1)
  call shdf5_orec(ndims,idims,trim(varn),rvara=accpg)

  varn = 'ACCPH'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1)
  call shdf5_orec(ndims,idims,trim(varn),rvara=accph)

  varn = 'ACCPP'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1)
  call shdf5_orec(ndims,idims,trim(varn),rvara=accpp)

  varn = 'ACCPR'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1)
  call shdf5_orec(ndims,idims,trim(varn),rvara=accpr)

  varn = 'ACCPS'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1)
  call shdf5_orec(ndims,idims,trim(varn),rvara=accps)

  varn = 'ACONPR'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1)
  call shdf5_orec(ndims,idims,trim(varn),rvara=aconpr)

  varn = 'LAND%CAN_SHV'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwl
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1)
  call shdf5_orec(ndims,idims,trim(varn),rvara=land%can_shv)

  varn = 'LAND%CAN_TEMP'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwl
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1)
  call shdf5_orec(ndims,idims,trim(varn),rvara=land%can_temp)

  varn = 'LAND%STOM_RESIST'
  write(io6,*)'writing:',varn
  ndims      = 1
  idims(1) = mwl
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1)
  call shdf5_orec(ndims,idims,trim(varn),rvara=land%stom_resist)

  varn = 'LAND%SOIL_ENERGY'
  write(io6,*)'writing:',varn
  ndims      = 2
  idims(1) = nzg
  idims(2) = mwl
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1)
  call shdf5_orec(ndims,idims,trim(varn),rvara=land%soil_energy)

  varn = 'LAND%SOIL_WATER'
  write(io6,*)'writing:',varn
  ndims      = 2
  idims(1) = nzg
  idims(2) = mwl
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1)
  call shdf5_orec(ndims,idims,trim(varn),rvara=land%soil_water)

  varn = 'UC'
  write(io6,*)'writing:',varn
  ndims      = 2
  idims(1) = nza
  idims(2) = mua
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=uc)

  varn = 'WC'
  write(io6,*)'writing:',varn
  ndims      = 2
  idims(1) = nza
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=wc)

  varn = 'SH_W'
  write(io6,*)'writing:',varn
  ndims      = 2
  idims(1) = nza
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=sh_w)

  varn = 'SH_V'
  write(io6,*)'writing:',varn
  ndims      = 2
  idims(1) = nza
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=sh_v)

  varn = 'SH_C'
  write(io6,*)'writing:',varn
  ndims      = 2
  idims(1) = nza
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=sh_c)

  varn = 'THETA'
  write(io6,*)'writing:',varn
  ndims      = 2
  idims(1) = nza
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  call shdf5_orec(ndims,idims,trim(varn),rvara=theta)

  varn = 'SPRESS'
  write(io6,*)'writing:',varn
  ndims      = 2
  idims(1) = nza
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  do i=1,idims(1)
     do j=1,idims(2)
        press_s(i,j) = sngl(press(i,j))
     enddo
  enddo
  call shdf5_orec(ndims,idims,trim(varn),rvara=press_s)

  varn = 'SRHO'
  write(io6,*)'writing:',varn
  ndims      = 2
  idims(1) = nza
  idims(2) = mwa
  write (ioaunt,'(a32,1x,i3,2i8)') varn,ndims,idims(1:2)
  do i=1,idims(1)
     do j=1,idims(2)
        rho_s(i,j) = sngl(rho(i,j))
     enddo
  enddo
  call shdf5_orec(ndims,idims,trim(varn),rvara=rho_s)

  call shdf5_close()
  
  print 12,time8,hnamelh
12 format(/,1X,79('*'),/, '  History write         ','  Time8 = ',F11.0  &
        ,/,'      Header file name - ',A60,/,1X,79('*'))
  
  return
end subroutine write_lite_files_3d
