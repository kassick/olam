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
subroutine ed_vegetation_dynamics()

  use ed_structure_defs
  use mem_leaf, only: land
  use misc_coms, only: simtime, current_time, dtlm
  use leaf_coms, only: dt_leaf
  use ed_options, only: include_fire, frq_phenology
  use rastro_evts

  implicit none

  type(site), pointer :: cs
  real :: tfact1
  real :: tfact2
  integer :: doy
  integer :: julday

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAM_ED_VEGETATION_DYNAMICS_IN,rst_buf)
#endif


  ! find the day of year
  doy = julday(current_time%month, current_time%date, current_time%year)
  
  ! Time factor for averaging FRQSTATE
  tfact1 = dt_leaf / frq_phenology

  ! Time factor for averaging dailies 
  tfact2 = frq_phenology / (365.25 * 86400.0)

  cs => land%first_site
  do while(associated(cs))

     call normalize_ed_daily_vars(cs)

     call phenology_driver(cs, doy, current_time%month, tfact1)

     call dbalive_dt(cs, tfact2)

     if(current_time%date == 1)then

        call structural_growth(cs, current_time%month)
        call reproduction(cs, current_time%month)
        cs%min_monthly_temp = 500.0

        if(include_fire == 1)call fire_frequency(current_time%month, cs)
        call site_disturbance_rates(current_time%month, current_time%year, cs)

        if(current_time%month == 1)then

           call fuse_patches(cs)
           call apply_disturbances(cs)

        endif

     endif

     call update_C_and_N_pools(cs)

     call reinitialize_ed_daily_vars(cs)

     cs => cs%next_site
  enddo

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_ED_VEGETATION_DYNAMICS_OUT,rst_buf)
#endif

  return
end subroutine ed_vegetation_dynamics

!========================================

subroutine update_model_time(ctime,dtlong)
  use misc_coms, only: simtime
  use rastro_evts

  implicit none

  type(simtime) :: ctime
  real :: dtlong

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAM_UPDATE_MODEL_TIME_IN,rst_buf)
#endif

  ctime%time = ctime%time + dtlong
  if(ctime%time.ge.86400.0)then
     ctime%time = ctime%time - 86400.0
     ctime%date = ctime%date + 1
     if(ctime%date.eq.29 .and. ctime%month .eq. 2   &
          .and. mod(ctime%year,4).ne.0)then
        ctime%date = 1
        ctime%month = ctime%month + 1
     elseif(ctime%date.eq.30 .and. ctime%month.eq.2)then
        ctime%date = 1
        ctime%month = ctime%month + 1
     elseif(ctime%date.eq.31)then
        if(ctime%month.eq.4 .or. ctime%month.eq.6 .or. ctime%month.eq.9  &
             .or. ctime%month.eq.11)then
           ctime%date = 1
           ctime%month = ctime%month + 1
        endif
     elseif(ctime%date.eq.32)then
        ctime%date = 1
        ctime%month = ctime%month + 1
        if(ctime%month.eq.13)then
           ctime%month = 1
           ctime%year = ctime%year + 1
        endif
     endif
  elseif(ctime%time.lt.0.0)then
     ctime%time = ctime%time + 86400.0
     ctime%date = ctime%date - 1
     if(ctime%date.eq.0)then
        ctime%month = ctime%month - 1
        if(ctime%month.eq.0)then
           ctime%month = 12
           ctime%year = ctime%year - 1
           ctime%date = 31
        else
           if(ctime%month.eq.2)then
              if(mod(ctime%year,4).eq.0)then
                 ctime%date = 29
              else
                 ctime%date = 28
              endif
           elseif(ctime%month.eq.4 .or. ctime%month.eq.6 .or.   &
                ctime%month.eq.9 .or. ctime%month.eq.11)then
              ctime%date = 30
           else
              ctime%date = 31
           endif
        endif
     endif
  endif

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_UPDATE_MODEL_TIME_OUT,rst_buf)
#endif

  return
end subroutine update_model_time
