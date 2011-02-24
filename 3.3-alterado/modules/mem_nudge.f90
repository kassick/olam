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
Module mem_nudge

   integer, allocatable :: itab_nudp(:,:)

   real, allocatable :: xenudp(:)
   real, allocatable :: yenudp(:)
   real, allocatable :: zenudp(:)

   real, allocatable ::    rho_obsp(:,:)
   real, allocatable ::  theta_obsp(:,:)
   real, allocatable ::    shw_obsp(:,:)
   real, allocatable :: uzonal_obsp(:,:)
   real, allocatable :: umerid_obsp(:,:)

   real, allocatable ::    rho_obsf(:,:)
   real, allocatable ::  theta_obsf(:,:)
   real, allocatable ::    shw_obsf(:,:)
   real, allocatable :: uzonal_obsf(:,:)
   real, allocatable :: umerid_obsf(:,:)

   real, allocatable ::    rho_obs(:,:)
   real, allocatable ::  theta_obs(:,:)
   real, allocatable ::    shw_obs(:,:)
   real, allocatable :: uzonal_obs(:,:)
   real, allocatable :: umerid_obs(:,:)

   real, allocatable ::    rho_sim(:,:)
   real, allocatable ::  theta_sim(:,:)
   real, allocatable ::    shw_sim(:,:)
   real, allocatable :: uzonal_sim(:,:)
   real, allocatable :: umerid_sim(:,:)

   integer :: nudflag
   integer :: nudrad
   integer :: nnudfiles
   integer :: nnudfl
   integer :: mnudp

   real :: tnudcent
   real :: wt_nudge_uv
   real :: wt_nudge_th
   real :: wt_nudge_pi
   real :: wt_nudge_rt
   real :: wt_nudge_grid

Contains

   subroutine alloc_nudge1()

   use misc_coms, only: io6

   implicit none

! Allocate arrays based on options (if necessary)

write(io6,*) 'allocating nudge1 ',mnudp

   allocate (itab_nudp(mnudp,6)) ; itab_nudp(1:mnudp,1:6) = 1

   allocate (xenudp(mnudp)) ;  xenudp(1:mnudp) = 0.
   allocate (yenudp(mnudp)) ;  yenudp(1:mnudp) = 0.
   allocate (zenudp(mnudp)) ;  zenudp(1:mnudp) = 0.

   return
   end subroutine alloc_nudge1

!=========================================================================

   subroutine alloc_nudge2(mza)

   use misc_coms, only: io6

   implicit none
   
   integer, intent(in) :: mza

! Allocate arrays based on options (if necessary)

write(io6,*) 'allocating nudge2 ',mza,mnudp

   allocate (   rho_obsp(mza,mnudp)) ;     rho_obsp(1:mza,1:mnudp) = 0.
   allocate ( theta_obsp(mza,mnudp)) ;   theta_obsp(1:mza,1:mnudp) = 0.
   allocate (   shw_obsp(mza,mnudp)) ;     shw_obsp(1:mza,1:mnudp) = 0.
   allocate (uzonal_obsp(mza,mnudp)) ;  uzonal_obsp(1:mza,1:mnudp) = 0.
   allocate (umerid_obsp(mza,mnudp)) ;  umerid_obsp(1:mza,1:mnudp) = 0.

   allocate (   rho_obsf(mza,mnudp)) ;     rho_obsf(1:mza,1:mnudp) = 0.
   allocate ( theta_obsf(mza,mnudp)) ;   theta_obsf(1:mza,1:mnudp) = 0.
   allocate (   shw_obsf(mza,mnudp)) ;     shw_obsf(1:mza,1:mnudp) = 0.
   allocate (uzonal_obsf(mza,mnudp)) ;  uzonal_obsf(1:mza,1:mnudp) = 0.
   allocate (umerid_obsf(mza,mnudp)) ;  umerid_obsf(1:mza,1:mnudp) = 0.

   allocate (   rho_obs(mza,mnudp)) ;     rho_obs(1:mza,1:mnudp) = 0.
   allocate ( theta_obs(mza,mnudp)) ;   theta_obs(1:mza,1:mnudp) = 0.
   allocate (   shw_obs(mza,mnudp)) ;     shw_obs(1:mza,1:mnudp) = 0.
   allocate (uzonal_obs(mza,mnudp)) ;  uzonal_obs(1:mza,1:mnudp) = 0.
   allocate (umerid_obs(mza,mnudp)) ;  umerid_obs(1:mza,1:mnudp) = 0.

   allocate (   rho_sim(mza,mnudp)) ;     rho_sim(1:mza,1:mnudp) = 0.
   allocate ( theta_sim(mza,mnudp)) ;   theta_sim(1:mza,1:mnudp) = 0.
   allocate (   shw_sim(mza,mnudp)) ;     shw_sim(1:mza,1:mnudp) = 0.
   allocate (uzonal_sim(mza,mnudp)) ;  uzonal_sim(1:mza,1:mnudp) = 0.
   allocate (umerid_sim(mza,mnudp)) ;  umerid_sim(1:mza,1:mnudp) = 0.

   return
   end subroutine alloc_nudge2

!=========================================================================

   subroutine filltab_nudge(mza)

   use var_tables, only: vtables

   implicit none
   
   integer, intent(in) :: mza
   
   integer :: ndims
   integer :: idims(2)

! NUDGE VARIABLES

   ndims = 2
   idims(1) = mza
   idims(2) = mnudp

!------------------------------------------------------------------------
   if (allocated(rho_obsp))    call vtables(ndims,idims,'RHO_OBSP   :hist'  &
         ,rvara2=rho_obsp)
!------------------------------------------------------------------------
   if (allocated(theta_obsp))  call vtables(ndims,idims,'THETA_OBSP :hist'  &
         ,rvara2=theta_obsp)
!------------------------------------------------------------------------
   if (allocated(shw_obsp))    call vtables(ndims,idims,'SHW_OBSP   :hist'  &
         ,rvara2=shw_obsp)
!------------------------------------------------------------------------
   if (allocated(uzonal_obsp)) call vtables(ndims,idims,'UZONAL_OBSP:hist'  &
         ,rvara2=uzonal_obsp)
!------------------------------------------------------------------------
   if (allocated(umerid_obsp)) call vtables(ndims,idims,'UMERID_OBSP:hist'  &
         ,rvara2=umerid_obsp)
!------------------------------------------------------------------------
   if (allocated(rho_obsf))    call vtables(ndims,idims,'RHO_OBSF   :hist'  &
         ,rvara2=rho_obsf)
!------------------------------------------------------------------------
   if (allocated(theta_obsf))  call vtables(ndims,idims,'THETA_OBSF :hist'  &
         ,rvara2=theta_obsf)
!------------------------------------------------------------------------
   if (allocated(shw_obsf))    call vtables(ndims,idims,'SHW_OBSF   :hist'  &
         ,rvara2=shw_obsf)
!------------------------------------------------------------------------
   if (allocated(uzonal_obsf)) call vtables(ndims,idims,'UZONAL_OBSF:hist'  &
         ,rvara2=uzonal_obsf)
!------------------------------------------------------------------------
   if (allocated(umerid_obsf)) call vtables(ndims,idims,'UMERID_OBSF:hist'  &
         ,rvara2=umerid_obsf)
!------------------------------------------------------------------------

   return
   end subroutine filltab_nudge

End Module mem_nudge
