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
Module mem_turb

   real, allocatable :: tkep    (:,:)
   real, allocatable :: epsp    (:,:)
   real, allocatable :: hkm     (:,:)
   real, allocatable :: vkm     (:,:)
   real, allocatable :: vkh     (:,:)
   real, allocatable :: sxfer_tk(:,:)
   real, allocatable :: sxfer_rk(:,:)

   real, allocatable :: vkm_sfc (:)
   real, allocatable :: sflux_w (:)
   real, allocatable :: sflux_t (:)
   real, allocatable :: sflux_r (:)
   real, allocatable :: vels    (:)
   real, allocatable :: ustar   (:)

Contains

!===============================================================================

   subroutine alloc_turb(mza,mwa,nsw_max,idiffk)

   use misc_coms, only: io6

   implicit none
   
   integer, intent(in) :: mza,mwa,nsw_max,idiffk

! Allocate arrays based on options (if necessary)
! Initialize arrays to zero

! CHECK FOR IDIFFK BASED ON GRID 1 FOR NOW
      
   if (idiffk == 1 .or. idiffk == 4 .or. idiffk == 5 .or. idiffk == 6) then
      allocate (tkep(mza,mwa));  tkep(1:mza,1:mwa) = 0.
   endif
   
   if (idiffk == 6) then
      allocate (epsp(mza,mwa));  epsp(1:mza,1:mwa) = 0.
   endif
   
   write(io6,*) idiffk,allocated(epsp)

   allocate (hkm(mza,mwa));  hkm(1:mza,1:mwa) = 0.
   allocate (vkm(mza,mwa));  vkm(1:mza,1:mwa) = 0.
   allocate (vkh(mza,mwa));  vkh(1:mza,1:mwa) = 0.

   allocate (sxfer_tk(nsw_max,mwa));  sxfer_tk(1:nsw_max,1:mwa) = 0.
   allocate (sxfer_rk(nsw_max,mwa));  sxfer_rk(1:nsw_max,1:mwa) = 0.

   allocate (vkm_sfc(mwa));  vkm_sfc(1:mwa) = 0.
   allocate (sflux_w(mwa));  sflux_w(1:mwa) = 0.
   allocate (sflux_t(mwa));  sflux_t(1:mwa) = 0.
   allocate (sflux_r(mwa));  sflux_r(1:mwa) = 0.
   allocate (vels   (mwa));  vels   (1:mwa) = 0.
   allocate (ustar  (mwa));  ustar  (1:mwa) = 0.
                          
   return
   end subroutine alloc_turb

!===============================================================================

   subroutine dealloc_turb()

   implicit none
   
   if (allocated(tkep))    deallocate (tkep)
   if (allocated(epsp))    deallocate (epsp)
   if (allocated(hkm))     deallocate (hkm)
   if (allocated(vkm))     deallocate (vkm)
   if (allocated(vkh))     deallocate (vkh)
   if (allocated(vkm_sfc)) deallocate (vkm_sfc)
   if (allocated(sflux_t)) deallocate (sflux_t)
   if (allocated(sflux_r)) deallocate (sflux_r)
   if (allocated(vels))    deallocate (vels)
   if (allocated(ustar))   deallocate (ustar)

   return
   end subroutine dealloc_turb

!===============================================================================

   subroutine filltab_turb(mza,mwa,nsw_max,idiffk)

   use var_tables, only: vtables
   use misc_coms,  only: io6
   
   implicit none

   integer, intent(in) :: mza,mwa,nsw_max,idiffk
   
   integer :: ndims
   integer, dimension(2) :: idims

! Fill pointers to arrays into variable tables

   ndims = 2
   idims(1) = mza
   idims(2) = mwa

   if (allocated(tkep)) call vtables(ndims,idims,'TKEP :hist:mpt1',rvara2=tkep)
   if (allocated(epsp)) call vtables(ndims,idims,'EPSP :hist:mpt1',rvara2=epsp)
   if (allocated(hkm))  call vtables(ndims,idims,'HKM  :hist:'    ,rvara2=hkm)
   if (allocated(vkm))  call vtables(ndims,idims,'VKM  :hist:'    ,rvara2=vkm)
   if (allocated(vkh))  call vtables(ndims,idims,'VKH  :hist'     ,rvara2=vkh)

   idims(1) = nsw_max

   if (allocated(sxfer_tk)) call vtables(ndims,idims,'SXFER_TK :hist'  &
         ,rvara2=sxfer_tk)
   if (allocated(sxfer_rk)) call vtables(ndims,idims,'SXFER_RK :hist'  &
         ,rvara2=sxfer_rk)

   ndims = 1
   idims(1) = mwa

   if (allocated(vkm_sfc)) call vtables(ndims,idims,'VKM_SFC :hist'  &
         ,rvara1=vkm_sfc)
   if (allocated(sflux_t)) call vtables(ndims,idims,'SFLUX_T :hist'  &
         ,rvara1=sflux_t)
   if (allocated(sflux_r)) call vtables(ndims,idims,'SFLUX_R :hist'  &
         ,rvara1=sflux_r)
   if (allocated(ustar))   call vtables(ndims,idims,'USTAR   :hist'  &
         ,rvara1=ustar)

   return
   end subroutine filltab_turb

End Module mem_turb
