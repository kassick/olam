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
Module mem_radiate

   integer :: jday ! Julian day
   
   real :: solfac  ! solar-constant coefficient for variable Earth-Sun dist
   real :: sunx    ! x-component of unit vector pointing to sun [m]
   real :: suny    ! y-component of unit vector pointing to sun [m]
   real :: sunz    ! z-component of unit vector pointing to sun [m]

   real, allocatable :: fthrd   (:,:)
   real, allocatable :: fthrd_lw(:,:)
                          
   real, allocatable :: rshort        (:)
   real, allocatable :: rshort_top    (:)
   real, allocatable :: rshortup_top  (:)
   real, allocatable :: rlongup_top   (:)
   real, allocatable :: rlong         (:)
   real, allocatable :: rlongup       (:)
   real, allocatable :: albedt        (:)
   real, allocatable :: albedt_beam   (:)
   real, allocatable :: albedt_diffuse(:)
   real, allocatable :: cosz          (:)
   real, allocatable :: rlong_albedo  (:)
   real, allocatable :: rshort_diffuse(:)
   real, allocatable :: aero_opt_depth(:)

   integer, allocatable ::rad_region(:)

!  These are used for adding extra levels at the top with the Mclatchy soundings
   integer, parameter :: maxadd_rad = 10 ! max allowed # of added rad levels
   integer            :: nadd_rad        ! actual # of added radiation levels
   real               :: zmrad = 30.e3   ! top of radiation grid

Contains

!===============================================================================

   subroutine alloc_radiate(mza,mwa,ilwrtyp,iswrtyp)

   use misc_coms, only: io6
   
   implicit none

   integer, intent(in) :: mza
   integer, intent(in) :: mwa
   integer, intent(in) :: ilwrtyp
   integer, intent(in) :: iswrtyp

! Allocate arrays based on options (if necessary)
! Initialize arrays to zero
      
   if (ilwrtyp + iswrtyp > 0)  then
   
      write(io6,*) 'allocating rad ',mwa,mza
   
      allocate (fthrd   (mza,mwa)) ; fthrd   (1:mza,1:mwa) = 0.
      allocate (fthrd_lw(mza,mwa)) ; fthrd_lw(1:mza,1:mwa) = 0.
      allocate (rshort      (mwa)) ; rshort        (1:mwa) = 0.

      if (iswrtyp == 3 .or. ilwrtyp == 3) then
         allocate (rad_region(mwa))
                   rad_region(1:mwa) = 0
      endif

      if (iswrtyp == 4 .or. ilwrtyp == 4) then
         allocate (aero_opt_depth   (mwa))
                   aero_opt_depth (1:mwa) = 0.
      endif
      
      allocate (rshort_top    (mwa)) ; rshort_top     (1:mwa) = 0.
      allocate (rshortup_top  (mwa)) ; rshortup_top   (1:mwa) = 0.
      allocate (rlongup_top   (mwa)) ; rlongup_top    (1:mwa) = 0.
      allocate (rshort_diffuse(mwa)) ; rshort_diffuse (1:mwa) = 0.
      allocate (rlong         (mwa)) ; rlong          (1:mwa) = 0.
      allocate (rlongup       (mwa)) ; rlongup        (1:mwa) = 0.
      allocate (rlong_albedo  (mwa)) ; rlong_albedo   (1:mwa) = 0.
      allocate (albedt        (mwa)) ; albedt         (1:mwa) = 0.
      allocate (albedt_beam   (mwa)) ; albedt_beam    (1:mwa) = 0.
      allocate (albedt_diffuse(mwa)) ; albedt_diffuse (1:mwa) = 0.
      allocate (cosz          (mwa)) ; cosz           (1:mwa) = 0.
   endif

   return
   end subroutine alloc_radiate

!===============================================================================

   subroutine dealloc_radiate()

   implicit none

   if (allocated(fthrd))          deallocate (fthrd)
   if (allocated(fthrd_lw))       deallocate (fthrd_lw)
   if (allocated(rshort))         deallocate (rshort)
   if (allocated(aero_opt_depth)) deallocate (aero_opt_depth)
   if (allocated(rad_region))     deallocate (rad_region)
   if (allocated(rshort_top))     deallocate (rshort_top)
   if (allocated(rshortup_top))   deallocate (rshortup_top)
   if (allocated(rlongup_top))    deallocate (rlongup_top)
   if (allocated(rshort_diffuse)) deallocate (rshort_diffuse)
   if (allocated(rlong))          deallocate (rlong)
   if (allocated(rlongup))        deallocate (rlongup)
   if (allocated(rlong_albedo))   deallocate (rlong_albedo)
   if (allocated(albedt))         deallocate (albedt)
   if (allocated(albedt_beam))    deallocate (albedt_beam)
   if (allocated(albedt_diffuse)) deallocate (albedt_diffuse)
   if (allocated(cosz))           deallocate (cosz)

   return
   end subroutine dealloc_radiate

!===============================================================================

   subroutine filltab_radiate(mza,mwa,ilwrtyp,iswrtyp)

   use var_tables, only: vtables

   implicit none

   integer, intent(in) :: mza,mwa,ilwrtyp,iswrtyp
   integer :: ndims
   integer :: idims(2)

   ndims = 2
   idims(1) = mza
   idims(2) = mwa

   if (allocated(fthrd))    call vtables(ndims,idims,'FTHRD    :hist'  &
         ,rvara2=fthrd)
   if (allocated(fthrd_lw)) call vtables(ndims,idims,'FTHRD_LW :hist'  &
         ,rvara2=fthrd_lw)

   ndims = 1
   idims(1) = mwa

   if (allocated(rshort))   call vtables(ndims,idims,'RSHORT   :hist'  &
         ,rvara1=rshort)
   if (allocated(rlong))    call vtables(ndims,idims,'RLONG    :hist'  &
         ,rvara1=rlong)
   if (allocated(rlongup))  call vtables(ndims,idims,'RLONGUP  :hist'  &
         ,rvara1=rlongup)
   if (allocated(albedt))   call vtables(ndims,idims,'ALBEDT   :hist'  &
         ,rvara1=albedt)
   if (allocated(cosz))     call vtables(ndims,idims,'COSZ     :hist'  &
         ,rvara1=cosz)

   return
   end subroutine filltab_radiate

End Module mem_radiate
