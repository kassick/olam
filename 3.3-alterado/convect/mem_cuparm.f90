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
Module mem_cuparm

   real, allocatable :: thsrc(:,:)
   real, allocatable :: rtsrc(:,:)
   real, allocatable :: thsrcsh(:,:)
   real, allocatable :: rtsrcsh(:,:)

   real, allocatable :: aconpr(:)
   real, allocatable :: conprr(:)

Contains

!===============================================================================

   subroutine alloc_cuparm(mza,mwa,nqparm,nqparm_sh)

   implicit none

   integer, intent(in) :: mza,mwa
   integer, intent(in) :: nqparm,nqparm_sh ! these are max values over all MRLs
   
! Allocate all cuparm arrays if either deep or shallow cuparm is activated
! on any MRL
      
   if (nqparm > 0 .or. nqparm_sh > 0 ) then
      allocate (thsrc(mza,mwa)) ; thsrc(1:mza,1:mwa) = 0.
      allocate (rtsrc(mza,mwa)) ; rtsrc(1:mza,1:mwa) = 0.
      allocate (aconpr(mwa))    ; aconpr(1:mwa)      = 0.
      allocate (conprr(mwa))    ; conprr(1:mwa)      = 0.

      allocate (thsrcsh(mza,mwa)) ; thsrcsh(1:mza,1:mwa) = 0.
      allocate (rtsrcsh(mza,mwa)) ; rtsrcsh(1:mza,1:mwa) = 0.
   endif

   return
   end subroutine alloc_cuparm

!===============================================================================

   subroutine dealloc_cuparm()

   implicit none

   if (allocated(thsrc))   deallocate (thsrc)
   if (allocated(rtsrc))   deallocate (rtsrc)
   if (allocated(thsrcsh)) deallocate (thsrcsh)
   if (allocated(rtsrcsh)) deallocate (rtsrcsh)
   if (allocated(aconpr))  deallocate (aconpr)
   if (allocated(conprr))  deallocate (conprr)

   return
   end subroutine dealloc_cuparm
   
!===============================================================================

   subroutine filltab_cuparm(mza,mwa)

   use var_tables, only: vtables

   implicit none

   integer, intent(in) :: mza,mwa

   integer :: ndims
   integer :: idims(2)

   ndims = 2
   idims(1) = mza
   idims(2) = mwa

   if (allocated(thsrc))   call vtables(ndims,idims,'THSRC  :hist',rvara2=thsrc)
   if (allocated(rtsrc))   call vtables(ndims,idims,'RTSRC  :hist',rvara2=rtsrc)

   if (allocated(thsrcsh)) call vtables(ndims,idims,'THSRCSH :hist',rvara2=thsrcsh)
   if (allocated(rtsrcsh)) call vtables(ndims,idims,'RTSRCSH :hist',rvara2=rtsrcsh)

   ndims = 1
   idims(1) = mwa

   if (allocated(aconpr))  call vtables(ndims,idims,'ACONPR :hist',rvara1=aconpr)
   if (allocated(conprr))  call vtables(ndims,idims,'CONPRR :hist',rvara1=conprr)

   return
   end subroutine filltab_cuparm

End Module mem_cuparm
