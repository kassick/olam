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
Module mem_addsc
   
   ! Added scalar variables and tendencies

   Type addsc_vars   
      real, allocatable :: sclp(:,:)
      real, allocatable :: sclt(:,:)
      real, allocatable :: drydep(:)
   End Type
   
   type (addsc_vars), allocatable :: addsc(:)
  
Contains

!===============================================================================

   subroutine alloc_addsc(mza,mwa,naddsc)

   use misc_coms, only: io6

   implicit none

   integer, intent(in) :: mza,mwa,naddsc
   
   integer :: iaddsc

! Allocate arrays based on options (if necessary).

   do iaddsc = 1,naddsc
   
      write(io6,*) 'alloc_addsc ',iaddsc,mza,mwa,naddsc
   
      allocate (addsc(iaddsc)%sclp(mza,mwa))
      allocate (addsc(iaddsc)%drydep(mwa))

      addsc(iaddsc)%sclp(1:mza,1:mwa) = 0.
      addsc(iaddsc)%drydep(1:mwa) = 0.

   enddo

   return
   end subroutine alloc_addsc

!===============================================================================

   subroutine dealloc_addsc(naddsc)

   implicit none

   integer, intent(in) :: naddsc
   integer :: iaddsc

!  Deallocate arrays

   do iaddsc = 1,naddsc
      if (allocated(addsc(iaddsc)%sclp))   deallocate (addsc(iaddsc)%sclp)
      if (allocated(addsc(iaddsc)%drydep)) deallocate (addsc(iaddsc)%drydep)
   enddo
           
   return
   end subroutine dealloc_addsc

!===============================================================================
               
   subroutine filltab_addsc(mza,mwa,naddsc)

   use var_tables, only: vtables
 
   implicit none

   integer, intent(in) :: mza,mwa,naddsc

   integer :: iaddsc
   character (len=7) :: sname

   integer :: ndims
   integer :: idims(2)

   do iaddsc = 1,naddsc
      if (allocated(addsc(iaddsc)%sclp)) then

         ndims = 2
         idims(1) = mza
         idims(2) = mwa

         write(sname,'(a4,i3.3)') 'SCLP',iaddsc
         call vtables(ndims,idims,sname//' :hist:mpt1'  &
            ,rvara2=addsc(iaddsc)%sclp)

         ndims = 1
         idims(1) = mwa

         write(sname,'(a4,i3.3)') 'SCDD',iaddsc
         call vtables(ndims,idims,sname//' :hist:mpt1'  &
            ,rvara1=addsc(iaddsc)%drydep)

      endif
   enddo
    
   return
   end subroutine filltab_addsc

End Module mem_addsc
