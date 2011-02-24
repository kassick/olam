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
Module mem_basic

   real, allocatable :: ump  (:,:) ! past horiz momentum [kg/(m^2 s)]
   real, allocatable :: umc  (:,:) ! current horiz momentum [kg/(m^2 s)]
   real, allocatable :: wmc  (:,:) ! current vert momentum [kg/(m^2 s)]
   real, allocatable :: uc   (:,:) ! current horiz velocity [m/s]
   real, allocatable :: wc   (:,:) ! current horiz velocity [m/s]
   real, allocatable :: sh_w (:,:) ! tot water spec dens [kg_wat/kg_air]
   real, allocatable :: sh_v (:,:) ! spec hum [kg_vap/kg_air]
   real, allocatable :: thil (:,:) ! ice-liquid pot temp [K]
   real, allocatable :: theta(:,:) ! pot temp [K]

   real(kind=8), allocatable :: press(:,:) ! air pressure [Pa]
   real(kind=8), allocatable :: rho  (:,:) ! total air density [kg/m^3]

   real, allocatable :: cputime(:)

Contains

!===============================================================================

   subroutine alloc_basic(mza,mua,mwa)

   implicit none

   integer, intent(in) :: mza,mua,mwa

! Allocate basic memory needed for 'INITIAL' or 'HISTORY' runs
! and initialize allocated arrays to zero

   allocate (ump   (mza,mua)) ; ump   (1:mza,1:mua) = 0.
   allocate (umc   (mza,mua)) ; umc   (1:mza,1:mua) = 0.
   allocate (uc    (mza,mua)) ; uc    (1:mza,1:mua) = 0.
   allocate (wmc   (mza,mwa)) ; wmc   (1:mza,1:mwa) = 0.
   allocate (wc    (mza,mwa)) ; wc    (1:mza,1:mwa) = 0.
   allocate (thil  (mza,mwa)) ; thil  (1:mza,1:mwa) = 0.
   allocate (theta (mza,mwa)) ; theta (1:mza,1:mwa) = 0.
   allocate (rho   (mza,mwa)) ; rho   (1:mza,1:mwa) = 0.
   allocate (press (mza,mwa)) ; press (1:mza,1:mwa) = 0.
   allocate (sh_w  (mza,mwa)) ; sh_w  (1:mza,1:mwa) = 0.
   allocate (sh_v  (mza,mwa)) ; sh_v  (1:mza,1:mwa) = 0.

   allocate (cputime(mwa))   ; cputime(1:mwa) = 0.

   return
   end subroutine alloc_basic

!===============================================================================

   subroutine dealloc_basic()

   implicit none

   if (allocated(ump))     deallocate (ump)
   if (allocated(umc))     deallocate (umc)
   if (allocated(wmc))     deallocate (wmc)
   if (allocated(uc))      deallocate (uc)
   if (allocated(wc))      deallocate (wc)
   if (allocated(rho))     deallocate (rho)
   if (allocated(sh_w))    deallocate (sh_w)
   if (allocated(sh_v))    deallocate (sh_v)
   if (allocated(press))   deallocate (press)
   if (allocated(thil))    deallocate (thil)
   if (allocated(theta))   deallocate (theta)
   if (allocated(cputime)) deallocate (cputime)

   return
   end subroutine dealloc_basic

!===============================================================================

   subroutine filltab_basic(mza,mua,mwa)

   use var_tables, only: vtables

   implicit none

   integer, intent(in) :: mza,mua,mwa

   integer :: ndims
   integer :: idims(2)

   ndims = 2
   idims(1) = mza
   idims(2) = mua

   if (allocated(ump))     call vtables(ndims,idims,'UMP   :hist'  &
         ,rvara2=ump)
   if (allocated(umc))     call vtables(ndims,idims,'UMC   :hist'  &
         ,rvara2=umc)
   if (allocated(uc))      call vtables(ndims,idims,'UC    :hist'  &
         ,rvara2=uc)

   idims(2) = mwa

   if (allocated(wmc))     call vtables(ndims,idims,'WMC   :hist'  &
         ,rvara2=wmc)
   if (allocated(wc))      call vtables(ndims,idims,'WC    :hist'  &
         ,rvara2=wc)
   if (allocated(sh_w))    call vtables(ndims,idims,'SH_W  :hist:mpt1'  &
         ,rvara2=sh_w)
   if (allocated(sh_v))    call vtables(ndims,idims,'SH_V  :hist:mpt1'  &
         ,rvara2=sh_v)
   if (allocated(thil))    call vtables(ndims,idims,'THIL  :hist'  &
         ,rvara2=thil)
   if (allocated(theta))   call vtables(ndims,idims,'THETA :hist:mpt1'  &
         ,rvara2=theta)
   if (allocated(rho))     call vtables(ndims,idims,'RHO   :hist'  &
         ,dvara2=rho)
   if (allocated(press))   call vtables(ndims,idims,'PRESS :hist'  &
         ,dvara2=press)

   ndims = 1
   idims(1) = mwa

   if (allocated(cputime)) call vtables(ndims,idims,'CPUTIME ',rvara1=cputime)

   return
   end subroutine filltab_basic

End Module mem_basic
