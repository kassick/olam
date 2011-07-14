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
subroutine sea_startup()

use sea_coms, only: mms, mus, mws, maxjms, isstflg, seatmp,  &
                    iupdsst, iseaiceflg, iupdseaice

use mem_sea,  only: sea, alloc_sea, filltab_sea

use hdf5_utils
use misc_coms, only: io6, runtype
use rastro_evts

implicit none

integer :: ierr
integer :: ndims, idims(1)
integer :: iws

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAM_SEA_STARTUP_IN,rst_buf)
#endif

! Subroutine SEA_STARTUP allocates some sea arrays and initializes sst

! THIS SUBROUTINE DOES NOT INITIALIZE canopy temperature and moisture
! values, which depend on atmospheric conditions.

!-------------------------------------------------------------------------------
! STEP 1: Call alloc_sea and filltab_sea (sea grid arrays already allocated)
!-------------------------------------------------------------------------------

call alloc_sea(mws)
call filltab_sea(mms,mus,mws)

!-------------------------------------------------------------------------------
! STEP 2a: Fill sst values
!-------------------------------------------------------------------------------

if (isstflg == 2) then

! Default initialization of SST

   do iws = 2,mws
      sea%seatp(iws) = seatmp
      sea%seatf(iws) = seatmp
   enddo

elseif (runtype /= 'PLOTONLY' .and. runtype /= 'PARCOMBINE') then

! Read SST database
! Not needed for a plotonly run

   write(io6,'(/,a)') 'calling sst_database_read(0)'

   call sst_database_read(0)

   write(io6,'(/,a)') 'calling sst_database_read(1)'

   call sst_database_read(1)

endif

!-------------------------------------------------------------------------------
! STEP 2b: Fill sea ice values
!-------------------------------------------------------------------------------

if (iseaiceflg == 2) then

! Default initialization of SEAICE

   do iws = 2,mws
      sea%seaicef(iws) = 0.
   enddo

elseif (runtype /= 'PLOTONLY' .and. runtype /= 'PARCOMBINE') then

! Read SEA ICE database
! Not needed for a plotonly run

   write(io6,'(/,a)') 'calling seaice_database_read(0)'

   call seaice_database_read(0)

   if (iupdseaice == 1) then

      write(io6,'(/,a)') 'calling seaice_database_read(1)'

      call seaice_database_read(1)

   endif
   
endif

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_SEA_STARTUP_OUT,rst_buf)
#endif


return
end subroutine sea_startup

