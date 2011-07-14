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
subroutine ndvi_database_read(iaction)

use mem_leaf,  only: land, itab_wl

use leaf_coms, only: mml, mwl, iupdndvi,   &
                     fnames_ndvi, ctotdate_ndvi, s1900_ndvi,  &
                     indvifile, ndvi_database

use misc_coms, only: io6, iyear1, imonth1, idate1, itime1, timmax8,  &
                     time8, runtype, s1900_init, s1900_sim

use leaf_db,   only: leaf_database_read
use mem_para,  only: myrank
use rastro_evts

implicit none

integer, intent(in) :: iaction

integer :: indviy,indvim,indvid,indvih
integer :: iyears, imonths, idates, ihours

character(len=14) :: totdate_start
character(len=14) :: totdate_init
character(len=14) :: totdate
character(len=14) :: totdatem
character(len=14) :: ctotdate_current
character(len=1)  :: dummy
character(len=128) :: flnm
character(len=80) :: line, line2, filemonth

integer :: nf
integer :: iwl
integer :: ntimes, jtime

real(kind=8) :: secs_init
real(kind=8) :: secs1, secs2

real :: datp(mwl)

integer, save :: indvicyclic, nndvifiles

logical :: there

#ifdef OLAM_RASTRO
call rst_event_i_f(OLAM_NDVI_DATABASE_READ_IN,iaction)
#endif


! Check type of call to ndvi_database_read

if (iaction == 0) then

! Convert current model time from s1900 to years, months, dates, hours

   call date_secs_ymdt(s1900_sim,iyears,imonths,idates,ihours)

! Initialize ndvi cyclic flag to zero

   indvicyclic = 0

   flnm = trim(ndvi_database)//'HEADER'

   write(io6,*)  'Checking ndvi database header file ',trim(flnm)

   inquire(file=flnm,exist=there)

   if (.not. there) then
      write(io6,*) 'NDVI database file header was not found - stopping run'
      stop 'stop: no ndvi file'
   endif

! Open & read NDVI_DATABASE HEADER file; make inventory of file names & times
! Assumes data in header file is chronologically ordered oldest to newest.

   call rams_f_open(29,flnm,'FORMATTED','OLD','READ',0)

   read(29,*) dummy
   read(29,*) ntimes

   nndvifiles = ntimes

   do jtime = 1,ntimes
      read(29,'(A80)') line                       ! Read filetime info
      call char_strip_var(line,filemonth,line2)   ! Extract 3-letter file month
      read (line2,*) indviy,indvim,indvid,indvih  ! Extract year, month, day, hr
      
! If file year is read as zero, ndvi data is expected to be cyclic over 1 year.
! Increment indvicyclic to indicate this and use current simulation year for
! ndvi database file times.

      if (indviy == 0) then
         indvicyclic = indvicyclic + 1
         indviy = iyears
      endif

! Compute and store ndvi database file times in arrays

      fnames_ndvi(jtime) = trim(ndvi_database)//trim(filemonth)
      call date_make_big (indviy,indvim,indvid,indvih*100,ctotdate_ndvi(jtime))
      call date_abs_secs2(indviy,indvim,indvid,indvih*100,s1900_ndvi(jtime))
   enddo                     

   close(29)

! If ndvi cyclic flag > 0, check its value against ntimes and stop if they 
! are unequal.  If they are equal, reset ndvi cyclic flag to 1 and augment
! ndvi file arrays by 1 at each end.

   if (indvicyclic > 0) then
      if (indvicyclic /= ntimes) then
         write(io6,'(/,a)') 'Some but not all ndvi database files do not have'
         write(io6,'(a)')   'year 0000, which is ambiguous.  Stopping model.'
         stop 'stop_ndvi_inv'
      endif
      indvicyclic = 1
      nndvifiles = ntimes + 2

! Shift ndvi data file names and times by one array element

      do jtime = ntimes,1,-1
         fnames_ndvi  (jtime+1) = fnames_ndvi  (jtime)
         ctotdate_ndvi(jtime+1) = ctotdate_ndvi(jtime)
         s1900_ndvi   (jtime+1) = s1900_ndvi   (jtime)
      enddo      

! Add new ndvi member at beginning of time sequence

      fnames_ndvi(1) = fnames_ndvi(ntimes+1)

      call date_unmake_big(indviy,indvim,indvid,indvih,ctotdate_ndvi(ntimes+1))
      call date_make_big(indviy-1,indvim,indvid,indvih,ctotdate_ndvi(1))
      call date_abs_secs2(indviy-1,indvim,indvid,indvih,s1900_ndvi(1))

! Add new ndvi member at end of time sequence

      fnames_ndvi(ntimes+2) = fnames_ndvi(2)

      call date_unmake_big(indviy,indvim,indvid,indvih,ctotdate_ndvi(2))
      call date_make_big(indviy+1,indvim,indvid,indvih,ctotdate_ndvi(ntimes+2))
      call date_abs_secs2(indviy+1,indvim,indvid,indvih,s1900_ndvi(ntimes+2))

   endif

! Loop over number of NDVI_DATABASE file times and search for the one that
! corresponds to current or most recent model time.

   indvifile = 0
   do nf = 1,nndvifiles

      write(io6,*) 'nndvif0 ',nf,s1900_ndvi(nf),' ',s1900_sim

      if (s1900_ndvi(nf) <= s1900_sim) then
         indvifile = nf
      endif
   enddo

   if (indvifile < 1) then
      write(io6,*) ' '
      write(io6,*) 'Unable to find previous or current ndvi file for current'
      write(io6,*) 'model time.  Stopping model '
      stop 'stop: no current ndvi file'
   endif

elseif (iaction == 1) then
    
! Processing next ndvi file (only called with iaction = 1 if iupdndvi = 1)

   indvifile = indvifile + 1
   
   if (indvifile > nndvifiles)then
      if(indvicyclic == 0) then
         write(io6,*) ' '
         write(io6,*) 'No future ndvi file is available for nudging '
         write(io6,*) 'Stopping model '
         stop 'stop: no future ndvi file'
      else
         indvifile = 3
         do jtime = 1, nndvifiles
            call date_unmake_big(indviy,indvim,indvid,indvih,ctotdate_ndvi(jtime))
            call date_make_big(indviy+1,indvim,indvid,indvih,ctotdate_ndvi(jtime))
            call date_abs_secs2(indviy+1,indvim,indvid,indvih,s1900_ndvi(jtime))
         enddo
      endif
   endif
      
   land%veg_ndvip(:) = land%veg_ndvif(:)

endif

! Read and interpolate current ndvi file
  
call leaf_database_read(mwl,                 &
                        mml,                 &
                        itab_wl,             &
                        land%xeml,           &
                        land%yeml,           &
                        land%zeml,           &
                        trim(ndvi_database), &
                        trim(fnames_ndvi(indvifile)), &
                        'ndvi',              &
                        datp=datp            )

! Loop over all land cells (already defined and filled with leaf_class)

do iwl = 2,mwl
   land%veg_ndvif(iwl) = max(.05,datp(iwl))
enddo

if (iaction == 0) then
   land%veg_ndvip(:) = land%veg_ndvif(:)
endif

#ifdef OLAM_RASTRO
call rst_event_i_f(OLAM_NDVI_DATABASE_READ_OUT,iaction)
#endif

return
end subroutine ndvi_database_read

