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
subroutine seaice_database_read(iaction)

use mem_sea,     only: sea, itab_ws

use sea_coms,    only: mms, mws, iupdseaice, iseaicecyclic, nseaicefiles, &
                       fnames_seaice, ctotdate_seaice, s1900_seaice,      &
                       iseaicefile, seaice_database

use misc_coms,   only: io6, iyear1, imonth1, idate1, itime1, timmax8,  &
                       time8, runtype, s1900_init, s1900_sim

use consts_coms, only: erad, piu180
use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec, shdf5_info
use max_dims,    only: maxsstfiles

implicit none

integer, intent(in) :: iaction

integer :: iseaicey,iseaicem,iseaiced,iseaiceh
integer :: iyears, imonths, idates, ihours

real, allocatable :: dato(:,:)

integer :: im
integer :: j
integer :: nio, njo, nperdeg
integer :: io1, io2, jo1, jo2
integer :: nf
integer :: iws
integer :: ntimes, jtime
integer :: totdateseaice
integer :: totdate_init
integer :: slen
integer :: ndims, idims(2)

real :: wio1, wio2, wjo1, wjo2
real :: offpix
real :: glat, glon
real :: rio, rjo

character(len=128) :: flnm
character(len=10)  :: sdate

#ifdef OLAM_RASTRO
call rst_event_i_f(OLAM_SEAICE_DATABASE_READ_IN,iaction)
#endif


! This subroutine is simpler than topm_database because it assumes that 
! each seaice_database file covers the entire geographic area of the model.
! If this ever changes, this subroutine must be modified.

! Check type of call to seaice_database_read

if (iaction == 0) then

! Convert current model time from s1900 to years, months, dates, hours

   call date_secs_ymdt(s1900_sim,iyears,imonths,idates,ihours)

! Initialize seaice cyclic flag to zero

   iseaicecyclic = 0
   nseaicefiles  = 0

   flnm = trim(seaice_database)//'??????????.h5'
   write(io6,*) 'Checking for seaice database files'
   call OLAM_filelist(fnames_seaice, maxsstfiles, flnm, nseaicefiles)

   if (nseaicefiles < 1) then
      write(io6,*) 'SEAICE database files '//flnm//' were not found.'
      write(io6,*) 'Stopping run.'
      stop 'stop: no seaice file'
   endif
   
   ntimes = nseaicefiles

   do jtime=1,ntimes

      ! Assume SEAICE file names should always end with YYYYMMDDHH.h5, and
      ! use the file name to infer the file date and time

      flnm  = fnames_seaice(jtime)
      slen = len_trim(flnm)
      sdate = flnm(slen-12:slen-3)

      read(sdate,'(i4,i2,i2,i2)') iseaicey, iseaicem, iseaiced, iseaiceh
      
      ! If file year is read as zero, seaice data is expected to be cyclic
      ! over 1 year. Increment iseaicecyclic to indicate this and use current
      ! simulation year for seaice database file times.

      if (iseaicey == 0) then
         iseaicecyclic = iseaicecyclic + 1
         iseaicey = iyears
      endif

      call date_make_big (iseaicey,iseaicem,iseaiced,iseaiceh*100, &
           ctotdate_seaice(jtime))
      call date_abs_secs2(iseaicey,iseaicem,iseaiced,iseaiceh*100, &
           s1900_seaice(jtime))
   enddo

   ! Make sure files are sorted by date
   call dintsort28(nseaicefiles,ctotdate_seaice,fnames_seaice,s1900_seaice)

! If seaice cyclic flag > 0, check its value against ntimes and stop if they 
! are unequal.  If they are equal, reset seaice cyclic flag to 1 and augment
! seaice file arrays by 1 at each end.

   if (iseaicecyclic > 0) then
      if (iseaicecyclic /= ntimes) then
         write(io6,'(/,a)') 'Some but not all seaice database files do not have'
         write(io6,'(a)')   'year 0000, which is ambiguous.  Stopping model.'
         stop 'stop_seaice_inv'
      endif
      iseaicecyclic = 1
      nseaicefiles = ntimes + 2

! Shift seaice data file names and times by one array element

      do jtime = ntimes,1,-1
         fnames_seaice  (jtime+1) = fnames_seaice  (jtime)
         ctotdate_seaice(jtime+1) = ctotdate_seaice(jtime)
         s1900_seaice   (jtime+1) = s1900_seaice   (jtime)
      enddo      

! Add new seaice member at beginning of time sequence

      fnames_seaice(1) = fnames_seaice(ntimes+1)

      call date_unmake_big(iseaicey,iseaicem,iseaiced,iseaiceh, &
           ctotdate_seaice(ntimes+1))
      call date_make_big(iseaicey-1,iseaicem,iseaiced,iseaiceh, &
           ctotdate_seaice(1))
      call date_abs_secs2(iseaicey-1,iseaicem,iseaiced,iseaiceh, &
           s1900_seaice(1))

! Add new seaice member at end of time sequence

      fnames_seaice(ntimes+2) = fnames_seaice(2)

      call date_unmake_big(iseaicey,iseaicem,iseaiced,iseaiceh, &
           ctotdate_seaice(2))
      call date_make_big(iseaicey+1,iseaicem,iseaiced,iseaiceh, &
           ctotdate_seaice(ntimes+2))
      call date_abs_secs2(iseaicey+1,iseaicem,iseaiced,iseaiceh, &
           s1900_seaice(ntimes+2))

   endif

! Loop over number of SEAICE_DATABASE file times and search for the one that
! corresponds to current or most recent model time.

   iseaicefile = 0
   do nf = 1,nseaicefiles

      write(io6,*) 'nseaicef0 ',nf,s1900_seaice(nf),' ',s1900_sim

      if (s1900_seaice(nf) <= s1900_sim) then
         iseaicefile = nf
      endif
   enddo

   if (iseaicefile < 1) then
      write(io6,*) ' '
      write(io6,*) 'Unable to find previous or current seaice file for current'
      write(io6,*) 'model time.  Stopping model.'
      stop 'stop: no current seaice file'
   endif

elseif (iaction == 1) then

! Processing next iseaice file (only called with iaction = 1 if iupdseaice = 1)

   iseaicefile = iseaicefile + 1
   
   if (iseaicefile > nseaicefiles) then
      if(iseaicecyclic == 0)then
         write(io6,*) ' '
         write(io6,*) 'No future seaice file is available for nudging '
         write(io6,*) 'Stopping model '
         stop 'stop: no future seaice file'
      else
         iseaicefile = 3
         do jtime = 1, nseaicefiles
            call date_unmake_big(iseaicey,iseaicem,iseaiced,iseaiceh,ctotdate_seaice(jtime))
            call date_make_big(iseaicey+1,iseaicem,iseaiced,iseaiceh,ctotdate_seaice(jtime))
            call date_abs_secs2(iseaicey+1,iseaicem,iseaiced,iseaiceh,s1900_seaice(jtime))
         enddo
      endif
   endif
      
   sea%seaicep(:) = sea%seaicef(:)   

endif

! Open, read, and close seaice dataset file

write(io6,*) 'seaice_database_read2 ', iseaicefile, &
     trim(fnames_seaice(iseaicefile))

call shdf5_open(fnames_seaice(iseaicefile),'R')
call shdf5_info('ice',ndims,idims)
nio = idims(1)
njo = idims(2)
allocate(dato(nio,njo))
call shdf5_irec(ndims,idims,'ice',rvara=dato)
call shdf5_close()

offpix = 0.
nperdeg = nio/360
if (mod(nio,nperdeg) == 2) offpix = .5

! Fill seaice array

do iws = 2,mws

   glat = asin(sea%zews(iws) / erad) * piu180
   glon = atan2(sea%yews(iws),sea%xews(iws)) * piu180

   glon = max(-179.999,min(179.999,glon))

   rio = 1. + (glon + 180.) * nperdeg + offpix
   rjo = 1. + (glat +  90.) * nperdeg + offpix

   io1 = int(rio)
   jo1 = int(rjo)
         
   wio2 = rio - float(io1)
   wjo2 = rjo - float(jo1)
           
   wio1 = 1. - wio2
   wjo1 = 1. - wjo2

   io2 = min(nio,io1 + 1)
   jo2 = min(njo,jo1 + 1)
      
   sea%seaicef(iws) =  &
      wio1 * (wjo1 * dato(io1,jo1) + wjo2 * dato(io1,jo2)) + &
      wio2 * (wjo1 * dato(io2,jo1) + wjo2 * dato(io2,jo2))

enddo

deallocate(dato)

if (iaction == 0) then
   sea%seaicep(:) = sea%seaicef(:)
endif

#ifdef OLAM_RASTRO
call rst_event_i_f(OLAM_SEAICE_DATABASE_READ_OUT,iaction)
#endif

return
end subroutine seaice_database_read
