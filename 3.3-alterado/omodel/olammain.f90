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
program main

use misc_coms, only: io6, iparallel
use mem_para,  only: myrank, mgroupsize

implicit none

character(len=128) :: name_name = '/mnt/pvfs2/OLAMIN'
character(len=128) :: cargv
character(len=30)  :: io6file
character(len=30)  :: logfile
integer :: numarg
integer :: i,n
integer :: rastro_id
integer :: bad = 0
real, external :: walltime
real, external :: etime
real :: totaltime, elapsed(2)
integer :: time_array_0(8), time_array_1(8)
integer(kind=4) :: start_time, end_time
real :: omp_time1, omp_time2
real, external :: omp_get_wtime
character(len=80) :: hostname
integer hn_status


! Determine if this run is parallel, and determine myrank and mgroupsize

call olam_mpi_init()

rastro_id=10
call rst_init_f(myrank,rastro_id)

hn_status = hostnm(hostname)
call rst_event_s_f(OLAM_INIT,hostname)

iparallel = 0
if (mgroupsize > 1) iparallel = 1

! If run is sequential, default choice is to set io6 to standard output unit 6.

! io6 = 6

! if (iparallel == 1) then
! .and. myrank > 0) then

! If run is parallel, default choice is to attach output unit io6 to separate files

   io6 = 20
   
   write (io6file,'(i10)') myrank
   io6file = '/tmp/o.io6_r'//trim(adjustl(io6file))

! First, remove file in case it exists


   call system('rm -f '//trim(io6file)//char(0))

   open(io6,file=io6file,status='new',form='formatted')

! endif

write(io6,'(/,a,i6)') ' myrank     = ',myrank
write(io6,'(  a,i6)') ' mgroupsize = ',mgroupsize
write(io6,'(  a,i6)') ' iparallel  = ',iparallel

numarg = iargc()
write(io6,*) 'numarg:', numarg

numarg = numarg + 1

! Parse the command line arguments

i = 1

do while (i <= numarg)
   call getarg(i,cargv)
   ! write(io6,*) 'args: ',i,cargv

   if (cargv(1:1) == '-') then
      if (cargv(2:2) == 'f') then
         call getarg(i+1,name_name)
         if (len_trim(name_name) < 1) bad = bad + 1
         i = i + 2
      else
         ! write(io6,*) 'OLAM unknown option: ', cargv
         i = i + 1
      endif
   else
      ! write(io6,*) 'OLAM unknown option: ', cargv
      i = i + 1
   endif
enddo

if (bad > 0) then
   write(io6,*) 'OLAM usage: ''exec name'' '
   write(io6,*) '  [-f ''Namelist file''] '
   stop 'bad command line arguments'
endif

write(io6,*) 'OLAM input namelist file: ',trim(name_name)

! Initialize, execute, and end olam run

call olam_run(name_name)

! If this run is parallel, finalize MPI and close io6 file

call olam_mpi_finalize()
if (iparallel == 1) then
   close(io6)
endif

call rst_finalize_f()
stop 'olam_end'
end program main
