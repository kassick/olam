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
Module hdf5_utils

Contains

subroutine shdf5_open(locfn,access,idelete)

use misc_coms, only: io6

implicit none

character(len=*) :: locfn     ! file name
character(len=*) :: access    ! File access ('R','W','RW')
integer, optional :: idelete  ! If W, delete/overwrite file if exists? 1=yes, 0=no
                              ! Only needed when access='W'

integer :: hdferr ! Error flag
integer :: iaccess ! int access flag
character(len=2) :: caccess ! File access ('R ','W ','RW')

logical :: exists ! File existence

real, external :: walltime
real :: wtime_start_shdf5open, inicio

wtime_start_shdf5open = walltime(0.)

caccess = access

! Check for existence of RAMS file.

inicio = walltime(wtime_start_shdf5open)

inquire(file=trim(locfn),exist=exists)

write(io6, *) '                             ==T== Ver se o arquivo existia: ',(walltime(wtime_start_shdf5open)-inicio)

! Create a new file or open an existing RAMS file.
if (access(1:1) == 'R') then
   if (.not.exists) then
      print*,'shdf5_open:'
      print*,'   Attempt to open a file for reading that does not exist.'
      print*,'   Filename: ',trim(locfn)
      stop 'shdf5_open: no file'
   else
      if (caccess == 'R ') iaccess = 1
      if (caccess == 'RW') iaccess = 2

    inicio = walltime(wtime_start_shdf5open)

      call fh5f_open(trim(locfn)//char(0), iaccess, hdferr)

      
    write(io6, *) '                             ==T== Tempo total gasto na fh5f_open(): ',(walltime(wtime_start_shdf5open)-inicio)
      
      if (hdferr < 0) then
         print*,'shdf5_open:'
         print*,'   Error opening hdf5 file - error -',hdferr
         print*,'   Filename: ',trim(locfn)
         stop 'shdf5_open: open error'      
      endif
   endif
elseif (access(1:1) == 'W') then
   if (.not.exists) then
      iaccess=2
    
      inicio = walltime(wtime_start_shdf5open)
      
    call fh5f_create(trim(locfn)//char(0), iaccess, hdferr)
    
    write(io6, *) '                             ==T== Tempo gasto na fh5f_create(): ',(walltime(wtime_start_shdf5open)-inicio)

    else
      if(.not.present(idelete) ) then
         print*,'shdf5_open: idelete not specified when access=W'
         stop 'shdf5_open: no idelete'
      endif
      
      if(idelete == 0) then
         print*,'In shdf5_open:'
         print*,'   Attempt to open an existing file for writing, '
         print*,'      but overwrite is disabled. idelete=',idelete
         print*,'   Filename: ',trim(locfn)
         stop 'shdf5_open'
      else
         call system('rm -f '//trim(locfn)//char(0))
         iaccess=1
      
         inicio = walltime(wtime_start_shdf5open)
         
    call fh5f_create(trim(locfn)//char(0), iaccess, hdferr)
    
    write(io6, *) '                             ==T== Tempo gasto na fh5f_create(): ',(walltime(wtime_start_shdf5open)-inicio)
    
    endif
   endif
   if(hdferr < 0) then
      print*,'HDF5 file create failed:',hdferr
      print*,'file name:',trim(locfn),' ',trim(access), idelete
      stop 'shdf5_open: bad create'
   endif
endif


return
end subroutine shdf5_open

!===============================================================================

subroutine shdf5_info(dsetname,ndims,dims)

implicit none

character(len=*) :: dsetname ! Dataset name
integer :: dims(*)

integer :: ndims ! Dataset rank (in file)

integer :: hdferr ! Error flag


! Open the dataset.

call fh5d_open( trim(dsetname)//char(0), hdferr)

if (hdferr < 0) then
   print*, 'In shdf5_info:'
   print*, 'Variable ', trim(dsetname), ' is not in the currently opened hdf5 file'
   ndims   = 0
   dims(1) = 0
   return
endif

! Get dataset's dimensions

call fh5s_get_ndims(ndims)
call fh5s_get_dims(dims)

!print*,'ndims: ',ndims
!print*,'dims: ',dims(1:ndims)

call fh5d_close(hdferr)

return
end subroutine shdf5_info


!===============================================================================

subroutine shdf5_orec(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                         ,ivars,rvars,cvars,dvars,lvars)
use misc_coms, only: io6

implicit none

character(len=*) :: dsetname ! Variable label
integer :: ndims             ! Number of dimensions or rank
integer, dimension(*) :: dims ! Dataset dimensions.

! Array and scalar arguments for different types. Only specify one in each call.
integer,          optional :: ivara(*),ivars
real,             optional :: rvara(*),rvars
character(len=*), optional :: cvara(*),cvars
real(kind=8),     optional :: dvara(*),dvars
logical,          optional :: lvara(*),lvars

integer:: h5_type   ! Local type designator

integer, dimension(4) :: dimsh ! Dataset dimensions.

character(len=2) :: ctype    ! Variable type: int, real, char
integer :: hdferr ! Error flag

real, external :: walltime
real :: wtime_start_shdf5orec, inicio

wtime_start_shdf5orec = walltime(0.)

! Find which data type is input
    if(present(ivars)) then ; ctype='is'
elseif(present(rvars)) then ; ctype='rs'
elseif(present(cvars)) then ; ctype='cs'
elseif(present(dvars)) then ; ctype='ds'
elseif(present(lvars)) then ; ctype='ls'
elseif(present(ivara)) then ; ctype='ia'
elseif(present(rvara)) then ; ctype='ra'
elseif(present(cvara)) then ; ctype='ca'
elseif(present(dvara)) then ; ctype='da'
elseif(present(lvara)) then ; ctype='la'
else
   print*,'Incorrect or missing data field argument in shdf5_orec'
   stop 'shdf5_orec: bad data field'
endif

! Check dimensions and set compression chunk size

if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
   print*,'Dimension error in shdf5_orec:',ndims,dims(1:ndims)
   stop 'shdf5_orec: bad dims'
endif
dimsh(1:ndims) = dims(1:ndims)
     
inicio = walltime(wtime_start_shdf5orec)

! Prepare memory and options for the write
call fh5_prepare_write(ndims, dimsh, hdferr)
write(io6, *) '                             ==T== Tempo gasto na fh5_prepare_write(): ',(walltime(wtime_start_shdf5orec)-inicio)

if (hdferr /= 0) then
   print*,'shdf5_orec: can''t prepare requested field:',trim(dsetname)
   return
endif

if (ctype(1:1) == 'i') h5_type=1
if (ctype(1:1) == 'r') h5_type=2
if (ctype(1:1) == 'c') h5_type=3
if (ctype(1:1) == 'd') h5_type=4  ! If native precision is 8 bytes, do h5_type=2
if (ctype(1:1) == 'l') h5_type=5

inicio = walltime(wtime_start_shdf5orec)

! Write the dataset.
if (ctype == 'is') then
   call fh5_write(h5_type, ivars, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'rs') then
   call fh5_write(h5_type, rvars, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'cs') then
   call fh5_write(h5_type, cvars, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'ds') then
   call fh5_write(h5_type, dvars, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'ls') then
   call fh5_write(h5_type, lvars, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'ia') then
   call fh5_write(h5_type, ivara, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'ra') then
   call fh5_write(h5_type, rvara, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'ca') then
   call fh5_write(h5_type, cvara, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'da') then
   call fh5_write(h5_type, dvara, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'la') then
   call fh5_write(h5_type, lvara, trim(dsetname)//char(0), hdferr)
endif

write(io6, *) '                             ==T== Tempo gasto na fh5_write(): ',(walltime(wtime_start_shdf5orec)-inicio)

if (hdferr /= 0) then
   print*,'In shdf5_orec: hdf5 write error =',hdferr
   stop 'shdf5_orec: hdf5 write error'
endif

! Close the dataset, the dataspace for the dataset, and the dataspace properties.

inicio = walltime(wtime_start_shdf5orec)

call fh5_close_write(hdferr)

write(io6, *) '                             ==T== Tempo gasto na fh5_close_write(): ',(walltime(wtime_start_shdf5orec)-inicio)

return
end subroutine

!===============================================================================

subroutine shdf5_irec(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                         ,ivars,rvars,cvars,dvars,lvars)

use misc_coms, only: io6
        
implicit none

character(len=*) :: dsetname ! Dataset name
integer :: ndims             ! Number of dimensions or rank
integer, dimension(*) :: dims ! Dataset dimensions.

! Array and scalar arguments for different types. Only specify one in each call.
integer,          optional :: ivara(*),ivars
real,             optional :: rvara(*),rvars
character(len=*), optional :: cvara(*),cvars
real(kind=8),     optional :: dvara(*),dvars
logical,          optional :: lvara(*),lvars

integer:: h5_type   ! Local type designator

integer, dimension(4) :: dimsh ! Dataset dimensions.

integer :: hdferr ! Error flag

character(len=2) :: ctype

real, external :: walltime
real :: wtime_start_shdf5irec, inicio

wtime_start_shdf5irec = walltime(0.)

! Find which data type will be read
    if(present(ivars)) then ; ctype='is'
elseif(present(rvars)) then ; ctype='rs'
elseif(present(cvars)) then ; ctype='cs'
elseif(present(dvars)) then ; ctype='ds'
elseif(present(lvars)) then ; ctype='ls'
elseif(present(ivara)) then ; ctype='ia'
elseif(present(rvara)) then ; ctype='ra'
elseif(present(cvara)) then ; ctype='ca'
elseif(present(dvara)) then ; ctype='da'
elseif(present(lvara)) then ; ctype='la'
else
   print*,'Incorrect or missing data field argument in shdf5_irec'
   stop 'shdf5_irec: bad data field'
endif

! Check dimensions and set compression chunk size

if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
   print*,'Dimension error in shdf5_irec:',ndims,dims(1:ndims)
   stop 'shdf5_irec: bad dims'
endif

dimsh(1:ndims) = dims(1:ndims)
!print*,'-----:',trim(dsetname)
    
inicio = walltime(wtime_start_shdf5irec)

! Prepare file and memory space for the read
call fh5_prepare_read(trim(dsetname)//char(0), ndims, dimsh, hdferr)

write(io6, *) '                                     ==T== Tempo gasto na fh5_prepare_read(): ',(walltime(wtime_start_shdf5irec)-inicio)

if (hdferr < 0) then
   print*,'shdf5_irec: can''t prepare requested field:',trim(dsetname)
   return
endif

! Read data from hyperslab in the file into the hyperslab in memory.
if (ctype(1:1) == 'i') h5_type=1
if (ctype(1:1) == 'r') h5_type=2
if (ctype(1:1) == 'c') h5_type=3
if (ctype(1:1) == 'd') h5_type=4  ! If native precision is 8 bytes, do h5_type=2
if (ctype(1:1) == 'l') h5_type=5

inicio = walltime(wtime_start_shdf5irec)

if (ctype == 'is') then
   call fh5d_read(h5_type,ivars,hdferr)
elseif (ctype == 'rs') then
   call fh5d_read(h5_type,rvars,hdferr)
elseif (ctype == 'cs') then
   call fh5d_read(h5_type,cvars,hdferr)
elseif (ctype == 'ds') then
   call fh5d_read(h5_type,dvars,hdferr)
elseif (ctype == 'ls') then
   call fh5d_read(h5_type,lvars,hdferr)
elseif (ctype == 'ia') then
   call fh5d_read(h5_type,ivara,hdferr)
elseif (ctype == 'ra') then
   call fh5d_read(h5_type,rvara,hdferr)
elseif (ctype == 'ca') then
   call fh5d_read(h5_type,cvara,hdferr)
elseif (ctype == 'da') then
   call fh5d_read(h5_type,dvara,hdferr)
elseif (ctype == 'la') then
   call fh5d_read(h5_type,lvara,hdferr)
endif

write(io6, *) '                                     ==T== Tempo gasto na fh5d_read(): ',(walltime(wtime_start_shdf5irec)-inicio)


if (hdferr /= 0) then
   print*,'shdf5_irec: call fh5d_read: hdf5 error =',hdferr
   stop
endif

! Close the dataset, the dataspace for the dataset, and the memory space.

inicio = walltime(wtime_start_shdf5irec)

call fh5_close_read(hdferr)

write(io6, *) '                                     ==T== Tempo gasto na fh5_close_read(): ',(walltime(wtime_start_shdf5irec)-inicio)


return
end subroutine

!===============================================================================

subroutine shdf5_close()
        
use misc_coms, only: io6, iparallel

implicit none

integer :: hdferr  ! Error flags
real, external :: walltime
real :: wtime_start_shdf5close, inicio

wtime_start_shdf5close = walltime(0.)

! Close RAMS hdf file.

call fh5f_close(hdferr)

inicio = walltime(wtime_start_shdf5close)
write(io6, *) '                                     ==T== Tempo gasto na fh5f_close(): ',(walltime(wtime_start_shdf5close)-inicio)

return
end  subroutine

!===============================================================================

subroutine shdf5_io(action,ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                              ,ivars,rvars,cvars,dvars,lvars)

  implicit none

  character(len=*)           :: dsetname, action
  integer                    :: ndims, dims(*)
  integer,          optional :: ivara(*), ivars
  real,             optional :: rvara(*), rvars
  character(len=*), optional :: cvara(*), cvars
  real(kind=8),     optional :: dvara(*), dvars
  logical,          optional :: lvara(*), lvars
 
  ! THIS ROUTINE CALLS SHDF5_IREC OR SHDF5_OREC TO READ OR WRITE A VARIABLE
  ! DEPENDING ON WHETHER 'ACTION' EQUALS 'READ' OR 'WRITE'

  if (trim(action) == 'READ') then
     
     call shdf5_irec(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
          ,ivars,rvars,cvars,dvars,lvars)
     
  elseif (trim(action) == 'WRITE') then
     
     call shdf5_orec(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
          ,ivars,rvars,cvars,dvars,lvars)
     
  else
     
     print *, "Illegal action in shdf5_io."
     print *, "Action should be 'READ' or 'WRITE'"
     stop     "Ending model run"

  endif
  
end subroutine shdf5_io

end module
