Module damaris_helper

use misc_coms,   only: io6

integer :: dummy

Contains

#ifdef DAMARIS_SUPPORT

subroutine damaris_io(action, iteration, ndims, sdims, dims, dsetname &
                                              ,ivara,rvara,cvara,dvara,lvara  &
                                              ,ivars,rvars,cvars,dvars,lvars)
  implicit none


  character(len=*)           :: dsetname, action
  integer                    :: ndims, dims(*), sdims(*), iteration
  integer,          optional :: ivara(*), ivars
  real,             optional :: rvara(*), rvars
  character(len=*), optional :: cvara(*), cvars
  real(kind=8),     optional :: dvara(*), dvars
  logical,          optional :: lvara(*), lvars

  integer*8                  :: chunk_handle
  integer                    :: err


  ! damaris requires a start and end dimensions; we only need to say it begins
  ! in 0 and ends in the values of dims vector


  call df_chunk_set(ndims,sdims,dims, chunk_handle)

  if ( chunk_handle == 0 ) then
    stop 'ERROR: Damaris Chunk Handle is 0'
  endif


  ! now write it out
  if(present(ivars)) then 
      call df_chunk_write(chunk_handle, dsetname, iteration, ivars, err)
  elseif(present(rvars)) then
      call df_chunk_write(chunk_handle, dsetname, iteration, rvars, err)
  elseif(present(cvars)) then
      call df_chunk_write(chunk_handle, dsetname, iteration, cvars, err)
  elseif(present(dvars)) then
      call df_chunk_write(chunk_handle, dsetname, iteration, dvars, err)
  elseif(present(lvars)) then
      call df_chunk_write(chunk_handle, dsetname, iteration, lvars, err)
  elseif(present(ivara)) then
      call df_chunk_write(chunk_handle, dsetname, iteration, ivara, err)
  elseif(present(rvara)) then
      call df_chunk_write(chunk_handle, dsetname, iteration, rvara, err)
  elseif(present(cvara)) then
      call df_chunk_write(chunk_handle, dsetname, iteration, cvara, err)
  elseif(present(dvara)) then
      call df_chunk_write(chunk_handle, dsetname, iteration, dvara, err)
  elseif(present(lvara)) then
      call df_chunk_write(chunk_handle, dsetname, iteration, lvara, err)
  else
     print*,'Incorrect or missing data field argument in damaris_io'
     stop 'damaris_io: bad data field'
  endif


  if (err < 0 ) then
    write (io6,*) 'error writing damaris variable ', dsetname
    stop 'error writing variable'
  endif


end subroutine damaris_io



end module

#endif
