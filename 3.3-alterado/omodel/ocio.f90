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

subroutine commio(action)

! THIS ROUTINE READS OR WRITES NON-HORIZONTALLY VARYING FIELDS
! THAT ARE COMMON TO ALL PROCESSES. CALL THIS ROUTINE WITH
! 'ACTION' EQUALS 'READ' OR 'WRITE' TO INPUT OR OUTPUT THE FIELDS

use misc_coms,  only: io6, itime1, idate1, imonth1, iyear1, nxp, ngrids,  &
                      centlat, centlon, grdlen, grdwid, grdaxis, nzp,  &
                      mdomain, deltaz, deltax, dzrat, dzmax,  &
                      itopoflg, time8, zbase
use mem_grid,   only: nza, mua, mwa, mma, mza, &
                      zm, zt, dzm ,dzt, dzim, dzit
use leaf_coms,  only: nzg, nzs, slz, ivegflg, isfcl
use hdf5_utils, only: shdf5_orec, shdf5_irec, shdf5_io

implicit none
integer          :: ndims, idims(2), k
character(len=*) :: action

if ((trim(action) /= "READ") .and. (trim(action) /= "WRITE")) then
   write(io6,*) 'Illegal action in routine commio'
   write(io6,*) 'action must be "READ" or "WRITE"'
   stop     'Stopping run'
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GROUP1: NAMELIST PARAMETERS THAT WE DON'T WANT CHANGED ON HISTORY RESTART.
!         THESE SHOULD CORRESPOND WITH THE VARIABLES IN THE NOT_HISTORY 
!         SECTION OF THE SUBROUTINE COPY_NL.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ndims = 1
idims(1) = 1
call shdf5_io(action, ndims, idims, 'nl%itime1',  ivars=itime1)
call shdf5_io(action, ndims, idims, 'nl%idate1' , ivars=idate1)
call shdf5_io(action, ndims, idims, 'nl%imonth1', ivars=imonth1)
call shdf5_io(action, ndims, idims, 'nl%iyear1',  ivars=iyear1)
call shdf5_io(action, ndims, idims, 'nl%ngrids',  ivars=ngrids)
call shdf5_io(action, ndims, idims, 'nl%nzp',     ivars=nzp)
call shdf5_io(action, ndims, idims, 'nl%nzg',     ivars=nzg)
call shdf5_io(action, ndims, idims, 'nl%nzs',     ivars=nzs)
call shdf5_io(action, ndims, idims, 'nl%nxp',     ivars=nxp)
call shdf5_io(action, ndims, idims, 'nl%mdomain', ivars=mdomain)
call shdf5_io(action, ndims, idims, 'nl%isfcl',   ivars=isfcl)
call shdf5_io(action, ndims, idims, 'nl%itopoflg',ivars=itopoflg)
call shdf5_io(action, ndims, idims, 'nl%ivegflg', ivars=ivegflg)
call shdf5_io(action, ndims, idims, 'nl%deltax',  rvars=deltax)

call shdf5_io(action, ndims, idims, 'nl%deltaz',  rvars=deltaz)
call shdf5_io(action, ndims, idims, 'nl%dzrat',   rvars=dzrat)
call shdf5_io(action, ndims, idims, 'nl%dzmax',   rvars=dzmax)
call shdf5_io(action, ndims, idims, 'nl%zbase',   rvars=zbase)

ndims = 1
idims(1) = ngrids
call shdf5_io(action, ndims, idims, 'nl%centlat', rvara=centlat)
call shdf5_io(action, ndims, idims, 'nl%centlon', rvara=centlon)
call shdf5_io(action, ndims, idims, 'nl%grdlen',  rvara=grdlen)
call shdf5_io(action, ndims, idims, 'nl%grdwid',  rvara=grdwid)
call shdf5_io(action, ndims, idims, 'nl%grdaxis', rvara=grdaxis)

ndims=1
idims(1) = nzg
call shdf5_io(action, ndims, idims, 'nl%slz', rvara=slz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GROUP 2: NON-HORIZONTALLY-VARYING (COMMOM) MODEL VARIABLES
!           THAT ARE NOT PART OF THE VTABLES I/O MECHANISM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ndims=1
idims(1) = 1
call shdf5_io(action, ndims, idims, 'time8', dvars=time8)

return
end subroutine commio
