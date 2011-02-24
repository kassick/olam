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
subroutine para_init_sea(seaflag)

use misc_coms,  only: io6

use mem_ijtabs, only: itabg_w, mrls

use mem_para,   only: mgroupsize, myrank,  &
!                     send_us, recv_us,  &
                      send_ws, recv_ws,  &
                      send_wsf, recv_wsf,  &
!                     nsends_us, nrecvs_us,  &
                      nsends_ws, nrecvs_ws,  &
                      nsends_wsf, nrecvs_wsf

use sea_coms,   only: mms, mus, mws, nms, nus, nws

use mem_sea,    only: itab_ms, itab_us, itab_ws, ltab_ms, ltab_us, ltab_ws,  &
                      itabg_ms, itabg_us, itabg_ws, alloc_sea_grid, sea,  &
                      jtab_ws_mpi
                        
use mem_sflux,  only: seaflux, nseaflux, jseaflux, seaflux_temp

implicit none

logical, intent(in) :: seaflag(nws)

integer :: j,jsend
integer :: iw,ims,ius,iws
integer :: isf
integer :: im1s,im2s, iw1s,iw2s

integer :: ims_myrank = 1 ! Counter for MS points to be included on this rank
integer :: ius_myrank = 1 ! Counter for US points to be included on this rank
integer :: iws_myrank = 1 ! Counter for WS points to be included on this rank

! Automatic arrays

logical :: myrankflag_ms(nms) ! Flag for MS points existing on this rank
logical :: myrankflag_us(nus) ! Flag for US points existing on this rank
logical :: myrankflag_ws(nws) ! Flag for WS points existing on this rank

real :: xems_0(nms),yems_0(nms),zems_0(nms)
real :: xews_0(nws),yews_0(nws),zews_0(nws)

real :: sea_area_0(nws)


!write(io6,*) 'pas1 ',allocated(seaflux(303)%send_sf)

! Allocate temporary ltab data structures

allocate (ltab_ms(nms)) ! Duplicate for itab_ms
allocate (ltab_us(nus)) ! Duplicate for itab_us
allocate (ltab_ws(nws)) ! Duplicate for itab_ws

! Allocate send & recv counter arrays and initialize to zero

allocate (nsends_ws (mrls)) ; nsends_ws (1:mrls) = 0
allocate (nsends_wsf(mrls)) ; nsends_wsf(1:mrls) = 0
allocate (nrecvs_ws (mrls)) ; nrecvs_ws (1:mrls) = 0
allocate (nrecvs_wsf(mrls)) ; nrecvs_wsf(1:mrls) = 0

! Initialize myrank flag arrays to .false.

myrankflag_ms(1:nms) = .false.
myrankflag_us(1:nus) = .false.
myrankflag_ws(1:nws) = .false.

! Copy itab data structures to ltab data structures

ltab_ms(1:nms) = itab_ms(1:nms)
ltab_us(1:nus) = itab_us(1:nus)
ltab_ws(1:nws) = itab_ws(1:nws)

! Copy grid coordinates and other quantities to temporary arrays

do ims = 2,nms
   xems_0(ims) = sea%xems(ims)
   yems_0(ims) = sea%yems(ims)
   zems_0(ims) = sea%zems(ims)
enddo

do iws = 2,nws
   xews_0(iws) = sea%xews(iws)
   yews_0(iws) = sea%yews(iws)
   zews_0(iws) = sea%zews(iws)

   sea_area_0(iws) = sea%area(iws)
enddo

! Deallocate itab data structures and main grid coordinate arrays

deallocate (itab_ms, itab_us, itab_ws)

deallocate (sea%area)

deallocate (sea%xems, sea%yems, sea%zems)      
deallocate (sea%xews, sea%yews, sea%zews)      


!write(io6,*) 'pas2 ',allocated(seaflux(303)%send_sf)

!------------------------------------------------------------------------------
! The following section (which would correspond to landcell-to-landcell exchange
! of information between parallel processors) is likely not required for sea 
! cells.  The prognostic ocean model will perform its own parallel exchange of
! information.
!------------------------------------------------------------------------------

!t Loop over all US points, and for each whose assigned irank is equal to myrank, 
!t flag its WS neighbors for inclusion on this rank.

!t do ius = 2,nus
!t    if (itabg_us(ius)%irank == myrank) then

!t       iw1s = ltab_us(ius)%iw1
!t       iw2s = ltab_us(ius)%iw2

!t       myrankflag_ws(iw1s) = .true.
!t       myrankflag_ws(iw2s) = .true.

!t    endif
!t enddo

! Loop over all WS points, and for each whose seaflag value is .true., flag
! the WS point for inclusion on this rank.  Then, for any WS point flagged, 
! flag its US and MS neighbors for inclusion on this rank.

do iws = 2,nws
   if (seaflag(iws)) then
      myrankflag_ws(iws) = .true.
   endif

   if (myrankflag_ws(iws)) then

      do j = 1,ltab_ws(iws)%jm
         ims = ltab_ws(iws)%im(j)
         ius = ltab_ws(iws)%iu(j)

         myrankflag_ms(ims) = .true.
         myrankflag_us(ius) = .true.
      enddo

   endif
enddo


!write(io6,*) 'pas3 '

! Loop over all MS, US, and WS points and count the ones that have been flagged
! for inclusion on this rank.

do ims = 2,nms
   if (myrankflag_ms(ims)) then
      ims_myrank = ims_myrank + 1
   endif
enddo

do ius = 2,nus
   if (myrankflag_us(ius)) then
      ius_myrank = ius_myrank + 1
   endif
enddo

do iws = 2,nws
   if (myrankflag_ws(iws)) then
      iws_myrank = iws_myrank + 1
   endif
enddo


!write(io6,*) 'pas4 '

! Set mms, mus, mws values for this rank

mms = ims_myrank
mus = ius_myrank
mws = iws_myrank

! Allocate itab data structures and main grid coordinate arrays

call alloc_sea_grid(mms,mus,mws,mgroupsize)

! Reset point counts to 1

ims_myrank = 1
ius_myrank = 1
iws_myrank = 1

! Store new myrank MS, US, WS indices in itabg data structures


!write(io6,*) 'pas5 ',allocated(seaflux(303)%send_sf)

do ims = 2,nms
   if (myrankflag_ms(ims)) then
      ims_myrank = ims_myrank + 1

      itabg_ms(ims)%ims_myrank = ims_myrank      
   endif
enddo

do ius = 2,nus
   if (myrankflag_us(ius)) then
      ius_myrank = ius_myrank + 1

      itabg_us(ius)%ius_myrank = ius_myrank
   endif
enddo

do iws = 2,nws
   if (myrankflag_ws(iws)) then
      iws_myrank = iws_myrank + 1

      itabg_ws(iws)%iws_myrank = iws_myrank
   endif
enddo

! Memory copy to main tables


!write(io6,*) 'pas6 '

do ims = 2,nms
   if (myrankflag_ms(ims)) then
      ims_myrank = itabg_ms(ims)%ims_myrank

      itab_ms(ims_myrank)%imglobe = ims

      sea%xems(ims_myrank) = xems_0(ims)
      sea%yems(ims_myrank) = yems_0(ims)
      sea%zems(ims_myrank) = zems_0(ims)
   endif
enddo

do ius = 2,nus
   if (myrankflag_us(ius)) then
      ius_myrank = itabg_us(ius)%ius_myrank

      itab_us(ius_myrank)%irank = itabg_us(ius)%irank
      itab_us(ius_myrank)%iuglobe = ius

      im1s = ltab_us(ius)%im1
      im2s = ltab_us(ius)%im2
      iw1s = ltab_us(ius)%iw1
      iw2s = ltab_us(ius)%iw2

!--------------------------------------
if (iw1s < 1) print*, 'iw1s ',ius,iw1s
if (iw2s < 1) print*, 'iw2s ',ius,iw2s
!--------------------------------------

      if (myrankflag_ms(im1s)) itab_us(ius_myrank)%im1 = itabg_ms(im1s)%ims_myrank
      if (myrankflag_ms(im2s)) itab_us(ius_myrank)%im2 = itabg_ms(im2s)%ims_myrank
      if (myrankflag_ws(iw1s)) itab_us(ius_myrank)%iw1 = itabg_ws(iw1s)%iws_myrank
      if (myrankflag_ws(iw2s)) itab_us(ius_myrank)%iw2 = itabg_ws(iw2s)%iws_myrank
   endif
enddo

do iws = 2,nws
   if (myrankflag_ws(iws)) then
      iws_myrank = itabg_ws(iws)%iws_myrank

      itab_ws(iws_myrank)%irank   = itabg_ws(iws)%irank
      itab_ws(iws_myrank)%iwglobe = iws
      itab_ws(iws_myrank)%jm      = ltab_ws(iws)%jm

      sea%xews(iws_myrank) = xews_0(iws)
      sea%yews(iws_myrank) = yews_0(iws)
      sea%zews(iws_myrank) = zews_0(iws)

      sea%area(iws_myrank) = sea_area_0(iws)

      do j = 1,ltab_ws(iws)%jm
         ims = ltab_ws(iws)%im(j)
         ius = ltab_ws(iws)%iu(j)

         if (myrankflag_ms(ims)) itab_ws(iws_myrank)%im(j) = itabg_ms(ims)%ims_myrank
         if (myrankflag_us(ius)) itab_ws(iws_myrank)%iu(j) = itabg_us(ius)%ius_myrank
      enddo

   endif
enddo


!write(io6,*) 'pas7 ',allocated(seaflux(303)%send_sf)

! Initialize iremote_u and iremote_w values to 0.

do j = 1,mgroupsize
!   send_us(j)%iremote = 0
   send_ws(j)%iremote = 0
   send_wsf(j)%iremote = 0

!   recv_us(j)%iremote = 0
   recv_ws(j)%iremote = 0
   recv_wsf(j)%iremote = 0
enddo

!------------------------------------------------------------------------------
! The following section (which would correspond to landcell-to-landcell exchange
! of information between parallel processors) is likely not required for sea 
! cells.  The prognostic ocean model will perform its own parallel exchange of
! information.
!------------------------------------------------------------------------------

!t Loop over all US points and get the 2 neighbor WS indices of each

!t do ius = 2,nus
!t    iw1s = ltab_us(ius)%iw1
!t    iw2s = ltab_us(ius)%iw2

!t    if (itabg_us(ius)%irank == myrank) then

! US point is on myrank.  If either WS neighbor is on remote rank, myrank must
! receive that WS cell.

!t       if (iw1s > 1 .and. itabg_ws(iw1s)%irank /= myrank) then
!t          call recv_table_ws(itabg_ws(iw1s)%irank)
!t       endif

!t       if (iw2s > 1 .and. itabg_ws(iw2s)%irank /= myrank) then
!t          call recv_table_ws(itabg_ws(iw2s)%irank)
!t       endif

!t    elseif (itabg_us(ius)%irank /= myrank) then

! US point is on remote rank.  If either WS neighbor is on myrank, myrank must
! send that WS cell.

!t       if (iw1s > 1 .and. itabg_ws(iw1s)%irank == myrank) then
!t          call send_table_ws(iw1s,itabg_us(ius)%irank)
!t       endif

!t       if (iw2s > 1 .and. itabg_ws(iw2s)%irank == myrank) then
!t          call send_table_ws(iw2s,itabg_us(ius)%irank)
!t       endif

!t    endif

!t enddo

! Loop over all seaflux cells
! (This section is required for computing seacell-atmosphere fluxes and 
! involves only surface and canopy properties.)


!write(io6,*) 'pas8 '

do isf = 2,nseaflux


!write(io6,*) 'pas8.1 ',isf,nseaflux,allocated(seaflux(303)%send_sf)

   iw  = seaflux_temp(isf)%iw   ! full-domain index
   iws = seaflux_temp(isf)%iws  ! full-domain index


!write(io6,*) 'pas8.2 ',isf,iw,iws,allocated(seaflux(303)%send_sf)

   if (itabg_w(iw)%irank /= myrank .and. itabg_ws(iws)%irank == myrank) then  


!write(io6,*) 'pas8.3 ',isf,iw,iws

! ATM cell is on remote rank while SEA cell is on myrank, so seaflux is computed
! on remote rank and must be received by myrank, and SEA cell fields must be
! sent to remote rank.

      call send_table_ws(iws,itabg_w(iw)%irank)

!write(io6,*) 'pas8.4 ',isf,iw,iws,allocated(seaflux(303)%send_sf)

      call recv_table_wsf(itabg_w(iw)%irank)


!write(io6,*) 'pas8.5 ',isf,iw,iws

   endif

   if (itabg_w(iw)%irank == myrank .and. itabg_ws(iws)%irank /= myrank) then  


!write(io6,*) 'pas8.6 ',isf,iw,iws,itabg_ws(iws)%irank,allocated(seaflux(303)%send_sf)

! ATM cell is on myrank while SEA cell is on remote rank, so seaflux is computed
! on myrank and must be sent to remote rank, and SEA cell fields must be
! received from remote rank.

      call send_table_wsf(isf,itabg_ws(iws)%irank)

!write(io6,*) 'pas8.7 ',isf,iw,iws

      call recv_table_ws(itabg_ws(iws)%irank)


!write(io6,*) 'pas8.8 ',isf,iw,iws

   endif


!write(io6,*) 'pas8.9 ',isf,iw,iws

enddo


!write(io6,*) 'pas9 '

! Deallocate temporary data structures and arrays

deallocate (ltab_ms,ltab_us,ltab_ws)

return
end subroutine para_init_sea

!===============================================================================

subroutine recv_table_ws(iremote)

use mem_para,  only: nrecvs_ws, recv_ws
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote is already in table of ranks to receive from

jrecv = 1
do while (jrecv <= nrecvs_ws(1) .and. recv_ws(jrecv)%iremote /= iremote)
   jrecv = jrecv + 1
enddo

! If jrecv exceeds nrecvs_ws(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_ws(1).

if (jrecv > nrecvs_ws(1)) nrecvs_ws(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_ws(jrecv)%iremote = iremote

return
end subroutine recv_table_ws

!===============================================================================

subroutine recv_table_wsf(iremote)

use mem_para,  only: nrecvs_wsf, recv_wsf
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote is already in table of ranks to receive from

jrecv = 1
do while (jrecv <= nrecvs_wsf(1) .and. recv_wsf(jrecv)%iremote /= iremote)
   jrecv = jrecv + 1
enddo

! If jrecv exceeds nrecvs_wsf(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_wsf(1).

if (jrecv > nrecvs_wsf(1)) nrecvs_wsf(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_wsf(jrecv)%iremote = iremote

return
end subroutine recv_table_wsf

!===============================================================================

subroutine send_table_ws(iws,iremote)

use mem_sea,   only: itab_ws, itabg_ws
use mem_para,  only: nsends_ws, send_ws
use max_dims,  only: maxremote
use misc_coms, only: io6

implicit none

integer, intent(in) :: iws
integer, intent(in) :: iremote

integer :: jsend
integer :: iws_myrank

! Check whether iremote is already in table of ranks to send to

jsend = 1
do while (jsend <= nsends_ws(1) .and. send_ws(jsend)%iremote /= iremote)
   jsend = jsend + 1
enddo

! If jsend exceeds nsends_ws(1), jsend represents a rank not yet entered in the
! table, so increase nsends_ws(1).

if (jsend > nsends_ws(1)) nsends_ws(1) = jsend

! If nsends_ws(1) exceeds maxremote, print error message and stop

if (jsend > maxremote) then
   write(io6,*) 'In subroutine send_table_ws, nsends_ws(1) exceeds maxremote'
   write(io6,*) 'nsends_ws(1) = ',nsends_ws(1)
   write(io6,*) 'maxremote = ',maxremote
   write(io6,*) 'Stopping model '
   stop 'stopping in send_table_ws '
endif

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

iws_myrank = itabg_ws(iws)%iws_myrank

itab_ws(iws_myrank)%send_s(jsend) = .true.

send_ws(jsend)%iremote = iremote

return
end subroutine send_table_ws

!===============================================================================

subroutine send_table_wsf(isf,iremote)

use mem_sflux, only: seaflux, seafluxg
use mem_para,  only: nsends_wsf, send_wsf
use max_dims,  only: maxremote
use misc_coms, only: io6

implicit none

integer, intent(in) :: isf
integer, intent(in) :: iremote

integer :: jsend
integer :: isf_myrank

! Check whether iremote is already in table of ranks to send to

!write(io6,*) 'pbs1 ',isf,iremote

jsend = 1
do while (jsend <= nsends_wsf(1) .and. send_wsf(jsend)%iremote /= iremote)
   jsend = jsend + 1
enddo


!write(io6,*) 'pbs2 ',isf,iremote,jsend

! If jsend exceeds nsends_wsf(1), jsend represents a rank not yet entered in the
! table, so increase nsends_wsf(1).

if (jsend > nsends_wsf(1)) nsends_wsf(1) = jsend


!write(io6,*) 'pbs3 ',isf,iremote,jsend

! If nsends_wsf(1) exceeds maxremote, print error message and stop

if (jsend > maxremote) then
   write(io6,*) 'In subroutine send_table_wsf, nsends_wsf(1) exceeds maxremote'
   write(io6,*) 'nsends_wsf(1) = ',nsends_wsf(1)
   write(io6,*) 'maxremote = ',maxremote
   write(io6,*) 'Stopping model '
   stop 'stopping in send_table_wsf '
endif

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

isf_myrank = seafluxg(isf)%isf_myrank


!write(io6,*) 'pbs4 ',isf,iremote,jsend,isf_myrank


!write(io6,*) 'pbs4.1 ',allocated(seaflux)
!write(io6,*) 'pbs4.2 ',size(seaflux)
!write(io6,*) 'pbs4.3 ',allocated(seaflux(isf_myrank)%send_sf)
!write(io6,*) 'pbs4.4 ',size(seaflux(isf_myrank)%send_sf)

seaflux(isf_myrank)%send_sf(jsend) = .true.


!write(io6,*) 'pbs5 ',isf,iremote,jsend

send_wsf(jsend)%iremote = iremote


!write(io6,*) 'pbs6 ',isf,iremote,jsend

return
end subroutine send_table_wsf

