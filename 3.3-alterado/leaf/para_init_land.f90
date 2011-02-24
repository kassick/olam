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
subroutine para_init_land(landflag)

use misc_coms,  only: io6

use mem_ijtabs, only: itabg_w, mrls

use mem_para,   only: mgroupsize, myrank,  &
!                     send_ul, recv_ul,  &
                      send_wl, recv_wl,  &
                      send_wlf, recv_wlf,  &
!                     nsends_ul, nrecvs_ul,  &
                      nsends_wl, nrecvs_wl,  &
                      nsends_wlf, nrecvs_wlf

use leaf_coms,  only: mml, mul, mwl, nml, nul, nwl, nzg

use mem_leaf,   only: itab_ml, itab_ul, itab_wl, ltab_ml, ltab_ul, ltab_wl,  &
                      itabg_ml, itabg_ul, itabg_wl, alloc_land_grid, land,  &
                      jtab_wl_mpi
                        
use mem_sflux,  only: landflux, nlandflux, jlandflux, landflux_temp

implicit none

logical, intent(in) :: landflag(nwl)

integer :: j,jsend
integer :: iw,iml,iul,iwl
integer :: ilf
integer :: im1l,im2l,iw1l,iw2l

integer :: iml_myrank = 1 ! Counter for ML points to be included on this rank
integer :: iul_myrank = 1 ! Counter for UL points to be included on this rank
integer :: iwl_myrank = 1 ! Counter for WL points to be included on this rank

! Automatic arrays

logical :: myrankflag_ml(nml) ! Flag for ML points existing on this rank
logical :: myrankflag_ul(nul) ! Flag for UL points existing on this rank
logical :: myrankflag_wl(nwl) ! Flag for WL points existing on this rank

real :: xeml_0(nml),yeml_0(nml),zeml_0(nml)
real :: xewl_0(nwl),yewl_0(nwl),zewl_0(nwl)

real :: land_area_0(nwl)

integer :: leaf_class_0(nwl)
integer :: ntext_soil_0(nzg,nwl)


!write(io6,*) 'pal1 '

! Allocate temporary ltab data structures

allocate (ltab_ml(nml)) ! Duplicate for itab_ml
allocate (ltab_ul(nul)) ! Duplicate for itab_ul
allocate (ltab_wl(nwl)) ! Duplicate for itab_wl

! Allocate send & recv counter arrays and initialize to zero

allocate (nsends_wl (mrls)) ; nsends_wl (1:mrls) = 0
allocate (nsends_wlf(mrls)) ; nsends_wlf(1:mrls) = 0
allocate (nrecvs_wl (mrls)) ; nrecvs_wl (1:mrls) = 0
allocate (nrecvs_wlf(mrls)) ; nrecvs_wlf(1:mrls) = 0

! Initialize myrank flag arrays to .false.

myrankflag_ml(1:nml) = .false.
myrankflag_ul(1:nul) = .false.
myrankflag_wl(1:nwl) = .false.

! Copy itab data structures to ltab data structures

ltab_ml(1:nml) = itab_ml(1:nml)
ltab_ul(1:nul) = itab_ul(1:nul)
ltab_wl(1:nwl) = itab_wl(1:nwl)

! Copy grid coordinates and other quantities to temporary arrays

do iml = 2,nml
   xeml_0(iml) = land%xeml(iml)
   yeml_0(iml) = land%yeml(iml)
   zeml_0(iml) = land%zeml(iml)
enddo

do iwl = 2,nwl
   xewl_0(iwl) = land%xewl(iwl)
   yewl_0(iwl) = land%yewl(iwl)
   zewl_0(iwl) = land%zewl(iwl)

   leaf_class_0(iwl) = land%leaf_class(iwl)
   land_area_0(iwl) = land%area(iwl)

   ntext_soil_0(1:nzg,iwl) = land%ntext_soil(1:nzg,iwl)
enddo


!write(io6,*) 'pal2 '

! Deallocate itab data structures and main grid coordinate arrays

deallocate (itab_ml, itab_ul, itab_wl)

deallocate (land%leaf_class)
deallocate (land%area)
deallocate (land%ntext_soil)

deallocate (land%xeml, land%yeml, land%zeml)      
deallocate (land%xewl, land%yewl, land%zewl)      


!write(io6,*) 'pal3 '

!------------------------------------------------------------------------------
! The following section will be required for landcell-to-landcell transport of 
! water via rivers and aquifers, which is not yet implemented in LEAF.  The MPI
! send/recv tables set up by this section will be distinct from those set up
! for landcell-atmosphere fluxes in the next section, and will exchange all
! prognostic LEAF quantities.)
!------------------------------------------------------------------------------

!t Loop over all UL points, and for each whose assigned irank is equal to myrank, 
!t flag its WL neighbors for inclusion on this rank.

!t do iul = 2,nul
!t    if (itabg_ul(iul)%irank == myrank) then

!t       iw1l = ltab_ul(iul)%iw1
!t       iw2l = ltab_ul(iul)%iw2

!t       myrankflag_wl(iw1l) = .true.
!t       myrankflag_wl(iw2l) = .true.

!t    endif
!t enddo

! Loop over all WL points, and for each whose landflag value is .true., flag 
! the WL point for inclusion on this rank.  Then, for any WL point flagged, 
! flag its UL and ML neighbors for inclusion on this rank.

do iwl = 2,nwl
   if (landflag(iwl)) then
      myrankflag_wl(iwl) = .true.
   endif

   if (myrankflag_wl(iwl)) then

      do j = 1,ltab_wl(iwl)%jm
         iml = ltab_wl(iwl)%im(j)
         iul = ltab_wl(iwl)%iu(j)

         myrankflag_ml(iml) = .true.
         myrankflag_ul(iul) = .true.
      enddo

   endif
enddo


!write(io6,*) 'pal4 '

! Loop over all ML, UL, and WL points and count the ones that have been flagged
! for inclusion on this rank.

do iml = 2,nml
   if (myrankflag_ml(iml)) then
      iml_myrank = iml_myrank + 1
   endif
enddo

do iul = 2,nul
   if (myrankflag_ul(iul)) then
      iul_myrank = iul_myrank + 1
   endif
enddo

do iwl = 2,nwl
   if (myrankflag_wl(iwl)) then
      iwl_myrank = iwl_myrank + 1
   endif
enddo


!write(io6,*) 'pal5 '

! Set mml, mul, mwl values for this rank

mml = iml_myrank
mul = iul_myrank
mwl = iwl_myrank

! Allocate itab data structures and main grid coordinate arrays

call alloc_land_grid(mml,mul,mwl,nzg,mgroupsize)


!write(io6,*) 'pal6 '

! Reset point counts to 1

iml_myrank = 1
iul_myrank = 1
iwl_myrank = 1

! Store new myrank ML, UL, WL indices in itabg data structures

do iml = 2,nml
   if (myrankflag_ml(iml)) then
      iml_myrank = iml_myrank + 1

      itabg_ml(iml)%iml_myrank = iml_myrank      
   endif
enddo

do iul = 2,nul
   if (myrankflag_ul(iul)) then
      iul_myrank = iul_myrank + 1

      itabg_ul(iul)%iul_myrank = iul_myrank
   endif
enddo

do iwl = 2,nwl
   if (myrankflag_wl(iwl)) then
      iwl_myrank = iwl_myrank + 1

      itabg_wl(iwl)%iwl_myrank = iwl_myrank      
   endif
enddo


!write(io6,*) 'pal7 '

! Memory copy to main tables

do iml = 2,nml
   if (myrankflag_ml(iml)) then
      iml_myrank = itabg_ml(iml)%iml_myrank

      itab_ml(iml_myrank)%imglobe = iml

      land%xeml(iml_myrank) = xeml_0(iml)
      land%yeml(iml_myrank) = yeml_0(iml)
      land%zeml(iml_myrank) = zeml_0(iml)
   endif
enddo

do iul = 2,nul
   if (myrankflag_ul(iul)) then
      iul_myrank = itabg_ul(iul)%iul_myrank

      itab_ul(iul_myrank)%irank = itabg_ul(iul)%irank
      itab_ul(iul_myrank)%iuglobe = iul

      im1l = ltab_ul(iul)%im1
      im2l = ltab_ul(iul)%im2
      iw1l = ltab_ul(iul)%iw1
      iw2l = ltab_ul(iul)%iw2

      if (myrankflag_ml(im1l)) itab_ul(iul_myrank)%im1 = itabg_ml(im1l)%iml_myrank
      if (myrankflag_ml(im2l)) itab_ul(iul_myrank)%im2 = itabg_ml(im2l)%iml_myrank
      if (myrankflag_wl(iw1l)) itab_ul(iul_myrank)%iw1 = itabg_wl(iw1l)%iwl_myrank
      if (myrankflag_wl(iw2l)) itab_ul(iul_myrank)%iw2 = itabg_wl(iw2l)%iwl_myrank
   endif
enddo


!write(io6,*) 'pal8 '

do iwl = 2,nwl
   if (myrankflag_wl(iwl)) then
      iwl_myrank = itabg_wl(iwl)%iwl_myrank

      itab_wl(iwl_myrank)%irank   = itabg_wl(iwl)%irank
      itab_wl(iwl_myrank)%iwglobe = iwl
      itab_wl(iwl_myrank)%jm      = ltab_wl(iwl)%jm

      land%xewl(iwl_myrank) = xewl_0(iwl)
      land%yewl(iwl_myrank) = yewl_0(iwl)
      land%zewl(iwl_myrank) = zewl_0(iwl)

      land%leaf_class      (iwl_myrank) = leaf_class_0      (iwl)
      land%area            (iwl_myrank) = land_area_0       (iwl) 
      land%ntext_soil(1:nzg,iwl_myrank) = ntext_soil_0(1:nzg,iwl)

      do j = 1,ltab_wl(iwl)%jm
         iml = ltab_wl(iwl)%im(j)
         iul = ltab_wl(iwl)%iu(j)
         
         if (myrankflag_ml(iml)) itab_wl(iwl_myrank)%im(j) = itabg_ml(iml)%iml_myrank
         if (myrankflag_ul(iul)) itab_wl(iwl_myrank)%iu(j) = itabg_ul(iul)%iul_myrank
      enddo

   endif
enddo

! Initialize iremote_u and iremote_w values to 0.


!write(io6,*) 'pal9 '

do j = 1,mgroupsize
!   send_ul(j)%iremote = 0
   send_wl(j)%iremote = 0
   send_wlf(j)%iremote = 0

!   recv_ul(j)%iremote = 0
   recv_wl(j)%iremote = 0
   recv_wlf(j)%iremote = 0
enddo

!------------------------------------------------------------------------------
! The following section will be required for landcell-to-landcell transport of 
! water via rivers and aquifers, which is not yet implemented in LEAF.  The MPI
! send/recv tables set up by this section will be distinct from those set up
! for landcell-atmosphere fluxes in the next section, and will exchange all
! prognostic LEAF quantities.)
!------------------------------------------------------------------------------

! Loop over all UL points and get the 2 neighbor WL indices of each

!t do iul = 2,nul
!t    iw1l = ltab_ul(iul)%iw1
!t    iw2l = ltab_ul(iul)%iw2

!t    if (itabg_ul(iul)%irank == myrank) then

! UL point is on myrank.  If either WL neighbor is on remote rank, myrank must
! receive that WL cell.

!t       if (iw1l > 1 .and. itabg_wl(iw1l)%irank /= myrank) then
!t          call recv_table_wl(itabg_wl(iw1l)%irank)
!t       endif

!t       if (iw2l > 1 .and. itabg_wl(iw2l)%irank /= myrank) then
!t          call recv_table_wl(itabg_wl(iw2l)%irank)
!t       endif

!t    elseif (itabg_ul(iul)%irank /= myrank) then

! UL point is on remote rank.  If either WL neighbor is on myrank, myrank must
! send that WL cell.

!t       if (iw1l > 1 .and. itabg_wl(iw1l)%irank == myrank) then
!t          call send_table_wl(iw1l,itabg_ul(iul)%irank)
!t       endif

!t       if (iw2l > 1 .and. itabg_wl(iw2l)%irank == myrank) then
!t          call send_table_wl(iw2l,itabg_ul(iul)%irank)
!t       endif

!t    endif

!t enddo

! Loop over all landflux cells
! (This section is required for computing landcell-atmosphere fluxes and 
! involves only surface and canopy properties.)


!write(io6,*) 'pal10 '

do ilf = 2,nlandflux


!write(io6,*) 'pal8.1 ',ilf,nlandflux

   iw  = landflux_temp(ilf)%iw   ! full-domain index
   iwl = landflux_temp(ilf)%iwl  ! full-domain index


!write(io6,*) 'pal8.2 ',ilf,iw,iwl

   if (itabg_w(iw)%irank /= myrank .and. itabg_wl(iwl)%irank == myrank) then  


!write(io6,*) 'pal8.3 ',ilf,iw,iwl

! ATM cell is on remote rank while LAND cell is on myrank, so landflux is computed
! on remote rank and must be received by myrank, and LAND cell fields must be
! sent to remote rank.

      call send_table_wl(iwl,itabg_w(iw)%irank)

!write(io6,*) 'pal8.4 ',ilf,iw,iwl

      call recv_table_wlf(itabg_w(iw)%irank)


!write(io6,*) 'pal8.5 ',ilf,iw,iwl

   endif

   if (itabg_w(iw)%irank == myrank .and. itabg_wl(iwl)%irank /= myrank) then  


!write(io6,*) 'pal8.6 ',ilf,iw,iwl

! ATM cell is on myrank while LAND cell is on remote rank, so landflux is computed
! on myrank and must be sent to remote rank, and LAND cell fields must be
! received from remote rank.

      call send_table_wlf(ilf,itabg_wl(iwl)%irank)

!write(io6,*) 'pal8.7 ',ilf,iw,iwl

      call recv_table_wl(itabg_wl(iwl)%irank)


!write(io6,*) 'pal8.8 ',ilf,iw,iwl

   endif


!write(io6,*) 'pal8.9 ',ilf,iw,iwl

enddo


!write(io6,*) 'pal11 '

! Deallocate temporary data structures and arrays

deallocate (ltab_ml,ltab_ul,ltab_wl)


!write(io6,*) 'pal12 '

return
end subroutine para_init_land

!===============================================================================

subroutine recv_table_wl(iremote)

use mem_para,  only: nrecvs_wl, recv_wl
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote is already in table of ranks to receive from

jrecv = 1
do while (jrecv <= nrecvs_wl(1) .and. recv_wl(jrecv)%iremote /= iremote)
   jrecv = jrecv + 1
enddo

! If jrecv exceeds nrecvs_wl(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_wl(1).

if (jrecv > nrecvs_wl(1)) nrecvs_wl(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_wl(jrecv)%iremote = iremote

return
end subroutine recv_table_wl

!===============================================================================

subroutine recv_table_wlf(iremote)

use mem_para,  only: nrecvs_wlf, recv_wlf, myrank
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote is already in table of ranks to receive from

jrecv = 1
do while (jrecv <= nrecvs_wlf(1) .and. recv_wlf(jrecv)%iremote /= iremote)
   jrecv = jrecv + 1
enddo

! If jrecv exceeds nrecvs_wlf(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_wlf(1).

if (jrecv > nrecvs_wlf(1)) nrecvs_wlf(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_wlf(jrecv)%iremote = iremote

return
end subroutine recv_table_wlf

!===============================================================================

subroutine send_table_wl(iwl,iremote)

use mem_leaf,  only: itab_wl, itabg_wl
use mem_para,  only: nsends_wl, send_wl
use max_dims,  only: maxremote
use misc_coms, only: io6

implicit none

integer, intent(in) :: iwl
integer, intent(in) :: iremote

integer :: jsend
integer :: iwl_myrank

! Check whether iremote is already in table of ranks to send to

jsend = 1
do while (jsend <= nsends_wl(1) .and. send_wl(jsend)%iremote /= iremote)
   jsend = jsend + 1
enddo

! If jsend exceeds nsends_wl(1), jsend represents a rank not yet entered in the
! table, so increase nsends_wl(1).

if (jsend > nsends_wl(1)) nsends_wl(1) = jsend

! If nsends_wl(1) exceeds maxremote, print error message and stop

if (jsend > maxremote) then
   write(io6,*) 'In subroutine send_table_wl, nsends_wl(1) exceeds maxremote'
   write(io6,*) 'nsends_wl(1) = ',nsends_wl(1)
   write(io6,*) 'maxremote = ',maxremote
   write(io6,*) 'Stopping model '
   stop 'stopping in send_table_wl '
endif

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

iwl_myrank = itabg_wl(iwl)%iwl_myrank

itab_wl(iwl_myrank)%send_l(jsend) = .true.

send_wl(jsend)%iremote = iremote

return
end subroutine send_table_wl

!===============================================================================

subroutine send_table_wlf(ilf,iremote)

use mem_sflux,  only: landflux, landfluxg
use mem_para,   only: nsends_wlf, send_wlf
use mem_ijtabs, only: nloops_w
use max_dims,   only: maxremote
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ilf
integer, intent(in) :: iremote

integer :: jsend
integer :: ilf_myrank

! Check whether iremote is already in table of ranks to send to

jsend = 1
do while (jsend <= nsends_wlf(1) .and. send_wlf(jsend)%iremote /= iremote)
   jsend = jsend + 1
enddo

! If jsend exceeds nsends_wlf(1), jsend represents a rank not yet entered in the
! table, so increase nsends_wlf(1).

if (jsend > nsends_wlf(1)) nsends_wlf(1) = jsend

! If nsends_wlf(1) exceeds maxremote, print error message and stop

if (jsend > maxremote) then
   write(io6,*) 'In subroutine send_table_wlf, nsends_wlf(1) exceeds maxremote'
   write(io6,*) 'nsends_wlf(1) = ',nsends_wlf(1)
   write(io6,*) 'maxremote = ',maxremote
   write(io6,*) 'Stopping model '
   stop 'stopping in send_table_wlf '
endif

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

ilf_myrank = landfluxg(ilf)%ilf_myrank

landflux(ilf_myrank)%send_lf(jsend) = .true.

send_wlf(jsend)%iremote = iremote

return
end subroutine send_table_wlf

