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
subroutine olam_mpi_init()

use mem_para,  only: mgroupsize, myrank,                    &
                     send_u, send_w, send_uf,               &
                     recv_u, recv_w, recv_uf,               &
                     send_wl, send_wlf, send_ws, send_wsf,  &
                     recv_wl, recv_wlf, recv_ws, recv_wsf

use misc_coms, only: io6, iparallel

implicit none

#ifdef OLAM_MPI

include 'mpif.h'
integer :: ierr

! Initialize MPI and determine process groupsize and myrank

call MPI_Init(ierr)
call MPI_Comm_size(MPI_COMM_WORLD,mgroupsize,ierr)
call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)

! Allocate send and recv tables

allocate(send_u (mgroupsize))
allocate(send_w (mgroupsize))
allocate(send_uf(mgroupsize))

allocate(recv_u (mgroupsize))
allocate(recv_w (mgroupsize))
allocate(recv_uf(mgroupsize))

allocate(send_wl (mgroupsize))
allocate(send_wlf(mgroupsize))
allocate(send_ws (mgroupsize))
allocate(send_wsf(mgroupsize))

allocate(recv_wl (mgroupsize))
allocate(recv_wlf(mgroupsize))
allocate(recv_ws (mgroupsize))
allocate(recv_wsf(mgroupsize))

#else

mgroupsize = 1
myrank = 0

#endif

return
end subroutine olam_mpi_init

!===============================================================================

subroutine olam_mpi_finalize()

use mem_para, only: send_u, send_w, recv_u,                &
                    recv_w, send_uf, recv_uf,              &
                    send_wl, send_wlf, send_ws, send_wsf,  &
                    recv_wl, recv_wlf, recv_ws, recv_wsf

implicit none

#ifdef OLAM_MPI

include 'mpif.h'
integer :: ierr

call MPI_Finalize(ierr)

if (allocated(send_u))  deallocate(send_u)
if (allocated(send_w))  deallocate(send_w)
if (allocated(send_uf)) deallocate(send_uf)

if (allocated(recv_u))  deallocate(recv_u)
if (allocated(recv_w))  deallocate(recv_w)
if (allocated(recv_uf)) deallocate(recv_uf)

if (allocated(send_wl))   deallocate(send_wl)
if (allocated(send_wlf))  deallocate(send_wlf)
if (allocated(send_ws))   deallocate(send_ws)
if (allocated(send_wsf))  deallocate(send_wsf)

if (allocated(recv_wl))  deallocate(recv_wl)
if (allocated(recv_wlf)) deallocate(recv_wlf)
if (allocated(recv_ws))  deallocate(recv_ws)
if (allocated(recv_wsf)) deallocate(recv_wsf)

#endif

return
end subroutine olam_mpi_finalize

!==================================================================

subroutine olam_alloc_mpi(mza,mrls)

use mem_ijtabs, only: jtab_u, jtab_w
use mem_para,   only: myrank, nrecvs_u, nrecvs_w, nsends_u, nsends_w,  &
                      recv_u, recv_w, recv_uf, send_u, send_w, send_uf
use misc_coms,  only: io6
use var_tables, only: nvar_par

implicit none

integer, intent(in) :: mza
integer, intent(in) :: mrls

#ifdef OLAM_MPI

include 'mpif.h'

integer :: nbytes_int
integer :: nbytes_real
integer :: nbytes_real8

integer :: nbytes_per_iu
integer :: nbytes_per_iuf
integer :: nbytes_per_iw

integer :: itag10 = 10
integer :: itag110 = 110

integer :: status(MPI_STATUS_SIZE)
integer :: ierr
integer :: jsend
integer :: jrecv

integer :: nupts, nwpts
integer :: mrl
integer :: ibuf(2)

! allocate send buffers

call MPI_Pack_size(1,MPI_INTEGER,MPI_COMM_WORLD,nbytes_int  ,ierr)
call MPI_Pack_size(1,MPI_REAL   ,MPI_COMM_WORLD,nbytes_real ,ierr)
call MPI_Pack_size(1,MPI_REAL8  ,MPI_COMM_WORLD,nbytes_real8,ierr)

! Determine number of bytes to send per IU column

nbytes_per_iu = nbytes_int  &
              + mza * 2 * nbytes_real

nbytes_per_iuf = nbytes_int  &
              + mza * 3 * nbytes_real

! Loop over all U sends for mrl = 1

do jsend = 1,nsends_u(1)

! Determine size of send_u buffer for mrl = 1

   send_u(jsend)%nbytes = nbytes_int   &
                        + nbytes_per_iu * jtab_u(25+jsend)%jend(1)

   send_uf(jsend)%nbytes = nbytes_int   &
                        + nbytes_per_iuf * jtab_u(25+jsend)%jend(1)
                      
! Allocate buffer

   allocate(send_u(jsend)%buff(send_u(jsend)%nbytes))
   allocate(send_uf(jsend)%buff(send_uf(jsend)%nbytes))

! Initialize send_u%irequest to 'null' value

   send_u(jsend)%irequest  = -999
   send_uf(jsend)%irequest = -999

! Send buffer sizes to receive ranks
   
   ibuf(1) = send_u(jsend)%nbytes
   ibuf(2) = send_uf(jsend)%nbytes

   call MPI_Send(ibuf,2,MPI_INTEGER  &
                ,send_u(jsend)%iremote,itag10,MPI_COMM_WORLD,ierr)

!----------------------------------------------------------
   write(io6,*) myrank, 'sent u/uf ',send_u(jsend)%iremote, itag10
!----------------------------------------------------------

! Loop over all mrl values greater than 1

   do mrl = 2,mrls

! Send jtab_u(25+jsend)%jend(mrl) to receive ranks

      call MPI_Send(jtab_u(25+jsend)%jend(mrl),1,MPI_INTEGER  &
                   ,send_u(jsend)%iremote,itag10+mrl,MPI_COMM_WORLD,ierr)


!----------------------------------------------------------
   write(io6,*) myrank, 'sent u2 ',send_u(jsend)%iremote, itag10+mrl
!----------------------------------------------------------

! If at least 1 U point needs to be sent to current remote rank for 
! current mrl, increase nsends_u(mrl) by 1.

      if (jtab_u(25+jsend)%jend(mrl) > 0) nsends_u(mrl) = nsends_u(mrl) + 1

   enddo

enddo

! Determine number of bytes to send per IW column

nbytes_per_iw = nbytes_int                                     &
              + mza * max(3 * nbytes_real8 + 2 * nbytes_real,  &
                                      nvar_par * nbytes_real)

! Loop over all W sends for mrl = 1

do jsend = 1,nsends_w(1)

! Determine size of send_w buffer for mrl = 1

   send_w(jsend)%nbytes = nbytes_int                               &
                        + nbytes_per_iw * jtab_w(35+jsend)%jend(1)

! Allocate buffer

   allocate(send_w(jsend)%buff(send_w(jsend)%nbytes))

! Initialize send_w%irequest to 'null' value

   send_w(jsend)%irequest = -999
   
! Send buffer sizes to receive ranks

   call MPI_Send(send_w(jsend)%nbytes,1,MPI_INTEGER  &
                ,send_w(jsend)%iremote,itag110,MPI_COMM_WORLD,ierr)


!----------------------------------------------------------
   write(io6,*) myrank, 'sent w ',send_w(jsend)%iremote, itag110
!----------------------------------------------------------

! Loop over all mrl values greater than 1

   do mrl = 2,mrls

! Send jtab_w(35+jsend)%jend(mrl) to receive ranks

      call MPI_Send(jtab_w(35+jsend)%jend(mrl),1,MPI_INTEGER  &
                   ,send_w(jsend)%iremote,itag110+mrl,MPI_COMM_WORLD,ierr)


!----------------------------------------------------------
   write(io6,*) myrank, 'sent w2 ',send_w(jsend)%iremote, itag110+mrl
!----------------------------------------------------------

! If at least 1 W point needs to be sent to current remote rank for 
! current mrl, increase nsends_w(mrl) by 1.

      if (jtab_w(35+jsend)%jend(mrl) > 0) nsends_w(mrl) = nsends_w(mrl) + 1

   enddo

enddo

! Loop over all U receives for mrl = 1

do jrecv = 1,nrecvs_u(1)

! Get recv_u buffer sizes


!----------------------------------------------------------
   write(io6,*) myrank, 'receiving u/uf ',recv_u(jrecv)%iremote,itag10
!----------------------------------------------------------

   call MPI_Recv(ibuf,2,MPI_INTEGER  &
                ,recv_u(jrecv)%iremote,itag10,MPI_COMM_WORLD,status,ierr)


!----------------------------------------------------------
   write(io6,*) myrank, 'received u/uf ',recv_u(jrecv)%iremote,itag10
!----------------------------------------------------------

   recv_u(jrecv)%nbytes  = ibuf(1)
   recv_uf(jrecv)%nbytes = ibuf(2)

! Allocate recv_u buffers

   allocate(recv_u(jrecv)%buff(recv_u(jrecv)%nbytes))
   allocate(recv_uf(jrecv)%buff(recv_uf(jrecv)%nbytes))

! Initialize recv_u%irequest to 'null' value

   recv_u(jrecv)%irequest = -999
   recv_uf(jrecv)%irequest = -999

! Loop over all mrl values greater than 1

   do mrl = 2,mrls

! Get number of receive points for this mrl and this receive rank


!----------------------------------------------------------
   write(io6,*) myrank, 'receiving u2 ',recv_u(jrecv)%iremote, itag10+mrl
!----------------------------------------------------------

      call MPI_Recv(nupts,1,MPI_INTEGER  &
                   ,recv_u(jrecv)%iremote,itag10+mrl,MPI_COMM_WORLD,status,ierr)


!----------------------------------------------------------
   write(io6,*) myrank, 'received u2 ',recv_u(jrecv)%iremote, itag10+mrl
!----------------------------------------------------------

! If at least 1 U point needs to be received from current remote rank for 
! current mrl, increase nrecvs_u(mrl) by 1.

      if (nupts > 0) nrecvs_u(mrl) = nrecvs_u(mrl) + 1

   enddo

enddo

! Loop over all W receives for mrl = 1

do jrecv = 1,nrecvs_w(1)

! Get recv_w buffer sizes


!----------------------------------------------------------
   write(io6,*) myrank, 'receiving w ',recv_w(jrecv)%iremote, itag110
!----------------------------------------------------------

   call MPI_Recv(recv_w(jrecv)%nbytes,1,MPI_INTEGER  &
                ,recv_w(jrecv)%iremote,itag110,MPI_COMM_WORLD,status,ierr)


!----------------------------------------------------------
   write(io6,*) myrank, 'received w ',recv_w(jrecv)%iremote, itag110
!----------------------------------------------------------

! Allocate recv_w buffers

   allocate(recv_w(jrecv)%buff(recv_w(jrecv)%nbytes))

! Initialize recv_w%irequest to 'null' value

   recv_w(jrecv)%irequest = -999

! Loop over all mrl values greater than 1

   do mrl = 2,mrls

! Get number of receive points for this mrl and this receive rank


!----------------------------------------------------------
   write(io6,*) myrank, 'receiving w2 ',recv_w(jrecv)%iremote, itag110+mrl
!----------------------------------------------------------

      call MPI_Recv(nwpts,1,MPI_INTEGER  &
                   ,recv_w(jrecv)%iremote,itag110+mrl,MPI_COMM_WORLD,status,ierr)


!----------------------------------------------------------
   write(io6,*) myrank, 'received w2 ',recv_w(jrecv)%iremote, itag110+mrl
!----------------------------------------------------------

! If at least 1 W point needs to be received from current remote rank for 
! current mrl, increase nrecvs_w(mrl) by 1.

      if (nwpts > 0) nrecvs_w(mrl) = nrecvs_w(mrl) + 1

   enddo

enddo

#endif

return
end subroutine olam_alloc_mpi

!===============================================================================

subroutine mpi_send_u(sendgroup)

! Subroutine to perform a parallel MPI send of a "U group"
! of field variables

use mem_basic,  only: umc,uc
use mem_para,   only: send_u, recv_u, nsends_u, nrecvs_u
use mem_ijtabs, only: itab_u, jtab_u, mrl_begs, istp
use mem_grid,   only: mza
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: sendgroup

#ifdef OLAM_MPI

include 'mpif.h'

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: itag1 = 1
integer :: mrl
integer :: j
integer :: iu
integer :: iuglobe

! Set MRL and return if mrl < 1 for this step

if (sendgroup == 'I') then
   mrl = 1
else
   mrl = mrl_begs(istp)
endif

if (mrl < 1) return

! Before we send anything, post the receives

do jrecv = 1,nrecvs_u(mrl)

   call MPI_Irecv(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,MPI_PACKED,  &
                  recv_u(jrecv)%iremote,itag1,MPI_COMM_WORLD,          &
                  recv_u(jrecv)%irequest,ierr                          )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_u(mrl)

   ipos = 0

   call MPI_Pack(jtab_u(25+jsend)%jend(mrl),1,MPI_INTEGER,  &
      send_u(jsend)%buff,send_u(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,jtab_u(25+jsend)%jend(mrl)
      iu = jtab_u(25+jsend)%iu(j)
      iuglobe = itab_u(iu)%iuglobe
!----------------------------------------------------------------
      call qsub('U',iu)

      call MPI_Pack(iuglobe,1,MPI_INTEGER,  &
         send_u(jsend)%buff,send_u(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      call MPI_Pack(umc(1,iu),mza,MPI_REAL,  &
         send_u(jsend)%buff,send_u(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      call MPI_Pack(uc(1,iu),mza,MPI_REAL,  &
         send_u(jsend)%buff,send_u(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   enddo
   call rsub('Usend',25+jsend)

   call MPI_Isend(send_u(jsend)%buff,ipos,MPI_PACKED,          &
                  send_u(jsend)%iremote,itag1,MPI_COMM_WORLD,  &
                  send_u(jsend)%irequest,ierr                  )

enddo

#endif

return
end subroutine mpi_send_u

!===============================================================================

subroutine mpi_send_uf(hcnum_u,hcnum_w,hflux_t)

! Subroutine to perform a parallel MPI send of a "U group"
! of field variables

use mem_para,   only: send_uf, recv_uf, nsends_u, nrecvs_u
use mem_ijtabs, only: jtab_u, itab_u, mrl_begs, istp
use mem_grid,   only: mza, mua
use misc_coms,  only: io6

implicit none

real, intent(in) :: hcnum_u(mza,mua)
real, intent(in) :: hcnum_w(mza,mua)
real, intent(in) :: hflux_t(mza,mua)

#ifdef OLAM_MPI

include 'mpif.h'

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: itag2 = 2
integer :: mrl
integer :: j
integer :: iu
integer :: iuglobe

! Return if mrl < 1 for this step

mrl = mrl_begs(istp)
if (mrl < 1) return

! Before we send anything, post the receives

do jrecv = 1,nrecvs_u(mrl)

   call MPI_Irecv(recv_uf(jrecv)%buff,recv_uf(jrecv)%nbytes,MPI_PACKED,  &
                  recv_uf(jrecv)%iremote,itag2,MPI_COMM_WORLD,          &
                  recv_uf(jrecv)%irequest,ierr                          )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_u(mrl)

   ipos = 0

   call MPI_Pack(jtab_u(25+jsend)%jend(mrl),1,MPI_INTEGER,  &
      send_uf(jsend)%buff,send_uf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,jtab_u(25+jsend)%jend(mrl)
      iu = jtab_u(25+jsend)%iu(j)
      iuglobe = itab_u(iu)%iuglobe
!----------------------------------------------------------------
      call qsub('U',iu)

      call MPI_Pack(iuglobe,1,MPI_INTEGER,  &
         send_uf(jsend)%buff,send_uf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      call MPI_Pack(hcnum_u(1,iu),mza,MPI_REAL,  &
         send_uf(jsend)%buff,send_uf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      call MPI_Pack(hcnum_w(1,iu),mza,MPI_REAL,  &
         send_uf(jsend)%buff,send_uf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      call MPI_Pack(hflux_t(1,iu),mza,MPI_REAL,  &
         send_uf(jsend)%buff,send_uf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   enddo

   call rsub('Usendf',25+jsend)

   call MPI_Isend(send_uf(jsend)%buff,ipos,MPI_PACKED,          &
                  send_uf(jsend)%iremote,itag2,MPI_COMM_WORLD,  &
                  send_uf(jsend)%irequest,ierr                  )

enddo

#endif

return
end subroutine mpi_send_uf

!=============================================================================

subroutine mpi_send_w(sendgroup,wmarw)

! Subroutine to perform a parallel MPI send of a "W group" or "S group"
! of field variables

use mem_basic,  only: wmc,wc,thil,rho,press
use mem_turb,   only: hkm,vkm,vkm_sfc
use var_tables, only: nvar_par, vtab_par
use mem_ijtabs, only: jtab_w, itab_w, mrl_begs, mrl_begl, istp
use mem_grid,   only: mza, mwa
use misc_coms,  only: io6
use micro_coms, only: level
use mem_para,   only: nrecvs_w, nsends_w, recv_w, send_w

implicit none

character(1), intent(in) :: sendgroup
real(kind=8), intent(in) :: wmarw(mza,mwa)

#ifdef OLAM_MPI

include 'mpif.h'

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: itag3 = 3
integer :: mrl
integer :: j
integer :: iw
integer :: iwglobe

! Set MRL and return if mrl < 1 for this step

if (sendgroup == 'I') then
   mrl = 1
elseif (sendgroup == 'S') then
   mrl = mrl_begl(istp)
else
   mrl = mrl_begs(istp)
endif

if (mrl < 1) return

! Before we send anything, post the receives

do jrecv = 1,nrecvs_w(mrl)

   call MPI_Irecv(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,MPI_PACKED,  &
                  recv_w(jrecv)%iremote,itag3,MPI_COMM_WORLD,          &
                  recv_w(jrecv)%irequest,ierr                          )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_w(mrl)

   ipos = 0

   call MPI_Pack(jtab_w(35+jsend)%jend(mrl),1,MPI_INTEGER,  &
      send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,jtab_w(35+jsend)%jend(mrl)
      iw = jtab_w(35+jsend)%iw(j)
      iwglobe = itab_w(iw)%iwglobe
!----------------------------------------------------------------
      call qsub('W',iw)

      call MPI_Pack(iwglobe,1,MPI_INTEGER,  &
         send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      if (sendgroup == 'S') then

         do ivar = 1,nvar_par

            call MPI_Pack(vtab_par(ivar)%rvar2_p(1,iw),mza,MPI_REAL,  &
               send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         enddo

      elseif (sendgroup == 'I') then

         call MPI_Pack(press(1,iw),mza,MPI_REAL8,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(rho(1,iw),mza,MPI_REAL8,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(wc(1,iw),mza,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(thil(1,iw),mza,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'K') then

         call MPI_Pack(hkm(1,iw),mza,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(vkm(1,iw),mza,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(vkm_sfc(iw),1,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'P') then

         call MPI_Pack(wmarw(1,iw),mza,MPI_REAL8,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(press(1,iw),mza,MPI_REAL8,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(rho(1,iw),mza,MPI_REAL8,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'T') then

         call MPI_Pack(wc(1,iw),mza,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(thil(1,iw),mza,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         if (level == 3) then
            call MPI_Pack(rho(1,iw),mza,MPI_REAL8,  &
               send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
         endif
 
      endif

   enddo
   call rsub('Wsend',35+jsend)

   call MPI_Isend(send_w(jsend)%buff,ipos,MPI_PACKED,          &
                  send_w(jsend)%iremote,itag3,MPI_COMM_WORLD,  &
                  send_w(jsend)%irequest,ierr                  )

enddo

#endif

return
end subroutine mpi_send_w

!=============================================================================

subroutine mpi_recv_u(recvgroup)

! Subroutine to perform a parallel MPI receive of a "U group"
! of field variables

use mem_basic,  only: umc,uc
use mem_para,   only: send_u, recv_u, nsends_u, nrecvs_u, mgroupsize
use mem_ijtabs, only: itabg_u, mrl_begs, istp
use mem_grid,   only: mza
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: recvgroup

#ifdef OLAM_MPI

include 'mpif.h'

integer :: ierr,ipos
integer :: jrecv,jsend,ivar,jtmp
integer :: nupts
integer :: mrl
integer :: j
integer :: iu
integer :: iuglobe
integer :: status(MPI_STATUS_SIZE)
integer :: statuses(MPI_STATUS_SIZE,mgroupsize)

! Set MRL and return if mrl < 1 for this step

if (recvgroup == 'I') then
   mrl = 1
else
   mrl = mrl_begs(istp)
endif

if (mrl < 1) return

!  Now, let's wait on our receives

do jtmp = 1,nrecvs_u(mrl)

   call MPI_Waitany(nrecvs_u(mrl), recv_u(1:nrecvs_u(mrl))%irequest, jrecv, status, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,ipos,  &
      nupts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,nupts
      call MPI_Unpack(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,ipos,  &
         iuglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      iu = itabg_u(iuglobe)%iu_myrank
!----------------------------------------------------------------
      call qsub('U',iu)

      call MPI_Unpack(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,ipos,  &
         umc(1,iu),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

      call MPI_Unpack(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,ipos,  &
         uc(1,iu),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

   enddo
   call rsub('Urecv',25+jrecv)

enddo

! Make sure all of our sends are finished and de-allocated
call MPI_Waitall(nsends_u(mrl), send_u(1:nsends_u(mrl))%irequest, statuses, ierr)

#endif

return
end subroutine mpi_recv_u

!=============================================================================

subroutine mpi_recv_uf(hcnum_u,hcnum_w,hflux_t)

! Subroutine to perform a parallel MPI receive of a "U group"
! of field variables

use mem_para,   only: send_uf, recv_uf, nsends_u, nrecvs_u, mgroupsize
use mem_ijtabs, only: mrl_begs, itabg_u, istp
use mem_grid,   only: mza, mua
use misc_coms,  only: io6

implicit none

real, intent(inout) :: hcnum_u(mza,mua)
real, intent(inout) :: hcnum_w(mza,mua)
real, intent(inout) :: hflux_t(mza,mua)

#ifdef OLAM_MPI

include 'mpif.h'

integer :: ierr,ipos
integer :: jrecv,jsend,ivar,jtmp
integer :: nupts
integer :: mrl
integer :: j
integer :: iu
integer :: iuglobe
integer :: status(MPI_STATUS_SIZE)
integer :: statuses(MPI_STATUS_SIZE,mgroupsize)

! Return if mrl < 1 for this step

mrl = mrl_begs(istp)
if (mrl < 1) return

!  Now, let's wait on our receives

do jtmp = 1,nrecvs_u(mrl)

   call MPI_Waitany(nrecvs_u(mrl), recv_uf(1:nrecvs_u(mrl))%irequest, jrecv, status, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_uf(jrecv)%buff,recv_uf(jrecv)%nbytes,ipos,  &
      nupts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,nupts
      call MPI_Unpack(recv_uf(jrecv)%buff,recv_uf(jrecv)%nbytes,ipos,  &
         iuglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      iu = itabg_u(iuglobe)%iu_myrank
!----------------------------------------------------------------
      call qsub('U',iu)

      call MPI_Unpack(recv_uf(jrecv)%buff,recv_uf(jrecv)%nbytes,ipos,  &
         hcnum_u(1,iu),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

      call MPI_Unpack(recv_uf(jrecv)%buff,recv_uf(jrecv)%nbytes,ipos,  &
         hcnum_w(1,iu),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

      call MPI_Unpack(recv_uf(jrecv)%buff,recv_uf(jrecv)%nbytes,ipos,  &
         hflux_t(1,iu),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

   enddo
   call rsub('Urecv',25+jrecv)

enddo

! Make sure all of our sends are finished and de-allocated
call MPI_Waitall(nsends_u(mrl), send_uf(1:nsends_u(mrl))%irequest, statuses, ierr)

#endif

return
end subroutine mpi_recv_uf

!=============================================================================

subroutine mpi_recv_w(recvgroup,wmarw)

! Subroutine to perform a parallel MPI receive of a "W group" or "S group"
! of field variables

use mem_basic,  only: wmc,wc,thil,rho,press
use mem_turb,   only: hkm,vkm,vkm_sfc
use var_tables, only: vtab_par, nvar_par
use mem_para,   only: nrecvs_w, nsends_w, recv_w, send_w, mgroupsize
use mem_ijtabs, only: itabg_w, mrl_begs, mrl_begl, istp
use mem_grid,   only: mza, mwa
use misc_coms,  only: io6
use micro_coms, only: level

implicit none

character(1), intent(in) :: recvgroup
real(kind=8), intent(inout) :: wmarw(mza,mwa)

#ifdef OLAM_MPI

include 'mpif.h'

integer :: ierr,ipos
integer :: jrecv,jsend,ivar,jtmp
integer :: nwpts
integer :: mrl
integer :: j
integer :: iw
integer :: iwglobe
integer :: status(MPI_STATUS_SIZE)
integer :: statuses(MPI_STATUS_SIZE,mgroupsize)

! Set MRL and return if mrl < 1 for this step

if (recvgroup == 'I') then
   mrl = 1
elseif (recvgroup == 'S') then
   mrl = mrl_begl(istp)
else
   mrl = mrl_begs(istp)
endif

if (mrl < 1) return

!  Now, let's wait on our receives

do jtmp = 1,nrecvs_w(mrl)

   call MPI_Waitany(nrecvs_w(mrl), recv_w(1:nrecvs_w(mrl))%irequest, jrecv, status, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
      nwpts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,nwpts
      call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
         iwglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      iw = itabg_w(iwglobe)%iw_myrank
!----------------------------------------------------------------
      call qsub('W',iw)

      if (recvgroup == 'S') then

         do ivar = 1,nvar_par

            call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
               vtab_par(ivar)%rvar2_p(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         enddo

      elseif (recvgroup == 'I') then

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            press(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            rho(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            wc(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            thil(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

      elseif (recvgroup == 'K') then

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            hkm(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            vkm(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            vkm_sfc(iw),1,MPI_REAL,MPI_COMM_WORLD,ierr)

      elseif (recvgroup == 'P') then
            
         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            wmarw(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            press(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            rho(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

      elseif (recvgroup == 'T') then

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            wc(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            thil(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         if (level == 3) then
            call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
               rho(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)
         endif

      endif

   enddo
   call rsub('Wrecv',35+jrecv)

enddo

! Make sure all of our sends are finished and de-allocated
call MPI_Waitall(nsends_w(mrl), send_w(1:nsends_w(mrl))%irequest, statuses, ierr)

#endif

return
end subroutine mpi_recv_w

!=============================================================================

subroutine olam_stop(message)
use mem_para,  only: myrank
implicit none

#ifdef OLAM_MPI
include 'mpif.h'
integer :: ierr
#endif

character(*), intent(in) :: message

#ifdef MPI  
  write(*,'(A,I0,A)') "Node ", myrank, ":"
  write(*,'(A)') "STOP "//message
  call mpi_abort(MPI_COMM_WORLD,1,ier)
  stop
#else
  write(*,*) "STOPPING: "//message
  stop
#endif

end subroutine olam_stop



