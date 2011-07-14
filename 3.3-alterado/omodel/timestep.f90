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
subroutine timestep()

use rastro_evts
use misc_coms,   only: io6, time8, time_istp8, nqparm, initial, ilwrtyp,   &
                       iswrtyp, dtsm, nqparm_sh, dtlm, iparallel,   &
                       s1900_init, s1900_sim
use mem_ijtabs,  only: nstp, istp, mrls, leafstep
use mem_nudge,   only: nudflag
use mem_grid,    only: mza, mua, mwa
use micro_coms,  only: level
use leaf_coms,   only: isfcl
use mem_para,    only: myrank

use mem_basic  ! needed only when print statements below are uncommented
use mem_tend   ! needed only when print statements below are uncommented
use mem_leaf   ! needed only when print statements below are uncommented
use ed_options, only: frq_phenology

implicit none

integer :: is,mvl,isl4,jstp
real :: t1,w1

! automatic arrays

real :: umarusc    (mza,mua) ! U mom density times U-face surface area [kg/s]
real :: wmarwsc    (mza,mwa) ! W mom density times W-face surface area [kg/s]
real :: rho_old    (mza,mwa) ! density from beginning of timestep [kg/m^3]
real :: alpha_press(mza,mwa) ! 
real :: rhot       (mza,mwa) ! grid-cell total mass tendency [kg/s]

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAM_TIMESTEP_IN,rst_buf)
#endif

! +----------------------------------------------------------------------------+
! |  Each call to subroutine timestep drives all steps in advancing by dtlong  |
! +----------------------------------------------------------------------------+

time_istp8 = time8

if (time8 < 1.e-3) then
!   call bubble()
endif

do jstp = 1,nstp  ! nstp = no. of finest-grid-level aco steps in dtlm(1)
#ifdef OLAM_RASTRO
   call rst_event_ii_f(OLAM_INNERSTEP_IN,istp,jstp)
#endif
   istp = jstp

   call tend0(rhot)

! call check_nans(1)

   if (ilwrtyp + iswrtyp > 0) then
      call radiate()
   endif

! call check_nans(2)

   call surface_turb_flux()

! call check_nans(3)

   call turb_k()  ! Compute K's

! call check_nans(4)

   if (iparallel == 1) then
      call mpi_send_w('K',wmarwsc)  ! Send K's; wmarwsc used as dummy array here
   endif

! call check_nans(5)

   if (any(nqparm   (1:mrls) > 0) .or.  &
       any(nqparm_sh(1:mrls) > 0)) then

      call cuparm_driver(rhot)

      if (isfcl == 1) then
         call surface_cuparm_flux()
      endif
   endif

! call check_nans(6)

   if (iparallel == 1) then
      call mpi_recv_w('K',wmarwsc)  ! Recv K's; wmarwsc used as dummy array here
   endif

! call check_nans(7)

   call thiltend_long(alpha_press,rhot)

! call check_nans(8)

   call veltend_long()

! call check_nans(9)

   if (initial == 2 .and. nudflag == 1)  &
   call obs_nudge(rhot)

! call check_nans(10)

   call zero_massflux(umarusc,wmarwsc,rho_old)

! call check_nans(11)

   call prog_wrtu(umarusc,wmarwsc,alpha_press,rhot)

! call check_nans(12)

   call timeavg_massflux(umarusc,wmarwsc)

! call check_nans(14)

   call scalar_transport(umarusc,wmarwsc,rho_old)

! call check_nans(15)

   call predtr(rho_old)

! call check_nans(16)

   if (iparallel == 1) then
      call mpi_recv_u('U')  ! Recv U group
   endif

! call check_nans(17)

   if (level /= 3) then
      call thermo()
   endif

! call check_nans(18)

   if (level == 3) then
      call micro()  ! maybe later make freq. uniform

      if (isfcl == 1) then
         call surface_precip_flux()
      endif
   endif

! call check_nans(19)

   call trsets()  

! call check_nans(20)

   if (iparallel == 1) then
      call mpi_send_w('T',wmarwsc)  ! Send W group; wmarwsc used as dummy array here
      call mpi_recv_w('T',wmarwsc)  ! Recv W group
   endif

! call check_nans(21)

   if (iparallel == 1) then
      call mpi_send_w('S',wmarwsc)  ! Send scalars; wmarwsc used as dummy array here
      call mpi_recv_w('S',wmarwsc)  ! Recv scalars; wmarwsc used as dummy array here
   endif

! call check_nans(22)

   if (leafstep(istp) > 0) then
      call leaf3()

      if (iparallel == 1) call mpi_send_wl('T')
      if (iparallel == 1) call mpi_recv_wl('T')

      call seacells()

      if (iparallel == 1) call mpi_send_ws('T')
      if (iparallel == 1) call mpi_recv_ws('T')
   endif

   time_istp8 = time8 + float(istp) * dtsm(mrls)  ! Update precise time
   s1900_sim = s1900_init + time_istp8

! call check_nans(23)

#ifdef OLAM_RASTRO
   call rst_event_ii_f(OLAM_INNERSTEP_OUT,istp,jstp)
#endif

enddo

! Call ED model if it is time to do vegetation dynamics

if (mod(real(time8)+dtlm(1),frq_phenology) < dtlm(1)) call ed_vegetation_dynamics()  

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_TIMESTEP_OUT,rst_buf)
#endif

return
end subroutine timestep

!==========================================================================

!!subroutine check_nans(icall)
!!
!!use mem_basic,  only: sh_w,rho,thil,sh_v,wc,wmc,press,umc,uc
!!use mem_micro,  only: sh_c,sh_r
!!use mem_grid,   only: mza,mwa,lpw,volti,mua
!!use mem_tend,   only: thilt, umt
!!use micro_coms, only: level
!!use misc_coms,  only: io6, iparallel
!!use mem_ijtabs, only: itab_u, itab_w
!!
!!implicit none
!!  
!!integer :: i,k,iu
!!integer, intent(in) :: icall
!!
!!!do iu = 1,mua
!!!   if (itab_u(iu)%iuglobe == 1206) then
!!!      write(io6,'(a,i6,5e15.7)') 'cns ',icall,umc(2,iu),uc(2,iu),umt(2,iu)
!!!   endif
!!!enddo
!!
!!return
!!
!!do i = 2,mwa
!!   do k = lpw(i),mza
!!      if (isnan(sh_w (k,i)) .or.  &
!!          isnan(rho  (k,i)) .or.  &
!!          isnan(thil (k,i)) .or.  &
!!          isnan(thilt(k,i)) .or.  &
!!          isnan(press(k,i)) .or.  &
!!          isnan(wc   (k,i)) .or.  &
!!          isnan(wmc  (k,i)) .or.  &
!!          isnan(sh_v (k,i)) .or.  &
!!          thil(k,i) < 100.0) then
!!
!!         write(io6,*) ''
!!         write(io6,*) 'check_nans',k,i,icall
!!         write(io6,*) 'sh_w,rho,thil',sh_w(k,i),rho(k,i),thil(k,i)
!!         write(io6,*) 'thilt,press',thilt(k,i),press(k,i)
!!         write(io6,*) 'wc, wmc, sh_v',wc(k,i),wmc(k,i),sh_v(k,i)
!!         stop
!!      endif
!!
!!!      if (level >= 3) then
!!!         if (isnan(sh_c(k,i)) .or.  &
!!!             isnan(sh_r(k,i))) then
!!               
!!!            write(io6,*) 'nan',k,i,icall
!!!            write(io6,*) 'sh_c,sh_r',sh_c(k,i),sh_r(k,i)
!!!            stop
!!!         endif
!!!      endif
!!
!!   enddo
!!enddo
!!
!!return
!!end subroutine check_nans

!==========================================================================

subroutine tend0(rhot)

use mem_ijtabs, only: jtab_w, jtab_u, istp, mrl_begl
use var_tables, only: scalar_tab, num_scalar
use mem_grid,   only: mza, mwa, mua, lcu, lpw
use mem_tend,   only: wmt, umt, thilt
use misc_coms,  only: io6
use rastro_evts

!$ use omp_lib

implicit none

real, intent(in) :: rhot

integer :: n,mrl,j,k,iw,iu

#ifdef OLAM_RASTRO
call rst_event_r_f(OLAM_TEND0_IN,rhot)
#endif

! SET SCALAR TENDENCIES TO ZERO

mrl = mrl_begl(istp)

if (mrl > 0) then
   do n = 1,num_scalar
      call tnd0(scalar_tab(n)%var_t)
   enddo
   call tnd0(thilt)
   call tnd0(rhot)
endif

! SET W MOMENTUM TENDENCY TO ZERO

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_IN,rst_buf)
#endif
!$omp parallel do private(iw,k)
do j = 1,jtab_w(14)%jend(mrl); iw = jtab_w(14)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = lpw(iw),mza-1
      wmt(k,iw) = 0.
   enddo
enddo
!$omp end parallel do
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_OUT,rst_buf)
#endif
endif
call rsub('W',14)


! SET U MOMENTUM TENDENCY TO ZERO

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_IN,rst_buf)
#endif
!$omp parallel do private(iu,k)
do j = 1,jtab_u(11)%jend(mrl); iu = jtab_u(11)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)
   do k = lcu(iu),mza-1
      umt(k,iu) = 0.
   enddo
enddo
!$omp end parallel do
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_OUT,rst_buf)
#endif
endif
call rsub('U',11)

#ifdef OLAM_RASTRO
call rst_event_r_f(OLAM_TEND0_OUT,rhot)
#endif

return
end subroutine tend0

!==========================================================================

subroutine tnd0(vart)

use mem_ijtabs, only: jtab_w, istp, mrl_begl
use mem_grid,   only: mza, mwa, lpw
use misc_coms,  only: io6
use rastro_evts
!$ use omp_lib

implicit none

real, intent(out) :: vart(mza,mwa)

integer :: j,iw,k,mrl

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_IN,rst_buf)
#endif
!$omp parallel do private(iw,k)
do j = 1,jtab_w(11)%jend(mrl); iw = jtab_w(11)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = lpw(iw)-1,mza
      vart(k,iw) = 0.
   enddo
enddo
!$omp end parallel do
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_OUT,rst_buf)
#endif
endif
call rsub('W',11)

return
end subroutine tnd0

!==========================================================================

subroutine predtr(rho_old)

use var_tables, only: num_scalar, scalar_tab
use mem_ijtabs, only: istp, jtab_w, mrl_endl
use mem_grid,   only: mza, mwa
use misc_coms,  only: io6
use rastro_evts

implicit none

real, intent(in) :: rho_old(mza,mwa)

integer :: n,mrl

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAM_PREDTR_IN,rst_buf)
#endif


!   -  Step thermodynamic variables from  t  to  t+1.
!   -  Set top, lateral and bottom boundary conditions on some variables
!        if needed.
!   -  Call adjustment to assure all positive definite quantities
!        remain positive.
!   -  Rediagnose some thermodynamic quantities for use on the small
!        timestep.

!     Update the scalars and apply lateral, top, and bottom boundary
!     conditions.

mrl = mrl_endl(istp)
if (mrl > 0) then
do n = 1,num_scalar
   call o_update(n,scalar_tab(n)%var_p,scalar_tab(n)%var_t,rho_old)
enddo
endif

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PREDTR_OUT,rst_buf)
#endif

return
end subroutine predtr

!==========================================================================

subroutine o_update(n,varp,vart,rho_old)

use mem_ijtabs, only: jtab_w, istp, itab_w, mrl_endl
use mem_basic,  only: rho
use misc_coms,  only: io6, dtlm
use mem_grid,   only: mza, mwa, lpw
use rastro_evts

!$ use omp_lib

implicit none

integer, intent(in) :: n
real, intent(in) :: vart(mza,mwa),rho_old(mza,mwa)
real, intent(out) :: varp(mza,mwa)

integer :: iw,j,k,mrl
real :: dtl

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_IN,rst_buf)
#endif
!$omp parallel do private(iw,k,dtl)
do j = 1,jtab_w(27)%jend(mrl); iw = jtab_w(27)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   dtl = dtlm(itab_w(iw)%mrlw)
   do k = lpw(iw),mza-1
      varp(k,iw) = (varp(k,iw) * rho_old(k,iw) + dtl * vart(k,iw)) / rho(k,iw) 
   enddo
enddo
!$omp end parallel do
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_OUT,rst_buf)
#endif
endif
call rsub('W',27)!RRR

return
end subroutine o_update

!==========================================================================

subroutine bubble()

use mem_basic, only: thil, theta
use misc_coms, only: io6
!use mem_grid
!use mem_ijtabs

implicit none

integer :: iw,i,j,k


   do k = 2,2
      thil(k,17328) = thil(k,17328) + 5. 
      theta(k,17328) = theta(k,17328) + 5. 

      thil(k,17329) = thil(k,17329) + 5. 
      theta(k,17329) = theta(k,17329) + 5. 

      thil(k,17333) = thil(k,17333) + 5. 
      theta(k,17333) = theta(k,17333) + 5. 

      thil(k,17334) = thil(k,17334) + 5. 
      theta(k,17334) = theta(k,17334) + 5. 

      thil(k,17335) = thil(k,17335) + 5. 
      theta(k,17335) = theta(k,17335) + 5. 

      thil(k,17336) = thil(k,17336) + 5. 
      theta(k,17336) = theta(k,17336) + 5. 
   enddo


return
end subroutine bubble
