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
subroutine zero_massflux(umarusc,wmarwsc,rho_old)

use mem_ijtabs, only: jtab_u, jtab_w, istp, mrl_begl
use mem_grid,   only: mza, mua, mwa, lpw
use mem_basic,  only: rho
use misc_coms,  only: io6
use rastro_evts

!$ use omp_lib

implicit none

real, intent(out) :: umarusc(mza,mua)
real, intent(out) :: wmarwsc(mza,mwa)
real, intent(out) :: rho_old(mza,mwa)

integer :: j,k,iu,iw,mrl

#ifdef OLAM_RASTRO
character*1 :: rst_buf = '_'
call rst_event_s_f(OLAM_ZERO_MASSFLUX_IN,rst_buf)
#endif


! Zero out long timestep mass flux components (used for scalar advective 
! transport) so they may be summed over small timesteps

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_IN,rst_buf)
#endif
!$omp parallel do private (iu,k)
do j = 1,jtab_u(14)%jend(mrl); iu = jtab_u(14)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)
   do k = 1,mza-1          ! begin at level 2 even if below ground
      umarusc(k,iu) = 0.   ! initialize horiz mass flux for scalar advection
   enddo
enddo
!$omp end parallel do
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_OUT,rst_buf)
#endif
endif
call rsub('U',14)

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_IN,rst_buf)
#endif
!$omp parallel do private (iw,k)
do j = 1,jtab_w(18)%jend(mrl); iw = jtab_w(18)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = 1,lpw(iw)-2
      rho_old(k,iw) = 0.0
   enddo
   do k = lpw(iw)-1,mza-1
      wmarwsc(k,iw) = 0.        ! initialize vert mass flux for scalar advection
      rho_old(k,iw) = rho(k,iw) ! Save DTL density for use with scalar updates
   enddo
enddo
!$omp end parallel do
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_OUT,rst_buf)
#endif
endif
call rsub('W',18)

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_ZERO_MASSFLUX_OUT,rst_buf)
#endif

return
end subroutine zero_massflux

!===========================================================================

subroutine timeavg_massflux(umarusc,wmarwsc)

use mem_ijtabs, only: jtab_u, itab_u, jtab_w, itab_w, istp, mrl_endl
use mem_grid,   only: mza, mua, mwa, lpu, lpw
use misc_coms,  only: io6, nacoust
use rastro_evts

!$ use omp_lib

implicit none

real, intent(inout) :: umarusc(mza,mua)
real, intent(inout) :: wmarwsc(mza,mwa)

integer :: j,k,iu,iw,mrl,mrlu,mrlw
real :: acoi,acoi2

#ifdef OLAM_RASTRO
character*1 :: rst_buf = '_'
call rst_event_s_f(OLAM_TIMEAVG_MASSFLUX_IN,rst_buf)
#endif


call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_IN,rst_buf)
#endif
!$omp parallel do private (iu,mrlu,k,acoi)
do j = 1,jtab_u(20)%jend(mrl); iu = jtab_u(20)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)
   mrlu = itab_u(iu)%mrlu
   acoi = 1. / float(nacoust(mrlu))
   do k = lpu(iu),mza-1
      umarusc(k,iu) = umarusc(k,iu) * acoi  ! upsum
   enddo
enddo
!$omp end parallel do
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_OUT,rst_buf)
#endif
endif
call rsub('U',20)

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_IN,rst_buf)
#endif
!$omp parallel do private (iw,mrlw,k,acoi2)
do j = 1,jtab_w(25)%jend(mrl); iw = jtab_w(25)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   mrlw = itab_w(iw)%mrlw
   acoi2 = 1. / float(nacoust(mrlw))
   do k = lpw(iw),mza-2
      wmarwsc(k,iw) = wmarwsc(k,iw) * acoi2  ! wpsum
   enddo
enddo
!$omp end parallel do
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_OUT,rst_buf)
#endif
endif
call rsub('W',25)

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_TIMEAVG_MASSFLUX_OUT,rst_buf)
#endif

return
end subroutine timeavg_massflux

!===========================================================================

subroutine tridiffo(m1,ka,kz,cim1,ci,cip1,rhs,soln)

implicit none

integer, intent(in) :: m1,ka,kz
real, intent(in) :: cim1(m1),ci(m1),cip1(m1),rhs(m1)
real, intent(out) :: soln(m1)
real :: scr1(m1)  ! automatic array
real :: scr2(m1)  ! automatic array

integer :: k
real cji

scr1(ka) = cip1(ka) / ci(ka)
scr2(ka) = rhs(ka) / ci(ka)

do k = ka+1,kz
   soln(k) = ci(k) - cim1(k) * scr1(k-1)
   cji = 1. / soln(k)
   scr1(k) = cip1(k) * cji
   scr2(k) = (rhs(k) - cim1(k) * scr2(k-1)) * cji
enddo

soln(kz) = scr2(kz)

do k = kz-1,ka,-1
   soln(k) = scr2(k) - scr1(k) * soln(k+1)
enddo

return
end subroutine tridiffo

