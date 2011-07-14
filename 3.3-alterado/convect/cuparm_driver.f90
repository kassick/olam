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
subroutine cuparm_driver(rhot)

use mem_grid,    only: mwa, mza, lpw, volt, mua
use grell_coms,  only: alloc_grell
use misc_coms,   only: io6, time_istp8, nqparm, nqparm_sh, confrq, dtlong,  &
                       initial, itime1
use mem_ijtabs,  only: itab_w, jtab_w, mrl_begl, istp
use mem_cuparm,  only: thsrc, rtsrc, thsrcsh, rtsrcsh, aconpr, conprr
use mem_tend,    only: thilt, sh_wt

use mem_basic,   only: uc, wc, theta, press, rho, sh_v
use consts_coms, only: tkmin, r8
use mem_micro,   only: sh_c
use rastro_evts

!$ use omp_lib

implicit none

real, intent(inout) :: rhot(mza,mwa)
real, parameter :: cptime = 7200.0

integer, save :: ialloc = 0
integer :: j
integer :: iw
integer :: k
integer :: mrl, mrlw

#ifdef OLAM_RASTRO
character*1 :: rst_buf = '_'
call rst_event_s_f(OLAM_CUPARM_DRIVER_IN,rst_buf)
#endif


! Memory allocation for Grell cumulus transport

if (ialloc == 0)  then
   call alloc_grell(mza, mwa)
   ialloc = 1 
endif

! If model run has been initialized from observational data, avoid cumulus
! parameterization during an initial period to allow gravity waves to settle.

if ((initial == 2) .and. (time_istp8 < real(cptime,r8))) return

! Check whether it is time to update cumulus parameterization tendencies

if ((istp == 1) .and. (mod(time_istp8+0.001_r8,real(confrq,r8)) < dtlong)) then

! Print message that cumulus parameterization is being computed

   write(io6, '(A,f0.2,A,f0.2)') &
        ' Convective tendencies updated at time = ', time_istp8, &
        ' Real time (hrs) = ', time_istp8/3600.+(itime1/100+mod(itime1,100)/60.)

! Initialize indices for search for maximum heating rate
! (commented out because incorrect in parallel operation)

!    dthmax = 0.
!    iwqmax = 0
!    kqmax = 0

! Loop over all IW grid cells where cumulus parameterization may be done

   call psub()
!----------------------------------------------------------------------
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_IN,rst_buf)
#endif

   !$omp parallel do private (iw,mrlw)
   do j = 1,jtab_w(15)%jend(1); iw = jtab_w(15)%iw(j) ! jend(1) for mrl = 1
!----------------------------------------------------------------------
   call qsub('W',iw)

! MRL for current IW column

      mrlw = itab_w(iw)%mrlw

! Select cumulus convection scheme based on MRL of current IW column

      if (nqparm(mrlw) == 1) then
   
! Kuo deep convection

         call cuparm_kuo(iw)

      elseif (nqparm(mrlw) == 2) then
   
! Grell deep convection

         call cuparth(iw,dtlong)

      endif
   
      if (nqparm_sh(mrlw) == 1) then
   
! Grell shallow convection

         call cuparth_shal(iw,dtlong)
      
      endif

! Update indices for maximum heating rate
! (commented out because incorrect in parallel operation)

!   do k = lpw(iw),mza-1
!      ftcon = thsrc(k,iw) / rho(k,iw)
!      if (ftcon > dthmax) then
!         dthmax = ftcon
!         iwqmax = iw
!         kqmax = k
!      endif
!   enddo

   enddo
   !$omp end parallel do
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_OUT,rst_buf)
#endif
   call rsub('Wa',15)

! Print maximum heating rate
! (commented out because incorrect in parallel operation)

!   write(io6, '(A,I0,A,I0,A,F0.3,A)') " MAX CONVECTIVE HEATING RATE AT IW=",  &
!        iwqmax, " K=", kqmax, " IS ", dthmax*86400., " K/DAY"


endif

! Add current value of convective tendencies to thilt and sh_wt arrays
! every long timestep (whether cumulus parameterization is updated 
! this timestep or not)

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_IN,rst_buf)
#endif
!$omp parallel do private (iw,k)
do j = 1,jtab_w(15)%jend(mrl); iw = jtab_w(15)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   aconpr(iw) = aconpr(iw) + conprr(iw) * dtlong

   do k = lpw(iw),mza-1
      thilt(k,iw) = thilt(k,iw)  &
                  + (thsrc(k,iw) + thsrcsh(k,iw)) * rho(k,iw) * volt(k,iw)

      sh_wt(k,iw) = sh_wt(k,iw)  &
                  + (rtsrc(k,iw) + rtsrcsh(k,iw)) * rho(k,iw)

      rhot (k,iw) = rhot (k,iw)  &
                  + (rtsrc(k,iw) + rtsrcsh(k,iw)) * rho(k,iw) * volt(k,iw)
   enddo

enddo
!$omp end parallel do
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_OUT,rst_buf)
#endif
endif
call rsub('Wb',15)

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_CUPARM_DRIVER_OUT,rst_buf)
#endif

return
end subroutine cuparm_driver
