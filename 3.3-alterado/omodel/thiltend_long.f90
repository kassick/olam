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
subroutine thiltend_long(alpha_press,rhot)

use mem_ijtabs, only: istp, jtab_w, mrl_begl
use mem_grid,   only: mza, mwa
use misc_coms,  only: io6
use rastro_evts

!$ use omp_lib

implicit none

real, intent(out) :: alpha_press(mza,mwa)
real, intent(inout) :: rhot(mza,mwa)

integer :: j,iw,mrl

#ifdef OLAM_RASTRO
character*1 :: rst_buf = '_'
call rst_event_s_f(OLAM_THILTEND_LONG_IN,rst_buf)
#endif


! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_IN,rst_buf)
#endif
!$omp parallel do private (iw)
do j = 1,jtab_w(17)%jend(mrl); iw = jtab_w(17)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call thiltend_long0(iw,alpha_press,rhot)

enddo
!$omp end parallel do
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PARBLOCK_OUT,rst_buf)
#endif
endif
call rsub('W',17)

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_THILTEND_LONG_OUT,rst_buf)
#endif

return
end subroutine thiltend_long

!===============================================================================

subroutine thiltend_long0(iw,alpha_press,rhot)

use mem_ijtabs,  only: itab_w
use misc_coms,   only: io6, dtlm
use mem_basic,   only: rho, thil, sh_w, sh_v, theta
use mem_turb,    only: vkh, hkm, sxfer_tk, sxfer_rk
use mem_tend,    only: thilt
use consts_coms, only: pc1, rdry, rvap, cpocv
use mem_grid,    only: mza, mwa, lpw, lsw, arw, dzim, volt, dniu, aru

implicit none

integer, intent(in) :: iw

real, intent(out) :: alpha_press(mza,mwa)
real, intent(inout) :: rhot(mza,mwa)

integer :: k,ka,ks
integer :: iu1,iu2,iu3
integer :: iw1,iw2,iw3

real :: df,hdxiu,c1,c1volti,tw,wt1,hdyiv,zmkf,flux,dtl,dtli
real :: hdniu1,hdniu2,hdniu3

! Automatic arrays:

real, dimension(mza) :: akodz,dtomass,vctr5  &
                       ,vctr6,vctr7,vctr8,vctr9,del_thil

iu1 = itab_w(iw)%iu1; iu2 = itab_w(iw)%iu2; iu3 = itab_w(iw)%iu3
iw1 = itab_w(iw)%iw1; iw2 = itab_w(iw)%iw2; iw3 = itab_w(iw)%iw3

dtl = dtlm(itab_w(iw)%mrlw)
dtli = 1. / dtl
ka = lpw(iw)

! Vertical loop over T levels 

do k = ka,mza-1
   dtomass(k) = dtl / (rho(k,iw) * volt(k,iw)) 
enddo

! Vertical loop over W levels: fill tri-diagonal matrix coefficients and r.h.s.

do k = ka,mza-2
   akodz(k) = arw(k,iw) * vkh(k,iw) * dzim(k)
   vctr5(k) = - akodz(k) * dtomass(k)
   vctr7(k) = - akodz(k) * dtomass(k+1)
   vctr6(k) = 1. - vctr5(k) - vctr7(k)
   vctr8(k) = akodz(k) * (thil(k,iw) - thil(k+1,iw))
enddo

! Vertical loop over T levels that are adjacent to surface

do ks = 1,lsw(iw)
   k = ka + ks - 1

! Apply surface heat xfer [kg_a K] directly to thilt [kg_a K / s]
   
   thilt(k,iw) = thilt(k,iw) + dtli * sxfer_tk(ks,iw)

! Apply surface vapor xfer [kg_vap] directly to rhot [kg_air / s]

   rhot(k,iw) = rhot(k,iw) + dtli * sxfer_rk(ks,iw)

! Change in thil from surface heat xfer

   del_thil(k) = sxfer_tk(ks,iw) / (rho(k,iw) * volt(k,iw))
   
! Zero out sxfer_tk(ks,iw) now that it has been transferred to the atm

   sxfer_tk(ks,iw) = 0.  
enddo

! Lowest T level that is not adjacent to surface

del_thil(ka+lsw(iw)) = 0.

! Vertical loop over W levels that are adjacent to surface

do ks = 1,lsw(iw)
   k = ka + ks - 1

! Change in vctr8 from surface heat xfer

   vctr8(k) = vctr8(k) + akodz(k) * (del_thil(k) - del_thil(k+1))
enddo

! Solve tri-diagonal matrix

call tridiffo(mza,ka,mza-2,vctr5,vctr6,vctr7,vctr8,vctr9)

! Set bottom and top internal turbulent fluxes to zero

vctr9(ka-1) = 0.
vctr9(mza-1) = 0.

hdniu1 = .5 * dniu(iu1)  ! use this 1/dx form now - it seems better than A/V
hdniu2 = .5 * dniu(iu2)  ! use this 1/dx form now - it seems better than A/V
hdniu3 = .5 * dniu(iu3)  ! use this 1/dx form now - it seems better than A/V

! Vertical loop over T levels

do k = ka,mza-1

! Update thil tendency from vertical and horizontal turbulent fluxes

   thilt(k,iw) = thilt(k,iw)                     &
      
      + vctr9(k-1) - vctr9(k)                    & ! vertical diffusive fluxes

      + aru(k,iu1) * (hkm(k,iw1) + hkm(k,iw))    & ! horiz diff across ARU1
          * hdniu1 * (thil(k,iw1) - thil(k,iw))  & ! horiz diff across ARU1

      + aru(k,iu2) * (hkm(k,iw2) + hkm(k,iw))    & ! horiz diff across ARU2
          * hdniu2 * (thil(k,iw2) - thil(k,iw))  & ! horiz diff across ARU2

      + aru(k,iu3) * (hkm(k,iw3) + hkm(k,iw))    & ! horiz diff across ARU3
          * hdniu3 * (thil(k,iw3) - thil(k,iw))    ! horiz diff across ARU3

! Compute the ALPHA term for pressure used in prdctw2

! old form:

!     alpha_press(k,iw) = pc1                              &
!        * ((1. - sh_w(k,iw)) * rdry + sh_v(k,iw) * rvap)  &
!        * theta(k,iw) / thil(k,iw)

   alpha_press(k,iw) = pc1                               &
      * (((1. - sh_w(k,iw)) * rdry + sh_v(k,iw) * rvap)  &
      * theta(k,iw) / thil(k,iw)) ** cpocv

enddo

return
end subroutine thiltend_long0
