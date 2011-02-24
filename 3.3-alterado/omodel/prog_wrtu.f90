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
subroutine prog_wrtu(umarusc,wmarwsc,alpha_press,rhot)

use mem_ijtabs, only: jtab_u, jtab_w, itab_u, istp, itab_w,  &
                      mrl_begs, mrl_ends, mrl_endl
use mem_basic,  only: rho, thil, wc, press, wmc, ump, umc, uc
use mem_grid,   only: zt, zm, dzim, lpw, mza, mua, mwa, aru, lpu, lcu,  &
                      volui, volt, unx, volwi, dzm, dzt
use misc_coms,  only: io6, iparallel, cnum_vel, cnum_sclr

!$ use omp_lib

!------------------------------------------
! Only for ncar test cases:
!use ncar_testcases_all, only: ncar_testcase
!------------------------------------------

implicit none

real, intent(inout) :: umarusc    (mza,mua)
real, intent(inout) :: wmarwsc    (mza,mwa)
real, intent(in   ) :: alpha_press(mza,mwa)
real, intent(in   ) :: rhot       (mza,mwa)

integer :: j,iu,iw,k,ka,kb,iwp,mrl,iup

integer :: iw1,iw2

real :: dts,dts2,wmarw2,cnvel,pgf

! automatic arrays

real(kind=8) :: umaru(mza,mua) ! U mom density times U-face surface area [kg/s]
real(kind=8) :: wmarw(mza,mwa) ! W mom density times W-face surface area [kg/s]
!real :: umaru(mza,mua) ! U mom density times U-face surface area [kg/s]
!real :: wmarw(mza,mwa) ! W mom density times W-face surface area [kg/s]

real :: hcnum_u(mza,mua)
real :: hcnum_w(mza,mua)
real :: hflux_t(mza,mua)

real(kind=8) :: rho_s(mza,mwa)

real :: vadvwt1(mza)
real :: vadvwt2(mza)

real, parameter :: chi = .1
real, parameter :: chic = 1.5 + chi
real, parameter :: chip = -.5 - chi

! Make copy of rho array

rho_s(:,:) = rho(:,:)

! Compute vertical advective weights for stretched grid

do k = 1,mza-1
   vadvwt1(k) = (zt(k+1) - zm(k)) * dzim(k)
   vadvwt2(k) =  (zm(k) - zt(k))  * dzim(k)
enddo

! Horizontal mass fluxes for acoustic and long timesteps, and ump update

!$omp parallel do private(k) 
do iu = 2,mua
   do k = 1,mza-1
      umaru  (k,iu) = (chic * umc(k,iu) + chip * ump(k,iu)) * aru(k,iu)
      umarusc(k,iu) = umarusc(k,iu) + umaru(k,iu)
      ump    (k,iu) = umc(k,iu)
   enddo
enddo
!$omp end parallel do 

! Horizontal fluxes or courant numbers for acoustic timestep quantities

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iu) 
do j = 1,jtab_u(15)%jend(mrl); iu = jtab_u(15)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)

   call hflux(iu,rho_s,umaru,hcnum_u,hcnum_w,hflux_t)

enddo
!$omp end parallel do 
endif
call rsub('U',15)

! Send/recv horizontal fluxes or cnums

if (iparallel == 1) then
   call mpi_send_uf(hcnum_u,hcnum_w,hflux_t)
   call mpi_recv_uf(hcnum_u,hcnum_w,hflux_t)
endif

!-------------------------------------
! Only for ncar test cases:
!   if (ncar_testcase == 3) return
!-------------------------------------

! Main loop over W columns for updating WM, WC, WMARW, RHO, THIL, and PRESS

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw) 
do j = 1,jtab_w(19)%jend(mrl); iw = jtab_w(19)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call prog_wrt(iw,umaru,wmarw,wmarwsc,alpha_press,rhot  &
      ,hflux_t,hcnum_w,rho_s,vadvwt1,vadvwt2)
      
enddo
!$omp end parallel do 
endif
call rsub('Wa',19)

! LATERAL BOUNDARY CONDITION ONLY FOR LIMITED-AREA-DOMAIN MODEL CONFIGURATION
! Copy LBC for PRESS, RHO, WMARW

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw,iwp,k) 
do j = 1,jtab_w(22)%jend(mrl); iw = jtab_w(22)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
call qsub('W',iw)

   do k = 1,mza
      press(k,iw) = press(k,iwp)
      rho(k,iw) = rho(k,iwp)
      wmarw(k,iw) = wmarw(k,iwp)
   enddo
   
enddo
!$omp end parallel do 
endif
call rsub('W',22)

if (iparallel == 1) then
   call mpi_send_w('P',wmarw)  ! Send P group
   call mpi_recv_w('P',wmarw)  ! Recv P group
endif

! Horizontal loop over U points to update UMC

call psub()
!----------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Uncomment for 10-meter mountain experiment only !!!!
!   call uwcomp(0,0,0,0.,0.,0.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp parallel do private(iu)
do j = 1,jtab_u(16)%jend(mrl); iu = jtab_u(16)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)

   call prog_u(iu,umaru,wmarw,hcnum_u,rho_s)

enddo
!$omp end parallel do 
endif
call rsub('Ua',16)

! WC diagnosis from WMC and RHO - Form from Wenneker    

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw,ka,k) 
do j = 1,jtab_w(19)%jend(mrl); iw = jtab_w(19)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   ka = lpw(iw)
   do k = ka,mza-2
      wc(k,iw) = 2. * wmc(k,iw) * dzm(k)  &
               / (dzt(k+1) * rho(k,iw) + dzt(k) * rho(k+1,iw))
   enddo

   wc(1:ka-1,iw) = wc(ka,iw)  ! WC topography BC
   wc(mza-1,iw) = 0.
enddo
!$omp end parallel do 
endif
call rsub('Wb',19)

! LATERAL BOUNDARY CONDITION ONLY FOR LIMITED-AREA-DOMAIN MODEL CONFIGURATION
! Copy LBC for WMC, WC, THIL

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw,iwp,k) 
do j = 1,jtab_w(24)%jend(mrl); iw = jtab_w(24)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
call qsub('W',iw)

   do k = 1,mza
      wmc (k,iw) = wmc (k,iwp)
      wc  (k,iw) = wc  (k,iwp)
      thil(k,iw) = thil(k,iwp)
   enddo

enddo
!$omp end parallel do 
endif
call rsub('W',24)

! UC diagnosis from UMC and RHO - Form from Wenneker      

call psub()
!----------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!$omp parallel do private(iu,iw1,iw2,ka,k)
do j = 1,jtab_u(16)%jend(mrl); iu = jtab_u(16)%iu(j)
iw1 = itab_u(iu)%iw1; iw2 = itab_u(iu)%iw2
!----------------------------------------------------------------------
call qsub('U',iu)

   ka = lcu(iu)
   do k = ka,mza-1
      uc(k,iu) = umc(k,iu) / (volui(k,iu)  &
         * (volt(k,iw2) * rho(k,iw1) + volt(k,iw1) * rho(k,iw2)))
   enddo

   uc(1:ka-1,iu) = uc(ka,iu)

enddo
!$omp end parallel do 
endif
call rsub('Ub',16)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Uncomment for 10-meter mountain experiment only !!!!
!   call uwcomp(2,0,0,0.,0.,0.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! LATERAL BOUNDARY CONDITION ONLY FOR LIMITED-AREA-DOMAIN MODEL CONFIGURATION
! Copy LBC for UMC, UC

call psub()
!----------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!$omp parallel do private(iu,iup,k) 
do j = 1,jtab_u(18)%jend(mrl); iu = jtab_u(18)%iu(j)
   iup = itab_u(iu)%iup
!----------------------------------------------------------------------
call qsub('U',iu)

   do k = 1,mza-1
      umc(k,iu) = umc(k,iup)
      uc(k,iu)  = uc(k,iup)
   enddo

enddo
!$omp end parallel do 
endif
call rsub('U',18)

if (iparallel == 1) then
   call mpi_send_u('U')  ! Send U group
endif

return
end subroutine prog_wrtu

!=========================================================================

subroutine hflux(iu,rho_s,umaru,hcnum_u,hcnum_w,hflux_t)

use mem_ijtabs,  only: itab_u
use mem_basic,   only: thil
use misc_coms,   only: io6, cnum_sclr, cnum_vel, dtsm
use mem_grid,    only: mza, mua, mwa, lpu, dts_w, volt

implicit none

integer, intent(in) :: iu

real(kind=8), intent(in ) :: umaru  (mza,mua)
!real, intent(in ) :: umaru  (mza,mua)
real, intent(out) :: hcnum_u(mza,mua)
real, intent(out) :: hcnum_w(mza,mua)
real, intent(out) :: hflux_t(mza,mua)

real(kind=8), intent(in) :: rho_s (mza,mwa)

integer :: iw1,iw2,iw3,iw4,iw5,iw6

integer :: k,kb

real :: dts,dts2
real :: hcnum_s
real :: cnum_w

real, parameter :: cnu1 = .4
real, parameter :: cnu2 = .2

real, parameter :: cnw1 = .4
real, parameter :: cnw2 = .2

real, parameter :: cns1 = .4
real, parameter :: cns2 = .2

! Automatic array

real :: umass(mza)

iw1 = itab_u(iu)%iw1; iw2 = itab_u(iu)%iw2; iw3 = itab_u(iu)%iw3
iw4 = itab_u(iu)%iw4; iw5 = itab_u(iu)%iw5; iw6 = itab_u(iu)%iw6

dts = dtsm(itab_u(iu)%mrlu)
dts2 = dts * 2.

kb = lpu(iu)

! Loop over T levels

hcnum_u(1:kb-1,iu) = 0.
hflux_t(1:kb-1,iu) = 0.

umass(kb-1) = 0.0
do k = kb,mza-1

! Mass in combined iw1 and iw2 cells

   umass(k) = rho_s(k,iw1) * volt(k,iw1) + rho_s(k,iw2) * volt(k,iw2)

! Horizontal courant numbers for U and S (double value for S 
! compensates for double mass in denominator since flux used for T cell) 

   hcnum_u(k,iu) = umaru(k,iu) * dts / umass(k)
   hcnum_s       = 2. * hcnum_u(k,iu)

! Impose minimum values on hcnum_u and hcnum_s

   if (cnum_vel > abs(hcnum_u(k,iu)))  &
      hcnum_u(k,iu) = sign(cnum_vel,hcnum_u(k,iu))

   if (cnum_sclr > abs(hcnum_s))  &
      hcnum_s = sign(cnum_sclr,hcnum_s) 

!------------------------------------------------------------------------
! Simple increase in hcnum_u and hcnum_s near topography.  Later, increase 
! more selectively using monotonic scheme based on local & nearby field values

   go to 1

   if (k <= kb) then
      hcnum_u(k,iu) = sign(cnu1,hcnum_u(k,iu)) 
      hcnum_s       = sign(cns1,hcnum_s) 
   endif
   
   if (k == kb+1) then
      hcnum_u(k,iu) = sign(cnu2,hcnum_u(k,iu)) 
      hcnum_s       = sign(cns2,hcnum_s) 
   endif
1  continue
!------------------------------------------------------------------------

! Horizontal flux of thil

   hflux_t(k,iu) = umaru(k,iu) * .5 * (thil(k,iw1) + thil(k,iw2)    & 
                          + hcnum_s * (thil(k,iw1) - thil(k,iw2)))

enddo

! Loop over W levels (of U column) for horizontal W Courant number

hcnum_w(1:kb-2,iu) = 0.

do k = kb-1,mza-2
   hcnum_w(k,iu) = (umaru(k,iu) + umaru(k+1,iu)) * dts2  &
                 / (umass(k) + umass(k+1))

! Impose minimum value on hcnum_w

   if (cnum_vel > abs(hcnum_w(k,iu)))  &
      hcnum_w(k,iu) = sign(cnum_vel,hcnum_w(k,iu)) 

!------------------------------------------------------------------------
! Simple increase in hcnum_w near topography.  Later, increase more selectively
! using monotonic scheme based on WC at (iw1,iw2,iw3,iw4,iw5,iw6).

   go to 2
   if (k <= kb) then
      hcnum_w(k,iu) = sign(cnw1,hcnum_w(k,iu)) 
   endif
   
   if (k == kb+1) then
      hcnum_w(k,iu) = sign(cnw2,hcnum_w(k,iu)) 
   endif
2  continue
!------------------------------------------------------------------------

enddo

return
end subroutine hflux

!=========================================================================

subroutine prog_wrt(iw,umaru,wmarw,wmarwsc,alpha_press,rhot  &
   ,hflux_t,hcnum_w,rho_s,vadvwt1,vadvwt2)

use mem_tend,    only: thilt, wmt
use mem_ijtabs,  only: itab_w
use mem_basic,   only: wmc, rho, thil, wc, uc, theta, press
use misc_coms,   only: io6, cnum_sclr, cnum_vel, initial, dn01d, th01d, &
                       deltax, nxp, mdomain, time8
use consts_coms, only: cpocv, gravo2, grav
use mem_grid,    only: mza, mua, mwa, lpw, dts_w, arw, volt, volti, volwi, dzt, &
                       dzim, xew, zm
use mem_rayf,    only: rayfw_distim, rayf_cofw, rayf_distim, rayf_cof

implicit none

integer, intent(in) :: iw

real(kind=8), intent(in   ) :: umaru      (mza,mua)
real(kind=8), intent(out  ) :: wmarw      (mza,mwa)
!real, intent(in   ) :: umaru      (mza,mua)
!real, intent(out  ) :: wmarw      (mza,mwa)
real, intent(inout) :: wmarwsc    (mza,mwa)
real, intent(in   ) :: alpha_press(mza,mwa)
real, intent(in   ) :: rhot       (mza,mwa)
real, intent(in   ) :: hflux_t    (mza,mua)
real, intent(in   ) :: hcnum_w    (mza,mua)

real(kind=8), intent(in) :: rho_s (mza,mwa)

real, intent(in) :: vadvwt1 (mza)
real, intent(in) :: vadvwt2 (mza)

integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9
integer :: iw1,iw2,iw3

real :: diru1,diru2,diru3
real :: fwu4,fwu5,fwu6,fwu7,fwu8,fwu9
real :: fww1,fww2,fww3

integer :: k,ka
integer :: k1,k2,k3

real :: dts,dtso2,dts2,dts8,flux_rhothil

real :: c6,c7,c8,c9,c10,wmarw2  &
   ,pfrac,sinlat,coslat,theq
real :: fracx, rayfx

real :: ak = 1. / (40. * 86400.)  ! needed only for HS experiment
real :: as = 1. / (4.  * 86400.)  ! needed only for HS experiment

! Dirichlet nudging inverse time scale for W

real :: bctimei = 1. / 1000.
   
real :: cnum_w
real :: vproj1,vproj2,vproj3
real :: del_rhothil
real :: hcnsclr

! Vertical implicit scheme weighting parameters

real, parameter :: fw = .55  & ! wmc   ( .5 was std)  later .85
                  ,fr = .55  & ! rho   ( .5 was std)  later .50
                  ,fp = .75    ! press (.75 was std)  later .85

real, parameter :: pc2 = fp * cpocv

real :: sc7,sc8,sc9,sc11,sc12,sc13,sc14,sc21,sc22,sc23,sc24,sc25,sc26
                  
! Automatic arrays

real :: delex_wm      (mza)
real(kind=8) :: delex_rho     (mza)
real :: delex_rhothil (mza)
real :: del_wm        (mza)
real(kind=8) :: del_wmarw     (mza)
real :: wmtharw       (mza)
real :: del_wmtharw   (mza)
real :: wvertflx      (mza)

real(kind=8) :: rhothil (mza)
real(kind=8) :: press_t (mza)

real, dimension(mza) :: vctr1,vctr2,vctr5,vctr6,vctr10
real, dimension(mza) :: vctr31,vctr32,vctr33,vctr34,vctr36

real, parameter :: hcns1 = .2
real, parameter :: hcns2 = .1

real, parameter :: cnw1 = .4
real, parameter :: cnw2 = .2


iu1 = itab_w(iw)%iu1; iu2 = itab_w(iw)%iu2; iu3 = itab_w(iw)%iu3
iu4 = itab_w(iw)%iu4; iu5 = itab_w(iw)%iu5; iu6 = itab_w(iw)%iu6
iu7 = itab_w(iw)%iu7; iu8 = itab_w(iw)%iu8; iu9 = itab_w(iw)%iu9

iw1 = itab_w(iw)%iw1; iw2 = itab_w(iw)%iw2; iw3 = itab_w(iw)%iw3

diru1 = itab_w(iw)%diru1; diru2 = itab_w(iw)%diru2; diru3 = itab_w(iw)%diru3

fwu4 = itab_w(iw)%fwu4;  fwu5 = itab_w(iw)%fwu5; fww1 = itab_w(iw)%fww1
fwu6 = itab_w(iw)%fwu6;  fwu7 = itab_w(iw)%fwu7; fww2 = itab_w(iw)%fww2
fwu8 = itab_w(iw)%fwu8;  fwu9 = itab_w(iw)%fwu9; fww3 = itab_w(iw)%fww3

ka = lpw(iw)

dts  = dts_w(iw)
dtso2 = .5 * dts
dts2 = 2. * dts
dts8 = 8. * dts

! Set bottom & top advective mass and heat fluxes to zero

wmarw(1:ka-1,iw) = 0.
wmarw(mza-1,iw) = 0.

wmtharw(ka-1)  = 0.
wmtharw(mza-1) = 0.

! Loop over W levels

do k = ka,mza-2

! vertical mass flux for time level t

   wmarw(k,iw) = wmc(k,iw) * arw(k,iw)

! half vertical courant number for scalars

   hcnsclr = wmarw(k,iw) * dts * volwi(k,iw) / (rho_s(k,iw) + rho_s(k+1,iw))
   if (.5*cnum_sclr > abs(hcnsclr)) hcnsclr = sign(.5*cnum_sclr,hcnsclr)

 go to 1
   if (k == ka)   hcnsclr = sign(hcns1,hcnsclr)
   if (k == ka+1) hcnsclr = sign(hcns2,hcnsclr)
1  continue

! vertical heat flux for time level t

   wmtharw(k) = wmarw(k,iw)      &
      * (vadvwt1(k) * thil(k,iw) + vadvwt2(k) * thil(k+1,iw)  &
      + hcnsclr * (thil(k,iw) - thil(k+1,iw)))
enddo

! Loop over T levels

do k = ka,mza-1

! Change in RHO from horizontal and EXPLICIT vertical advection

   delex_rho(k) = dts * volti(k,iw)  &
      * (rhot(k,iw)                  & ! long timestep tend (from nudging)
      + diru1 * umaru(k,iu1)         &
      + diru2 * umaru(k,iu2)         &
      + diru3 * umaru(k,iu3)         &
      + wmarw(k-1,iw) - wmarw(k,iw)  )
      
! Change in RHOTHIL from long timestep tendencies plus
!   horizontal and EXPLICIT vertical advection

   delex_rhothil(k) = dts * volti(k,iw)  &
      * (thilt(k,iw)                     & ! long timestep tend
      +  diru1 * hflux_t(k,iu1)          &
      +  diru2 * hflux_t(k,iu2)          &
      +  diru3 * hflux_t(k,iu3)          &
      + wmtharw(k-1) - wmtharw(k)        ) ! vertical advection

! RHOTHIL(t) and PRESS(t)
! alpha_press is evaluated on long timestep in thiltend_long

   rhothil(k) = rho_s(k,iw) * thil(k,iw)
   press_t(k) = alpha_press(k,iw) * rhothil(k) ** cpocv

! double vertical mass flux averaged vertically to T points

   wmarw2 = wmarw(k,iw) + wmarw(k-1,iw) ! Use unequal wts on stretched grid?

! vertical courant number for vertical WC advection (double mass flux in numerator
! compensates for double mass in denominator since flux used for W cell)

   cnum_w = wmarw2 * dtso2 / (rho_s(k,iw) * volt(k,iw))
   if (cnum_vel > abs(cnum_w)) cnum_w = sign(cnum_vel,cnum_w)

 go to 2
   if (k == ka)   cnum_w = sign(cnw1,cnum_w)
   if (k == ka+1) cnum_w = sign(cnw2,cnum_w)
2  continue

! vertical advective flux of WC (use /= wts?)

   wvertflx(k) = .25 * wmarw2 * (wc(k-1,iw) + wc(k,iw)  &
                     + cnum_w * (wc(k-1,iw) - wc(k,iw)))

enddo

   wvertflx(ka) = 0.

! Loop over W levels for update of DELEX_WM (downward loop to handle ka-1 level)

do k = mza-2,ka-1,-1

   k1 = max(k,lpw(iw1)-1)
   k2 = max(k,lpw(iw2)-1)
   k3 = max(k,lpw(iw3)-1)

! Projected velocity transported by IU1

   vproj1 = fwu4 * (uc(k1,iu4) + uc(k1+1,iu4))  &
          + fwu5 * (uc(k1,iu5) + uc(k1+1,iu5))  &
          + fww1 * wc(k1,iw1)

! Projected velocity transported by IU2

   vproj2 = fwu6 * (uc(k2,iu6) + uc(k2+1,iu6))  &
          + fwu7 * (uc(k2,iu7) + uc(k2+1,iu7))  &
          + fww2 * wc(k2,iw2)

! Projected velocity transported by IU3

   vproj3 = fwu8 * (uc(k3,iu8) + uc(k3+1,iu8))  &
          + fwu9 * (uc(k3,iu9) + uc(k3+1,iu9))  &
          + fww3 * wc(k3,iw3)

   if (k >= ka) then

! Change in WM from EXPLICIT terms (long timestep tendency, 3 horizontal
! advective fluxes, 2 vertical advective fluxes, vertical pgf, gravity)

      delex_wm(k) = dts * (volwi(k,iw) * (  &

           wmt(k,iw)  + .25 *              &

          ((umaru(k,iu1) + umaru(k+1,iu1)) * (diru1 * (vproj1 + wc(k,iw))   &
                                   + hcnum_w(k,iu1) * (vproj1 - wc(k,iw)))  &

         + (umaru(k,iu2) + umaru(k+1,iu2)) * (diru2 * (vproj2 + wc(k,iw))   &
                                   + hcnum_w(k,iu2) * (vproj2 - wc(k,iw)))  &
             
         + (umaru(k,iu3) + umaru(k+1,iu3)) * (diru3 * (vproj3 + wc(k,iw))   &
                                   + hcnum_w(k,iu3) * (vproj3 - wc(k,iw)))) &

         + wvertflx(k) - wvertflx(k+1))          &

         + dzim(k) * (press_t(k) - press_t(k+1)  &

         - gravo2 * (dzt(k) * rho_s(k,iw) + dzt(k+1) * rho_s(k+1,iw))))

   else
   
!!!!!!!!!!!!!!!!!!!!!special   
   go to 55
!!!!!!!!!!!!!!!!!!!!!! end special

! Change in WM(ka) from EXPLICIT terms (3 horizontal advective fluxes at ka-1)

      delex_wm(ka) = delex_wm(ka) + dts * volwi(ka,iw) * .25 *              &

          ((umaru(k,iu1) + umaru(k+1,iu1)) * (diru1 * (vproj1 + wc(k,iw))   &
                                   + hcnum_w(k,iu1) * (vproj1 - wc(k,iw)))  &

         + (umaru(k,iu2) + umaru(k+1,iu2)) * (diru2 * (vproj2 + wc(k,iw))   &
                                   + hcnum_w(k,iu2) * (vproj2 - wc(k,iw)))  &
             
         + (umaru(k,iu3) + umaru(k+1,iu3)) * (diru3 * (vproj3 + wc(k,iw))   &
                                   + hcnum_w(k,iu3) * (vproj3 - wc(k,iw)))  )

!!!!!!!!!!!!!! special
55 continue
!!!!!!!!!!!!!!!!end special
   endif
enddo

c6  = dts * .5 * fw
c7  = dts * .25 * fw
c8  = dts * fp * cpocv
c9  = dts * (-.5) * fr * grav
c10 = dts * fw

! SPECIAL - EXTRA RAYF AT ENDS OF CYCLIC DOMAIN

fracx = abs(xew(iw)) / (real(nxp-1) * .866 * deltax)

! END SPECIAL

! RAYLEIGH FRICTION ON WM

if (rayfw_distim > 1.e-6) then
    rayfx = .2 * (-2. + 3. * fracx) * rayf_cofw(mza-2)
    rayfx = 0.   ! Default: no extra RAYF
   do k = ka,mza-2
      wmc(k,iw) = wmc(k,iw) - dts * max(rayf_cofw(k),rayfx) * wmc(k,iw)
   enddo
endif

! RAYLEIGH FRICTION ON THIL

if (rayf_distim > 1.e-6) then
   if (initial == 1) then   ! HHI case
      rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza-1)
      rayfx = 0.   ! Default: no extra RAYF
      do k = ka,mza-1
! Form based on theta alone
         delex_rhothil(k) = delex_rhothil(k) + max(rayf_cof(k),rayfx)  & 
            * dn01d(k) * (th01d(k) - theta(k,iw))
      enddo
   else                     ! LHI/VARI case
      ! Need implementation for LHI/VARI (use vartp for merid. variation?)
   endif
endif

!!!!!!!!!!!!!!! Special - uncomment only for Held-Suarez experiment !!!!!!!

!hs   do k = ka,mza-1                                                     !hs
                     
!hs      pfrac = 1.e-5 * press(k,iw)                                      !hs
!hs      sinlat = sin(glatw(iw)*pio180)                                   !hs
!hs      coslat = cos(glatw(iw)*pio180)                                   !hs
!hs      theq = max(200.*pfrac**(-.285714)  &                             !hs
!hs             ,315. - 60. * sinlat**2 - 10. * log(pfrac) * coslat**2)   !hs

!hs      rayf_cof(k) = ak  &                                              !hs
!hs         + (as - ak) * max(0.,(press(k,iw)-7.e4)/3.e4) * coslat**4     !hs
            
!hs      delex_rhothil(k) = delex_rhothil(k) + dts * rayf_cof(k)  &       !hs
!hs         * rho_s(k,iw) * (theq - theta(k,iw))                            !hs
               
!hs   enddo                                                               !hs
!!!!!!!!!!!!!!!!!!!!!! End special !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Fill matrix coefficients for implicit update of WM

do k = ka,mza-1
   vctr1(k)  = wc(k,iw) + wc(k-1,iw)     ! T pts
   vctr2(k)  = thil(k,iw) + thil(k+1,iw) ! W pts
   vctr5(k)  = press_t(k) / rhothil(k)   ! T pts
   vctr6(k)  = c6 * volti(k,iw)          ! T pts
   vctr10(k) = c10 * volti(k,iw)         ! T pts
enddo
vctr2(ka-1) = vctr2(ka)

do k = ka,mza-2
         
   sc7  = c7   * volwi(k,iw) ! W pts
   sc8  = c8   * dzim(k)     ! W pts
   sc9  = c9   * dzim(k)     ! W pts
   sc11 = sc8  * vctr5(k)    ! W pts
   sc12 = sc8  * vctr5(k+1)  ! W pts
   sc13 = sc9  * dzt(k)      ! W pts
   sc14 = sc9  * dzt(k+1)    ! W pts
      
   sc21 = sc7  * vctr1(k)    ! W pts
   sc22 = sc7  * vctr1(k+1)  ! W pts
   sc23 = sc11 * vctr6(k)    ! W pts
   sc24 = sc12 * vctr6(k+1)  ! W pts
   sc25 = sc13 * vctr10(k)   ! W pts
   sc26 = sc14 * vctr10(k+1) ! W pts
      
   vctr32(k) = 1. + arw(k,iw)  &
             * (sc22 - sc21 + vctr2(k) * (sc23 + sc24) + sc25 - sc26)
         
   vctr31(k) = - arw(k-1,iw) * (sc21 + sc23 * vctr2(k-1) + sc25)
      
   vctr33(k) =   arw(k+1,iw) * (sc22 - sc24 * vctr2(k+1) + sc26) 
     
   vctr34(k) = delex_wm(k)  &
             + sc11 * delex_rhothil(k) - sc12 * delex_rhothil(k+1)  &
             + sc13 * delex_rho(k)     + sc14 * delex_rho(k+1)

enddo

! Solve implicit tri-diagonal matrix equation for delta WM (del_wm)

if (ka <= mza-2) then
   call tridiffo(mza,ka,mza-2,vctr31,vctr32,vctr33,vctr34,del_wm)
endif
       
! Vertical loop over W points

do k = ka,mza-2

! Update WMC from (t) to (t+1) due to change in WM (del_wm)

   wmc(k,iw) = wmc(k,iw) + del_wm(k)

! Change in vertical mass and heat fluxes from (t) to (t+fw)

   del_wmarw(k)   = fw * del_wm(k) * arw(k,iw)
   del_wmtharw(k) = del_wmarw(k) * .5 * (thil(k,iw) + thil(k+1,iw))

! Add vertical mass flux at (t+fw) to arrays for vertical UC and scalar transport

   wmarw(k,iw)   = wmarw(k,iw)   + del_wmarw(k)
   wmarwsc(k,iw) = wmarwsc(k,iw) + wmarw(k,iw)

enddo

! Set top and bottom values in mass-flux-change and heat-flux-change arrays to zero

del_wmarw(ka-1)  = 0.
del_wmarw(mza-1) = 0.

del_wmtharw(ka-1)  = 0.
del_wmtharw(mza-1) = 0.

! Vertical loop over T points

do k = ka,mza-1

! Change of rho from (t) to (t+1)

   rho(k,iw) = rho(k,iw) + delex_rho(k)  &
      + dts * volti(k,iw) * (del_wmarw(k-1) - del_wmarw(k))

! Change of rhothil from (t) to (t+1)

   del_rhothil = delex_rhothil(k)  &
      + dts * volti(k,iw) * (del_wmtharw(k-1) - del_wmtharw(k))

! Update pressure from (t) to (t+tp)

   press(k,iw) = press_t(k) + pc2 * vctr5(k) * del_rhothil
! Update thil from (t) to (t+1)

   thil(k,iw) = (rhothil(k) + del_rhothil) / rho(k,iw)

enddo

return
end subroutine prog_wrt

!============================================================================

subroutine prog_u(iu,umaru,wmarw,hcnum_u,rho_s)

use mem_tend,    only: umt
use mem_ijtabs,  only: itab_u
use mem_basic,   only: uc, wc, press, umc
use misc_coms,   only: io6, dtsm, cnum_vel, initial, mdomain, u01d, v01d, dn01d, &
                       deltax, nxp
use consts_coms, only: erad
use mem_grid,    only: lpu, lcu, volt, aru, volui, xeu, yeu, zeu,  &
                       unx, uny, unz, mza, mua, mwa
use mem_rayf,    only: rayf_distim, rayf_cof

implicit none

integer, intent(in) :: iu

real(kind=8), intent(in) :: umaru  (mza,mua)
real(kind=8), intent(in) :: wmarw  (mza,mwa)
!real, intent(in) :: umaru  (mza,mua)
!real, intent(in) :: wmarw  (mza,mwa)
real, intent(in) :: hcnum_u(mza,mua)

real(kind=8), intent(in) :: rho_s(mza,mwa)

integer :: k,ka,kb

integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9,iu10,iu11,iu12
integer :: iw1,iw2,iw3,iw4,iw5,iw6

real :: fuu5,fuu6,fuu7,fuu8,fuu9,fuu10,fuu11,fuu12
real :: fuw3,fuw4,fuw5,fuw6
real :: diru1,diru2,diru3,diru4
real :: vproj1,vproj2,vproj3,vproj4
real :: pgc12,pgc45,pgc63,pgc12b,pgc45b,pgc12c,pgc63c,pgc12d

real :: dts,dts2,wmarw2,cnvel,raxis,uv01dr,uv01dx,uv01dy,uv01dz,ucref,uc2
!hs real :: pressloc
real :: fracx, rayfx

! Automatic arrays

real :: fadvz(mza)
real :: umass(mza)
real :: umt_rayf(mza)
real :: pgf(mza)

! Extract neighbor indices and coefficients for this point in the U stencil

iu1  = itab_u(iu)%iu1;  iu2  = itab_u(iu)%iu2;  iu3  = itab_u(iu)%iu3
iu4  = itab_u(iu)%iu4;  iu5  = itab_u(iu)%iu5;  iu6  = itab_u(iu)%iu6
iu7  = itab_u(iu)%iu7;  iu8  = itab_u(iu)%iu8;  iu9  = itab_u(iu)%iu9
iu10 = itab_u(iu)%iu10; iu11 = itab_u(iu)%iu11; iu12 = itab_u(iu)%iu12

iw1 = itab_u(iu)%iw1; iw2 = itab_u(iu)%iw2; iw3 = itab_u(iu)%iw3
iw4 = itab_u(iu)%iw4; iw5 = itab_u(iu)%iw5; iw6 = itab_u(iu)%iw6

fuu5  = itab_u(iu)%fuu5;  fuu6  = itab_u(iu)%fuu6;  fuw3 = itab_u(iu)%fuw3
fuu7  = itab_u(iu)%fuu7;  fuu8  = itab_u(iu)%fuu8;  fuw4 = itab_u(iu)%fuw4
fuu9  = itab_u(iu)%fuu9;  fuu10 = itab_u(iu)%fuu10; fuw5 = itab_u(iu)%fuw5
fuu11 = itab_u(iu)%fuu11; fuu12 = itab_u(iu)%fuu12; fuw6 = itab_u(iu)%fuw6

diru1 = itab_u(iu)%diru1; diru2 = itab_u(iu)%diru2
diru3 = itab_u(iu)%diru3; diru4 = itab_u(iu)%diru4

pgc12  = itab_u(iu)%pgc12;  pgc45  = itab_u(iu)%pgc45
pgc63  = itab_u(iu)%pgc63
pgc12b = itab_u(iu)%pgc12b; pgc45b = itab_u(iu)%pgc45b
pgc12c = itab_u(iu)%pgc12c; pgc63c = itab_u(iu)%pgc63c
pgc12d = itab_u(iu)%pgc12d

dts = dtsm(itab_u(iu)%mrlu)
dts2 = dts * 2.

ka = lcu(iu)
kb = lpu(iu)

! RAYLEIGH FRICTION ON UMC

umt_rayf(:) = 0.

if (rayf_distim > 1.e-6) then

! Vertical loop over U points

   do k = kb,mza-1

      if (initial == 1) then      ! HHI case

! Must rotate reference wind to local UC orientation

         if (mdomain <= 1) then  ! Model uses "earth" coordinates
            raxis = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)  ! dist from earth axis
            uv01dr = -v01d(k) * zeu(iu) / erad  ! radially outward from axis

            uv01dx = (-u01d(k) * yeu(iu) + uv01dr * xeu(iu)) / raxis 
            uv01dy = ( u01d(k) * xeu(iu) + uv01dr * yeu(iu)) / raxis 
            uv01dz =   v01d(k) * raxis / erad 

            ucref = uv01dx * unx(iu) + uv01dy * uny(iu) + uv01dz * unz(iu)
         else
            ucref = u01d(k) * unx(iu) + v01d(k) * uny(iu)
         endif

! SPECIAL - EXTRA RAYF AT ENDS OF CYCLIC DOMAIN

         fracx = abs(xeu(iu)) / (real(nxp-1) * .866 * deltax)
         rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza-1)
         rayfx = 0.   ! Default: no extra RAYF
         
! END SPECIAL

         umt_rayf(k) = max(rayf_cof(k),rayfx) * dn01d(k) * (ucref - uc(k,iu))

      else                     ! LHI/VARI case
         ! Need implementation for LHI/VARI (use vartp for merid. variation?)
         ! Ok to do this in veltend_long???
      endif

   enddo

endif

! PRESSURE GRADIENT FORCE ON UMC

pgf(ka:kb-1) = 0.

! Vertical loop over T levels

do k = kb,mza-1

   if ((lpu(iu1) > k .or. lpu(iu4) > k) .and.  &
       (lpu(iu2) > k .or. lpu(iu3) > k)) then

! 2-point stencil

      pgf(k) = pgc12d * (press(k,iw1) - press(k,iw2))

   elseif (lpu(iu2) > k .or. lpu(iu3) > k) then

! 4-point stencil with pts 6-3

      pgf(k) = (pgc12c * (press(k,iw1) - press(k,iw2))  &
             +  pgc63c * (press(k,iw6) - press(k,iw3))  )

   elseif (lpu(iu4) > k .or. lpu(iu1) > k) then

! 4-point stencil with pts 4-5

      pgf(k) = (pgc12b * (press(k,iw1) - press(k,iw2))  &
             +  pgc45b * (press(k,iw4) - press(k,iw5))  )

   else

! 6-point stencil

      pgf(k) = (pgc12 * (press(k,iw1) - press(k,iw2))  &
             +  pgc45 * (press(k,iw4) - press(k,iw5))  &
             +  pgc63 * (press(k,iw6) - press(k,iw3))  )

   endif

enddo

! Vertical loop over T levels

do k = ka,mza-1

! Mass in UM control volume

   umass(k) = rho_s(k,iw1) * volt(k,iw1) + rho_s(k,iw2) * volt(k,iw2)

enddo

! Vertical loop over W levels

do k = ka,mza-2
   
! Vertical mass flux at top of U control volume

   wmarw2 = wmarw(k,iw1) + wmarw(k,iw2)  ! U ctrl vol top mass flux
      
! Vertical courant number (double time in numerator, double mass in denom)

   cnvel = wmarw2 * dts2 / (umass(k) + umass(k+1))
   if (cnum_vel > abs(cnvel)) cnvel = sign(cnum_vel,cnvel)

go to 3
   if (k == ka)   cnvel = sign(.9,cnvel)
   if (k == ka+1) cnvel = sign(.45,cnvel)
3  continue

! Vertical advective flux of UC

   fadvz(k) = .5 * wmarw2 * (uc(k,iu) + uc(k+1,iu)  &
                  + cnvel * (uc(k,iu) - uc(k+1,iu)) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Uncomment for 10-meter mountain experiment only !!!!
!      uc2 = uc(k,iu) + uc(k+1,iu) + cnvel * (uc(k,iu) - uc(k+1,iu))
!      call uwcomp(1,iu,k,fadvz(k),wmarw2,uc2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

enddo

! Set bottom & top advective fluxes to zero

fadvz(ka-1) = 0.
fadvz(mza-1) = 0.

! Vertical loop over U points

do k = ka,mza-1

! Projected velocities

   vproj1 = fuu5  * uc(k,iu5)               &
          + fuu6  * uc(k,iu6)               &
          + fuw3  * (wc(k,iw3) + wc(k-1,iw3))

   vproj2 = fuu7  * uc(k,iu7)               &
          + fuu8  * uc(k,iu8)               &
          + fuw4  * (wc(k,iw4) + wc(k-1,iw4))

   vproj3 = fuu9  * uc(k,iu9)               &
          + fuu10 * uc(k,iu10)              &
          + fuw5  * (wc(k,iw5) + wc(k-1,iw5))

   vproj4 = fuu11 * uc(k,iu11)              &
          + fuu12 * uc(k,iu12)              &
          + fuw6  * (wc(k,iw6) + wc(k-1,iw6))

! Update UM from long timestep tendencies, advection, and pressure gradient force

   umc(k,iu) = umc(k,iu) + dts * (pgf(k) + umt_rayf(k) + volui(k,iu)  &
      * (umt(k,iu) + .5 *                    &

       (umaru(k,iu1) * (diru1 * (vproj1 + uc(k,iu))    &
             + hcnum_u(k,iu1) * (vproj1 - uc(k,iu)))   &

      + umaru(k,iu2) * (diru2 * (vproj2 + uc(k,iu))    &
             + hcnum_u(k,iu2) * (vproj2 - uc(k,iu)))   &

      + umaru(k,iu3) * (diru3 * (vproj3 + uc(k,iu))    &
             + hcnum_u(k,iu3) * (vproj3 - uc(k,iu)))   &

      + umaru(k,iu4) * (diru4 * (vproj4 + uc(k,iu))    &
             + hcnum_u(k,iu4) * (vproj4 - uc(k,iu))))  &

      + fadvz(k-1) - fadvz(k)))

enddo

!!!!!!!!!!!! Special section (uncomment only for Held-Suarez experiment)
!hs   do k = lpu(iu),mza-1                                         !hs
!hs      pressloc = .5 * (press(k,iw1) + press(k,iw2))             !hs
!hs      rayf_cof(k) = max(0.,(pressloc-7.e4)/3.e4) / 86400.       !hs

!hs      umc(k,iu) = umc(k,iu) + dts * rayf_cof(k)  &              !hs
!hs         * .5 * (rho_s(k,iw1) + rho_s(k,iw2)) * (0. - uc(k,iu)) !hs
!hs   enddo                                                        !hs
!!!!!!!!!!!! End special section (uncomment only for Held-Suarez experiment)

return
end subroutine prog_u
