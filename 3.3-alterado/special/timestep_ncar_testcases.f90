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

!------------------------------------------
! Only for ncar test cases:
use ncar_testcases_all, only: ncar_testcase
!------------------------------------------

implicit none

integer :: is,mvl,isl4,jstp
real :: t1,w1

! automatic arrays

real :: umarusc    (mza,mua) ! U mom density times U-face surface area [kg/s]
real :: wmarwsc    (mza,mwa) ! W mom density times W-face surface area [kg/s]
real :: rho_old    (mza,mwa) ! density from beginning of timestep [kg/m^3]
real :: alpha_press(mza,mwa) ! 
real :: rhot       (mza,mwa) ! grid-cell total mass tendency [kg/s]

! +----------------------------------------------------------------------------+
! |  Each call to subroutine timestep drives all steps in advancing by dtlong  |
! +----------------------------------------------------------------------------+

time_istp8 = time8

if (time8 < 1.e-3) then
!   call bubble()
endif

do jstp = 1,nstp  ! nstp = no. of finest-grid-level aco steps in dtlm(1)
   istp = jstp

   call tend0(rhot)

!-------------------------------------
! Only for ncar test cases:
      if (jstp == 1) call diagn_global()
!-------------------------------------

!-------------------------------------
! Only for ncar test cases:
   if (ncar_testcase == 3) go to 31
!-------------------------------------

call check_nans(1)

   if (ilwrtyp + iswrtyp > 0) then
      call radiate()
   endif

call check_nans(2)

   call surface_turb_flux()

call check_nans(3)

   if (any(nqparm   (1:mrls) > 0) .or.  &
       any(nqparm_sh(1:mrls) > 0)) then

      call cuparm_driver(rhot)

      if (isfcl == 1) then
         call surface_cuparm_flux()
      endif
   endif

call check_nans(4)

   call veltend_long()

call check_nans(5)

   call thiltend_long(alpha_press,rhot)

call check_nans(6)

   if (initial == 2 .and. nudflag == 1)  &
   call obs_nudge(rhot)

call check_nans(7)

!-------------------------------------
! Only for ncar test cases:
31 continue
!-------------------------------------

   call zero_massflux(umarusc,wmarwsc,rho_old)

call check_nans(8)

!-------------------------------------------------------------------
! Most of the following subroutine is bypassed for ncar test case 3:
   call prog_wrtu(umarusc,wmarwsc,alpha_press,rhot)
!-------------------------------------------------------------------

call check_nans(9)

   call timeavg_massflux(umarusc,wmarwsc)

   call scalar_transport(umarusc,wmarwsc,rho_old)

call check_nans(10)

   call predtr(rho_old)

call check_nans(11)

   if (level /= 3) then
      call thermo()
   endif

call check_nans(12)

   if (level == 3) then
      call micro()  ! maybe later make freq. uniform

      if (isfcl == 1) then
         call surface_precip_flux()
      endif
   endif

call check_nans(13)

   call trsets()  

call check_nans(14)

   if (iparallel == 1) then
      call mpi_recv_w('W',wmarwsc)  ! Recv W group; wmarwsc used as dummy array here
      call mpi_recv_u('U')  ! Recv U group
   endif

call check_nans(15)

   call turb_k()  ! Get K's for next timestep

call check_nans(16)

   if (iparallel == 1) then
      call mpi_send_w('S',wmarwsc)  ! Send scalars; wmarwsc used as dummy array here
      call mpi_recv_w('S',wmarwsc)  ! Recv scalars; wmarwsc used as dummy array here
   endif

call check_nans(17)

   if (leafstep(istp) > 0) then
      call leaf3()

      call mpi_send_wl('T')
      call mpi_recv_wl('T')

      call seacells()

      call mpi_send_ws('T')
      call mpi_recv_ws('T')
   endif

   time_istp8 = time8 + float(istp) * dtsm(mrls)  ! Update precise time
   s1900_sim = s1900_init + time_istp8

call check_nans(18)

enddo

! Call ED model if it is time to do vegetation dynamics

!if (mod(real(time8)+dtlm(1),frq_phenology) < dtlm(1)) call ed_vegetation_dynamics()  

return
end subroutine timestep

!==========================================================================

subroutine check_nans(icall)

use mem_basic,  only: sh_w,rho,thil,sh_v,wc,wmc,press,umc,uc
use mem_micro,  only: sh_c,sh_r
use mem_grid,   only: mza,mwa,lpw,volti,mua
use mem_tend,   only: thilt, umt
use micro_coms, only: level
use misc_coms,  only: io6, iparallel
use mem_ijtabs, only: itab_u, itab_w

implicit none
  
integer :: i,k,iu
integer, intent(in) :: icall

!do iu = 1,mua
!   if (itab_u(iu)%iuglobe == 1206) then
!      write(io6,'(a,i6,5e15.7)') 'cns ',icall,umc(2,iu),uc(2,iu),umt(2,iu)
!   endif
!enddo

return

do i = 2,mwa
   do k = lpw(i),mza
!      if (isnan(sh_w (k,i)) .or.  &
!          isnan(rho  (k,i)) .or.  &
!          isnan(thil (k,i)) .or.  &
!          isnan(thilt(k,i)) .or.  &
!          isnan(press(k,i)) .or.  &
!          isnan(wc   (k,i)) .or.  &
!          isnan(wmc  (k,i)) .or.  &
!          isnan(sh_v (k,i)) .or.  &
!          thil(k,i) < 100.0) then

!         write(io6,*) ''
!         write(io6,*) 'check_nans',k,i,icall
!         write(io6,*) 'sh_w,rho,thil',sh_w(k,i),rho(k,i),thil(k,i)
!         write(io6,*) 'thilt,press',thilt(k,i),press(k,i)
!         write(io6,*) 'wc, wmc, sh_v',wc(k,i),wmc(k,i),sh_v(k,i)
!         stop
!      endif

!      if (level >= 3) then
!         if (isnan(sh_c(k,i)) .or.  &
!             isnan(sh_r(k,i))) then
               
!            write(io6,*) 'nan',k,i,icall
!            write(io6,*) 'sh_c,sh_r',sh_c(k,i),sh_r(k,i)
!            stop
!         endif
!      endif

   enddo
enddo

return
end subroutine check_nans

!==========================================================================

subroutine tend0(rhot)

use mem_ijtabs, only: jtab_w, jtab_u, istp, mrl_begl
use var_tables, only: scalar_tab, num_scalar
use mem_grid,   only: mza, mwa, mua, lcu, lpw
use mem_tend,   only: wmt, umt, thilt
use misc_coms,  only: io6

implicit none

real, intent(in) :: rhot

integer :: n,mrl,j,k,iw,iu

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
do j = 1,jtab_w(14)%jend(mrl); iw = jtab_w(14)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = lpw(iw),mza-1
      wmt(k,iw) = 0.
   enddo
enddo
endif
call rsub('W',14)


! SET U MOMENTUM TENDENCY TO ZERO

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
do j = 1,jtab_u(11)%jend(mrl); iu = jtab_u(11)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)
   do k = lcu(iu),mza-1
      umt(k,iu) = 0.
   enddo
enddo
endif
call rsub('U',11)

return
end subroutine tend0

!==========================================================================

subroutine tnd0(vart)

use mem_ijtabs, only: jtab_w, istp, mrl_begl
use mem_grid,   only: mza, mwa, lpw
use misc_coms,  only: io6

implicit none

real, intent(out) :: vart(mza,mwa)

integer :: j,iw,k,mrl

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
do j = 1,jtab_w(11)%jend(mrl); iw = jtab_w(11)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = lpw(iw)-1,mza
      vart(k,iw) = 0.
   enddo
enddo
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

implicit none

real, intent(in) :: rho_old(mza,mwa)

integer :: n,mrl

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

return
end subroutine predtr

!==========================================================================

subroutine o_update(n,varp,vart,rho_old)

use mem_ijtabs, only: jtab_w, istp, itab_w, mrl_endl
use mem_basic,  only: rho
use misc_coms,  only: io6, dtlm
use mem_grid,   only: mza, mwa, lpw

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
do j = 1,jtab_w(27)%jend(mrl); iw = jtab_w(27)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   dtl = dtlm(itab_w(iw)%mrlw)
   do k = lpw(iw),mza-1
      varp(k,iw) = (varp(k,iw) * rho_old(k,iw) + dtl * vart(k,iw)) / rho(k,iw) 
   enddo
enddo
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

!==========================================================================

subroutine diagn_global()

! subroutine for diagnosis of williamson et al expt 1
     
use mem_basic
use mem_ijtabs
use misc_coms
use consts_coms
use mem_grid
use oplot_coms,  only: op
use mem_addsc,   only: addsc
!------------------------------------------
! Only for ncar test cases:
use ncar_testcases_all, only: ncar_testcase
!------------------------------------------


implicit none

integer, save :: ncall = 0

integer :: j,iw,kk,k
integer :: iu1,iu2,iu3
real :: vxu1,vxu2,vxu3,vyu1,vyu2,vyu3,vzu1,vzu2,vzu3
real :: vx,vy,vz,raxis,u,v,w

real, save, dimension(36000) :: ge1,ge2,ge3,ge4,ge5,ge6,ge7,vctr18,  &
                                q1l1,q1l2,q1li,   &
                                q2l1,q2l2,q2li

character(len=2) :: title

real, save :: aspect = .7
real, save :: scalelab = .012
real, save :: timebeg,timeend,timedif,timeinc
real :: exner, temp
real(kind=8) :: enk, ent, tmass
real(kind=8) :: enk_sum, ent_sum, tmass_sum
real(kind=8) :: q1mass_sum, q2mass_sum, q3mass_sum, q4mass_sum

real(kind=8), save :: enk_sum_init, ent_sum_init, tmass_sum_init
real(kind=8), save :: q1mass_sum_init, q2mass_sum_init
real(kind=8), save :: q3mass_sum_init, q4mass_sum_init

real(kind=8), save :: sum_absq1dif,sum_absq1tr
real(kind=8), save :: sum_q1dif2,sum_q1tr2
real(kind=8), save :: absq1dif_max,absq1tr_max

real(kind=8), save :: sum_absq2dif,sum_absq2tr
real(kind=8), save :: sum_q2dif2,sum_q2tr2
real(kind=8), save :: absq2dif_max,absq2tr_max

real(kind=8) :: q1,q2,q3,q4

real(kind=8), save, allocatable :: q1_tr(:,:),q2_tr(:,:)

ncall = ncall + 1

! On the first call to this subroutine, compute the plotting time increment

if (ncall == 1) then
   if (ncar_testcase == 3) then
      allocate(q1_tr(mza,mwa),q2_tr(mza,mwa))
   endif
   
   timebeg = time8   / 86400.
   timeend = timmax8 / 86400.
   timedif = timeend - timebeg

   if (timedif < .03) then
      timeinc = .001
   elseif (timedif < .06) then
      timeinc = .002
   elseif (timedif < .1) then
      timeinc = .004

   elseif (timedif < .3) then
      timeinc = .01
   elseif (timedif < .6) then
      timeinc = .02
   elseif (timedif < 1.) then
      timeinc = .04

   elseif (timedif < 3.) then
      timeinc = .1
   elseif (timedif < 6.) then
      timeinc = .2
   elseif (timedif < 10.) then
      timeinc = .4

   elseif (timedif < 30.) then
      timeinc = 1.
   elseif (timedif < 60.) then
      timeinc = 2.
   elseif (timedif < 100.) then
      timeinc = 4.

   elseif (timedif < 300.) then
      timeinc = 10.
   elseif (timedif < 600.) then
      timeinc = 20.
   elseif (timedif < 1000.) then
      timeinc = 40.
   endif
      
endif

vctr18(ncall) = time8 / 86400.

! Initialize summation and max/min quantities to zero

enk_sum = 0.
ent_sum = 0.
tmass_sum = 0.
q1mass_sum = 0.
q2mass_sum = 0.
q3mass_sum = 0.
q4mass_sum = 0.

sum_absq1dif = 0.
sum_absq1tr  = 0.

sum_q1dif2   = 0.
sum_q1tr2    = 0.

absq1dif_max = 0.
absq1tr_max  = 0.

sum_absq2dif = 0.
sum_absq2tr  = 0.

sum_q2dif2   = 0.
sum_q2tr2    = 0.

absq2dif_max = 0.
absq2tr_max  = 0.

! Horizontal loop over all IW points

!----------------------------------------------------------------------
do j = 1,jtab_w(8)%jend(1); iw = jtab_w(8)%iw(j)
!---------------------------------------------------------------------

   iu1 = itab_w(iw)%iu1      	
   iu2 = itab_w(iw)%iu2      	
   iu3 = itab_w(iw)%iu3      	

   vxu1 = itab_w(iw)%vxu1; vyu1 = itab_w(iw)%vyu1; vzu1 = itab_w(iw)%vzu1
   vxu2 = itab_w(iw)%vxu2; vyu2 = itab_w(iw)%vyu2; vzu2 = itab_w(iw)%vzu2
   vxu3 = itab_w(iw)%vxu3; vyu3 = itab_w(iw)%vyu3; vzu3 = itab_w(iw)%vzu3

! Vertical loop over all active T levels

   do k = lpw(iw),mza-1

! Zonal and meridional velocity components (U,V) at current (K,IW) point

      vx = vxu1 * uc(k,iu1) + vxu2 * uc(k,iu2) + vxu3 * uc(k,iu3)
      vy = vyu1 * uc(k,iu1) + vyu2 * uc(k,iu2) + vyu3 * uc(k,iu3)
      vz = vzu1 * uc(k,iu1) + vzu2 * uc(k,iu2) + vzu3 * uc(k,iu3)

      raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis
      u = (vy * xew(iw) - vx * yew(iw)) / raxis
      v = vz * raxis / erad  &
        - (vx * xew(iw) + vy * yew(iw)) * zew(iw) / (raxis * erad) 

! Vertical velocity W at current (K,IW) point

      w = .5 * (wc(k-1,iw) + wc(k,iw))

! Air temperature at current (K,IW) point

      exner = (press(k,iw) / p00) ** rocp   ! Defined WITHOUT CP factor
      temp = exner * theta(k,iw)
      
! Kinetic energy per unit mass at current (K,IW) point

      enk = sqrt(u**2 + v**2 + w**2)
      
! Total energy per unit mass at current (K,IW) point

      ent = enk + cv * temp + grav * zt(k)


! Sums of kinetic energy, total energy, mass

      enk_sum   = enk_sum   + enk * rho(k,iw) * volt(k,iw)
      ent_sum   = ent_sum   + ent * rho(k,iw) * volt(k,iw)
      tmass_sum = tmass_sum +       rho(k,iw) * volt(k,iw)

! Sums and max values for q1

      if (ncar_testcase <= 2) then

         q1 = addsc(1)%sclp(k,iw)
         q2 = addsc(2)%sclp(k,iw)
         q3 = addsc(3)%sclp(k,iw)
         q4 = addsc(4)%sclp(k,iw)

         q1mass_sum = q1mass_sum + q1 * rho(k,iw) * volt(k,iw)
         q2mass_sum = q2mass_sum + q2 * rho(k,iw) * volt(k,iw)
         q3mass_sum = q3mass_sum + q3 * rho(k,iw) * volt(k,iw)
         q4mass_sum = q4mass_sum + q4 * rho(k,iw) * volt(k,iw)

      endif

      if (ncar_testcase == 3) then

         q1 = addsc(1)%sclp(k,iw)
         q2 = addsc(2)%sclp(k,iw)

         q1mass_sum = q1mass_sum + q1 * rho(k,iw) * volt(k,iw)
         q2mass_sum = q2mass_sum + q2 * rho(k,iw) * volt(k,iw)

! On the first call to this subroutine, save initial values

         if (ncall == 1) then
            q1_tr(k,iw) = q1
            q2_tr(k,iw) = q2
         endif

         sum_absq1dif = sum_absq1dif + abs(q1 - q1_tr(k,iw)) * rho(k,iw) * volt(k,iw)
         sum_absq1tr  = sum_absq1tr  + abs(q1_tr(k,iw))      * volt(k,iw)

         sum_q1dif2   = sum_q1dif2   + (q1 - q1_tr(k,iw))**2 * volt(k,iw)
         sum_q1tr2    = sum_q1tr2    + q1_tr(k,iw)**2        * volt(k,iw)

         absq1dif_max = max(absq1dif_max, abs(q1 - q1_tr(k,iw)))
         absq1tr_max  = max(absq1tr_max, q1_tr(k,iw))

         sum_absq2dif = sum_absq2dif + abs(q2 - q2_tr(k,iw)) * rho(k,iw) * volt(k,iw)
         sum_absq2tr  = sum_absq2tr  + abs(q2_tr(k,iw))      * volt(k,iw)

         sum_q2dif2   = sum_q2dif2   + (q2 - q2_tr(k,iw))**2 * volt(k,iw)
         sum_q2tr2    = sum_q2tr2    + q2_tr(k,iw)**2        * volt(k,iw)

         absq2dif_max = max(absq2dif_max, abs(q2 - q2_tr(k,iw)))
         absq2tr_max  = max(absq2tr_max, q2_tr(k,iw))

      endif

   enddo
   
enddo

! first time in this routine, save initial values

if (ncall == 1) then

   enk_sum_init   = enk_sum
   ent_sum_init   = ent_sum
   tmass_sum_init = tmass_sum

   if (ncar_testcase <= 3) then
      q1mass_sum_init = q1mass_sum
      q2mass_sum_init = q2mass_sum
   endif
   
   if (ncar_testcase <= 2) then
      q3mass_sum_init = q3mass_sum
      q4mass_sum_init = q4mass_sum
   endif
   
endif

! Time series of normalized global integrals

if (abs(enk_sum_init) < 1.e-14_8) then
   ge1(ncall) = enk_sum
else
   ge1(ncall) = (enk_sum - enk_sum_init) / enk_sum_init
endif
ge2(ncall) = (ent_sum - ent_sum_init) / ent_sum_init
ge3(ncall) = (tmass_sum - tmass_sum_init) / tmass_sum_init

if (ncar_testcase == 3) then

   ge4(ncall) = (q1mass_sum - q1mass_sum_init) / q1mass_sum_init
   ge5(ncall) = (q2mass_sum - q2mass_sum_init) / q2mass_sum_init

! First 3 normalized global errors (NCAR cases Eqs. 54-56)

   q1l1(ncall) = sum_absq1dif / sum_absq2tr
   q1l2(ncall) = sqrt(sum_q1dif2) / sqrt(sum_q2tr2)
   q1li(ncall) = absq1dif_max / absq2tr_max

   q2l1(ncall) = sum_absq2dif / sum_absq2tr
   q2l2(ncall) = sqrt(sum_q2dif2) / sqrt(sum_q2tr2)
   q2li(ncall) = absq2dif_max / absq2tr_max

   write(6,'(a,6e13.4)') 'q12norms ',q1l1(ncall),q1l2(ncall),q1li(ncall),  &
                                     q2l1(ncall),q2l2(ncall),q2li(ncall)
endif

if (ncar_testcase <= 2) then

   ge4(ncall) = (q1mass_sum - q1mass_sum_init) / q1mass_sum_init
   ge5(ncall) = (q2mass_sum - q2mass_sum_init) / q2mass_sum_init
   ge6(ncall) = (q3mass_sum - q3mass_sum_init) / q3mass_sum_init
   ge7(ncall) = (q4mass_sum - q4mass_sum_init) / q4mass_sum_init

endif

!print*, 'ge3-mass ',ge3(ncall)

if (time8 + 1.5 * dtlong > timmax8) then

! Reopen the current graphics output workstation if it is closed

   call o_reopnwk()

! Plot time series

!----------------------------------------------------------------------
   call plotback()
   call oplot_xy2('0','N',aspect,scalelab   &
                 ,ncall,  vctr18,ge1            &
                 ,'time(days)','kinetic energy' &
                 ,timebeg,timeend,timeinc,5  ,-1.e-1,1.e-1,0.1e-1,5  )
   call o_frame()
!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2('0','N',aspect,scalelab   &
                 ,ncall,  vctr18,ge2            &
                 ,'time(days)','total energy' &
                 ,timebeg,timeend,timeinc,5  ,-1.e-3,1.e-3,0.1e-3,5  )
   call o_frame()
!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2('0','N',aspect,scalelab        &
                 ,ncall,  vctr18,ge3                 &
                 ,'time(days)','total mass deviation' &
                 ,timebeg,timeend,timeinc,5  ,-1.e-12,1.e-12,0.1e-12,10  )
   call o_frame()
!-------------------------------------------------------------------

   if (ncar_testcase == 3) then

!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2('0','N',aspect,scalelab        &
                 ,ncall,  vctr18,ge4                 &
                 ,'time(days)','q5 mass deviation' &
                 ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
   call o_frame()
!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2('0','N',aspect,scalelab        &
                 ,ncall,  vctr18,ge5                 &
                 ,'time(days)','q6 mass deviation' &
                 ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
   call o_frame()
!-------------------------------------------------------------------

! Tracer mass norms

!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2log10('0','N',aspect,scalelab      &
                 ,ncall,  vctr18,q1l1               &
                 ,'time(days)','q5 mass l1 norm'    &
                 ,timebeg,timeend,timeinc,5  ,-4,1  )
   call o_frame()
!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2log10('0','N',aspect,scalelab      &
                 ,ncall,  vctr18,q1l2               &
                 ,'time(days)','q5 mass l2 norm'    &
                 ,timebeg,timeend,timeinc,5  ,-4,1  )
   call o_frame()
!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2log10('0','N',aspect,scalelab      &
                 ,ncall,  vctr18,q1li               &
                 ,'time(days)','q5 mass l_inf norm' &
                 ,timebeg,timeend,timeinc,5  ,-4,1  )
   call o_frame()
!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2log10('0','N',aspect,scalelab      &
                 ,ncall,  vctr18,q2l1               &
                 ,'time(days)','q6 mass l1 norm'    &
                 ,timebeg,timeend,timeinc,5  ,-4,1  )
   call o_frame()
!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2log10('0','N',aspect,scalelab      &
                 ,ncall,  vctr18,q2l2               &
                 ,'time(days)','q6 mass l2 norm'    &
                 ,timebeg,timeend,timeinc,5  ,-4,1  )
   call o_frame()
!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2log10('0','N',aspect,scalelab      &
                 ,ncall,  vctr18,q2li               &
                 ,'time(days)','q6 mass linf norm'  &
                 ,timebeg,timeend,timeinc,5  ,-4,1  )
   call o_frame()
!-------------------------------------------------------------------

   endif

   if (ncar_testcase <= 2) then

!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2('0','N',aspect,scalelab        &
                 ,ncall,  vctr18,ge4                 &
                 ,'time(days)','q1 mass deviation' &
                 ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
   call o_frame()
!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2('0','N',aspect,scalelab        &
                 ,ncall,  vctr18,ge5                 &
                 ,'time(days)','q2 mass deviation' &
                 ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
   call o_frame()
!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2('0','N',aspect,scalelab        &
                 ,ncall,  vctr18,ge6                 &
                 ,'time(days)','q3 mass deviation' &
                 ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
   call o_frame()
!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2('0','N',aspect,scalelab        &
                 ,ncall,  vctr18,ge7                 &
                 ,'time(days)','q4 mass deviation' &
                 ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
   call o_frame()
!-------------------------------------------------------------------

   endif

! Close the current workstation if output
! is to a NCAR graphics meta file. This allows viewing the complete
! meta file (including the last frame) during a run and in case the
! simulation crashes.

   call o_clswk()

endif

return
end subroutine diagn_global
