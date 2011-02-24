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
subroutine inithh()
! +--------------------------------------------------------+
! _  Initialization of model grid, topography, and fields  _
! _    for horizontally homogeneous fields.                _
! +--------------------------------------------------------+

implicit none

call arrsnd()  !  Arrange the input sounding.
call refs1d()  !  Compute the 1-D reference state variables 
call fldshhi()  !  Compute the 3-D model fields.

return
end subroutine inithh

!===============================================================================

subroutine arrsnd()

use consts_coms, only: pio180, eps_virt, grav, rdry, p00k, rocp, cp, cpor,  &
                       p00i, p00, eps_vap
use misc_coms,   only: io6, nsndg, ipsflg, itsflg, irtsflg, iusflg,  &
                       ps, ts, rts, us, vs, hs, thds, p_sfc

implicit none

integer :: ksndg
real :: toffset,dir,spd,zold2,zold1,tavg,rtss,wt,ssh_sndg,vapor_press
real, external :: eslf
real :: tvirt(nsndg)  ! automatic array

! +----------------------------------------------------------------------+
! RTS and RT01D now defined as specific humidity instead of mixing ratio
! +----------------------------------------------------------------------+

!     Arrange the input sounding

if (itsflg == 0) then
   toffset = 273.15
elseif (itsflg == 1) then
   toffset = 0.
endif

do ksndg = 1,nsndg

   if (iusflg /= 0) then
      dir = us(ksndg)
      spd = vs(ksndg)
      us(ksndg) = -spd * sin(pio180*dir)
      vs(ksndg) = -spd * cos(pio180*dir)
   endif

   if (ipsflg == 0) then
!     PS array is pressure in millibars with HS(1) = height of lowest sounding level
      ps(ksndg) = ps(ksndg) * 100.
      
   elseif (ipsflg == 1) then
!     PS array is height in meters with P_SFC = surface pressure

!     If sounding moisture is expressed as a mixing ratio (IRTSFLG=2),
!     take advantage of knowing virtual temperature effect when
!     integrating hydrostatically to get sounding pressure.
!
      if (irtsflg == 2) then
         tvirt(ksndg) = (1. + eps_virt * rts(ksndg) * 1.e-3)
      else
         tvirt(ksndg) = 1.
      endif

      if (ksndg == 1) then
         zold2 = ps(ksndg)
         ps(ksndg) = p_sfc * 100.
      else
         zold1 = zold2
         zold2 = ps(ksndg)
         if (itsflg == 0 .or. itsflg == 1) then
            tavg = .5 * ( (ts(ksndg) + toffset) * tvirt(ksndg)  &
                 +        ts(ksndg-1) * tvirt(ksndg-1) )
            ps(ksndg) = ps(ksndg-1) * exp(-grav * (zold2 - zold1) / (rdry * tavg))
         elseif (itsflg == 2) then
            tavg = (ts(ksndg) * tvirt(ksndg)  &
                + ts(ksndg-1) * tvirt(ksndg-1) * p00k / ps(ksndg-1)**rocp) * .5
            ps(ksndg) = (ps(ksndg-1)**rocp - grav * (zold2 - zold1) * p00k  &
               / (cp * tavg))**cpor
         endif
      endif
   endif

   if (itsflg == 0) then
!     Temperature in degrees celsius
      ts(ksndg) = ts(ksndg) + 273.15
   elseif (itsflg == 1) then
!     Temperature in degrees kelvin
   elseif (itsflg == 2) then
!     Temperature is potential temperature in kelvin
      ts(ksndg) = (ps(ksndg) * p00i)**rocp * ts(ksndg)
   endif

   if (irtsflg == 0) then

! RTS is given as dew point in degrees celsius
! Compute ambient vapor pressure

      vapor_press = eslf(rts(ksndg))

! Do not allow vapor pressure to exceed ambient pressure

      vapor_press = min(ps(ksndg),vapor_press)

! Compute specific humidity (pcol_rt) from vapor pressure and ambient pressure

      rts(ksndg) = eps_vap * vapor_press  &
                 / (ps(ksndg) + vapor_press * (eps_vap - 1.))

   elseif (irtsflg == 1) then

! RTS is given as dew point in degrees kelvin
! Compute ambient vapor pressure

      vapor_press = eslf(rts(ksndg)-273.15)

! Do not allow vapor pressure to exceed ambient pressure

      vapor_press = min(ps(ksndg),vapor_press)

! Compute specific humidity (pcol_rt) from vapor pressure and ambient pressure

      rts(ksndg) = eps_vap * vapor_press  &
                 / (ps(ksndg) + vapor_press * (eps_vap - 1.))

   elseif (irtsflg == 2) then

! RTS is given as specific humidity in g/kg

      rts(ksndg) = rts(ksndg) * 1.e-3

   elseif (irtsflg == 3) then

! RTS is given as relative humidity in percent
! Compute ambient vapor pressure based on relative humidity (pcol_r)
! and saturation vapor pressure

      vapor_press = .01 * rts(ksndg) * eslf(ts(ksndg)-273.15)

! Do not allow vapor pressure to exceed ambient pressure

      vapor_press = min(ps(ksndg),vapor_press)

! Compute specific humidity from vapor pressure and ambient pressure

      rts(ksndg) = eps_vap * vapor_press  &
                 / (ps(ksndg) + vapor_press * (eps_vap - 1.))

   elseif (irtsflg == 4) then

! RTS is given as dew point depression in kelvin
! Compute ambient vapor pressure

      vapor_press = eslf(ts(ksndg)-rts(ksndg)-273.15)

! Do not allow vapor pressure to exceed ambient pressure

      vapor_press = min(ps(ksndg),vapor_press)

! Compute specific humidity (pcol_rt) from vapor pressure and ambient pressure

      rts(ksndg) = eps_vap * vapor_press  &
                 / (ps(ksndg) + vapor_press * (eps_vap - 1.))

   endif
enddo

!     compute height levels of input sounding.

do ksndg = 2,nsndg
   hs(ksndg) = hs(ksndg-1) - rdry * .5  &
      * (ts(ksndg) * (1. + eps_virt * rts(ksndg))   &
      + ts(ksndg-1) * (1. + eps_virt * rts(ksndg-1)))  &
      * (log(ps(ksndg)) - log(ps(ksndg-1))) / grav

!   write(io6,*) 'hs_1 ',ksndg,hs(ksndg-1),hs(ksndg),ts(ksndg),ps(ksndg)  &
!                        ,ts(ksndg-1),ps(ksndg-1),rts(ksndg),rts(ksndg-1)
enddo

do ksndg = 1,nsndg
   thds(ksndg) = ts(ksndg) * (p00 / ps(ksndg))**rocp
enddo

return
end subroutine arrsnd

!===============================================================================

subroutine refs1d()
! +---------------------------------------------------------------------
! \   This routine computes the reference state sounding on the model
! \     levels from input sounding defined on pressure levels.
! +---------------------------------------------------------------------

use misc_coms,   only: io6, nsndg, hs, thds, us, vs, rts, ps,  &
                       pr01d, dn01d, rt01d, th01d, u01d, v01d
use consts_coms, only: cvocp, p00k, rdry, eps_virt, gravo2
use mem_grid,    only: mza, zm, zt, dzt
use micro_coms,  only: level

implicit none

integer :: k,kk,iter
real :: dens1

write(io6,*) 'Beginning refs1d '

if (zt(mza) > hs(nsndg)) then
   write(io6,*) ' !!! Input sounding is not high enough !!!'
   write(io6,*) ' !!! Sounding top (m): ',hs(nsndg)
   write(io6,*) ' !!! Model top (m):    ',zt(mza)
   stop 'refs1d'
endif

call htint(nsndg,thds,hs,mza,th01d,zt)
call htint(nsndg,  us,hs,mza, u01d,zt)
call htint(nsndg,  vs,hs,mza, v01d,zt)
if (level >= 1) then
   call htint(nsndg,rts,hs,mza,rt01d,zt)
else
   rt01d(1:mza) = 0.
endif

! New interpolation of pressure prior to iterative hydrostatic integration

call htint(nsndg,ps,hs,mza,pr01d,zt)
dens1 = ps(1) ** cvocp * p00k / (rdry * thds(1) * (1. + eps_virt * rts(1)))

! Use iterative method for hydrostatic integration (assume no condensate)

do iter = 1,100
   dn01d(1) = pr01d(1) ** cvocp * p00k  &
            / (rdry * th01d(1) * (1. + eps_virt * rt01d(1)))
   pr01d(1) = ps(1) - gravo2 * (dens1 + dn01d(1)) * (zt(1) - hs(1)) 

   do k = 2,mza
      dn01d(k) = pr01d(k) ** cvocp * p00k  &
               / (rdry * th01d(k) * (1. + eps_virt * rt01d(k)))
! Impose minimum value of 1 Pa to avoid overshoot to negative values 
! during iteration
      pr01d(k) = max(1.,  &
         pr01d(k-1) - gravo2 * (dn01d(k-1) * dzt(k-1) + dn01d(k) * dzt(k)))
   enddo
enddo

rt01d(1) = rt01d(2)
th01d(1) = th01d(2)

! Print out initial state column 

write(io6,*) ' '
write(io6,*) '========================================================================='
write(io6,*) '                    OLAM INITIAL STATE COLUMN (hhi)'
write(io6,*) '========================================================================='
write(io6,*) '   zm(m)    k     zt(m)   pr01d(Pa)  dn01d(kg/m3)  th01d(K)  rt01d(g/kg)'
write(io6,*) '========================================================================='
write(io6,*) ' '

do k = mza-1,2,-1
   write(io6, '(f10.2,1x,9(''-------''))') zm(k)
   write(io6, '(10x,i5,f10.2,f11.2,f11.4,f12.2,f11.4)')  &
       k,zt(k),pr01d(k),dn01d(k),th01d(k),rt01d(k)*1.e3
enddo
   
write(io6, '(f10.2,1x,9(''-------''))') zm(1)
write(io6,*) ' '

return
end subroutine refs1d

!===============================================================================

subroutine fldshhi()
   
use mem_basic,   only: theta, thil, press, rho, wc, wmc, uc, ump, umc,  &
                       sh_w, sh_v
use mem_micro,   only: sh_c
use micro_coms,  only: level
use mem_ijtabs,  only: jtab_w, jtab_u, itab_w, itab_u, rotation_angle
use misc_coms,   only: io6, mdomain, th01d, pr01d, dn01d, rt01d, u01d, v01d,  &
                       iparallel, naddsc
use consts_coms, only: cvocp, p00k, rdry, rvap, p00, rocp, alvlocp,  &
                       gravo2, erad, grav, &
                       pi1, pio180, pi2, omega, piu180, pio2
use mem_grid,    only: mza, dzt, mua, lpu, lpw, unx, uny, unz, xeu, yeu, zeu,  &
                       aru, volt, volui,  &
                       glatw, glonw, xem, yem, zem, zew, zt   
                       
use mem_addsc,   only: addsc
                       
!nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn
! Module input for NCAR test cases

use leaf_coms, only: nslcon,nvgcon
use sea_coms,  only: seatmp

! In the NCAR test cases, the above 2 variables (nvgcon and seatmp) are not
! required for their usual purpose in OLAM, so they are used here as a means 
! for specifying ncar_testcase number and rotation angle of earth's axis of 
! rotation.  Their values are specified in the namelist file OLAMIN.

use ncar_testcases_all                  

! end ncar
!nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn

implicit none   

integer :: j,iw,k,iu,iv,iter,iw1,iw2,iup,ivp,iuv,im1,im2,iwp
real :: rcloud,dummy,temp,rvls,exner  &
   ,uv01dx,uv01dy,uv01dz,uv01dr,raxis,rhovs

real, external :: rhovsl

!nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn
! NCAR test cases

real(r8) :: uuu,vvv,phyd
real(r8) :: eta, tempn, slon, slat, phis, eta_c, psfc, q5, q6, pressu

ncar_testcase = nslcon
ncar_choice   = nvgcon
rotation_angle = seatmp

! end ncar
!nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(8)%jend(1); iw = jtab_w(8)%iw(j)
!---------------------------------------------------------------------
call qsub('W',iw)

! Fill arrays prior to iteration

   do k = 1,mza
      theta(k,iw) = th01d(k)
      thil(k,iw)  = theta(k,iw)
      press(k,iw) = pr01d(k)
      rho(k,iw)   = dn01d(k)
      wc(k,iw)    = 0.
      wmc(k,iw)   = 0.
      
      if (level == 0) then
         sh_w(k,iw) = 0.
         sh_v(k,iw) = 0.
      elseif (level == 1) then
         sh_w(k,iw) = rt01d(k) 
         sh_v(k,iw) = rt01d(k)
      else
         sh_w(k,iw) = rt01d(k)
         sh_v(k,iw) = rt01d(k)
      endif
   enddo

   do iter = 1,100
   
! Diagnose density in all levels

      do k = 1,mza
      
         if (level == 0) then
            rho(k,iw) = press(k,iw) ** cvocp * p00k / (rdry * theta(k,iw))
         elseif (level == 1) then
            rho(k,iw) = press(k,iw) ** cvocp * p00k  &
               / (theta(k,iw) * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))
         else
            exner = (press(k,iw) / p00) ** rocp
            temp = exner * theta(k,iw)
            rhovs = rhovsl(temp-273.15)
            sh_c(k,iw) = max(0.,sh_w(k,iw)-rhovs/real(rho(k,iw)))
            sh_v(k,iw) = sh_w(k,iw) - sh_c(k,iw)
! As for iteration in subroutine satadjst, use (0.3,0.7) weighting to damp iteration
!            rho(k,iw) = .5 * rho(k,iw)  &
!                      + .5 * press(k,iw) ** cvocp * p00k  &
            rho(k,iw) = press(k,iw) ** cvocp * p00k  &
               / (theta(k,iw) * (1. - sh_c(k,iw))  &
               * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))
            thil(k,iw) = theta(k,iw)  &
               / (1. + alvlocp * sh_c(k,iw) / max(temp,253.))
         endif
         
      enddo
         
      k = lpw(iw)

      do while (k <= mza)

!nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn
! NCAR test cases: Use temperature field for 1 & 2

         if (ncar_testcase == 1 .or.  &
             ncar_testcase == 2) then

            eta = 1.d-5 * press(k,iw)

            tempn = temperature (real(pio180*glonw(iw),r8),  &
                                 real(pio180*glatw(iw),r8),  &
                                 eta,rotation_angle)

            if (k == lpw(iw)) then
               phis = surface_geopotential(dble(pio180*glonw(iw)),  &
                                           dble(pio180*glatw(iw)),  &
                                           rotation_angle)

               psfc = 1.d5
            endif
            
         elseif (ncar_testcase == 3) then

            call advection ('56',                            &
                            real(pio180*glonw(iw),r8),       &
                            real(pio180*glatw(iw),r8),       &
                            real(zt(k),r8), rotation_angle,  &
                            uuu, vvv, tempn,                 &
                            phis,psfc,                       &
                            q5,q6                            )

            addsc(1)%sclp(k,iw) = q5
            addsc(2)%sclp(k,iw) = q6
            
         elseif (ncar_testcase == 4) then

            call Rossby_Haurwitz (real(pio180*glonw(iw),r8),  &
                                  real(pio180*glatw(iw),r8),  &
                                  press(k,iw),                &
                                  uuu, vvv, tempn,            &
                                  phis,psfc                   )

         elseif (ncar_testcase == 5) then

            call mountain_Rossby (real(pio180*glonw(iw),r8),  &
                                  real(pio180*glatw(iw),r8),  &
                                  press(k,iw),                &
                                  uuu, vvv, tempn,            &
                                  phis,psfc                   )

         elseif (ncar_testcase == 6) then

            call gravity_wave (ncar_choice,                &
                               real(pio180*glonw(iw),r8),  &
                               real(pio180*glatw(iw),r8),  &
                               real(zt(k),r8),             &
                               uuu, vvv, tempn,            &
                               phis,psfc                   )

         endif

         exner = (press(k,iw) / p00) ** rocp
         theta(k,iw) = tempn / exner

         thil(k,iw) = theta(k,iw)

! Hydrostatic integration

         if (k == lpw(iw)) then

! If k = lpw(iw), hydrostatic integration from (phis, psfc) surface using 
! density in k = lpw(iw) level

            phyd = psfc - rho(k,iw) * (grav * zt(k) - phis)
               
            press(k,iw) = phyd

         elseif (k < lpw(iw)) then

! If k < lpw(iw), hydrostatic integration from higher level to get current level

            phyd = press(k+1,iw)  &
                 + gravo2 * (rho(k+1,iw) * dzt(k+1) + rho(k,iw) * dzt(k))

            press(k,iw) = phyd

         else

! If k > lpw(iw), hydrostatic integration from lower level to get current level

            phyd = press(k-1,iw)  &
                 - gravo2 * (rho(k-1,iw) * dzt(k-1) + rho(k,iw) * dzt(k))

! Impose minimum value of 1 Pa to avoid overshoot to negative values 
! during iteration.

            phyd = max(1.0_r8,phyd)

! Use weighting to damp oscillations

            press(k,iw) = .05_r8 * press(k,iw)  &
                        + .95_r8 * phyd

         endif
             
   
! end ncar
!nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn

         if (k <= lpw(iw)) k = k - 1
         if (k > lpw(iw))  k = k + 1
         if (k < 1)        k = lpw(iw) + 1
         
      enddo
   enddo

!nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn
! NCAR tracer variables q1-q4

   if (ncar_testcase == 1 .or.  &
       ncar_testcase == 2) then

      do k = 1,mza

         eta = 1.d-5 * press(k,iw)
      
         if (naddsc >= 1) then
            eta_c = 0.6d0

            addsc(1)%sclp(k,iw) = tracer_q1_q2(dble(pio180*glonw(iw)),  &
                                               dble(pio180*glatw(iw)),  &
                                               eta,rotation_angle,eta_c ) 
         endif
                                         
         if (naddsc >= 2) then
            eta_c = 1.d0

            addsc(2)%sclp(k,iw) = tracer_q1_q2(dble(pio180*glonw(iw)),  &
                                               dble(pio180*glatw(iw)),  &
                                               eta,rotation_angle,eta_c ) 
         endif
                                         
         if (naddsc >= 3) then
            addsc(3)%sclp(k,iw) = tracer_q3(dble(pio180*glonw(iw)),  &
                                            dble(pio180*glatw(iw)),  &
                                            eta,rotation_angle       ) 
         endif
                                         
         if (naddsc >= 4) then
            addsc(4)%sclp(k,iw) = 1.
         endif
                                         
      enddo
      
   endif

! end ncar
!nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn


   if (ncar_testcase == 3) then
      do k = 1,mza
         rho(k,iw) = 1.
      enddo
   endif

enddo

! Set reference state to 3d fields at IW = 2.

pr01d(1:mza) = press(1:mza,2)
dn01d(1:mza) = rho(1:mza,2)
th01d(1:mza) = theta(1:mza,2)

call rsub('Wa',8)

! Copy WMC, etc to open lateral boundary points

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(9)%jend(1); iw = jtab_w(9)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
call qsub('W',iw)

   do k = 1,mza
      wmc(k,iw)   = wmc(k,iwp) 
      wc(k,iw)    = wc(k,iwp) 
      rho(k,iw)   = rho(k,iwp) 
      press(k,iw) = press(k,iwp)
      thil(k,iw)  = thil(k,iwp)
      theta(k,iw) = theta(k,iwp)
      sh_w(k,iw)  = sh_w(k,iwp)
      sh_v(k,iw)  = sh_v(k,iwp)
   enddo

   if (level > 1) then
      do k = 1,mza
         sh_c(k,iw) = sh_c(k,iwp)
      enddo
   endif

enddo
call rsub('Wa',9)

if (iparallel == 1) then
   call mpi_send_w('I',wc)  ! Send W group; wc used as dummy array here
   call mpi_recv_w('I',wc)  ! Recv W group; wc used as dummy array here
endif

! UMC, UC

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_u(8)%jend(1); iu = jtab_u(8)%iu(j)
   iw1 = itab_u(iu)%iw1; iw2 = itab_u(iu)%iw2
!----------------------------------------------------------------------
call qsub('U',iu)

   if (iw1 == 1) iw1 = iw2
   if (iw2 == 1) iw2 = iw1

! If sounding winds are to be interpreted as eastward (U) and 
! northward (V) components, rotate winds from geographic to
! polar stereographic orientation

! U point coordinates and normal vector components

   do k = 1,mza

      if (mdomain <= 1) then  ! Model uses "earth" coordinates

! Standard section

         raxis = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)  ! dist from earth axis
         uv01dr = -v01d(k) * zeu(iu) / erad  ! radially outward from axis

         uv01dx = (-u01d(k) * yeu(iu) + uv01dr * xeu(iu)) / raxis 
         uv01dy = ( u01d(k) * xeu(iu) + uv01dr * yeu(iu)) / raxis 
         uv01dz =   v01d(k) * raxis / erad 

         uc(k,iu) = uv01dx * unx(iu) + uv01dy * uny(iu) + uv01dz * unz(iu)

! End standard section

!nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn
! NCAR test cases: Velocity field for cases 1-6

         slon = pio180 * .5 * (glonw(iw1) + glonw(iw2))
         slat = pio180 * .5 * (glatw(iw1) + glatw(iw2))
         pressu = .5 * (press(k,iw1) + press(k,iw2))

! Correction for points near +/- 180 degrees longitude

         if (abs(glonw(iw1) - glonw(iw2)) > 180.) slon = slon + pi1

         if     (ncar_testcase == 1) then

            eta = 1.d-5 * pressu
            uuu = u_wind(slon,slat,eta,.false.,rotation_angle)
            vvv = v_wind(slon,slat,eta,.false.,rotation_angle)                  

         elseif (ncar_testcase == 2) then

            eta = 1.d-5 * pressu
            uuu = u_wind(slon,slat,eta,.true. ,rotation_angle)
            vvv = v_wind(slon,slat,eta,.true. ,rotation_angle)                  

         elseif (ncar_testcase == 3) then

            call advection ('00',                             &
                            slon,                            &
                            slat,                            &
                            real(zt(k),r8), rotation_angle,  &
                            uuu, vvv, tempn,                 &
                            phis,psfc,                       &
                            q5,q6                            )

         elseif (ncar_testcase == 4) then

            call Rossby_Haurwitz (slon,             &
                                  slat,             &
                                  pressu,           &
                                  uuu, vvv, tempn,  &
                                  phis,psfc         )

         elseif (ncar_testcase == 5) then

            call mountain_Rossby (slon,             &
                                  slat,             &
                                  pressu,           &
                                  uuu, vvv, tempn,  &
                                  phis,psfc         )

         elseif (ncar_testcase == 6) then

            call gravity_wave (ncar_choice,      &
                               slon,             &
                               slat,             &
                               real(zt(k),r8),   &
                               uuu, vvv, tempn,  &
                               phis,psfc         )

         endif

         raxis = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)  ! dist from earth axis
         uv01dr = -vvv * zeu(iu) / erad  ! radially outward from axis

         uv01dx = (-uuu * yeu(iu) + uv01dr * xeu(iu)) / raxis 
         uv01dy = ( uuu * xeu(iu) + uv01dr * yeu(iu)) / raxis 
         uv01dz =   vvv * raxis / erad 

         uc(k,iu) = uv01dx * unx(iu) + uv01dy * uny(iu) + uv01dz * unz(iu)

! end ncar
!nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn

      else
         uc(k,iu) = u01d(k) * unx(iu) + v01d(k) * uny(iu)
      endif

      umc(k,iu) = uc(k,iu) * volui(k,iu)  &
         * (volt(k,iw2) * rho(k,iw1) + volt(k,iw1) * rho(k,iw2))

   enddo      

enddo
call rsub('Ua',8)

! Copy UMC, UC to inner and outer lateral boundary points

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_u(9)%jend(1); iu = jtab_u(9)%iu(j)
   iup = itab_u(iu)%iup
!----------------------------------------------------------------------
call qsub('U',iu)
   do k = 1,mza-1
      umc(k,iu) = umc(k,iup)
      uc(k,iu)  = uc(k,iup)
   enddo
enddo
call rsub('Ua',9)!RRR

! Set UMC, UC to zero wherever ARU = 0.

call psub()!PPPPPPPP
!----------------------------------------------------------------------
do iu = 2,mua
!----------------------------------------------------------------------
call qsub('U',iu) !QQQQQ
   do k = 1,mza-1
      if (aru(k,iu) < 1.e-9) then
         umc(k,iu) = 0.
         uc(k,iu)  = 0.
      endif
   enddo
enddo
call rsub('Ua',0)

if (iparallel == 1) then
   call mpi_send_u('I')  ! Send U group
   call mpi_recv_u('I')  ! Recv U group
endif

! Set UMP to UMC

ump(:,:) = umc(:,:)

return
end subroutine fldshhi
