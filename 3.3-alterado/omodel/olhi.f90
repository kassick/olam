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
subroutine fldslhi()

use mem_basic,   only: theta, thil, rho, press, sh_w, sh_v,  &
                       wc, wmc, uc, umc, ump
use mem_micro,   only: sh_c
use micro_coms,  only: level
use mem_ijtabs,  only: jtab_w, jtab_u, itab_u
use consts_coms, only: p00, rocp, cvocp, p00k, rdry, rvap, alvlocp, gravo2
use mem_grid,    only: zt, dzt, xeu, zeu, yeu, unx, uny, mua, zm,  &
                       mza, mwa, glatw, aru
use mem_zonavg,  only: zonz_vect, zonu_vect, zont_vect, zonr_vect,  &
                       zonp_vect, zonz, zonu, zont, zonr
use misc_coms,   only: io6, iparallel
use rastro_evts
    

implicit none   

integer :: j,iw,k,iu,iv,iter,iw1,iw2,iup,ivp,iuv,lf,lv,llat  &
   ,im,iplev,ilat,im1,im2,ilatn,ilats
real :: rcloud,qlatu,qlonu,dummy,qlatv,qlonv,temp,rvls,exner  &
   ,qlatuv,qlonuv,ugx,ugy,raxis,pnorth,psouth,ug,ug1,ug2,wt1,pkhyd  &
   ,alat,dyo2g,fcorn,fcors,pilo,pihi,thetavlo,thetavhi,rlat,wt2,cpo2g,rhovs
character(len=1) :: line

real, dimension(mza) :: vctr1,vctr2,vctr3  ! automatic arrays
real, dimension(mza,mwa) :: uzonal  ! automatic array

real, external :: rhovsl

#ifdef OLAM_RASTRO
character*1 :: rst_buf = '_'
call rst_event_s_f(OLAM_FLDSLHI_IN,rst_buf)
#endif

call zonavg_init()  ! Read in 'ZONAVG_CLIMATE' data and fill zonavg arrarys

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(8)%jend(1); iw = jtab_w(8)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

! Linearly interpolate zonavg arrays by latitude to current IW column 
! and K level

   rlat = .4 * (glatw(iw) + 93.75)
   ilat = int(rlat)
   wt2 = rlat - float(ilat)
   
   zonz_vect(1:22) = (1. - wt2) * zonz(ilat,1:22) + wt2 * zonz(ilat+1,1:22)
   zont_vect(1:22) = (1. - wt2) * zont(ilat,1:22) + wt2 * zont(ilat+1,1:22)
   zonu_vect(1:22) = (1. - wt2) * zonu(ilat,1:22) + wt2 * zonu(ilat+1,1:22)
   zonr_vect(1:22) = (1. - wt2) * zonr(ilat,1:22) + wt2 * zonr(ilat+1,1:22)
   
! Interpolate zonavg vector arrays in height to model levels

! Fill pressure, theta, air density, and vapor density arrays from zonavg
! vector arrays prior to iteration

   call hintrp_ee(22,zonp_vect,zonz_vect,mza,vctr1,zt) ! pressure
   call hintrp_ee(22,zont_vect,zonz_vect,mza,vctr2,zt) ! temp
   call hintrp_ee(22,zonr_vect,zonz_vect,mza,vctr3,zt) ! specific humidity
   call hintrp_ee(22,zonu_vect,zonz_vect,mza,uzonal(1,iw),zt) ! uzonal

   press(1:mza,iw) = vctr1(1:mza)
   theta(1:mza,iw) = vctr2(1:mza) * (p00 / vctr1(1:mza)) ** rocp  ! temp to theta
   thil(1:mza,iw) = theta(1:mza,iw)
   rho(1:mza,iw)  = press(1:mza,iw) ** cvocp * p00k / (rdry * theta(1:mza,iw))
   sh_v(1:mza,iw) = vctr3(1:mza)  ! spec hum
   sh_w(1:mza,iw) = sh_v(1:mza,iw)
   wc(1:mza,iw)   = 0.
   wmc(1:mza,iw)  = 0.

! Iterative hydrostatic integration

   do iter = 1,100
   
      do k = 1,mza

! Try this: hold Mclatchy temp (vctr2) constant during iterations
      
         theta(k,iw) = vctr2(k) * (p00 / press(k,iw)) ** rocp
         thil(k,iw) = theta(k,iw)

         if (level == 0) then
            rho(k,iw) = press(k,iw) ** cvocp * p00k / (rdry * theta(k,iw))
         elseif (level == 1) then
            rho(k,iw) = press(k,iw) ** cvocp * p00k  &
               / (theta(k,iw) * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))
         else
            exner = (press(k,iw) / p00) ** rocp   ! Defined WITHOUT CP factor
            temp = exner * theta(k,iw)
            rhovs = rhovsl(temp-273.15)
            sh_c(k,iw) = max(0.,sh_w(k,iw)-rhovs/real(rho(k,iw)))
            sh_v(k,iw) = sh_w(k,iw) - sh_c(k,iw)
! As is done for iteration in sub satadjst, use (0.3,0.7) weights to damp oscil.
!            rho(k,iw) = .5 * rho(k,iw)  &
!                      + .5 * press(k,iw) ** cvocp * p00k  &
            rho(k,iw) = press(k,iw) ** cvocp * p00k  &
               / (theta(k,iw) * (1. - sh_c(k,iw))  &
               * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))
            thil(k,iw) = theta(k,iw)  &
               / (1. + alvlocp * sh_c(k,iw) / max(temp,253.))
         endif

         if (k >= 2) then
            pkhyd = press(k-1,iw)  &
               - gravo2 * (rho(k-1,iw) * dzt(k-1) + rho(k,iw) * dzt(k))
! Impose minimum value of .01 Pa to avoid overshoot to negative values 
! during iteration.  Use weighting to damp oscillations
            press(k,iw) = .15 * press(k,iw) + .85 * max(.01,pkhyd)

         endif

      enddo
      
   enddo

enddo
call rsub('Wb',8)

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

! Average winds to U point and rotate at U point (assumed to be global simulation)

   do k = 1,mza
      ug = .5 * (uzonal(k,iw1) + uzonal(k,iw2))

      raxis = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)  ! dist from earth axis

      ugx = -ug * yeu(iu) / raxis 
      ugy =  ug * xeu(iu) / raxis 

      uc(k,iu) = ugx * unx(iu) + ugy * uny(iu)

      umc(k,iu) = uc(k,iu) * .5 * (rho(k,iw1) + rho(k,iw2))
   enddo

enddo
call rsub('Ub',8)

! Set UMC, UC to zero wherever ARU = 0.

call psub()
!----------------------------------------------------------------------
do iu = 2,mua
!----------------------------------------------------------------------
call qsub('U',iu)
   do k = 1,mza-1
      if (aru(k,iu) < 1.e-9) then
         umc(k,iu) = 0.
         uc(k,iu)  = 0.
      endif
   enddo
enddo
call rsub('Ub',0)

if (iparallel == 1) then
   call mpi_send_u('I')  ! Send U group
   call mpi_recv_u('I')  ! Recv U group
endif

! Set UMP to UMC

ump(:,:) = umc(:,:)

! print out initial state column from column 2

iw = 2

write(io6,*) ' '
write(io6,*) '========================================================================='
write(io6,*) '                    OLAM INITIAL STATE COLUMN (lhi)'
write(io6,*) '========================================================================='
write(io6,*) '   zm(m)    k     zt(m)   press(Pa)   rho(kg/m3)   theta(K)   sh_w(g/kg)'
write(io6,*) '========================================================================='
write(io6,*) ' '

do k = mza-1,2,-1
   write(io6, '(f10.2,1x,9(''-------''))') zm(k)
   write(io6, '(10x,i5,f10.2,f11.2,f11.4,f12.2,f11.4)')  &
       k,zt(k),press(k,iw),rho(k,iw),theta(k,iw),sh_w(k,iw)*1.e3
enddo
   
write(io6, '(f10.2,1x,9(''-------''))') zm(1)
write(io6,*) ' '

deallocate(zont,zonu,zonz,zonr  &
   ,zonp_vect,zont_vect,zonu_vect,zonz_vect,zonr_vect)

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_FLDSLHI_OUT,rst_buf)
#endif

return
end subroutine fldslhi

!===============================================================================

subroutine zonavg_init()

use misc_coms,   only: io6, idate1, imonth1, zonclim                    
use consts_coms, only: pio180, cp, rocp, eps_virt, omega, grav2, dlat
use mem_zonavg,  only: zonz, zonu, zont, zonr, zonp_vect,  &
                       alloc_zonavg, fill_zonavg

implicit none   

integer :: iplev,ilat,ilatn,ilats
real :: dyo2g,fcorn,fcors,pilo,pihi,thetavlo,thetavhi,cpo2g

! Read in 'ZONAVG_CLIMATE' dataset

call alloc_zonavg()
call fill_zonavg(zonclim,idate1,imonth1)

! In order to obtain zonally averaged water vapor mixing ratio, interpolate
! Mclatchy soundings in time to mclat array and obtain coefficients for 
! latitudinal spline interpolation

call fill_zonr_mclat()

! Compute zonal height field from hydrostatic and geostrophic balances

zonz(37,1) = 113. ! Std hgt of 1000 mb sfc in tropics (from Mclatchy sounding)
zonz(38,1) = 113. ! Std hgt of 1000 mb sfc in tropics (from Mclatchy sounding)

! Use geostrophic balance eqn to get 1000 mb heights at all latitudes

dyo2g = 2.5 * dlat / grav2  ! spacing (m) between zonavg values divided by 2g

do ilatn = 39,74
   ilats = 75 - ilatn 

! Coriolis parameter at midpoint between zonavg values (north and south hemispheres)
   
   fcorn = 2. * omega * sin(2.5 * (ilatn-38) * pio180)
   fcors = 2. * omega * sin(2.5 * (ilats-37) * pio180)

! Compute and apply delta zp using midpoint average of two zonu values

   zonz(ilatn,1) = zonz(ilatn-1,1)   &
      - fcorn * dyo2g * (zonu(ilatn-1,1) + zonu(ilatn,1))
   zonz(ilats,1) = zonz(ilats+1,1)  &
      - fcorn * dyo2g * (zonu(ilats+1,1) + zonu(ilats,1))
enddo

! Use hydrostatic equation to get all heights above 1000 mb surface

cpo2g = cp / grav2

do ilat = 1,74
   do iplev = 2,22
      pilo = (1.e-5 * zonp_vect(iplev-1)) ** rocp
      pihi = (1.e-5 * zonp_vect(iplev)) ** rocp
      thetavlo = zont(ilat,iplev-1) * (1. + eps_virt * zonr(ilat,iplev-1))  &
               / pilo
      thetavhi = zont(ilat,iplev) * (1. + eps_virt * zonr(ilat,iplev))  &
               / pihi
      zonz(ilat,iplev) = zonz(ilat,iplev-1) + cpo2g * (pilo - pihi)  &
         * (thetavlo + thetavhi)

!write(io6,3302) ilat,iplev,zonp_vect(iplev),pihi,thetavhi,zont(ilat,iplev)  &
!    ,zonr(ilat,iplev),zonz(ilat,iplev)
!3302 format('zonz1 ',2i6,12e15.6)
         
   enddo
enddo

return
end subroutine zonavg_init

!===============================================================================

subroutine fill_zonr_mclat()

use misc_coms,   only: io6, imonth1, idate1, iyear1
use mem_mclat,   only: slat, mclat, ypp_mclat, mcol, mclat_spline
use mem_zonavg,  only: zonr_vect, zonp_vect, zonr
use mem_radiate, only: jday

implicit none   

integer :: ilat,lv
real :: alat
integer, external :: julday

! Subroutine fill_zonr_mclat fills the ZONR array with specific humidity values 
! by means of interpolation from the McLatchy sounding.  The ZONR array is 2-D, 
! defined at 22 vertical pressure levels and 74 latitudes (2.5 degree spacing),
! and is required to supplement the zonu and zont arrays which are filled with
! monthly climatological data from the U.K. Met office.  

! Get current Julian day

jday = julday(imonth1,idate1,iyear1)

! Compute spline coefficients in preparation for interpolation

call mclat_spline(jday)  

! Loop over latitude columns of zonavg data

do ilat = 1,74
   alat = -93.75 + 2.5 * ilat

! Loop over vertical levels in McLatchy sounding

   do lv = 1,33

! Spline-interpolate (by latitude) Mclatchy pressure, vapor density, and 
! total density to current latitude
 
      call spline2(13,slat,mclat(1,lv,2),ypp_mclat(1,lv,2),alat,mcol(lv,2))
      call spline2(13,slat,mclat(1,lv,4),ypp_mclat(1,lv,4),alat,mcol(lv,4))
      call spline2(13,slat,mclat(1,lv,6),ypp_mclat(1,lv,6),alat,mcol(lv,6))

! Compute specific humidity from quotient of vapor and total density

      mcol(lv,4) = mcol(lv,4) / mcol(lv,6)
   enddo

! Vertically interpolate Mclatchy vapor specific humidity BY PRESSURE to zonavg 
! data levels

   call pintrp_ee(33,mcol(1,4),mcol(1,2),22,zonr_vect,zonp_vect)
   zonr(ilat,1:22) = zonr_vect(1:22)
enddo

return
end subroutine fill_zonr_mclat



