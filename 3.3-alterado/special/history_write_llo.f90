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
subroutine history_write_ll(nlon,nlat)

use mem_ijtabs,  only: itab_w, itab_u, rotation_angle
use mem_basic,   only: uc, wc, rho, press, theta

use mem_grid,    only: xew, yew, zew, &
                       topm, glatw, glonw, lpw, mza, mua, mwa, zm, zt, lpu,  &
                       xem, yem, zem, xeu, yeu, zeu,  &
                       unx, uny, unz, utx, uty, utz, dzt

use mem_addsc,   only: addsc

use misc_coms,   only: io6, time8, naddsc, timmax8, dtlong

use consts_coms, only: p00, rocp, piu180, grav, erad, pio180, cp, r8

!------------------------------------------
! Only for ncar test cases:
use ncar_testcases_all, only: ncar_testcase, ncar_choice
!------------------------------------------

implicit none

integer, intent(in) :: nlon,nlat
integer, save :: ncall = 0

!include 'netcdf.inc'
character(60), save :: nc_fname
character(20) :: zmodel = 'olam_dycore'
character(20) :: zversion  = ''
character(20) :: zres
character(1)  :: zcase, zchoice
character(4)  :: zscalar, zhgt
integer       :: ncid, iret
real(r8)      :: lat(nlat)
real(r8)      :: lon(nlon)
real(r8)      :: lev(mza-2)
real(r8)      :: cur_time

integer :: k,ka
real :: pressw
real :: vx,vy,vz,raxis,u,v
real :: tempk
real :: tuu1,tuu2,tuu3,tuu4,vc

real :: alon, alat
real :: fldval, phis, theta_bar

integer :: iu,iu1,iu2,iu3,iu4,im1,im2,im3,iw1,iw2,iw
integer :: ilat,ilon,ilatlon
integer :: kll

real :: scr(mza,mwa), scr1(mwa), scr2(mwa)

! Interpolate some model fields to uniform lat-lon grid with spacing 
! comparable to global grid (at equator).

real :: phis_ll(nlon,nlat)
real :: ps_ll  (nlon,nlat)

real :: u_ll (nlon,nlat,mza-2)
real :: v_ll (nlon,nlat,mza-2)
real :: w_ll (nlon,nlat,mza-2)
real :: t_ll (nlon,nlat,mza-2)
real :: p_ll (nlon,nlat,mza-2)
real :: r_ll (nlon,nlat,mza-2)

real :: z3_ll(nlon,nlat,mza-2)
real :: q1_ll(nlon,nlat,mza-2)
real :: q2_ll(nlon,nlat,mza-2)
real :: q3_ll(nlon,nlat,mza-2)
real :: q4_ll(nlon,nlat,mza-2)
real :: theta_pert_ll(nlon,nlat,mza-2)

real :: uzonal(mza,mua), umerid(mza,mua)

real, save, allocatable :: ubar_init(:,:)
real :: alat1,alat2,dnorm,dnorm1,dnorm2,anorm1,anorm2,ubar

real :: unorm1(1000),unorm2(1000),vctr18(1000)
real, save :: aspect = .7
real, save :: scalelab = .012
real, save :: timebeg,timeend,timedif,timeinc

! netcdf variable ids
integer, save :: time_id, PHIS_id, PS_id, U_id, V_id, W_id, T_id, P_id
integer, save :: R_id, Z3_id, Q1_id, Q2_id, Q3_id, Q4_id

do ilat = 1,nlat
   lat(ilat) = -90.0_r8 + 180.0_r8 * real(ilat-1,r8) / real(nlat-1,r8)
enddo

do ilon = 1,nlon
   lon(ilon) = 360.0_r8 * real(ilon-1,r8) / real(nlon,r8)
enddo

do k=2,mza-1
   lev(k-1) = zt(k)
enddo

cur_time = time8 / 86400.0_r8

scr  = 0.0
scr1 = 0.0
scr2 = 0.0

if (ncall == 0) then
   allocate(ubar_init(nlat,mza-1))

   timebeg = time8   / 86400.0_r8
   timeend = timmax8 / 86400.0_r8
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
      
   write(zcase,'(I1)') ncar_testcase
   write(zhgt, '(I0)') mza-2

   if (ncar_testcase <= 3) then

      if (rotation_angle < 1.0_r8) then
         zchoice = '0'
      elseif (rotation_angle < 46.0_r8) then
         zchoice = '3'
      else
         zchoice = '6'
      endif

   elseif (ncar_testcase <= 5) then
      zchoice = '0'
   else
      write(zchoice,'(I1)') ncar_choice
   endif

   if (ncar_testcase <= 2) then
      zscalar = '1234'
   elseif (ncar_testcase <= 3) then
      zscalar = '56'
   else
      zscalar = '0'
   endif
   
   if (nlon <= 92) then
      zres = 'low'
   elseif (nlon <= 182) then
      zres = 'medium'
   elseif (nlon <= 362) then
      zres = 'medium_high'
   elseif (nlon <= 722) then
      zres = 'high'
   else
      zres = 'ultra_high'
   endif

   if (len_trim(zversion) > 0) then
      nc_fname = trim(zmodel)//'_'//zcase//'-'//zchoice//'-'//trim(zscalar)// &
           '_'//trim(zres)//'_L'//trim(zhgt)//'_'//trim(zversion)//'.nc'
   else
      nc_fname = trim(zmodel)//'_'//zcase//'-'//zchoice//'-'//trim(zscalar)// &
           '_'//trim(zres)//'_L'//trim(zhgt)//'.nc'
   endif

!   call create_ncdf()

else

!   iret = nf_open(nc_fname, NF_WRITE, ncid)
!   iret = nf_enddef(ncid)

endif

!------------------------------------------------------------
! UZONAL and UMERID wind components at U point
!------------------------------------------------------------

do iu = 2,mua
   iu1 = itab_u(iu)%iu1
   iu2 = itab_u(iu)%iu2
   iu3 = itab_u(iu)%iu3
   iu4 = itab_u(iu)%iu4

   tuu1 = itab_u(iu)%tuu1
   tuu2 = itab_u(iu)%tuu2
   tuu3 = itab_u(iu)%tuu3
   tuu4 = itab_u(iu)%tuu4

   do k = 2,mza - 1
      kll = k - 1

      fldval = uc(k,iu1) * tuu1  &
             + uc(k,iu2) * tuu2  &
             + uc(k,iu3) * tuu3  &
             + uc(k,iu4) * tuu4

      vx = unx(iu) * uc(k,iu) + utx(iu) * fldval
      vy = uny(iu) * uc(k,iu) + uty(iu) * fldval
      vz = unz(iu) * uc(k,iu) + utz(iu) * fldval

      raxis = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)  ! dist from earth axis

      uzonal(k,iu) = (vy * xeu(iu) - vx * yeu(iu)) / raxis
      umerid(k,iu) = vz * raxis / erad  &
         - (vx * xeu(iu) + vy * yeu(iu)) * zeu(iu) / (raxis * erad) 

   enddo
enddo

print*, 'interpolating uzonal to lat-lon'

call interp_hnu_ll(nlon,nlat,mza,mza-2,uzonal,u_ll)

print*, 'interpolating umerid to lat-lon'

call interp_hnu_ll(nlon,nlat,mza,mza-2,umerid,v_ll)

!------------------------------------------------------------
! PHIS (topography height) and PS (surface pressure)
!------------------------------------------------------------

print*, 'interpolating uzonal and umerid to lat-lon'

do iw = 2,mwa
   im1 = itab_w(iw)%im1
   im2 = itab_w(iw)%im2
   im3 = itab_w(iw)%im3
   
   ka = lpw(iw)
   
   scr1(iw) = (topm(im1) + topm(im2) + topm(im3)) / 3.  ! PHIS

   scr2(iw) = press(ka,iw) + grav * rho(ka,iw) * (zt(ka) - scr1(iw))  ! PS
enddo

call interp_htw_ll(nlon,nlat,1,1,scr1,phis_ll)

call interp_htw_ll(nlon,nlat,1,1,scr2,ps_ll)

!------------------------------------------------------------
! Z3 height of model levels
!------------------------------------------------------------

print*, 'filling z3_ll values'

do k = 2,mza-1
   kll = k - 1
   do ilat = 1,nlat
      do ilon = 1,nlon
         z3_ll(ilon,ilat,kll) = zt(k)
      enddo
   enddo
enddo

!------------------------------------------------------------
! W vertical velocity at each K level
!------------------------------------------------------------

print*, 'interpolating w to lat-lon'

do iw = 2,mwa
   do k = 2,mza
      scr(k,iw) = .5 * (wc(k-1,iw) + wc(k,iw))
   enddo
enddo

call interp_htw_ll(nlon,nlat,mza,mza-2,scr,w_ll)

!------------------------------------------------------------
! Pressure
!------------------------------------------------------------

print*, 'interpolating pressure to lat-lon'

do iw = 2,mwa
   do k = 1,mza
      scr(k,iw) = press(k,iw)
   enddo
enddo

call interp_htw_ll(nlon,nlat,mza,mza-2,scr,p_ll)

!------------------------------------------------------------
! Density
!------------------------------------------------------------

print*, 'interpolating density to lat-lon'

do iw = 2,mwa
   do k = 1,mza
      scr(k,iw) = rho(k,iw)
   enddo
enddo

call interp_htw_ll(nlon,nlat,mza,mza-2,scr,r_ll)

!------------------------------------------------------------
! Temperature
!------------------------------------------------------------

print*, 'interpolating temperature to lat-lon'

do iw = 2,mwa
   do k = 1,mza
      scr(k,iw) = theta(k,iw) * (press(k,iw) / p00) ** rocp
   enddo
enddo

call interp_htw_ll(nlon,nlat,mza,mza-2,scr,t_ll)

!------------------------------------------------------------
! ADDSC_1
!------------------------------------------------------------

if (naddsc >= 1) then

   print*, 'interpolating addsc1 to lat-lon'

   call interp_htw_ll(nlon,nlat,mza,mza-2,addsc(1)%sclp,q1_ll)

endif

!------------------------------------------------------------
! ADDSC_2
!------------------------------------------------------------

if (naddsc >= 2) then

   print*, 'interpolating addsc2 to lat-lon'

   call interp_htw_ll(nlon,nlat,mza,mza-2,addsc(2)%sclp,q2_ll)

endif

!------------------------------------------------------------
! ADDSC_3
!------------------------------------------------------------

if (naddsc >= 3) then

   print*, 'interpolating addsc3 to lat-lon'

   call interp_htw_ll(nlon,nlat,mza,mza-2,addsc(3)%sclp,q3_ll)

endif

!------------------------------------------------------------
! ADDSC_4
!------------------------------------------------------------

if (naddsc >= 4) then

   print*, 'interpolating addsc4 to lat-lon'

   call interp_htw_ll(nlon,nlat,mza,mza-2,addsc(4)%sclp,q4_ll)

endif

!------------------------------------------------------------
! THETA perturbation
!------------------------------------------------------------

if (ncar_testcase == 6) then

   print*, 'interpolating theta perturbation to lat-lon'

   do k = 1,mza

      if (ncar_choice == 0) then
         theta_bar = 300. * exp(.0001 * zt(k) / grav)
      else
         theta_bar = 300. * exp(grav * zt(k) / (cp * 300.))
      endif

      do iw = 2,mwa
         scr(k,iw) = theta(k,iw) - theta_bar
      enddo
   enddo
   
   call interp_htw_ll(nlon,nlat,mza,mza-2,scr,theta_pert_ll)

endif

!------------------------------------------------------------
! L2 error norms of zonal wind
!------------------------------------------------------------

if (ncar_testcase == 1) then

   anorm1 = 0.
   anorm2 = 0.
   dnorm1 = 0.

   do ilat = 1,nlat
      alat = -90. + 180. * real(ilat-1) / real(nlat-1)
      alat1 = (max(-90.,alat - 90. / real(nlat-1))) * pio180
      alat2 = (min( 90.,alat + 90. / real(nlat-1))) * pio180

      do k = 2,mza-1
         kll = k - 1
         ubar = 0.

         do ilon = 1,nlon
            ubar = ubar + u_ll(ilon,ilat,kll)
         enddo

         ubar = ubar / real(nlon)

         if (ncall == 0) then
            ubar_init(ilat,kll) = ubar
         endif

         dnorm2 = 0.

         do ilon = 1,nlon
            dnorm = r_ll(ilon,ilat,kll) * dzt(k) * (sin(alat2) - sin(alat1))

            anorm1 = anorm1 + (u_ll(ilon,ilat,kll) - ubar)**2 * dnorm

            dnorm1 = dnorm1 + dnorm
            dnorm2 = dnorm2 + dnorm
         enddo
      
         anorm2 = anorm2 + (ubar - ubar_init(ilat,kll))**2 * dnorm2

      enddo
   enddo
   
endif

ncall = ncall + 1

if (ncar_testcase == 1) then
   unorm1(ncall) = sqrt(anorm1/dnorm1)
   unorm2(ncall) = sqrt(anorm2/dnorm1)
   vctr18(ncall) = time8 / 86400.

   print*, 'anorms_ll ',ncall,unorm1(ncall),unorm2(ncall)
endif

!------------------------------------------------------------
! Plot test
!------------------------------------------------------------

! Reopen the current graphics output workstation if it is closed

call o_reopnwk()

if (ncar_testcase <= 2) then
   call contplot_ll(nlon,nlat,mza-2,109,   u_ll,'Z',5 )
   call contplot_ll(nlon,nlat,mza-2,109,   v_ll,'Z',5 )
   call contplot_ll(nlon,nlat,mza-2, 26,   w_ll,'Z',5 )
   call contplot_ll(nlon,nlat,mza-2, 45,   t_ll,'Z',5 )
   call contplot_ll(nlon,nlat,mza-2, 41,   p_ll,'Z',5 )
   call contplot_ll(nlon,nlat,mza-2, 10,   r_ll,'Z',5 )
   call contplot_ll(nlon,nlat,1    , 41,  ps_ll,'Z',1 )
   call contplot_ll(nlon,nlat,1    , 84,phis_ll,'Z',1 )
endif

if (ncar_testcase == 3) then
!   if (naddsc >= 1) call contplot_ll(nlon,nlat,mza-2,95,q1_ll,'X',nlon/2+1 )
   if (naddsc >= 1) call contplot_ll(nlon,nlat,mza-2,107,q1_ll,'Y',(nlat+1)/2 )
   if (naddsc >= 2) call contplot_ll(nlon,nlat,mza-2,107,q2_ll,'Y',(nlat+1)/2 )
   if (naddsc >= 1) call contplot_ll(nlon,nlat,mza-2,107,q1_ll,'Z',24 )
   if (naddsc >= 2) call contplot_ll(nlon,nlat,mza-2,107,q2_ll,'Z',24 )
endif

if (ncar_testcase == 5) then
   call contplot_ll(nlon,nlat,mza-2,100,w_ll,'Y',(3*nlat+1)/4 )
endif

if (ncar_testcase == 6) then
   if (ncar_choice <= 2) then
      call contplot_ll(nlon,nlat,mza-2,104,theta_pert_ll,'Y',(nlat+1)/2 )
      call contplot_ll(nlon,nlat,mza-2,104,theta_pert_ll,'Z',mza/2 )
   else
      call contplot_ll(nlon,nlat,mza-2,104,theta_pert_ll,'Y',(3*nlat+1)/4 )
   endif
endif

if (time8 + 1.5 * dtlong > timmax8) then

! Plot time series

   if (ncar_testcase == 1) then
!----------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab   &
                    ,ncall,  vctr18,unorm1            &
                    ,'time(days)','UNORM1' &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!----------------------------------------------------------------------
      call plotback()
      call oplot_xy2('0','N',aspect,scalelab   &
                    ,ncall,  vctr18,unorm2            &
                    ,'time(days)','UNORM2' &
                    ,timebeg,timeend,timeinc,5  ,0.,12.,.5,4  )
      call o_frame()
!-------------------------------------------------------------------
   endif

endif

call o_clswk()

! CALL NetCDF output routine here
!call netcdf_output()

! close NetCDF file
!iret = nf_close(ncid)

return
end subroutine history_write_ll

!==========================================================================

subroutine tileplot_ll(nlon,nlat,itab,fld)

use plotcolors,  only: clrtab
use oplot_coms, only: op

implicit none

integer, intent(in) :: nlon,nlat,itab
real, intent(in) :: fld(nlon,nlat)

real :: alon, alat
real :: fldval

real :: htpn(4)
real :: vtpn(4)
integer :: ilat, ilon, ival, icolor

call plotback()

call o_set(0.,1.,0.,1.,-5.,365.,-95.,95.,1)

op%xmin =  -5.
op%xmax = 365.
op%ymin = -95.
op%ymax =  95.

do ilat = 1,nlat
   alat = -90. + 180. * real(ilat-1) / real(nlat-1)

   vtpn(1) = alat - 90. / real(nlat-1)
   vtpn(2) = vtpn(1)
   vtpn(3) = alat + 90. / real(nlat-1)
   vtpn(4) = vtpn(3)

   do ilon = 1,nlon
      alon = 0. + 360. * real(ilon-1) / real(nlon)

      htpn(1) = alon - 180. / real(nlon)
      htpn(2) = alon + 180. / real(nlon)
      htpn(3) = htpn(2)
      htpn(4) = htpn(1)

! Set field value and extract contour color from color table

      fldval = fld(ilon,ilat)

      ival = 1
      do while (fldval > clrtab(itab)%vals(ival) .and.  &
                  ival < clrtab(itab)%nvals             )
         ival = ival + 1
      enddo
      icolor = clrtab(itab)%ipal(ival)
      
      call fillpolyg(4,htpn,vtpn,icolor)

   enddo
enddo

call o_frame()

return
end subroutine tileplot_ll

!==========================================================================

subroutine contplot_ll(nlon,nlat,nlev,itab,fld,slab,islab)

use plotcolors, only: clrtab
use oplot_coms, only: op
use mem_grid,   only: zm,zt,mza

implicit none

integer, intent(in) :: nlon,nlat,nlev,itab,islab
real, intent(in) :: fld(nlon,nlat,nlev)
character, intent(in) :: slab*1

real :: fldvals(4)

real :: htpn(4)
real :: vtpn(4)
integer :: ilat, ilon, ilon1, ilon2, ival, icolor, ilev

if (trim(slab) == 'Y') then

do ilev = 1,nlev
do ilat = 1,nlat
do ilon = 1,nlon

!write(6,442) 'ccpgA ',ilev,ilat,ilon,fld(ilon,ilat,ilev)
!442 format(a,3i5,f10.4)

enddo
enddo
enddo

endif

call plotback()

if (trim(slab) == 'X') then

   op%xmin = -95.
   op%xmax =  95.
   op%ymin =   0.
   op%ymax = zm(mza-1)

   call o_set(.05,.95,.05,.95,op%xmin,op%xmax,op%ymin,op%ymax,1)

   do ilat = 1,nlat-1
      htpn(1) = -90. + 180. * real(ilat-1) / real(nlat-1)
      htpn(2) = -90. + 180. * real(ilat  ) / real(nlat-1)
      htpn(3) = htpn(2)
      htpn(4) = htpn(1)

      do ilev = 1,nlev-1

         vtpn(1) = zt(ilev+1)
         vtpn(2) = vtpn(1)
         vtpn(3) = zt(ilev+2)
         vtpn(4) = vtpn(3)

! Set field value and extract contour color from color table

         fldvals(1) = fld(islab,ilat  ,ilev  )
         fldvals(2) = fld(islab,ilat+1,ilev  )
         fldvals(3) = fld(islab,ilat+1,ilev+1)
         fldvals(4) = fld(islab,ilat  ,ilev+1)

         call contpolyg(itab,1,4,htpn,vtpn,fldvals)
         call contpolyg(itab,0,4,htpn,vtpn,fldvals)

      enddo
   enddo

elseif (trim(slab) == 'Y') then

   op%xmin =  -5.
   op%xmax = 365.
   op%ymin =   0.
   op%ymax = zm(mza-1)

   call o_set(.05,.95,.05,.95,op%xmin,op%xmax,op%ymin,op%ymax,1)

   do ilev = 1,nlev-1
      vtpn(1) = zt(ilev+1)
      vtpn(2) = vtpn(1)
      vtpn(3) = zt(ilev+2)
      vtpn(4) = vtpn(3)

      do ilon1 = 1,nlon
         ilon2 = ilon1 + 1
         if (ilon1 == nlon) ilon2 = 1

         htpn(1) = 0. + 360. * real(ilon1-1) / real(nlon)
         htpn(2) = 0. + 360. * real(ilon1  ) / real(nlon)
         htpn(3) = htpn(2)
         htpn(4) = htpn(1)

! Set field value and extract contour color from color table

         fldvals(1) = fld(ilon1,islab,ilev  )
         fldvals(2) = fld(ilon2,islab,ilev  )
         fldvals(3) = fld(ilon2,islab,ilev+1)
         fldvals(4) = fld(ilon1,islab,ilev+1)


!write(6,443) 'ccpgY ',ilev,ilon1,htpn(1:4),vtpn(1:4),fldvals(1:4)
!443 format(a,2i5,12f9.3)

         call contpolyg(itab,1,4,htpn,vtpn,fldvals)
         call contpolyg(itab,0,4,htpn,vtpn,fldvals)

      enddo
   enddo

elseif (trim(slab) == 'Z') then

   op%xmin =  -5.
   op%xmax = 365.
   op%ymin = -95.
   op%ymax =  95.

   call o_set(.05,.95,.05,.95,op%xmin,op%xmax,op%ymin,op%ymax,1)

   do ilat = 1,nlat-1
      vtpn(1) = -90. + 180. * real(ilat-1) / real(nlat-1)
      vtpn(2) = vtpn(1)
      vtpn(3) = -90. + 180. * real(ilat  ) / real(nlat-1)
      vtpn(4) = vtpn(3)

      do ilon1 = 1,nlon
         ilon2 = ilon1 + 1
         if (ilon1 == nlon) ilon2 = 1
      
         htpn(1) = 0. + 360. * real(ilon1-1) / real(nlon)
         htpn(2) = 0. + 360. * real(ilon1  ) / real(nlon)
         htpn(3) = htpn(2)
         htpn(4) = htpn(1)

! Set field value and extract contour color from color table

         fldvals(1) = fld(ilon1,ilat  ,islab)
         fldvals(2) = fld(ilon2,ilat  ,islab)
         fldvals(3) = fld(ilon2,ilat+1,islab)
         fldvals(4) = fld(ilon1,ilat+1,islab)

         call contpolyg(itab,1,4,htpn,vtpn,fldvals)
         call contpolyg(itab,0,4,htpn,vtpn,fldvals)

      enddo
   enddo

endif
   
call o_frame()

return
end subroutine contplot_ll
