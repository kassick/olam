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
subroutine gridinit()

use misc_coms,   only: io6, runtype, mdomain, itopoflg, ngrids, initial,  &
                       nzp, timmax8, alloc_misc, iparallel,               &
                       iyear1, imonth1, idate1, itime1, s1900_init, s1900_sim

use leaf_coms,   only: nzg, nzs, isfcl, nwl
use sea_coms,    only: nws

use mem_ijtabs,  only: istp, mrls, fill_jtabs
use oplot_coms,  only: op
use mem_grid,    only: nma, nua, nwa, mma, mua, mwa, mza, nwa, zm, zt,  &
                       alloc_grid1, alloc_grid2
use mem_nudge,   only: nudflag
use mem_sflux,   only: init_fluxcells, nlandflux, nseaflux

implicit none

real, external :: walltime
real :: wtime_start_gridinit,inicio

wtime_start_gridinit = walltime(0.)

! Read LAND and SEA files

if (runtype /= 'MAKESFC' .and. isfcl == 1) then

   write(io6,'(/,a)') 'gridinit calling landfile_read'

   inicio = walltime(wtime_start_gridinit)

   call landfile_read()
   write(io6, *) '      ==T== Tempo total gasto na landfile_read(): ',(walltime(wtime_start_gridinit)-inicio)

   write(io6,'(/,a)') 'gridinit calling seafile_read'
  
   inicio = walltime(wtime_start_gridinit)
   
   call seafile_read()
   write(io6, *) '      ==T== Tempo total gasto na seafile_read(): ',(walltime(wtime_start_gridinit)-inicio)

endif

! Generate OLAM grid structure for 'MAKESFC' or 'MAKEGRID' runtype

if (runtype == 'MAKESFC' .or. runtype == 'MAKEGRID') then

! Vertical grid coordinate setup

      write(io6,'(/,a)') 'gridinit calling gridset'
      call gridset()

! Horizontal grid setup

   if (mdomain == 0) then

      write(io6,'(/,a)') 'gridinit calling icosahedron'
      call icosahedron()     ! global spherical domain; calls 2 allocs

   elseif (mdomain == 2 .or. mdomain == 3) then

      write(io6,'(/,a)') 'gridinit calling cartesian_2d'
      call cartesian_2d()    ! 2D cartesian channel domain; calls 2 allocs

   elseif (mdomain == 4) then

      write(io6,'(/,a)') 'gridinit calling cartesian_3d'
      call cartesian_3d()    ! 3D cartesian channel domain; calls 2 allocs

   endif

   write(io6,'(/,a)') 'gridinit after icosahedron or cartesian'
   write(io6,'(a,i8)')    ' mma = ',mma
   write(io6,'(a,i8)')    ' mua = ',mua
   write(io6,'(a,i8)')    ' mwa = ',mwa

   if (ngrids > 1 .and. mdomain /= 2 .and. mdomain /= 3) then

      write(io6,'(/,a,i5)') 'gridinit calling spawn_nest; ngrids = ',ngrids

      call spawn_nest()  ! calls 2 allocs

      write(io6,'(/,a)') 'olam_run finished spawn_nest'
      write(io6,'(a,i8)')   ' mma = ',mma
      write(io6,'(a,i8)')   ' mua = ',mua
      write(io6,'(a,i8)')   ' mwa = ',mwa

   endif

! Allocate (xeu, yeu, zeu, xew, yew, zew, glatw, glonw) for full domain

   write(io6,'(/,a)') 'gridinit calling alloc_grid1 for full domain'

   call alloc_grid1()

! special - perform cvt grid adjustment
!   call pcvt()
! end special

! Fill (xeu, yeu, zeu, xew, yew, zew, glatw, glonw) for full domain

   write(io6,'(/,a)') 'gridinit calling triangle_geometry for full domain'

   call triangle_geometry()

! 'MAKESFC' run uses OLAM atmos grid configuration for land/sea cells

   if (runtype == 'MAKESFC') then
      write(io6,'(/,a)') 'gridinit calling makesfc'
      call makesfc()
      return
   endif

! Initialize dtlm, dtsm, ndtrat, and nacoust, 
! and compute the timestep schedule for all grid operations.

   write(io6,'(/,a)') 'gridinit calling modsched'

   call modsched()

   write(io6,'(/,a)') 'gridinit calling fill_jtabs'

   call fill_jtabs(mma,mua,mwa)

! Allocate remaining unstructured grid geometry arrays

   write(io6,'(/,a)') 'gridinit calling alloc_grid2'

   call alloc_grid2()

! Fill topography values

   write(io6,'(/,a)') 'gridinit calling mksfc_topo'

   call topo_init()  ! Initialize TOPM from default value (may customize here)
   
   if (itopoflg == 1) then
      call topo_database_read() ! Read TOPM from dataset
   endif

! Set up control volumes

   write(io6,'(/,a)') 'gridinit calling ctrlvols'

   call ctrlvols()        ! Fill unstructured grid control volumes

! Check whether LEAF and SEA are being used

   if (isfcl == 1) then

! Determine and initialize flux cells for entire model domain

      write(io6,'(/,a)') 'gridinit after alloc(landflag,seaflag)'

      write(io6,'(/,a)') 'gridinit calling init_fluxcells'
      call init_fluxcells()

   endif

! Write GRIDFILE

   write(io6,'(/,a)') 'gridinit calling gridfile_write'
   call gridfile_write()

else

! Read atmos grid for INITIAL/HISTORY/PLOTONLY/PARCOMBINE run

   write(io6,'(/,a)') 'gridinit calling gridfile_read'

   inicio = walltime(wtime_start_gridinit)

   call gridfile_read()
   write(io6, *) '      ==T== Tempo total gasto na gridfile_read(): ',(walltime(wtime_start_gridinit)-inicio)

endif

if (isfcl == 1) then
   write(io6,'(/,a)') 'gridinit before return to olam_run'
   write(io6,'(a,i8)')   ' nwl       = ',nwl
   write(io6,'(a,i8)')   ' nws       = ',nws
   write(io6,'(a,i8)')   ' nlandflux = ',nlandflux
   write(io6,'(a,i8)')   ' nseaflux  = ',nseaflux
endif

return
end subroutine gridinit

!===============================================================================

subroutine gridset()

use mem_grid,   only: nza, mza,  &
                      zm, zt, dzm, dzt, dzim, dzit, alloc_gridz
use misc_coms,  only: io6, nzp, deltaz, dzrat, dzmax, ztop, zbase
use oname_coms, only: nl

implicit none

integer :: ifm,icm,k,iinc,icnt,if1,jinc  &
   ,jcnt,jf,kinc,kcnt,kf,nrat,i,j,kcy,kcw,kk,ng,npg
integer :: nidiag,ijcorner,numd
real :: centx1,centy1,centx,centy,dzr,dsum,dzrcm,dzrfm,tsum,dzrati
real :: dxmax
real, allocatable, dimension(:) :: zmvec,ztvec

nza = nzp
mza = nzp

call alloc_gridz()
allocate (zmvec(-1:mza+1),ztvec(-1:mza+1))

! calculate zm

if ( deltaz < spacing(0.) ) then
   zmvec(1:nzp) = nl%zz(1:nzp)
   zmvec(nzp+1) = 2. * zmvec(nzp) - zmvec(nzp-1)
else
   zmvec(1) = zbase
   zmvec(2) = zbase + deltaz
   do k = 3,nzp+1
      zmvec(k) = zmvec(k-1)  &
               + min(dzrat * (zmvec(k-1) - zmvec(k-2)),max(deltaz,dzmax))
   enddo
endif
dzrati = (zmvec(2) - zmvec(1)) / (zmvec(3) - zmvec(2))
zmvec(0) = zmvec(1) - (zmvec(2) - zmvec(1)) * dzrati
zmvec(-1) = zmvec(0) - (zmvec(1) - zmvec(0)) * dzrati

ztop = zmvec(nzp-1)

! compute zt values by geometric interpolation.

do k = 1,nza
   dzrfm = sqrt(sqrt((zmvec(k+1) - zmvec(k)) / (zmvec(k-1) - zmvec(k-2))))
   ztvec(k) = zmvec(k-1) + (zmvec(k) - zmvec(k-1)) / (1. + dzrfm)
enddo
ztvec(nza+1) = .5 * (zmvec(nza) + zmvec(nza+1))

! Other vertical coordinate values

do k = 1,nza
   zm(k)   = zmvec(k)
   zt(k)   = ztvec(k)
   dzm(k)  = ztvec(k+1) - ztvec(k)
   dzt(k)  = zmvec(k) - zmvec(k-1)
   dzim(k) = 1. / dzm(k)
   dzit(k) = 1. / dzt(k)
enddo

deallocate (zmvec,ztvec)

return
end subroutine gridset

!===============================================================================

subroutine ctrlvols()

use mem_ijtabs,  only: jtab_m, jtab_u, jtab_w, itab_u, itab_w
use misc_coms,   only: io6, mdomain
use consts_coms, only: erad
use mem_grid,    only: nsw_max, nza, nma, nua, nwa, lpu, lcu, lpw, lsw,  &
                       topm, zm, dzt, zt,  &
                       dtu, arw0, aru, arw, volt, volti, volui, volwi

implicit none

integer :: j,iw,iwp,iter,iu,iup,im1,im2,k,km,i,im  &
   ,im11,im21,im12,im22,iu1,iu2,iw1,iw2,kp  &
   ,iu1a,iu1b,iu2a,iu2b  &
   ,iuo1a,iuo1b,iuo2a,iuo2b  &
   ,im3,iu3,ka

integer, parameter :: npass = 2
integer :: ipass

real :: hmin,hmax,arwo4,sum1,sumk  &
   ,arwo3,w1,w2,zfac,t1,t2,t3,hm

!!!!!!!!!!!!! special quadrature parameters

integer, parameter :: np = 10
integer :: ip,jp,ipair

real :: top1,top2,top3,top4
real :: topp,topm13,topm32
real :: facj1,facj2,faci1,faci2
real :: zfacm,zfact,del_arw

!!!!!!!!!!!!! end special

! Define face areas and volumes for all cells
! This algorithm applies to regular input topography defined at M points,
!    and must be replaced or modified for more complex boundary topologies
! Small corrections to topography are applied wherever necessary to 
!    prevent very small apertures

! Copy topm to lateral boundary points

!call psub()
!----------------------------------------------------------------------
!do j = 1,jtab_m(2)%jend(1); im = jtab_m(2)%im(j)
!   itopm = itab_m(im)%itopm
!----------------------------------------------------------------------
!call qsub('M',im)

!   topm(im) = topm(itopm)

!enddo
!call rsub('M',2)

! ARU, LPU, and TOPM adjustment

topm(2:nma) = max(topm(2:nma),zm(1))

do ipass = 1,npass

   write(io6,*) 'Defining contol volume areas: ipass = ',ipass

   call psub()
!----------------------------------------------------------------------
   do j = 1,jtab_u(1)%jend(1); iu = jtab_u(1)%iu(j)
      im1 = itab_u(iu)%im1; im2 = itab_u(iu)%im2
!----------------------------------------------------------------------
   call qsub('U',iu)

      aru(1,iu) = 0.
      lpu(iu) = nza

      if (dtu(iu) < 1.e-6) then

         do k = 2,nza
            aru(k,iu) = 0.
         enddo

      else

         if (topm(im2) > topm(im1)) then
            hmin = topm(im1)
            hmax = topm(im2)
         else
            hmin = topm(im2)
            hmax = topm(im1)
         endif

         do k = nza,2,-1
            km = k - 1

            if (zm(k) <= hmin) then

               aru(k,iu) = 0.

            elseif (zm(km) >= hmax) then

               aru(k,iu) = dtu(iu) * dzt(k)

            elseif(zm(k) < hmax .and. zm(km) < hmin) then

               aru(k,iu) = dtu(iu) * .5 * (zm(k) - hmin)**2 / (hmax - hmin)
                              
               if (ipass < npass .and. aru(k,iu) < .01 * dtu(iu) * dzt(k)) then
                  topm(im1) = max(topm(im1),zm(k) + .003)
                  topm(im2) = max(topm(im2),zm(k) + .003)
                  aru(k,iu) = 0.
               endif
               
            elseif(zm(k) <  hmax .and. zm(km) >=  hmin) then

               aru(k,iu) = dtu(iu) * dzt(k)  &
                  * (.5 * (zm(k) + zm(km)) - hmin) / (hmax - hmin)

            elseif(zm(k) >= hmax .and. zm(km) < hmin) then

               aru(k,iu) = dtu(iu) * (zm(k) - .5 * (hmax + hmin))
               
               if (ipass < npass .and. aru(k,iu) < .01 * dtu(iu) * dzt(k)) then
                  topm(im1) = max(topm(im1),zm(k))
                  topm(im2) = max(topm(im2),zm(k))
                  aru(k,iu) = 0.
               endif
               
            elseif(zm(k) >= hmax .and. zm(km) >=  hmin) then

               aru(k,iu) = dtu(iu)  &
                  * (dzt(k) - .5 * (hmax - zm(km)) ** 2 / (hmax - hmin))

            else

               write(io6,*) 'aru option not reached ',k,i,j,  &
                  zm(k),zm(km),hmax,hmin
               stop 'stop aru defn'   

            endif

            if (aru(k,iu) > .005 * dtu(iu) * dzt(k)) lpu(iu) = k

         enddo

      endif

! Expand ARU with height for spherical geometry

      if (mdomain < 2) then
         do k = 2,nza
            zfac = (erad + zt(k)) / erad
            aru(k,iu) = aru(k,iu) * zfac
         enddo
      endif

   enddo
   call rsub('U',1)

enddo  ! end ipass loop  
      
! Set ARU to zero at non-topo walls

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_u(3)%jend(1); iu = jtab_u(3)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)

   do k = 2,nza
      aru(k,iu) = 0.
   enddo
   lpu(iu) = nza 
      
enddo
call rsub('U',3)

! Initialize nsw_max

nsw_max = 1

! Fill arw(k,iw) from topography

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(1)%jend(1); iw = jtab_w(1)%iw(j)
   im1 = itab_w(iw)%im1; im2 = itab_w(iw)%im2
   im3 = itab_w(iw)%im3
!----------------------------------------------------------------------
call qsub('W',iw)

   t1 = topm(im1)
   t2 = topm(im2)
   t3 = topm(im3)

   lsw(iw) = 1 ! # levels of surface W

   do k = 2,nza-2
      hm = zm(k)

      if     (t1 >= hm .and. t2 >= hm .and. t3 >= hm) then
         arw(k,iw) = 0.
      elseif (t1 <  hm .and. t2 <  hm .and. t3 <  hm) then
         arw(k,iw) = arw0(iw) 
      elseif (t1 >= hm .and. t2 >= hm .and. t3 <  hm) then
         arw(k,iw) = arw0(iw) * (hm - t3)**2 / ((t1 - t3) * (t2 - t3))
         lsw(iw) = lsw(iw) + 1
      elseif (t1 >= hm .and. t2 <  hm .and. t3 >= hm) then
         arw(k,iw) = arw0(iw) * (hm - t2)**2 / ((t1 - t2) * (t3 - t2))
         lsw(iw) = lsw(iw) + 1
      elseif (t1 <  hm .and. t2 >= hm .and. t3 >= hm) then
         arw(k,iw) = arw0(iw) * (hm - t1)**2 / ((t2 - t1) * (t3 - t1))
         lsw(iw) = lsw(iw) + 1
      elseif (t1 >= hm .and. t2 <  hm .and. t3 <  hm) then
         arw(k,iw) = arw0(iw) * (1. - (t1 - hm)**2 / ((t1 - t2) * (t1 - t3)))
         lsw(iw) = lsw(iw) + 1
      elseif (t1 <  hm .and. t2 >= hm .and. t3 <  hm) then
         arw(k,iw) = arw0(iw) * (1. - (t2 - hm)**2 / ((t2 - t1) * (t2 - t3)))
         lsw(iw) = lsw(iw) + 1
      elseif (t1 <  hm .and. t2 <  hm .and. t3 >= hm) then
         arw(k,iw) = arw0(iw) * (1. - (t3 - hm)**2 / ((t3 - t1) * (t3 - t2)))
         lsw(iw) = lsw(iw) + 1
         cycle
      endif

   enddo

! Set arw = 0 for bottom (k = 1), wall-on-top (k = nza-1), 
! and top (k = nza) levels

   arw(1,iw) = 0.   
   arw(nza-1,iw) = 0.   
   arw(nza,iw) = 0.   

! Increase nsw_max if necessary

   if (lsw(iw) > nsw_max) nsw_max = lsw(iw)
   
! Expand ARW with height for spherical geometry

   if (mdomain < 2) then
      do k = 2,nza-2
         zfac = (erad + zm(k)) / erad
         arw(k,iw) = arw(k,iw) * zfac ** 2
      enddo
   endif

enddo
call rsub('W',1)

! VOLT and VOLTI from ARU/ARW

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(3)%jend(1); iw = jtab_w(3)%iw(j)
   im1 = itab_w(iw)%im1; im2 = itab_w(iw)%im2; im3 = itab_w(iw)%im3
   iu1 = itab_w(iw)%iu1; iu2 = itab_w(iw)%iu2; iu3 = itab_w(iw)%iu3
!----------------------------------------------------------------------
call qsub('W',iw)
 
   t1 = topm(im1)
   t2 = topm(im2)
   t3 = topm(im3)
   
   do k = 2,nza

      zfacm = 1.
      zfact = 1.

! Expand ARW with height for spherical geometry

      if (mdomain < 2) then
         zfacm = (erad + zm(k)) / erad
         zfact = (erad + zt(k)) / erad
      endif

      if (t1 < zm(k-1) .and. t2 < zm(k-1) .and. t3 < zm(k-1)) then
!         arw(k,iw) = arw0(iw) * zfacm**2
         volt(k,iw) = arw0(iw) * zfact**2 * dzt(k)
      else
!         arw(k,iw) = 0.
         volt(k,iw) = 1.e-9
 
! USE NUMERICAL QUADRATURE APPROACH

! FIND INTERPOLATED TOPOGRAPHY HEIGHT AT MULTIPLE POINTS WITHIN IW TRIANGLE   

! special - find xe,ye,ze of center of np*np sub-triangles in IW

         topm13 = topm(im3) - topm(im1)
         topm32 = topm(im2) - topm(im3)

         del_arw = arw0(iw) / real(np*np)

         do jp = 1,np

            facj1 = real(jp-1) / real(np)
            facj2 = real(jp)   / real(np)

            do ip = 1,jp

               faci1 = real(ip-1) / real(np)
               faci2 = real(ip)   / real(np)

               top1 = topm(im1) + facj1 * topm13 + faci1 * topm32
               top2 = topm(im1) + facj1 * topm13 + faci2 * topm32
               top3 = topm(im1) + facj2 * topm13 + faci1 * topm32
               top4 = topm(im1) + facj2 * topm13 + faci2 * topm32

               do ipair = 1,2

                  if (ipair == 1) then
                     topp = (top1 + top3 + top4) / 3.
                  else
                     topp = (top1 + top2 + top4) / 3.
                  endif

                  if (topp < zm(k)) then
!                     arw(k,iw) = arw0(iw) + del_arw * zfacm**2
                     volt(k,iw) = volt(k,iw) + del_arw * zfact**2  &
                                * (zm(k) - max(topp,zm(k-1)))
                  endif

                  if (ip == jp) exit 

               enddo

            enddo

         enddo

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!go to 1
! Option for stability: expand volt if too small relative to any grid cell face

         volt(k,iw) = max(volt(k,iw)                             &
            ,real(arw(k,iw),8)   * dzt(k)                                &
            ,real(arw(k-1,iw),8) * dzt(k)                                &
            ,real(arw0(iw),8) * zfac * aru(k,iu1) / max(1.e-3,dtu(iu1))  &
            ,real(arw0(iw),8) * zfac * aru(k,iu2) / max(1.e-3,dtu(iu2))  &
            ,real(arw0(iw),8) * zfac * aru(k,iu3) / max(1.e-3,dtu(iu3))  )
1 continue
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      endif

      volti(k,iw) = 1. / volt(k,iw)
   enddo

enddo
call rsub('W',3)

! LPW

lpw(1:nwa) = nza

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(6)%jend(1); iw = jtab_w(6)%iw(j)
   iu1 = itab_w(iw)%iu1;  iu2 = itab_w(iw)%iu2;  iu3 = itab_w(iw)%iu3
!----------------------------------------------------------------------
call qsub('W',iw)

   lpw(iw) = nza-1   	
   do k = nza-2,2,-1
      if (lpu(iu1) <= k .or. lpu(iu2) <= k .or. lpu(iu3) <= k) then
         lpw(iw) = k
      else
         arw(k,iw) = 0.  ! close area if all 3 surrounding U's are closed
         volt(k,iw) = 1.e-9  ! close volume if all 3 surrounding U's are closed
         volti(k,iw) = 1.e9  ! close volume if all 3 surrounding U's are closed
      endif
   enddo

enddo
call rsub('W',6)

! VOLUI from VOLT

lcu(1:nua) = nza

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_u(4)%jend(1); iu = jtab_u(4)%iu(j)
   iw1 = itab_u(iu)%iw1; iw2 = itab_u(iu)%iw2
!----------------------------------------------------------------------
call qsub('U',iu)

   do k = nza,2,-1
      volui(k,iu) = 1. / (volt(k,iw1) + volt(k,iw2))      

      if (volt(k,iw1) + volt(k,iw2) > 1.e-8) then
         lcu(iu) = k
      endif
   enddo

   volui(1,iu) = 1.e9

enddo
call rsub('U',4)

! VOLWI from VOLT

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(5)%jend(1); iw = jtab_w(5)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   do k = 2,nza
      kp = min(k+1,nza)
      volwi(k,iw) = 2. / (volt(k,iw) + volt(kp,iw)) 
   enddo

! modify volwi for lpw and lpw-1 levels

   ka = lpw(iw)
   volwi(ka,iw)   = 1. / (volt(ka,iw) + .5 * volt(ka+1,iw))
   volwi(ka-1,iw) = 1. / (.5 * volt(ka,iw))

enddo
call rsub('W',5)

return
end subroutine ctrlvols

!===============================================================================

subroutine topo_init()

use mem_grid, only: nma, topm, xem
use misc_coms, only: io6, deltax

implicit none

integer :: im

real :: hfwid
real :: hgt
real :: hfwid2

! Fill the TOPM array with a default value of 0 or modify it as desired.  
! If itopoflg is set to 1, these values will be overridden in the call to
! topo_database, which inputs a standard OLAM topography dataset.

hfwid = 10000.

! dudhia expts
! hfwid = 5. * .866 * deltax

! hgt = 405.
! hgt = 1012.
! end dudhia expts

hfwid2 = hfwid**2

do im = 2,nma
   topm(im) = 0.

!   topm(im) = 200. * mod(im,4)

! SPECIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   topm(im) = max(0.,hgt * hfwid2 / (hfwid2 + xem(im)**2) - 1.)  
!   write(io6,*) 'topm ',im,xem(im),topm(im)
! TOPM = 0 AT LARGE DISTANCE FROM HILL
! END SPECIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

enddo

return
end subroutine topo_init

!===============================================================================

subroutine gridfile_write()

use misc_coms,  only: io6, ngrids, gridfile, mdomain, nzp, nxp,  &
                      iclobber, itopoflg,  &
                      deltax, deltaz, dzmax, dzrat, zbase,  &
                      centlat, centlon, grdlen, grdwid, grdaxis
use mem_ijtabs, only: nloops_m, nloops_u, nloops_w, maxtpn, mrls,  &
                      itab_m, itab_u, itab_w
use mem_grid,   only: nza, nma, nua, nwa, nsw_max,  &
                      zm, zt, dzm, dzt, dzim, dzit,  &
                      lpu, lcu, lpw, lsw,  &
                      topm, xem, yem, zem, xeu, yeu, zeu, xew, yew, zew,  &
                      unx, uny, unz, utx, uty, utz, wnx, wny, wnz,  &
                      dtu, dnu, dniu, arw0, glatw, glonw, glatm, glonm,  &
                      aru, volui, arw, volwi, volt, volti,  &
                      alloc_grid1, alloc_grid2
use leaf_coms,  only: isfcl
use mem_sflux,  only: nseaflux, nlandflux, seaflux, landflux, ntraps, ntrapl,  &
                      xemstrap, yemstrap, zemstrap, xemltrap, yemltrap, zemltrap

use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
use mem_nudge,  only: mnudp

implicit none

! This routine writes the grid variables to the grid file.

integer :: im, iu, iw

integer :: ndims, idims(2)

! Scratch arrays for copying output

integer, allocatable :: iscr(:,:)
real,    allocatable :: rscr(:,:)
logical, allocatable :: lscr(:,:)

write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
write(io6,*) 'grid_write: opening file:', trim(gridfile)
write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

call shdf5_open(gridfile,'W',iclobber)

! Write the gridfile information that exists in namelist

ndims = 1
idims(1) = 1
idims(2) = 1

call shdf5_orec(ndims, idims, 'NZP'     , ivars=nzp)
call shdf5_orec(ndims, idims, 'NXP'     , ivars=nxp)
call shdf5_orec(ndims, idims, 'MDOMAIN' , ivars=mdomain)
call shdf5_orec(ndims, idims, 'NGRIDS'  , ivars=ngrids)
call shdf5_orec(ndims, idims, 'ISFCL'   , ivars=isfcl)
call shdf5_orec(ndims, idims, 'ITOPOFLG', ivars=itopoflg)
call shdf5_orec(ndims, idims, 'DELTAX'  , rvars=deltax)
call shdf5_orec(ndims, idims, 'DELTAZ'  , rvars=deltaz)
call shdf5_orec(ndims, idims, 'DZRAT'   , rvars=dzrat)
call shdf5_orec(ndims, idims, 'DZMAX'   , rvars=dzmax)
call shdf5_orec(ndims, idims, 'ZBASE'   , rvars=zbase)

idims(1) = ngrids

call shdf5_orec(ndims, idims, 'CENTLAT', rvara=centlat)
call shdf5_orec(ndims, idims, 'CENTLON', rvara=centlon)
call shdf5_orec(ndims, idims, 'GRDLEN' , rvara=grdlen)
call shdf5_orec(ndims, idims, 'GRDWID' , rvara=grdwid)
call shdf5_orec(ndims, idims, 'GRDAXIS', rvara=grdaxis)

! Write the grid dimensions

idims(1) = 1

call shdf5_orec(ndims, idims, 'NZA'    , ivars=nza)
call shdf5_orec(ndims, idims, 'NMA'    , ivars=nma)
call shdf5_orec(ndims, idims, 'NUA'    , ivars=nua)
call shdf5_orec(ndims, idims, 'NWA'    , ivars=nwa)
call shdf5_orec(ndims, idims, 'NSW_MAX', ivars=nsw_max)
call shdf5_orec(ndims, idims, 'MRLS'   , ivars=mrls)
call shdf5_orec(ndims, idims, 'MNUDP'  , ivars=mnudp)

! Write grid structure variables

idims(1) = nza

call shdf5_orec(ndims, idims, 'ZM'  , rvara=zm)
call shdf5_orec(ndims, idims, 'ZT'  , rvara=zt)
call shdf5_orec(ndims, idims, 'DZM' , rvara=dzm)
call shdf5_orec(ndims, idims, 'DZT' , rvara=dzt)
call shdf5_orec(ndims, idims, 'DZIM', rvara=dzim)
call shdf5_orec(ndims, idims, 'DZIT', rvara=dzit)

idims(1) = nma

call shdf5_orec(ndims, idims, 'TOPM' , rvara=topm)
call shdf5_orec(ndims, idims, 'XEM'  , rvara=xem)
call shdf5_orec(ndims, idims, 'YEM'  , rvara=yem)
call shdf5_orec(ndims, idims, 'ZEM'  , rvara=zem)
call shdf5_orec(ndims, idims, 'GLATM', rvara=glatm)
call shdf5_orec(ndims, idims, 'GLONM', rvara=glonm)

idims(1) = nua

call shdf5_orec(ndims, idims, 'LPU' , ivara=lpu)
call shdf5_orec(ndims, idims, 'LCU' , ivara=lcu)
call shdf5_orec(ndims, idims, 'XEU' , rvara=xeu)
call shdf5_orec(ndims, idims, 'YEU' , rvara=yeu)
call shdf5_orec(ndims, idims, 'ZEU' , rvara=zeu)
call shdf5_orec(ndims, idims, 'DTU' , rvara=dtu)
call shdf5_orec(ndims, idims, 'DNU' , rvara=dnu)
call shdf5_orec(ndims, idims, 'DNIU', rvara=dniu)
call shdf5_orec(ndims, idims, 'UNX' , rvara=unx)
call shdf5_orec(ndims, idims, 'UNY' , rvara=uny)
call shdf5_orec(ndims, idims, 'UNZ' , rvara=unz)
call shdf5_orec(ndims, idims, 'UTX' , rvara=utx)
call shdf5_orec(ndims, idims, 'UTY' , rvara=uty)
call shdf5_orec(ndims, idims, 'UTZ' , rvara=utz)

idims(1) = nwa

call shdf5_orec(ndims, idims, 'LPW'  , ivara=lpw)
call shdf5_orec(ndims, idims, 'LSW'  , ivara=lsw)
call shdf5_orec(ndims, idims, 'XEW'  , rvara=xew)
call shdf5_orec(ndims, idims, 'YEW'  , rvara=yew)
call shdf5_orec(ndims, idims, 'ZEW'  , rvara=zew)
call shdf5_orec(ndims, idims, 'ARW0' , rvara=arw0)
call shdf5_orec(ndims, idims, 'GLATW', rvara=glatw)
call shdf5_orec(ndims, idims, 'GLONW', rvara=glonw)
call shdf5_orec(ndims, idims, 'WNX'  , rvara=wnx)
call shdf5_orec(ndims, idims, 'WNY'  , rvara=wny)
call shdf5_orec(ndims, idims, 'WNZ'  , rvara=wnz)

ndims = 2
idims(1) = nza
idims(2) = nua

call shdf5_orec(ndims, idims, 'ARU'  , rvara=aru)
call shdf5_orec(ndims, idims, 'VOLUI', rvara=volui)

idims(2) = nwa

call shdf5_orec(ndims, idims, 'ARW'  , rvara=arw)
call shdf5_orec(ndims, idims, 'VOLWI', rvara=volwi)
call shdf5_orec(ndims, idims, 'VOLT' , dvara=volt)
call shdf5_orec(ndims, idims, 'VOLTI', dvara=volti)

! Write ITAB_M SCALARS

ndims = 1
idims(1) = nma
idims(2) = 1

call shdf5_orec(ndims,idims,'itab_m%ntpn'   ,ivara=itab_m(:)%ntpn)
call shdf5_orec(ndims,idims,'itab_m%itopm'  ,ivara=itab_m(:)%itopm)
call shdf5_orec(ndims,idims,'itab_m%imglobe',ivara=itab_m(:)%imglobe)
call shdf5_orec(ndims,idims,'itab_m%arm'    ,rvara=itab_m(:)%arm)

! Write ITAB_M LOOP ARRAY using memory copy

allocate (lscr(nloops_m,nma))

ndims = 2
idims(1) = nloops_m
idims(2) = nma

do im = 1,nma
   lscr(1:nloops_m,im) = itab_m(im)%loop(1:nloops_m)
enddo

call shdf5_orec(ndims,idims,'itab_m%loop',lvara=lscr)

deallocate (lscr)

! Write ITAB_M INTEGER ARRAYS using memory copy

allocate (iscr(maxtpn,nma))

ndims = 2
idims(1) = maxtpn
idims(2) = nma

do im = 1,nma
   iscr(1:maxtpn,im) = itab_m(im)%iw(1:maxtpn)
enddo

call shdf5_orec(ndims,idims,'itab_m%iw',ivara=iscr)

do im = 1,nma
   iscr(1:maxtpn,im) = itab_m(im)%iu(1:maxtpn)
enddo

call shdf5_orec(ndims,idims,'itab_m%iu',ivara=iscr)

deallocate (iscr)

! Write ITAB_U SCALARS

ndims = 1
idims(1) = nua
idims(2) = 1

call shdf5_orec(ndims,idims,'itab_u%iup'    ,ivara=itab_u(:)%iup)

call shdf5_orec(ndims,idims,'itab_u%im1'    ,ivara=itab_u(:)%im1)
call shdf5_orec(ndims,idims,'itab_u%im2'    ,ivara=itab_u(:)%im2)

call shdf5_orec(ndims,idims,'itab_u%iu1'    ,ivara=itab_u(:)%iu1)
call shdf5_orec(ndims,idims,'itab_u%iu2'    ,ivara=itab_u(:)%iu2)
call shdf5_orec(ndims,idims,'itab_u%iu3'    ,ivara=itab_u(:)%iu3)
call shdf5_orec(ndims,idims,'itab_u%iu4'    ,ivara=itab_u(:)%iu4)
call shdf5_orec(ndims,idims,'itab_u%iu5'    ,ivara=itab_u(:)%iu5)
call shdf5_orec(ndims,idims,'itab_u%iu6'    ,ivara=itab_u(:)%iu6)
call shdf5_orec(ndims,idims,'itab_u%iu7'    ,ivara=itab_u(:)%iu7)
call shdf5_orec(ndims,idims,'itab_u%iu8'    ,ivara=itab_u(:)%iu8)
call shdf5_orec(ndims,idims,'itab_u%iu9'    ,ivara=itab_u(:)%iu9)
call shdf5_orec(ndims,idims,'itab_u%iu10'   ,ivara=itab_u(:)%iu10)
call shdf5_orec(ndims,idims,'itab_u%iu11'   ,ivara=itab_u(:)%iu11)
call shdf5_orec(ndims,idims,'itab_u%iu12'   ,ivara=itab_u(:)%iu12)

call shdf5_orec(ndims,idims,'itab_u%iw1'    ,ivara=itab_u(:)%iw1)
call shdf5_orec(ndims,idims,'itab_u%iw2'    ,ivara=itab_u(:)%iw2)
call shdf5_orec(ndims,idims,'itab_u%iw3'    ,ivara=itab_u(:)%iw3)
call shdf5_orec(ndims,idims,'itab_u%iw4'    ,ivara=itab_u(:)%iw4)
call shdf5_orec(ndims,idims,'itab_u%iw5'    ,ivara=itab_u(:)%iw5)
call shdf5_orec(ndims,idims,'itab_u%iw6'    ,ivara=itab_u(:)%iw6)

call shdf5_orec(ndims,idims,'itab_u%mrlu'   ,ivara=itab_u(:)%mrlu)
call shdf5_orec(ndims,idims,'itab_u%iuglobe',ivara=itab_u(:)%iuglobe)

call shdf5_orec(ndims,idims,'itab_u%diru1'  ,rvara=itab_u(:)%diru1)
call shdf5_orec(ndims,idims,'itab_u%diru2'  ,rvara=itab_u(:)%diru2)
call shdf5_orec(ndims,idims,'itab_u%diru3'  ,rvara=itab_u(:)%diru3)
call shdf5_orec(ndims,idims,'itab_u%diru4'  ,rvara=itab_u(:)%diru4)

call shdf5_orec(ndims,idims,'itab_u%fuu5'   ,rvara=itab_u(:)%fuu5)
call shdf5_orec(ndims,idims,'itab_u%fuu6'   ,rvara=itab_u(:)%fuu6)
call shdf5_orec(ndims,idims,'itab_u%fuu7'   ,rvara=itab_u(:)%fuu7)
call shdf5_orec(ndims,idims,'itab_u%fuu8'   ,rvara=itab_u(:)%fuu8)
call shdf5_orec(ndims,idims,'itab_u%fuu9'   ,rvara=itab_u(:)%fuu9)
call shdf5_orec(ndims,idims,'itab_u%fuu10'  ,rvara=itab_u(:)%fuu10)
call shdf5_orec(ndims,idims,'itab_u%fuu11'  ,rvara=itab_u(:)%fuu11)
call shdf5_orec(ndims,idims,'itab_u%fuu12'  ,rvara=itab_u(:)%fuu12)

call shdf5_orec(ndims,idims,'itab_u%fuw3'   ,rvara=itab_u(:)%fuw3)
call shdf5_orec(ndims,idims,'itab_u%fuw4'   ,rvara=itab_u(:)%fuw4)
call shdf5_orec(ndims,idims,'itab_u%fuw5'   ,rvara=itab_u(:)%fuw5)
call shdf5_orec(ndims,idims,'itab_u%fuw6'   ,rvara=itab_u(:)%fuw6)

call shdf5_orec(ndims,idims,'itab_u%tuu1'   ,rvara=itab_u(:)%tuu1)
call shdf5_orec(ndims,idims,'itab_u%tuu2'   ,rvara=itab_u(:)%tuu2)
call shdf5_orec(ndims,idims,'itab_u%tuu3'   ,rvara=itab_u(:)%tuu3)
call shdf5_orec(ndims,idims,'itab_u%tuu4'   ,rvara=itab_u(:)%tuu4)

call shdf5_orec(ndims,idims,'itab_u%pgc12'  ,rvara=itab_u(:)%pgc12)
call shdf5_orec(ndims,idims,'itab_u%pgc45'  ,rvara=itab_u(:)%pgc45)
call shdf5_orec(ndims,idims,'itab_u%pgc63'  ,rvara=itab_u(:)%pgc63)
call shdf5_orec(ndims,idims,'itab_u%pgc12b' ,rvara=itab_u(:)%pgc12b)
call shdf5_orec(ndims,idims,'itab_u%pgc45b' ,rvara=itab_u(:)%pgc45b)
call shdf5_orec(ndims,idims,'itab_u%pgc12c' ,rvara=itab_u(:)%pgc12c)
call shdf5_orec(ndims,idims,'itab_u%pgc63c' ,rvara=itab_u(:)%pgc63c)
call shdf5_orec(ndims,idims,'itab_u%pgc12d' ,rvara=itab_u(:)%pgc12d)

call shdf5_orec(ndims,idims,'itab_u%vxu1_u' ,rvara=itab_u(:)%vxu1_u)
call shdf5_orec(ndims,idims,'itab_u%vxu2_u' ,rvara=itab_u(:)%vxu2_u)
call shdf5_orec(ndims,idims,'itab_u%vxu3_u' ,rvara=itab_u(:)%vxu3_u)
call shdf5_orec(ndims,idims,'itab_u%vxu4_u' ,rvara=itab_u(:)%vxu4_u)

call shdf5_orec(ndims,idims,'itab_u%vxw1_u' ,rvara=itab_u(:)%vxw1_u)
call shdf5_orec(ndims,idims,'itab_u%vxw2_u' ,rvara=itab_u(:)%vxw2_u)

call shdf5_orec(ndims,idims,'itab_u%vyu1_u' ,rvara=itab_u(:)%vyu1_u)
call shdf5_orec(ndims,idims,'itab_u%vyu2_u' ,rvara=itab_u(:)%vyu2_u)
call shdf5_orec(ndims,idims,'itab_u%vyu3_u' ,rvara=itab_u(:)%vyu3_u)
call shdf5_orec(ndims,idims,'itab_u%vyu4_u' ,rvara=itab_u(:)%vyu4_u)

call shdf5_orec(ndims,idims,'itab_u%vyw1_u' ,rvara=itab_u(:)%vyw1_u)
call shdf5_orec(ndims,idims,'itab_u%vyw2_u' ,rvara=itab_u(:)%vyw2_u)

! Write ITAB_U LOOP ARRAY using memory copy

allocate (lscr(nloops_u,nua))

ndims = 2
idims(1) = nloops_u
idims(2) = nua

do iu = 1,nua
   lscr(1:nloops_u,iu) = itab_u(iu)%loop(1:nloops_u)
enddo

call shdf5_orec(ndims,idims,'itab_u%loop',lvara=lscr)

deallocate (lscr)

! Write ITAB_W SCALARS

ndims = 1
idims(1) = nwa
idims(2) = 1

call shdf5_orec(ndims,idims,'itab_w%iwp'      ,ivara=itab_w(:)%iwp)

call shdf5_orec(ndims,idims,'itab_w%im1'    ,ivara=itab_w(:)%im1)
call shdf5_orec(ndims,idims,'itab_w%im2'    ,ivara=itab_w(:)%im2)
call shdf5_orec(ndims,idims,'itab_w%im3'    ,ivara=itab_w(:)%im3)

call shdf5_orec(ndims,idims,'itab_w%iu1'    ,ivara=itab_w(:)%iu1)
call shdf5_orec(ndims,idims,'itab_w%iu2'    ,ivara=itab_w(:)%iu2)
call shdf5_orec(ndims,idims,'itab_w%iu3'    ,ivara=itab_w(:)%iu3)
call shdf5_orec(ndims,idims,'itab_w%iu4'    ,ivara=itab_w(:)%iu4)
call shdf5_orec(ndims,idims,'itab_w%iu5'    ,ivara=itab_w(:)%iu5)
call shdf5_orec(ndims,idims,'itab_w%iu6'    ,ivara=itab_w(:)%iu6)
call shdf5_orec(ndims,idims,'itab_w%iu7'    ,ivara=itab_w(:)%iu7)
call shdf5_orec(ndims,idims,'itab_w%iu8'    ,ivara=itab_w(:)%iu8)
call shdf5_orec(ndims,idims,'itab_w%iu9'    ,ivara=itab_w(:)%iu9)

call shdf5_orec(ndims,idims,'itab_w%iw1'    ,ivara=itab_w(:)%iw1)
call shdf5_orec(ndims,idims,'itab_w%iw2'    ,ivara=itab_w(:)%iw2)
call shdf5_orec(ndims,idims,'itab_w%iw3'    ,ivara=itab_w(:)%iw3)

call shdf5_orec(ndims,idims,'itab_w%iwglobe'  ,ivara=itab_w(:)%iwglobe)
call shdf5_orec(ndims,idims,'itab_w%mrlw'     ,ivara=itab_w(:)%mrlw)
call shdf5_orec(ndims,idims,'itab_w%mrlw_orig',ivara=itab_w(:)%mrlw_orig)
call shdf5_orec(ndims,idims,'itab_w%mrow'     ,ivara=itab_w(:)%mrow)
call shdf5_orec(ndims,idims,'itab_w%mrowh'    ,ivara=itab_w(:)%mrowh)

call shdf5_orec(ndims,idims,'itab_w%diru1'  ,rvara=itab_w(:)%diru1)
call shdf5_orec(ndims,idims,'itab_w%diru2'  ,rvara=itab_w(:)%diru2)
call shdf5_orec(ndims,idims,'itab_w%diru3'  ,rvara=itab_w(:)%diru3)

call shdf5_orec(ndims,idims,'itab_w%fwu4'   ,rvara=itab_w(:)%fwu4)
call shdf5_orec(ndims,idims,'itab_w%fwu5'   ,rvara=itab_w(:)%fwu5)
call shdf5_orec(ndims,idims,'itab_w%fwu6'   ,rvara=itab_w(:)%fwu6)
call shdf5_orec(ndims,idims,'itab_w%fwu7'   ,rvara=itab_w(:)%fwu7)
call shdf5_orec(ndims,idims,'itab_w%fwu8'   ,rvara=itab_w(:)%fwu8)
call shdf5_orec(ndims,idims,'itab_w%fwu9'   ,rvara=itab_w(:)%fwu9)

call shdf5_orec(ndims,idims,'itab_w%fww1'   ,rvara=itab_w(:)%fww1)
call shdf5_orec(ndims,idims,'itab_w%fww2'   ,rvara=itab_w(:)%fww2)
call shdf5_orec(ndims,idims,'itab_w%fww3'   ,rvara=itab_w(:)%fww3)

call shdf5_orec(ndims,idims,'itab_w%vxu1' ,rvara=itab_w(:)%vxu1)
call shdf5_orec(ndims,idims,'itab_w%vxu2' ,rvara=itab_w(:)%vxu2)
call shdf5_orec(ndims,idims,'itab_w%vxu3' ,rvara=itab_w(:)%vxu3)

call shdf5_orec(ndims,idims,'itab_w%vxw' ,rvara=itab_w(:)%vxw)

call shdf5_orec(ndims,idims,'itab_w%vyu1' ,rvara=itab_w(:)%vyu1)
call shdf5_orec(ndims,idims,'itab_w%vyu2' ,rvara=itab_w(:)%vyu2)
call shdf5_orec(ndims,idims,'itab_w%vyu3' ,rvara=itab_w(:)%vyu3)

call shdf5_orec(ndims,idims,'itab_w%vyw' ,rvara=itab_w(:)%vyw)

call shdf5_orec(ndims,idims,'itab_w%vzu1' ,rvara=itab_w(:)%vzu1)
call shdf5_orec(ndims,idims,'itab_w%vzu2' ,rvara=itab_w(:)%vzu2)
call shdf5_orec(ndims,idims,'itab_w%vzu3' ,rvara=itab_w(:)%vzu3)

call shdf5_orec(ndims,idims,'itab_w%vzw' ,rvara=itab_w(:)%vzw)

call shdf5_orec(ndims,idims,'itab_w%vxu1_w' ,rvara=itab_w(:)%vxu1_w)
call shdf5_orec(ndims,idims,'itab_w%vxu2_w' ,rvara=itab_w(:)%vxu2_w)
call shdf5_orec(ndims,idims,'itab_w%vxu3_w' ,rvara=itab_w(:)%vxu3_w)

call shdf5_orec(ndims,idims,'itab_w%vyu1_w' ,rvara=itab_w(:)%vyu1_w)
call shdf5_orec(ndims,idims,'itab_w%vyu2_w' ,rvara=itab_w(:)%vyu2_w)
call shdf5_orec(ndims,idims,'itab_w%vyu3_w' ,rvara=itab_w(:)%vyu3_w)

! Write ITAB_W LOOP ARRAY using memory copy

allocate (lscr(nloops_w,nwa))

ndims = 2
idims(1) = nloops_w
idims(2) = nwa

do iw = 1,nwa
   lscr(1:nloops_w,iw) = itab_w(iw)%loop(1:nloops_w)
enddo

call shdf5_orec(ndims,idims,'itab_w%loop',lvara=lscr)

deallocate (lscr)

! Write ITAB_W INTEGER AND REAL ARRAYS using memory copy

allocate (iscr(3,nwa))
allocate (rscr(3,nwa))

ndims = 2
idims(1) = 3
idims(2) = nwa

do iw = 1,nwa
   iscr(1:3,iw) = itab_w(iw)%inudp(1:3)
   rscr(1:3,iw) = itab_w(iw)%fnudp(1:3)
enddo

call shdf5_orec(ndims,idims,'itab_w%inudp',ivara=iscr)
call shdf5_orec(ndims,idims,'itab_w%fnudp',rvara=rscr)

deallocate (iscr,rscr)

! Check whether LAND/SEA models are used

if (isfcl == 1) then

! Write SEA VALUES

   ndims = 1
   idims(1) = 1
   idims(2) = 1

   call shdf5_orec(ndims, idims, 'NSEAFLUX' ,ivars=nseaflux)
   call shdf5_orec(ndims, idims, 'NTRAPS'   ,ivars=ntraps)

   idims(1) = nseaflux

   call shdf5_orec(ndims,idims,'seaflux%isfglobe',ivara=seaflux(:)%isfglobe)
   call shdf5_orec(ndims,idims,'seaflux%iw'      ,ivara=seaflux(:)%iw)
   call shdf5_orec(ndims,idims,'seaflux%kw'      ,ivara=seaflux(:)%kw)
   call shdf5_orec(ndims,idims,'seaflux%iws'     ,ivara=seaflux(:)%iws)
   call shdf5_orec(ndims,idims,'seaflux%jtrap'   ,ivara=seaflux(:)%jtrap)
   call shdf5_orec(ndims,idims,'seaflux%itrap'   ,ivara=seaflux(:)%itrap)
   call shdf5_orec(ndims,idims,'seaflux%area'    ,rvara=seaflux(:)%area)
   call shdf5_orec(ndims,idims,'seaflux%xef'     ,rvara=seaflux(:)%xef)
   call shdf5_orec(ndims,idims,'seaflux%yef'     ,rvara=seaflux(:)%yef)
   call shdf5_orec(ndims,idims,'seaflux%zef'     ,rvara=seaflux(:)%zef)
   call shdf5_orec(ndims,idims,'seaflux%arf_atm' ,rvara=seaflux(:)%arf_atm)
   call shdf5_orec(ndims,idims,'seaflux%arf_sea' ,rvara=seaflux(:)%arf_sea)

   if (ntraps > 0) then

      ndims = 2
      idims(1) = 4
      idims(2) = ntraps

      call shdf5_orec(ndims,idims,'xemstrap' ,rvara=xemstrap)
      call shdf5_orec(ndims,idims,'yemstrap' ,rvara=yemstrap)
      call shdf5_orec(ndims,idims,'zemstrap' ,rvara=zemstrap)

   endif
   
! Write LAND VALUES

   ndims = 1
   idims(1) = 1
   idims(2) = 1

   call shdf5_orec(ndims, idims, 'NLANDFLUX' ,ivars=nlandflux)
   call shdf5_orec(ndims, idims, 'NTRAPL'    ,ivars=ntrapl)

   idims(1) = nlandflux

   call shdf5_orec(ndims,idims,'landflux%ilfglobe',ivara=landflux(:)%ilfglobe)
   call shdf5_orec(ndims,idims,'landflux%iw'      ,ivara=landflux(:)%iw)
   call shdf5_orec(ndims,idims,'landflux%kw'      ,ivara=landflux(:)%kw)
   call shdf5_orec(ndims,idims,'landflux%iwl'     ,ivara=landflux(:)%iwl)
   call shdf5_orec(ndims,idims,'landflux%jtrap'   ,ivara=landflux(:)%jtrap)
   call shdf5_orec(ndims,idims,'landflux%itrap'   ,ivara=landflux(:)%itrap)
   call shdf5_orec(ndims,idims,'landflux%area'    ,rvara=landflux(:)%area)
   call shdf5_orec(ndims,idims,'landflux%xef'     ,rvara=landflux(:)%xef)
   call shdf5_orec(ndims,idims,'landflux%yef'     ,rvara=landflux(:)%yef)
   call shdf5_orec(ndims,idims,'landflux%zef'     ,rvara=landflux(:)%zef)
   call shdf5_orec(ndims,idims,'landflux%arf_atm' ,rvara=landflux(:)%arf_atm)
   call shdf5_orec(ndims,idims,'landflux%arf_land',rvara=landflux(:)%arf_land)

   if (ntrapl > 0) then

      ndims = 2
      idims(1) = 4
      idims(2) = ntrapl

      call shdf5_orec(ndims,idims,'xemltrap' ,rvara=xemltrap)
      call shdf5_orec(ndims,idims,'yemltrap' ,rvara=yemltrap)
      call shdf5_orec(ndims,idims,'zemltrap' ,rvara=zemltrap)

   endif

endif

! Close GRIDFILE

call shdf5_close()

return
end subroutine gridfile_write

!===============================================================================

subroutine gridfile_read()

use misc_coms,  only: io6, ngrids, gridfile, mdomain, nzp, nxp,  &
                      itopoflg,  &
                      deltax, deltaz, dzmax, dzrat, zbase,  &
                      centlat, centlon, grdlen, grdwid, grdaxis
use mem_ijtabs, only: nloops_m, nloops_u, nloops_w, maxtpn, mrls,  &
                      itab_m, itab_u, itab_w,  &
                      alloc_itabs
use mem_grid,   only: nza, nma, nua, nwa, mza, mma, mua, mwa, nsw_max,  &
                      zm, zt, dzm, dzt, dzim, dzit,  &
                      lpu, lcu, lpw, lsw,  &
                      topm, xem, yem, zem, xeu, yeu, zeu, xew, yew, zew,  &
                      unx, uny, unz, utx, uty, utz, wnx, wny, wnz,  &
                      dtu, dnu, dniu, arw0, glatw, glonw, glatm, glonm,  &
                      aru, volui, arw, volwi, volt, volti,  &
                      alloc_gridz, alloc_xyzem, alloc_grid1, alloc_grid2
use leaf_coms,  only: isfcl
use mem_sflux,  only: nseaflux, nlandflux, mseaflux, mlandflux,  &
                      ntraps, ntrapl, mtraps, mtrapl,  &
                      seaflux, landflux,  &
                      xemstrap, yemstrap, zemstrap, xemltrap, yemltrap, zemltrap

use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
use mem_nudge,  only: mnudp
use mem_para,   only: myrank

! This subroutine checks for the existence of a gridfile, and if it exists, 
! also checks for agreement of grid configuration between the file and the 
! current model run.  If the file does not exist or does not match grid
! configuration, the run is stopped.

implicit none

integer :: im, iu, iw

integer :: ierr

integer :: ngr
integer :: ndims, idims(2)

integer :: ngrids0, mdomain0, nxp0, nzp0, itopoflg0, isfcl0

real :: deltax0, deltaz0, dzrat0, dzmax0, zbase0

character(len=128) :: flnm
character(len=10) :: number

logical  :: there, exans

real :: centlat0(ngrids)
real :: centlon0(ngrids)
real :: grdlen0 (ngrids)
real :: grdwid0 (ngrids)
real :: grdaxis0(ngrids)

! Scratch arrays for copying intput

integer, allocatable :: iscr(:,:)
real,    allocatable :: rscr(:,:)
logical, allocatable :: lscr(:,:)

real, external :: walltime
real :: wtime_start_gridfileread, inicio

wtime_start_gridfileread = walltime(0.)


! Check if grid file exists
inicio = walltime(wtime_start_gridfileread)

inquire(file=gridfile, exist=exans)


if (exans) then

! Grid file exists.  Open, read, and close file.

   write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   write(io6,*) 'Opening grid file ', trim(gridfile)
   write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

   call shdf5_open(trim(gridfile),'R')

   write(io6, *) '              ==T== Abrir o gridfile: ',(walltime(wtime_start_gridfileread)-inicio)

   inicio = walltime(wtime_start_gridfileread)

! Read the grid information that exists in namelist

   ndims = 1
   idims(1) = 1
   idims(2) = 1

   call shdf5_irec(ndims, idims, 'NZP'     , ivars=nzp0)
   call shdf5_irec(ndims, idims, 'NXP'     , ivars=nxp0)
   call shdf5_irec(ndims, idims, 'MDOMAIN' , ivars=mdomain0)
   call shdf5_irec(ndims, idims, 'NGRIDS'  , ivars=ngrids0)
   call shdf5_irec(ndims, idims, 'ISFCL'   , ivars=isfcl0)
   call shdf5_irec(ndims, idims, 'ITOPOFLG', ivars=itopoflg0)
   call shdf5_irec(ndims, idims, 'DELTAX'  , rvars=deltax0)
   call shdf5_irec(ndims, idims, 'DELTAZ'  , rvars=deltaz0)
   call shdf5_irec(ndims, idims, 'DZRAT'   , rvars=dzrat0)
   call shdf5_irec(ndims, idims, 'DZMAX'   , rvars=dzmax0)
   call shdf5_irec(ndims, idims, 'ZBASE'   , rvars=zbase0)

   idims(1) = ngrids0

   call shdf5_irec(ndims, idims, 'CENTLAT', rvara=centlat0)
   call shdf5_irec(ndims, idims, 'CENTLON', rvara=centlon0)
   call shdf5_irec(ndims, idims, 'GRDLEN' , rvara=grdlen0)
   call shdf5_irec(ndims, idims, 'GRDWID' , rvara=grdwid0)
   call shdf5_irec(ndims, idims, 'GRDAXIS', rvara=grdaxis0)

! Check equality between grid file information and namelist variables

   ierr = 0

   if (nzp0      /= nzp     ) ierr = 1 
   if (nxp0      /= nxp     ) ierr = 1 
   if (mdomain0  /= mdomain ) ierr = 1 
   if (ngrids0   /= ngrids  ) ierr = 1 
   if (isfcl0    /= isfcl   ) ierr = 1 
   if (itopoflg0 /= itopoflg) ierr = 1 

   if (abs(deltax0 - deltax) > 1.e-3) ierr = 1 
   if (abs(deltaz0 - deltaz) > 1.e-3) ierr = 1 
   if (abs(dzrat0  - dzrat ) > 1.e-3) ierr = 1 
   if (abs(dzmax0  - dzmax ) > 1.e-3) ierr = 1 
   if (abs(zbase0  - zbase ) > 1.e-3) ierr = 1 
   
   do ngr = 1,ngrids0
      if (abs(centlat0(ngr) - centlat(ngr)) > 1.e-3) ierr = 1
      if (abs(centlon0(ngr) - centlon(ngr)) > 1.e-3) ierr = 1
      if (abs(grdlen0 (ngr) - grdlen (ngr)) > 1.e1 ) ierr = 1
      if (abs(grdwid0 (ngr) - grdwid (ngr)) > 1.e1 ) ierr = 1
      if (abs(grdaxis0(ngr) - grdaxis(ngr)) > 1.e-3) ierr = 1
   enddo

   if (ierr == 1) then

      write(io6,*) 'GRIDFILE mismatch with OLAMIN namelist: Stopping model run'
      write(io6,*) 'Values: gridfile, namelist'
      write(io6,*) '-----------------------------------------------'
      write(io6,*)              'nzp:      ',nzp0     ,nzp
      write(io6,*)              'nxp:      ',nxp0     ,nxp
      write(io6,*)              'mdomain:  ',mdomain0 ,mdomain
      write(io6,*)              'ngrids:   ',ngrids0  ,ngrids
      write(io6,*)              'isfcl:    ',isfcl0   ,isfcl
      write(io6,*)              'itopoflg: ',itopoflg0,itopoflg
      write(io6,*)              'deltax:   ',deltax0  ,deltax
      write(io6,*)              'deltaz:   ',deltaz0  ,deltaz
      write(io6,*)              'dzrat:    ',dzrat0   ,dzrat
      write(io6,*)              'dzmax:    ',dzmax0   ,dzmax
      write(io6,*)              'zbase:    ',zbase0   ,zbase
      write(io6,*) ' '
      write(io6, '(a,20f10.3)') 'centlat0: ',centlat0(1:ngrids)
      write(io6, '(a,20f10.3)') 'centlat:  ',centlat (1:ngrids)
      write(io6,*) ' '
      write(io6, '(a,20f10.3)') 'centlon0: ',centlon0(1:ngrids)
      write(io6, '(a,20f10.3)') 'centlon:  ',centlon (1:ngrids)
      write(io6,*) ' '
      write(io6, '(a,20f12.1)') 'grdlen0:  ',grdlen0 (1:ngrids)
      write(io6, '(a,20f12.1)') 'grdlen:   ',grdlen  (1:ngrids)
      write(io6,*) ' '
      write(io6, '(a,20f12.1)') 'grdwid0:  ',grdwid0 (1:ngrids)
      write(io6, '(a,20f12.1)') 'grdwid:   ',grdwid  (1:ngrids)
      write(io6,*) ' '
      write(io6, '(a,20f12.1)') 'grdaxis0: ',grdaxis0(1:ngrids)
      write(io6, '(a,20f12.1)') 'grdaxis:  ',grdaxis (1:ngrids)
      write(io6,*) '-----------------------------------------------'

      stop 'stop - gridfile mismatch'
   
   endif

! Read the grid dimensions

   call shdf5_irec(ndims, idims, 'NZA'    , ivars=nza)
   call shdf5_irec(ndims, idims, 'NMA'    , ivars=nma)
   call shdf5_irec(ndims, idims, 'NUA'    , ivars=nua)
   call shdf5_irec(ndims, idims, 'NWA'    , ivars=nwa)
   call shdf5_irec(ndims, idims, 'NSW_MAX', ivars=nsw_max)
   call shdf5_irec(ndims, idims, 'MRLS'   , ivars=mrls)
   call shdf5_irec(ndims, idims, 'MNUDP'  , ivars=mnudp)

! Copy grid dimensions

   mza = nza
   mma = nma
   mua = nua
   mwa = nwa

! Allocate and read grid structure variables

   call alloc_gridz()
   call alloc_itabs(nma,nua,nwa)
   call alloc_xyzem()
   call alloc_grid1()
   call alloc_grid2()

   idims(1) = nza

   call shdf5_irec(ndims, idims, 'ZM'  , rvara=zm)
   call shdf5_irec(ndims, idims, 'ZT'  , rvara=zt)
   call shdf5_irec(ndims, idims, 'DZM' , rvara=dzm)
   call shdf5_irec(ndims, idims, 'DZT' , rvara=dzt)
   call shdf5_irec(ndims, idims, 'DZIM', rvara=dzim)
   call shdf5_irec(ndims, idims, 'DZIT', rvara=dzit)

   idims(1) = nma

   call shdf5_irec(ndims, idims, 'TOPM' , rvara=topm)
   call shdf5_irec(ndims, idims, 'XEM'  , rvara=xem)
   call shdf5_irec(ndims, idims, 'YEM'  , rvara=yem)
   call shdf5_irec(ndims, idims, 'ZEM'  , rvara=zem)
   call shdf5_irec(ndims, idims, 'GLATM', rvara=glatm)
   call shdf5_irec(ndims, idims, 'GLONM', rvara=glonm)

   idims(1) = nua

   call shdf5_irec(ndims, idims, 'LPU' , ivara=lpu)
   call shdf5_irec(ndims, idims, 'LCU' , ivara=lcu)
   call shdf5_irec(ndims, idims, 'XEU' , rvara=xeu)
   call shdf5_irec(ndims, idims, 'YEU' , rvara=yeu)
   call shdf5_irec(ndims, idims, 'ZEU' , rvara=zeu)
   call shdf5_irec(ndims, idims, 'DTU' , rvara=dtu)
   call shdf5_irec(ndims, idims, 'DNU' , rvara=dnu)
   call shdf5_irec(ndims, idims, 'DNIU', rvara=dniu)
   call shdf5_irec(ndims, idims, 'UNX' , rvara=unx)
   call shdf5_irec(ndims, idims, 'UNY' , rvara=uny)
   call shdf5_irec(ndims, idims, 'UNZ' , rvara=unz)
   call shdf5_irec(ndims, idims, 'UTX' , rvara=utx)
   call shdf5_irec(ndims, idims, 'UTY' , rvara=uty)
   call shdf5_irec(ndims, idims, 'UTZ' , rvara=utz)

   idims(1) = nwa

   call shdf5_irec(ndims, idims, 'LPW'  , ivara=lpw)
   call shdf5_irec(ndims, idims, 'LSW'  , ivara=lsw)
   call shdf5_irec(ndims, idims, 'XEW'  , rvara=xew)
   call shdf5_irec(ndims, idims, 'YEW'  , rvara=yew)
   call shdf5_irec(ndims, idims, 'ZEW'  , rvara=zew)
   call shdf5_irec(ndims, idims, 'ARW0' , rvara=arw0)
   call shdf5_irec(ndims, idims, 'GLATW', rvara=glatw)
   call shdf5_irec(ndims, idims, 'GLONW', rvara=glonw)
   call shdf5_irec(ndims, idims, 'WNX'  , rvara=wnx)
   call shdf5_irec(ndims, idims, 'WNY'  , rvara=wny)
   call shdf5_irec(ndims, idims, 'WNZ'  , rvara=wnz)

   ndims = 2
   idims(1) = nza
   idims(2) = nua

   call shdf5_irec(ndims, idims, 'ARU'  , rvara=aru)
   call shdf5_irec(ndims, idims, 'VOLUI', rvara=volui)
   
   idims(2) = nwa

   call shdf5_irec(ndims, idims, 'ARW'  , rvara=arw)
   call shdf5_irec(ndims, idims, 'VOLWI', rvara=volwi)
   call shdf5_irec(ndims, idims, 'VOLT' , dvara=volt)
   call shdf5_irec(ndims, idims, 'VOLTI', dvara=volti)
   
! Read ITAB_M SCALARS

   ndims = 1
   idims(1) = nma
   idims(2) = 1

   call shdf5_irec(ndims,idims,'itab_m%ntpn'   ,ivara=itab_m(:)%ntpn)
   call shdf5_irec(ndims,idims,'itab_m%itopm'  ,ivara=itab_m(:)%itopm)
   call shdf5_irec(ndims,idims,'itab_m%imglobe',ivara=itab_m(:)%imglobe)
   call shdf5_irec(ndims,idims,'itab_m%arm'    ,rvara=itab_m(:)%arm)

! Read ITAB_M LOOP ARRAY using memory copy

   allocate (lscr(nloops_m,nma))

   ndims = 2
   idims(1) = nloops_m
   idims(2) = nma

   call shdf5_irec(ndims,idims,'itab_m%loop',lvara=lscr)

   do im = 1,nma
      itab_m(im)%loop(1:nloops_m) = lscr(1:nloops_m,im)
   enddo

   deallocate (lscr)

! Read ITAB_M INTEGER ARRAYS using memory copy

   allocate (iscr(maxtpn,nma))

   ndims = 2
   idims(1) = maxtpn
   idims(2) = nma

   call shdf5_irec(ndims,idims,'itab_m%iw',ivara=iscr)

   do im = 1,nma
      itab_m(im)%iw(1:maxtpn) = iscr(1:maxtpn,im)
   enddo

   call shdf5_irec(ndims,idims,'itab_m%iu',ivara=iscr)

   do im = 1,nma
      itab_m(im)%iu(1:maxtpn) = iscr(1:maxtpn,im)
   enddo

   deallocate (iscr)

! Read ITAB_U SCALARS

   ndims = 1
   idims(1) = nua
   idims(2) = 1

   call shdf5_irec(ndims,idims,'itab_u%iup'    ,ivara=itab_u(:)%iup)

   call shdf5_irec(ndims,idims,'itab_u%im1'    ,ivara=itab_u(:)%im1)
   call shdf5_irec(ndims,idims,'itab_u%im2'    ,ivara=itab_u(:)%im2)

   call shdf5_irec(ndims,idims,'itab_u%iu1'    ,ivara=itab_u(:)%iu1)
   call shdf5_irec(ndims,idims,'itab_u%iu2'    ,ivara=itab_u(:)%iu2)
   call shdf5_irec(ndims,idims,'itab_u%iu3'    ,ivara=itab_u(:)%iu3)
   call shdf5_irec(ndims,idims,'itab_u%iu4'    ,ivara=itab_u(:)%iu4)
   call shdf5_irec(ndims,idims,'itab_u%iu5'    ,ivara=itab_u(:)%iu5)
   call shdf5_irec(ndims,idims,'itab_u%iu6'    ,ivara=itab_u(:)%iu6)
   call shdf5_irec(ndims,idims,'itab_u%iu7'    ,ivara=itab_u(:)%iu7)
   call shdf5_irec(ndims,idims,'itab_u%iu8'    ,ivara=itab_u(:)%iu8)
   call shdf5_irec(ndims,idims,'itab_u%iu9'    ,ivara=itab_u(:)%iu9)
   call shdf5_irec(ndims,idims,'itab_u%iu10'   ,ivara=itab_u(:)%iu10)
   call shdf5_irec(ndims,idims,'itab_u%iu11'   ,ivara=itab_u(:)%iu11)
   call shdf5_irec(ndims,idims,'itab_u%iu12'   ,ivara=itab_u(:)%iu12)

   call shdf5_irec(ndims,idims,'itab_u%iw1'    ,ivara=itab_u(:)%iw1)
   call shdf5_irec(ndims,idims,'itab_u%iw2'    ,ivara=itab_u(:)%iw2)
   call shdf5_irec(ndims,idims,'itab_u%iw3'    ,ivara=itab_u(:)%iw3)
   call shdf5_irec(ndims,idims,'itab_u%iw4'    ,ivara=itab_u(:)%iw4)
   call shdf5_irec(ndims,idims,'itab_u%iw5'    ,ivara=itab_u(:)%iw5)
   call shdf5_irec(ndims,idims,'itab_u%iw6'    ,ivara=itab_u(:)%iw6)

   call shdf5_irec(ndims,idims,'itab_u%iuglobe',ivara=itab_u(:)%iuglobe)
   call shdf5_irec(ndims,idims,'itab_u%mrlu'   ,ivara=itab_u(:)%mrlu)

   call shdf5_irec(ndims,idims,'itab_u%diru1'  ,rvara=itab_u(:)%diru1)
   call shdf5_irec(ndims,idims,'itab_u%diru2'  ,rvara=itab_u(:)%diru2)
   call shdf5_irec(ndims,idims,'itab_u%diru3'  ,rvara=itab_u(:)%diru3)
   call shdf5_irec(ndims,idims,'itab_u%diru4'  ,rvara=itab_u(:)%diru4)

   call shdf5_irec(ndims,idims,'itab_u%fuu5'   ,rvara=itab_u(:)%fuu5)
   call shdf5_irec(ndims,idims,'itab_u%fuu6'   ,rvara=itab_u(:)%fuu6)
   call shdf5_irec(ndims,idims,'itab_u%fuu7'   ,rvara=itab_u(:)%fuu7)
   call shdf5_irec(ndims,idims,'itab_u%fuu8'   ,rvara=itab_u(:)%fuu8)
   call shdf5_irec(ndims,idims,'itab_u%fuu9'   ,rvara=itab_u(:)%fuu9)
   call shdf5_irec(ndims,idims,'itab_u%fuu10'  ,rvara=itab_u(:)%fuu10)
   call shdf5_irec(ndims,idims,'itab_u%fuu11'  ,rvara=itab_u(:)%fuu11)
   call shdf5_irec(ndims,idims,'itab_u%fuu12'  ,rvara=itab_u(:)%fuu12)

   call shdf5_irec(ndims,idims,'itab_u%fuw3'   ,rvara=itab_u(:)%fuw3)
   call shdf5_irec(ndims,idims,'itab_u%fuw4'   ,rvara=itab_u(:)%fuw4)
   call shdf5_irec(ndims,idims,'itab_u%fuw5'   ,rvara=itab_u(:)%fuw5)
   call shdf5_irec(ndims,idims,'itab_u%fuw6'   ,rvara=itab_u(:)%fuw6)

   call shdf5_irec(ndims,idims,'itab_u%tuu1'   ,rvara=itab_u(:)%tuu1)
   call shdf5_irec(ndims,idims,'itab_u%tuu2'   ,rvara=itab_u(:)%tuu2)
   call shdf5_irec(ndims,idims,'itab_u%tuu3'   ,rvara=itab_u(:)%tuu3)
   call shdf5_irec(ndims,idims,'itab_u%tuu4'   ,rvara=itab_u(:)%tuu4)

   call shdf5_irec(ndims,idims,'itab_u%pgc12'  ,rvara=itab_u(:)%pgc12)
   call shdf5_irec(ndims,idims,'itab_u%pgc45'  ,rvara=itab_u(:)%pgc45)
   call shdf5_irec(ndims,idims,'itab_u%pgc63'  ,rvara=itab_u(:)%pgc63)
   call shdf5_irec(ndims,idims,'itab_u%pgc12b' ,rvara=itab_u(:)%pgc12b)
   call shdf5_irec(ndims,idims,'itab_u%pgc45b' ,rvara=itab_u(:)%pgc45b)
   call shdf5_irec(ndims,idims,'itab_u%pgc12c' ,rvara=itab_u(:)%pgc12c)
   call shdf5_irec(ndims,idims,'itab_u%pgc63c' ,rvara=itab_u(:)%pgc63c)
   call shdf5_irec(ndims,idims,'itab_u%pgc12d' ,rvara=itab_u(:)%pgc12d)

   call shdf5_irec(ndims,idims,'itab_u%vxu1_u' ,rvara=itab_u(:)%vxu1_u)
   call shdf5_irec(ndims,idims,'itab_u%vxu2_u' ,rvara=itab_u(:)%vxu2_u)
   call shdf5_irec(ndims,idims,'itab_u%vxu3_u' ,rvara=itab_u(:)%vxu3_u)
   call shdf5_irec(ndims,idims,'itab_u%vxu4_u' ,rvara=itab_u(:)%vxu4_u)

   call shdf5_irec(ndims,idims,'itab_u%vxw1_u' ,rvara=itab_u(:)%vxw1_u)
   call shdf5_irec(ndims,idims,'itab_u%vxw2_u' ,rvara=itab_u(:)%vxw2_u)

   call shdf5_irec(ndims,idims,'itab_u%vyu1_u' ,rvara=itab_u(:)%vyu1_u)
   call shdf5_irec(ndims,idims,'itab_u%vyu2_u' ,rvara=itab_u(:)%vyu2_u)
   call shdf5_irec(ndims,idims,'itab_u%vyu3_u' ,rvara=itab_u(:)%vyu3_u)
   call shdf5_irec(ndims,idims,'itab_u%vyu4_u' ,rvara=itab_u(:)%vyu4_u)

   call shdf5_irec(ndims,idims,'itab_u%vyw1_u' ,rvara=itab_u(:)%vyw1_u)
   call shdf5_irec(ndims,idims,'itab_u%vyw2_u' ,rvara=itab_u(:)%vyw2_u)

! Read ITAB_U LOOP ARRAY using memory copy

   allocate (lscr(nloops_u,nua))

   ndims = 2
   idims(1) = nloops_u
   idims(2) = nua

   call shdf5_irec(ndims,idims,'itab_u%loop',lvara=lscr)

   do iu = 1,nua
      itab_u(iu)%loop(1:nloops_u) = lscr(1:nloops_u,iu)
   enddo

   deallocate (lscr)

! Read ITAB_W SCALARS

   ndims = 1
   idims(1) = nwa
   idims(2) = 1

   call shdf5_irec(ndims,idims,'itab_w%iwp'      ,ivara=itab_w(:)%iwp)

   call shdf5_irec(ndims,idims,'itab_w%im1'    ,ivara=itab_w(:)%im1)
   call shdf5_irec(ndims,idims,'itab_w%im2'    ,ivara=itab_w(:)%im2)
   call shdf5_irec(ndims,idims,'itab_w%im3'    ,ivara=itab_w(:)%im3)

   call shdf5_irec(ndims,idims,'itab_w%iu1'    ,ivara=itab_w(:)%iu1)
   call shdf5_irec(ndims,idims,'itab_w%iu2'    ,ivara=itab_w(:)%iu2)
   call shdf5_irec(ndims,idims,'itab_w%iu3'    ,ivara=itab_w(:)%iu3)
   call shdf5_irec(ndims,idims,'itab_w%iu4'    ,ivara=itab_w(:)%iu4)
   call shdf5_irec(ndims,idims,'itab_w%iu5'    ,ivara=itab_w(:)%iu5)
   call shdf5_irec(ndims,idims,'itab_w%iu6'    ,ivara=itab_w(:)%iu6)
   call shdf5_irec(ndims,idims,'itab_w%iu7'    ,ivara=itab_w(:)%iu7)
   call shdf5_irec(ndims,idims,'itab_w%iu8'    ,ivara=itab_w(:)%iu8)
   call shdf5_irec(ndims,idims,'itab_w%iu9'    ,ivara=itab_w(:)%iu9)

   call shdf5_irec(ndims,idims,'itab_w%iw1'    ,ivara=itab_w(:)%iw1)
   call shdf5_irec(ndims,idims,'itab_w%iw2'    ,ivara=itab_w(:)%iw2)
   call shdf5_irec(ndims,idims,'itab_w%iw3'    ,ivara=itab_w(:)%iw3)

   call shdf5_irec(ndims,idims,'itab_w%iwglobe'  ,ivara=itab_w(:)%iwglobe)
   call shdf5_irec(ndims,idims,'itab_w%mrlw'     ,ivara=itab_w(:)%mrlw)
   call shdf5_irec(ndims,idims,'itab_w%mrlw_orig',ivara=itab_w(:)%mrlw_orig)
   call shdf5_irec(ndims,idims,'itab_w%mrow'     ,ivara=itab_w(:)%mrow)
   call shdf5_irec(ndims,idims,'itab_w%mrowh'    ,ivara=itab_w(:)%mrowh)

   call shdf5_irec(ndims,idims,'itab_w%diru1'  ,rvara=itab_w(:)%diru1)
   call shdf5_irec(ndims,idims,'itab_w%diru2'  ,rvara=itab_w(:)%diru2)
   call shdf5_irec(ndims,idims,'itab_w%diru3'  ,rvara=itab_w(:)%diru3)

   call shdf5_irec(ndims,idims,'itab_w%fwu4'   ,rvara=itab_w(:)%fwu4)
   call shdf5_irec(ndims,idims,'itab_w%fwu5'   ,rvara=itab_w(:)%fwu5)
   call shdf5_irec(ndims,idims,'itab_w%fwu6'   ,rvara=itab_w(:)%fwu6)
   call shdf5_irec(ndims,idims,'itab_w%fwu7'   ,rvara=itab_w(:)%fwu7)
   call shdf5_irec(ndims,idims,'itab_w%fwu8'   ,rvara=itab_w(:)%fwu8)
   call shdf5_irec(ndims,idims,'itab_w%fwu9'   ,rvara=itab_w(:)%fwu9)

   call shdf5_irec(ndims,idims,'itab_w%fww1'   ,rvara=itab_w(:)%fww1)
   call shdf5_irec(ndims,idims,'itab_w%fww2'   ,rvara=itab_w(:)%fww2)
   call shdf5_irec(ndims,idims,'itab_w%fww3'   ,rvara=itab_w(:)%fww3)

   call shdf5_irec(ndims,idims,'itab_w%vxu1' ,rvara=itab_w(:)%vxu1)
   call shdf5_irec(ndims,idims,'itab_w%vxu2' ,rvara=itab_w(:)%vxu2)
   call shdf5_irec(ndims,idims,'itab_w%vxu3' ,rvara=itab_w(:)%vxu3)

   call shdf5_irec(ndims,idims,'itab_w%vxw' ,rvara=itab_w(:)%vxw)

   call shdf5_irec(ndims,idims,'itab_w%vyu1' ,rvara=itab_w(:)%vyu1)
   call shdf5_irec(ndims,idims,'itab_w%vyu2' ,rvara=itab_w(:)%vyu2)
   call shdf5_irec(ndims,idims,'itab_w%vyu3' ,rvara=itab_w(:)%vyu3)

   call shdf5_irec(ndims,idims,'itab_w%vyw' ,rvara=itab_w(:)%vyw)

   call shdf5_irec(ndims,idims,'itab_w%vzu1' ,rvara=itab_w(:)%vzu1)
   call shdf5_irec(ndims,idims,'itab_w%vzu2' ,rvara=itab_w(:)%vzu2)
   call shdf5_irec(ndims,idims,'itab_w%vzu3' ,rvara=itab_w(:)%vzu3)

   call shdf5_irec(ndims,idims,'itab_w%vzw' ,rvara=itab_w(:)%vzw)

   call shdf5_irec(ndims,idims,'itab_w%vxu1_w' ,rvara=itab_w(:)%vxu1_w)
   call shdf5_irec(ndims,idims,'itab_w%vxu2_w' ,rvara=itab_w(:)%vxu2_w)
   call shdf5_irec(ndims,idims,'itab_w%vxu3_w' ,rvara=itab_w(:)%vxu3_w)

   call shdf5_irec(ndims,idims,'itab_w%vyu1_w' ,rvara=itab_w(:)%vyu1_w)
   call shdf5_irec(ndims,idims,'itab_w%vyu2_w' ,rvara=itab_w(:)%vyu2_w)
   call shdf5_irec(ndims,idims,'itab_w%vyu3_w' ,rvara=itab_w(:)%vyu3_w)

! Read itab_w LOOP ARRAY using memory copy

   allocate (lscr(nloops_w,nwa))

   ndims = 2
   idims(1) = nloops_w
   idims(2) = nwa

   call shdf5_irec(ndims,idims,'itab_w%loop',lvara=lscr)

   do iw = 1,nwa
      itab_w(iw)%loop(1:nloops_w) = lscr(1:nloops_w,iw)
   enddo

   deallocate (lscr)

! Write ITAB_W INTEGER AND REAL ARRAYS using memory copy

   allocate (iscr(3,nwa))
   allocate (rscr(3,nwa))

   ndims = 2
   idims(1) = 3
   idims(2) = nwa

   call shdf5_irec(ndims,idims,'itab_w%inudp',ivara=iscr)
   call shdf5_irec(ndims,idims,'itab_w%fnudp',rvara=rscr)

   do iw = 1,nwa
      itab_w(iw)%inudp(1:3) = iscr(1:3,iw)
      itab_w(iw)%fnudp(1:3) = rscr(1:3,iw)
   enddo

   deallocate (iscr,rscr)

! Check whether LAND/SEA models are used

   if (isfcl == 1) then

! Read SEAFLUX VALUES

      ndims = 1
      idims(1) = 1
      idims(2) = 1

      call shdf5_irec(ndims, idims, 'NSEAFLUX' ,ivars=nseaflux)
      call shdf5_irec(ndims, idims, 'NTRAPS'   ,ivars=ntraps)

      mseaflux = nseaflux
      mtraps = ntraps

      allocate (seaflux(nseaflux))

      allocate (xemstrap(4,ntraps))
      allocate (yemstrap(4,ntraps))
      allocate (zemstrap(4,ntraps))


      idims(1) = nseaflux

      call shdf5_irec(ndims,idims,'seaflux%isfglobe',ivara=seaflux(:)%isfglobe)
      call shdf5_irec(ndims,idims,'seaflux%iw'      ,ivara=seaflux(:)%iw)
      call shdf5_irec(ndims,idims,'seaflux%kw'      ,ivara=seaflux(:)%kw)
      call shdf5_irec(ndims,idims,'seaflux%iws'     ,ivara=seaflux(:)%iws)
      call shdf5_irec(ndims,idims,'seaflux%jtrap'   ,ivara=seaflux(:)%jtrap)
      call shdf5_irec(ndims,idims,'seaflux%itrap'   ,ivara=seaflux(:)%itrap)
      call shdf5_irec(ndims,idims,'seaflux%area'    ,rvara=seaflux(:)%area)
      call shdf5_irec(ndims,idims,'seaflux%xef'     ,rvara=seaflux(:)%xef)
      call shdf5_irec(ndims,idims,'seaflux%yef'     ,rvara=seaflux(:)%yef)
      call shdf5_irec(ndims,idims,'seaflux%zef'     ,rvara=seaflux(:)%zef)
      call shdf5_irec(ndims,idims,'seaflux%arf_atm' ,rvara=seaflux(:)%arf_atm)
      call shdf5_irec(ndims,idims,'seaflux%arf_sea' ,rvara=seaflux(:)%arf_sea)

      if (ntraps > 0) then

         ndims = 2
         idims(1) = 4
         idims(2) = ntraps

         call shdf5_irec(ndims,idims,'xemstrap' ,rvara=xemstrap)
         call shdf5_irec(ndims,idims,'yemstrap' ,rvara=yemstrap)
         call shdf5_irec(ndims,idims,'zemstrap' ,rvara=zemstrap)

      endif

! Read LANDFLUX VALUES

      ndims = 1
      idims(1) = 1
      idims(2) = 1

      call shdf5_irec(ndims, idims, 'NLANDFLUX' ,ivars=nlandflux)
      call shdf5_irec(ndims, idims, 'NTRAPL'    ,ivars=ntrapl)

      mlandflux = nlandflux
      mtrapl = ntrapl

      allocate (landflux(nlandflux))

      allocate (xemltrap(4,ntrapl))
      allocate (yemltrap(4,ntrapl))
      allocate (zemltrap(4,ntrapl))

      idims(1) = nlandflux

      call shdf5_irec(ndims,idims,'landflux%ilfglobe',ivara=landflux(:)%ilfglobe)
      call shdf5_irec(ndims,idims,'landflux%iw'      ,ivara=landflux(:)%iw)
      call shdf5_irec(ndims,idims,'landflux%kw'      ,ivara=landflux(:)%kw)
      call shdf5_irec(ndims,idims,'landflux%iwl'     ,ivara=landflux(:)%iwl)
      call shdf5_irec(ndims,idims,'landflux%jtrap'   ,ivara=landflux(:)%jtrap)
      call shdf5_irec(ndims,idims,'landflux%itrap'   ,ivara=landflux(:)%itrap)
      call shdf5_irec(ndims,idims,'landflux%area'    ,rvara=landflux(:)%area)
      call shdf5_irec(ndims,idims,'landflux%xef'     ,rvara=landflux(:)%xef)
      call shdf5_irec(ndims,idims,'landflux%yef'     ,rvara=landflux(:)%yef)
      call shdf5_irec(ndims,idims,'landflux%zef'     ,rvara=landflux(:)%zef)
      call shdf5_irec(ndims,idims,'landflux%arf_atm' ,rvara=landflux(:)%arf_atm)
      call shdf5_irec(ndims,idims,'landflux%arf_land',rvara=landflux(:)%arf_land)

      if (ntrapl > 0) then

         ndims = 2
         idims(1) = 4
         idims(2) = ntrapl

         call shdf5_irec(ndims,idims,'xemltrap' ,rvara=xemltrap)
         call shdf5_irec(ndims,idims,'yemltrap' ,rvara=yemltrap)
         call shdf5_irec(ndims,idims,'zemltrap' ,rvara=zemltrap)

      endif

   endif

! Close the GRIDFILE
   write(io6, *) '              ==T== Ler o gridfile: ',(walltime(wtime_start_gridfileread)-inicio)

   inicio = walltime(wtime_start_gridfileread)

   call shdf5_close()

   write(io6, *) '              ==T== Fechar o gridfile: ',(walltime(wtime_start_gridfileread)-inicio)
else

! Grid file does not exist.

   write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(io6,*) '!!!  Gridfile does not exist:'
   write(io6,*) '!!!  '//trim(gridfile)
   write(io6,*) '!!!  Stopping run'
   write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   
   stop 'stop - no gridfile'
   
endif

return
end subroutine gridfile_read



