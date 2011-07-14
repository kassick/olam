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
Module mem_sflux

! Total number of sea and land flux cells

   integer :: nseaflux, mseaflux
   integer :: nlandflux, mlandflux

   integer :: ntrapl, mtrapl ! number of trapezoids for all land flux cells
   integer :: ntraps, mtraps ! number of trapezoids for all sea flux cells

   real, allocatable :: xemltrap(:,:)
   real, allocatable :: yemltrap(:,:)
   real, allocatable :: zemltrap(:,:)

   real, allocatable :: xemstrap(:,:)
   real, allocatable :: yemstrap(:,:)
   real, allocatable :: zemstrap(:,:)

! Sea flux variables

   Type seaflux_vars
      logical, allocatable :: send_sf(:)

      integer :: isfglobe = 1
      integer :: iwrank = -1
      integer :: iw
      integer :: kw
      integer :: iws
      integer :: jtrap = 0  ! number of trapezoids in this flux cell
      integer :: itrap = 1  ! index of first trapezoid in this flux cell

      real :: dtf
      real :: area
      real :: xef,yef,zef
      real :: arf_atm
      real :: arf_sea
      real :: sfluxt  = 0.0
      real :: sfluxr  = 0.0
      real :: ustar   = 0.0
      real :: rhos    = 0.0
      real :: airtemp = 0.0
      real :: airshv  = 0.0
      real :: sxfer_t = 0.0
      real :: sxfer_r = 0.0
      real :: pcpg    = 0.0
      real :: qpcpg   = 0.0
      real :: dpcpg   = 0.0
      real :: rlong   = 0.0
      real :: rshort  = 0.0
      real :: rshort_diffuse = 0.0
   End Type

! Land flux variables

   Type landflux_vars
      logical, allocatable :: send_lf(:)

      integer :: ilfglobe = 1
      integer :: iwrank = -1
      integer :: iw
      integer :: kw
      integer :: iwl
      integer :: jtrap = 0  ! number of trapezoids in this flux cell
      integer :: itrap = 1  ! index of first trapezoid in this flux cell

      real :: dtf
      real :: area
      real :: xef,yef,zef
      real :: arf_atm
      real :: arf_land
      real :: sfluxt  = 0.0
      real :: sfluxr  = 0.0
      real :: ustar   = 0.0
      real :: rhos    = 0.0
      real :: prss    = 0.0
      real :: vels    = 0.0
      real :: airtemp = 0.0
      real :: airshv  = 0.0
      real :: sxfer_t = 0.0
      real :: sxfer_r = 0.0
      real :: pcpg    = 0.0
      real :: qpcpg   = 0.0
      real :: dpcpg   = 0.0
      real :: rlong   = 0.0
      real :: rshort  = 0.0
      real :: rshort_diffuse = 0.0
   End Type

   type (seaflux_vars), allocatable :: seaflux(:), seaflux_temp(:)
   type (landflux_vars), allocatable :: landflux(:), landflux_temp(:)

!----------------------------------------------------------------------------

   Type seafluxg_vars
      integer :: isf_myrank = -1
      integer :: irank = -1
   End Type seafluxg_vars

   Type landfluxg_vars
      integer :: ilf_myrank = -1
      integer :: irank = -1
   End Type landfluxg_vars

   type (seafluxg_vars),  allocatable, target :: seafluxg(:)
   type (landfluxg_vars), allocatable, target :: landfluxg(:)

!----------------------------------------------------------------------------

   Type jseaflux_vars
      integer, allocatable :: iseaflux(:)
      integer, allocatable :: jend(:)
   End Type

   Type jlandflux_vars
      integer, allocatable :: ilandflux(:)
      integer, allocatable :: jend(:)
   End Type

   type (jseaflux_vars)  :: jseaflux(12)  ! 2 + 10 remote nodes
   type (jlandflux_vars) :: jlandflux(12) ! 2 + 10 remote nodes

Contains

!===============================================================================

   subroutine filltab_sflux(mseaflux,mlandflux)

   use var_tables, only: vtables
   use misc_coms,  only: io6, iparallel
   
   implicit none

   integer, intent(in) :: mseaflux,mlandflux
   
   integer :: ndims
   integer, dimension(2) :: idims

! Fill pointers to arrays into variable tables

   ndims = 1
   idims(1) = mseaflux
   idims(2) = 1

   if (allocated(seaflux)) then
!------------------------------------------------------------------------
      if (iparallel == 1) call vtables(ndims,idims  &
           ,'SEAFLUX%ISFGLOBE :hist:nohist',ivara1=seaflux(:)%isfglobe)
!------------------------------------------------------------------------
      if (iparallel == 1) call vtables(ndims,idims  &
           ,'SEAFLUX%IWRANK :hist:nohist',ivara1=seaflux(:)%iwrank)
!------------------------------------------------------------------------
      call vtables(ndims,idims  &
           ,'SEAFLUX%SFLUXT :hist',rvara1=seaflux(:)%sfluxt)
!------------------------------------------------------------------------
      call vtables(ndims,idims  &
           ,'SEAFLUX%SFLUXR :hist',rvara1=seaflux(:)%sfluxr)
!------------------------------------------------------------------------
   endif

   idims(1) = mlandflux

   if (allocated(landflux)) then
!------------------------------------------------------------------------
      if (iparallel == 1)  call vtables(ndims,idims  &
           ,'LANDFLUX%ILFGLOBE :hist:nohist',ivara1=landflux(:)%ilfglobe)
!------------------------------------------------------------------------
      if (iparallel == 1)  call vtables(ndims,idims  &
           ,'LANDFLUX%IWRANK :hist:nohist',ivara1=landflux(:)%iwrank)
!------------------------------------------------------------------------
      call vtables(ndims,idims  &
           ,'LANDFLUX%SFLUXT :hist',rvara1=landflux(:)%sfluxt)
!------------------------------------------------------------------------
      call vtables(ndims,idims  &
           ,'LANDFLUX%SFLUXR :hist',rvara1=landflux(:)%sfluxr)
!------------------------------------------------------------------------
   endif

   return
   end subroutine filltab_sflux

!============================================================================

   subroutine fill_jflux()

   use mem_grid,   only: arw0, lpw, lsw, xem, yem, zem,  &
                         glatw, glonw, zm, topm
   use mem_ijtabs, only: itab_w, mrls, itabg_w
   use misc_coms,  only: io6, dtlm, iparallel
   use leaf_coms,  only: dt_leaf
   use mem_leaf,   only: itab_wl, itabg_wl
   use mem_sea,    only: itab_ws, itabg_ws
   use mem_para,   only: myrank, nsends_wsf, nsends_wlf
   use rastro_evts

   implicit none

   integer :: mrl,mrlw,isf,ilf,iw,iws,iwl,iloop,jsend,jend

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAM_FILL_JFLUX_IN,rst_buf)
#endif
! Allocate and initialize JSEAFLUX%JEND

   do iloop = 1,12
      allocate (jseaflux(iloop)%jend(mrls))
                jseaflux(iloop)%jend(1:mrls) = 0
   enddo

! Allocate ISEAFLUX member for first 2 JSEAFLUX data structures

   allocate (jseaflux(1)%iseaflux(mseaflux))
   allocate (jseaflux(2)%iseaflux(mseaflux))

   do mrl = mrls,1,-1
      do isf = 2,mseaflux
         iw  = seaflux(isf)%iw   ! global index
         iws = seaflux(isf)%iws  ! global index
      
! If run is parallel, convert iw to local domain

         if (iparallel == 1) then
            iw  = itabg_w (iw )%iw_myrank
            iws = itabg_ws(iws)%iws_myrank
         endif

         mrlw = itab_w(iw)%mrlw

         seaflux(isf)%dtf = min(dtlm(mrlw),dt_leaf)

         if (itab_w(iw)%mrlw == mrl) then
      
            if (iparallel == 0 .or.  &
               (iparallel == 1 .and. itab_w(iw)%irank == myrank)) then

               jseaflux(1)%jend(1:mrl) = jseaflux(1)%jend(1:mrl) + 1
               jseaflux(1)%iseaflux(jseaflux(1)%jend(1)) = isf
            endif

            if (iparallel == 0 .or.  &
               (iparallel == 1 .and. itab_ws(iws)%irank == myrank)) then

               jseaflux(2)%jend(1:mrl) = jseaflux(2)%jend(1:mrl) + 1
               jseaflux(2)%iseaflux(jseaflux(2)%jend(1)) = isf
            endif

         endif
      enddo
   enddo

! Allocate and initialize JLANDFLUX%JEND

   do iloop = 1,12
      allocate (jlandflux(iloop)%jend(mrls))
                jlandflux(iloop)%jend(1:mrls) = 0
   enddo

! Allocate ILANDFLUX member for first 2 JLANDFLUX structures

   allocate (jlandflux(1)%ilandflux(mlandflux))
   allocate (jlandflux(2)%ilandflux(mlandflux))

   do mrl = mrls,1,-1
      do ilf = 2,mlandflux
         iw  = landflux(ilf)%iw
         iwl = landflux(ilf)%iwl

! If run is parallel, convert iw to local domain

         if (iparallel == 1) then
            iw  = itabg_w (iw )%iw_myrank
            iwl = itabg_wl(iwl)%iwl_myrank
         endif

         mrlw = itab_w(iw)%mrlw

         landflux(ilf)%dtf = min(dtlm(mrlw),dt_leaf)

         if (itab_w(iw)%mrlw == mrl) then

            if (iparallel == 0 .or.  &
               (iparallel == 1 .and. itab_w(iw)%irank == myrank)) then

               jlandflux(1)%jend(1:mrl) = jlandflux(1)%jend(1:mrl) + 1
               jlandflux(1)%ilandflux(jlandflux(1)%jend(1)) = ilf
            endif

            if (iparallel == 0 .or.  &
               (iparallel == 1 .and. itab_wl(iwl)%irank == myrank)) then

               jlandflux(2)%jend(1:mrl) = jlandflux(2)%jend(1:mrl) + 1
               jlandflux(2)%ilandflux(jlandflux(2)%jend(1)) = ilf
            endif

         endif
      enddo
   enddo

! Return if run is not parallel

   if (iparallel == 0) return

! Compute and store JSEAFLUX%JEND(1)

   do jsend = 1,nsends_wsf(1)
      jseaflux(2+jsend)%jend(1) = 0
      do isf = 2,mseaflux
         if (seaflux(isf)%send_sf(jsend)) then
            jseaflux(2+jsend)%jend(1) = jseaflux(2+jsend)%jend(1) + 1
         endif
      enddo
      jseaflux(2+jsend)%jend(1) = max(1,jseaflux(2+jsend)%jend(1))
   enddo

! Allocate and zero-fill JSEAFLUX%ISEAFLUX

   do jsend = 1,nsends_wsf(1)
      jend = jseaflux(2+jsend)%jend(1)
      allocate (jseaflux(2+jsend)%iseaflux(jend))
                jseaflux(2+jsend)%iseaflux(1:jend) = 0
   enddo

! Initialize JSEAFLUX%JEND counters to zero

   do jsend = 1,nsends_wsf(1)
      jseaflux(2+jsend)%jend(1:mrls) = 0
   enddo

! Compute JSEAFLUX%ISEAFLUX

   do mrl = mrls,1,-1
      do isf = 2,mseaflux
         do jsend = 1,nsends_wsf(1)

            if (seaflux(isf)%send_sf(jsend)) then
               jseaflux(2+jsend)%jend(1:mrl) = jseaflux(2+jsend)%jend(1:mrl) + 1
               jseaflux(2+jsend)%iseaflux(jseaflux(2+jsend)%jend(1)) = isf
            endif

         enddo
      enddo
   enddo

! Compute and store JLANDFLUX%JEND(1)

   do jsend = 1,nsends_wlf(1)
      jlandflux(2+jsend)%jend(1) = 0
      do ilf = 2,mlandflux
         if (landflux(ilf)%send_lf(jsend)) then
            jlandflux(2+jsend)%jend(1) = jlandflux(2+jsend)%jend(1) + 1
         endif
      enddo
      jlandflux(2+jsend)%jend(1) = max(1,jlandflux(2+jsend)%jend(1))
   enddo

! Allocate and zero-fill JLANDFLUX%ILANDFLUX

   do jsend = 1,nsends_wlf(1)
      jend = jlandflux(2+jsend)%jend(1)
      allocate (jlandflux(2+jsend)%ilandflux(jend))
                jlandflux(2+jsend)%ilandflux(1:jend) = 0
   enddo

! Initialize JLANDFLUX%JEND counters to zero

   do jsend = 1,nsends_wlf(1)
      jlandflux(2+jsend)%jend(1:mrls) = 0
   enddo

! Compute JLANDFLUX%IWS

   do mrl = mrls,1,-1
      do ilf = 2,mlandflux
         do jsend = 1,nsends_wlf(1)

            if (landflux(ilf)%send_lf(jsend)) then
               jlandflux(2+jsend)%jend(1:mrl) = jlandflux(2+jsend)%jend(1:mrl) + 1
               jlandflux(2+jsend)%ilandflux(jlandflux(2+jsend)%jend(1)) = ilf
            endif

         enddo
      enddo
   enddo

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_FILL_JFLUX_OUT,rst_buf)
#endif

   return
   end subroutine fill_jflux

!============================================================================

subroutine init_fluxcells()

use mem_grid,   only: nwa, arw0, lpw, lsw, xem, yem, zem, xew, yew, zew, &
                      glatw, glonw, zm, topm
use mem_leaf,   only: land, itab_wl
use mem_sea,    only: sea, itab_ws
use mem_ijtabs, only: itab_w, mrls
use leaf_coms,  only: nwl, dt_leaf
use sea_coms,   only: nws
use misc_coms,  only: io6

implicit none

! local variables

! for now, only nqmax = 3 needed but later, land polygons will be more complex

integer, parameter :: npmax = 5, nqmax = 5

integer :: iws, iwl
integer :: isf, ilf
integer :: iw, im
integer :: j
integer :: im1, im2, im3
integer :: ks, k
integer :: mrl
integer :: np, nq
integer :: nseaflux_est, nlandflux_est, ntrapl_est, ntraps_est
integer :: incr_flux
integer :: incr_trap
integer :: jtrap, jt, it, imt

real :: x1,x2,x3,y1,y2,y3
real :: xp(npmax),yp(npmax)
real :: xq(nqmax),yq(nqmax)
real :: area
real :: xtrap(4,npmax+nqmax+npmax*nqmax)  ! trapezoid x coordinates
real :: ytrap(4,npmax+nqmax+npmax*nqmax)  ! trapezoid x coordinates
real :: traparea(npmax+nqmax+npmax*nqmax) ! trapezoid area

integer :: itmax, jtmax
integer :: iter
real :: trapareamax

real :: aatmin, aatmax, astmin, astmax, altmin, altmax
real :: dists, distn, wts, wtn

! Temporary scratch space during initialization only:

real, allocatable :: scrtrap(:,:)

real :: arf_atm_test(nwa)
real :: arf_sea_test(nws)
real :: arf_land_test(nwl)

real :: xeamin(nwa),xelmin(nwl),xesmin(nws)
real :: yeamin(nwa),yelmin(nwl),yesmin(nws)
real :: zeamin(nwa),zelmin(nwl),zesmin(nws)

real :: xeamax(nwa),xelmax(nwl),xesmax(nws)
real :: yeamax(nwa),yelmax(nwl),yesmax(nws)
real :: zeamax(nwa),zelmax(nwl),zesmax(nws)

integer :: imv(3), jm

! Estimate required array sizes and allocate arrays

incr_flux = 2 * nwa
incr_trap = 5 * nwa

nseaflux_est = incr_flux
nlandflux_est = incr_flux

ntrapl_est = incr_trap
ntraps_est = incr_trap

allocate (seaflux(nseaflux_est))
allocate (landflux(nlandflux_est))

allocate (xemltrap(4,ntrapl_est))
allocate (yemltrap(4,ntrapl_est))
allocate (zemltrap(4,ntrapl_est))

allocate (xemstrap(4,ntraps_est))
allocate (yemstrap(4,ntraps_est))
allocate (zemstrap(4,ntraps_est))

! new - find maximum extents of all atm, land, and sea cells

do iw = 2,nwa
   imv = (/ itab_w(iw)%im1, itab_w(iw)%im2, itab_w(iw)%im3 /)

   xeamin(iw) = minval(xem(imv))
   yeamin(iw) = minval(yem(imv))
   zeamin(iw) = minval(zem(imv))

   xeamax(iw) = maxval(xem(imv))
   yeamax(iw) = maxval(yem(imv))
   zeamax(iw) = maxval(zem(imv))
enddo

do iwl = 2,nwl
   jm = itab_wl(iwl)%jm

   xelmin(iwl) = minval(land%xeml(itab_wl(iwl)%im(1:jm)))
   yelmin(iwl) = minval(land%yeml(itab_wl(iwl)%im(1:jm)))
   zelmin(iwl) = minval(land%zeml(itab_wl(iwl)%im(1:jm)))

   xelmax(iwl) = maxval(land%xeml(itab_wl(iwl)%im(1:jm)))
   yelmax(iwl) = maxval(land%yeml(itab_wl(iwl)%im(1:jm)))
   zelmax(iwl) = maxval(land%zeml(itab_wl(iwl)%im(1:jm)))
enddo

do iws = 2,nws
   jm = itab_ws(iws)%jm

   xesmin(iws) = minval(sea%xems(itab_ws(iws)%im(1:jm)))
   yesmin(iws) = minval(sea%yems(itab_ws(iws)%im(1:jm)))
   zesmin(iws) = minval(sea%zems(itab_ws(iws)%im(1:jm)))

   xesmax(iws) = maxval(sea%xems(itab_ws(iws)%im(1:jm)))
   yesmax(iws) = maxval(sea%yems(itab_ws(iws)%im(1:jm)))
   zesmax(iws) = maxval(sea%zems(itab_ws(iws)%im(1:jm)))
enddo

np = 3

! Initialize sea and land flux cell counters to 1

isf = 1
ilf = 1

ntrapl = 0
ntraps = 0

landflux(1)%ilfglobe = 1
landflux(1)%iw = 1
landflux(1)%iwl = 1
seaflux(1)%isfglobe = 1
seaflux(1)%iw = 1
seaflux(1)%iws = 1

! Loop over all W points

do iw = 2,nwa

   if (.not. itab_w(iw)%loop(20)) cycle  ! If sfc fluxes not done at this IW

   im1 = itab_w(iw)%im1
   im2 = itab_w(iw)%im2
   im3 = itab_w(iw)%im3

! Evaluate x,y coordinates of ATM cell M points on polar stereographic plane
! tangent at IW

   call e_ps(xem(im1),yem(im1),zem(im1),glatw(iw),glonw(iw),x1,y1)
   call e_ps(xem(im2),yem(im2),zem(im2),glatw(iw),glonw(iw),x2,y2)
   call e_ps(xem(im3),yem(im3),zem(im3),glatw(iw),glonw(iw),x3,y3)

! Loop over vertical T levels in IW that are intersected by topography

   do ks = 1,lsw(iw)
      k = lpw(iw) + ks - 1

! Find polygon of intersection between topography and IW triangle for current
! K level

      call sfc_polygon(zm(k-1),zm(k),x1,x2,x3,y1,y2,y3,  &
         topm(im1),topm(im2),topm(im3),np,xp,yp)

! Loop over all WL points

      do iwl = 2,nwl

! new - skip interaction using new non-overlap check

         if (abs(xew(iw) + land%xewl(iwl)) < 12.e6) then
            if (xeamin(iw ) > xelmax(iwl)) cycle
            if (xelmin(iwl) > xeamax(iw )) cycle
         endif
         
         if (abs(yew(iw) + land%yewl(iwl)) < 12.e6) then
            if (yeamin(iw ) > yelmax(iwl)) cycle
            if (yelmin(iwl) > yeamax(iw )) cycle
         endif

         if (abs(zew(iw) + land%zewl(iwl)) < 12.e6) then
            if (zeamin(iw ) > zelmax(iwl)) cycle
            if (zelmin(iwl) > zeamax(iw )) cycle
         endif

         nq = itab_wl(iwl)%jm

         do j = 1,itab_wl(iwl)%jm
            im = itab_wl(iwl)%im(j)

! Evaluate x,y coordinates of LAND cell M points on polar stereographic plane
! tangent at IW

            call e_ps(land%xeml(im),land%yeml(im),land%zeml(im),  &
                      glatw(iw),glonw(iw),xq(j),yq(j))
         enddo

! Evaluate possible overlap of ATM and LAND polygons

         call polygon_overlap(iwl,np,nq,xp,yp,xq,yq,area,jtrap,xtrap,ytrap,traparea)

         if (area > 1.e-6 * min(arw0(iw),land%area(iwl))) then

! This IW/K trapezoid overlaps with land cell IWL. Loop over all trapezoids 
! of current overlap

            trapareamax = 0.

            do jt = 1,jtrap
               it = ntrapl + jt

! Find trapezoid of largest area 

               if (traparea(jt) > trapareamax) then
                  trapareamax = traparea(jt)
                  jtmax = jt
                  itmax = it
               endif   

! Convert xtrap,ytrap from ps to earth coordinates

               do imt = 1,4
                  call ps_e(xemltrap(imt,it),           &
                            yemltrap(imt,it),           &
                            zemltrap(imt,it),           &
                            glatw(iw),glonw(iw),        &
                            xtrap(imt,jt),ytrap(imt,jt) )

               enddo
            enddo

! Define 2 weights for positioning xef, yef, zef within largest trapezoid
     
            dists = xtrap(2,jtmax) - xtrap(1,jtmax)     
            distn = xtrap(3,jtmax) - xtrap(4,jtmax)

            wtn = .333 + .333 * distn / (dists + distn)
            wts = 1. - wtn

! Increment landflux cell index and fill cell values

            ilf = ilf + 1

            landflux(ilf)%jtrap = jtrap
            landflux(ilf)%itrap = ntrapl + 1

            landflux(ilf)%ilfglobe = ilf ! full domain index
            landflux(ilf)%iw       = iw  ! full domain index
            landflux(ilf)%kw       = k
            landflux(ilf)%iwl      = iwl ! full domain index
            landflux(ilf)%area     = area

            landflux(ilf)%xef      = .5 * (wts * xemltrap(1,itmax)  &
                                        +  wts * xemltrap(2,itmax)  &
                                        +  wtn * xemltrap(3,itmax)  &
                                        +  wtn * xemltrap(4,itmax))
            landflux(ilf)%yef      = .5 * (wts * yemltrap(1,itmax)  &
                                        +  wts * yemltrap(2,itmax)  &
                                        +  wtn * yemltrap(3,itmax)  &
                                        +  wtn * yemltrap(4,itmax))
            landflux(ilf)%zef      = .5 * (wts * zemltrap(1,itmax)  &
                                        +  wts * zemltrap(2,itmax)  &
                                        +  wtn * zemltrap(3,itmax)  &
                                        +  wtn * zemltrap(4,itmax))
            landflux(ilf)%arf_atm  = landflux(ilf)%area / arw0(iw)
            landflux(ilf)%arf_land = landflux(ilf)%area / land%area(iwl)

! Increment ntrapl

            ntrapl = ntrapl + jtrap

         endif

      enddo ! iwl loop
      
! Loop over all WS points

      do iws = 2,nws

! new - skip interaction using new non-overlap check

         if (abs(xew(iw) + sea%xews(iws)) < 12.e6) then
            if (xeamin(iw ) > xesmax(iws)) cycle
            if (xesmin(iws) > xeamax(iw )) cycle
         endif

         if (abs(yew(iw) + sea%yews(iws)) < 12.e6) then
            if (yeamin(iw ) > yesmax(iws)) cycle
            if (yesmin(iws) > yeamax(iw )) cycle
         endif

         if (abs(zew(iw) + sea%zews(iws)) < 12.e6) then         
            if (zeamin(iw ) > zesmax(iws)) cycle
            if (zesmin(iws) > zeamax(iw )) cycle
         endif

         nq = itab_ws(iws)%jm

         do j = 1,itab_ws(iws)%jm
            im = itab_ws(iws)%im(j)

! Evaluate x,y coordinates of SEA cell M points on polar stereographic plane
! tangent at IW

            call e_ps(sea%xems(im),sea%yems(im),sea%zems(im),  &
                      glatw(iw),glonw(iw),xq(j),yq(j))

         enddo

! Evaluate possible overlap of ATM and SEA polygons

         call polygon_overlap(iws,np,nq,xp,yp,xq,yq,area,  &
             jtrap,xtrap,ytrap,traparea)

         if (area > 1.e-6 * min(arw0(iw),sea%area(iws))) then

! This IW/K trapezoid overlaps with sea cell IWS.  Loop over all trapezoids 
! of current overlap

            trapareamax = 0.

            do jt = 1,jtrap
               it = ntraps + jt
               
! Find trapezoid of largest area 

               if (traparea(jt) > trapareamax) then
                  trapareamax = traparea(jt)
                  itmax = it
                  jtmax = jt
               endif   

! Convert xtrap,ytrap from ps to earth coordinates

               do imt = 1,4
                  call ps_e(xemstrap(imt,it),           &
                            yemstrap(imt,it),           &
                            zemstrap(imt,it),           &
                            glatw(iw),glonw(iw),        &
                            xtrap(imt,jt),ytrap(imt,jt) )
               enddo

            enddo

! Define 2 weights for positioning xef, yef, zef within largest trapezoid
     
            dists = xtrap(2,jtmax) - xtrap(1,jtmax)     
            distn = xtrap(3,jtmax) - xtrap(4,jtmax)

            wtn = .333 + .333 * distn / (dists + distn)
            wts = 1. - wtn

! Increment seaflux cell index and fill cell values

            isf = isf + 1

            seaflux(isf)%jtrap = jtrap
            seaflux(isf)%itrap = ntraps + 1

            seaflux(isf)%isfglobe = isf  ! full domain index
            seaflux(isf)%iw       = iw  ! full domain index
            seaflux(isf)%kw       = k
            seaflux(isf)%iws      = iws ! full domain index
            seaflux(isf)%area     = area
            seaflux(isf)%xef      = .5 * (wts * xemstrap(1,itmax)  &
                                       +  wts * xemstrap(2,itmax)  &
                                       +  wtn * xemstrap(3,itmax)  &
                                       +  wtn * xemstrap(4,itmax))
            seaflux(isf)%yef      = .5 * (wts * yemstrap(1,itmax)  &
                                       +  wts * yemstrap(2,itmax)  &
                                       +  wtn * yemstrap(3,itmax)  &
                                       +  wtn * yemstrap(4,itmax))
            seaflux(isf)%zef      = .5 * (wts * zemstrap(1,itmax)  &
                                       +  wts * zemstrap(2,itmax)  &
                                       +  wtn * zemstrap(3,itmax)  &
                                       +  wtn * zemstrap(4,itmax))
            seaflux(isf)%arf_atm  = seaflux(isf)%area / arw0(iw)
            seaflux(isf)%arf_sea  = seaflux(isf)%area / sea%area(iws)

! Increment ntraps

            ntraps = ntraps + jtrap

         endif
         
      enddo ! iws loop

   enddo ! ks loop
   
! If number of flux cells or trapezoids is getting close to allocated 
! quantity, allocate more space

   if (nseaflux_est - isf < 1000) then 

      nseaflux_est = nseaflux_est + incr_flux
      allocate (seaflux_temp(isf))
      seaflux_temp(1:isf) = seaflux(1:isf)

      deallocate (seaflux)
      allocate (seaflux(nseaflux_est))
      seaflux(1:isf) = seaflux_temp(1:isf)
      
      deallocate (seaflux_temp)
      
   endif

   if (nlandflux_est - ilf < 1000) then 

      nlandflux_est = nlandflux_est + incr_flux
      allocate (landflux_temp(ilf))
      landflux_temp(1:ilf) = landflux(1:ilf)

      deallocate (landflux)
      allocate (landflux(nlandflux_est))
      landflux(1:ilf) = landflux_temp(1:ilf)
      
      deallocate (landflux_temp)
      
   endif

   if (ntraps_est - ntraps < 1000) then 

      ntraps_est = ntraps_est + incr_trap
      allocate (scrtrap(4,ntraps))

      scrtrap(1:4,1:ntraps) = xemstrap(1:4,1:ntraps)
      deallocate (xemstrap)
      allocate (xemstrap(4,ntraps_est))
      xemstrap(1:4,1:ntraps) = scrtrap(1:4,1:ntraps)
      
      scrtrap(1:4,1:ntraps) = yemstrap(1:4,1:ntraps)
      deallocate (yemstrap)
      allocate (yemstrap(4,ntraps_est))
      yemstrap(1:4,1:ntraps) = scrtrap(1:4,1:ntraps)
      
      scrtrap(1:4,1:ntraps) = zemstrap(1:4,1:ntraps)
      deallocate (zemstrap)
      allocate (zemstrap(4,ntraps_est))
      zemstrap(1:4,1:ntraps) = scrtrap(1:4,1:ntraps)
      
      deallocate (scrtrap)
      
   endif      

   if (ntrapl_est - ntrapl < 1000) then 

      ntrapl_est = ntrapl_est + incr_trap
      allocate (scrtrap(4,ntrapl))

      scrtrap(1:4,1:ntrapl) = xemltrap(1:4,1:ntrapl)
      deallocate (xemltrap)
      allocate (xemltrap(4,ntrapl_est))
      xemltrap(1:4,1:ntrapl) = scrtrap(1:4,1:ntrapl)
      
      scrtrap(1:4,1:ntrapl) = yemltrap(1:4,1:ntrapl)
      deallocate (yemltrap)
      allocate (yemltrap(4,ntrapl_est))
      yemltrap(1:4,1:ntrapl) = scrtrap(1:4,1:ntrapl)
      
      scrtrap(1:4,1:ntrapl) = zemltrap(1:4,1:ntrapl)
      deallocate (zemltrap)
      allocate (zemltrap(4,ntrapl_est))
      zemltrap(1:4,1:ntrapl) = scrtrap(1:4,1:ntrapl)
      
      deallocate (scrtrap)
      
   endif
      
enddo ! iw loop

nseaflux = isf
nlandflux = ilf

mseaflux = isf
mlandflux = ilf

! Check sums of arf values to make sure they are close to 1.0

do iter = 1,2

   arf_atm_test(1:nwa) = 0.
   arf_sea_test(1:nws) = 0.
   arf_land_test(1:nwl) = 0.

   do isf = 2,nseaflux
      iw  = seaflux(isf)%iw
      iws = seaflux(isf)%iws

      arf_atm_test(iw)  = arf_atm_test(iw)  + seaflux(isf)%arf_atm
      arf_sea_test(iws) = arf_sea_test(iws) + seaflux(isf)%arf_sea
   enddo

   do ilf = 2,nlandflux
      iw  = landflux(ilf)%iw
      iwl = landflux(ilf)%iwl

      arf_atm_test(iw)   = arf_atm_test(iw)   + landflux(ilf)%arf_atm
      arf_land_test(iwl) = arf_land_test(iwl) + landflux(ilf)%arf_land
   enddo

   aatmin = 1.e9
   aatmax = -1.e9

   astmin = 1.e9
   astmax = -1.e9

   altmin = 1.e9
   altmax = -1.e9

   do iw = 2,nwa
      if (aatmin > arf_atm_test(iw)) aatmin = arf_atm_test(iw)
      if (aatmax < arf_atm_test(iw)) aatmax = arf_atm_test(iw)
   enddo

   do iws = 2,nws
      if (astmin > arf_sea_test(iws)) astmin = arf_sea_test(iws)
      if (astmax < arf_sea_test(iws)) astmax = arf_sea_test(iws)
   enddo

   do iwl = 2,nwl
      if (altmin > arf_land_test(iwl)) altmin = arf_land_test(iwl)
      if (altmax < arf_land_test(iwl)) altmax = arf_land_test(iwl)
   enddo

   if (iter == 1) then

      do isf = 2,nseaflux
         iw  = seaflux(isf)%iw
         iws = seaflux(isf)%iws

         seaflux(isf)%arf_atm = seaflux(isf)%arf_atm / arf_atm_test(iw) 
         seaflux(isf)%arf_sea = seaflux(isf)%arf_sea / arf_sea_test(iws)
      enddo

      do ilf = 2,nlandflux
         iw  = landflux(ilf)%iw
         iwl = landflux(ilf)%iwl

         landflux(ilf)%arf_atm  = landflux(ilf)%arf_atm  / arf_atm_test(iw) 
         landflux(ilf)%arf_land = landflux(ilf)%arf_land / arf_land_test(iwl)
      enddo
      
   endif

enddo

return
end subroutine init_fluxcells

!============================================================================

subroutine sfc_polygon(zm1,zm2,x1,x2,x3,y1,y2,y3,z1,z2,z3,np,xp,yp)

implicit none

real, intent(in) :: zm1,zm2  ! Bottom, top of current model T level
real, intent(in) :: z1,z2,z3 ! Topography heights at 3 M pts
real, intent(in) :: x1,x2,x3 ! x ps coords of 3 M pts
real, intent(in) :: y1,y2,y3 ! y ps coords of 3 M pts

integer, intent(out) :: np ! Number of polygon pts
real, intent(out) :: xp(5),yp(5) ! x,y ps coords of polygon pts
  
if     (z1 <= z2 .and. z1 <= z3) then
   call sfc_polyg(zm1,zm2,z1,z2,z3,x1,x2,x3,y1,y2,y3,np,xp,yp)
elseif (z2 <= z3 .and. z2 <= z1) then
   call sfc_polyg(zm1,zm2,z2,z3,z1,x2,x3,x1,y2,y3,y1,np,xp,yp)
else
   call sfc_polyg(zm1,zm2,z3,z1,z2,x3,x1,x2,y3,y1,y2,np,xp,yp)
endif

return
end subroutine sfc_polygon

!============================================================================

subroutine sfc_polyg(zm1,zm2,z1,z2,z3,x1,x2,x3,y1,y2,y3,np,xp,yp)

use misc_coms,  only: io6

implicit none

real, intent(in) :: zm1,zm2  ! Bottom, top of current model T level
real, intent(in) :: z1,z2,z3 ! Topo heights at 3 M pts, cyclic order, z1 lowest
real, intent(in) :: x1,x2,x3 ! x ps coords of 3 M pts
real, intent(in) :: y1,y2,y3 ! y ps coords of 3 M pts

integer, intent(out) :: np ! Number of polygon pts
real, intent(out) :: xp(5),yp(5) ! x,y ps coords of polygon pts

if (zm1 <= z1) then

   xp(1) = x1
   yp(1) = y1

   if (zm2 <= z2) then

      xp(2) = x1 + (x2-x1) * (zm2-z1) / (z2-z1)
      yp(2) = y1 + (y2-y1) * (zm2-z1) / (z2-z1)

      if (zm2 <= z3) then

         xp(3) = x1 + (x3-x1) * (zm2-z1) / (z3-z1)
         yp(3) = y1 + (y3-y1) * (zm2-z1) / (z3-z1)

         np = 3

      else
      
         xp(3) = x3 + (x2-x3) * (zm2-z3) / (z2-z3)
         yp(3) = y3 + (y2-y3) * (zm2-z3) / (z2-z3)

         xp(4) = x3
         yp(4) = y3

         np = 4
         
      endif         
      
   else
      
      xp(2) = x2
      yp(2) = y2

      if (zm2 <= z3) then

         xp(3) = x2 + (x3-x2) * (zm2-z2) / (z3-z2)
         yp(3) = y2 + (y3-y2) * (zm2-z2) / (z3-z2)

         xp(4) = x1 + (x3-x1) * (zm2-z1) / (z3-z1)
         yp(4) = y1 + (y3-y1) * (zm2-z1) / (z3-z1)

         np = 4

      else
      
         xp(3) = x3
         yp(3) = y3

         np = 3
         
      endif         
      
   endif

elseif (zm1 <= z2) then

   xp(1) = x1 + (x2-x1) * (zm1-z1) / (z2-z1)
   yp(1) = y1 + (y2-y1) * (zm1-z1) / (z2-z1)


   if (zm2 <= z2) then

      xp(2) = x1 + (x2-x1) * (zm2-z1) / (z2-z1)
      yp(2) = y1 + (y2-y1) * (zm2-z1) / (z2-z1)

      if (zm2 <= z3) then

         xp(3) = x1 + (x3-x1) * (zm2-z1) / (z3-z1)
         yp(3) = y1 + (y3-y1) * (zm2-z1) / (z3-z1)

         xp(4) = x1 + (x3-x1) * (zm1-z1) / (z3-z1)
         yp(4) = y1 + (y3-y1) * (zm1-z1) / (z3-z1)

         np = 4
            
      else

         xp(3) = x3 + (x2-x3) * (zm2-z3) / (z2-z3)
         yp(3) = y3 + (y2-y3) * (zm2-z3) / (z2-z3)

         if (zm1 <= z3) then

            xp(4) = x3
            yp(4) = y3

            xp(5) = x1 + (x3-x1) * (zm1-z1) / (z3-z1)
            yp(5) = y1 + (y3-y1) * (zm1-z1) / (z3-z1)

            np = 5
               
         else
            
            xp(4) = x3 + (x2-x3) * (zm1-z3) / (z2-z3)
            yp(4) = y3 + (y2-y3) * (zm1-z3) / (z2-z3)

            np = 4
            
         endif
            
      endif
          
   else

      xp(2) = x2
      yp(2) = y2

      if (zm2 <= z3) then

         xp(3) = x2 + (x3-x2) * (zm2-z2) / (z3-z2)
         yp(3) = y2 + (y3-y2) * (zm2-z2) / (z3-z2)

         xp(4) = x1 + (x3-x1) * (zm2-z1) / (z3-z1)
         yp(4) = y1 + (y3-y1) * (zm2-z1) / (z3-z1)

         xp(5) = x1 + (x3-x1) * (zm1-z1) / (z3-z1)
         yp(5) = y1 + (y3-y1) * (zm1-z1) / (z3-z1)

         np = 5
               
      else
         
         if (zm1 <= z3) then

            xp(3) = x3
            yp(3) = y3

            xp(4) = x1 + (x3-x1) * (zm1-z1) / (z3-z1)
            yp(4) = y1 + (y3-y1) * (zm1-z1) / (z3-z1)

            np = 4
               
         else
            
            xp(3) = x3 + (x2-x3) * (zm1-z3) / (z2-z3)
            yp(3) = y3 + (y2-y3) * (zm1-z3) / (z2-z3)

            np = 3
            
         endif

      endif
      
   endif

else

   xp(1) = x2 + (x3-x2) * (zm1-z2) / (z3-z2)
   yp(1) = y2 + (y3-y2) * (zm1-z2) / (z3-z2)

   if (zm2 <= z3) then

      xp(2) = x2 + (x3-x2) * (zm2-z2) / (z3-z2)
      yp(2) = y2 + (y3-y2) * (zm2-z2) / (z3-z2)

      xp(3) = x1 + (x3-x1) * (zm2-z1) / (z3-z1)
      yp(3) = y1 + (y3-y1) * (zm2-z1) / (z3-z1)

      xp(4) = x1 + (x3-x1) * (zm1-z1) / (z3-z1)
      yp(4) = y1 + (y3-y1) * (zm1-z1) / (z3-z1)

      np = 4
               
   else   

      xp(2) = x3
      yp(2) = y3

      xp(3) = x1 + (x3-x1) * (zm1-z1) / (z3-z1)
      yp(3) = y1 + (y3-y1) * (zm1-z1) / (z3-z1)

      np = 3

   endif
            
endif

return
end subroutine sfc_polyg

!============================================================================

subroutine polygon_overlap(iwls,np,nq,xp,yp,xq,yq,area,ntrap,xtrap,ytrap,traparea)

! Given x,y coordinates of the vertices of polygons p and q, compute the area
! of overlap between the polygons using a sweepline algorithm.

! Method adapted from:
! Zerzan, J., 1989, Computers & Geosciences, Vol. 15, No. 7, pp. 1109-1114.

implicit none

integer, intent(in) :: np,nq,iwls ! Number of vertices in p and q
real, intent(in) :: xp(np),yp(np) ! x,y coordinates of p vertices
real, intent(in) :: xq(nq),yq(nq) ! x,y coordinates of q vertices

integer, intent(out) :: ntrap              ! number of trapezoids
real, intent(out) :: xtrap(4,np+nq+np*nq)  ! trapezoid x coordinates
real, intent(out) :: ytrap(4,np+nq+np*nq)  ! trapezoid y coordinates
real, intent(out) :: traparea(np+nq+np*nq) ! trapezoid area
real, intent(out) :: area                  ! area of overlap of p and q

real :: yev(np+nq+np*nq)  ! y-coordinates of event

integer :: nev      ! # of events
integer :: nsect    ! # of intersections between strip centerline and p,q edges
real :: ymid        ! y-coord of centerline of strip between consecutive events
real :: xmid(np+nq) ! x-coords where strip centerline intersects p and q edges 
real :: x1(np+nq)   ! x-coords where strip .1 line intersects p and q edges 
real :: xcent       ! x-coord of midpoint of segment between xmid values

integer :: ip,iq,ipa,ipb,iqa,iqb,iflag,iev,i,ia,ib,is
real :: p0,q0,dx,dy,y1,dxtrap

real :: alpha  ! angle swept by ray from point to polygon perimeter during circuit

! Find vertices in p that are not outside q

nev = 0
area = 0.
ntrap = 0

do ip = 1,np
   call inout_check(nq,xq,yq,xp(ip),yp(ip),alpha)

   if (abs(alpha) > .2) then 
      nev = nev + 1
      yev(nev) = yp(ip)
   endif
enddo

! Find vertices in q that are not outside p

do iq = 1,nq
   call inout_check(np,xp,yp,xq(iq),yq(iq),alpha)

   if (abs(alpha) > .2) then 
      nev = nev + 1
      yev(nev) = yq(iq)
   endif
enddo

! Find intersecting edges of polygons p and q

do ipa = 1,np
   ipb = ipa + 1
   if (ipa == np) ipb = 1

   do iqa = 1,nq
      iqb = iqa + 1
      if (iqa == nq) iqb = 1

      call intersect(xp(ipa),yp(ipa),xp(ipb),yp(ipb)  &
                    ,xq(iqa),yq(iqa),xq(iqb),yq(iqb),p0,q0,iflag)

      if (iflag == -1) cycle
      if (p0 < -0.000001 .or. p0 > 1.000001) cycle
      if (q0 < -0.000001 .or. q0 > 1.000001) cycle

! Line segments pa-pb and qa-qb intersect; find y-coordinate of intersection

      nev = nev + 1
      yev(nev) = (1. - p0) * yp(ipa) + p0 * yp(ipb)

   enddo
enddo

if (nev == 0) then
!   print*, 'no overlap'
   return
endif

! Sort event list to obtain increasing order of yev

call fltsort(nev,yev)

! Loop over event points

do iev = 1,nev-1

! dy = width of strip between current and next event

   dy = yev(iev+1) - yev(iev)
   
! Reject overlap event if dy is less than 1 meter (threshold ok for single precision)
   
   if (dy < 1.) cycle

! ymid = y-coordinate of strip centerline.
! Initialize dx (centerline length sum) to 0.

   ymid = yev(iev) + .5 * dy
   y1   = yev(iev) + .1 * dy
   dx = 0.

! Find x-coordinate of intersections of strip centerline with edges of p and q

   nsect = 0

   do ia = 1,np
      ib = ia + 1
      if (ia == np) ib = 1

      if (ymid < min(yp(ia),yp(ib)) .or. ymid > max(yp(ia),yp(ib))) cycle

      nsect = nsect + 1
      xmid(nsect) = xp(ia)  &
                  + (xp(ib) - xp(ia)) * (ymid - yp(ia)) / (yp(ib) - yp(ia))
      x1(nsect) = xp(ia)  &
                + (xp(ib) - xp(ia)) * (y1 - yp(ia)) / (yp(ib) - yp(ia))
   enddo

   do ia = 1,nq
      ib = ia + 1
      if (ia == nq) ib = 1

      if (ymid < min(yq(ia),yq(ib)) .or. ymid > max(yq(ia),yq(ib))) cycle

      nsect = nsect + 1
      xmid(nsect) = xq(ia)  &
                  + (xq(ib) - xq(ia)) * (ymid - yq(ia)) / (yq(ib) - yq(ia))
      x1(nsect) = xq(ia)  &
                + (xq(ib) - xq(ia)) * (y1 - yq(ia)) / (yq(ib) - yq(ia))
   enddo

! Sort xmid values into increasing order

   call fltsort2(nsect,xmid,x1)

! see if the segment is inside both polygons

   do is = 1,nsect - 1
      xcent = .5 * (xmid(is) + xmid(is+1))

      if (xcent == xmid(is)) cycle          ! if segment length is 0

      call inout_check(np,xp,yp,xcent,ymid,alpha)

      if (abs(alpha) < 1.) cycle

      call inout_check(nq,xq,yq,xcent,ymid,alpha)

      if (abs(alpha) < 1.) cycle

      dxtrap = xmid(is+1) - xmid(is)

      dx = dx + dxtrap

! Find x at 4 corners of current trapezoidal region

      ntrap = ntrap + 1
      
      xtrap(1,ntrap) =  1.25 * x1(is)   -  .25 * xmid(is) 
      xtrap(2,ntrap) =  1.25 * x1(is+1) -  .25 * xmid(is+1) 
      xtrap(3,ntrap) = -1.25 * x1(is+1) + 2.25 * xmid(is+1) 
      xtrap(4,ntrap) = -1.25 * x1(is)   + 2.25 * xmid(is)
      
      ytrap(1,ntrap) = yev(iev) 
      ytrap(2,ntrap) = yev(iev) 
      ytrap(3,ntrap) = yev(iev+1) 
      ytrap(4,ntrap) = yev(iev+1) 

      traparea(ntrap) = dxtrap * dy

   enddo

   area = area + dx * dy

enddo

return
end subroutine polygon_overlap

!============================================================================

subroutine inout_check(n,x,y,x0,y0,th1)

! Given planar Cartesian coordinates of the vertices of a simple closed polygon,
! x(1),y(1),...,x(n),y(n), and of an additional point, x0,y0, determine
! whether the point is inside, outside, or on the border of the polygon.

! Set:
! th1 = 0 if outside
! th1 = pi if on the border
! th1 = 2 * pi if inside

implicit none

integer, intent(in) :: n  ! number of points in polygon

real, intent(in) :: x(n),y(n)  ! x,y coordinates of polygon points      
real, intent(in) :: x0,y0      ! x,y coordinates of additional point
real, intent(out) :: th1       ! output value

! Local variables

integer :: i
real :: theta
real :: xh1, xh2
real :: edgemax
real :: x1(n+1), y1(n+1)

real, parameter :: pi1 = 3.1415926536, pi2 = 2. * pi1

do i = 1,n
   x1(i) = x(i) - x0
   y1(i) = y(i) - y0

   ! Check whether x0,y0 is exactly on a vertex
   if ((abs(x1(i)) < 1.e-6) .and. (abs(y1(i)) < 1.e-6)) then
      th1 = pi1
      return
   endif
enddo

x1(n+1) = x1(1)
y1(n+1) = y1(1)
th1 = 0.

xh2 = atan2(y1(1),x1(1))
if (xh2 < 0.) xh2 = xh2 + pi2

do i = 1,n
   xh1 = atan2(y1(i+1),x1(i+1))
   if (xh1 < 0.) xh1 = xh1 + pi2

   theta = xh1 - xh2
   if (theta < -pi1) theta = theta + pi2

   if ((abs(abs(theta) - pi1) < .00001)) then
      th1 = pi1
      return
   endif

   if (theta > pi1) theta = theta - pi2
   th1 = th1 + theta
   
   xh2 = xh1
enddo

th1 = abs(th1)

return
end subroutine inout_check

!============================================================================

subroutine intersect(xpa,ypa,xpb,ypb,xqa,yqa,xqb,yqb,p0,q0,iflag)

! Given x,y coordinates of points pa and pb on line p and points qa and qb on
! line q, find parameteric values p0 and q0 where p and q intersect.
! If no intersection, set iflag = -1.

implicit none

real, intent(in) :: xpa,ypa,xpb,ypb  ! x,y coordinates of pa and pb
real, intent(in) :: xqa,yqa,xqb,yqb  ! x,y coordinates of qa and qb

integer, intent(out) :: iflag
real, intent(out) :: p0,q0

real :: pabxqab  ! cross product of vectors pa-pb and qa-qb 

iflag = -1

pabxqab = (xpb - xpa) * (yqb - yqa) - (ypb - ypa) * (xqb - xqa)

if (pabxqab /= 0.) then  ! Infinite lines intersect
   iflag = 0

   p0 = ((yqb - yqa) * (xqa - xpa) + (yqa - ypa) * (xqa - xqb)) / pabxqab
   q0 = ((ypb - ypa) * (xqa - xpa) + (yqa - ypa) * (xpa - xpb)) / pabxqab
endif

!print*, 'pabxqab ',pabxqab,p0,q0,iflag

return
end subroutine intersect

!============================================================================

subroutine fltsort(n,f)

! Sort n floating point numbers f into ascending order

implicit none

integer, intent(in) :: n
real, intent(inout) :: f(n)

integer :: i,j

real :: f0

do i = 1,n-1
   do j = i+1,n
      if (f(j) < f(i)) then
         f0 = f(i)
         f(i) = f(j)
         f(j) = f0
      endif
   enddo
enddo

return
end subroutine fltsort
   
!============================================================================

subroutine fltsort2(n,f1,f2)

! Sort n floating point numbers in each of f1 and f2 into ascending order by f1

implicit none

integer, intent(in) :: n
real, intent(inout) :: f1(n),f2(n)

integer :: i,j

real :: f0

do i = 1,n-1
   do j = i+1,n
      if (f1(j) < f1(i)) then
         f0 = f1(i)
         f1(i) = f1(j)
         f1(j) = f0

         f0 = f2(i)
         f2(i) = f2(j)
         f2(j) = f0
      endif
   enddo
enddo

return
end subroutine fltsort2

End Module mem_sflux
