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
subroutine para_decomp()

! Decompose global grid into multiple subdomains for parallel computation

use mem_para,   only: mgroupsize, myrank
use misc_coms,  only: io6

use mem_ijtabs, only: itab_u, itab_w, itabg_m, itabg_u, itabg_w
use mem_grid,   only: nma, nua, nwa, xem, yem, zem

use leaf_coms,  only: nml, nul, nwl, isfcl
use mem_leaf,   only: itab_ul, itab_wl, itabg_ml, itabg_ul, itabg_wl, land

use sea_coms,   only: nms, nus, nws
use mem_sea,    only: itab_us, itab_ws, itabg_ms, itabg_us, itabg_ws, sea

use mem_sflux,  only: nseaflux, nlandflux, seafluxg, landfluxg

implicit none

integer :: im,iu,iw
integer :: im1,im2,im3,iw1,iw2
integer :: igp,jgp
integer :: i, j, ii, jj, iil, jjl, iis, jjs

integer :: iter,ibin
integer :: ngroups
integer :: numtot, numcut

real :: xewm (nwa),yewm (nwa),zewm (nwa)
real :: xewml(nwl),yewml(nwl),zewml(nwl)
real :: xewms(nws),yewms(nws),zewms(nws)

real :: xmin,ymin,zmin
real :: xmax,ymax,zmax
real :: xrange, yrange, zrange
real :: cmin, cmax, cmin0, cut

integer :: igsize(mgroupsize)
integer :: nwg   (mgroupsize)
integer :: nwgl  (mgroupsize)
integer :: nwgs  (mgroupsize)

integer :: iwtemp (nwa-1), jwtemp (nwa-1)
integer :: iwltemp(nwl-1), jwltemp(nwl-1)
integer :: iwstemp(nws-1), jwstemp(nws-1)

integer :: num(1002)

real :: val (nwa)
real :: vall(nwl)
real :: vals(nws)

Type grp_var
   integer, allocatable :: iw(:)
   integer, allocatable :: iwl(:)
   integer, allocatable :: iws(:)
End type

type (grp_var) :: grp(mgroupsize)

! Allocate permanent itabg data structures

allocate (itabg_m(nma))
allocate (itabg_u(nua))
allocate (itabg_w(nwa))

if (isfcl == 1) then
   allocate (itabg_ml(nml))
   allocate (itabg_ul(nul))
   allocate (itabg_wl(nwl))

   allocate (itabg_ms(nms))
   allocate (itabg_us(nus))
   allocate (itabg_ws(nws))

   allocate (seafluxg(nseaflux))
   allocate (landfluxg(nlandflux))
endif

! Define temp variables for W

do iw = 2,nwa
   im1 = itab_w(iw)%im1
   im2 = itab_w(iw)%im2
   im3 = itab_w(iw)%im3

   xewm(iw) = max(xem(im1),xem(im2),xem(im3))
   yewm(iw) = max(yem(im1),yem(im2),yem(im3))
   zewm(iw) = max(zem(im1),zem(im2),zem(im3))
enddo

! Allocate and fill grp%iw, grp%iwl, and grp%iws for group 1

allocate (grp(1)%iw(nwa-1))

do iw = 2,nwa
   grp(1)%iw(iw-1) = iw
enddo

! Check if LEAF and SEA models are being used

if (isfcl == 1) then

! Define temp variables for WL

   do iw = 2,nwl
      xewml(iw) = -1.e9
      yewml(iw) = -1.e9
      zewml(iw) = -1.e9

      do j = 1,itab_wl(iw)%jm
         im = itab_wl(iw)%im(j)
      
         if (xewml(iw) < land%xeml(im)) xewml(iw) = land%xeml(im)
         if (yewml(iw) < land%yeml(im)) yewml(iw) = land%yeml(im)
         if (zewml(iw) < land%zeml(im)) zewml(iw) = land%zeml(im)
      enddo
   enddo

! Define temp variables for WS

   do iw = 2,nws
      xewms(iw) = -1.e9
      yewms(iw) = -1.e9
      zewms(iw) = -1.e9

      do j = 1,itab_ws(iw)%jm
         im = itab_ws(iw)%im(j)
      
         if (xewms(iw) < sea%xems(im)) xewms(iw) = sea%xems(im)
         if (yewms(iw) < sea%yems(im)) yewms(iw) = sea%yems(im)
         if (zewms(iw) < sea%zems(im)) zewms(iw) = sea%zems(im)
      enddo
   enddo

! Allocate and fill grp%iwl, and grp%iws for group 1

   allocate (grp(1)%iwl(nwl-1))
   allocate (grp(1)%iws(nws-1))

   do iw = 2,nwl
      grp(1)%iwl(iw-1) = iw
   enddo

   do iw = 2,nws
      grp(1)%iws(iw-1) = iw
   enddo

endif

ngroups = 1
jgp = 1
igsize(1) = mgroupsize

nwg (1) = nwa - 1
nwgl(1) = nwl - 1
nwgs(1) = nws - 1

do while (ngroups < mgroupsize)

   do igp = 1,ngroups

      if (igsize(igp) > 1) then

         jgp = jgp + 1
         igsize(jgp) = igsize(igp) / 2
         igsize(igp) = igsize(igp) - igsize(jgp)
         numcut = (nwg(igp) * igsize(igp)) / (igsize(igp) + igsize(jgp))

         xmin = 1.e9
         ymin = 1.e9
         zmin = 1.e9

         xmax = -1.e9
         ymax = -1.e9
         zmax = -1.e9

! Find max and min x,y,z of current group

         do i = 1,nwg(igp)
            iw = grp(igp)%iw(i)

            if (xmin > xewm(iw)) xmin = xewm(iw)
            if (ymin > yewm(iw)) ymin = yewm(iw)
            if (zmin > zewm(iw)) zmin = zewm(iw)

            if (xmax < xewm(iw)) xmax = xewm(iw)
            if (ymax < yewm(iw)) ymax = yewm(iw)
            if (zmax < zewm(iw)) zmax = zewm(iw)
         enddo

! Determine whether to cut in x, y, or z direction

         if (1.1 * (zmax - zmin) > xmax - xmin  .and.  &
             1.1 * (zmax - zmin) > ymax - ymin) then

! ATM cells - z direction

            do i = 1,nwg(igp)
               iw = grp(igp)%iw(i)
               val(iw) = zewm(iw)
            enddo

            cmin = zmin
            cmax = zmax

           if (isfcl == 1) then

! LAND cells - z direction

               do i = 1,nwgl(igp)
                  iw = grp(igp)%iwl(i)
                  vall(iw) = zewml(iw)
               enddo

! SEA cells - z direction

               do i = 1,nwgs(igp)
                  iw = grp(igp)%iws(i)
                  vals(iw) = zewms(iw)
               enddo
               
            endif

         elseif (1.1 * (xmax - xmin) > ymax - ymin) then

! ATM cells - x direction

            do i = 1,nwg(igp)
               iw = grp(igp)%iw(i)
               val(iw) = xewm(iw)
            enddo

            cmin = xmin
            cmax = xmax

            if (isfcl == 1) then

! LAND cells - x direction

               do i = 1,nwgl(igp)
                  iw = grp(igp)%iwl(i)
                  vall(iw) = xewml(iw)
               enddo

! SEA cells - x direction

               do i = 1,nwgs(igp)
                  iw = grp(igp)%iws(i)
                  vals(iw) = xewms(iw)
               enddo
               
            endif

         else

! ATM cells - y direction

            do i = 1,nwg(igp)
               iw = grp(igp)%iw(i)
               val(iw) = yewm(iw)
            enddo

            cmin = ymin
            cmax = ymax

            if (isfcl == 1) then

! LAND cells - y direction

               do i = 1,nwgl(igp)
                  iw = grp(igp)%iwl(i)
                  vall(iw) = yewml(iw)
               enddo

! SEA cells - y direction

               do i = 1,nwgs(igp)
                  iw = grp(igp)%iws(i)
                  vals(iw) = yewms(iw)
               enddo
               
            endif

         endif

! Determine cut value, iterating 3 times

         do iter = 1,3

            num(:) = 0

            do i = 1,nwg(igp)
               iw = grp(igp)%iw(i)
               if (val(iw) <= cmin) then
                  ibin = 1
               elseif (val(iw) >= cmax) then
                  ibin = 1002
               else
                  ibin = int(1000. * (val(iw) - cmin) / (cmax - cmin)) + 2
               endif
               num(ibin) = num(ibin) + 1
            enddo

! Sum number in each bin until reaching half of total

            numtot = num(1)
            ibin = 1
            do while (numtot + num(ibin+1) <= numcut)
               numtot = numtot + num(ibin+1)
               ibin = ibin + 1
            enddo
            
            cmin0 = cmin + (cmax - cmin) * .001 * (real(ibin) - 1.1)
            cmax = cmin + (cmax - cmin) * .001 * (real(ibin) +  .1)
            cmin = cmin0
         enddo
         cut = .5 * (cmin + cmax)

! Transfer a number of IW points from igp to jgp

         jj = 0
         ii = 0
         do i = 1,nwg(igp)
            iw = grp(igp)%iw(i)

            if (val(iw) > cut) then
               jj = jj + 1
               jwtemp(jj) = iw
            else
               ii = ii + 1
               iwtemp(ii) = iw
            endif
         enddo

         nwg(igp) = ii
         nwg(jgp) = jj

! Deallocate 1 old group and allocate 2 new groups

         deallocate(grp(igp)%iw)
         allocate(grp(igp)%iw(ii))
         allocate(grp(jgp)%iw(jj))

! Fill 2 new groups of ATM cellsfrom temporary arrays

         do j = 1,jj
            grp(jgp)%iw(j)=jwtemp(j)
         enddo

         do i = 1,ii
            grp(igp)%iw(i)=iwtemp(i)
         enddo

         if (isfcl == 1) then

! Transfer a number of IWL points from igp to jgp

            jjl = 0
            iil = 0
            do i = 1,nwgl(igp)
               iw = grp(igp)%iwl(i)

               if (vall(iw) > cut) then
                  jjl = jjl + 1
                  jwltemp(jjl) = iw
               else
                  iil = iil + 1
                  iwltemp(iil) = iw
               endif
            enddo

            nwgl(igp) = iil
            nwgl(jgp) = jjl

! Deallocate 1 old group and allocate 2 new groups

            deallocate(grp(igp)%iwl)
            allocate(grp(igp)%iwl(iil))
            allocate(grp(jgp)%iwl(jjl))

! Fill 2 new groups of LAND cells from temporary arrays

            do j = 1,jjl
               grp(jgp)%iwl(j)=jwltemp(j)
            enddo

            do i = 1,iil
               grp(igp)%iwl(i)=iwltemp(i)
            enddo

! Transfer a number of IWS points from igp to jgp

            jjs = 0
            iis = 0
            do i = 1,nwgs(igp)
               iw = grp(igp)%iws(i)

               if (vals(iw) > cut) then
                  jjs = jjs + 1
                  jwstemp(jjs) = iw
               else
                  iis = iis + 1
                  iwstemp(iis) = iw
               endif
            enddo

            nwgs(igp) = iis
            nwgs(jgp) = jjs

! Deallocate 1 old group and allocate 2 new groups

            deallocate(grp(igp)%iws)
            allocate(grp(igp)%iws(iis))
            allocate(grp(jgp)%iws(jjs))

! Fill 2 new groups of SEA cells from temporary arrays

            do j = 1,jjs
               grp(jgp)%iws(j)=jwstemp(j)
            enddo

            do i = 1,iis
               grp(igp)%iws(i)=iwstemp(i)
            enddo

         endif

      endif
   enddo

   ngroups = jgp

enddo

! Fill irank for each IW point from group array

do igp = 1,ngroups

! ATM cells

   do i = 1,nwg(igp)
      iw = grp(igp)%iw(i)
      itabg_w(iw)%irank = igp - 1
   enddo

   if (isfcl == 1) then

! LAND cells

      do i = 1,nwgl(igp)
         iw = grp(igp)%iwl(i)
         itabg_wl(iw)%irank = igp - 1
      enddo

! SEA cells

      do i = 1,nwgs(igp)
         iw = grp(igp)%iws(i)
         itabg_ws(iw)%irank = igp - 1
      enddo

   endif

enddo

! Loop over all U points and assign its rank to the higher rank of its
! two IW neighbors

! ATM cells

do iu = 2,nua
   iw1 = itab_u(iu)%iw1
   iw2 = itab_u(iu)%iw2

   itabg_u(iu)%irank = max(itabg_w(iw1)%irank,itabg_w(iw2)%irank)
enddo

if (isfcl == 1) then

! LAND cells

   do iu = 2,nul
      iw1 = itab_ul(iu)%iw1
      iw2 = itab_ul(iu)%iw2

! iw1 or iw2 may be zero at edge of land grid

      if (iw1 < 2) then
         itabg_ul(iu)%irank = itabg_wl(iw2)%irank
      elseif (iw2 < 2) then
         itabg_ul(iu)%irank = itabg_wl(iw1)%irank
      else
         itabg_ul(iu)%irank = max(itabg_wl(iw1)%irank,itabg_wl(iw2)%irank)
      endif
   enddo

! SEA cells

   do iu = 2,nus
      iw1 = itab_us(iu)%iw1
      iw2 = itab_us(iu)%iw2

! iw1 or iw2 may be zero at edge of sea grid

      if (iw1 < 2) then
         itabg_us(iu)%irank = itabg_ws(iw2)%irank
      elseif (iw2 < 2) then
         itabg_us(iu)%irank = itabg_ws(iw1)%irank
      else
         itabg_us(iu)%irank = max(itabg_ws(iw1)%irank,itabg_ws(iw2)%irank)
      endif
      
   enddo

endif

return
end subroutine para_decomp

!===============================================================================

subroutine para_init()

use misc_coms,  only: io6

use mem_ijtabs, only: itab_m, itab_u, itab_w, ltab_m, ltab_u, ltab_w,  &
                      itabg_m, itabg_u, itabg_w, nloops_m, nloops_u, nloops_w, &
                      alloc_itabs, mrls, maxtpn

use mem_grid,   only: nza, nma, nua, nwa, mma, mua, mwa, lpu, lcu, lpw, lsw,  &
                      xem, yem, zem, xeu, yeu, zeu, xew, yew, zew,  &
                      unx, uny, unz, utx, uty, utz, wnx, wny, wnz,  &
                      dnu, dniu, dtu, arw0, topm, glatw, glonw, glatm, glonm, &
                      aru, arw, volui, volwi, volt, volti,  &
                      alloc_xyzem, alloc_grid1, alloc_grid2

use mem_para,   only: mgroupsize, myrank, send_u, recv_u, send_w, recv_w,  &
                      nsends_u, nsends_w, nrecvs_u, nrecvs_w, send_uf, recv_uf

use mem_sflux,  only: nseaflux, mseaflux, seaflux, seaflux_temp, seafluxg,  &
                      nlandflux,  mlandflux, landflux, landflux_temp, landfluxg

use sea_coms,   only: nws

use leaf_coms,  only: nwl, isfcl

use mem_sea,    only: itabg_ws

use mem_leaf,   only: itabg_wl

implicit none

integer :: j,k
integer :: im,iu,iw
integer :: itopm,iup,iwp
integer :: im1,im2,im3
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9,iu10,iu11,iu12
integer :: iw1,iw2,iw3,iw4,iw5,iw6
integer :: isf,ilf,iws,iwl
integer :: itpn

integer :: im_myrank = 1 ! Counter for M points to be included on this rank
integer :: iu_myrank = 1 ! Counter for U points to be included on this rank
integer :: iw_myrank = 1 ! Counter for W points to be included on this rank

! Automatic arrays

logical :: myrankflag_m(nma) ! Flag for M points existing on this rank
logical :: myrankflag_u(nua) ! Flag for U points existing on this rank
logical :: myrankflag_w(nwa) ! Flag for W points existing on this rank

logical :: seaflag(nws)
logical :: landflag(nwl)

integer :: lpu_temp(nua),lcu_temp(nua)
integer :: lpw_temp(nwa),lsw_temp(nwa)

real :: topm_temp(nma), glatm_temp(nma), glonm_temp(nma)

real :: xem_temp(nma),yem_temp(nma),zem_temp(nma)
real :: xeu_temp(nua),yeu_temp(nua),zeu_temp(nua)
real :: xew_temp(nwa),yew_temp(nwa),zew_temp(nwa)
real :: unx_temp(nua),uny_temp(nua),unz_temp(nua)
real :: utx_temp(nua),uty_temp(nua),utz_temp(nua)
real :: wnx_temp(nwa),wny_temp(nwa),wnz_temp(nwa)

real :: glatw_temp(nwa),glonw_temp(nwa),arw0_temp(nwa)
real :: dtu_temp(nua),dnu_temp(nua),dniu_temp(nua)

real :: aru_temp(nza,nua),arw_temp(nza,nwa)

real :: volui_temp(nza,nua),volwi_temp(nza,nwa)

real(kind=8) :: volt_temp(nza,nwa),volti_temp(nza,nwa)

! Allocate temporary ltab data structures

allocate (ltab_m(nma))  ! Duplicate for itab_m
allocate (ltab_u(nua))  ! Duplicate for itab_u
allocate (ltab_w(nwa))  ! Duplicate for itab_w

! Allocate temporary sflux data structures

!write(io6,*) 'pin1 '

if (isfcl == 1) then
   allocate (seaflux_temp(nseaflux))    ! Duplicate for seaflux
   allocate (landflux_temp(nlandflux))  ! Duplicate for landflux
   
   seaflag(1:nws) = .false.
   landflag(1:nwl) = .false.
endif

! Allocate send & recv counter arrays and initialize to zero


!write(io6,*) 'pin2 '

allocate (nsends_u(mrls)) ; nsends_u(1:mrls) = 0
allocate (nsends_w(mrls)) ; nsends_w(1:mrls) = 0
allocate (nrecvs_u(mrls)) ; nrecvs_u(1:mrls) = 0
allocate (nrecvs_w(mrls)) ; nrecvs_w(1:mrls) = 0

! Initialize myrank flag arrays to .false.


!write(io6,*) 'pin3 '

myrankflag_m(1:nma) = .false.
myrankflag_u(1:nua) = .false.
myrankflag_w(1:nwa) = .false.

! Copy itab data structures to ltab data structures

ltab_m(1:nma) = itab_m(1:nma)
ltab_u(1:nua) = itab_u(1:nua)
ltab_w(1:nwa) = itab_w(1:nwa)

! Copy sflux data structures to sflux_temp data structures

if (isfcl == 1) then
   seaflux_temp (1:nseaflux)  = seaflux (1:nseaflux)
   landflux_temp(1:nlandflux) = landflux(1:nlandflux)
endif

! Copy grid coordinates to temporary arrays


!write(io6,*) 'pin4 '

do im = 1,nma
   topm_temp(im) = topm(im)

   xem_temp(im)  = xem(im)
   yem_temp(im)  = yem(im)
   zem_temp(im)  = zem(im)

   glatm_temp(im) = glatm(im)
   glonm_temp(im) = glonm(im)
enddo

do iu = 1,nua
   lpu_temp(iu) = lpu(iu)
   lcu_temp(iu) = lcu(iu)

   xeu_temp(iu) = xeu(iu)
   yeu_temp(iu) = yeu(iu)
   zeu_temp(iu) = zeu(iu)
   
   unx_temp(iu) = unx(iu)
   uny_temp(iu) = uny(iu)
   unz_temp(iu) = unz(iu)

   utx_temp(iu) = utx(iu)
   uty_temp(iu) = uty(iu)
   utz_temp(iu) = utz(iu)

   dtu_temp(iu)  = dtu(iu)
   dnu_temp(iu)  = dnu(iu)
   dniu_temp(iu) = dniu(iu)
   
   do k = 1,nza
      aru_temp(k,iu) = aru(k,iu)
      volui_temp(k,iu) = volui(k,iu)
   enddo
enddo

do iw = 1,nwa
   lpw_temp(iw) = lpw(iw)
   lsw_temp(iw) = lsw(iw)

   xew_temp(iw) = xew(iw)
   yew_temp(iw) = yew(iw)
   zew_temp(iw) = zew(iw)

   wnx_temp(iw) = wnx(iw)
   wny_temp(iw) = wny(iw)
   wnz_temp(iw) = wnz(iw)

   glatw_temp(iw) = glatw(iw)
   glonw_temp(iw) = glonw(iw)

   arw0_temp(iw) = arw0(iw)
   
   do k = 1,nza
      arw_temp(k,iw) = arw(k,iw)
      volwi_temp(k,iw) = volwi(k,iw)
      volt_temp(k,iw) = volt(k,iw)
      volti_temp(k,iw) = volti(k,iw)
   enddo
enddo

! Deallocate itab data structures and main grid coordinate arrays


!write(io6,*) 'pin5 '

deallocate (itab_m, itab_u, itab_w)
if (isfcl == 1) deallocate (seaflux, landflux)

deallocate (lpu, lcu, lpw, lsw)
deallocate (xem, yem, zem, topm, glatm, glonm)
deallocate (xeu, yeu, zeu, unx, uny, unz, utx, uty, utz, dtu, dnu, dniu)
deallocate (aru, volui)
deallocate (xew, yew, zew, wnx, wny, wnz, glatw, glonw, arw0)
deallocate (arw, volwi, volt, volti)

! Loop over all U points, and for each whose assigned irank is equal to myrank,
! flag all U and W points in its computational stencil for inclusion on this
! rank, excluding IUP and IWP.


!write(io6,*) 'pin6 '

do iu = 2,nua

   if (itabg_u(iu)%irank == myrank) then

      iu1  = ltab_u(iu)%iu1 
      iu2  = ltab_u(iu)%iu2 
      iu3  = ltab_u(iu)%iu3 
      iu4  = ltab_u(iu)%iu4
      iu5  = ltab_u(iu)%iu5
      iu6  = ltab_u(iu)%iu6
      iu7  = ltab_u(iu)%iu7
      iu8  = ltab_u(iu)%iu8
      iu9  = ltab_u(iu)%iu9
      iu10 = ltab_u(iu)%iu10
      iu11 = ltab_u(iu)%iu11
      iu12 = ltab_u(iu)%iu12

      iw1 = ltab_u(iu)%iw1
      iw2 = ltab_u(iu)%iw2
      iw3 = ltab_u(iu)%iw3
      iw4 = ltab_u(iu)%iw4
      iw5 = ltab_u(iu)%iw5
      iw6 = ltab_u(iu)%iw6

      myrankflag_u(iu) = .true.

      myrankflag_u(iu1) = .true.
      myrankflag_u(iu2) = .true.
      myrankflag_u(iu3) = .true.
      myrankflag_u(iu4) = .true.
      myrankflag_u(iu5) = .true.
      myrankflag_u(iu6) = .true.
      myrankflag_u(iu7) = .true.
      myrankflag_u(iu8) = .true.
      myrankflag_u(iu9) = .true.
      myrankflag_u(iu10) = .true.
      myrankflag_u(iu11) = .true.
      myrankflag_u(iu12) = .true.

      myrankflag_w(iw1) = .true.
      myrankflag_w(iw2) = .true.
      myrankflag_w(iw3) = .true.
      myrankflag_w(iw4) = .true.
      myrankflag_w(iw5) = .true.
      myrankflag_w(iw6) = .true.

   endif
enddo

! Loop over all U points, and for each that has been flagged for inclusion
! on this rank, flag both its M points for inclusion on this rank.
! Count U points also.

do iu = 2,nua

   if (myrankflag_u(iu)) then

      im1 = ltab_u(iu)%im1 
      im2 = ltab_u(iu)%im2

      myrankflag_m(im1) = .true.
      myrankflag_m(im2) = .true.

      iu_myrank = iu_myrank + 1

   endif
enddo

! Loop over all M and W points and count the ones that have been flagged
! for inclusion on this rank.

do im = 2,nma
   if (myrankflag_m(im)) then
      im_myrank = im_myrank + 1
   endif
enddo

do iw = 2,nwa
   if (myrankflag_w(iw)) then
      iw_myrank = iw_myrank + 1
   endif
enddo

! Set mma, mua, mwa values for this rank


!write(io6,*) 'pin7 '

mma = im_myrank
mua = iu_myrank
mwa = iw_myrank

! Allocate itab data structures and main grid coordinate arrays

call alloc_itabs(mma,mua,mwa)
call alloc_xyzem()
call alloc_grid1()
call alloc_grid2()

! Reset point counts to 1

im_myrank = 1
iu_myrank = 1
iw_myrank = 1

! Store new myrank M, U, W indices in itabg data structures

do im = 1,nma
   if (myrankflag_m(im)) then
      im_myrank = im_myrank + 1
      
      itabg_m(im)%im_myrank = im_myrank      
   endif
enddo

do iu = 1,nua
   if (myrankflag_u(iu)) then
      iu_myrank = iu_myrank + 1
      
      itabg_u(iu)%iu_myrank = iu_myrank
      
! Fill itabg_u(iu)%iu_myrank value for IUP point

      iup = ltab_u(iu)%iup
      itabg_u(iup)%iu_myrank = iu_myrank
   endif
enddo

do iw = 1,nwa
   if (myrankflag_w(iw)) then
      iw_myrank = iw_myrank + 1

      itabg_w(iw)%iw_myrank = iw_myrank      

! Fill itabg_w(iw)%iw_myrank value for IWP point

      iwp = ltab_w(iw)%iwp
      itabg_w(iwp)%iw_myrank = iw_myrank
   endif
enddo

! M point memory copy

do im = 1,nma
   if (myrankflag_m(im)) then
      im_myrank = itabg_m(im)%im_myrank

! First, copy entire data type from global index to subdomain index

      itab_m(im_myrank) = ltab_m(im)

! Next, redefine some individual itab_m members

      itab_m(im_myrank)%imglobe = im

! Reset IM neighbor indices to 1

      itab_m(im_myrank)%itopm = 1
      itab_m(im_myrank)%iw(1:maxtpn) = 1
      itab_m(im_myrank)%iu(1:maxtpn) = 1

! Global indices of neighbors of IM
! Set indices of neighbors of IM that are present on this rank

      itopm = ltab_m(im)%itopm
      if (myrankflag_m(itopm)) itab_m(im_myrank)%itopm = itabg_m(itopm)%im_myrank

      do itpn = 1,ltab_m(im)%ntpn
         iu = ltab_m(im)%iu(itpn)
         iw = ltab_m(im)%iw(itpn)
      
         if (myrankflag_u(iu)) itab_m(im_myrank)%iu(itpn) = itabg_u(iu)%iu_myrank
         if (myrankflag_w(iw)) itab_m(im_myrank)%iw(itpn) = itabg_w(iw)%iw_myrank
      enddo
      
! Copy M point grid values

      xem(im_myrank) = xem_temp(im)
      yem(im_myrank) = yem_temp(im)
      zem(im_myrank) = zem_temp(im)

      topm(im_myrank) = topm_temp(im)

      glatm(im_myrank) = glatm_temp(im)
      glonm(im_myrank) = glonm_temp(im)
   endif
enddo

! U point memory copy


!write(io6,*) 'pin8 '

do iu = 1,nua
   if (myrankflag_u(iu)) then
      iu_myrank = itabg_u(iu)%iu_myrank

! First, copy entire data type from global index to subdomain index

      itab_u(iu_myrank) = ltab_u(iu)

! Next, redefine some individual itab_u members

      itab_u(iu_myrank)%irank   = itabg_u(iu)%irank
      itab_u(iu_myrank)%iuglobe = iu

! Reset some loop flags if this IU point is primary on a remote rank

      if (itab_u(iu_myrank)%irank /= myrank) then      

         call uloops('n',iu_myrank,-4,-7,-8,-12,-13,-15,-16, 0, 0, 0)

! Turn off LBC copy if IUP point is on remote node

         if (.not. myrankflag_u(iup))  &
            call uloops('n',iu_myrank,-2,-9,-18,0,0,0,0,0,0,0)

      endif

! Reset IU neighbor indices to 1

      itab_u(iu_myrank)%iup = 1

      itab_u(iu_myrank)%im1 = 1
      itab_u(iu_myrank)%im2 = 1

      itab_u(iu_myrank)%iu1 = 1
      itab_u(iu_myrank)%iu2 = 1
      itab_u(iu_myrank)%iu3 = 1
      itab_u(iu_myrank)%iu4 = 1
      itab_u(iu_myrank)%iu5 = 1
      itab_u(iu_myrank)%iu6 = 1
      itab_u(iu_myrank)%iu7 = 1
      itab_u(iu_myrank)%iu8 = 1
      itab_u(iu_myrank)%iu9 = 1
      itab_u(iu_myrank)%iu10 = 1
      itab_u(iu_myrank)%iu11 = 1
      itab_u(iu_myrank)%iu12 = 1

      itab_u(iu_myrank)%iw1 = 1
      itab_u(iu_myrank)%iw2 = 1
      itab_u(iu_myrank)%iw3 = 1
      itab_u(iu_myrank)%iw4 = 1
      itab_u(iu_myrank)%iw5 = 1
      itab_u(iu_myrank)%iw6 = 1

! Global indices of neighbors of IU

      iup = ltab_u(iu)%iup

      im1 = ltab_u(iu)%im1
      im2 = ltab_u(iu)%im2

      iu1  = ltab_u(iu)%iu1
      iu2  = ltab_u(iu)%iu2
      iu3  = ltab_u(iu)%iu3
      iu4  = ltab_u(iu)%iu4
      iu5  = ltab_u(iu)%iu5
      iu6  = ltab_u(iu)%iu6
      iu7  = ltab_u(iu)%iu7
      iu8  = ltab_u(iu)%iu8
      iu9  = ltab_u(iu)%iu9
      iu10 = ltab_u(iu)%iu10
      iu11 = ltab_u(iu)%iu11
      iu12 = ltab_u(iu)%iu12

      iw1 = ltab_u(iu)%iw1
      iw2 = ltab_u(iu)%iw2
      iw3 = ltab_u(iu)%iw3
      iw4 = ltab_u(iu)%iw4
      iw5 = ltab_u(iu)%iw5
      iw6 = ltab_u(iu)%iw6

! Set indices of neighbors of IU that are present on this rank

      if (myrankflag_u(iup)) itab_u(iu_myrank)%iup = itabg_u(iup)%iu_myrank

      if (myrankflag_m(im1)) itab_u(iu_myrank)%im1 = itabg_m(im1)%im_myrank
      if (myrankflag_m(im2)) itab_u(iu_myrank)%im2 = itabg_m(im2)%im_myrank

      if (myrankflag_u(iu1))  itab_u(iu_myrank)%iu1  = itabg_u(iu1)%iu_myrank
      if (myrankflag_u(iu2))  itab_u(iu_myrank)%iu2  = itabg_u(iu2)%iu_myrank
      if (myrankflag_u(iu3))  itab_u(iu_myrank)%iu3  = itabg_u(iu3)%iu_myrank
      if (myrankflag_u(iu4))  itab_u(iu_myrank)%iu4  = itabg_u(iu4)%iu_myrank
      if (myrankflag_u(iu5))  itab_u(iu_myrank)%iu5  = itabg_u(iu5)%iu_myrank
      if (myrankflag_u(iu6))  itab_u(iu_myrank)%iu6  = itabg_u(iu6)%iu_myrank
      if (myrankflag_u(iu7))  itab_u(iu_myrank)%iu7  = itabg_u(iu7)%iu_myrank
      if (myrankflag_u(iu8))  itab_u(iu_myrank)%iu8  = itabg_u(iu8)%iu_myrank
      if (myrankflag_u(iu9))  itab_u(iu_myrank)%iu9  = itabg_u(iu9)%iu_myrank
      if (myrankflag_u(iu10)) itab_u(iu_myrank)%iu10 = itabg_u(iu10)%iu_myrank
      if (myrankflag_u(iu11)) itab_u(iu_myrank)%iu11 = itabg_u(iu11)%iu_myrank
      if (myrankflag_u(iu12)) itab_u(iu_myrank)%iu12 = itabg_u(iu12)%iu_myrank

      if (myrankflag_w(iw1)) itab_u(iu_myrank)%iw1 = itabg_w(iw1)%iw_myrank
      if (myrankflag_w(iw2)) itab_u(iu_myrank)%iw2 = itabg_w(iw2)%iw_myrank
      if (myrankflag_w(iw3)) itab_u(iu_myrank)%iw3 = itabg_w(iw3)%iw_myrank
      if (myrankflag_w(iw4)) itab_u(iu_myrank)%iw4 = itabg_w(iw4)%iw_myrank
      if (myrankflag_w(iw5)) itab_u(iu_myrank)%iw5 = itabg_w(iw5)%iw_myrank
      if (myrankflag_w(iw6)) itab_u(iu_myrank)%iw6 = itabg_w(iw6)%iw_myrank

! Copy U point grid values

      lpu(iu_myrank) = lpu_temp(iu)
      lcu(iu_myrank) = lcu_temp(iu)

      xeu(iu_myrank) = xeu_temp(iu)
      yeu(iu_myrank) = yeu_temp(iu)
      zeu(iu_myrank) = zeu_temp(iu)

      unx(iu_myrank) = unx_temp(iu)
      uny(iu_myrank) = uny_temp(iu)
      unz(iu_myrank) = unz_temp(iu)

      utx(iu_myrank) = utx_temp(iu)
      uty(iu_myrank) = uty_temp(iu)
      utz(iu_myrank) = utz_temp(iu)

      dtu(iu_myrank) = dtu_temp(iu)
      dnu(iu_myrank) = dnu_temp(iu)
      dniu(iu_myrank) = dniu_temp(iu)
      
      do k = 1,nza
         aru  (k,iu_myrank) = aru_temp (k,iu)
         volui(k,iu_myrank) = volui_temp(k,iu)
      enddo
   endif
enddo


!write(io6,*) 'pin9 '

! W point memory copy

do iw = 1,nwa
   if (myrankflag_w(iw)) then
      iw_myrank = itabg_w(iw)%iw_myrank

! First, copy entire data type from global index to subdomain index

      itab_w(iw_myrank) = ltab_w(iw)

! Next, redefine some individual itab_w members

      itab_w(iw_myrank)%irank   = itabg_w(iw)%irank
      itab_w(iw_myrank)%iwglobe = iw

! Reset some loop flags if this IW point is primary on a remote rank

      if (itab_w(iw_myrank)%irank /= myrank) then      

         call wloops('n',iw_myrank, -5,-12,-13,-15,-16,-17,-19,-20,-23,-26)
         call wloops('n',iw_myrank,-27,-29,-30,-34,  0,  0,  0,  0,  0,  0)

! Turn off LBC copy if IWP point is on remote node

         if (.not. myrankflag_w(iwp))  &
            call wloops('n',iw_myrank,-2,-4,-9,-22,-24,-31,-35,0,0,0)

      endif

! Reset IW neighbor indices to 1

      itab_w(iw_myrank)%iwp = 1

      itab_w(iw_myrank)%im1 = 1
      itab_w(iw_myrank)%im2 = 1
      itab_w(iw_myrank)%im3 = 1

      itab_w(iw_myrank)%iu1 = 1
      itab_w(iw_myrank)%iu2 = 1
      itab_w(iw_myrank)%iu3 = 1
      itab_w(iw_myrank)%iu4 = 1
      itab_w(iw_myrank)%iu5 = 1
      itab_w(iw_myrank)%iu6 = 1
      itab_w(iw_myrank)%iu7 = 1
      itab_w(iw_myrank)%iu8 = 1
      itab_w(iw_myrank)%iu9 = 1

      itab_w(iw_myrank)%iw1 = 1
      itab_w(iw_myrank)%iw2 = 1
      itab_w(iw_myrank)%iw3 = 1

! Global indices of neighbors of IW

      iwp = ltab_w(iw)%iwp

      im1 = ltab_w(iw)%im1
      im2 = ltab_w(iw)%im2
      im3 = ltab_w(iw)%im3

      iu1 = ltab_w(iw)%iu1
      iu2 = ltab_w(iw)%iu2
      iu3 = ltab_w(iw)%iu3
      iu4 = ltab_w(iw)%iu4
      iu5 = ltab_w(iw)%iu5
      iu6 = ltab_w(iw)%iu6
      iu7 = ltab_w(iw)%iu7
      iu8 = ltab_w(iw)%iu8
      iu9 = ltab_w(iw)%iu9

      iw1 = ltab_w(iw)%iw1
      iw2 = ltab_w(iw)%iw2
      iw3 = ltab_w(iw)%iw3

! Set indices of neighbors of IW that are present on this rank

      if (myrankflag_w(iwp)) itab_w(iw_myrank)%iwp = itabg_w(iwp)%iw_myrank

      if (myrankflag_m(im1)) itab_w(iw_myrank)%im1 = itabg_m(im1)%im_myrank
      if (myrankflag_m(im2)) itab_w(iw_myrank)%im2 = itabg_m(im2)%im_myrank
      if (myrankflag_m(im3)) itab_w(iw_myrank)%im3 = itabg_m(im3)%im_myrank

      if (myrankflag_u(iu1)) itab_w(iw_myrank)%iu1 = itabg_u(iu1)%iu_myrank
      if (myrankflag_u(iu2)) itab_w(iw_myrank)%iu2 = itabg_u(iu2)%iu_myrank
      if (myrankflag_u(iu3)) itab_w(iw_myrank)%iu3 = itabg_u(iu3)%iu_myrank
      if (myrankflag_u(iu4)) itab_w(iw_myrank)%iu4 = itabg_u(iu4)%iu_myrank
      if (myrankflag_u(iu5)) itab_w(iw_myrank)%iu5 = itabg_u(iu5)%iu_myrank
      if (myrankflag_u(iu6)) itab_w(iw_myrank)%iu6 = itabg_u(iu6)%iu_myrank
      if (myrankflag_u(iu7)) itab_w(iw_myrank)%iu7 = itabg_u(iu7)%iu_myrank
      if (myrankflag_u(iu8)) itab_w(iw_myrank)%iu8 = itabg_u(iu8)%iu_myrank
      if (myrankflag_u(iu9)) itab_w(iw_myrank)%iu9 = itabg_u(iu9)%iu_myrank

      if (myrankflag_w(iw1)) itab_w(iw_myrank)%iw1 = itabg_w(iw1)%iw_myrank
      if (myrankflag_w(iw2)) itab_w(iw_myrank)%iw2 = itabg_w(iw2)%iw_myrank
      if (myrankflag_w(iw3)) itab_w(iw_myrank)%iw3 = itabg_w(iw3)%iw_myrank

! Copy W point grid values

      lpw(iw_myrank) = lpw_temp(iw)
      lsw(iw_myrank) = lsw_temp(iw)

      xew(iw_myrank) = xew_temp(iw)
      yew(iw_myrank) = yew_temp(iw)
      zew(iw_myrank) = zew_temp(iw)

      wnx(iw_myrank) = wnx_temp(iw)
      wny(iw_myrank) = wny_temp(iw)
      wnz(iw_myrank) = wnz_temp(iw)

      glatw(iw_myrank) = glatw_temp(iw)
      glonw(iw_myrank) = glonw_temp(iw)

      arw0(iw_myrank) = arw0_temp(iw)

      do k = 1,nza
         arw  (k,iw_myrank) = arw_temp  (k,iw)
         volwi(k,iw_myrank) = volwi_temp(k,iw)
         volt (k,iw_myrank) = volt_temp (k,iw)
         volti(k,iw_myrank) = volti_temp(k,iw)
      enddo
      
   endif
enddo

! Check whether LAND/SEA models are used


!write(io6,*) 'pin10 '

if (isfcl == 1) then

! Copy SEAFLUX values


!write(io6,*) 'pin1 '

   mseaflux = 1

   do isf = 2,nseaflux
      iw  = seaflux_temp(isf)%iw
      iws = seaflux_temp(isf)%iws

      if (itabg_w (iw )%irank == myrank .or.  &
          itabg_ws(iws)%irank == myrank) then

         mseaflux = mseaflux + 1
         seaflag(iws) = .true.
      endif
   enddo


!write(io6,*) 'pin11 '

   allocate (seaflux(mseaflux))


!write(io6,*) 'pin12 '

   mseaflux = 1

   seaflux(1)%isfglobe = 1
   seaflux(1)%iw = 1
   seaflux(1)%iws = 1

   do isf = 2,nseaflux
      iw  = seaflux_temp(isf)%iw   ! full-domain index
      iws = seaflux_temp(isf)%iws  ! full-domain index

      if (itabg_w (iw )%irank == myrank .or.  &
          itabg_ws(iws)%irank == myrank) then

         mseaflux = mseaflux + 1

         seaflux(mseaflux) = seaflux_temp(isf) ! retain global indices

         seafluxg(isf)%isf_myrank = mseaflux
      endif
   enddo

!write(io6,*) 'pin13 ',mgroupsize,mseaflux

! Allocate seaflux%send_sf AFTER the preceding copy from seaflux_temp to seaflux
! because seaflux_temp%send_sf is not (and does not need to be) allocated.

   do isf = 1,mseaflux
      iw = seaflux(isf)%iw
      seaflux(isf)%iwrank = itabg_w(iw)%irank
      allocate (seaflux(isf)%send_sf(mgroupsize))
      seaflux(isf)%send_sf(1:mgroupsize) = .false.
   enddo

! Copy LANDLUX values

!write(io6,*) 'pin14 ',allocated(seaflux(303)%send_sf)

   mlandflux = 1

   do ilf = 2,nlandflux
      iw  = landflux_temp(ilf)%iw
      iwl = landflux_temp(ilf)%iwl

      if (itabg_w (iw )%irank == myrank .or.  &
          itabg_wl(iwl)%irank == myrank) then

         mlandflux = mlandflux + 1
         landflag(iwl) = .true.
      endif
   enddo


!write(io6,*) 'pin15 ',allocated(seaflux(303)%send_sf)

   allocate (landflux(mlandflux))


!write(io6,*) 'pin16 ',allocated(seaflux(303)%send_sf)


   mlandflux = 1

   landflux(1)%ilfglobe = 1
   landflux(1)%iw = 1
   landflux(1)%iwl = 1

   do ilf = 2,nlandflux
      iw  = landflux_temp(ilf)%iw   ! full-domain index
      iwl = landflux_temp(ilf)%iwl  ! full-domain index

      if (itabg_w (iw )%irank == myrank .or.  &
          itabg_wl(iwl)%irank == myrank) then

         mlandflux = mlandflux + 1

         landflux(mlandflux) = landflux_temp(ilf) ! retain global indices

         landfluxg(ilf)%ilf_myrank = mlandflux
      endif
   enddo

!write(io6,*) 'pin17 ',allocated(seaflux(303)%send_sf)

! Allocate landflux%send_lf AFTER the preceding copy from landflux_temp to landflux
! because landflux_temp%send_lf is not (and does not need to be) allocated.

   do ilf = 1,mlandflux
      iw = landflux(ilf)%iw
      landflux(ilf)%iwrank = itabg_w(iw)%irank
      allocate (landflux(ilf)%send_lf(mgroupsize))
      landflux(ilf)%send_lf(1:mgroupsize) = .false.
   enddo


!write(io6,*) 'pin18 ',allocated(seaflux(303)%send_sf)

   call para_init_sea(seaflag)

!write(io6,*) 'pin18.5 '

   call para_init_land(landflag)


!write(io6,*) 'pin19 '

endif

! Initialize iremote_u and iremote_w values to 0.


!write(io6,*) 'pin20 '

do j = 1,mgroupsize
   send_u(j)%iremote = 0
   send_w(j)%iremote = 0
   send_uf(j)%iremote = 0

   recv_u(j)%iremote = 0
   recv_w(j)%iremote = 0
   recv_uf(j)%iremote = 0
enddo

! Loop over all U points and for each that is primary on a remote rank, 
! access all U and W points in its stencil.

do iu = 2,nua
   if (itabg_u(iu)%irank /= myrank) then
   
      iup  = ltab_u(iu)%iup 

      iu1  = ltab_u(iu)%iu1 
      iu2  = ltab_u(iu)%iu2 
      iu3  = ltab_u(iu)%iu3 
      iu4  = ltab_u(iu)%iu4
      iu5  = ltab_u(iu)%iu5
      iu6  = ltab_u(iu)%iu6
      iu7  = ltab_u(iu)%iu7
      iu8  = ltab_u(iu)%iu8
      iu9  = ltab_u(iu)%iu9
      iu10 = ltab_u(iu)%iu10
      iu11 = ltab_u(iu)%iu11
      iu12 = ltab_u(iu)%iu12

      iw1 = ltab_u(iu)%iw1
      iw2 = ltab_u(iu)%iw2
      iw3 = ltab_u(iu)%iw3
      iw4 = ltab_u(iu)%iw4
      iw5 = ltab_u(iu)%iw5
      iw6 = ltab_u(iu)%iw6

! IU point is primary on remote rank.  

! If IU point is in memory of myrank, its value must be received from 
! remote rank.  Add that remote rank to receive table.

      if (myrankflag_u(iu)) call recv_table_u(itabg_u(iu)%irank)

! Add to send table any U or W point in the stencil of IU that is primary
! on myrank.  

      if (itabg_u(iup )%irank == myrank) call send_table_u(iup ,itabg_u(iu)%irank)

      if (itabg_u(iu1 )%irank == myrank) call send_table_u(iu1 ,itabg_u(iu)%irank)
      if (itabg_u(iu2 )%irank == myrank) call send_table_u(iu2 ,itabg_u(iu)%irank)
      if (itabg_u(iu3 )%irank == myrank) call send_table_u(iu3 ,itabg_u(iu)%irank)
      if (itabg_u(iu4 )%irank == myrank) call send_table_u(iu4 ,itabg_u(iu)%irank)
      if (itabg_u(iu5 )%irank == myrank) call send_table_u(iu5 ,itabg_u(iu)%irank)
      if (itabg_u(iu6 )%irank == myrank) call send_table_u(iu6 ,itabg_u(iu)%irank)
      if (itabg_u(iu7 )%irank == myrank) call send_table_u(iu7 ,itabg_u(iu)%irank)
      if (itabg_u(iu8 )%irank == myrank) call send_table_u(iu8 ,itabg_u(iu)%irank)
      if (itabg_u(iu9 )%irank == myrank) call send_table_u(iu9 ,itabg_u(iu)%irank)
      if (itabg_u(iu10)%irank == myrank) call send_table_u(iu10,itabg_u(iu)%irank)
      if (itabg_u(iu11)%irank == myrank) call send_table_u(iu11,itabg_u(iu)%irank)
      if (itabg_u(iu12)%irank == myrank) call send_table_u(iu12,itabg_u(iu)%irank)

      if (itabg_w(iw1 )%irank == myrank) call send_table_w(iw1,itabg_u(iu)%irank)
      if (itabg_w(iw2 )%irank == myrank) call send_table_w(iw2,itabg_u(iu)%irank)
      if (itabg_w(iw3 )%irank == myrank) call send_table_w(iw3,itabg_u(iu)%irank)
      if (itabg_w(iw4 )%irank == myrank) call send_table_w(iw4,itabg_u(iu)%irank)
      if (itabg_w(iw5 )%irank == myrank) call send_table_w(iw5,itabg_u(iu)%irank)
      if (itabg_w(iw6 )%irank == myrank) call send_table_w(iw6,itabg_u(iu)%irank)

   endif
enddo


!write(io6,*) 'pin21 '

! Loop over all W points and for each that is primary on a remote rank,
! access all U and W points in its stencil.

do iw = 2,nwa
   if (itabg_w(iw)%irank /= myrank) then

      iu1  = ltab_w(iw)%iu1 
      iu2  = ltab_w(iw)%iu2 
      iu3  = ltab_w(iw)%iu3 
      iu4  = ltab_w(iw)%iu4
      iu5  = ltab_w(iw)%iu5
      iu6  = ltab_w(iw)%iu6
      iu7  = ltab_w(iw)%iu7
      iu8  = ltab_w(iw)%iu8
      iu9  = ltab_w(iw)%iu9

      iwp = ltab_w(iw)%iwp

      iw1 = ltab_w(iw)%iw1
      iw2 = ltab_w(iw)%iw2
      iw3 = ltab_w(iw)%iw3

! If IW point is in memory of myrank, it must be received from remote rank.
! Add that remote rank to receive table.

      if (myrankflag_w(iw)) call recv_table_w(itabg_w(iw)%irank)

! Add to send table any U or W point in the stencil of IW that is primary
! on myrank.  

      if (itabg_u(iu1)%irank == myrank) call send_table_u(iu1,itabg_w(iw)%irank)
      if (itabg_u(iu2)%irank == myrank) call send_table_u(iu2,itabg_w(iw)%irank)
      if (itabg_u(iu3)%irank == myrank) call send_table_u(iu3,itabg_w(iw)%irank)
      if (itabg_u(iu4)%irank == myrank) call send_table_u(iu4,itabg_w(iw)%irank)
      if (itabg_u(iu5)%irank == myrank) call send_table_u(iu5,itabg_w(iw)%irank)
      if (itabg_u(iu6)%irank == myrank) call send_table_u(iu6,itabg_w(iw)%irank)
      if (itabg_u(iu7)%irank == myrank) call send_table_u(iu7,itabg_w(iw)%irank)
      if (itabg_u(iu8)%irank == myrank) call send_table_u(iu8,itabg_w(iw)%irank)
      if (itabg_u(iu9)%irank == myrank) call send_table_u(iu9,itabg_w(iw)%irank)

      if (itabg_w(iwp)%irank == myrank) call send_table_w(iwp,itabg_w(iw)%irank)

      if (itabg_w(iw1)%irank == myrank) call send_table_w(iw1,itabg_w(iw)%irank)
      if (itabg_w(iw2)%irank == myrank) call send_table_w(iw2,itabg_w(iw)%irank)
      if (itabg_w(iw3)%irank == myrank) call send_table_w(iw3,itabg_w(iw)%irank)

   endif
enddo


!write(io6,*) 'pin22 '

! Deallocate temporary data structures and arrays

deallocate (ltab_m,ltab_u,ltab_w)

return
end subroutine para_init

!===============================================================================

subroutine recv_table_u(iremote)

use mem_para,  only: nrecvs_u, recv_u, recv_uf
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote is already in table of ranks to receive from

do jrecv=1,nrecvs_u(1)
   if (recv_u(jrecv)%iremote == iremote) exit
enddo

! If jrecv exceeds nrecvs_u(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_u(1).

if (jrecv > nrecvs_u(1)) nrecvs_u(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_u(jrecv)%iremote = iremote
recv_uf(jrecv)%iremote = iremote

return
end subroutine recv_table_u

!===============================================================================

subroutine recv_table_w(iremote)

use mem_para,  only: nrecvs_w, recv_w
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote_w is already in table of ranks to receive from

do jrecv = 1,nrecvs_w(1)
   if (recv_w(jrecv)%iremote == iremote) exit
enddo

! If jrecv exceeds nrecvs_w(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_w(1).

if (jrecv > nrecvs_w(1)) nrecvs_w(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_w(jrecv)%iremote = iremote

return
end subroutine recv_table_w

!===============================================================================

subroutine send_table_u(iu,iremote)

use mem_ijtabs, only: itab_u, itabg_u
use mem_para,   only: nsends_u, send_u, send_uf, mgroupsize
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iu
integer, intent(in) :: iremote

integer :: jsend
integer :: iu_myrank

! Check whether iremote_u is already in table of ranks to send to

do jsend=1,nsends_u(1)
   if (send_u(jsend)%iremote == iremote) exit
enddo

! If jsend exceeds nsends_u, jsend represents a rank not yet entered in the
! table, so increase nsends_u.

if (jsend > nsends_u(1)) nsends_u(1) = jsend

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

iu_myrank = itabg_u(iu)%iu_myrank

itab_u(iu_myrank)%loop(25+jsend) = .true.
send_u(jsend)%iremote = iremote
send_uf(jsend)%iremote = iremote

return
end subroutine send_table_u

!===============================================================================

subroutine send_table_w(iw,iremote)

use mem_ijtabs, only: itab_w, itabg_w
use mem_para,   only: nsends_w, send_w
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iw
integer, intent(in) :: iremote

integer :: jsend
integer :: iw_myrank

! Check whether iremote_w is already in table of ranks to send to

do jsend=1,nsends_w(1)
   if (send_w(jsend)%iremote == iremote) exit
enddo

! If jsend exceeds nsends_w, jsend represents a rank not yet entered in the
! table, so increase nsends_w.

if (jsend > nsends_w(1)) nsends_w = jsend

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

iw_myrank = itabg_w(iw)%iw_myrank

itab_w(iw_myrank)%loop(35+jsend) = .true.
send_w(jsend)%iremote = iremote

return
end subroutine send_table_w

