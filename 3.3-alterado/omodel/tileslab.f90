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
subroutine tileslab_hmp(iplt,action)

use oplot_coms, only: op
use mem_grid,   only: mza, mma, zm, zt, xem, yem, zem, xew, yew, zew, lpw
use mem_ijtabs, only: itab_m
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: ntpn
integer :: kt, k, ko
integer :: im
integer :: itpn
integer :: ici
integer :: ng
integer :: iw
integer :: notavail

real :: hpt
real :: vpt
real :: fldval

real :: htpn(8)
real :: vtpn(8)

! Set integer ifill flag and color table value for underground cells

! Find KT or KM level to plot

kt = 2
do while (kt < mza .and. zm(kt) < op%slabloc(iplt))
   kt = kt + 1
enddo
k = kt

! If plotting M point, set k for correct M level

if (op%stagpt == 'M' .and. zt(kt) > op%slabloc(iplt)) k = kt - 1

! If field is 3d, first plot underground points with underground color

if (action == 'T' .and. op%dimens == '3') then
   call plot_underground_m(iplt,kt)
endif

do im = 2,mma

! Check for points to be skipped over

   if (.not. itab_m(im)%loop(1)) cycle  ! For now, skip pts that don't read in topm

! Get tile plot coordinates.

   call oplot_transform(iplt,xem(im),yem(im),zem(im),hpt,vpt)

   ntpn = itab_m(im)%ntpn

   do itpn = 1,ntpn

      iw = itab_m(im)%iw(itpn)
      
! Skip this M point if current IW point index < 2 (which occurs at lateral boundary
! of limited-area domain or parallel subdomain)      

      if (iw < 2) go to 9
      
      call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),htpn(itpn),vtpn(itpn))

! Avoid wrap-around for lat-lon plot

      if (op%projectn(iplt) == 'L') call ll_unwrap(hpt,htpn(itpn))

   enddo

! Jump out of loop if any cell corner is on other side of earth

   if (any(htpn(1:ntpn) > 1.e11)) cycle

! Jump out of loop if entire cell is outside plot window. 
   
   if ( all(htpn(1:ntpn) < op%xmin) .or. all(htpn(1:ntpn) > op%xmax) .or.  &
        all(vtpn(1:ntpn) < op%ymin) .or. all(vtpn(1:ntpn) > op%ymax) ) cycle

! Jump out of loop if field is 3-D and any surrounding T cell is underground

   if (op%dimens == '3' .and. any(kt < lpw(itab_m(im)%iw(1:ntpn)))) cycle

! If field is 3-D plot point that is not adjacent to any underground T cells,
! or if field is 2-D, plot each point

   call oplot_lib(iplt,k,im,'VALUE',op%fldname(iplt),fldval,notavail,ko)
   if (notavail > 0) cycle 
   
   call celltile(iplt,im,ntpn,htpn,vtpn,hpt,vpt,fldval,action)
   	
   9 continue     
        
enddo

return
end subroutine tileslab_hmp

!===============================================================================

subroutine tileslab_htw(iplt,action)

use oplot_coms, only: op, xepc, yepc, zepc
use mem_grid,   only: mza, mwa, zm, zt, xew, yew, zew, xem, yem, zem, lpw
use mem_ijtabs, only: itab_w
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: ntpn
integer :: kt, k, ko
integer :: iw
integer :: itpn
integer :: ici
integer :: ng
integer :: im1
integer :: im2
integer :: im3
integer :: iok
integer :: notavail

real :: hpt
real :: vpt
real :: fldval

real :: htpn(3)
real :: vtpn(3)
real :: topo1, topo2

! Set integer ifill flag and color table value for underground cells

! Find KT or KM level to plot (not used if op%pltlev = 'p' or 's')

kt = 2
do while (kt < mza .and. zm(kt) < op%slabloc(iplt))
   kt = kt + 1
enddo
k = kt

! If plotting W point, set k for correct W level

if (op%stagpt == 'W' .and. zt(kt) > op%slabloc(iplt)) k = kt - 1

do iw = 2,mwa

! Get tile plot coordinates.  

   call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),hpt,vpt)

   im1 = itab_w(iw)%im1
   im2 = itab_w(iw)%im2
   im3 = itab_w(iw)%im3

   call oplot_transform(iplt,xem(im1),yem(im1),zem(im1),htpn(1),vtpn(1))
   call oplot_transform(iplt,xem(im2),yem(im2),zem(im2),htpn(2),vtpn(2))
   call oplot_transform(iplt,xem(im3),yem(im3),zem(im3),htpn(3),vtpn(3))

! Avoid wrap-around for lat-lon plot

   if (op%projectn(iplt) == 'L') then
      do itpn = 1,3
         call ll_unwrap(hpt,htpn(itpn))
      enddo
   endif

! Jump out of loop if any cell corner is on other side of earth

   if (any(htpn(1:3) > 1.e11)) cycle

! Jump out of loop if entire cell is outside plot window

   if ( all(htpn(1:3) < op%xmin) .or. all(htpn(1:3) > op%xmax) .or.  &
        all(vtpn(1:3) < op%ymin) .or. all(vtpn(1:3) > op%ymax) ) cycle

! Get cell value and plot if 'available'

   call oplot_lib(iplt,k,iw,'VALUE',op%fldname(iplt),fldval,notavail,ko)

!print*, 'ts1 ',op%fldname(iplt),notavail,op%dimens

   if (notavail > 0) then
      if (notavail == 1 .and. action == 'T') then

! Tile-plot cell with underground color

         call fillpolyg(3,htpn,vtpn,op%icigrnd)
      endif

      cycle
   endif
   
   call celltile(iplt,iw,3,htpn,vtpn,hpt,vpt,fldval,action)

! Plot cone circle if so specified

   if (op%pltcone(iplt) == 'C') then
   
      call coneplot_w(iw,topo1,topo2,iok,htpn)
      
      if (iok == 1) then
      
         call oplot_transform(iplt,xepc(1),yepc(1),zepc(1),htpn(1),vtpn(1))
         call oplot_transform(iplt,xepc(2),yepc(2),zepc(2),htpn(2),vtpn(2))
   
         call o_frstpt(htpn(1),vtpn(1))
         call o_vector(htpn(2),vtpn(2))

      endif
 
   endif
      
enddo

return
end subroutine tileslab_htw

!===============================================================================

subroutine tileslab_hu(iplt,action)

use oplot_coms, only: op
use mem_grid,   only: mza, mua, zm, xeu, yeu, zeu, xem, yem, zem,  &
                    xew, yew, zew, lpw
use mem_ijtabs, only: itab_u
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: kt, k, ko
integer :: iu
integer :: iw1
integer :: iw2
integer :: itpn
integer :: ici
integer :: ng
integer :: im1
integer :: im2
integer :: notavail

real :: fldval
real :: hpt
real :: vpt
real :: x1
real :: x4
real :: x2
real :: x3

real :: htpn(4)
real :: vtpn(4)

real :: htpn2,htpn4
real :: vtpn2,vtpn4

! Set integer ifill flag and color table value for underground cells

! Find KT or KM level to plot (not used if op%pltlev = 'p' or 's')

kt = 2
do while (kt < mza .and. zm(kt) < op%slabloc(iplt))
   kt = kt + 1
enddo
k = kt

do iu = 2,mua

! Transform tile plot X and Y coordinates.  

   call oplot_transform(iplt,xeu(iu),yeu(iu),zeu(iu),hpt,vpt)

   im1 = itab_u(iu)%im1
   im2 = itab_u(iu)%im2
   iw1 = itab_u(iu)%iw1
   iw2 = itab_u(iu)%iw2

   call oplot_transform(iplt,xem(im1),yem(im1),zem(im1),htpn(1),vtpn(1))
   call oplot_transform(iplt,xem(im2),yem(im2),zem(im2),htpn(3),vtpn(3))

   if (iw1 > 1) then
      call oplot_transform(iplt,xew(iw1),yew(iw1),zew(iw1),htpn(2),vtpn(2))
   else
      htpn(2) = htpn(3)
      vtpn(2) = vtpn(3)
   endif

   if (iw2 > 1) then
      call oplot_transform(iplt,xew(iw2),yew(iw2),zew(iw2),htpn(4),vtpn(4))
   else
      htpn(4) = htpn(3)
      vtpn(4) = vtpn(3)
   endif

!  Avoid wrap-around for lat-lon plot

   if (op%projectn(iplt) == 'L') then
      do itpn = 1,4
         call ll_unwrap(hpt,htpn(itpn))
      enddo
   endif

! Jump out of loop if any cell corner is on other side of earth

   if (any(htpn(1:4) > 1.e11)) cycle

! Jump out of loop if entire cell is outside plot window. 

   if (all(htpn(1:4) < op%xmin) .or. all(htpn(1:4) > op%xmax) .or.  &
       all(vtpn(1:4) < op%ymin) .or. all(vtpn(1:4) > op%ymax)) cycle

! Save copy of 2nd and 4th htpn,vtpn pts

   htpn2 = htpn(2)
   htpn4 = htpn(4)
   vtpn2 = vtpn(2)
   vtpn4 = vtpn(4)
   
! Get cell value and 'available' flag

   call oplot_lib(iplt,k,iu,'VALUE',op%fldname(iplt),fldval,notavail,ko)    

! Check if both IW1 and IW2 are above ground

   if (op%dimens == '2' .or.  &
      (ko >= lpw(iw1) .and. ko >= lpw(iw2))) then

! Plot full cell value unless not available

      if (notavail > 0) cycle 
      call celltile(iplt,iu,4,htpn,vtpn,hpt,vpt,fldval,action)
      
! Else, check if IW1 is above ground

   elseif (ko >= lpw(iw1)) then

! Tile-plot IW2 half of cell with underground color

      if (action == 'T') then
         htpn(2) = htpn(3)
         vtpn(2) = vtpn(3)

         call fillpolyg(4,htpn,vtpn,op%icigrnd)
      endif

! Plot IW1 half of cell with tile color unless not available

      htpn(2) = htpn2
      vtpn(2) = vtpn2

      htpn(4) = htpn(3)
      vtpn(4) = vtpn(3)

      if (notavail > 0) cycle 
      call celltile(iplt,iu,4,htpn,vtpn,hpt,vpt,fldval,action)
      
! Else, check if IW2 is above ground

   elseif (ko >= lpw(iw2)) then

! Tile-plot IW1 half of cell with underground color

      if (action == 'T') then
         htpn(4) = htpn(3)
         vtpn(4) = vtpn(3)

         call fillpolyg(4,htpn,vtpn,op%icigrnd)
      endif

! Plot IW2 half of cell with tile color unless not available

      htpn(4) = htpn2
      vtpn(4) = vtpn2

      htpn(2) = htpn(3)
      vtpn(2) = vtpn(3)

      if (notavail > 0) cycle 
      call celltile(iplt,iu,4,htpn,vtpn,hpt,vpt,fldval,action)
      
! Else, tile-plot both full IU cell with underground color

   else

! Tile-plot cell with underground color

      call fillpolyg(4,htpn,vtpn,op%icigrnd)
   endif

enddo

return
end subroutine tileslab_hu

!===============================================================================

subroutine tileslab_hs(iplt,action)

use oplot_coms, only: op
use mem_sea,    only: sea, itab_ws
use sea_coms,   only: mws, maxjms
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: ntpn
integer :: itpn
integer :: ici
integer :: ng
integer :: ipm
integer :: jpm
integer :: notavail
integer :: iws
integer :: j
integer :: ims
integer :: jms
integer :: ko

real :: hpt
real :: vpt
real :: fldval

real :: htpn(maxjms)
real :: vtpn(maxjms)

do iws = 2,mws

   jms = itab_ws(iws)%jm

! Initialize hpt and vpt to zero

   hpt = 0.
   vpt = 0.

! Get tile plot coordinates.  

   do j = 1,jms

      ims = itab_ws(iws)%im(j)
      
      call oplot_transform(iplt           &
                          ,sea%xems(ims)  &
                          ,sea%yems(ims)  &
                          ,sea%zems(ims)  &
                          ,htpn(j)        &
                          ,vtpn(j)        )

      hpt = hpt + htpn(j)
      vpt = vpt + vtpn(j)
   enddo

   hpt = hpt / real(jms)
   vpt = vpt / real(jms)

! Avoid wrap-around for lat-lon plot

   if (op%projectn(iplt) == 'L') then
      do j = 1, jms
         call ll_unwrap(hpt,htpn(j))
      enddo
   endif
   
! Jump out of loop if any cell corner is on other side of earth

   if (any(htpn(1:jms) > 1.e11)) cycle

! Jump out of loop if any is cell corner is outside plot window. 

   if ( all(htpn(1:jms) < op%xmin) .or. all(htpn(1:jms) > op%xmax) .or. &
        all(vtpn(1:jms) < op%ymin) .or. all(vtpn(1:jms) > op%ymax) ) cycle

! Plot cell
   	
   call oplot_lib(iplt,1,iws,'VALUE',op%fldname(iplt),fldval,notavail,ko)    
   if (notavail > 0) cycle 
   call celltile(iplt,iws,itab_ws(iws)%jm,htpn,vtpn,hpt,vpt,fldval,action)

enddo

return
end subroutine tileslab_hs

!===============================================================================

subroutine tileslab_hl(iplt,action)

use oplot_coms, only: op
use mem_leaf,   only: land, itab_wl
use leaf_coms,  only: maxjml, mwl, nzg, nzs
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: ntpn
integer :: ici
integer :: ng
integer :: ipm
integer :: k, ko
integer :: notavail
integer :: iwl
integer :: j
integer :: iml
integer :: jml

real :: hpt
real :: vpt
real :: fldval

real :: htpn(maxjml)
real :: vtpn(maxjml)

! Find K level to plot if field is 3d

if (op%dimens == '3G') then
   k = min(nzg,max(1,nint(op%slabloc(iplt))))
elseif (op%dimens == '3S') then
   k = min(nzs,max(1,nint(op%slabloc(iplt))))
else
   k = 1
endif

do iwl = 2,mwl

   jml = itab_wl(iwl)%jm

! Initialize hpt and vpt to zero

   hpt = 0.
   vpt = 0.

! Get tile plot coordinates.  
   do j = 1,jml

      iml = itab_wl(iwl)%im(j)

      call oplot_transform(iplt            &
                          ,land%xeml(iml)  &
                          ,land%yeml(iml)  &
                          ,land%zeml(iml)  &
                          ,htpn(j)         &
                          ,vtpn(j)         )

      hpt = hpt + htpn(j)
      vpt = vpt + vtpn(j)

   enddo

   hpt = hpt / real(jml)
   vpt = vpt / real(jml)
   
! Avoid wrap-around for lat-lon plots

   if (op%projectn(iplt) == 'L') then
      do j = 1,jml
         call ll_unwrap(hpt,htpn(j))
      enddo
   endif
   
! Jump out of loop if any cell corner is on other side of earth

   if (any(htpn(1:jml) > 1.e11)) cycle

! Jump out of loop if entire cell is outside plot window. 

   if ( all(htpn(1:jml) < op%xmin) .or. all(htpn(1:jml) > op%xmax) .or.  &
        all(vtpn(1:jml) < op%ymin) .or. all(vtpn(1:jml) > op%ymax) ) cycle

! Plot cell
   	
   call oplot_lib(iplt,k,iwl,'VALUE',op%fldname(iplt),fldval,notavail,ko)    
   
   call celltile(iplt,iwl,jml,htpn,vtpn,hpt,vpt,fldval,action)

enddo

return
end subroutine tileslab_hl

!===============================================================================

subroutine tileslab_hfs(iplt,action)

use oplot_coms, only: op
use mem_sea,    only: sea, itab_ws
use leaf_coms,  only: mwl, nzg, nzs
use misc_coms,  only: io6
use mem_sflux,  only: mseaflux,seaflux,xemstrap,yemstrap,zemstrap

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: ntpn
integer :: ici
integer :: ng
integer :: ipm
integer :: k, ko
integer :: notavail
integer :: iws
integer :: j
integer :: ims,jms
integer :: isf,itrap,jtrap,ktrap
integer :: isfprev

real :: hpt
real :: vpt
real :: fldval

real :: htpn(4)
real :: vtpn(4)

k = 1
isfprev = 0

do isf = 2,mseaflux

   jtrap = seaflux(isf)%jtrap
   itrap = seaflux(isf)%itrap

   call oplot_transform(iplt              &
                       ,seaflux(isf)%xef  &
                       ,seaflux(isf)%yef  &
                       ,seaflux(isf)%zef  &
                       ,hpt               &
                       ,vpt               )

   do ktrap = itrap,itrap + jtrap - 1

      do j = 1,4
         call oplot_transform(iplt               &
                             ,xemstrap(j,ktrap)  &
                             ,yemstrap(j,ktrap)  &
                             ,zemstrap(j,ktrap)  &
                             ,htpn(j)            &
                             ,vtpn(j)            )

! Avoid wrap-around for lat-lon plots

         if (op%projectn(iplt) == 'L') call ll_unwrap(hpt,htpn(j))
      enddo

! Jump out of loop if any cell corner is on other side of earth

      if (any(htpn(1:4) > 1.e11)) cycle

! Jump out of loop if entire cell is outside plot window. 

      if ( all(htpn(1:4) < op%xmin) .or. all(htpn(1:4) > op%xmax) .or.  &
           all(vtpn(1:4) < op%ymin) .or. all(vtpn(1:4) > op%ymax) ) cycle

      if (isfprev /= isf) then
         isfprev = isf
         
         call oplot_lib(iplt,k,isf,'VALUE',op%fldname(iplt),fldval,notavail,ko)   
      endif
      
      call celltile(iplt,isf,4,htpn,vtpn,hpt,vpt,fldval,action)
      if (action == 'P') cycle

   enddo

enddo

return
end subroutine tileslab_hfs

!===============================================================================

subroutine tileslab_hfl(iplt,action)

use oplot_coms, only: op
use mem_leaf,   only: land, itab_wl
use leaf_coms,  only: mwl, nzg, nzs
use misc_coms,  only: io6
use mem_sflux,  only: mlandflux,landflux,xemltrap,yemltrap,zemltrap

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: ntpn
integer :: ici
integer :: ng
integer :: ipm
integer :: k, ko
integer :: notavail
integer :: iwl
integer :: j
integer :: iml,jml
integer :: ilf,itrap,jtrap,ktrap
integer :: ilfprev

real :: hpt
real :: vpt
real :: fldval

real :: htpn(4)
real :: vtpn(4)

k = 1
ilfprev = 0

do ilf = 2,mlandflux

   jtrap = landflux(ilf)%jtrap
   itrap = landflux(ilf)%itrap
   
   call oplot_transform(iplt               &
                       ,landflux(ilf)%xef  &
                       ,landflux(ilf)%yef  &
                       ,landflux(ilf)%zef  &
                       ,hpt                &
                       ,vpt                )


   do ktrap = itrap,itrap + jtrap - 1

      do j = 1,4
         call oplot_transform(iplt               &
                             ,xemltrap(j,ktrap)  &
                             ,yemltrap(j,ktrap)  &
                             ,zemltrap(j,ktrap)  &
                             ,htpn(j)            &
                             ,vtpn(j)            )

          if (op%projectn(iplt) == 'L') call ll_unwrap(hpt,htpn(j))
      enddo

! Jump out of loop if any cell corner is on other side of earth

      if (any(htpn(1:4) > 1.e11)) cycle

! Jump out of loop if entire cell is outside plot window. 

      if ( all(htpn(1:4) < op%xmin) .or. all(htpn(1:4) > op%xmax) .or.  &
           all(vtpn(1:4) < op%ymin) .or. all(vtpn(1:4) > op%ymax) ) cycle

      if (ilfprev /= ilf) then
         ilfprev = ilf
         call oplot_lib(iplt,k,ilf,'VALUE',op%fldname(iplt),fldval,notavail,ko)
      endif

      call celltile(iplt,ilf,4,htpn,vtpn,hpt,vpt,fldval,action)
      if (action == 'P') cycle

   enddo

enddo

return
end subroutine tileslab_hfl

!===============================================================================

subroutine tileslab_vt(iplt,action)

use oplot_coms, only: op
use mem_grid,   only: mwa, mza, zm, zt, lpw
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: k, ko
integer :: iw
integer :: ju1
integer :: ju2
integer :: iok
integer :: notavail

real :: hpt
real :: vpt
real :: fldval

real :: htpn(4)
real :: vtpn(4)
real :: topo1,topo2

! Loop over W points

do iw = 2,mwa

! Get horizontal plot coordinates for this W point

   if (op%projectn(iplt) == 'C') then
      call coneplot_w(iw,topo1,topo2,iok,htpn)
   elseif (op%projectn(iplt)(1:1) == 'V') then
      call xyplot_w(iplt,iw,ju1,ju2,topo1,topo2,iok,htpn)
   endif

   hpt = .5 * (htpn(1) + htpn(2))

! Jump out of loop if this W point does not intersect plot cone

   if (iok /= 1) cycle

! Jump out of loop if either cell side is outside plot window. 

   if (htpn(1) < op%xmin .or. htpn(1) > op%xmax .or.  &
       htpn(2) < op%xmin .or. htpn(2) > op%xmax) cycle
   
   do k = 2,mza-1

! Get T-cell vertical coordinates

      vtpn(1) = zm(k-1)
      vtpn(2) = vtpn(1)
      vtpn(3) = zm(k)
      vtpn(4) = vtpn(3)
      vpt = zt(k)

! Check if cell is above ground

      if (k >= lpw(iw)) then 

! Cell is above ground

         call oplot_lib(iplt,k,iw,'VALUE',op%fldname(iplt),fldval,notavail,ko)
         if (notavail > 0) cycle 

         call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

      else
           
! Cell is underground - Tile-plot it with underground color
  
         call fillpolyg(4,htpn(1),vtpn(1),op%icigrnd)
      endif

   enddo

enddo
   
return
end subroutine tileslab_vt

!===============================================================================

subroutine tileslab_vu(iplt,action)

use oplot_coms, only: op
use mem_grid,   only: mua, mza, zm, zt, lpw
use mem_ijtabs, only: itab_u
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: k, ko
integer :: iu
integer :: ng
integer :: im1
integer :: im2
integer :: iw1
integer :: iw2
integer :: ju1
integer :: ju2
integer :: iok
integer :: notavail

real :: hpt
real :: vpt
real :: fldval

real :: htpn(4)
real :: vtpn(4)

real :: htpn1(4)
real :: htpn2(4)


real :: hpt1
real :: hpt2
real :: topo1,topo2

! Loop over U points

do iu = 2,mua

   iw1 = itab_u(iu)%iw1
   iw2 = itab_u(iu)%iw2
   
   if (iw1 < 2 .or. iw2 < 2) cycle

! Get horizontal plot coordinates for IW1 point

   if (op%projectn(iplt) == 'C') then
      call coneplot_w(iw1,topo1,topo2,iok,htpn1)
   elseif (op%projectn(iplt)(1:1) == 'V') then
      call xyplot_w(iplt,iw1,ju1,ju2,topo1,topo2,iok,htpn1)
   endif

! Jump out of loop if this IW1 point does not intersect plot cone

   if (iok /= 1) cycle

   hpt1 = .5 * (htpn1(1) + htpn1(2))

! Get horizontal plot coordinates for IW2 point

   if (op%projectn(iplt) == 'C') then
      call coneplot_w(iw2,topo1,topo2,iok,htpn2)
   elseif (op%projectn(iplt)(1:1) == 'V') then
      call xyplot_w(iplt,iw2,ju1,ju2,topo1,topo2,iok,htpn2)
   endif

! Jump out of loop if this IW2 point does not intersect plot cone

   if (iok /= 1) cycle

   hpt2 = .5 * (htpn2(1) + htpn2(2))

! Avoid wrap_around - use htpn that is between hpt1 and hpt2

   if (op%projectn(iplt) == 'C') then

      if ((hpt1 < htpn(1) .and. htpn(1) < hpt2) .or.  &
          (hpt1 > htpn(1) .and. htpn(1) > hpt2)) then
         call ll_unwrap(htpn(1),hpt1)
         call ll_unwrap(htpn(1),hpt2)
      else
         call ll_unwrap(htpn(2),hpt1)
         call ll_unwrap(htpn(2),hpt2)
      endif
   endif

! Jump out of loop if either cell side is outside plot window. 

   if (hpt1 < op%xmin .or. hpt1 > op%xmax .or.  &
       hpt2 < op%xmin .or. hpt2 > op%xmax) cycle
   
   do k = 2,mza-1  ! Loop is over T levels

! Get U-cell vertical coordinates

      vtpn(1) = zm(k-1)
      vtpn(2) = vtpn(1)
      vtpn(3) = zm(k)
      vtpn(4) = vtpn(3)
      vpt = zt(k)

! Check if both IW1 and IW2 are above ground

      if (iw1 > 1 .and. k >= lpw(iw1) .and.  &
          iw2 > 1 .and. k >= lpw(iw2)) then 

! Yes - plot full cell with tile color

         htpn(1) = hpt1
         htpn(2) = hpt2
         htpn(3) = htpn(2)
         htpn(4) = htpn(1)
         hpt = .5 * (hpt1 + hpt2)
   
         call oplot_lib(iplt,k,iu,'VALUE',op%fldname(iplt),fldval,notavail,ko)
         if (notavail > 0) cycle 
         
         call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)

! Else, check if only IW1 is above ground

      elseif (iw1 > 1 .and. k >= lpw(iw1)) then 

! Yes - Make sure IW2 cell is not at eastern boundary

         if (iw2 /= 1) then 

! OK - plot IW2 cell with underground color

            htpn(1) = .5 * (hpt1 + hpt2)
            htpn(2) = hpt2
            htpn(3) = htpn(2)
            htpn(4) = htpn(1)

            call fillpolyg(4,htpn,vtpn,op%icigrnd)

         endif

! Plot IW1 half cell with tile color

         htpn(1) = .5 * (hpt1 + hpt2)
         htpn(2) = hpt1
         htpn(3) = htpn(2)
         htpn(4) = htpn(1)
         hpt = .5 * (hpt1 + hpt2)

         call oplot_lib(iplt,k,iu,'VALUE',op%fldname(iplt),fldval,notavail,ko)
         if (notavail > 0) cycle 
         call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)

! Else, check if only IW2 is above ground

      elseif (iw2 > 1 .and. k >= lpw(iw2)) then 

! Yes - make sure IW1 cell is not at western boundary

         if (iw1 /= 1) then 

! OK - plot IW1 cell with underground color

            htpn(1) = .5 * (hpt1 + hpt2)
            htpn(2) = hpt1
            htpn(3) = htpn(2)
            htpn(4) = htpn(1)

            call fillpolyg(4,htpn,vtpn,op%icigrnd)

         endif

! Plot IW2 half cell with tile color

         htpn(1) = .5 * (hpt1 + hpt2)
         htpn(2) = hpt2
         htpn(3) = htpn(2)
         htpn(4) = htpn(1)
         hpt = .5 * (hpt1 + hpt2)

         call oplot_lib(iplt,k,iu,'VALUE',op%fldname(iplt),fldval,notavail,ko)
         if (notavail > 0) cycle 
         call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)

      else

! At this point, both IW1 and IW2 are below ground or at E/W boundary -
! Plot one or both with underground color as appropriate     

! Check if IW1 cell is at western boundary

         htpn(1) = hpt1
         htpn(2) = hpt2
         htpn(3) = htpn(2)
         htpn(4) = htpn(1)

         call fillpolyg(4,htpn,vtpn,op%icigrnd)

      endif

   enddo

enddo

return
end subroutine tileslab_vu

!===============================================================================

subroutine tileslab_vw(iplt,action)

use oplot_coms, only: op
use mem_grid,   only: mwa, mza, zm, zt, lpw
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: k, ko
integer :: iw
integer :: ng
integer :: iv1
integer :: iv2
integer :: iu1
integer :: iu2
integer :: ju1
integer :: ju2
integer :: iok
integer :: notavail

real :: hpt
real :: vpt
real :: fldval

real :: htpn(4)
real :: vtpn(4)
real :: topo1,topo2

! Loop over W points

do iw = 2,mwa

! Get horizontal plot coordinates for this W point

   if (op%projectn(iplt) == 'C') then
      call coneplot_w(iw,topo1,topo2,iok,htpn)
   elseif (op%projectn(iplt) == 'V') then
      call xyplot_w(iplt,iw,ju1,ju2,topo1,topo2,iok,htpn)
   endif

   hpt = .5 * (htpn(1) + htpn(2))

! Jump out of loop if this W point does not intersect plot cone

   if (iok /= 1) cycle

! Jump out of loop if either cell side is outside plot window. 

   if (htpn(1) < op%xmin .or. htpn(1) > op%xmax .or.  &
       htpn(2) < op%xmin .or. htpn(2) > op%xmax) cycle
   
   do k = 1,mza-1
      vpt = zm(k)
      
! Check if lower T cell is above ground

      if (k >= lpw(iw)) then 

! Yes - plot full cell with tile color

         vtpn(1) = zt(k)
         vtpn(2) = vtpn(1)
         vtpn(3) = min(zt(k+1),zm(mza-1))
         vtpn(4) = vtpn(3)
         
         call oplot_lib(iplt,k,iw,'VALUE',op%fldname(iplt),fldval,notavail,ko)
         if (notavail > 0) cycle 
         call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)

! Else, check if only upper T cell is above ground

      elseif (k+1 >= lpw(iw)) then 

! Yes - plot upper half W cell with tile color

         vtpn(1) = zm(k)
         vtpn(2) = vtpn(1)
         vtpn(3) = zt(k+1)
         vtpn(4) = vtpn(3)

         call oplot_lib(iplt,k,iw,'VALUE',op%fldname(iplt),fldval,notavail,ko)
         if (notavail > 0) cycle 
         call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)
         
      else
      
! At this point, upper T cell is below ground - Plot with underground color

         vtpn(1) = zm(k)
         vtpn(2) = vtpn(1)
         vtpn(3) = zm(k+1)
         vtpn(4) = vtpn(3)

         call fillpolyg(4,htpn(1),vtpn(1),op%icigrnd)

      endif

   enddo

! Plot top half W cell with tile color

!   vtpn(1) = zt(mza-1)
!   vtpn(2) = vtpn(1)
!   vtpn(3) = zm(mza-1)
!   vtpn(4) = vtpn(3)
!   vpt = zm(mza-1)

!   call oplot_lib(iplt,mza-1,iw,'VALUE',op%fldname(iplt),fldval,notavail,ko)
!   call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval)
      
enddo
   
return
end subroutine tileslab_vw

!===============================================================================

subroutine tileslab_vl(iplt,action)

use oplot_coms, only: op
use leaf_coms,  only: mwl, nzg, nzs
use mem_leaf,   only: land
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: k, ko
integer :: iw
integer :: ng
integer :: iv1
integer :: iv2
integer :: iu1
integer :: iu2
integer :: ip
integer :: itpn
integer :: notavail

real :: hpt
real :: vpt
real :: fldval
real :: patwidth
real :: botk
real :: delzk

real :: htpn(4)
real :: vtpn(4)

!!!!!!!! VERTICAL XSECTION NEEDS WORK

do ip = 2,mwl
!   iw = land%iw(ip)

! Skip iw column if not intersected by slabloc or if outside window bounds


! Get horizontal plot coordinates for cells in this column

   op%stagpt = 'LA'  ! Get land cell fractional area
   call oplot_lib(iplt,0,ip,'VALUE',op%fldname(iplt),fldval,notavail,ko)

!   patwidth = (xu(iu2) - xu(iu1)) * fldval

!   htpn(1) = vt2da(iw)
!   htpn(2) = htpn(1) + patwidth
!   htpn(3) = htpn(2)
!   htpn(4) = htpn(1)

!   vt2da(iw) = htpn(2)

!   hpt = .5 * (htpn(1) + htpn(2))

   botk = - real(nzg+nzs+4)   ! level of bottom of bottom soil layer (negative)
   delzk = op%ymin / botk  ! This is a positive number

   do k = 1,nzg  ! Loop over soil layers

! Get vertical coordinates for soil layers

      vtpn(1) = delzk * (float(k-1) + botk) 
      vtpn(2) = vtpn(1)
      vtpn(3) = delzk * (float(k)   + botk)
      vtpn(4) = vtpn(3)
      vpt = .5 * (vtpn(1) + vtpn(3))

! plot soil layers

      op%stagpt = 'L'  ! Get soil value
      call oplot_lib(iplt,k,ip,'VALUE',op%fldname(iplt),fldval,notavail,ko)
      call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

      if (op%pltgrid(iplt) == 'G') then
         call o_frstpt(htpn(4),vtpn(4))
         do itpn = 1,4
            call o_vector(htpn(itpn),vtpn(itpn))
         enddo
      endif

   enddo

   do k = 1,nzs  ! Loop over sfcwater layers

! check for existence of sfcwater in current level k

      call oplot_lib(iplt,k,ip,'VALUE','SFCWATER_MASS',fldval,notavail,ko)

      if (fldval > 0.) then
      
! Get vertical coordinates for sfcwater layer k

         vtpn(1) = delzk * (float(nzg+k-1) + .5 + botk) 
         vtpn(2) = vtpn(1)
         vtpn(3) = delzk * (float(nzg+k)   + .5 + botk) 
         vtpn(4) = vtpn(3)
         vpt = .5 * (vtpn(1) + vtpn(3))

! plot sfcwater layer k

         op%stagpt = 'LW'  ! Get sfcwater value
         call oplot_lib(iplt,k,ip,'VALUE',op%fldname(iplt),fldval,notavail,ko)
         call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

         if (op%pltgrid(iplt) == 'G') then
            call o_frstpt(htpn(4),vtpn(4))
            do itpn = 1,4
               call o_vector(htpn(itpn),vtpn(itpn))
            enddo
         endif
         
      endif
      
   enddo

! plot vegetation layer

   vtpn(1) = delzk * (float(nzg+nzs) + 1. + botk) 
   vtpn(2) = vtpn(1)
   vtpn(3) = delzk * (float(nzg+nzs+1) + 1. + botk) 
   vtpn(4) = vtpn(3)
   vpt = .5 * (vtpn(1) + vtpn(3))

   op%stagpt = 'LV'  ! Get vegetation value
   call oplot_lib(iplt,0,ip,'VALUE',op%fldname(iplt),fldval,notavail,ko)
   call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

   if (op%pltgrid(iplt) == 'G') then
      call o_frstpt(htpn(4),vtpn(4))
      do itpn = 1,4
         call o_vector(htpn(itpn),vtpn(itpn))
      enddo
   endif

! plot canopy air layer

   vtpn(1) = delzk * (float(nzg+nzs+1) + 1.5 + botk) 
   vtpn(2) = vtpn(1)
   vtpn(3) = delzk * (float(nzg+nzs+2) + 1.5 + botk) 
   vtpn(4) = vtpn(3)
   vpt = .5 * (vtpn(1) + vtpn(3))

   op%stagpt = 'LC'  ! Get vegetation value
   call oplot_lib(iplt,0,ip,'VALUE',op%fldname(iplt),fldval,notavail,ko)
   call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

   if (op%pltgrid(iplt) == 'G') then
      call o_frstpt(htpn(4),vtpn(4))
      do itpn = 1,4
         call o_vector(htpn(itpn),vtpn(itpn))
      enddo
   endif

enddo
   
return
end subroutine tileslab_vl

!===============================================================================

subroutine tileslab_vs(iplt,action)

use oplot_coms, only: op
use leaf_coms,  only: nzg, nzs
use sea_coms,   only: mws
use mem_sea,    only: sea
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: k, ko
integer :: iw
integer :: ng
integer :: iv1
integer :: iv2
integer :: iu1
integer :: iu2
integer :: ip
integer :: itpn
integer :: notavail

real :: hpt
real :: vpt
real :: fldval
real :: patwidth
real :: botk
real :: delzk

real :: htpn(4)
real :: vtpn(4)

!!!!!!!! VERTICAL XSECTION NEEDS WORK

do ip = 2,mws
!   iw = sea%iw(ip)

! Skip iw column if not intersected by slabloc or if outside window bounds

! Get horizontal plot coordinates for cells in this column

   op%stagpt = 'SA'  ! Get sea cell fractional area
   call oplot_lib(iplt,0,ip,'VALUE',op%fldname(iplt),fldval,notavail,ko)

   botk = - real(nzg+nzs+4)   ! level of bottom of bottom soil layer (negative)
   delzk = op%ymin / botk  ! This is a positive number

! plot (top) sea layer

   vtpn(1) = delzk * (float(nzg-1) + botk) 
   vtpn(2) = vtpn(1)
   vtpn(3) = delzk * (float(nzg) + botk) 
   vtpn(4) = vtpn(3)
   vpt = .5 * (vtpn(1) + vtpn(3))

   op%stagpt = 'S'  ! Get vegetation value
   call oplot_lib(iplt,0,ip,'VALUE',op%fldname(iplt),fldval,notavail,ko)
   call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

   if (op%pltgrid(iplt) == 'G') then
      call o_frstpt(htpn(4),vtpn(4))
      do itpn = 1,4
         call o_vector(htpn(itpn),vtpn(itpn))
      enddo
   endif

! plot canopy air layer

   vtpn(1) = delzk * (float(nzg+nzs+1) + 1.5 + botk) 
   vtpn(2) = vtpn(1)
   vtpn(3) = delzk * (float(nzg+nzs+2) + 1.5 + botk) 
   vtpn(4) = vtpn(3)
   vpt = .5 * (vtpn(1) + vtpn(3))

   op%stagpt = 'SC'  ! Get canopy air value
   call oplot_lib(iplt,0,ip,'VALUE',op%fldname(iplt),fldval,notavail,ko)
   call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

   if (op%pltgrid(iplt) == 'G') then
      call o_frstpt(htpn(4),vtpn(4))
      do itpn = 1,4
         call o_vector(htpn(itpn),vtpn(itpn))
      enddo
   endif

enddo
   
return
end subroutine tileslab_vs

!===============================================================================

subroutine celltile(iplt,i,ntpn,htpn,vtpn,hpt,vpt,fldval,action)

use oplot_coms, only: op
use plotcolors, only: clrtab
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt
integer, intent(in) :: i
integer, intent(in) :: ntpn

real, intent(in) :: htpn(*)
real, intent(in) :: vtpn(*)
real, intent(in) :: hpt
real, intent(in) :: vpt
real, intent(in) :: fldval

character(1), intent(in) :: action

real    :: fldval1
integer :: icolor
integer :: itab
integer :: ival

if (action == 'T') then

   itab = op%icolortab(iplt)

! Cyclic treatment of color palette (used with integer-type data)

   fldval1 = fldval
   if (clrtab(itab)%ifmt(1) == 20)  &
      fldval1 = mod(fldval-1.,real(clrtab(itab)%nvals-2)) + 1.

! Extract contour color from color table

   ival = 1
   do while (fldval1 > clrtab(itab)%vals(ival) .and.  &
               ival < clrtab(itab)%nvals             )
      ival = ival + 1
   enddo
   icolor = clrtab(itab)%ipal(ival)
   
   call fillpolyg(ntpn,htpn,vtpn,icolor)
endif

if (action == 'P') then
   if ( (hpt > op%xmin) .and. (hpt < op%xmax) .and. &
        (vpt > op%ymin) .and. (vpt < op%ymax) ) then
      call oplot_prtvalue(fldval,hpt,vpt,op%vsprd,.7*op%psiz,op%icolortab(iplt))
   endif
endif

return
end subroutine celltile

!===============================================================================

subroutine coneplot_w(iw,topo1,topo2,iok,htpn)

use oplot_coms,  only: op, xepc, yepc, zepc
use mem_grid,    only: xem, yem, zem, topm
use mem_ijtabs,  only: itab_w
use consts_coms, only: pio180, erad, piu180, pi2
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iw
integer, intent(out) :: iok
real, intent(out) :: htpn(4)
real, intent(out) :: topo1,topo2

integer :: im1,im2,im3

real :: unxec3,unyec3,unzec3
real :: xcent,ycent,zcent
real :: xwincent,ywincent,zwincent
real :: radcone
real :: raxis
real :: vecx,vecy,vecz
real :: vecleftx,veclefty,vecleftz
real :: valx,valy
real :: angm1,angm2,angm3
real :: angmin,angmax
real :: ang0
real :: wt1,wt2
real :: xrad1,yrad1,zrad1
real :: xrad2,yrad2,zrad2
real :: ange1,ange2
real :: sinplat,cosplat
real :: sinvaz,cosvaz
real :: sinconang,cosconang
real :: sinconlat,cosconlat
real :: sinconlon,cosconlon

real :: xemc(3),yemc(3),zemc(3)  ! Coords of 3 M points
real :: topmc(3) ! Topo height of 3 M points

sinplat = sin(op%plat3 * pio180)
cosplat = cos(op%plat3 * pio180)

! cone is viewed from INSIDE, so cone axis is in opposite direction from viewazim 

sinvaz = -sin(op%viewazim * pio180)
cosvaz = -cos(op%viewazim * pio180)

sinconang = sin(op%coneang * pio180)
cosconang = cos(op%coneang * pio180)

! Use formulas (5-5) and (5-6) of USGS Map Projections Manual to get lat/lon
! of plot-cone axis given lat/lon of plot center, plot-cone angle, and view azimuth
! (Avoid special case of where cone center is +/- 90 deg longitude from plot center)

op%conelat = asin(sinplat * cosconang + cosplat * sinconang * cosvaz)

if (abs(cosplat * cosconang - sinplat * sinconang * cosvaz) > 1.e-6) then

   op%conelon = op%plon3 * pio180                                   &
              + atan2(sinconang * sinvaz,                           &
               (cosplat * cosconang - sinplat * sinconang * cosvaz))

elseif (op%viewazim < 180.) then
   op%conelon = (op%plon3 - 90.) * pio180
else
   op%conelon = (op%plon3 + 90.) * pio180
endif

sinconlat = sin(op%conelat)
cosconlat = cos(op%conelat)

sinconlon = sin(op%conelon)
cosconlon = cos(op%conelon)

! Earth components of unit vector outward along cone axis

unxec3 = cosconlat * cosconlon
unyec3 = cosconlat * sinconlon
unzec3 = sinconlat

! Intersection of cone and earth is a circle - find earth coords of circle center

xcent = unxec3 * erad * cosconang 
ycent = unyec3 * erad * cosconang 
zcent = unzec3 * erad * cosconang 

! Cone radius

radcone = erad * sinconang

! Earth coordinates of plot window center

zwincent = erad * sinplat
raxis = erad * cosplat
xwincent = raxis * cos(op%plon3 * pio180) 
ywincent = raxis * sin(op%plon3 * pio180) 

! Indexes of 3 M points for this IW triangle

im1 = itab_w(iw)%im1
im2 = itab_w(iw)%im2
im3 = itab_w(iw)%im3

! Compute angle of each of 3 triangle corners with plot cone center using
! dot products (Adequate precision may require cross product for cone angles 
! close to 0 or 180)

angm1 = acos((xem(im1)*unxec3 + yem(im1)*unyec3 + zem(im1)*unzec3) / erad) * piu180
angm2 = acos((xem(im2)*unxec3 + yem(im2)*unyec3 + zem(im2)*unzec3) / erad) * piu180
angm3 = acos((xem(im3)*unxec3 + yem(im3)*unyec3 + zem(im3)*unzec3) / erad) * piu180

angmin = min(angm1,angm2,angm3)
angmax = max(angm1,angm2,angm3)

! Return with iok = 0 if iw column is not intersected by plot slab

iok = 0

if (angmin > op%coneang .or. angmax < op%coneang) return

iok = 1  ! Since we got here, IW column is intersected by plot cone

! Find IM point that is closest to cone axis and copy to temporary M point #1
! Fill other 2 points in cyclic order

if (angm1 <= angm2 .and. angm1 <= angm3) then  ! m1 is inner point
   xemc(1) = xem(im1)
   yemc(1) = yem(im1)
   zemc(1) = zem(im1)
   topmc(1) = topm(im1)

   xemc(2) = xem(im2)
   yemc(2) = yem(im2)
   zemc(2) = zem(im2)
   topmc(2) = topm(im2)

   xemc(3) = xem(im3)
   yemc(3) = yem(im3)
   zemc(3) = zem(im3)
   topmc(3) = topm(im3)
elseif (angm2 <= angm1 .and. angm2 <= angm3) then  ! m2 is inner point
   xemc(1) = xem(im2)
   yemc(1) = yem(im2)
   zemc(1) = zem(im2)
   topmc(1) = topm(im2)

   xemc(2) = xem(im3)
   yemc(2) = yem(im3)
   zemc(2) = zem(im3)
   topmc(2) = topm(im3)

   xemc(3) = xem(im1)
   yemc(3) = yem(im1)
   zemc(3) = zem(im1)
   topmc(3) = topm(im1)
      
   ang0 = angm1
   angm1 = angm2
   angm2 = angm3
   angm3 = ang0
elseif (angm3 <= angm1 .and. angm3 <= angm2) then  ! m3 is inner point
   xemc(1) = xem(im3)
   yemc(1) = yem(im3)
   zemc(1) = zem(im3)
   topmc(1) = topm(im3)

   xemc(2) = xem(im1)
   yemc(2) = yem(im1)
   zemc(2) = zem(im1)
   topmc(2) = topm(im1)

   xemc(3) = xem(im2)
   yemc(3) = yem(im2)
   zemc(3) = zem(im2)
   topmc(3) = topm(im2)

   ang0 = angm1
   angm1 = angm3
   angm3 = angm2
   angm2 = ang0
endif
   
! Find two points of intersection between current IW triangle and cone

if (angm2 > op%coneang .or. angm2 > angm3) then
   wt2 = (op%coneang - angm1) / (angm2 - angm1)
   wt1 = 1. - wt2
   xepc(2) = wt1 * xemc(1) + wt2 * xemc(2)
   yepc(2) = wt1 * yemc(1) + wt2 * yemc(2)
   zepc(2) = wt1 * zemc(1) + wt2 * zemc(2)
   topo2   = wt1 * topmc(1) + wt2 * topmc(2) ! topo height(2)

   if (angm3 > op%coneang) then
      wt2 = (op%coneang - angm1) / (angm3 - angm1)
      wt1 = 1. - wt2
      xepc(1) = wt1 * xemc(1) + wt2 * xemc(3)
      yepc(1) = wt1 * yemc(1) + wt2 * yemc(3)
      zepc(1) = wt1 * zemc(1) + wt2 * zemc(3)
      topo1   = wt1 * topmc(1) + wt2 * topmc(3) ! topo height(1)
   else
      wt2 = (op%coneang - angm3) / (angm2 - angm3)
      wt1 = 1. - wt2
      xepc(1) = wt1 * xemc(3) + wt2 * xemc(2)
      yepc(1) = wt1 * yemc(3) + wt2 * yemc(2)
      zepc(1) = wt1 * zemc(3) + wt2 * zemc(2)
      topo1   = wt1 * topmc(3) + wt2 * topmc(2) ! topo height(1)
   endif
   
elseif (angm2 < angm3) then
   
   wt2 = (op%coneang - angm2) / (angm3 - angm2)
   wt1 = 1. - wt2
   xepc(2) = wt1 * xemc(2) + wt2 * xemc(3)
   yepc(2) = wt1 * yemc(2) + wt2 * yemc(3)
   zepc(2) = wt1 * zemc(2) + wt2 * zemc(3)
   topo2   = wt1 * topmc(2) + wt2 * topmc(3) ! topo height(2)

   wt2 = (op%coneang - angm1) / (angm3 - angm1)
   wt1 = 1. - wt2
   xepc(1) = wt1 * xemc(1) + wt2 * xemc(3)
   yepc(1) = wt1 * yemc(1) + wt2 * yemc(3)
   zepc(1) = wt1 * zemc(1) + wt2 * zemc(3)
   topo1   = wt1 * topmc(1) + wt2 * topmc(3) ! topo height(1)

else ! angm2 = angm3

   wt2 = (op%coneang - angm1) / (angm3 - angm1)
   wt1 = 1. - wt2
   xepc(1) = wt1 * xemc(1) + wt2 * xemc(3)
   yepc(1) = wt1 * yemc(1) + wt2 * yemc(3)
   zepc(1) = wt1 * zemc(1) + wt2 * zemc(3)
   topo1   = wt1 * topmc(1) + wt2 * topmc(3) ! topo height(1)

   wt2 = (op%coneang - angm1) / (angm2 - angm1)
   wt1 = 1. - wt2
   xepc(2) = wt1 * xemc(1) + wt2 * xemc(2)
   yepc(2) = wt1 * yemc(1) + wt2 * yemc(2)
   zepc(2) = wt1 * zemc(1) + wt2 * zemc(2)
   topo2   = wt1 * topmc(1) + wt2 * topmc(2) ! topo height(1)

endif

! Transform horizontal point coordinates

! Components of vector from circle center to plot window center

vecx = xwincent - xcent
vecy = ywincent - ycent
vecz = zwincent - zcent

! Components of vector 90 degrees to the left (in azimuth) from preceding vector

vecleftx = unyec3 * vecz - unzec3 * vecy
veclefty = unzec3 * vecx - unxec3 * vecz
vecleftz = unxec3 * vecy - unyec3 * vecx

! Compute dot product between vector from circle center to plot center
! and vector from circle center to current point

valx = vecx * (xepc(1) - xcent)  &
     + vecy * (yepc(1) - ycent)  &
     + vecz * (zepc(1) - zcent)

! Repeat with 90-left vector

valy = vecleftx * (xepc(1) - xcent)  &
     + veclefty * (yepc(1) - ycent)  &
     + vecleftz * (zepc(1) - zcent)

ange1 = atan2(-valy,valx)  ! Angle increases clockwise

! Repeat dot product for second point

valx = vecx * (xepc(2) - xcent)  &
     + vecy * (yepc(2) - ycent)  &
     + vecz * (zepc(2) - zcent)

valy = vecleftx * (xepc(2) - xcent)  &
     + veclefty * (yepc(2) - ycent)  &
     + vecleftz * (zepc(2) - zcent)

ange2 = atan2(-valy,valx)  ! Angle increases clockwise

! Avoid wrap_around

if (ange2 < ange1) ange2 = ange2 + pi2

! Scale angles to htpn coordinates (in meters along cone circle)

htpn(1) = ange1 * radcone 
htpn(2) = ange2 * radcone
htpn(3) = htpn(2)
htpn(4) = htpn(1)

return
end subroutine coneplot_w

!===============================================================================

subroutine xplot_w(iplt,iw,iok,htpn)

! This subroutine finds points of intersection between vertical plot slab and
! triangular (prism) grid cells

use oplot_coms, only: op
use mem_grid, only: xem, yem
use mem_ijtabs, only: itab_w
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt
integer, intent(in) :: iw

integer, intent(out) :: iok
real, intent(out) :: htpn(4)

integer :: im1
integer :: im2
integer :: im3

real :: s1
real :: s2
real :: wt1
real :: wt2

real :: xem3(3)  ! Coords of 3 M points
real :: yem3(3)  ! Coords of 3 M points

im1 = itab_w(iw)%im1
im2 = itab_w(iw)%im2
im3 = itab_w(iw)%im3

xem3(1) = xem(im1)
xem3(2) = xem(im2)
xem3(3) = xem(im3)

yem3(1) = yem(im1)
yem3(2) = yem(im2)
yem3(3) = yem(im3)

s1 = min(xem3(1),xem3(2),xem3(3))
s2 = max(xem3(1),xem3(2),xem3(3))

! Return with iok = 0 if iw column is not intersected by plot cone

iok = 0

if (s1 > op%slabloc(iplt) .or. s2 < op%slabloc(iplt)) return

iok = 1  ! Since we got here, IW column is intersected by plot cone

! Find IM point that has lowest X coordinate and copy to temporary M point #1
! Fill other 2 points in cyclic order

if (xem(im1) <= xem(im2) .and. xem(im1) <= xem(im3)) then  ! m1 is inner point
   xem3(1) = xem(im1)
   yem3(1) = yem(im1)

   xem3(2) = xem(im2)
   yem3(2) = yem(im2)

   xem3(3) = xem(im3)
   yem3(3) = yem(im3)
elseif (xem(im2) <= xem(im1) .and. xem(im2) <= xem(im3)) then  ! m2 is inner point
   xem3(1) = xem(im2)
   yem3(1) = yem(im2)

   xem3(2) = xem(im3)
   yem3(2) = yem(im3)

   xem3(3) = xem(im1)
   yem3(3) = yem(im1)
elseif (xem(im3) <= xem(im1) .and. xem(im3) <= xem(im2)) then  ! m3 is inner point
   xem3(1) = xem(im3)
   yem3(1) = yem(im3)

   xem3(2) = xem(im1)
   yem3(2) = yem(im1)

   xem3(3) = xem(im2)
   yem3(3) = yem(im2)
endif
   
! Find two points of intersection between current IW triangle and X slab

if (xem3(2) > op%slabloc(iplt)) then
   wt2 = (op%slabloc(iplt) - xem3(1)) / (xem3(2) - xem3(1))
   wt1 = 1. - wt2
   htpn(1) = wt1 * yem3(1) + wt2 * yem3(2)

   if (xem3(3) > op%slabloc(iplt)) then
      wt2 = (op%slabloc(iplt) - xem3(3)) / (xem3(1) - xem3(3))
      wt1 = 1. - wt2
      htpn(2) = wt1 * yem3(3) + wt2 * yem3(1)
   else
      wt2 = (op%slabloc(iplt) - xem3(3)) / (xem3(2) - xem3(3))
      wt1 = 1. - wt2
      htpn(2) = wt1 * yem3(3) + wt2 * yem3(2)
   endif
else
   wt2 = (op%slabloc(iplt) - xem3(2)) / (xem3(3) - xem3(2))
   wt1 = 1. - wt2
   htpn(1) = wt1 * yem3(2) + wt2 * yem3(3)

   wt2 = (op%slabloc(iplt) - xem3(3)) / (xem3(1) - xem3(3))
   wt1 = 1. - wt2
   htpn(2) = wt1 * yem3(3) + wt2 * yem3(1)
endif

htpn(3) = htpn(2)
htpn(4) = htpn(1)

return
end subroutine xplot_w

!===============================================================================

subroutine yplot_w(iplt,iw,iok,htpn)

! This subroutine finds points of intersection between vertical plot slab and
! triangular (prism) grid cells

use oplot_coms, only: op
use mem_grid,   only: xem, yem
use mem_ijtabs, only: itab_w
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt
integer, intent(in) :: iw
integer, intent(out) :: iok
real, intent(out) :: htpn(4)

integer :: im1
integer :: im2
integer :: im3

real :: s1
real :: s2
real :: wt1
real :: wt2

real :: xem3(3)  ! Coords of 3 M points
real :: yem3(3)  ! Coords of 3 M points

im1 = itab_w(iw)%im1
im2 = itab_w(iw)%im2
im3 = itab_w(iw)%im3

xem3(1) = xem(im1)
xem3(2) = xem(im2)
xem3(3) = xem(im3)

yem3(1) = yem(im1)
yem3(2) = yem(im2)
yem3(3) = yem(im3)

s1 = min(yem3(1),yem3(2),yem3(3))
s2 = max(yem3(1),yem3(2),yem3(3))

! Return with iok = 0 if iw column is not intersected by plot cone

iok = 0

if (s1 > op%slabloc(iplt) .or. s2 < op%slabloc(iplt)) return

iok = 1  ! Since we got here, IW column is intersected by plot cone

! Find IM point that has lowest X coordinate and copy to temporary M point #1
! Fill other 2 points in cyclic order

if (yem(im1) <= yem(im2) .and. yem(im1) <= yem(im3)) then  ! m1 is inner point
   xem3(1) = xem(im1)
   yem3(1) = yem(im1)

   xem3(2) = xem(im2)
   yem3(2) = yem(im2)

   xem3(3) = xem(im3)
   yem3(3) = yem(im3)
elseif (yem(im2) <= yem(im1) .and. yem(im2) <= yem(im3)) then  ! m2 is inner point
   xem3(1) = xem(im2)
   yem3(1) = yem(im2)

   xem3(2) = xem(im3)
   yem3(2) = yem(im3)

   xem3(3) = xem(im1)
   yem3(3) = yem(im1)
elseif (yem(im3) <= yem(im1) .and. yem(im3) <= yem(im2)) then  ! m3 is inner point
   xem3(1) = xem(im3)
   yem3(1) = yem(im3)

   xem3(2) = xem(im1)
   yem3(2) = yem(im1)

   xem3(3) = xem(im2)
   yem3(3) = yem(im2)
endif

! Find two points of intersection between current IW triangle and X slab

if (yem3(2) > op%slabloc(iplt)) then
   wt2 = (op%slabloc(iplt) - yem3(1)) / (yem3(2) - yem3(1))
   wt1 = 1. - wt2
   htpn(1) = wt1 * xem3(1) + wt2 * xem3(2)

   if (yem3(3) > op%slabloc(iplt)) then
      wt2 = (op%slabloc(iplt) - yem3(3)) / (yem3(1) - yem3(3))
      wt1 = 1. - wt2
      htpn(2) = wt1 * xem3(3) + wt2 * xem3(1)
   else
      wt2 = (op%slabloc(iplt) - yem3(3)) / (yem3(2) - yem3(3))
      wt1 = 1. - wt2
      htpn(2) = wt1 * xem3(3) + wt2 * xem3(2)
   endif
else
   wt2 = (op%slabloc(iplt) - yem3(2)) / (yem3(3) - yem3(2))
   wt1 = 1. - wt2
   htpn(1) = wt1 * xem3(2) + wt2 * xem3(3)

   wt2 = (op%slabloc(iplt) - yem3(3)) / (yem3(1) - yem3(3))
   wt1 = 1. - wt2
   htpn(2) = wt1 * xem3(3) + wt2 * xem3(1)
endif

htpn(3) = htpn(2)
htpn(4) = htpn(1)

return
end subroutine yplot_w

!===============================================================================

subroutine xyplot_w(iplt,iw,ju1,ju2,topo1,topo2,iok,htpn)

! This subroutine finds points of intersection between vertical plot slab and
! triangular (prism) grid cells

use oplot_coms,  only: op
use mem_grid,    only: xem, yem, topm
use mem_ijtabs,  only: itab_w
use consts_coms, only: pio180
use oname_coms,  only: nl
use misc_coms,  only: io6

implicit none

integer, intent(in)  :: iplt
integer, intent(in)  :: iw
integer, intent(out) :: ju1
integer, intent(out) :: ju2
integer, intent(out) :: iok

real, intent(out) :: htpn(4)
real, intent(out) :: topo1,topo2

integer :: im1
integer :: im2
integer :: im3
integer :: iu1
integer :: iu2
integer :: iu3
integer :: iuc1
integer :: iuc2
integer :: iuc3

real :: s1
real :: s2
real :: s3
real :: sc1
real :: sc2
real :: sc3
real :: smin
real :: smax
real :: wt1
real :: wt2
real :: x1
real :: x2
real :: y1
real :: y2
real :: sinvaz
real :: cosvaz

real :: xemc(3)  ! Coords of 3 M points
real :: yemc(3)  ! Coords of 3 M points
real :: topmc(3) ! Topo height of 3 M points

im1 = itab_w(iw)%im1
im2 = itab_w(iw)%im2
im3 = itab_w(iw)%im3

iu1 = itab_w(iw)%iu1
iu2 = itab_w(iw)%iu2
iu3 = itab_w(iw)%iu3

sinvaz = sin((90. - op%viewazim) * pio180)
cosvaz = cos((90. - op%viewazim) * pio180)

! Location of 3 M points along line perpendicular to plot slab

s1 = (xem(im1) - nl%plotcoord1(iplt)) * cosvaz  &
   + (yem(im1) - nl%plotcoord2(iplt)) * sinvaz

s2 = (xem(im2) - nl%plotcoord1(iplt)) * cosvaz  &
   + (yem(im2) - nl%plotcoord2(iplt)) * sinvaz

s3 = (xem(im3) - nl%plotcoord1(iplt)) * cosvaz  &
   + (yem(im3) - nl%plotcoord2(iplt)) * sinvaz

smin = min(s1,s2,s3)
smax = max(s1,s2,s3)

! Return with iok = 0 if iw column is not intersected by plot slab

iok = 0

if (smin > op%slabloc(iplt) .or. smax < op%slabloc(iplt)) return

iok = 1  ! Since we got here, IW column is intersected by plot cone

! Find IM point that has lowest S coordinate and copy to temporary M point #1
! Fill other 2 points in cyclic order

if (s1 <= s2 .and. s1 <= s3) then  ! m1 has lowest S value
   xemc(1)  = xem(im1)
   yemc(1)  = yem(im1)
   topmc(1) = topm(im1)
   
   xemc(2)  = xem(im2)
   yemc(2)  = yem(im2)
   topmc(2) = topm(im2)

   xemc(3)  = xem(im3)
   yemc(3)  = yem(im3)
   topmc(3) = topm(im3)

   sc1     = s1
   sc2     = s2
   sc3     = s3

   iuc1    = iu1
   iuc2    = iu2
   iuc3    = iu3
elseif (s2 <= s1 .and. s2 <= s3) then  ! m2 has lowest S value
   xemc(1)  = xem(im2)
   yemc(1)  = yem(im2)
   topmc(1) = topm(im2)

   xemc(2)  = xem(im3)
   yemc(2)  = yem(im3)
   topmc(2) = topm(im3)

   xemc(3)  = xem(im1)
   yemc(3)  = yem(im1)
   topmc(3) = topm(im1)

   sc1     = s2
   sc2     = s3
   sc3     = s1

   iuc1    = iu2
   iuc2    = iu3
   iuc3    = iu1
elseif (s3 <= s1 .and. s3 <= s2) then  ! m3 has lowest S value
   xemc(1)  = xem(im3)
   yemc(1)  = yem(im3)
   topmc(1) = topm(im3)

   xemc(2)  = xem(im1)
   yemc(2)  = yem(im1)
   topmc(2) = topm(im1)

   xemc(3)  = xem(im2)
   yemc(3)  = yem(im2)
   topmc(3) = topm(im2)

   sc1     = s3
   sc2     = s1
   sc3     = s2

   iuc1    = iu3
   iuc2    = iu1
   iuc3    = iu2
endif
   
! Find two points of intersection between current IW triangle and X slab

if (sc2 > op%slabloc(iplt)) then

   wt2 = (op%slabloc(iplt) - sc1) / (sc2 - sc1)
   wt1 = 1. - wt2

   x2    = wt1 * xemc(1)  + wt2 * xemc(2)  ! x coord of htpn(2)
   y2    = wt1 * yemc(1)  + wt2 * yemc(2)  ! y coord of htpn(2)
   topo2 = wt1 * topmc(1) + wt2 * topmc(2) ! topo height(2)

   htpn(2) = (x2 - nl%plotcoord1(iplt)) * sinvaz  &
           - (y2 - nl%plotcoord2(iplt)) * cosvaz
           
   ju2 = iuc3
   
   if (sc3 > op%slabloc(iplt)) then

      wt2 = (op%slabloc(iplt) - sc1) / (sc3 - sc1)
      wt1 = 1. - wt2

      x1    = wt1 * xemc(1)  + wt2 * xemc(3)  ! x coord of htpn(1)
      y1    = wt1 * yemc(1)  + wt2 * yemc(3)  ! y coord of htpn(1)
      topo1 = wt1 * topmc(1) + wt2 * topmc(3) ! topo height(1)

      htpn(1) = (x1 - nl%plotcoord1(iplt)) * sinvaz  &
              - (y1 - nl%plotcoord2(iplt)) * cosvaz
   
      ju1 = iuc2
      
   else

      wt2 = (op%slabloc(iplt) - sc3) / (sc2 - sc3)
      wt1 = 1. - wt2

      x1    = wt1 * xemc(3)  + wt2 * xemc(2)  ! x coord of htpn(1)
      y1    = wt1 * yemc(3)  + wt2 * yemc(2)  ! y coord of htpn(1)
      topo1 = wt1 * topmc(3) + wt2 * topmc(2) ! topo height(1)

      htpn(1) = (x1 - nl%plotcoord1(iplt)) * sinvaz  &
              - (y1 - nl%plotcoord2(iplt)) * cosvaz
   
      ju1 = iuc1
      
   endif

else

   wt2 = (op%slabloc(iplt) - sc2) / (sc3 - sc2)
   wt1 = 1. - wt2

   x2    = wt1 * xemc(2)  + wt2 * xemc(3)  ! x coord of htpn(2)
   y2    = wt1 * yemc(2)  + wt2 * yemc(3)  ! y coord of htpn(2)
   topo2 = wt1 * topmc(2) + wt2 * topmc(3) ! topo height(2)

   htpn(2) = (x2 - nl%plotcoord1(iplt)) * sinvaz  &
           - (y2 - nl%plotcoord2(iplt)) * cosvaz
   
   ju2 = iuc1

   wt2 = (op%slabloc(iplt) - sc1) / (sc3 - sc1)
   wt1 = 1. - wt2

   x1    = wt1 * xemc(1)  + wt2 * xemc(3)  ! x coord of htpn(1)
   y1    = wt1 * yemc(1)  + wt2 * yemc(3)  ! y coord of htpn(1)
   topo1 = wt1 * topmc(1) + wt2 * topmc(3) ! topo height(1)

   htpn(1) = (x1 - nl%plotcoord1(iplt)) * sinvaz  &
           - (y1 - nl%plotcoord2(iplt)) * cosvaz

   ju1 = iuc2

endif

htpn(3) = htpn(2)
htpn(4) = htpn(1)

return
end subroutine xyplot_w

