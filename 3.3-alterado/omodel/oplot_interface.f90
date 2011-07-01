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
subroutine oplot_init()

use oplot_coms, only: op
use plotcolors, only: gks_colors
use misc_coms,  only: io6

implicit none

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAM_OPLOT_INIT_IN,rst_buf)
#endif

call o_opngks()
call gks_colors(1)

! Initialize some oplot parameters (not set in namelist).

op%loopplot = 0

op%dualpts  = 'PROG'
op%icigrnd = 13
op%iplotback = 0

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_OPLOT_INIT_OUT,rst_buf)
#endif

return
end subroutine oplot_init

!===============================================================================

subroutine plot_fields()

use misc_coms,  only: io6, iyear1, imonth1, idate1, itime1, time_istp8
use oplot_coms, only: op
use mem_grid,   only: mza
use misc_coms,  only: runtype

implicit none

integer :: iplt,labincx,labincy,notavail,ko
real :: delxmin,fldval
integer :: outyear,outmonth,outdate,outhour
real, save :: dummy(1)=0.,xinc,yinc
character(len=30) :: ylabel
real, external :: walltime

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAM_PLOT_FIELDS_IN,rst_buf)
#endif

! Reopen the current graphics output workstation if it is closed
call o_reopnwk()

do iplt = 1,op%nplt

! Plot a white background for this frame

   if (op%iplotback /= 1) then
      call plotback()
   endif

! Draw filled-in map if so specified

   if ((op%projectn(iplt) == 'L'  .or.   &
        op%projectn(iplt) == 'P'  .or.   & 
        op%projectn(iplt) == 'O') .and.  & 
       (op%maptyp(iplt)   == 'M'  .and.  &
        op%contrtyp(iplt) /= 'T'  .and.  &
        op%contrtyp(iplt) /= 'F'  .and.  &
        op%contrtyp(iplt) /= 'O')) then
!   call mkmap(iplt,'FILL')
   endif

! Get units and stagpoint information for this field

   call oplot_lib(iplt,1,2,'UNITS',trim(op%fldname(iplt)),fldval,notavail,ko)
   
! If notavail = 3 is returned from above call, current plot field is unavailable

   if (notavail > 2) then
      write(io6,*) 'FIELD ',trim(op%fldname(iplt)),' NOT AVAILABLE IN THIS RUN.'
      call oplot_set(iplt)  ! can remove this?
      call o_set (op%hp1,op%hp2,op%vp1,op%vp2,0.,1.,0.,1.,1)
      call o_plchhq(op%fx1,op%fnamey  &
         ,trim(op%fldname(iplt))//' NOT AVAILABLE IN THIS RUN'  &
         ,.012 * (op%hp2 - op%hp1),0.,-1.)
      go to 10
   endif 
   
! Plot current field values with tiles or filled contours

   if (op%contrtyp(iplt) == 'T' .or. op%contrtyp(iplt) == 'F' .or.  &
       op%contrtyp(iplt) == 'O') then

      op%ifill = 1
      call slab(iplt)
   endif

! Plot current field values with contour lines

   if (op%contrtyp(iplt) == 'L' .or. op%contrtyp(iplt) == 'O') then
      op%ifill = 0
      call slab(iplt)
   endif

! Plot grid cell indices if so specified

   if (op%pltindx(iplt) == 'I' .or. op%pltindx(iplt) == 'J') then
      call plot_index(iplt)
   endif

! Plot current field values with printed numbers if so specified

   if (op%prtval(iplt) == 'P') then
      call slab_val(iplt)
   endif

! Plot T-point or U-point vectors if so specified

   call vectslab(iplt)

! SPECIAL HURRICANE FRANCES LOCATION PLOT

!   call plot_frances(iplt)

! Draw land/sea cell boundaries if so specified

   if (op%pltgrid_landsea(iplt) == 'g') then
      call plot_grid_landsea(iplt)
   endif

! Draw dual grid cell boundaries if so specified

   if (op%pltdualgrid(iplt) == 'D' .or. op%pltdualgrid(iplt) == 'd') then
      call plot_dualgrid(iplt)
   endif

! Draw grid cell boundaries if so specified

   if (op%pltgrid(iplt) == 'G') then
      call plot_grid(iplt)
   endif

! Draw outer grid (frame) boundary if so specified

   if (op%pltborder(iplt) == 'b') then
      call plot_grid_frame()
   endif

! Draw plot frame, tick marks, X and Y tick labels, X and Y axis labels
!   if so specified.

   if (op%pltborder(iplt) == 't') then

      if (op%projectn(iplt) == 'L') then
         call niceinc20(op%xmin,op%xmax,xinc,labincx)
         call niceinc20(op%ymin,op%ymax,yinc,labincy)
         
         if (op%windowin(iplt) /= 'W') then
            xinc = 10.
            yinc = 10.
            labincx = 3
            labincy = 3
         endif

         call oplot_xy2(op%panel(iplt),op%colorbar(iplt),0.,.016  &
            ,1,dummy,dummy                                       &
            ,'LONGITUDE (deg)','LATITUDE (deg)'                  &
            ,op%xmin,op%xmax,xinc,labincx                        &
            ,op%ymin,op%ymax,yinc,labincy                        )
      else
         call niceinc20(.001*op%xmin,.001*op%xmax,xinc,labincx)
         call niceinc20(.001*op%ymin,.001*op%ymax,yinc,labincy)

         if (op%projectn(iplt) == 'C' .or. op%projectn(iplt) == 'V') then
            ylabel = 'Z (km)'
         else
            ylabel = 'Y (km)'
         endif

         call oplot_xy2(op%panel(iplt),op%colorbar(iplt),1.0,.016   &
            ,1,dummy,dummy                                         &
            ,'X (km)',ylabel                                       &
            ,.001*op%xmin,.001*op%xmax,xinc,labincx                &
            ,.001*op%ymin,.001*op%ymax,yinc,labincy                )
      endif

   endif

! Draw line map if so specified

   if ((op%projectn(iplt) == 'L'  .or.   &
        op%projectn(iplt) == 'P'  .or.   & 
        op%projectn(iplt) == 'O') .and.  & 
       (op%maptyp(iplt)   /= 'N') .and.  &
       (op%maptyp(iplt)   == 'm'  .or.   &
        op%contrtyp(iplt) == 'T'  .or.   &
        op%contrtyp(iplt) == 'F'  .or.   &
        op%contrtyp(iplt) == 'O')) then

      call mkmap(iplt,'LINE')

   endif

! Draw colorbar if so specified

   if (op%colorbar(iplt) == 'c') call plot_colorbar(iplt,op%icolortab(iplt))

! Draw label bar if so specified

   if (op%labelbar(iplt) == 'n' .or. op%labelbar(iplt) == 'i') then
      call date_add_to8 (iyear1,imonth1,idate1,itime1*100  &
         ,time_istp8,'s',outyear,outmonth,outdate,outhour)

      call plot_labelbar (iplt,op%label  &
         ,op%units,op%projectn(iplt),op%slabloc(iplt)      &
         ,op%fldval_min,op%fldval_max  &
         ,time_istp8,outhour,outdate,outmonth,outyear)
   endif

10 continue

! Complete current frame if so specified

   if (op%frameoff(iplt) /= 'f') then
      call o_frame()
   endif

enddo

! Close the current workstation if not a plotonly run and if output
! is to a NCAR graphics meta file. This allows viewing the complete
! meta file (including the last frame) during a run and in case the
! simulation crashes.

if ((trim(runtype) /= 'PLOTONLY') .and. (op%plttype == 0)) then


   call o_clswk()

endif
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_PLOT_FIELDS_OUT,rst_buf)
#endif

return
end subroutine plot_fields

!===============================================================================

subroutine slab(iplt)

use oplot_coms, only: op
use mem_grid,   only: mwa
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt

integer :: iw
real :: opvsprd

call oplot_set(iplt)

! Set plot line color (black)

call o_gsplci(10)
call o_gsfaci(10)
call o_gstxci(10)
call o_sflush()

if (op%projectn(iplt) == 'L' .or.  &
    op%projectn(iplt) == 'P' .or.  &
    op%projectn(iplt) == 'O' .or.  &
    op%projectn(iplt) == 'Z') then  ! For horizontal plotting...

   if (op%contrtyp(iplt) /= 'F' .and. op%contrtyp(iplt) /= 'L' .and.  &
       op%contrtyp(iplt) /= 'O') then

      if (op%stagpt == 'T' .or. op%stagpt == 'W') then
         call tileslab_htw(iplt,'T')
      elseif (op%stagpt == 'M' .or. op%stagpt == 'P') then
         call tileslab_hmp(iplt,'T')
      elseif (op%stagpt == 'U') then
         call tileslab_hu(iplt,'T')
      elseif (op%stagpt == 'S') then  ! for sea cells
         call tileslab_hs(iplt,'T')
      elseif (op%stagpt == 'L') then  ! for land cells
         call tileslab_hl(iplt,'T') 
      elseif (op%stagpt == 'B') then  ! for both sea and land cells
         op%stagpt = 'S'
         call tileslab_hs(iplt,'T')
         op%stagpt = 'L'
         call tileslab_hl(iplt,'T') 
         op%stagpt = 'B'
      elseif (op%stagpt == 'F') then  ! for flux cells (land + sea)
         op%stagpt = 'S'
         call tileslab_hfs(iplt,'T')
         op%stagpt = 'L'
         call tileslab_hfl(iplt,'T') 
         op%stagpt = 'F'
      endif

   else

      if     (op%stagpt == 'T' .or. op%stagpt == 'W') then
         call contslab_htw(iplt)
      elseif (op%stagpt == 'U' .or. op%stagpt == 'N') then
         call contslab_hnu(iplt)
      elseif (op%stagpt == 'M' .or. op%stagpt == 'P') then
         call contslab_hmp(iplt)
      endif

   endif

else  ! Vertical plots

   if (op%contrtyp(iplt) /= 'F' .and. op%contrtyp(iplt) /= 'L' .and.  &
       op%contrtyp(iplt) /= 'O') then

      if (op%stagpt == 'T') then
      	 call tileslab_vt(iplt,'T')
      elseif (op%stagpt == 'U') then
      	 call tileslab_vu(iplt,'T')
      elseif (op%stagpt == 'W') then
      	 call tileslab_vw(iplt,'T')
      elseif (op%stagpt == 'A') then  ! combined atmos, sea, and land vals

! Initialize cell location scratch array VT2DA

         do iw = 2,mwa
!            vt2da(iw) = xu(itab_w(iw)%iu1)  ! NEEDS WORK
         enddo

         op%stagpt = 'T'
       	 call tileslab_vt(iplt,'T')
         call tileslab_vs(iplt,'T')
         call tileslab_vl(iplt,'T') 

      endif

   else

      if (op%stagpt == 'T') then
         call contslab_vt(iplt)
      elseif (op%stagpt == 'U') then
         call contslab_vu(iplt)
      elseif (op%stagpt == 'W') then
         call contslab_vw(iplt)
      endif

   endif

endif   ! needs elseif sections for 'X' and 'Y'

return
end subroutine slab

!===============================================================================

subroutine slab_val(iplt)

use oplot_coms, only: op
use mem_grid,   only: mwa
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt

integer :: iw
real :: opvsprd

call oplot_set(iplt)

if (op%projectn(iplt) == 'L' .or.  &
    op%projectn(iplt) == 'P' .or.  &
    op%projectn(iplt) == 'O' .or.  &
    op%projectn(iplt) == 'Z') then  ! For horizontal plotting...

      if (op%stagpt == 'T' .or. op%stagpt == 'W') then
         call tileslab_htw(iplt,'P')
      elseif (op%stagpt == 'M' .or. op%stagpt == 'P') then
         call tileslab_hmp(iplt,'P')
      elseif (op%stagpt == 'U') then
         call tileslab_hu(iplt,'P')
      elseif (op%stagpt == 'S') then  ! for sea cells
         call tileslab_hs(iplt,'P')
      elseif (op%stagpt == 'L') then  ! for land cells
         call tileslab_hl(iplt,'P') 
      elseif (op%stagpt == 'B') then  ! for both sea and land cells
         op%stagpt = 'S'
         call tileslab_hs(iplt,'P')
         op%stagpt = 'L'
         call tileslab_hl(iplt,'P') 
         op%stagpt = 'B'
      elseif (op%stagpt == 'F') then  ! for flux cells (land + sea)
         op%stagpt = 'S'
         call tileslab_hfs(iplt,'P')
         op%stagpt = 'L'
         call tileslab_hfl(iplt,'P') 
         op%stagpt = 'F'
      endif

else  ! Vertical plots

   if (op%stagpt == 'T') then
      call tileslab_vt(iplt,'P')
   elseif (op%stagpt == 'U') then
      call tileslab_vu(iplt,'P')
   elseif (op%stagpt == 'W') then
      call tileslab_vw(iplt,'P')
   elseif (op%stagpt == 'A') then  ! combined atmos, sea, and land vals

! Initialize cell location scratch array VT2DA

      do iw = 2,mwa
!         vt2da(iw) = xu(itab_w(iw)%iu1)  ! NEEDS WORK
      enddo

      op%stagpt = 'T'
      call tileslab_vt(iplt,'P')
      call tileslab_vs(iplt,'P')
      call tileslab_vl(iplt,'P') 

   endif

endif   ! needs elseif sections for 'X' and 'Y'

return
end subroutine slab_val

!===============================================================================

subroutine plot_index(iplt)

use oplot_coms, only: op
use mem_grid,   only: mma, mua, mwa, xem, yem, zem, xeu, yeu, zeu,  &
                      xew, yew, zew
use sea_coms,   only: mws, mus
use mem_sea,    only: itab_ws, itab_us, sea
use leaf_coms,  only: mwl, mul
use mem_leaf,   only: itab_wl, itab_ul, land
use misc_coms,  only: io6
use mem_sflux,  only: mlandflux, landflux, mseaflux, seaflux,  &
                      xemltrap, yemltrap, zemltrap,  &
                      xemstrap, yemstrap, zemstrap

implicit none

integer, intent(in) :: iplt

integer :: im,iu,iw,iwl,iws,jm,ilf,isf,jtrap,itrap,ktrap
real :: hpt,vpt
real :: xewf1, yewf1, zewf1
real :: xewf2, yewf2, zewf2

! This subroutine plots grid point indices for staggered points not covered
! by the current field
call oplot_set(iplt)  ! not needed here?

if (op%projectn(iplt) /= 'L' .and.  &
    op%projectn(iplt) /= 'P' .and.  &
    op%projectn(iplt) /= 'O' .and.  &
    op%projectn(iplt) /= 'Z') return

! Plot M point indices

   call o_sflush()
   call o_gsplci(10)
   call o_gstxci(10)
   call o_sflush()

if (op%pltindx(iplt) == 'J' .or.  &
     trim(op%stagpt) == 'M' .or.  &
     trim(op%stagpt) == 'P') then
   do im = 2,mma
      call oplot_transform(iplt,xem(im),yem(im),zem(im),hpt,vpt)
      if (hpt < op%xmin .or. hpt > op%xmax .or.  &
          vpt < op%ymin .or. vpt > op%ymax) cycle
      call oplot_locindex(im,hpt,vpt,op%psiz)
   enddo
endif

! Plot U point indices

if (op%pltindx(iplt) == 'J' .or.  &
     trim(op%stagpt) == 'U') then
   do iu = 2,mua
      call oplot_transform(iplt,xeu(iu),yeu(iu),zeu(iu),hpt,vpt)
      if (hpt < op%xmin .or. hpt > op%xmax .or.  &
          vpt < op%ymin .or. vpt > op%ymax) cycle
      call oplot_locindex(iu,hpt,vpt,op%psiz) 
   enddo   
endif   

! Plot W point indices

if (op%pltindx(iplt) == 'J' .or.  &
     trim(op%stagpt) == 'W' .or.  &
     trim(op%stagpt) == 'T') then
   do iw = 2,mwa
      call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),hpt,vpt)
      if (hpt < op%xmin .or. hpt > op%xmax .or.  &
          vpt < op%ymin .or. vpt > op%ymax) cycle
      call oplot_locindex(iw,hpt,vpt,op%psiz) 
   enddo   
endif   

! Plot WL point indices

if (trim(op%stagpt) == 'L' .or.  &
    trim(op%stagpt) == 'B') then

   call o_sflush()
   call o_gsplci(6)
   call o_gstxci(6)
   call o_sflush()

   do iwl = 2,mwl

      call oplot_transform(iplt,  &
         land%xewl(iwl),land%yewl(iwl),land%zewl(iwl),hpt,vpt)

      if (hpt < op%xmin .or. hpt > op%xmax .or.  &
          vpt < op%ymin .or. vpt > op%ymax) cycle
      call oplot_locindex(iwl,hpt,vpt,.7 * op%psiz) 
   enddo

   call o_sflush()

endif   

! Plot WS point indices

if (trim(op%stagpt) == 'S' .or.  &
    trim(op%stagpt) == 'B') then

   call o_sflush()
   call o_gsplci(3)
   call o_gstxci(3)
   call o_sflush()

   do iws = 2,mws

      call oplot_transform(iplt,  &
         sea%xews(iws),sea%yews(iws),sea%zews(iws),hpt,vpt)

      if (hpt < op%xmin .or. hpt > op%xmax .or.  &
          vpt < op%ymin .or. vpt > op%ymax) cycle
      call oplot_locindex(iws,hpt,vpt,.7 * op%psiz) 
   enddo

   call o_sflush()

endif   

! Plot Flux cell indices

if (trim(op%stagpt) == 'F') then

! LANDFLUX CELLS

   do ilf = 2,mlandflux

      call oplot_transform(iplt,landflux(ilf)%xef,  &
                                landflux(ilf)%yef,  &
                                landflux(ilf)%zef,hpt,vpt)

      if (hpt < op%xmin .or. hpt > op%xmax .or.  &
          vpt < op%ymin .or. vpt > op%ymax) cycle

      call oplot_locindex(ilf,hpt,vpt,.7*op%psiz) 

   enddo

! SEAFLUX CELLS

   do isf = 2,mseaflux

      call oplot_transform(iplt,seaflux(isf)%xef,  &
                                seaflux(isf)%yef,  &
                                seaflux(isf)%zef,hpt,vpt)

      if (hpt < op%xmin .or. hpt > op%xmax .or.  &
          vpt < op%ymin .or. vpt > op%ymax) cycle

      call oplot_locindex(isf,hpt,vpt,.7*op%psiz) 

   enddo

endif   

return
end subroutine plot_index

!===============================================================================

subroutine vectslab(iplt)

use oplot_coms, only: op
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt

call oplot_set(iplt)

call o_gsplci(10)
call o_gstxci(10)

if (op%projectn(iplt) == 'L' .or.  &
    op%projectn(iplt) == 'P' .or.  &
    op%projectn(iplt) == 'O' .or.  &
    op%projectn(iplt) == 'Z') then

   if (op%vectbarb(iplt) == 'U' .or. &
       op%vectbarb(iplt) == 'u') call vecuslab_h(iplt)

   if (op%vectbarb(iplt) == 'V') call vectslab_h(iplt)
   
elseif (op%projectn(iplt) == 'X') then  ! etc for 'Y' and 'C'...


endif   

return
end subroutine vectslab

!===============================================================================

subroutine plot_underground(iplt,kt)

use oplot_coms, only: op
use mem_grid,   only: mwa, mza, zm, zt, lpw, xem, yem, zem, xew, yew, zew
use mem_ijtabs, only: itab_w
use misc_coms,  only: io6
use mem_basic,  only: press

implicit none

integer, intent(in) :: iplt,kt

real :: xpt,ypt,plev
integer :: k
integer :: ntpn,ifill,iw,itpn,iu1,iu2,iv1,iv2,im1,im2,im3,iok,ju1,ju2
real, dimension(8) :: htpn,vtpn
real :: topo1,topo2

! Set vertical index for points that may be interpolated vertically to pressure level

k = kt

if (op%projectn(iplt) == 'L' .or.  &
    op%projectn(iplt) == 'P' .or.  &
    op%projectn(iplt) == 'O' .or.  &
    op%projectn(iplt) == 'Z') then  ! Horizontal cross section

   do iw = 2,mwa

! Reset k value if plotting on pressure level

      if (op%projectn(iplt) /= 'Z' .and.  &
          op%pltlev(iplt) == 'p') then

         plev = op%slabloc(iplt) * 100.  ! convert from mb to Pa

         k = lpw(iw) - 1

         do while (press(k+1,iw) > plev .and. k < mza-1)
            k = k + 1
         enddo
         
      endif

! Skip this iw point if k is not underground

      if (k >= lpw(iw)) cycle

! Get tile plot coordinates

      if (op%projectn(iplt) == 'L') then    ! For determining wrap-around direction
         call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),xpt,ypt)
      endif
      
      im1 = itab_w(iw)%im1
      im2 = itab_w(iw)%im2
      im3 = itab_w(iw)%im3

      call oplot_transform(iplt,xem(im1),yem(im1),zem(im1),htpn(1),vtpn(1))
      call oplot_transform(iplt,xem(im2),yem(im2),zem(im2),htpn(2),vtpn(2))
      call oplot_transform(iplt,xem(im3),yem(im3),zem(im3),htpn(3),vtpn(3))

      do itpn = 1,3

! Avoid wrap-around

         if (op%projectn(iplt) == 'L') call ll_unwrap(xpt,htpn(itpn))
         
      enddo
 
! Jump out of loop if any cell corner is on other side of earth

   if (any(htpn(1:3) > 1.e11)) cycle

! Skip current iw cell if entire cell is outside plot window.

      if ( all(htpn(1:3) < op%xmin) .or. all(htpn(1:3) > op%xmax) .or.  &
           all(vtpn(1:3) < op%ymin) .or. all(vtpn(1:3) > op%ymax) ) cycle

! Tile-plot cell with underground color

      call fillpolyg(3,htpn,vtpn,op%icigrnd)

   enddo

else  ! Vertical cross section

   do iw = 2,mwa

! Get horizontal plot coordinates for cells in this column

      if (op%projectn(iplt) == 'C') then
         call coneplot_w(iw,topo1,topo2,iok,htpn)
      elseif (op%projectn(iplt) == 'X') then
         call xplot_w(iplt,iw,iok,htpn)
      elseif (op%projectn(iplt) == 'Y') then
         call yplot_w(iplt,iw,iok,htpn)
      elseif (op%projectn(iplt) == 'V') then
         call xyplot_w(iplt,iw,ju1,ju2,topo1,topo2,iok,htpn)
      endif

! Jump out of loop if this W point does not intersect plot cone

      if (iok /= 1) go to 9

! Jump out of loop if either cell side is outside plot window. 

      if (htpn(1) < op%xmin .or. htpn(1) > op%xmax .or.  &
          htpn(2) < op%xmin .or. htpn(2) > op%xmax) go to 9
   
      do k = 2,mza-1

! If cell is totally underground...

         if (k < lpw(iw)) then 

! Get T-cell vertical coordinates

            vtpn(1) = zm(k-1)
            vtpn(2) = vtpn(1)
            vtpn(3) = zm(k)
            vtpn(4) = vtpn(3)

! Plot cell with underground color

            call fillpolyg(4,htpn(1),vtpn(1),op%icigrnd)

! If cell is partially underground...

         elseif (k == lpw(iw)) then 

! Get T-cell vertical coordinates

            vtpn(1) = zm(k-1)
            vtpn(2) = vtpn(1)
            vtpn(3) = topo2
            vtpn(4) = topo1

            call fillpolyg(4,htpn(1),vtpn(1),op%icigrnd)


         endif

      enddo

9     continue

   enddo

endif

return
end subroutine plot_underground

!===============================================================================

subroutine plot_underground_m(iplt,kt)

use oplot_coms, only: op
use mem_grid,   only: mma, mza, zm, zt, lpw, xem, yem, zem, xew, yew, zew
use mem_ijtabs, only: itab_m
use misc_coms,  only: io6
use mem_basic,  only: press

implicit none

integer, intent(in) :: iplt,kt

real :: xpt,ypt,plev
integer :: ntpn,ifill,iw,itpn,iu1,iu2,iv1,iv2,k,im,iok,ju1,ju2
real, dimension(8) :: htpn,vtpn
integer :: kw(8)

kw(1:8) = kt

if (op%projectn(iplt) == 'L' .or.  &
    op%projectn(iplt) == 'P' .or.  &
    op%projectn(iplt) == 'O' .or.  &
    op%projectn(iplt) == 'Z') then  ! Horizontal cross section

   do im = 2,mma

!     Get tile plot coordinates

      if (op%projectn(iplt) == 'L') then    ! For determining wrap-around direction
         call oplot_transform(iplt,xem(im),yem(im),zem(im),xpt,ypt)
      endif
      
      ntpn = itab_m(im)%ntpn
      do itpn = 1,ntpn

         iw = itab_m(im)%iw(itpn)

! Reset kw values if plotting on pressure level

         if (op%projectn(iplt) /= 'Z' .and.  &
             op%pltlev(iplt) == 'p') then

            plev = op%slabloc(iplt) * 100.  ! convert from mb to Pa

            kw(itpn) = lpw(iw) - 1

            do while (press(kw(itpn)+1,iw) > plev .and. kw(itpn) < mza-1)
               kw(itpn) = kw(itpn) + 1
            enddo
         
         endif

         call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),htpn(itpn),vtpn(itpn))

         ! Avoid wrap-around
         if (op%projectn(iplt) == 'L') call ll_unwrap(xpt,htpn(itpn))
      enddo

! Jump out of loop if any cell corner is on other side of earth

      if (any(htpn(1:ntpn) > 1.e11)) cycle

!     Jump out of loop if entire cell is outside plot window.

      if ( all(htpn(1:ntpn) < op%xmin) .or. all(htpn(1:ntpn) > op%xmax) .or.  &
           all(vtpn(1:ntpn) < op%ymin) .or. all(vtpn(1:ntpn) > op%ymax) ) cycle
      
!     If cell is underground, tile-plot it with underground color

      if (any(kw(1:ntpn) < lpw(itab_m(im)%iw(1:ntpn)))) then
         call fillpolyg(ntpn,htpn,vtpn,op%icigrnd)
      endif

   enddo

endif
return
end subroutine plot_underground_m

!===============================================================================

subroutine plot_grid(iplt)

use oplot_coms, only: op
use mem_grid, only: mwa, mua, mza, xem, yem, zem, xew, yew, zew, zm
use mem_ijtabs, only: itab_w, itab_u
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt

integer :: iflag180,iskip

real :: xpt,ypt
integer :: ntpn,ifill,iw,itpn,iu1,iu2,iv1,iv2,k,im1,im2,im3,iok,ju1,ju2,iu
real, dimension(4) :: htpn,vtpn

real :: xp1,xp2,yp1,yp2
real :: xq1,xq2,yq1,yq2
real :: topo1,topo2

!call oplot_set ()

call o_sflush()
call o_gsplci(10)
call o_gstxci(10)

if (op%projectn(iplt) == 'L'  .or.  &
    op%projectn(iplt) == 'P'  .or.  &
    op%projectn(iplt) == 'O'  .or.  &
    op%projectn(iplt) == 'Z') then   ! Horizontal cross section

   do iu = 2,mua

! Get tile plot coordinates.

      im1 = itab_u(iu)%im1
      im2 = itab_u(iu)%im2

      call oplot_transform(iplt,xem(im1),yem(im1),zem(im1),xp1,yp1)
      call oplot_transform(iplt,xem(im2),yem(im2),zem(im2),xp2,yp2)

! Avoid wrap-around and set iflag180

      iflag180 = 0

      if (op%projectn(iplt) == 'L') then
         call ll_unwrap(xp1,xp2)

         if (xp1 < -180. .or. xp2 < -180.) iflag180 =  1
         if (xp1 >  180. .or. xp2 >  180.) iflag180 = -1
      endif

! Truncate segment if it crosses plot window boundary or skip if both endpoints
! are outside plot window

      call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

      if (iskip == 1) cycle

! Plot line segment

      call o_frstpt (xq1,yq1)
      call o_vector (xq2,yq2)

! If this segment crosses +/- 180 degrees longitude in lat/lon plot, re-plot
! at other end

      if (iflag180 /= 0) then

         xp1 = xp1 + 360. * iflag180
         xp2 = xp2 + 360. * iflag180

         call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

! Plot line segment

         call o_frstpt (xq1,yq1)
         call o_vector (xq2,yq2)
      
      endif

   enddo

else  ! Vertical cross section

   do iw = 2,mwa

! Get horizontal plot coordinates for cells in this column

      if (op%projectn(iplt) == 'C') then
         call coneplot_w(iw,topo1,topo2,iok,htpn)
      elseif (op%projectn(iplt) == 'X') then
         call xplot_w(iplt,iw,iok,htpn)
      elseif (op%projectn(iplt) == 'Y') then
         call yplot_w(iplt,iw,iok,htpn)
      elseif (op%projectn(iplt) == 'V') then
         call xyplot_w(iplt,iw,ju1,ju2,topo1,topo2,iok,htpn)
      endif

! Jump out of loop if this W point does not intersect plot cone or plane

      if (iok /= 1) go to 9

! Jump out of loop if either cell side is outside plot window. 

      if (htpn(1) < op%xmin .or. htpn(1) > op%xmax .or.  &
          htpn(2) < op%xmin .or. htpn(2) > op%xmax) go to 9

      do k = 2,mza-1

! Get T-cell vertical coordinates

         vtpn(1) = zm(k-1)
         vtpn(2) = vtpn(1)
         vtpn(3) = zm(k)
         vtpn(4) = vtpn(3)

! Plot cell

         call o_frstpt(htpn(4),vtpn(4))
         do itpn = 1,4
            call o_vector(htpn(itpn),vtpn(itpn))
         enddo
      enddo

9  continue

   enddo

endif

return
end subroutine plot_grid

!===============================================================================

subroutine plot_dualgrid(iplt)

use oplot_coms,  only: op
use mem_grid,    only: mwa, mua, mza, xem, yem, zem, xew, yew, zew, zm,  &
                     wnx, wny, wnz
use mem_ijtabs,  only: itab_w, itab_u
use misc_coms,   only: io6
use consts_coms, only: erad

implicit none

integer, intent(in) :: iplt

integer :: iflag180,iskip

real :: xpt,ypt
integer :: ntpn,ifill,iw,itpn,iu1,iu2,iv1,iv2,k,im1,im2,im3,iok,ju1,ju2,iu,  &
           iw1,iw2
real, dimension(4) :: htpn,vtpn

real :: xp1,xp2,yp1,yp2
real :: xq1,xq2,yq1,yq2
real :: topo1,topo2

!call oplot_set ()

call o_sflush()
call o_gsplci(6)
call o_gstxci(6)
call o_sflush()

if (op%projectn(iplt) == 'L'  .or.  &
    op%projectn(iplt) == 'P'  .or.  &
    op%projectn(iplt) == 'O'  .or.  &
    op%projectn(iplt) == 'Z') then   ! Horizontal cross section

   do iu = 2,mua

! Get tile plot coordinates.

      iw1 = itab_u(iu)%iw1
      iw2 = itab_u(iu)%iw2
      
      if (iw1 < 2 .or. iw2 < 2) cycle

      if (op%pltdualgrid(iplt) == 'D') then

         call oplot_transform(iplt,xew(iw1),yew(iw1),zew(iw1),xp1,yp1)
         call oplot_transform(iplt,xew(iw2),yew(iw2),zew(iw2),xp2,yp2)

      else

         call oplot_transform(iplt,wnx(iw1)*erad,  &
                                   wny(iw1)*erad,  &
                                   wnz(iw1)*erad,xp1,yp1)
         call oplot_transform(iplt,wnx(iw2)*erad,  &
                                   wny(iw2)*erad,  &
                                   wnz(iw2)*erad,xp2,yp2)

      endif
   
! Avoid wrap-around and set iflag180

      iflag180 = 0

      if (op%projectn(iplt) == 'L') then
         call ll_unwrap(xp1,xp2)

         if (xp1 < -180. .or. xp2 < -180.) iflag180 =  1
         if (xp1 >  180. .or. xp2 >  180.) iflag180 = -1
      endif

! Truncate segment if it crosses plot window boundary or skip if both endpoints
! are outside plot window

      call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

      if (iskip == 1) cycle

! Plot line segment

      call o_frstpt (xq1,yq1)
      call o_vector (xq2,yq2)

! If this segment crosses +/- 180 degrees longitude in lat/lon plot, re-plot
! at other end

      if (iflag180 /= 0) then

         xp1 = xp1 + 360. * iflag180
         xp2 = xp2 + 360. * iflag180

         call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

! Plot line segment

         call o_frstpt (xq1,yq1)
         call o_vector (xq2,yq2)
      
      endif

   enddo

else  ! Vertical cross section

endif

call o_sflush()

return
end subroutine plot_dualgrid

!===============================================================================

subroutine plot_grid_landsea(iplt)

use oplot_coms, only: op
use sea_coms,   only: mws, mus
use mem_sea,    only: sea, itab_ws, itab_us
use leaf_coms,  only: mwl, mul
use mem_leaf,   only: land, itab_wl, itab_ul
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt

integer :: iflag180,iskip

integer :: ntpn,ifill,iw,itpn,iu1,iu2,iv1,iv2,k,im1,im2,im3,iok,ju1,ju2,iu,  &
           ius,iul

real :: xp1,xp2,yp1,yp2
real :: xq1,xq2,yq1,yq2
real :: xpt,ypt

real :: htpn(4),vtpn(4)

!call oplot_set ()

if (op%projectn(iplt) == 'L'  .or.  &
    op%projectn(iplt) == 'P'  .or.  &
    op%projectn(iplt) == 'O'  .or.  &
    op%projectn(iplt) == 'Z') then   ! Horizontal cross section

!-------------------------------------------------
! Plot sea grid
!-------------------------------------------------

   call o_sflush()
   call o_gsplci(3)
   call o_gstxci(3)

   do ius = 2,mus

! Get tile plot coordinates.

      im1 = itab_us(ius)%im1
      im2 = itab_us(ius)%im2

      call oplot_transform  &
         (iplt,sea%xems(im1),sea%yems(im1),sea%zems(im1),xp1,yp1)
      call oplot_transform  &
         (iplt,sea%xems(im2),sea%yems(im2),sea%zems(im2),xp2,yp2)

! Avoid wrap-around and set iflag180

      iflag180 = 0

      if (op%projectn(iplt) == 'L') then
         call ll_unwrap(xp1,xp2)

         if (xp1 < -180. .or. xp2 < -180.) iflag180 =  1
         if (xp1 >  180. .or. xp2 >  180.) iflag180 = -1
      endif

! Truncate segment if it crosses plot window boundary or skip if both endpoints
! are outside plot window

      call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

      if (iskip == 1) cycle

! Plot line segment

      call o_frstpt (xq1,yq1)
      call o_vector (xq2,yq2)

! If this segment crosses +/- 180 degrees longitude in lat/lon plot, re-plot
! at other end

      if (iflag180 /= 0) then

         xp1 = xp1 + 360. * iflag180
         xp2 = xp2 + 360. * iflag180

         call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

! Plot line segment

         call o_frstpt (xq1,yq1)
         call o_vector (xq2,yq2)
      
      endif

   enddo

!-------------------------------------------------
! Plot land grid
!-------------------------------------------------

   call o_sflush()
   call o_gsplci(6)
   call o_gstxci(6)

   do iul = 2,mul

! Get tile plot coordinates.

      im1 = itab_ul(iul)%im1
      im2 = itab_ul(iul)%im2

      call oplot_transform  &
         (iplt,land%xeml(im1),land%yeml(im1),land%zeml(im1),xp1,yp1)
      call oplot_transform  &
         (iplt,land%xeml(im2),land%yeml(im2),land%zeml(im2),xp2,yp2)

! Avoid wrap-around and set iflag180

      iflag180 = 0

      if (op%projectn(iplt) == 'L') then
         call ll_unwrap(xp1,xp2)

         if (xp1 < -180. .or. xp2 < -180.) iflag180 =  1
         if (xp1 >  180. .or. xp2 >  180.) iflag180 = -1
      endif

! Truncate segment if it crosses plot window boundary or skip if both endpoints
! are outside plot window

      call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

      if (iskip == 1) cycle

! Plot line segment

      call o_frstpt (xq1,yq1)
      call o_vector (xq2,yq2)

! If this segment crosses +/- 180 degrees longitude in lat/lon plot, re-plot
! at other end

      if (iflag180 /= 0) then

         xp1 = xp1 + 360. * iflag180
         xp2 = xp2 + 360. * iflag180

         call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

! Plot line segment

         call o_frstpt (xq1,yq1)
         call o_vector (xq2,yq2)
      
      endif

   enddo
   call o_sflush()

else  ! Vertical cross section

   write(io6,*) 'Land/sea grid plot not available in vertical cross section.'

endif

return
end subroutine plot_grid_landsea

!===============================================================================

subroutine plot_grid_frame()

use oplot_coms, only: op
use misc_coms,  only: io6

implicit none

!call oplot_set ()

call o_sflush()
call o_gsplci(10)
call o_gstxci(10)

call o_frstpt(op%xmin,op%ymin)
call o_vector(op%xmax,op%ymin)
call o_vector(op%xmax,op%ymax)
call o_vector(op%xmin,op%ymax)
call o_vector(op%xmin,op%ymin)
call o_sflush()

return
end subroutine plot_grid_frame

!===============================================================================

subroutine mkmap(iplt,filltype)

use oplot_coms,  only: op
use consts_coms, only: erad, erad2
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt
character(len=*), intent(in) :: filltype

! Set plot color (black)

call o_gsplci(10)
call o_gsfaci(10)
call o_gstxci(10)
call o_sflush()

call oplot_set(iplt)

if (trim(filltype) == 'LINE') then

   call o_mapint()
   call o_mappos(op%h1,op%h2,op%v1,op%v2)
   call o_mapsti('GR',10)    ! For 10-degree spacing of plotted parallels and meridians
   call o_mapsti('DA',65535) ! To plot parallels and meridians with solid lines
   call o_mapstc('OU','PS')  ! To plot bnds of continents, countries, and US  states
!   call o_mapstc('OU','NO')  ! To plot NO geographic information

   if (op%projectn(iplt) == 'L') then
      call o_maproj('CE',0.,op%plon3,0.)  ! Force plat to be zero for CE projection
      call o_mapset('LI',op%xmin,op%xmax,op%ymin,op%ymax)
   elseif (op%projectn(iplt) == 'P') then
      call o_maproj('ST',op%plat3,op%plon3,0.)
      call o_mapset('LI',op%xmin/erad2,op%xmax/erad2  &
                      ,op%ymin/erad2,op%ymax/erad2)
   elseif (op%projectn(iplt) == 'O') then
      call o_maproj('OR',op%plat3,op%plon3,0.)
      call o_mapset('LI',op%xmin/erad,op%xmax/erad  &
                      ,op%ymin/erad,op%ymax/erad)
   endif

   call o_mapint()             ! Initialize the above parameters
!   call mapgrd()
   call o_maplot()
   call o_sflush()

endif

return
end subroutine mkmap

!===============================================================================

subroutine psub()

use oplot_coms, only: op
use misc_coms,  only: io6

implicit none

if (op%loopplot == 1) then

! Reopen the current graphics output workstation if it is closed

   call o_reopnwk()

   call plotback()  ! White background
   call oplot_set(1)  ! Use projection of FIRST plot field in OLAMIN 
endif

return
end subroutine psub

!===============================================================================

subroutine qsub(var,i)

use oplot_coms, only: op
use mem_grid,   only: xem, yem, zem, xeu, yeu, zeu, xew, yew, zew
use misc_coms,  only: io6

implicit none

integer, intent(in) :: i
character(len=*), intent(in) :: var

real :: hpt,vpt
integer, save :: jplt=1  ! transform to projection of FIRST field in OLAMIN

if (op%loopplot /= 1) return

if (var == 'M') then

   call oplot_transform(jplt,xem(i),yem(i),zem(i),hpt,vpt)
   if (hpt > op%xmin .and. hpt < op%xmax .and.  &
       vpt > op%ymin .and. vpt < op%ymax) then
      call oplot_locindex(i,hpt,vpt,op%psiz)
   endif

elseif (var == 'U') then

   call oplot_transform(jplt,xeu(i),yeu(i),zeu(i),hpt,vpt)
   
   if (hpt > op%xmin .and. hpt < op%xmax .and.  &
       vpt > op%ymin .and. vpt < op%ymax) then
      call oplot_locindex(i,hpt,vpt,op%psiz)
   endif

elseif (var == 'W') then

   call oplot_transform(jplt,xew(i),yew(i),zew(i),hpt,vpt)
   if (hpt > op%xmin .and. hpt < op%xmax .and.  &
       vpt > op%ymin .and. vpt < op%ymax) then
      call oplot_locindex(i,hpt,vpt,op%psiz)
   endif

endif

! Add new sections for 'S' and 'L'  (and 'SF', and 'LF' ??)

! For this, need xeland, yeland, zeland, xesea, yesea, zesea
!    (and xeseaflux, etc.)

return
end subroutine qsub

!===============================================================================

subroutine rsub(fldname,num)

use misc_coms,  only: io6, iyear1, imonth1, idate1, itime1, time_istp8
use oplot_coms, only: op

implicit none

character(len=*), intent(in) :: fldname
integer, intent(in) :: num

integer :: outyear,outmonth,outdate,outhour
character :: title*20

if (op%loopplot == 1) then
   write (title,'(i2)') num
   title = fldname//trim(title)

   call plot_grid(1)  ! Uses projection from FIRST plot field in OLAMIN

   call date_add_to8 (iyear1,imonth1,idate1,itime1*100  &
      ,time_istp8,'s',outyear,outmonth,outdate,outhour)
      
   call plot_labelbar (1,trim(title),' ','H',0.,0.,0.  &
      ,time_istp8,outhour,outdate,outmonth,outyear)

   call o_frame()

! Close the current workstation if not a plotonly run and if output
! is to a NCAR graphics meta file. This allows viewing the complete
! meta file (including the last frame) during a run and in case the
! simulation crashes.

   call o_clswk()

endif

return
end subroutine rsub

!===============================================================================

subroutine oplot_transform(iplt,xeq,yeq,zeq,xout,yout)

use oplot_coms, only: op
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt
real, intent(in) :: xeq,yeq,zeq
real, intent(out) :: xout,yout

! Transform to plot frame coordinates

if (op%projectn(iplt) == 'L') then
   call e_ec(xeq,yeq,zeq,xout,yout)
! xout is now true geographic longitude; shift according to specified window center
   xout = xout - op%plon3
   if (xout > 180.) then
      xout = xout - 360.
   elseif (xout < -180.) then
      xout = xout + 360.
   endif
elseif (op%projectn(iplt) == 'P') then
   call e_ps(xeq,yeq,zeq,op%plat3,op%plon3,xout,yout)
elseif (op%projectn(iplt) == 'O') then
   call e_or(xeq,yeq,zeq,op%plat3,op%plon3,xout,yout)
else  ! Cartesian domain
   xout = xeq
   yout = yeq
endif

return
end subroutine oplot_transform

!===============================================================================

subroutine oplot_panel(panel,colorbar0,aspect)

use oplot_coms, only: op
use misc_coms,  only: io6

implicit none

character(len=1), intent(in) :: panel,colorbar0
real, intent(in) :: aspect

real :: pscale, poffset, asp

! Relative positions within panel (if colorbar is used)

op%fx1 = .12  ! plot frame left side x coord
op%fx2 = .86  ! plot frame right side x coord
op%fy1 = .20  ! plot frame bottom y coord

op%fy2 = .94  ! plot frame top y coord

op%cbx1 = .88 ! colorbar left side x coord
op%cbx2 = .91 ! colorbar right side x coord
op%cblx = .92 ! colorbar label left end x coord

op%xlaby = .14 ! x-axis label y coord
op%xtlaby = .175 ! x-axis tick label y coord

op%ylabx = .03 ! y-axis label x coord
op%ytlabx = .11 ! y-axis tick label right end x coord

op%timex = .78 ! elapsed time (sec) left end x coord

! Information block coordinates

op%slabx = .78 ! slab left end x coord

op%timsy = .11 ! elapsed time (sec) y coord
op%timdy = .09 ! elapsed time (day) y coord
op%slaby = .07 ! slab y coord
op%sminy = .05 ! slab min y coord
op%smaxy = .03 ! slab max y coord

! Expand plot area if no colorbar

if (colorbar0 /= 'c') then

   op%fx1 = .14  ! plot frame left side x coord
   op%fx2 = .97  ! plot frame right side x coord
   op%fy1 = .11  ! plot frame bottom y coord
   op%fy2 = .94  ! plot frame top y coord

   op%xlaby = .04 ! x-axis label y coord
   op%xtlaby = .08 ! x-axis tick label y coord

   op%ylabx = .03 ! y-axis label x coord
   op%ytlabx = .13 ! y-axis tick label right end x coord

   op%timex = .03  ! elapsed time (sec) left end x coord
   op%timsy = .035 ! elapsed time (sec) y coord
   op%timdy = .015 ! elapsed time (day) y coord

   op%slabx = .78  ! slab left end x coord
   op%slaby = .055 ! slab y coord
   op%sminy = .035 ! slab min y coord
   op%smaxy = .015 ! slab max y coord

endif

! Reduce plot height under some circumstances (make plot window rectangular
! instead of square)

if (aspect < 1.e-3) then

! Aspect set to 0 scales plot window to physical coordinate limits

   asp = (op%ymax - op%ymin) / (op%xmax - op%xmin)
   op%fy2 = op%fy1 + (op%fx2 - op%fx1) * asp

elseif (aspect < .999) then

! Aspect between 0 and 1 is used as specified plot window height/width ratio

   op%fy2 = op%fy1 + (op%fx2 - op%fx1) * aspect
   
endif

op%fnamey = op%fy2 + .025 ! field name y coord

! Panel borders

op%hp1 = 0.
op%hp2 = 1.
op%vp1 = 0.
op%vp2 = 1.

!pscale = .48
!poffset = .52

pscale = .50
poffset = .50



!!!!!!!!!!!!!!!!!!!!!!!!!! special - modify plot coordinates

!op%fx2 = .97  ! plot frame right side x coord
!op%fy1 = .14  ! plot frame bottom y coord

!op%xlaby = .06 ! x-axis label y coord
!op%xtlaby = .10 ! x-axis tick label y coord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



if     (panel == '1') then
   op%hp1 = pscale * op%hp1
   op%hp2 = pscale * op%hp2
   op%vp1 = pscale * op%vp1 + poffset * aspect
   op%vp2 = pscale * op%vp2 + poffset * aspect
elseif (panel == '2') then
   op%hp1 = pscale * op%hp1 + poffset
   op%hp2 = pscale * op%hp2 + poffset
   op%vp1 = pscale * op%vp1 + poffset * aspect
   op%vp2 = pscale * op%vp2 + poffset * aspect
elseif (panel == '3') then
   op%hp1 = pscale * op%hp1
   op%hp2 = pscale * op%hp2
   op%vp1 = pscale * op%vp1
   op%vp2 = pscale * op%vp2
elseif (panel == '4') then
   op%hp1 = pscale * op%hp1 + poffset
   op%hp2 = pscale * op%hp2 + poffset
   op%vp1 = pscale * op%vp1
   op%vp2 = pscale * op%vp2
endif

! Frame borders

op%h1 = op%hp1 + op%fx1 * (op%hp2 - op%hp1)
op%h2 = op%hp1 + op%fx2 * (op%hp2 - op%hp1)
op%v1 = op%vp1 + op%fy1 * (op%vp2 - op%vp1)
op%v2 = op%vp1 + op%fy2 * (op%vp2 - op%vp1)

return
end subroutine oplot_panel

!===============================================================================

subroutine oplot_set(iplt)

use misc_coms,   only: io6, mdomain
use oplot_coms,  only: op
use oname_coms,  only: nl
use mem_grid,    only: xem, yem, xew, yew, zew, zm, mma, mwa, arw0, nza, mza
use consts_coms, only: erad, pio180

implicit none

integer, intent(in) :: iplt

integer :: iu,iw,iok,im,ju1,ju2
real :: h1,h2,v1,v2,dist,radcone,s1,s2,xmid,ymid
real :: delxmin,arwmin,xpt,ypt
real :: aspect
real :: topo1,topo2

real, dimension(8) :: htpn,vtpn

! Set plot window center or cone center coordinates (whichever is relevant),
! plus cone angle

op%plon3    = nl%plotcoord1(iplt)
op%plat3    = nl%plotcoord2(iplt)
op%coneang  = nl%slabloc(iplt)
op%viewazim = nl%viewazim(iplt)

! Set limits of plot window in model domain coordinates

if (op%windowin(iplt) == 'W') then

! If windowing in, use specified window size

   op%xmin = -.5 * nl%plotwid(iplt)
   op%xmax =  .5 * nl%plotwid(iplt)
   op%ymin = -.5 * nl%plotwid(iplt)
   op%ymax =  .5 * nl%plotwid(iplt)

! Special treatment for EC (lat/lon) plot: use plot center lat as window center
! (Plot center longitude is already accounted for in transformation) 

   if (op%projectn(iplt) == 'L') then
      op%xmin = op%plon3 - .5 * nl%plotwid(iplt)
      op%xmax = op%plon3 + .5 * nl%plotwid(iplt)
      op%ymin = op%plat3 - .5 * nl%plotwid(iplt)
      op%ymax = op%plat3 + .5 * nl%plotwid(iplt)
   endif

! Special treatment for vertical plot: ymin/ymax related to model Z coordinate

   if (op%projectn(iplt) == 'C' .or. op%projectn(iplt) == 'V') then
      op%ymin = zm(1)
      op%ymax = zm(nza-1)
      
      ! special for dudhia expts
      !op%ymax = 10.e3
      ! end special
      
      if (op%stagpt == 'A' .or. op%stagpt == 'L') op%ymin = -.2 * op%ymax
   endif

else

! If not windowing in, use following defaults

   if (op%projectn(iplt) == 'L') then

      op%xmin = -180.
      op%xmax =  180.
      op%ymin =  -90.
      op%ymax =   90.

   elseif (op%projectn(iplt) == 'P') then

      op%xmin = - 5. * erad
      op%xmax =   5. * erad
      op%ymin = - 5. * erad
      op%ymax =   5. * erad

   elseif (op%projectn(iplt) == 'O') then

      op%xmin = - 1.1 * erad
      op%xmax =   1.1 * erad
      op%ymin = - 1.1 * erad
      op%ymax =   1.1 * erad

   elseif (op%projectn(iplt) == 'C') then

      op%xmin = -3.4 * erad * sin(op%coneang * pio180)
      op%xmax =  3.4 * erad * sin(op%coneang * pio180)
      op%ymin = zm(1)
      op%ymax = zm(nza-1)

      if (op%stagpt == 'A' .or. op%stagpt == 'L') op%ymin = -.2 * op%ymax  

   elseif (op%projectn(iplt) == 'X') then

      op%xmin =  1.e9
      op%xmax = -1.e9
      do im = 2,mma
         op%xmin = min(op%xmin,yem(im))
         op%xmax = max(op%xmax,yem(im))
      enddo
      op%xmin = 1.01 * op%xmin
      op%xmax = 1.01 * op%xmax
      op%ymin = zm(1)
      op%ymax = zm(nza-1)

      if (op%stagpt == 'A' .or. op%stagpt == 'L') op%ymin = -.2 * op%ymax  

   elseif (op%projectn(iplt) == 'Y') then

      op%xmin =  1.e9
      op%xmax = -1.e9
      do im = 2,mma
         op%xmin = min(op%xmin,xem(im))
         op%xmax = max(op%xmax,xem(im))
      enddo
      op%xmin = 1.01 * op%xmin
      op%xmax = 1.01 * op%xmax
      op%ymin = zm(1)
      op%ymax = zm(nza-1)

      if (op%stagpt == 'A' .or. op%stagpt == 'L') op%ymin = -.2 * op%ymax  

   elseif (op%projectn(iplt) == 'V') then

      op%xmin =  1.e9
      op%xmax = -1.e9
      do im = 2,mma
         op%xmin = min(op%xmin,xem(im),yem(im))
         op%xmax = max(op%xmax,xem(im),yem(im))
      enddo
      op%xmin = 1.01 * op%xmin
      op%xmax = 1.01 * op%xmax
      op%ymin = zm(1)
      op%ymax = zm(nza-1)

      if (op%stagpt == 'A' .or. op%stagpt == 'L') op%ymin = -.2 * op%ymax  

   elseif (op%projectn(iplt) == 'Z') then

      op%xmin =  1.e9
      op%xmax = -1.e9
      op%ymin =  1.e9
      op%ymax = -1.e9
      do im = 2,mma
         op%xmin = min(op%xmin,xem(im))
         op%xmax = max(op%xmax,xem(im))
         op%ymin = min(op%ymin,yem(im))
         op%ymax = max(op%ymax,yem(im))
      enddo
      dist = max(op%xmax - op%xmin,op%ymax - op%ymin)
      xmid = .5 * (op%xmin + op%xmax)
      ymid = .5 * (op%ymin + op%ymax)

      op%xmin = xmid - .505 * dist
      op%xmax = xmid + .505 * dist
      op%ymin = ymid - .505 * dist
      op%ymax = ymid + .505 * dist

   endif

endif   

if (op%projectn(iplt) == 'L' .or.  &
    op%projectn(iplt) == 'P' .or.  &
    op%projectn(iplt) == 'O' .or.  &
    op%projectn(iplt) == 'Z') then   ! Horizontal cross section

! Find deltax of finest grid cell IN CURRENT PLOT WINDOW

   arwmin = 1.e13
   do iw = 2,mwa
      call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),xpt,ypt)

      if (xpt < op%xmin .or. xpt > op%xmax .or.  &
          ypt < op%ymin .or. ypt > op%ymax) cycle

      arwmin = min(arwmin,arw0(iw))
   enddo
   delxmin = sqrt(arwmin)

   if (op%projectn(iplt) == 'L') then
!      op%psiz = .05 * (delxmin / (6. * erad)) * (400. / (op%xmax - op%xmin))
      op%psiz = .20 * (delxmin / (6. * erad)) * (400. / (op%xmax - op%xmin))
      op%vsprd = .3 * (delxmin / (6. * erad)) * 400.
   else
   !   op%psiz = .20 * delxmin / (op%xmax - op%xmin)
   !    op%psiz = .10 * delxmin / (op%xmax - op%xmin)
      op%psiz = .07 * delxmin / (op%xmax - op%xmin)
   !   op%psiz = .02 * delxmin / (op%xmax - op%xmin)
      op%vsprd = .3 * delxmin
   endif

else                           ! Vertical cross section

! Get horizontal plot coordinates for cells in this column

   arwmin = 1.e13
   do iw = 2,mwa

      if (op%projectn(iplt) == 'C') then
         call coneplot_w(iw,topo1,topo2,iok,htpn)
      elseif (op%projectn(iplt) == 'X') then
         call xplot_w(iplt,iw,iok,htpn)
      elseif (op%projectn(iplt) == 'Y') then
         call yplot_w(iplt,iw,iok,htpn)
      elseif (op%projectn(iplt) == 'V') then
         call xyplot_w(iplt,iw,ju1,ju2,topo1,topo2,iok,htpn)
      endif

      if (iok /= 1) cycle ! If this W pt does not intersect plot cone

! Jump out of loop if either cell side is outside plot window. 

      if (htpn(1) < op%xmin .or. htpn(1) > op%xmax .or.  &
          htpn(2) < op%xmin .or. htpn(2) > op%xmax) cycle

      arwmin = min(arwmin,arw0(iw))
   enddo
   delxmin = sqrt(arwmin)

   if (mdomain <= 1) then
      radcone = erad * sin(op%coneang * pio180)
      op%psiz = .05 * delxmin / (radcone * pio180 * (op%xmax - op%xmin))
   else
      op%psiz = .06 * delxmin / (op%xmax - op%xmin)
   endif
   op%vsprd = 50.  ! for vert xsectn

endif

! Set limits of panel and frame window in plotter coordinates

if (op%projectn(iplt) == 'L' .and.  &
    op%windowin(iplt) /= 'W') then
    
   aspect = 0.
else
   aspect = 1.
endif
    
call oplot_panel(op%panel(iplt),op%colorbar(iplt),aspect)

! Scale plot window to selected model domain

call o_set(op%h1,op%h2,op%v1,op%v2,op%xmin,op%xmax,op%ymin,op%ymax,1)

return
end subroutine oplot_set

