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
subroutine contslab_hmp(iplt)

use oplot_coms, only: op
use mem_grid,   only: mza, mwa, zm, zt, lpw, xem, yem, zem, xew, yew, zew
use mem_ijtabs, only: itab_w
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt

integer :: ncpn,kt,k,ko,iw,itpn,ici,icpn,im1,im2,im3,notavail
integer :: iflag180
integer :: ipwx1,ipwx2,ipwy1,ipwy2

real :: hpt,vpt
real, dimension(4) :: hcpn,vcpn,fldvals

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

if (op%dimens == '3') then
   call plot_underground(iplt,kt)
endif

! Loop over W points for contouring M points

do iw = 2,mwa

! Get plot coordinates of current W point.  

   call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),hpt,vpt)

! Get plot coordinates of 3 M points.  

   im1 = itab_w(iw)%im1
   im2 = itab_w(iw)%im2
   im3 = itab_w(iw)%im3
   
   call oplot_transform(iplt,xem(im1),yem(im1),zem(im1),hcpn(1),vcpn(1))
   call oplot_transform(iplt,xem(im2),yem(im2),zem(im2),hcpn(2),vcpn(2))
   call oplot_transform(iplt,xem(im3),yem(im3),zem(im3),hcpn(3),vcpn(3))

! Set iflag180, and avoid wrap-around if lat/lon plot

   iflag180 = 0

   if (op%projectn(iplt) == 'L') then
      do icpn = 1,3
         call ll_unwrap(hpt,hcpn(icpn))
         if (hcpn(icpn) < -180.) iflag180 =  1
         if (hcpn(icpn) >  180.) iflag180 = -1
      enddo
   endif

! Initialize plot window flags to zero

   ipwx1 = 0
   ipwx2 = 0
   ipwy1 = 0
   ipwy2 = 0

! Loop over the 3 M points surrounding current W point

   do icpn = 1,3

! Skip this W point if any M point is far outside plot window

      if (abs(hcpn(icpn)) > 1.e11) go to 9
      
! Set plot window flag to 1 if any M point is on window side of 
! respective boundary

      if (hcpn(icpn) >= op%xmin) ipwx1 = 1
      if (hcpn(icpn) <= op%xmax) ipwx2 = 1
      if (vcpn(icpn) >= op%ymin) ipwy1 = 1
      if (vcpn(icpn) <= op%ymax) ipwy2 = 1

   enddo
   
! If any window flag is zero, all M points for this W point are outside
! the same window boundary, so skip this W point

   if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) go to 9

! Fill field value M points around current W point

   call oplot_lib(iplt,k,im1,'VALUE',op%fldname(iplt),fldvals(1),notavail,ko)
   if (notavail > 0) cycle 
   call oplot_lib(iplt,k,im2,'VALUE',op%fldname(iplt),fldvals(2),notavail,ko)
   if (notavail > 0) cycle 
   call oplot_lib(iplt,k,im3,'VALUE',op%fldname(iplt),fldvals(3),notavail,ko)
   if (notavail > 0) cycle 

! Contour plot cell of 2-D or 3-D field

   call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn,vcpn,fldvals)

! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
! at other end

   if (iflag180 /= 0) then

      do icpn = 1,3
         hcpn(icpn) = hcpn(icpn) + 360. * iflag180
      enddo
      
      call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn,vcpn,fldvals)
            
   endif

9  continue
            
enddo

return
end subroutine contslab_hmp

!===============================================================================

subroutine contslab_hnu(iplt)

use oplot_coms, only: op
use mem_grid,   only: mma, mwa, mza, zm, zt, lpu, lpw, xem, yem, zem,  &
                      xeu, yeu, zeu, xew, yew, zew
use mem_ijtabs, only: itab_m, itab_w, itab_u
use misc_coms,  only: io6, iparallel
use mem_para,   only: myrank

implicit none

integer, intent(in) :: iplt

integer :: ncpn,kt,k,ko,im,icpn,ici,iu,notavail,iu1,iu2,iu3,iw
integer :: iflag180
integer :: ipwx1,ipwx2,ipwy1,ipwy2

real :: hpt,vpt
real, dimension(8) :: hcpn,vcpn,fldvals

! Set integer ifill flag and color table value for underground cells

! Find KT or KM level to plot

kt = 2
do while (kt < mza .and. zm(kt) < op%slabloc(iplt))
   kt = kt + 1
enddo
k = kt

! If plotting N point, set k for correct N level

if (op%stagpt == 'N' .and. zt(kt) > op%slabloc(iplt)) k = kt - 1

! If field is 3d, first plot underground points with underground color

if (op%dimens == '3') then
   call plot_underground_m(iplt,kt)
endif

!------------------------------------------------------
! FIRST LOOP is over M points for contouring U points
!------------------------------------------------------

do im = 2,mma

! Get plot coordinates of current M point.

   call oplot_transform(iplt,xem(im),yem(im),zem(im),hpt,vpt)

   ncpn = itab_m(im)%ntpn

! Initialize iflag180 and plot window flags to zero

   iflag180 = 0

   ipwx1 = 0
   ipwx2 = 0
   ipwy1 = 0
   ipwy2 = 0

! Loop over all U points that surround current M point

   do icpn = 1,ncpn

! Current U point index   

      iu = itab_m(im)%iu(icpn)
      
! TEMPORARY FIX TO AVOID ACCESSING UNDEFINED VALUES IN PARALLEL
      if (iparallel == 1) then
         if (itab_u(iu)%irank /= myrank) goto 8
      endif

! Skip current U cell if index < 2

      if (iu < 2) go to 8

! Get plot coordinates of current U point

      call oplot_transform(iplt,xeu(iu),yeu(iu),zeu(iu),hcpn(icpn),vcpn(icpn))

! Skip this M point if current U point is far outside plot window 
! (which means that orthographic projection returned large value that
! indicates that point is on other side of Earth)

      if (abs(hcpn(icpn)) > 1.e11) go to 8
      
! Avoid wrap-around and set iflag180

      if (op%projectn(iplt)== 'L') then
         call ll_unwrap(hpt,hcpn(icpn))
         if (hcpn(icpn) < -180.001) iflag180 =  1
         if (hcpn(icpn) >  180.001) iflag180 = -1
      endif

! Set plot window flag to 1 if any M point is on window side of 
! respective boundary

      if (hcpn(icpn) >= op%xmin) ipwx1 = 1 
      if (hcpn(icpn) <= op%xmax) ipwx2 = 1 
      if (vcpn(icpn) >= op%ymin) ipwy1 = 1 
      if (vcpn(icpn) <= op%ymax) ipwy2 = 1 

   enddo
   
! If any window flag is zero, all M points for this U point are outside
! the same window boundary, so skip this U point

   if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) go to 8

! Loop over all U points that surround current M point and fill field values

   do icpn = 1,ncpn
      iu = itab_m(im)%iu(icpn)
      
      call oplot_lib(iplt,k,iu,'VALUE',op%fldname(iplt),fldvals(icpn),notavail,ko)
      if (notavail > 0) cycle 
   enddo

! Contour plot cell of 2-D or 3-D field

   call contpolyg(op%icolortab(iplt),op%ifill,ncpn,hcpn,vcpn,fldvals)
            
! If lat/lon plot and this polygon crosses +/- 180 degrees longitude, plot
! again at other end of plot

   if (iflag180 /= 0) then

      do icpn = 1,ncpn
         hcpn(icpn) = hcpn(icpn) + 360. * iflag180
      enddo
      
      call contpolyg(op%icolortab(iplt),op%ifill,ncpn,hcpn,vcpn,fldvals)
            
   endif

8  continue

enddo

!------------------------------------------------------
! SECOND LOOP is over W points for contouring U points
!------------------------------------------------------

do iw = 2,mwa

! Get plot coordinates of current W point.  

   call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),hpt,vpt)

! Get plot coordinates of 3 U points.  

   iu1 = itab_w(iw)%iu1
   iu2 = itab_w(iw)%iu2
   iu3 = itab_w(iw)%iu3
   
   call oplot_transform(iplt,xeu(iu1),yeu(iu1),zeu(iu1),hcpn(1),vcpn(1))
   call oplot_transform(iplt,xeu(iu2),yeu(iu2),zeu(iu2),hcpn(2),vcpn(2))
   call oplot_transform(iplt,xeu(iu3),yeu(iu3),zeu(iu3),hcpn(3),vcpn(3))

! Set iflag180, and avoid wrap-around if lat/lon plot

   iflag180 = 0

   if (op%projectn(iplt) == 'L') then
      do icpn = 1,3
         call ll_unwrap(hpt,hcpn(icpn))
         if (hcpn(icpn) < -180.) iflag180 =  1
         if (hcpn(icpn) >  180.) iflag180 = -1
      enddo
   endif

! Initialize plot window flags to zero

   ipwx1 = 0
   ipwx2 = 0
   ipwy1 = 0
   ipwy2 = 0

! Loop over the 3 M points surrounding current W point

   do icpn = 1,3

! Skip this W point if any M point is far outside plot window

      if (abs(hcpn(icpn)) > 1.e11) go to 9
      
! Set plot window flag to 1 if any M point is on window side of 
! respective boundary

      if (hcpn(icpn) >= op%xmin) ipwx1 = 1
      if (hcpn(icpn) <= op%xmax) ipwx2 = 1
      if (vcpn(icpn) >= op%ymin) ipwy1 = 1
      if (vcpn(icpn) <= op%ymax) ipwy2 = 1

   enddo
   
! If any window flag is zero, all M points for this W point are outside
! the same window boundary, so skip this W point

   if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) go to 9

! Fill field value M points around current W point

   call oplot_lib(iplt,k,iu1,'VALUE',op%fldname(iplt),fldvals(1),notavail,ko)
   if (notavail > 0) cycle 
   call oplot_lib(iplt,k,iu2,'VALUE',op%fldname(iplt),fldvals(2),notavail,ko)
   if (notavail > 0) cycle 
   call oplot_lib(iplt,k,iu3,'VALUE',op%fldname(iplt),fldvals(3),notavail,ko)
   if (notavail > 0) cycle 

! Contour plot cell of 2-D or 3-D field

   call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn,vcpn,fldvals)

! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
! at other end

   if (iflag180 /= 0) then

      do icpn = 1,3
         hcpn(icpn) = hcpn(icpn) + 360. * iflag180
      enddo
      
      call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn,vcpn,fldvals)
            
   endif

9  continue
            
enddo

return
end subroutine contslab_hnu

!===============================================================================

subroutine contslab_htw(iplt)

use oplot_coms, only: op
use mem_grid,   only: mma, mza, zm, zt, lpw, xem, yem, zem, xew, yew, zew
use mem_ijtabs, only: itab_m
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt

integer :: ncpn,kt,k,ko,im,icpn,ici,iwcpn,iw,notavail
integer :: iflag180
integer :: ipwx1,ipwx2,ipwy1,ipwy2

real :: hpt,vpt
real, dimension(8) :: hcpn,vcpn,fldvals

! Set integer ifill flag and color table value for underground cells

! Find KT or KM level to plot

kt = 2
do while (kt < mza .and. zm(kt) < op%slabloc(iplt))
   kt = kt + 1
enddo
k = kt

! If plotting W point, set k for correct W level

if (op%stagpt == 'W' .and. zt(kt) > op%slabloc(iplt)) k = kt - 1

! If field is 3d, first plot underground points with underground color

if (op%dimens == '3') then
   call plot_underground_m(iplt,kt)
endif

! Loop over M points for contouring W points

do im = 2,mma

! Get plot coordinates of current M point.

   call oplot_transform(iplt,xem(im),yem(im),zem(im),hpt,vpt)

   ncpn = itab_m(im)%ntpn

! Initialize iflag180 and plot window flags to zero

   iflag180 = 0

   ipwx1 = 0
   ipwx2 = 0
   ipwy1 = 0
   ipwy2 = 0

! Loop over all W points that surround current M point

   do icpn = 1,ncpn

! Current W point index   

      iw = itab_m(im)%iw(icpn)

! Skip current W cell if index < 2

      if (iw < 2) go to 9

! Get plot coordinates of current W point

      call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),hcpn(icpn),vcpn(icpn))

! Skip this M point if current W point is far outside plot window 
! (which means that orthographic projection returned large value that
! indicates that point is on other side of Earth)

      if (abs(hcpn(icpn)) > 1.e11) go to 9
      
! Avoid wrap-around and set iflag180

      if (op%projectn(iplt)== 'L') then
         call ll_unwrap(hpt,hcpn(icpn))
         if (hcpn(icpn) < -180.001) iflag180 =  1
         if (hcpn(icpn) >  180.001) iflag180 = -1
      endif

! Set plot window flag to 1 if any M point is on window side of 
! respective boundary

      if (hcpn(icpn) >= op%xmin) ipwx1 = 1 
      if (hcpn(icpn) <= op%xmax) ipwx2 = 1 
      if (vcpn(icpn) >= op%ymin) ipwy1 = 1 
      if (vcpn(icpn) <= op%ymax) ipwy2 = 1 

   enddo
   
! If any window flag is zero, all M points for this W point are outside
! the same window boundary, so skip this W point

   if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) go to 9

! Loop over all W points that surround current M point and fill field values

   do icpn = 1,ncpn
      iw = itab_m(im)%iw(icpn)
      call oplot_lib(iplt,k,iw,'VALUE',op%fldname(iplt),fldvals(icpn),notavail,ko)
      if (notavail > 0) cycle 
   enddo

! Contour plot cell of 2-D or 3-D field

   call contpolyg(op%icolortab(iplt),op%ifill,ncpn,hcpn,vcpn,fldvals)
            
! If lat/lon plot and this polygon crosses +/- 180 degrees longitude, plot
! again at other end of plot

   if (iflag180 /= 0) then

      do icpn = 1,ncpn
         hcpn(icpn) = hcpn(icpn) + 360. * iflag180
      enddo
      
      call contpolyg(op%icolortab(iplt),op%ifill,ncpn,hcpn,vcpn,fldvals)
            
   endif

9  continue

enddo

return
end subroutine contslab_htw

!===============================================================================

subroutine contslab_vt(iplt)

use oplot_coms,  only: op
use mem_grid,    only: mua, mza, lpw, zt
use mem_ijtabs,  only: itab_u
use misc_coms,   only: io6
use consts_coms, only: erad, pio180

implicit none

integer, intent(in) :: iplt

integer :: k,ko,iu,iw1,iw2,ju1,ju2,iok,notavail
real :: hpt,hpt1,hpt2
real, dimension(4) :: hcpn,hcpn1,hcpn2,vcpn,fldvals
real :: topo1,topo2,radcone

! First plot underground T cells with underground color

call plot_underground(iplt,0)

do iu = 2,mua  ! Loop is over U for contouring T points

   iw1 = itab_u(iu)%iw1
   iw2 = itab_u(iu)%iw2
   
! Jump to end of loop if either iw1 or iw2 is less than 1

   if (iw1 < 2 .or. iw2 < 2) cycle

! Get horizontal plot coordinates for IW1 point

   if (op%projectn(iplt) == 'C') then
      call coneplot_w(iw1,topo1,topo2,iok,hcpn1)
   elseif (op%projectn(iplt)(1:1) == 'V') then
      call xyplot_w(iplt,iw1,ju1,ju2,topo1,topo2,iok,hcpn1)
   endif

! Skip current IU point if this IW1 point does not intersect plot cone

   if (iok /= 1) cycle

   hpt1 = .5 * (hcpn1(1) + hcpn1(2))

! Get horizontal plot coordinates for IW2 point

   if (op%projectn(iplt) == 'C') then
      call coneplot_w(iw2,topo1,topo2,iok,hcpn2)
   elseif (op%projectn(iplt)(1:1) == 'V') then
      call xyplot_w(iplt,iw2,ju1,ju2,topo1,topo2,iok,hcpn2)
   endif

! Skip current IU point if this IW2 point does not intersect plot cone

   if (iok /= 1) cycle

   hpt2 = .5 * (hcpn2(1) + hcpn2(2))

! For now, skip point if +/- 180 degree point is crossed.  Later, truncate cells.

   if (op%projectn(iplt) == 'C') then
      radcone = erad * sin(op%coneang * pio180)

      if (abs(hpt2 - hpt1) > 3. * radcone) cycle
   endif

! Skip this IU point if either cell side is outside plot window. 

   if (hpt1 < op%xmin .or. hpt1 > op%xmax .or.  &
       hpt2 < op%xmin .or. hpt2 > op%xmax) cycle
   
   hcpn(1) = hpt1
   hcpn(2) = hpt2
   hcpn(3) = hcpn(2)
   hcpn(4) = hcpn(1)
   
   do k = 2,mza-2   ! Loop is over M levels

! Skip plot if any T cell around current M point is below ground 
! or is cell #1 (for now; later maybe draw contours 
! across partial cells).

      if (iw1 < 2 .or. iw2 < 2)           cycle
      if (k < lpw(iw1) .or. k < lpw(iw2)) cycle
   	
! Jump to end loop if either upper or lower cell center is outside plot window. 

   if (zt(k) < op%ymin .or. zt(k+1) > op%ymax) cycle
   
! Get T-cell vertical coordinates

      vcpn(1) = zt(k)
      vcpn(2) = vcpn(1)
      vcpn(3) = zt(k+1)
      vcpn(4) = vcpn(3)

! special for dudhia expts
!if (vcpn(4) > op%ymax) cycle
! end special

! Fill field values of 4 T points around current M point

      call oplot_lib(iplt,k,iw1,'VALUE',op%fldname(iplt),fldvals(1),notavail,ko)
      if (notavail > 0) cycle 
      call oplot_lib(iplt,k,iw2,'VALUE',op%fldname(iplt),fldvals(2),notavail,ko)
      if (notavail > 0) cycle 
      call oplot_lib(iplt,k+1,iw2,'VALUE',op%fldname(iplt),fldvals(3),notavail,ko)
      if (notavail > 0) cycle 
      call oplot_lib(iplt,k+1,iw1,'VALUE',op%fldname(iplt),fldvals(4),notavail,ko)
      if (notavail > 0) cycle 

! Contour plot cell around current M point

      call contpolyg(op%icolortab(iplt),op%ifill,4,hcpn,vcpn,fldvals)
      
   enddo

enddo
   
return
end subroutine contslab_vt

!===============================================================================

subroutine contslab_vu(iplt)

use oplot_coms, only: op
use mem_grid,   only: mwa, mza, lpw, zt
use misc_coms,  only: io6
use consts_coms, only: erad, pio180

implicit none

integer, intent(in) :: iplt

integer :: k,ko,iu,iw,ju1,ju2,iok,notavail
real :: hpt
real, dimension(4) :: hcpn,vcpn,fldvals
real :: topo1,topo2

! First plot underground T cells with underground color

call plot_underground(iplt,0)

do iw = 2,mwa  ! Loop is over W for contouring U points

! Get horizontal plot coordinates for IW point

   if (op%projectn(iplt) == 'C') then
!!needs ju's      call coneplot_w(iw,topo1,topo,iok,hcpn)
   elseif (op%projectn(iplt)(1:1) == 'V') then
      call xyplot_w(iplt,iw,ju1,ju2,topo1,topo2,iok,hcpn)
   endif

! Skip current IW point if it does not intersect plot cone

   if (iok /= 1) cycle

! Avoid wrap_around

   if (op%projectn(iplt) == 'C') then
      call ll_unwrap(hcpn(1),hcpn(2))
   endif

   hpt = .5 * (hcpn(1) + hcpn(2))

! Skip current IU point if either cell side is outside plot window. 

   if (hcpn(1) < op%xmin .or. hcpn(1) > op%xmax .or.  &
       hcpn(2) < op%xmin .or. hcpn(2) > op%xmax) cycle
   
   hcpn(3) = hcpn(2)
   hcpn(4) = hcpn(1)
   
   do k = lpw(iw),mza-2   ! Loop is over W levels

! Skip this K point if either upper or lower cell center is outside plot window. 

   if (zt(k) < op%ymin .or. zt(k+1) > op%ymax) cycle
   
! Get T-cell vertical coordinates

      vcpn(1) = zt(k)
      vcpn(2) = vcpn(1)
      vcpn(3) = zt(k+1)
      vcpn(4) = vcpn(3)

! special for dudhia expts
!if (vcpn(4) > op%ymax) cycle
! end special

! Fill field values of 4 T points around current M point

      call oplot_lib(iplt,k,ju1,'VALUE',op%fldname(iplt),fldvals(1),notavail,ko)
      if (notavail > 0) cycle 
      call oplot_lib(iplt,k,ju2,'VALUE',op%fldname(iplt),fldvals(2),notavail,ko)
      if (notavail > 0) cycle 
      call oplot_lib(iplt,k+1,ju2,'VALUE',op%fldname(iplt),fldvals(3),notavail,ko)
      if (notavail > 0) cycle 
      call oplot_lib(iplt,k+1,ju1,'VALUE',op%fldname(iplt),fldvals(4),notavail,ko)
      if (notavail > 0) cycle 

! Contour plot cell around current M point

      call contpolyg(op%icolortab(iplt),op%ifill,4,hcpn,vcpn,fldvals)
      
   enddo

enddo
   
return
end subroutine contslab_vu

!===============================================================================

subroutine contslab_vw(iplt)

use oplot_coms, only: op
use mem_grid,   only: mua, mza, lpw, zm
use mem_ijtabs, only: itab_u
use misc_coms,  only: io6
use consts_coms, only: erad, pio180

implicit none

integer, intent(in) :: iplt

integer :: k,ko,iu,iw1,iw2,ju1,ju2,iok,notavail
real :: hpt,hpt1,hpt2
real, dimension(4) :: hcpn,hcpn1,hcpn2,vcpn,fldvals
real :: topo1,topo2,radcone

! First plot underground T cells with underground color

call plot_underground(iplt,0)

do iu = 2,mua  ! Loop is over U for contouring W points

   iw1 = itab_u(iu)%iw1
   iw2 = itab_u(iu)%iw2

! Jump to end of loop if either iw1 or iw2 is less than 1

   if (iw1 < 2 .or. iw2 < 2) cycle

! Get horizontal plot coordinates for IW1 point

   if (op%projectn(iplt) == 'C') then
      call coneplot_w(iw1,topo1,topo2,iok,hcpn1)
   elseif (op%projectn(iplt)(1:1) == 'V') then
      call xyplot_w(iplt,iw1,ju1,ju2,topo1,topo2,iok,hcpn1)
   endif

! Skip current IU point if IW1 point does not intersect plot cone

   if (iok /= 1) cycle

   hpt1 = .5 * (hcpn1(1) + hcpn1(2))

! Get horizontal plot coordinates for IW2 point

   if (op%projectn(iplt) == 'C') then
      call coneplot_w(iw2,topo1,topo2,iok,hcpn2)
   elseif (op%projectn(iplt)(1:1) == 'V') then
      call xyplot_w(iplt,iw2,ju1,ju2,topo1,topo2,iok,hcpn2)
   endif

! Skip current IU point if IW2 point does not intersect plot cone

   if (iok /= 1) cycle

   hpt2 = .5 * (hcpn2(1) + hcpn2(2))

! For now, skip point if +/- 180 degree point is crossed.  Later, truncate cells.

   if (op%projectn(iplt) == 'C') then
      radcone = erad * sin(op%coneang * pio180)

      if (abs(hpt2 - hpt1) > 3. * radcone) cycle
   endif

! Skip this IU point if either cell side is outside plot window. 

   if (hpt1 < op%xmin .or. hpt1 > op%xmax .or.  &
       hpt2 < op%xmin .or. hpt2 > op%xmax) cycle

   hcpn(1) = hpt1
   hcpn(2) = hpt2
   hcpn(3) = hcpn(2)
   hcpn(4) = hcpn(1)
   
   do k = 2,mza-1   ! Loop is over T levels

! Skip plot if either T cell around current U point is below ground 
! or is cell #1 (for now; later maybe draw contours 
! across partial cells).

      if (iw1 < 2 .or. iw2 < 2)           cycle
      if (k < lpw(iw1) .or. k < lpw(iw2)) cycle

! Get W-cell vertical coordinates

      vcpn(1) = zm(k-1)
      vcpn(2) = vcpn(1)
      vcpn(3) = zm(k)
      vcpn(4) = vcpn(3)

! special for dudhia expts
!if (vcpn(4) > op%ymax) cycle
! end special

! Fill field values of 4 W points around current U point

      call oplot_lib(iplt,k-1,iw1,'VALUE',op%fldname(iplt),fldvals(1),notavail,ko)
      if (notavail > 0) cycle 
      call oplot_lib(iplt,k-1,iw2,'VALUE',op%fldname(iplt),fldvals(2),notavail,ko)
      if (notavail > 0) cycle 
      call oplot_lib(iplt,k,iw2,'VALUE',op%fldname(iplt),fldvals(3),notavail,ko)
      if (notavail > 0) cycle 
      call oplot_lib(iplt,k,iw1,'VALUE',op%fldname(iplt),fldvals(4),notavail,ko)
      if (notavail > 0) cycle 

! Contour plot cell around current M point

      call contpolyg(op%icolortab(iplt),op%ifill,4,hcpn,vcpn,fldvals)

   enddo

enddo

return
end subroutine contslab_vw

