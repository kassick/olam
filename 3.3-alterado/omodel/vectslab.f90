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
subroutine vectslab_h(iplt)

use oplot_coms, only: op
use mem_grid,   only: mza, mwa, lpw, zm, xew, yew, zew
use mem_ijtabs, only: itab_w
use mem_basic,  only: press, uc

use misc_coms,  only: io6, mdomain

use consts_coms

implicit none

integer, intent(in) :: iplt

integer :: kt,k,iw,ng,iu1,iu2,iu3,notavail

real :: uc1,uc2,uc3,wtbot  &
   ,uavg,vavg,pointx  &
   ,pointy,tailx,taily,stemangle,headangle1,headangle2   &
   ,headlen,head1x  &
   ,head1y,head2x,head2y,ptx,pty,tlx,tly,h1x,h1y,h2x,h2y  &
   ,fac,headlent  &
   ,vxu1,vyu1,vzu1,vxu2,vyu2,vzu2,vxu3,vyu3,vzu3,vx,vy,vz
real :: stemx,stemy,stemz,stemlen,snx,sny,snz,rnx,rny,rnz  &
   ,tailxe,tailye,tailze,head1xe,head1ye,head1ze,head2xe,head2ye,head2ze
real :: plev

! In case vector plot is to be made on pressure surface, define pressure value
! of that surface

plev = 0.

if (trim(op%fldname(iplt)) == 'Z1000MB') plev = 100000.
if (trim(op%fldname(iplt)) == 'Z850MB' ) plev = 85000.
if (trim(op%fldname(iplt)) == 'Z700MB' ) plev = 70000.
if (trim(op%fldname(iplt)) == 'Z500MB' ) plev = 50000.
if (trim(op%fldname(iplt)) == 'Z300MB' ) plev = 30000.
if (trim(op%fldname(iplt)) == 'Z200MB' ) plev = 20000.
if (trim(op%fldname(iplt)) == 'Z100MB' ) plev = 10000.
if (trim(op%fldname(iplt)) == 'Z50MB'  ) plev = 5000.
if (trim(op%fldname(iplt)) == 'Z30MB'  ) plev = 3000.
if (trim(op%fldname(iplt)) == 'Z10MB'  ) plev = 1000.
   
! Loop over all W points   

do iw = 2,mwa

! Skip this point if we want to plot vectors on a coarser mesh level

   if (itab_w(iw)%mrlw_orig > op%vec_maxmrl) cycle

   iu1 = itab_w(iw)%iu1      	
   iu2 = itab_w(iw)%iu2      	
   iu3 = itab_w(iw)%iu3      	

   vxu1 = itab_w(iw)%vxu1; vyu1 = itab_w(iw)%vyu1; vzu1 = itab_w(iw)%vzu1
   vxu2 = itab_w(iw)%vxu2; vyu2 = itab_w(iw)%vyu2; vzu2 = itab_w(iw)%vzu2
   vxu3 = itab_w(iw)%vxu3; vyu3 = itab_w(iw)%vyu3; vzu3 = itab_w(iw)%vzu3

   if (plev < .1) then

! If plev = 0, plot is on constant-height surface.  Find KT level of that surface.

      kt = 2
      do while (kt < mza .and. zm(kt) < op%slabloc(iplt))
         kt = kt + 1
      enddo

! If kt level is below ground, cycle

      if (kt < lpw(iw)) cycle

! Cell is above ground.  Get 3 wind components at KT level. 

      uc1 = uc(kt,iu1)
      uc2 = uc(kt,iu2)
      uc3 = uc(kt,iu3)

   else

! If plev > 0, plot is on constant-pressure surface.  Find KT level just below that surface.

      kt = 2
      do while (press(kt+1,iw) > plev .and. kt < mza-1)
         kt = kt + 1
      enddo  

! If kt level is below ground, cycle

      if (kt < lpw(iw)) cycle

! Interpolate 3 wind components vertically to pressure surface

       wtbot = (plev - press(kt+1,iw)) / (press(kt,iw) - press(kt+1,iw))    

       uc1 = wtbot * uc(kt,iu1) + (1. - wtbot) * uc(kt+1,iu1)
       uc2 = wtbot * uc(kt,iu2) + (1. - wtbot) * uc(kt+1,iu2)
       uc3 = wtbot * uc(kt,iu3) + (1. - wtbot) * uc(kt+1,iu3)

   endif

! 3D vector displacement (in time interval op%dtvec)

   stemx = (vxu1 * uc1 + vxu2 * uc2 + vxu3 * uc3) * op%dtvec
   stemy = (vyu1 * uc1 + vyu2 * uc2 + vyu3 * uc3) * op%dtvec
   stemz = (vzu1 * uc1 + vzu2 * uc2 + vzu3 * uc3) * op%dtvec

! Vector length and unit components

   stemlen = max(1.e-6,sqrt(stemx**2 + stemy**2 + stemz**2))      
   snx = stemx / stemlen
   sny = stemy / stemlen
   snz = stemz / stemlen

! "Right" unit components

   if (mdomain <= 1) then  ! Spherical geometry case
      rnx = (sny * zew(iw) - snz * yew(iw)) / erad
      rny = (snz * xew(iw) - snx * zew(iw)) / erad
      rnz = (snx * yew(iw) - sny * xew(iw)) / erad
   else                    ! Cartesian case
      rnx = sny
      rny = - snx
      rnz = 0.
   endif

! Earth coordinates of tail

   tailxe = xew(iw) - stemx
   tailye = yew(iw) - stemy
   tailze = zew(iw) - stemz

! Earth coordinates of left and right head tips

   headlen = op%headspeed * op%dtvec

   head1xe = xew(iw) + rnx * .42 * headlen - snx * .91 * headlen
   head1ye = yew(iw) + rny * .42 * headlen - sny * .91 * headlen
   head1ze = zew(iw) + rnz * .42 * headlen - snz * .91 * headlen

   head2xe = xew(iw) - rnx * .42 * headlen - snx * .91 * headlen
   head2ye = yew(iw) - rny * .42 * headlen - sny * .91 * headlen
   head2ze = zew(iw) - rnz * .42 * headlen - snz * .91 * headlen

! Transform all coordinates

   call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),pointx,pointy)
   call oplot_transform(iplt,tailxe,tailye,tailze,tailx,taily)
   call oplot_transform(iplt,head1xe,head1ye,head1ze,head1x,head1y)
   call oplot_transform(iplt,head2xe,head2ye,head2ze,head2x,head2y)

! Avoid wrap-around

   if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,tailx)
   if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,head1x)
   if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,head2x)

! Jump out of loop if vector head or tail is outside plot window. 

   if (pointx < op%xmin .or. pointx > op%xmax .or.  &
       pointy < op%ymin .or. pointy > op%ymax .or.  &
       tailx  < op%xmin .or. tailx  > op%xmax .or.  &
       taily  < op%ymin .or. taily  > op%ymax .or.  &
       head1x < op%xmin .or. head1x > op%xmax .or.  &
       head1y < op%ymin .or. head1y > op%ymax .or.  &
       head2x < op%xmin .or. head2x > op%xmax .or.  &
       head2y < op%ymin .or. head2y > op%ymax) cycle

! Compute re-scaling factor from transformation- STILL NEED THIS??????????

   fac = 1.0

! Draw vector

   call o_frstpt(pointx+fac*(tailx-pointx),pointy+fac*(taily-pointy))
   call o_vector(pointx,pointy)
   call o_frstpt(pointx+fac*(head1x-pointx),pointy+fac*(head1y-pointy))
   call o_vector(pointx,pointy)
   call o_vector(pointx+fac*(head2x-pointx),pointy+fac*(head2y-pointy))

enddo

return
end subroutine vectslab_h

!===============================================================================

subroutine vecuslab_h(iplt)

use oplot_coms,  only: op
use mem_grid,    only: mza, mua, zm, lpw, unx, uny, unz, xeu, yeu, zeu,  &
                       utx, uty, utz
use mem_ijtabs,  only: itab_u
use misc_coms,   only: io6,mdomain
use consts_coms, only: erad

implicit none

integer, intent(in) :: iplt

integer :: kt,k,ko,iu,iw1,iw2,notavail,im1,im2

real :: fldval,fldval2,uavg,vavg,pointx  &
   ,pointy,tailx,taily,stemangle,headangle1,headangle2   &
   ,headlen,head1x  &
   ,head1y,head2x,head2y,ptx,pty,tlx,tly,h1x,h1y,h2x,h2y  &
   ,fac,headlent  &
   ,vxu1,vyu1,vzu1,vxu2,vyu2,vzu2,vxu3,vyu3,vzu3,vx,vy,vz  &
   ,tailxe,tailye,tailze,stemlen
real :: stemx,stemy,stemz,snx,sny,snz,rnx,rny,rnz
real :: head1xe,head1ye,head1ze,head2xe,head2ye,head2ze

! Find KT level to plot

kt = 2
do while (kt < mza .and. zm(kt) < op%slabloc(iplt))
   kt = kt + 1
enddo
k = kt

do iu = 2,mua

   im1 = itab_u(iu)%im1
   im2 = itab_u(iu)%im2

   iw1 = itab_u(iu)%iw1
   iw2 = itab_u(iu)%iw2

! Check if both neighboring W cells are above ground

   if (kt >= lpw(iw1) .and. kt >= lpw(iw2)) then

! Cell is above ground 

      call oplot_lib(iplt,kt,iu,'VALUV','UC',fldval,notavail,ko)

! 3D vector displacement (in time interval op%dtvec)...
 
      if (op%vectbarb(iplt) == 'U') then

! Plot normal component (UC) only at U point

         stemx = unx(iu) * fldval * op%dtvec
         stemy = uny(iu) * fldval * op%dtvec
         stemz = unz(iu) * fldval * op%dtvec

      else

! Plot total horizontal vector at U point

         call oplot_lib(iplt,kt,iu,'VALUV','VC',fldval2,notavail,ko)

         stemx = (unx(iu) * fldval + utx(iu) * fldval2) * op%dtvec
         stemy = (uny(iu) * fldval + uty(iu) * fldval2) * op%dtvec
         stemz = (unz(iu) * fldval + utz(iu) * fldval2) * op%dtvec

      endif

! Vector length and unit components

      stemlen = max(1.e-6,sqrt(stemx**2 + stemy**2 + stemz**2))      
      snx = stemx / stemlen
      sny = stemy / stemlen
      snz = stemz / stemlen

! "Right" unit components

      if (mdomain <= 1) then  ! Spherical geometry case
         rnx = (sny * zeu(iu) - snz * yeu(iu)) / erad
         rny = (snz * xeu(iu) - snx * zeu(iu)) / erad
         rnz = (snx * yeu(iu) - sny * xeu(iu)) / erad
      else                    ! Cartesian case
         rnx = sny
         rny = - snx
         rnz = 0.
      endif

! Earth coordinates of tail

      tailxe = xeu(iu) - stemx
      tailye = yeu(iu) - stemy
      tailze = zeu(iu) - stemz

! Earth coordinates of left and right head tips

      headlen = op%headspeed * op%dtvec

      head1xe = xeu(iu) + rnx * .42 * headlen - snx * .91 * headlen
      head1ye = yeu(iu) + rny * .42 * headlen - sny * .91 * headlen
      head1ze = zeu(iu) + rnz * .42 * headlen - snz * .91 * headlen

      head2xe = xeu(iu) - rnx * .42 * headlen - snx * .91 * headlen
      head2ye = yeu(iu) - rny * .42 * headlen - sny * .91 * headlen
      head2ze = zeu(iu) - rnz * .42 * headlen - snz * .91 * headlen

! Transform all coordinates

      call oplot_transform(iplt,xeu(iu),yeu(iu),zeu(iu),pointx,pointy)
      call oplot_transform(iplt,tailxe,tailye,tailze,tailx,taily)
      call oplot_transform(iplt,head1xe,head1ye,head1ze,head1x,head1y)
      call oplot_transform(iplt,head2xe,head2ye,head2ze,head2x,head2y)

! Avoid wrap-around

      if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,tailx)
      if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,head1x)
      if (op%projectn(iplt) == 'L') call ll_unwrap(pointx,head2x)

! Jump out of loop if vector head or tail is outside plot window. 

      if (pointx < op%xmin .or. pointx > op%xmax .or.  &
          pointy < op%ymin .or. pointy > op%ymax .or.  &
          tailx  < op%xmin .or. tailx  > op%xmax .or.  &
          taily  < op%ymin .or. taily  > op%ymax .or.  &
          head1x < op%xmin .or. head1x > op%xmax .or.  &
          head1y < op%ymin .or. head1y > op%ymax .or.  &
          head2x < op%xmin .or. head2x > op%xmax .or.  &
          head2y < op%ymin .or. head2y > op%ymax) cycle
          
          
! Compute re-scaling factor from transformation- STILL NEED THIS??????????

       fac = 1.0

! Draw vector

      call o_frstpt(pointx+fac*(tailx-pointx),pointy+fac*(taily-pointy))
      call o_vector(pointx,pointy)
      call o_frstpt(pointx+fac*(head1x-pointx),pointy+fac*(head1y-pointy))
      call o_vector(pointx,pointy)
      call o_vector(pointx+fac*(head2x-pointx),pointy+fac*(head2y-pointy))

   endif
   
enddo

return
end subroutine vecuslab_h
