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
subroutine icosahedron()

use mem_ijtabs,  only: mrls, itab_m, itab_u, itab_w, alloc_itabs
use mem_grid,    only: mza, nma, nua, nwa, mma, mua, mwa, xem, yem, zem,  &
                       alloc_xyzem, impent
use misc_coms,   only: io6, nxp
use consts_coms, only: pi2, erad, erador5
use mem_nudge, only: nudflag, nudrad, mnudp, xenudp, yenudp, zenudp,  &
                       itab_nudp, alloc_nudge1

implicit none

real, parameter :: pwrd = 1.0  ! 0.9 is close to making uniform-sized triangles
                               ! 1.0 is original value

integer :: ibigd,i,j,idiamond,im_left,iu0,iu1,iu2,iu3,iu4,iw1,iw2,im  &
   ,idiamond_top,im_top,im_right,im_bot,nn10,idiamond_right,idiamond_bot  &
   ,iu,iw,itpn
integer :: id,iter,inudp,inudp1,inudp2

real :: wts,wtn,wtw,wte,expansion,anglen,anglew,anglee,angles,wtw0,wte0,sumwt

integer, save, dimension(10) :: ibigd_ne = (/6,7,8,9,10,7,8,9,10,6/)  &
                               ,ibigd_se = (/2,3,4,5,1,2,3,4,5,1/)

real, save, dimension(10) :: xed_s,xed_n,xed_w,xed_e  &
                            ,yed_s,yed_n,yed_w,yed_e  &
                            ,zed_s,zed_n,zed_w,zed_e

! Temporary scratch arrays for setting up nudging polygons

integer, allocatable :: nudflgm(:),iwnud(:)

! Define triangles, edges, and vertices for icosahedral faces and subdivisions

mrls = 1  ! Default value

! For now, use nxp to divide each face

nn10 = nxp * nxp * 10

! ADD 1 to total number of points needed

nma =     nn10 + 2 + 1  ! ADDING 1 for reference point (index = 1)
nua = 3 * nn10     + 1  ! ADDING 1 for reference point (index = 1)
nwa = 2 * nn10     + 1  ! ADDING 1 for reference point (index = 1)

mma = nma
mua = nua
mwa = nwa

! Allocate temporary arrays if using nudging polygons

if (nudflag > 0 .and. nudrad > 0) then
   allocate (nudflgm(nma),iwnud(nwa))
   nudflgm(1:nma) = 0
endif

! Allocate memory for itabs and M earth coords
! Initialize all neighbor indices to zero

call alloc_itabs(nma,nua,nwa)

call alloc_xyzem()

do im = 2,nma
   itab_m(im)%itopm = im
   call mloops('f',im,1,0,0,0)
enddo

do iu = 2,nua
   itab_u(iu)%iup = iu
   call uloops('f',iu, 1, 4, 5, 7, 8,11,12,13,14,15)
   call uloops('n',iu,16,20, 0, 0, 0, 0, 0, 0, 0, 0)
enddo

do iw = 2,nwa
   itab_w(iw)%iwp = iw
   call wloops('f',iw, 1, 3, 5, 6, 7, 8,11,12,13,14)
   call wloops('n',iw,15,16,17,18,19,20,23,25,26,27)
   call wloops('n',iw,28,29,30,33,34, 0, 0, 0, 0, 0)
enddo

! Fill big diamond corner coordinates

do id = 1,5

   anglen = .2 * (id-1) * pi2
   anglew = anglen - .1 * pi2
   anglee = anglen + .1 * pi2

   zed_s(id) = -erad
   xed_s(id) = 0.
   yed_s(id) = 0.

   zed_n(id) = erador5
   xed_n(id) = erador5 * 2. * cos(anglen)
   yed_n(id) = erador5 * 2. * sin(anglen)

   zed_w(id) = -erador5
   xed_w(id) = erador5 * 2. * cos(anglew)
   yed_w(id) = erador5 * 2. * sin(anglew)

   zed_e(id) = -erador5
   xed_e(id) = erador5 * 2. * cos(anglee)
   yed_e(id) = erador5 * 2. * sin(anglee)

enddo

do id = 6,10

   angles = .2 * (id-6) * pi2 + .1 * pi2
   anglew = angles - .1 * pi2
   anglee = angles + .1 * pi2

   zed_s(id) = -erador5
   xed_s(id) = erador5 * 2. * cos(angles)
   yed_s(id) = erador5 * 2. * sin(angles)

   zed_n(id) = erad
   xed_n(id) = 0.
   yed_n(id) = 0.

   zed_w(id) = erador5
   xed_w(id) = erador5 * 2. * cos(anglew)
   yed_w(id) = erador5 * 2. * sin(anglew)

   zed_e(id) = erador5
   xed_e(id) = erador5 * 2. * cos(anglee)
   yed_e(id) = erador5 * 2. * sin(anglee)

enddo

! Store IM index of south-pole and north-pole pentagonal points

impent(1) = 2
impent(12) = nma

do ibigd = 1,10

   do j = 1,nxp
      do i = 1,nxp

         idiamond = (ibigd - 1) * nxp * nxp  &
                  + (j - 1)     * nxp        &
                  +  i

! Indices that are "attached" to this diamond

         im_left  = idiamond + 2
         
         iu0 = 3 * idiamond
         iu1 = 3 * idiamond - 1
         iu3 = 3 * idiamond + 1
      
         iw1 = 2 * idiamond
         iw2 = 2 * idiamond + 1
        
! Store IM index of 10 out of 12 pentagonal points

         if (i == 1 .and. j == nxp) impent(ibigd+1) = im_left

! Indices that are "attached" to another diamond

         if (ibigd < 6) then   ! Southern 5 diamonds

            ! Top diamond indices

            if (i < nxp) then
               idiamond_top = idiamond + 1
            else
               idiamond_top = (ibigd_ne(ibigd) - 1) * nxp * nxp  &
                            + (j - 1)               * nxp        &
                            +  1
            endif

            im_top   = idiamond_top + 2
            iu4 = 3 * idiamond_top - 1   ! (it's the iu1 for id_top)

            ! Right diamond indices

            if (j > 1 .and. i < nxp) then
               idiamond_right = idiamond - nxp + 1
            elseif (j == 1) then
               idiamond_right = (ibigd_se(ibigd)-1) * nxp * nxp  &
                              + (i - 1)             * nxp        &
                              +  1
               iu2 = 3 * idiamond_right - 1 ! (it's the iu1 for id_right)
            else            ! case for i = nxpand j > 1
               idiamond_right = (ibigd_ne(ibigd)-1) * nxp * nxp  &
                              + (j - 2)             * nxp        &
                              +  1
            endif

            im_right = idiamond_right + 2

            ! Bottom diamond indices
            
            if (j > 1) then
               idiamond_bot = idiamond - nxp
               iu2 = 3 * idiamond_bot + 1 ! (it's the iu3 for id_bot)
            else
               idiamond_bot = (ibigd_se(ibigd)-1) * nxp * nxp  &
                            + (i - 2)             * nxp        &
                            +  1
            endif
            im_bot = idiamond_bot + 2
            
            if (i == 1 .and. j == 1) im_bot = 2

         else                  ! Northern 5 diamonds

           ! Top diamond indices

            if (i < nxp) then
               idiamond_top = idiamond + 1
               iu4 = 3 * idiamond_top - 1   ! (it's the iu1 for id_top)
            else
               idiamond_top = (ibigd_ne(ibigd)-1) * nxp * nxp  &
                            + (nxp - 1)           * nxp        &
                            +  j + 1
            endif

            im_top = idiamond_top + 2

            ! Right diamond indices

            if (j > 1 .and. i < nxp) then
               idiamond_right = idiamond - nxp + 1
            elseif (j == 1 .and. i < nxp) then
               idiamond_right = (ibigd_se(ibigd)-1)  * nxp * nxp  &
                              + (nxp - 1)            * nxp        &
                              + i + 1
            else            ! case for i = nxp
               idiamond_right = (ibigd_ne(ibigd)-1) * nxp * nxp  &
                              + (nxp - 1)           * nxp        &
                              +  j
               iu4 = 3 * idiamond_right + 1 ! (it's the iu3 for id_right)
            endif

            im_right = idiamond_right + 2

            ! Bottom diamond indices
            
            if (j > 1) then
               idiamond_bot = idiamond - nxp
            else
               idiamond_bot = (ibigd_se(ibigd)-1) * nxp * nxp  &
                            + (nxp - 1)           * nxp        &
                            + i 
            endif

            im_bot = idiamond_bot + 2
            iu2 = 3 * idiamond_bot + 1 ! (it's the iu3 for id_bot)
            
            if (i == nxp .and. j == nxp)  &
               im_top = 10 * nxp * nxp + 3

         endif

         call fill_diamond(im_left,im_right,im_top,im_bot  &
            ,iu0,iu1,iu2,iu3,iu4,iw1,iw2)
            
!!!!!!!!!!!!!!!!!!!!! new section for nudging points !!!!!!!!!!!!!!!!!
         if (nudflag > 0 .and. nudrad > 0) then
            if (mod(abs(j-i),3*nudrad) == 0) then
               if (mod(j-1,nudrad) == 0)    nudflgm(im_bot) = 1
               if (j == nxp .or. i == nxp)  nudflgm(im_top) = 1
            endif
         endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

! M point (xem,yem,zem) coordinates
         
         if (i + j <= nxp) then
            wts  = max(0.,min(1.,real(nxp + 1 - i - j) / real(nxp)))
            wtn  = 0.
            wtw0 = max(0.,min(1.,real(j) / real(i + j - 1)))
            wte0 = 1. - wtw0	
         else
            wts  = 0.
            wtn  = max(0.,min(1.,real(i + j - nxp - 1) / real(nxp)))
            wte0 = max(0.,min(1.,real(nxp - j)  &
                 / real(2 * nxp + 1 - i - j)))
            wtw0 = 1. - wte0	
         endif

! Experimental adjustment in spacing
! Compute sum of weights raised to pwrd

         wtw = (1. - wts - wtn) * wtw0
         wte = (1. - wts - wtn) * wte0

         sumwt = wts**pwrd + wtn**pwrd + wtw**pwrd + wte**pwrd

         wts = wts**pwrd / sumwt
         wtn = wtn**pwrd / sumwt
         wtw = wtw**pwrd / sumwt
         wte = wte**pwrd / sumwt

         
         xem(im_left) = wts * xed_s(ibigd)  &
                      + wtn * xed_n(ibigd)  &
                      + wtw * xed_w(ibigd)  &
                      + wte * xed_e(ibigd)

         yem(im_left) = wts * yed_s(ibigd)  &
                      + wtn * yed_n(ibigd)  &
                      + wtw * yed_w(ibigd)  &
                      + wte * yed_e(ibigd)

         zem(im_left) = wts * zed_s(ibigd)  &
                      + wtn * zed_n(ibigd)  &
                      + wtw * zed_w(ibigd)  &
                      + wte * zed_e(ibigd)

! Push M point coordinates out to earth radius

         expansion = erad / sqrt(xem(im_left) ** 2  &
                               + yem(im_left) ** 2  &
                               + zem(im_left) ** 2  )

         xem(im_left) = xem(im_left) * expansion
         yem(im_left) = yem(im_left) * expansion
         zem(im_left) = zem(im_left) * expansion
         
      enddo  ! end i loop
   enddo     ! end j loop

enddo        ! end idbig loop

xem(2) = 0.
yem(2) = 0.
zem(2) = -erad

xem(nma) = 0.
yem(nma) = 0.
zem(nma) = erad

!call twist()

call tri_neighbors()

! This is the place to do spring dynamics

call spring_dynamics()

! Begin counter for total number of nudging polygons (centered at M points)
! (start with 1 for dummy point)

mnudp = 1

! Define various quantities for polygon nudging if option is activated

if (nudflag > 0 .and. nudrad > 0) then

! Count number of M points that are flagged as nudge points

   inudp = 1
   do im = 2,mma
      if (nudflgm(im) == 1) inudp = inudp + 1
   enddo

! Allocate arrays for nudging groups

   mnudp = inudp
   call alloc_nudge1()

! Loop again over M points and check which are flagged as nudge points

   inudp = 1
   do im = 2,mma
      if (nudflgm(im) == 1) then
         inudp = inudp + 1

! Assign xe,ye,ze coordinates of nudging groups

         xenudp(inudp) = xem(im)
         yenudp(inudp) = yem(im)
         zenudp(inudp) = zem(im)

!      write(io6,40) im,inudp,xenudp(inudp), yenudp(inudp),zenudp(inudp)
!   40 format('im,inudp,xyz ',2i6,3e11.3)

! Assign nudge-point index to all W points that surround a nudge-flagged M point

         do itpn = 1,itab_m(im)%ntpn
            iw = itab_m(im)%iw(itpn)
            itab_w(iw)%inudp(1) = inudp
         enddo

      endif
   enddo

! Copy nudge-point indices to adjacent W points for nudrad > 1 cases

   do iw = 2,mwa
      iwnud(iw) = itab_w(iw)%inudp(1)
   enddo

   do iter = 1,2*(nudrad-1)
      do iu = 2,mua
         iw1 = itab_u(iu)%iw1
         iw2 = itab_u(iu)%iw2

         if (itab_w(iw1)%inudp(1) > 1) iwnud(iw2) = itab_w(iw1)%inudp(1)
         if (itab_w(iw2)%inudp(1) > 1) iwnud(iw1) = itab_w(iw2)%inudp(1)
      enddo

      do iw = 2,mwa
         itab_w(iw)%inudp(1) = iwnud(iw)
      enddo
   enddo

! We no longer need nudflgm or iwnud

   deallocate (nudflgm,iwnud)

! Now, the inudp(1) index is defined for all W points on global grid 1

! Search through all U points for those that have a different itab_w()%inudp(1)
! value on either side

   do iu = 2,mua
      iw1 = itab_u(iu)%iw1
      iw2 = itab_u(iu)%iw2

      inudp1 = itab_w(iw1)%inudp(1)
      inudp2 = itab_w(iw2)%inudp(1)

      if (inudp1 /= inudp2) then

!  If this condition is satisfied, exchange index information betwen nudp
!  points and add to itab_nudp array

         do i = 1,6
            if (itab_nudp(inudp1,i) == 1    ) itab_nudp(inudp1,i) = inudp2
            if (itab_nudp(inudp1,i) == inudp2) exit
         enddo

         do i = 1,6
            if (itab_nudp(inudp2,i) == 1    ) itab_nudp(inudp2,i) = inudp1
            if (itab_nudp(inudp2,i) == inudp1) exit
         enddo

      endif

   enddo

endif

return
end subroutine icosahedron

!===============================================================================

subroutine fill_diamond(im_left,im_right,im_top,im_bot  &
   ,iu0,iu1,iu2,iu3,iu4,iw1,iw2)

use mem_ijtabs, only: itab_u, itab_w
use misc_coms,  only: io6

implicit none

integer, intent(in) :: im_left,im_right,im_top,im_bot
integer, intent(in) :: iu0,iu1,iu2,iu3,iu4,iw1,iw2     

itab_u(iu0)%im1 = im_left
itab_u(iu0)%im2 = im_right
itab_u(iu0)%iw1 = iw1
itab_u(iu0)%iw2 = iw2
itab_u(iu0)%mrlu = 1
         
itab_u(iu1)%im1 = im_left
itab_u(iu1)%im2 = im_bot
itab_u(iu1)%iw2 = iw1
itab_u(iu1)%mrlu = 1
         
itab_u(iu2)%iw1 = iw1

itab_u(iu3)%im1 = im_top
itab_u(iu3)%im2 = im_left
itab_u(iu3)%iw2 = iw2
itab_u(iu3)%mrlu = 1

itab_u(iu4)%iw1 = iw2
         
itab_w(iw1)%iu1 = iu0 
itab_w(iw1)%iu2 = iu1 
itab_w(iw1)%iu3 = iu2 
itab_w(iw1)%mrlw = 1
itab_w(iw1)%mrlw_orig = 1
      
itab_w(iw2)%iu1 = iu0 
itab_w(iw2)%iu2 = iu4 
itab_w(iw2)%iu3 = iu3 
itab_w(iw2)%mrlw = 1
itab_w(iw2)%mrlw_orig = 1

return
end subroutine fill_diamond
   
!===============================================================================

subroutine twist()

use mem_ijtabs,  only: itab_m, itab_u, itab_w
use mem_grid,    only: nma, nua, xem, yem, zem
use misc_coms,   only: io6, nxp
use consts_coms, only: pi1, pi2

implicit none

integer :: im,im1,im2,im1_east,im2_east
integer :: iu,iu1,iu2,iu3,iu_east,iueq,iueq_east,iueq_fill
integer :: iw,iw1

integer :: nxpo2  ! half of nxp (nxp must be even if twist is done)

integer :: iueq_tab(10*nxp)
integer :: iweq_tab(10*nxp)

real :: radm  ! Distance of M point from Earth axis [m]
real :: angm  ! longitude (radians) of M point before shift is applied

nxpo2 = nxp / 2

! Initialize iueq_tab and iweq_tab to zero

iueq_tab(:) = 0
iweq_tab(:) = 0

! Shift all M points in southern hemisphere (excluding equator) 
! 36 degrees to the east

do im = 2,nma
   if (zem(im) < -100.) then

      radm = sqrt(xem(im) ** 2 + yem(im) ** 2)
      angm = atan2(yem(im),xem(im))

      xem(im) = radm * cos(angm + .2 * pi1)
      yem(im) = radm * sin(angm + .2 * pi1)

   endif
enddo

! Loop over all U points and search for those that are on the equator

do iu = 1,nua
   im1 = itab_u(iu)%im1
   im2 = itab_u(iu)%im2

   if (abs(zem(im1)) < 100. .and. abs(zem(im2)) < 100.) then

! Compute special index for equatorial U points that increases with longitude.
! This uses knowledge of icosahedron subroutine that im2 is always east of iu.
! IUEQ skips approximately every other integer and therefore needs to be collapsed.

      iueq = int(10 * nxp * (atan2(yem(im2),xem(im2)) + pi1) / pi2) + 1   

! Fill table of iu and iw1 indices for each iueq

      iueq_tab(iueq) = iu
      iweq_tab(iueq) = itab_u(iu)%iw1

   endif
enddo

! Collapse IUEQ tables in order to use consecutive integer indices

iueq_fill = 0

do iueq = 1,nxp * 10

   if (iueq_tab(iueq) > 1) then
      iueq_fill = iueq_fill + 1

      iueq_tab(iueq_fill) = iueq_tab(iueq)
      iweq_tab(iueq_fill) = iweq_tab(iueq)
   endif

enddo

! Loop over equatorial U points

do iueq = 1,nxp * 5

! Find equatorial index of U point that is 36 degrees to the east

   iueq_east = iueq + nxpo2
   if (iueq_east > nxp * 5) iueq_east = iueq_east - nxp * 5

! Get U point indices for iueq and iueq_east points

   iu      = iueq_tab(iueq)
   iu_east = iueq_tab(iueq_east)

! Get index for adjacent W point to the south of iu

   iw = iweq_tab(iueq)

! Get indices for M endpoints of iu and iu_east

   im1 = itab_u(iu)%im1
   im2 = itab_u(iu)%im2

   im1_east = itab_u(iu_east)%im1
   im2_east = itab_u(iu_east)%im2

! Get U neighbor indices for adjacent W point to the south

   iu1 = itab_w(iw)%iu1
   iu2 = itab_w(iw)%iu2
   iu3 = itab_w(iw)%iu3

! For U points on equator: shift their south W point neighbor index 

   itab_u(iu_east)%iw1 = iw

! For W points bordering equator on south: shift their equatorial U point
! neighbor index

   if     (iu1 == iu) then
      itab_w(iw)%iu1 = iu_east
   elseif (iu2 == iu) then
      itab_w(iw)%iu2 = iu_east
   elseif (iu3 == iu) then
      itab_w(iw)%iu3 = iu_east
   endif

! For U points bordering equator on south: shift their equatorial M point 
! neighbor index

   if (iu1 /= iu) then
      if     (itab_u(iu1)%im1 == im1) then
              itab_u(iu1)%im1 =  im1_east
      elseif (itab_u(iu1)%im2 == im1) then
              itab_u(iu1)%im2 =  im1_east
      elseif (itab_u(iu1)%im1 == im2) then
              itab_u(iu1)%im1 =  im2_east
      elseif (itab_u(iu1)%im2 == im2) then
              itab_u(iu1)%im2 =  im2_east
      endif
   endif

   if (iu2 /= iu) then
      if     (itab_u(iu2)%im1 == im1) then
              itab_u(iu2)%im1 =  im1_east
      elseif (itab_u(iu2)%im2 == im1) then
              itab_u(iu2)%im2 =  im1_east
      elseif (itab_u(iu2)%im1 == im2) then
              itab_u(iu2)%im1 =  im2_east
      elseif (itab_u(iu2)%im2 == im2) then
              itab_u(iu2)%im2 =  im2_east
      endif
   endif

   if (iu3 /= iu) then
      if     (itab_u(iu3)%im1 == im1) then
              itab_u(iu3)%im1 =  im1_east
      elseif (itab_u(iu3)%im2 == im1) then
              itab_u(iu3)%im2 =  im1_east
      elseif (itab_u(iu3)%im1 == im2) then
              itab_u(iu3)%im1 =  im2_east
      elseif (itab_u(iu3)%im2 == im2) then
              itab_u(iu3)%im2 =  im2_east
      endif
   endif

enddo

return
end subroutine twist

