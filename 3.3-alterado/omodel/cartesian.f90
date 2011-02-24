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
subroutine cartesian_2d()

use mem_ijtabs, only: itab_m, itab_u, itab_w, mrls, alloc_itabs
use mem_grid,   only: nma, nua, nwa, mma, mua, mwa, xem, yem, zem, alloc_xyzem
use misc_coms,  only: io6, nxp, mdomain, deltax

implicit none

integer :: i,im,iu,iw,im1,im2,im3,im4,iu1,iu2,iu3,iu4,iu5,iw0,iw1,iw2,iw3
real :: unit_dist,diamond_centx

! Define triangles, edges, and vertices for 2D cartesian channel domain
! THIS CASE APPLIES FOR MDOMAIN = 2 (open x bnds) or 3 (cyclic x bnds)

mrls = 1  ! Default value

! Use nxp to count triangles

nma =     nxp + nxp + 2 + 1 ! ADDING 1 for reference point (index = 1)
nua = 3 * nxp + nxp + 1 + 1 ! ADDING 1 for reference point (index = 1)
nwa = 2 * nxp           + 1 ! ADDING 1 for reference point (index = 1)

mma = nma
mua = nua
mwa = nwa

call alloc_itabs(nma,nua,nwa)

call alloc_xyzem()

do im = 2,nma
   itab_m(im)%itopm = im
   call mloops('f',im,1,0,0,0)
enddo

do i = 1,nxp

   im1 = 2 * i
   im2 = 2 * i + 1
   im3 = 2 * i + 2
   im4 = 2 * i + 3

   iu1 = 4 * i - 2
   iu2 = 4 * i - 1
   iu3 = 4 * i
   iu4 = 4 * i + 1
   iu5 = 4 * i + 2

   iw0 = 2 * i - 1
   iw1 = 2 * i
   iw2 = 2 * i + 1
   iw3 = 2 * i + 2

   if (i == 1  ) iw0 = 1
   if (i == nxp) iw3 = 1

! IU1

   itab_u(iu1)%im1 = im2
   itab_u(iu1)%im2 = im1
   itab_u(iu1)%iw1 = iw0
   itab_u(iu1)%iw2 = iw1

   if (i == 1) then
      if (mdomain == 2) then
         itab_u(iu1)%iup = iu1+4
         call uloops('f',iu1,1,5,7,8,11,14,15, 0,0,0)
      elseif (mdomain == 3) then
         itab_u(iu1)%iup = 4*nxp-6
         call uloops('f',iu1,1,5,7,8,11,16,17,29,0,0)
      endif
   elseif (i == nxp) then
      if (mdomain == 2) then
         itab_u(iu1)%iup = iu1
         call uloops('f',iu1, 1, 4, 5, 7, 8,11,12,13,14,15)
         call uloops('n',iu1,16,20, 0, 0, 0, 0, 0, 0, 0, 0)
      elseif (mdomain == 3) then
         itab_u(iu1)%iup = 6
         call uloops('f',iu1, 1, 4, 5, 7, 8,11,12,14,15,16)
         call uloops('n',iu1,18,20, 0, 0, 0, 0, 0, 0, 0, 0)
      endif
   else
      itab_u(iu1)%iup = iu1
      call uloops('f',iu1, 1, 4, 5, 7, 8,11,12,14,15,16)
      call uloops('n',iu1,20, 0, 0, 0, 0, 0, 0, 0, 0, 0)
   endif

! IU2

   itab_u(iu2)%im1 = im1
   itab_u(iu2)%im2 = im3
   itab_u(iu2)%iw1 = 1
   itab_u(iu2)%iw2 = iw1
   itab_u(iu2)%iup = iu2

   call uloops('f',iu2,3,11,14,15,0,0,0,0,0,0)

! IU3

   itab_u(iu3)%im1 = im2
   itab_u(iu3)%im2 = im3
   itab_u(iu3)%iw1 = iw1
   itab_u(iu3)%iw2 = iw2

   if (i == 1) then
      if (mdomain == 2) then
         itab_u(iu3)%iup = iu3+4
         call uloops('f',iu3, 1, 4, 5, 7, 8,11,13,14,15,17)
         call uloops('n',iu3,20, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      elseif (mdomain == 3) then
         itab_u(iu3)%iup = 4*nxp-4
         call uloops('f',iu3, 1, 4, 5, 7, 8,11,14,15,18, 0)
      endif
   elseif (i == nxp) then
      if (mdomain == 2) then
         itab_u(iu3)%iup = iu3-4
         call uloops('f',iu3, 1, 4, 5, 7, 8,11,13,14,15,19)
         call uloops('n',iu3,20, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      elseif (mdomain == 3) then
         itab_u(iu3)%iup = 8
         call uloops('f',iu3, 1, 4, 5, 7, 8,11,14,15,18, 0)
      endif
   else
      itab_u(iu3)%iup = iu3
      call uloops('f',iu3, 1, 4, 5, 7, 8,11,12,14,15,16)
      call uloops('n',iu3,20, 0, 0, 0, 0, 0, 0, 0, 0, 0)
   endif

! IU4

   itab_u(iu4)%im1 = im2
   itab_u(iu4)%im2 = im4
   itab_u(iu4)%iw1 = iw2
   itab_u(iu4)%iw2 = 1
   itab_u(iu4)%iup = iu4

   call uloops('f',iu4, 3,11,14,15, 0, 0, 0, 0, 0, 0)

! IU5

   if (i == nxp) then
      itab_u(iu5)%im1 = im4
      itab_u(iu5)%im2 = im3
      itab_u(iu5)%iw1 = iw2
      itab_u(iu5)%iw2 = iw3

      if (mdomain == 2) then
         itab_u(iu5)%iup = iu5-4
         call uloops('f',iu5, 1, 5, 7, 8,11,14,15,18, 0, 0)
      elseif (mdomain == 3) then
         itab_u(iu5)%iup = 10
         call uloops('f',iu5, 1, 5, 7, 8,11,14,15,18, 0, 0)
      endif
   endif

! IW1

   itab_w(iw1)%iu1 = iu1
   itab_w(iw1)%iu2 = iu2
   itab_w(iw1)%iu3 = iu3
   itab_w(iw1)%mrlw = 1
   itab_w(iw1)%mrlw_orig = 1

   if (i == 1) then
      if (mdomain == 2) then
         itab_w(iw1)%iwp = iw2
         call wloops('f',iw1, 1, 3, 5, 6, 7, 8,11,14,18,22)
         call wloops('n',iw1,23,32,33, 0, 0, 0, 0, 0, 0, 0)
      elseif (mdomain == 3) then
         itab_w(iw1)%iwp = 2*nxp-2
         call wloops('f',iw1, 1, 3, 5, 6, 7, 8,11,14,18,22)
         call wloops('n',iw1,24,31,33,35, 0, 0, 0, 0, 0, 0)
      endif
   elseif (i == nxp) then
      if (mdomain == 2) then
         itab_w(iw1)%iwp = iw1
         call wloops('f',iw1, 1, 3, 5, 6, 7, 8,11,12,13,14)
         call wloops('n',iw1,15,16,17,18,19,20,25,26,27,28)
         call wloops('n',iw1,29,30,33,34, 0, 0, 0, 0, 0, 0)
      elseif (mdomain == 3) then
         itab_w(iw1)%iwp = 4
         call wloops('f',iw1, 1, 3, 5, 6, 7, 8,11,14,18,22)
         call wloops('n',iw1,24,31,33,35, 0, 0, 0, 0, 0, 0)
      endif
   else
      itab_w(iw1)%iwp = iw1
      call wloops('f',iw1, 1, 3, 5, 6, 7, 8,11,12,13,14)
      call wloops('n',iw1,16,17,18,19,25,26,27,28,29,30)
      call wloops('n',iw1,33,34, 0, 0, 0, 0, 0, 0, 0, 0)
   endif

! IW2

   itab_w(iw2)%iu1 = iu3
   itab_w(iw2)%iu2 = iu5
   itab_w(iw2)%iu3 = iu4
   itab_w(iw2)%mrlw = 1
   itab_w(iw2)%mrlw_orig = 1
   if (i == 1) then
      if (mdomain == 2) then
         itab_w(iw2)%iwp = iw2
         call wloops('f',iw2, 1, 3, 5, 6, 7, 8,11,12,13,14)
         call wloops('n',iw2,15,16,17,18,19,20,25,26,27,28)
         call wloops('n',iw2,29,30,33,34, 0, 0, 0, 0, 0, 0)
      elseif (mdomain == 3) then
         itab_w(iw2)%iwp = 2*nxp-1
         call wloops('f',iw2, 1, 3, 5, 6, 7, 8,11,14,18,22)
         call wloops('n',iw2,24,31,33,35, 0, 0, 0, 0, 0, 0)
      endif
   elseif (i == nxp) then
      if (mdomain == 2) then
         itab_w(iw2)%iwp = iw1
         call wloops('f',iw2, 1, 3, 5, 6, 7, 8,11,14,18,19)
         call wloops('n',iw2,21,24,31,35, 0, 0, 0, 0, 0, 0)
      elseif (mdomain == 3) then
         itab_w(iw2)%iwp = 5
         call wloops('f',iw2, 1, 3, 5, 6, 7, 8,11,14,18,22)
         call wloops('n',iw2,24,31,33,35, 0, 0, 0, 0, 0, 0)
      endif
   else
      itab_w(iw2)%iwp = iw2
      call wloops('f',iw2, 1, 3, 5, 6, 7, 8,11,12,13,14)
      call wloops('n',iw2,16,17,18,19,25,26,27,28,29,30)
      call wloops('n',iw2,33,34, 0, 0, 0, 0, 0, 0, 0, 0)
   endif

   unit_dist = .5 * sqrt(3.) * deltax ! This is 1/2 of triangle face width
   diamond_centx = (2 * i - 1 - nxp) * unit_dist

   xem(im1) = diamond_centx - 1.5 * unit_dist
   xem(im2) = diamond_centx -  .5 * unit_dist
   xem(im3) = diamond_centx +  .5 * unit_dist
   xem(im4) = diamond_centx + 1.5 * unit_dist

   yem(im1) = - .75 * deltax
   yem(im2) =   .75 * deltax
   yem(im3) = - .75 * deltax
   yem(im4) =   .75 * deltax

enddo  ! end i loop

zem(2:nma) = 0.

call tri_neighbors()

return
end subroutine cartesian_2d

!=============================================================================

subroutine cartesian_3d()

use mem_ijtabs, only: itab_m, itab_u, itab_w, mrls, alloc_itabs
use mem_grid,   only: nma, nua, nwa, mma, mua, mwa, xem, yem, zem, alloc_xyzem
use misc_coms,  only: io6, nxp, mdomain, deltax

implicit none

integer :: i,im,iu,iw,iucyc,iwcyc
integer :: im1,im2,im3,im4,im5,im6,im7,im8,im9,im10
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9,iu10,iu11,iu12,iu13,iu14,iu15
integer :: iw1,iw2,iw3,iw4,iw5,iw6,iw7,iw8,iw9,iw10,iw11,iw12,iw13,iw14,iw15
integer :: iu16,iu17,iw16

real :: unit_dist,diamond_centx

! Define triangles, edges, and vertices for 3D cartesian channel domain
! Channel width = 4 allocated rows with cyclic repeat distance 2 rows wide
! THIS CASE APPLIES FOR MDOMAIN = 4 (cyclic x bnds & cyclic y bnds)

mrls = 1  ! Default value

! Use nxp to count triangles

nma =  5 * nxp + 5 + 1 ! ADDING 1 for reference point (index = 1)
nua = 13 * nxp + 4 + 1 ! ADDING 1 for reference point (index = 1)
nwa =  8 * nxp     + 1 ! ADDING 1 for reference point (index = 1)

mma = nma
mua = nua
mwa = nwa

call alloc_itabs(nma,nua,nwa)

call alloc_xyzem()

do im = 2,nma
   itab_m(im)%itopm = im
   call mloops('f',im,1,0,0,0)
enddo

do i = 1,nxp

   im1  = 5 * i - 3
   im2  = 5 * i - 2
   im3  = 5 * i - 1
   im4  = 5 * i
   im5  = 5 * i + 1
   im6  = 5 * i + 2
   im7  = 5 * i + 3
   im8  = 5 * i + 4
   im9  = 5 * i + 5
   im10 = 5 * i + 6

   iu1  = 13 * i - 11
   iu2  = 13 * i - 10
   iu3  = 13 * i - 9
   iu4  = 13 * i - 8
   iu5  = 13 * i - 7
   iu6  = 13 * i - 6
   iu7  = 13 * i - 5
   iu8  = 13 * i - 4
   iu9  = 13 * i - 3
   iu10 = 13 * i - 2
   iu11 = 13 * i - 1
   iu12 = 13 * i
   iu13 = 13 * i + 1
   iu14 = 13 * i + 2
   iu15 = 13 * i + 3
   iu16 = 13 * i + 4
   iu17 = 13 * i + 5

   iw1  = 8 * i - 10
   iw2  = 8 * i - 9
   iw3  = 8 * i - 8
   iw4  = 8 * i - 7
   iw5  = 8 * i - 6
   iw6  = 8 * i - 5
   iw7  = 8 * i - 4
   iw8  = 8 * i - 3
   iw9  = 8 * i - 2
   iw10 = 8 * i - 1
   iw11 = 8 * i
   iw12 = 8 * i + 1
   iw13 = 8 * i + 2
   iw14 = 8 * i + 3
   iw15 = 8 * i + 4
   iw16 = 8 * i + 5

   if (i == 1) then
      iw1 = 1
      iw2 = 1
      iw3 = 1
      iw4 = 1
   endif
   
   if (i == nxp) then
      iw13 = 1
      iw14 = 1
      iw15 = 1
      iw16 = 1
   endif

! Cyclic repeat distance along channel

   iucyc = 13 * (nxp - 2)
   iwcyc =  8 * (nxp - 2)

! IU1

   itab_u(iu1)%im1 = im4
   itab_u(iu1)%im2 = im1
   itab_u(iu1)%iw1 = iw1
   itab_u(iu1)%iw2 = iw5

   if (i == 1) then
      itab_u(iu1)%iup = iu3+iucyc
      call uloops('f',iu1, 1, 5, 7, 8,11,14,15,18, 0, 0)
   elseif (i == nxp) then
      itab_u(iu1)%iup = iu3-iucyc
      call uloops('f',iu1, 1, 4, 5, 7, 8,11,14,15,18, 0)
   else
      itab_u(iu1)%iup = iu3
      call uloops('f',iu1, 1, 4, 5, 7, 8,11,14,15,18, 0)
   endif

! IU2

   itab_u(iu2)%im1 = im2
   itab_u(iu2)%im2 = im4
   itab_u(iu2)%iw1 = iw2
   itab_u(iu2)%iw2 = iw6

   if (i == 1) then
      itab_u(iu2)%iup = iu2+iucyc
      call uloops('f',iu2, 1, 5, 7, 8,11,14,15,18, 0, 0)
   elseif (i == nxp) then
      itab_u(iu2)%iup = iu2-iucyc
      call uloops('f',iu2, 1, 4, 5, 7, 8,11,12,14,15,16)
      call uloops('n',iu2,18,20, 0, 0, 0, 0, 0, 0, 0, 0)
   else
      itab_u(iu2)%iup = iu2
      call uloops('f',iu2, 1, 4, 5, 7, 8,11,12,14,15,16)
      call uloops('n',iu2,20, 0, 0, 0, 0, 0, 0, 0, 0, 0)
   endif

! IU3

   itab_u(iu3)%im1 = im5
   itab_u(iu3)%im2 = im2
   itab_u(iu3)%iw1 = iw3
   itab_u(iu3)%iw2 = iw7

   if (i == 1) then
      itab_u(iu3)%iup = iu3+iucyc
      call uloops('f',iu3, 1, 5, 7, 8,11,14,15,18, 0, 0)
   elseif (i == nxp) then
      itab_u(iu3)%iup = iu3-iucyc
      call uloops('f',iu3, 1, 4, 5, 7, 8,11,12,14,15,16)
      call uloops('n',iu3,18,20, 0, 0, 0, 0, 0, 0, 0, 0)
   else
      itab_u(iu3)%iup = iu3
      call uloops('f',iu3, 1, 4, 5, 7, 8,11,12,14,15,16)
      call uloops('n',iu3,20, 0, 0, 0, 0, 0, 0, 0, 0, 0)
   endif

! IU4

   itab_u(iu4)%im1 = im3
   itab_u(iu4)%im2 = im5
   itab_u(iu4)%iw1 = iw4
   itab_u(iu4)%iw2 = iw8

   if (i == 1) then
      itab_u(iu4)%iup = iu2+iucyc
      call uloops('f',iu4, 1, 5, 7, 8,11,14,15,18, 0, 0)
   elseif (i == nxp) then
      itab_u(iu4)%iup = iu2-iucyc
      call uloops('f',iu4, 1, 4, 5, 7, 8,11,14,15,18, 0)
   else
      itab_u(iu4)%iup = iu2
      call uloops('f',iu4, 1, 4, 5, 7, 8,11,14,15,18, 0)
   endif

! IU5

   itab_u(iu5)%im1 = im1
   itab_u(iu5)%im2 = im6
   itab_u(iu5)%iw1 = 1
   itab_u(iu5)%iw2 = iw5

   if (i == 1) then
      itab_u(iu5)%iup = iu6+iucyc
      call uloops('f',iu5, 1, 4, 5, 7, 8,11,14,15,18, 0)
   elseif (i == nxp) then
      itab_u(iu5)%iup = iu6-iucyc
      call uloops('f',iu5, 1, 4, 5, 7, 8,11,14,15,18, 0)
   else
      itab_u(iu5)%iup = iu6
      call uloops('f',iu5, 1, 4, 5, 7, 8,11,14,15,18, 0)
   endif

! IU6

   itab_u(iu6)%im1 = im2
   itab_u(iu6)%im2 = im7
   itab_u(iu6)%iw1 = iw6
   itab_u(iu6)%iw2 = iw7

   if (i == 1) then
      itab_u(iu6)%iup = iu6+iucyc
      call uloops('f',iu6, 1, 4, 5, 7, 8,11,14,15,18, 0)
   elseif (i == nxp) then
      itab_u(iu6)%iup = iu6-iucyc
      call uloops('f',iu6, 1, 4, 5, 7, 8,11,14,15,18, 0)
   else
      itab_u(iu6)%iup = iu6
      call uloops('f',iu6, 1, 4, 5, 7, 8,11,12,14,15,16)
      call uloops('n',iu6,20, 0, 0, 0, 0, 0, 0, 0, 0, 0)
   endif

! IU7

   itab_u(iu7)%im1 = im3
   itab_u(iu7)%im2 = im8
   itab_u(iu7)%iw1 = iw8
   itab_u(iu7)%iw2 = 1

   if (i == 1) then
      itab_u(iu7)%iup = iu6+iucyc
      call uloops('f',iu7, 1, 4, 5, 7, 8,11,14,15,18, 0)
   elseif (i == nxp) then
      itab_u(iu7)%iup = iu6-iucyc
      call uloops('f',iu7, 1, 4, 5, 7, 8,11,14,15,18, 0)
   else
      itab_u(iu7)%iup = iu6
      call uloops('f',iu7, 1, 4, 5, 7, 8,11,14,15,18, 0)
   endif

! IU8

   itab_u(iu8)%im1 = im4
   itab_u(iu8)%im2 = im6
   itab_u(iu8)%iw1 = iw5
   itab_u(iu8)%iw2 = iw9

   if (i == 1) then
      itab_u(iu8)%iup = iu10+iucyc
      call uloops('f',iu8, 1, 4, 5, 7, 8,11,14,15,18, 0)
   elseif (i == nxp) then
      itab_u(iu8)%iup = iu10-iucyc
      call uloops('f',iu8, 1, 4, 5, 7, 8,11,14,15,18, 0)
   else
      itab_u(iu8)%iup = iu10
      call uloops('f',iu8, 1, 4, 5, 7, 8,11,14,15,18, 0)
   endif

! IU9

   itab_u(iu9)%im1 = im7
   itab_u(iu9)%im2 = im4
   itab_u(iu9)%iw1 = iw6
   itab_u(iu9)%iw2 = iw10

   if (i == 1) then
      itab_u(iu9)%iup = iu9+iucyc
      call uloops('f',iu9, 1, 4, 5, 7, 8,11,14,15,18, 0)
   elseif (i == nxp) then
      itab_u(iu9)%iup = iu9-iucyc
      call uloops('f',iu9, 1, 4, 5, 7, 8,11,14,15,18, 0)
   else
      itab_u(iu9)%iup = iu9
      call uloops('f',iu9, 1, 4, 5, 7, 8,11,12,14,15,16)
      call uloops('n',iu9,20, 0, 0, 0, 0, 0, 0, 0, 0, 0)
   endif 

! IU10

   itab_u(iu10)%im1 = im5
   itab_u(iu10)%im2 = im7
   itab_u(iu10)%iw1 = iw7
   itab_u(iu10)%iw2 = iw11

   if (i == 1) then
      itab_u(iu10)%iup = iu10+iucyc
      call uloops('f',iu10, 1, 4, 5, 7, 8,11,14,15,18, 0)
   elseif (i == nxp) then
      itab_u(iu10)%iup = iu10-iucyc
      call uloops('f',iu10, 1, 4, 5, 7, 8,11,14,15,18, 0)
   else
      itab_u(iu10)%iup = iu10
      call uloops('f',iu10, 1, 4, 5, 7, 8,11,12,14,15,16)
      call uloops('n',iu10,20, 0, 0, 0, 0, 0, 0, 0, 0, 0)
   endif

! IU11

   itab_u(iu11)%im1 = im8
   itab_u(iu11)%im2 = im5
   itab_u(iu11)%iw1 = iw8
   itab_u(iu11)%iw2 = iw12

   if (i == 1) then
      itab_u(iu11)%iup = iu9+iucyc
      call uloops('f',iu11, 1, 4, 5, 7, 8,11,14,15,18, 0)
   elseif (i == nxp) then
      itab_u(iu11)%iup = iu9-iucyc
      call uloops('f',iu11, 1, 4, 5, 7, 8,11,14,15,18, 0)
   else
      itab_u(iu11)%iup = iu9
      call uloops('f',iu11, 1, 4, 5, 7, 8,11,14,15,18, 0)
   endif

! IU12

   itab_u(iu12)%im1 = im4
   itab_u(iu12)%im2 = im9
   itab_u(iu12)%iw1 = iw9
   itab_u(iu12)%iw2 = iw10

   if (i == 1) then
      itab_u(iu12)%iup = iu12+iucyc
      call uloops('f',iu12, 1, 4, 5, 7, 8,11,14,15,18, 0)
   elseif (i == nxp) then
      itab_u(iu12)%iup = iu12-iucyc
      call uloops('f',iu12, 1, 4, 5, 7, 8,11,14,15,18, 0)
   else
      itab_u(iu12)%iup = iu12
      call uloops('f',iu12, 1, 4, 5, 7, 8,11,12,14,15,16)
      call uloops('n',iu12,20, 0, 0, 0, 0, 0, 0, 0, 0, 0)
   endif

! IU13

   itab_u(iu13)%im1 = im5
   itab_u(iu13)%im2 = im10
   itab_u(iu13)%iw1 = iw11
   itab_u(iu13)%iw2 = iw12

   if (i == 1) then
      itab_u(iu13)%iup = iu12+iucyc
      call uloops('f',iu13, 1, 4, 5, 7, 8,11,14,15,18, 0)
   elseif (i == nxp) then
      itab_u(iu13)%iup = iu12-iucyc
      call uloops('f',iu13, 1, 4, 5, 7, 8,11,14,15,18, 0)
   else
      itab_u(iu13)%iup = iu12
      call uloops('f',iu13, 1, 4, 5, 7, 8,11,14,15,18, 0)
   endif

! IU14

   if (i == nxp) then
      itab_u(iu14)%im1 = im9
      itab_u(iu14)%im2 = im6
      itab_u(iu14)%iw1 = iw9
      itab_u(iu14)%iw2 = iw13

      itab_u(iu14)%iup = iu16-iucyc
      call uloops('f',iu14, 1, 5, 7, 8,11,14,15,18, 0, 0)
   endif

! IU15

   if (i == nxp) then
      itab_u(iu15)%im1 = im7
      itab_u(iu15)%im2 = im9
      itab_u(iu15)%iw1 = iw10
      itab_u(iu15)%iw2 = iw14

      itab_u(iu15)%iup = iu15-iucyc
      call uloops('f',iu15, 1, 5, 7, 8,11,14,15,18, 0, 0)
   endif

! IU16

   if (i == nxp) then
      itab_u(iu16)%im1 = im10
      itab_u(iu16)%im2 = im7
      itab_u(iu16)%iw1 = iw11
      itab_u(iu16)%iw2 = iw15

      itab_u(iu16)%iup = iu13-iucyc
      call uloops('f',iu16, 1, 5, 7, 8,11,14,15,18, 0, 0)
   endif

! IU17

   if (i == nxp) then
      itab_u(iu17)%im1 = im8
      itab_u(iu17)%im2 = im10
      itab_u(iu17)%iw1 = iw12
      itab_u(iu17)%iw2 = iw16

      itab_u(iu17)%iup = iu15-iucyc
      call uloops('f',iu17, 1, 5, 7, 8,11,14,15,18, 0, 0)
   endif

! IW5

   itab_w(iw5)%iu1 = iu1
   itab_w(iw5)%iu2 = iu5
   itab_w(iw5)%iu3 = iu8
   itab_w(iw5)%mrlw = 1
   itab_w(iw5)%mrlw_orig = 1

   if (i == 1) then
      itab_w(iw5)%iwp = iw7+iwcyc
      call wloops('f',iw5, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw5,24,31,33,35, 0, 0, 0, 0, 0, 0)
   elseif (i == nxp) then
      itab_w(iw5)%iwp = iw7-iwcyc
      call wloops('f',iw5, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw5,24,31,33,35, 0, 0, 0, 0, 0, 0)
   else
      itab_w(iw5)%iwp = iw7
      call wloops('f',iw5, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw5,24,31,33,35, 0, 0, 0, 0, 0, 0)
   endif

! IW6

   itab_w(iw6)%iu1 = iu2
   itab_w(iw6)%iu2 = iu9
   itab_w(iw6)%iu3 = iu6
   itab_w(iw6)%mrlw = 1
   itab_w(iw6)%mrlw_orig = 1

   if (i == 1) then
      itab_w(iw6)%iwp = iw6+iwcyc
      call wloops('f',iw6, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw6,24,31,33,35, 0, 0, 0, 0, 0, 0)
   elseif (i == nxp) then
      itab_w(iw6)%iwp = iw6-iwcyc
      call wloops('f',iw6, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw6,24,31,33,35, 0, 0, 0, 0, 0, 0)
   else
      itab_w(iw6)%iwp = iw6
      call wloops('f',iw6, 1, 3, 5, 6, 7, 8,11,12,13,14)
      call wloops('n',iw6,16,17,18,19,25,26,27,28,29,30)
      call wloops('n',iw6,33,34, 0, 0, 0, 0 ,0 ,0 ,0 ,0)
   endif

! IW7

   itab_w(iw7)%iu1 = iu3
   itab_w(iw7)%iu2 = iu6
   itab_w(iw7)%iu3 = iu10
   itab_w(iw7)%mrlw = 1
   itab_w(iw7)%mrlw_orig = 1

   if (i == 1) then
      itab_w(iw7)%iwp = iw7+iwcyc
      call wloops('f',iw7, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw7,24,31,33,35, 0, 0, 0, 0, 0, 0)
   elseif (i == nxp) then
      itab_w(iw7)%iwp = iw7-iwcyc
      call wloops('f',iw7, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw7,24,31,33,35, 0, 0, 0, 0, 0, 0)
   else
      itab_w(iw7)%iwp = iw7
      call wloops('f',iw7, 1, 3, 5, 6, 7, 8,11,12,13,14)
      call wloops('n',iw7,16,17,18,19,25,26,27,28,29,30)
      call wloops('n',iw7,33,34, 0, 0, 0, 0, 0, 0, 0, 0)
   endif

! IW8

   itab_w(iw8)%iu1 = iu4
   itab_w(iw8)%iu2 = iu11
   itab_w(iw8)%iu3 = iu7
   itab_w(iw8)%mrlw = 1
   itab_w(iw8)%mrlw_orig = 1

   if (i == 1) then
      itab_w(iw8)%iwp = iw6+iwcyc
      call wloops('f',iw8, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw8,24,31,33,35, 0, 0, 0, 0, 0, 0)
   elseif (i == nxp) then
      itab_w(iw8)%iwp = iw6-iwcyc
      call wloops('f',iw8, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw8,24,31,33,35, 0, 0, 0, 0, 0, 0)
   else
      itab_w(iw8)%iwp = iw6
      call wloops('f',iw8, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw8,24,31,33,35, 0, 0, 0, 0, 0, 0)
   endif

! IW9

   itab_w(iw9)%iu1 = iu8
   itab_w(iw9)%iu2 = iu14
   itab_w(iw9)%iu3 = iu12
   itab_w(iw9)%mrlw = 1
   itab_w(iw9)%mrlw_orig = 1

   if (i == 1) then
      itab_w(iw9)%iwp = iw11+iwcyc
      call wloops('f',iw9, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw9,24,31,33,35, 0, 0, 0, 0, 0, 0)
   elseif (i == nxp) then
      itab_w(iw9)%iwp = iw11-iwcyc
      call wloops('f',iw9, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw9,24,31,33,35, 0, 0, 0, 0, 0, 0)
   else
      itab_w(iw9)%iwp = iw11
      call wloops('f',iw9, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw9,24,31,33,35, 0, 0, 0, 0, 0, 0)
   endif

! IW10

   itab_w(iw10)%iu1 = iu9
   itab_w(iw10)%iu2 = iu12
   itab_w(iw10)%iu3 = iu15
   itab_w(iw10)%mrlw = 1
   itab_w(iw10)%mrlw_orig = 1

   if (i == 1) then
      itab_w(iw10)%iwp = iw10+iwcyc
      call wloops('f',iw10, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw10,24,31,33,35, 0, 0, 0, 0, 0, 0)
   elseif (i == nxp) then
      itab_w(iw10)%iwp = iw10-iwcyc
      call wloops('f',iw10, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw10,24,31,33,35, 0, 0, 0, 0, 0, 0)
   else
      itab_w(iw10)%iwp = iw10
      call wloops('f',iw10, 1, 3, 5, 6, 7, 8,11,12,13,14)
      call wloops('n',iw10,16,17,18,19,25,26,27,28,29,30)
      call wloops('n',iw10,33,34, 0, 0, 0, 0, 0, 0, 0, 0)
   endif

! IW11

   itab_w(iw11)%iu1 = iu10
   itab_w(iw11)%iu2 = iu16
   itab_w(iw11)%iu3 = iu13
   itab_w(iw11)%mrlw = 1
   itab_w(iw11)%mrlw_orig = 1

   if (i == 1) then
      itab_w(iw11)%iwp = iw11+iwcyc
      call wloops('f',iw11, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw11,24,31,33,35, 0, 0, 0, 0, 0, 0)
   elseif (i == nxp) then
      itab_w(iw11)%iwp = iw11-iwcyc
      call wloops('f',iw11, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw11,24,31,33,35, 0, 0, 0, 0, 0, 0)
   else
      itab_w(iw11)%iwp = iw11
      call wloops('f',iw11, 1, 3, 5, 6, 7, 8,11,12,13,14)
      call wloops('n',iw11,16,17,18,19,25,26,27,28,29,30)
      call wloops('n',iw11,33,34, 0, 0, 0, 0, 0, 0, 0, 0)
   endif

! IW12

   itab_w(iw12)%iu1 = iu11
   itab_w(iw12)%iu2 = iu13
   itab_w(iw12)%iu3 = iu17
   itab_w(iw12)%mrlw = 1
   itab_w(iw12)%mrlw_orig = 1

   if (i == 1) then
      itab_w(iw12)%iwp = iw10+iwcyc
      call wloops('f',iw12, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw12,24,31,33,35, 0, 0, 0, 0, 0, 0)
   elseif (i == nxp) then
      itab_w(iw12)%iwp = iw10-iwcyc
      call wloops('f',iw12, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw12,24,31,33,35, 0, 0, 0, 0, 0, 0)
   else
      itab_w(iw12)%iwp = iw10
      call wloops('f',iw12, 1, 3, 5, 6, 7, 8,11,14,18,22)
      call wloops('n',iw12,24,31,33,35, 0, 0, 0, 0, 0, 0)
   endif

   unit_dist = .5 * sqrt(3.) * deltax ! This is 1/2 of triangle face width
   diamond_centx = (2 * i - 1 - nxp) * unit_dist

   xem(im1)  = diamond_centx - 1.5 * unit_dist
   xem(im2)  = xem(im1)
   xem(im3)  = xem(im1)
   xem(im4)  = diamond_centx -  .5 * unit_dist
   xem(im5)  = xem(im4)
   xem(im6)  = diamond_centx +  .5 * unit_dist
   xem(im7)  = xem(im6)
   xem(im8)  = xem(im6)
   xem(im9)  = diamond_centx + 1.5 * unit_dist
   xem(im10) = xem(im9)

   yem(im1)  = - 3.0 * deltax
   yem(im2)  =    .0 * deltax
   yem(im3)  =   3.0 * deltax
   yem(im4)  = - 1.5 * deltax
   yem(im5)  =   1.5 * deltax
   yem(im6)  = - 3.0 * deltax
   yem(im7)  =    .0 * deltax
   yem(im8)  =   3.0 * deltax
   yem(im9)  = - 1.5 * deltax
   yem(im10) =   1.5 * deltax

enddo  ! end i loop

zem(2:nma) = 0.

call tri_neighbors()

return
end subroutine cartesian_3d

