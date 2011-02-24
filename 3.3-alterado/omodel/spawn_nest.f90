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
subroutine spawn_nest()

! This subroutine adds nested grid regions at the beginning of a simulation.
! Later will make modified version to add nested grid region(s) during a
!   simulation.

use mem_ijtabs,  only: itab_m, itab_u, itab_w, ltab_m, ltab_u, ltab_w,      &
                       nest_u, nest_w, nloops_m, nloops_u, nloops_w, mrls,  &
                       alloc_itabs
use mem_grid,    only: nma, nua, nwa, mma, mua, mwa, xem, yem, zem,  &
                       alloc_xyzem, nrows, mrows
use misc_coms,   only: io6, ngrids, mdomain, grdlen, grdwid, grdaxis,  &
                       centlat, centlon, nxp
use consts_coms, only: pio180, erad, pi1, pi2

implicit none

integer :: iu,iw,iregion,im,iw1,iw2,im1,im2,im3 &
   ,iu1,iu2,iu3,ndiv,iu1o,iu2o,iu3o,iu1o_iw1,iu2o_iw1,iu3o_iw1  &
   ,iu4,iu5,iu6,iw3,ngr,mrlo,mrloo,mwa0,itpn,ntpn,nw,inudp

real :: focdist,xf1,xf2,yf1,yf2,widfac,xw,yw,xm1,ym1,xm2,ym2,xm3,ym3
real :: expansion

integer :: nside    ! Number of sides of nested grid polygon
integer :: nper(16) ! Value of perimeter side counter at end of each side
integer, allocatable :: imper(:) ! Ouside IW index at each perimeter index
integer, allocatable :: iuper(:) ! Boundary IU index at each perimeter index
integer :: kma,kua,kwa   ! New M,U,W indices while constructing nest perimeter
integer :: ngrp   ! Number of perimeter groups
integer :: npts   ! Estimated number of CM pts along FM perimeter

integer, allocatable :: jm(:,:),ju(:,:),iurow_pent(:),igsize(:)

real, allocatable :: xem_temp(:),yem_temp(:),zem_temp(:)

! Set number of transition rows to be used; may put this in OLAMIN eventually

mrows = 3

do ngr = 2,ngrids  ! Loop over nested grids

! Estimate number of CM points along FM perimeter

   npts = nxp * 2**(ngr-2) * grdlen(ngr) / 1.e6 + 10

! Allocate temporary tables

   allocate (ltab_m(nma))  ! Duplicates for itab_m
   allocate (ltab_u(nua))  ! Duplicates for itab_u
   allocate (ltab_w(nwa))  ! Duplicates for itab_w

   allocate (nest_u(nua), nest_w(nwa))  ! Nest relations
   allocate (xem_temp(nma), yem_temp(nma), zem_temp(nma))
   allocate (iurow_pent(nua))

   allocate (jm(nrows+1,npts), ju(nrows,npts))
   allocate (imper(npts), iuper(npts), igsize(npts))

   call pent_urows(iurow_pent)

! Copy ITAB information and M-point coordinates to temporary tables

   ltab_m(1:nma) = itab_m(1:nma)
   ltab_u(1:nua) = itab_u(1:nua)
   ltab_w(1:nwa) = itab_w(1:nwa)

   do im = 1,nma
      xem_temp(im) = xem(im)
      yem_temp(im) = yem(im)
      zem_temp(im) = zem(im)
   enddo

! Deallocate main tables

   deallocate (itab_m, itab_u, itab_w)
   deallocate (xem, yem, zem)      

! Distance of foci from grid center; Local PS coordinates of foci

   focdist = .5 * sqrt(abs(grdlen(ngr)**2 - grdwid(ngr)**2))
   xf2 =  focdist * cos((90. - grdaxis(ngr)) * pio180)   
   yf2 =  focdist * sin((90. - grdaxis(ngr)) * pio180)   
   xf1 = -xf2   
   yf1 = -yf2   

! Locate and flag all W triangles to be subdivided
! Define 3 new W indices and 3 new U indices for each W - attach to W.

   do iw = 2,nwa  ! Loop over all previously-existing W points

      iu1 = ltab_w(iw)%iu1
      iu2 = ltab_w(iw)%iu2
      
      im1 = ltab_u(iu1)%im1
      im2 = ltab_u(iu1)%im2
      im3 = ltab_u(iu2)%im1
      if (im3 == im1 .or. im3 == im2) im3 = ltab_u(iu2)%im2

! If using spherical geometry, transform IM1, IM2, IM3 locations to 
! local PS space wrt grid center

      if (mdomain < 2) then

         call e_ps(xem_temp(im1),yem_temp(im1),zem_temp(im1)  &
            ,centlat(ngr),centlon(ngr),xm1,ym1)
         call e_ps(xem_temp(im2),yem_temp(im2),zem_temp(im2)  &
            ,centlat(ngr),centlon(ngr),xm2,ym2)
         call e_ps(xem_temp(im3),yem_temp(im3),zem_temp(im3)  &
            ,centlat(ngr),centlon(ngr),xm3,ym3)

       else

! If NOT using spherical geometry, copy IM1, IM2, IM3 horizontal locations 
! to local cartesian space wrt grid center

         xm1 = xem_temp(im1) - centlon(ngr)
         xm2 = xem_temp(im2) - centlon(ngr)
         xm3 = xem_temp(im3) - centlon(ngr)

         ym1 = yem_temp(im1) - centlat(ngr)
         ym2 = yem_temp(im2) - centlat(ngr)
         ym3 = yem_temp(im3) - centlat(ngr)
         
      endif

! Check whether W location is within elliptical current nested region

      if (sqrt((xm1-xf1)**2 + (ym1-yf1)**2) +                    &
          sqrt((xm1-xf2)**2 + (ym1-yf2)**2) < grdlen(ngr) .and.  &
          sqrt((xm2-xf1)**2 + (ym2-yf1)**2) +                    &
          sqrt((xm2-xf2)**2 + (ym2-yf2)**2) < grdlen(ngr) .and.  &
          sqrt((xm3-xf1)**2 + (ym3-yf1)**2) +                    &
          sqrt((xm3-xf2)**2 + (ym3-yf2)**2) < grdlen(ngr)) then

! W cell is to be subdivided

!        write(io6,*) 'subdividing W cell ',iw

         nest_w(iw)%iu1 = mua + 1
         nest_w(iw)%iu2 = mua + 2
         nest_w(iw)%iu3 = mua + 3

         nest_w(iw)%iw1 = mwa + 1
         nest_w(iw)%iw2 = mwa + 2
         nest_w(iw)%iw3 = mwa + 3

         mua = mua + 3
         mwa = mwa + 3

      endif
      
   enddo

! Add W points to nested grid region in order to eliminate concavities.
! This requires iterative procedure

   mwa0 = 0  ! Counter of already-existing W points - initialize to zero to
             ! force at least one pass through the following DO WHILE loop

   do while (mwa > mwa0)

      mwa0 = mwa

      do im = 2,nma

         ntpn = ltab_m(im)%ntpn
         nw = 0  ! Initialize counter for subdivided W points around this M point


!   write(io6,*) 'expanding ',im,ntpn

         do itpn = 1,ntpn
            iw = ltab_m(im)%iw(itpn)
            
!     write(io6,*) 'counting ',itpn,iw
            
            if (nest_w(iw)%iw1 > 0) then
               nw = nw + 1  ! Count up subdivided W points around this M point
            endif
         enddo

! Check value of nw for illegal values

         if (ntpn > 4 .and. (ntpn - nw == 1 .or. ntpn - nw == 2)) then

! This M point is a concavity so activate remaining unactivated W points around it

            do itpn = 1,ntpn
               iw = ltab_m(im)%iw(itpn)

               if (nest_w(iw)%iw1 == 0) then
                  nest_w(iw)%iu1 = mua + 1
                  nest_w(iw)%iu2 = mua + 2
                  nest_w(iw)%iu3 = mua + 3

                  nest_w(iw)%iw1 = mwa + 1
                  nest_w(iw)%iw2 = mwa + 2
                  nest_w(iw)%iw3 = mwa + 3

                  mua = mua + 3
                  mwa = mwa + 3
                  
!               write(io6,*) 'Activiting W point ',iw,' to prevent concavity'  
               
               endif
            enddo

         endif

      enddo
   enddo

! Nested region should be fully expanded without concavities now

! Define new vertex index for midpoint of each original U edge that is adjacent
! to an original triangle that is being subdivided.  Attach new vertex
! index to old U edge.  Also, define new U index for second half of U.
! Attach to U.  [Make adjacent to U%m2.]

   do iu = 2,nua

      ! Check whether this U is adjacent to a W that is being subdivided

      iw1 = ltab_u(iu)%iw1
      iw2 = ltab_u(iu)%iw2

      if (nest_w(iw1)%iw3 > 0 .or. nest_w(iw2)%iw3 > 0) then

         nest_u(iu)%im = mma + 1
         nest_u(iu)%iu = mua + 1

         mma = mma + 1
         mua = mua + 1

      endif
   enddo

! Map out perimeter of new mesh refined region

   call perim_map(nper,nside,npts,imper,iuper)

! Save current values of mma, mua, mwa prior to adding boundary points

   kma = mma
   kua = mua
   kwa = mwa

! Form groups of original M and U indices on perimeter and increase mma, mua, 
! and mwa according to what groups will require

   call perim_add(jm,ju,npts,nside  &
                    ,nper,ngrp,igsize,imper,iuper,iurow_pent)

! Allocate main tables to expanded size
! Initialize all neighbor indices to zero

   call alloc_itabs(mma,mua,mwa)

   call alloc_xyzem()

! Memory copy to main tables

   do im = 1,nma
      itab_m(im)%loop(1:nloops_m) = ltab_m(im)%loop(1:nloops_m)
      itab_m(im)%itopm = ltab_m(im)%itopm
      xem(im) = xem_temp(im)
      yem(im) = yem_temp(im)
      zem(im) = zem_temp(im)
   enddo

   itab_u(1:nua) = ltab_u(1:nua)
   itab_w(1:nwa) = ltab_w(1:nwa)
  
! Average coordinates to new M points

   do iu = 2,nua
      if (nest_u(iu)%im > 0) then
      	 im = nest_u(iu)%im
      	 im1 = itab_u(iu)%im1
      	 im2 = itab_u(iu)%im2

      	 xem(im) = .5 * (xem(im1) + xem(im2))
         yem(im) = .5 * (yem(im1) + yem(im2))
         zem(im) = .5 * (zem(im1) + zem(im2))
      endif
   enddo

! Contruct tables for new fully subdivided triangles

   mrloo = 0  ! Initialize check variable for uniform mrlw over current nested grid

   do iw = 2,nwa

! Check if IW is fully subdivided cell

      if (nest_w(iw)%iw3 > 0) then

! This is fully subdivided W cell

! Mapping of original undivided triangles

         iu1o = ltab_w(iw)%iu1
         iu2o = ltab_w(iw)%iu2
         iu3o = ltab_w(iw)%iu3
         mrlo = ltab_w(iw)%mrlw
         inudp = ltab_w(iw)%inudp(1)

! Check of mrlw value for current nested grid

         if (mrloo == 0) mrloo = mrlo  ! Set to first nonzero mrlw encountered
         if (mrlo /= mrloo) then
            write(io6,*) 'Current nested grid ',ngr
            write(io6,*) 'crosses pre-existing grid boundary.'
            write(io6,*) 'iw = ',iw
            write(io6,*) 'stopping model'
            stop 'stop - nested grid out of bounds'
         endif

         iu1o_iw1 = ltab_u(iu1o)%iw1
         iu2o_iw1 = ltab_u(iu2o)%iw1
         iu3o_iw1 = ltab_u(iu3o)%iw1

         ! Mapping of new divided triangles

         iu1 = nest_w(iw)%iu1        
         iu2 = nest_w(iw)%iu2        
         iu3 = nest_w(iw)%iu3        

         itab_u(iu1)%iup = iu1
         itab_u(iu2)%iup = iu2
         itab_u(iu3)%iup = iu3
         
         call uloops('f',iu1, 1, 4, 5, 7, 8,11,12,13,14,15)
         call uloops('n',iu1,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

         call uloops('f',iu2, 1, 4, 5, 7, 8,11,12,13,14,15)
         call uloops('n',iu2,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

         call uloops('f',iu3, 1, 4, 5, 7, 8,11,12,13,14,15)
         call uloops('n',iu3,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

         iu4 = nest_u(iu1o)%iu       
         iu5 = nest_u(iu2o)%iu       
         iu6 = nest_u(iu3o)%iu

         itab_u(iu4)%loop(1:nloops_u) = itab_u(iu1o)%loop(1:nloops_u)
         itab_u(iu5)%loop(1:nloops_u) = itab_u(iu2o)%loop(1:nloops_u)
         itab_u(iu6)%loop(1:nloops_u) = itab_u(iu3o)%loop(1:nloops_u)

         itab_u(iu4)%iup = iu4
         itab_u(iu5)%iup = iu5
         itab_u(iu6)%iup = iu6

         iw1 = nest_w(iw)%iw1        
         iw2 = nest_w(iw)%iw2        
         iw3 = nest_w(iw)%iw3 

         ! Fill tables with new values

         itab_w(iw)%iu1 = iu1       
         itab_w(iw)%iu2 = iu2       
         itab_w(iw)%iu3 = iu3   

         itab_w(iw1)%iu1 = iu1    
         itab_w(iw2)%iu1 = iu2    
         itab_w(iw3)%iu1 = iu3

         itab_w(iw)%mrlw  = mrlo + 1
         itab_w(iw1)%mrlw = mrlo + 1
         itab_w(iw2)%mrlw = mrlo + 1
         itab_w(iw3)%mrlw = mrlo + 1

         itab_w(iw1)%mrlw_orig = mrlo + 1
         itab_w(iw2)%mrlw_orig = mrlo + 1
         itab_w(iw3)%mrlw_orig = mrlo + 1

         itab_w(iw)%inudp(1)  = inudp
         itab_w(iw1)%inudp(1) = inudp
         itab_w(iw2)%inudp(1) = inudp
         itab_w(iw3)%inudp(1) = inudp

         itab_u(iu1o)%im2 = nest_u(iu1o)%im
         itab_u(iu2o)%im2 = nest_u(iu2o)%im
         itab_u(iu3o)%im2 = nest_u(iu3o)%im

         itab_u(iu4)%im1 = nest_u(iu1o)%im
         itab_u(iu4)%im2 = ltab_u(iu1o)%im2
         itab_u(iu5)%im1 = nest_u(iu2o)%im
         itab_u(iu5)%im2 = ltab_u(iu2o)%im2
         itab_u(iu6)%im1 = nest_u(iu3o)%im
         itab_u(iu6)%im2 = ltab_u(iu3o)%im2

         if (iw == iu1o_iw1) then
            itab_w(iw3)%iu2 = iu1o
            itab_w(iw2)%iu3 = iu4

            itab_u(iu1)%im1 = nest_u(iu2o)%im            
            itab_u(iu1)%im2 = nest_u(iu3o)%im            
            itab_u(iu1)%iw1 = iw1
            itab_u(iu1)%iw2 = iw

            itab_u(iu1o)%iw1 = iw3
            itab_u(iu4)%iw1 = iw2
         else
            itab_w(iw3)%iu2 = iu4
            itab_w(iw2)%iu3 = iu1o

            itab_u(iu1)%im1 = nest_u(iu3o)%im            
            itab_u(iu1)%im2 = nest_u(iu2o)%im            
            itab_u(iu1)%iw1 = iw
            itab_u(iu1)%iw2 = iw1

            itab_u(iu1o)%iw2 = iw2
            itab_u(iu4)%iw2 = iw3
         endif    

         if (iw == iu2o_iw1) then
            itab_w(iw1)%iu2 = iu2o
            itab_w(iw3)%iu3 = iu5

            itab_u(iu2)%im1 = nest_u(iu3o)%im            
            itab_u(iu2)%im2 = nest_u(iu1o)%im            
            itab_u(iu2)%iw1 = iw2
            itab_u(iu2)%iw2 = iw

            itab_u(iu2o)%iw1 = iw1
            itab_u(iu5)%iw1 = iw3
         else
            itab_w(iw1)%iu2 = iu5
            itab_w(iw3)%iu3 = iu2o

            itab_u(iu2)%im1 = nest_u(iu1o)%im            
            itab_u(iu2)%im2 = nest_u(iu3o)%im            
            itab_u(iu2)%iw1 = iw
            itab_u(iu2)%iw2 = iw2
            
            itab_u(iu2o)%iw2 = iw3
            itab_u(iu5)%iw2 = iw1
         endif    

         if (iw == iu3o_iw1) then
            itab_w(iw2)%iu2 = iu3o
            itab_w(iw1)%iu3 = iu6

            itab_u(iu3)%im1 = nest_u(iu1o)%im            
            itab_u(iu3)%im2 = nest_u(iu2o)%im            
            itab_u(iu3)%iw1 = iw3
            itab_u(iu3)%iw2 = iw

            itab_u(iu3o)%iw1 = iw2
            itab_u(iu6)%iw1 = iw1
         else
            itab_w(iw2)%iu2 = iu6
            itab_w(iw1)%iu3 = iu3o

            itab_u(iu3)%im1 = nest_u(iu2o)%im            
            itab_u(iu3)%im2 = nest_u(iu1o)%im            
            itab_u(iu3)%iw1 = iw
            itab_u(iu3)%iw2 = iw3

            itab_u(iu3o)%iw2 = iw1
            itab_u(iu6)%iw2 = iw2
         endif    

      endif

   enddo    ! end of iw loop

! Fill transition zone

   call perim_fill(ngr,mrloo,kma,kua,kwa,jm,ju,npts,ngrp,igsize)

! Fill itabs loop tables for newly spawned points (U pts already done above)

   do im = nma+1,mma
      itab_m(im)%itopm = im
      call mloops('f',im,1,0,0,0)
   enddo

   do iw = nwa+1,mwa
      itab_w(iw)%iwp = iw
      call wloops('f',iw, 1, 3, 5, 6, 7, 8,11,12,13,14)
      call wloops('n',iw,15,16,17,18,19,20,25,26,27,28)
      call wloops('n',iw,29,30,33,34, 0, 0, 0 ,0 ,0 ,0)
   enddo

   nma = mma    ! Copy new counter value to nma
   nua = mua    ! Copy new counter value to nua
   nwa = mwa    ! Copy new counter value to nwa

   deallocate (ltab_m,ltab_u,ltab_w)
   deallocate (nest_u,nest_w)
   deallocate (xem_temp,yem_temp,zem_temp)

   deallocate (jm,ju)
   deallocate (imper,iuper,iurow_pent,igsize)

   call tri_neighbors()

! Call subroutine to ID W cells just outside and just inside current NGR
! border.  This is permanent ID, used in spring dynamics even when new
! grids are added.

   call perim_mrow()

! This is the place to do spring dynamics

   call spring_dynamics()

   write(io6,'(/,a,i2)') 'Finished spawning grid number ',ngr
   write(io6,'(a,i8)')   ' mma = ',mma
   write(io6,'(a,i8)')   ' mua = ',mua
   write(io6,'(a,i8)')   ' mwa = ',mwa

enddo   ! end of ngr loop

! Set mrls equal to maximum mrlw value

do iw = 2,nwa
   if (mrls < itab_w(iw)%mrlw) mrls = itab_w(iw)%mrlw
enddo

! Push all M point coordinates out to earth radius if not already there

do im = 2,nma
   expansion = erad / sqrt(xem(im) ** 2  &
                         + yem(im) ** 2  &
                         + zem(im) ** 2  )

   xem(im) = xem(im) * expansion
   yem(im) = yem(im) * expansion
   zem(im) = zem(im) * expansion
enddo

return
end subroutine spawn_nest

!==============================================================================

subroutine pent_urows(iurow_pent)

use mem_ijtabs, only: itab_m,itab_u,itab_w
use mem_grid,   only: nua, nwa, impent, mrows
use misc_coms,  only: io6

implicit none

integer, intent(out) :: iurow_pent(nua)

integer :: ipent,im,itpn,iw,irow,jrow,iw1,iw2,iw3,iwrow,iu

! automatic arrays

integer :: iwrow_temp1(nwa)
integer :: iwrow_temp2(nwa)

! Initialize temporary arrays to zero before loop over pentagon points

iurow_pent(1:nua) = 0
iwrow_temp1(1:nwa) = 0
iwrow_temp2(1:nwa) = 0

! Loop over all 12 pentagon M points

do ipent = 1,12

   im = impent(ipent)
! Loop over W points that surround current M point; set row flag to 1

   do itpn = 1,5
      iw = itab_m(im)%iw(itpn)
      iwrow_temp1(iw) = 1
      iwrow_temp2(iw) = 1
   enddo
   
enddo
   
! Advance outward and flag each row

do irow = 1,2*mrows-1
   jrow = mod(irow,2)

   do iw = 2,nwa

      if (iwrow_temp1(iw) == 0) then

! If IW is adjacent to any other IW cell with nonzero mrow, 
! set mrow for IW cell 

         iw1 = itab_w(iw)%iw1
         iw2 = itab_w(iw)%iw2
         iw3 = itab_w(iw)%iw3

! Check for positive mrow values

         iwrow = max(iwrow_temp1(iw1)  &
                    ,iwrow_temp1(iw2)  &
                    ,iwrow_temp1(iw3))

         if (iwrow > 0) iwrow_temp2(iw) = iwrow + jrow

      endif

   enddo

   do iw = 2,nwa
      iwrow_temp1(iw) = iwrow_temp2(iw)
   enddo

enddo

! Loop over all U points and flag those between unequal nonzero iwrow values

do iu = 2,nua
   iw1 = itab_u(iu)%iw1            
   iw2 = itab_u(iu)%iw2            

   if (iwrow_temp1(iw1) > 0 .and. iwrow_temp1(iw2) > 0 .and.  &
       iwrow_temp1(iw1) /= iwrow_temp1(iw2)) then
          
       iurow_pent(iu) = min(iwrow_temp1(iw1),iwrow_temp1(iw2))
   endif
enddo

return
end subroutine pent_urows

!==============================================================================

subroutine perim_map(nper,nside,npts,imper,iuper)

! Perim_map maps the perimeter points of a nested grid that is being spawned

use mem_grid,  only: nma
use misc_coms, only: io6

implicit none

integer, intent(out) :: nper(16) ! Value of perimeter side counter at end of each side
integer, intent(out) :: nside    ! Number of sides of nested grid polygon
integer, intent(in)  :: npts
integer, intent(out) :: imper(npts) ! Boundary IM index at each perimeter index
integer, intent(out) :: iuper(npts) ! Boundary IU index at each perimeter index

integer :: imstart  ! IM index of starting M point (at a corner of ngr perimeter)
integer :: im       ! dummy im index
integer :: nwdiv    ! number of subdivided W pts adjacent to current M pt
integer :: ima      ! Current M pt in counterclockwise path around perimeter
integer :: imb      ! Next M pt in counterclockwise path around perimeter
integer :: iua      ! Current U pt in counterclockwise path around perimeter
integer :: iper     ! Current perimeter point counter
integer :: iside    ! Current perimeter side counter

! Set IM starting point to zero

imstart = 0

! Loop over all ORIGINAL M points and find the first one that has exactly 2 
! fully-divided W neighbors in the grid being spawned.

do im = 2,nma   
   call count_wdiv(im,nwdiv)
   if (nwdiv == 2) then
      imstart = im
      exit
   endif
enddo

! Bug check:  make sure that imstart is not zero now.

if (imstart == 0) then
   write(io6,*) 'imstart is zero - stopping model'
   stop 'stop imstart'
endif

! March around NGR boundary in a counterclockwise direction beginning at IMSTART

ima = imstart
imper(1) = ima
iper = 0
iside = 1

nper(:) = 0

do  ! Loop over all original M points on ngr perimeter 

! Find next M and U points (imb, iua) on perimeter and outside W point (iwout)

   call perim_ngr(ima, imb, iua)

! Increment perimeter point counter and store point indices

   nper(iside) = nper(iside) + 1
   iper = iper + 1

   iuper(iper) = iua
   imper(iper+1) = imb

! Check if imb is a corner point.

   call count_wdiv(imb,nwdiv)

   if (nwdiv == 2) then

! imb is corner point.  Store current iper index in nper array.

      nper(iside) = iper

! Check if imb equals istart.  If it does, exit loop

      if (imb == imstart) then

! imb equals istart.  Store total number of sides in nside.  Exit Do loop
      
         nside = iside
         exit
      endif
   
! imb does not equal istart.  Advance to next side

      iside = iside + 1
         
   endif

   ima = imb

enddo

return
end subroutine perim_map

!==============================================================================

subroutine count_wdiv(im,nwdiv)

! Subroutine count_wdiv is to be used during the process of spawning a nested grid,
! after temporary arrays nest_u and nest_w have been filled for interior points.

! Given any M point IM, determine how many fully-divided W points are adjacent to it.

use mem_ijtabs, only: ltab_m, nest_w
use misc_coms,  only: io6

implicit none

integer, intent(in) :: im
integer, intent(out) :: nwdiv

integer :: itpn, iw

nwdiv = 0

! Loop over all W points that are adjacent to this M point

do itpn = 1,ltab_m(im)%ntpn
   iw = ltab_m(im)%iw(itpn)

   if (nest_w(iw)%iw3 > 0) then
      nwdiv = nwdiv + 1
   endif
enddo

return
end subroutine count_wdiv

!==============================================================================

subroutine perim_ngr(imstart, imnext, iunext)

! Subroutine perim_ngr is to be used during the process of spawning a nested grid,
! after temporary arrays nest_u and nest_w have been filled for interior points.

! Given any M point on the nested grid boundary that existed before the spawn 
! process began, and proceeding along the boundary in a counterclockwise direction,
! find the adjacent U and M points that existed before the spawn process began.

use mem_ijtabs, only: ltab_m, ltab_u, ltab_w, nest_u, nest_w
use misc_coms,  only: io6

implicit none

integer, intent(in) :: imstart  ! starting original M point

integer, intent(out) :: imnext  ! next original M pt (counterclockwise)
integer, intent(out) :: iunext  ! next original U pt (counterclockwise)

integer :: itpn, im1, im2, iw1, iw2, iu

imnext = 0
iunext = 0

! Loop over all U points connected to current M point

do itpn = 1,ltab_m(imstart)%ntpn
   iu = ltab_m(imstart)%iu(itpn)
   
   im1 = ltab_u(iu)%im1
   im2 = ltab_u(iu)%im2

   iw1 = ltab_u(iu)%iw1
   iw2 = ltab_u(iu)%iw2

! Check if current IU point is along NGR boundary

   if (nest_u(iu)%im > 0) then

! Current IU point is either along or inside NGR boundary

      if (im1 == imstart .and. nest_w(iw1)%iw3 == 0) then

! Current IU point is on boundary and in clockwise direction
! and next M point is im2.          

         iunext = iu
         imnext = im2

         exit

      elseif (im2 == imstart .and. nest_w(iw2)%iw3 == 0) then

! Current IU point is on boundary and in clockwise direction
! and next M point is im1.          

         iunext = iu
         imnext = im1

         exit

      endif

   endif

enddo

if (iunext < 2) then
   write(io6,*) 'iunext < 2 in subroutine perim_ngr - stopping model'
   stop 'perim_ngr1'
elseif (imnext < 2) then
   write(io6,*) 'imnext < 2 in subroutine perim_ngr - stopping model'
   stop 'perim_ngr2'
endif

return
end subroutine perim_ngr

!==============================================================================

subroutine perim_add(jm,ju,npts,nside  &
                    ,nper,ngrp,igsize,imper,iuper,iurow_pent)

use mem_grid,   only: nrows, mrows, nua, mma, mua, mwa
use mem_ijtabs, only: itab_m,itab_u,itab_w
use misc_coms,  only: io6

implicit none

integer, intent(inout) :: npts,nside
integer, intent(out) :: jm(nrows+1,npts)
integer, intent(out) :: ju(nrows,npts)
integer, intent(out) :: ngrp
integer, intent(out) :: igsize(npts)
integer, intent(in)  :: imper(npts)
integer, intent(in)  :: iuper(npts)
integer, intent(inout)  :: nper(16)
integer, intent(in) :: iurow_pent(nua)

integer :: iper,ig,iside,ileft,iwid
integer :: igs   ! counter for group number on current side

integer :: iu,jside

! Loop over all U points on perimeter of new grid

iper = 1
igs = 0  ! group size
iside = 1

do while (iside <= nside)

   do while (iper <= nper(iside))

      iu = iuper(iper)

! Check if current U point is influenced by pentagon M point.
! If so, increase group size.

      if (iurow_pent(iu) > 0) then
         igs = igs + 1
         
! Check if group size equals iurow_pent(iu)

         if (igs == iurow_pent(iu)) then         
         
! This group of points is influenced by pentagon M point.  Place them
! Into a new "side" of the nest.  Check if there are any remaining points
! on current ORIGINAL side.  Shift remaining iper points to new sides
! (a shift of 1 or 2 sides).
            
            if (iper == nper(iside)) then
               do jside = nside,iside,-1
                  nper(jside+1) = nper(jside)
               enddo
               nside = nside + 1
            else
               do jside = nside,iside,-1
                  nper(jside+2) = nper(jside)
               enddo
               nside = nside + 2
               nper(iside + 1) = iper
            endif
            
            nper(iside) = iper - igs
            igs = 0

         endif 
         
      endif

      iper = iper + 1

   enddo
   
   iside = iside + 1
   
enddo

! Determine width of transition zone at each point on perimeter

iper = 1
ig = 0  ! group number

do iside = 1,nside

   igs = 0

   do while (iper <= nper(iside))

      ig = ig + 1
      igs = igs + 1

      ileft = nper(iside) - iper + 1

! FIRST METHOD: START WITH MAX AND DECREASE AS END OF SIDE IS APPROACHED

      if     (ileft == 16 .and. mrows == 5) then
         iwid = 4
      elseif (ileft == 12 .and. mrows == 5) then
         iwid = 4
      elseif (ileft == 11 .and. mrows == 5) then
         iwid = 4
      elseif (ileft == 9 .and. mrows == 4) then
         iwid = 3
      elseif (ileft == 8 .and. mrows == 5) then
         iwid = 4
      elseif (ileft == 7 .and. mrows == 5) then
         iwid = 4
      elseif (ileft == 6 .and. mrows >= 4) then
         iwid = 3
      elseif (ileft == 5 .and. mrows == 4) then
         iwid = 3
      elseif (ileft == 4 .and. mrows == 3) then
         iwid = 2
      else
         iwid = min(ileft,mrows)
      endif

! SECOND METHOD: USE 2 AT BOTH ENDS

!      if (igs == 1) then
!         iwid = min(ileft,2)
!      else
!         iwid = min(mrows,max(2,ileft-2),ileft)
!      endif

      igsize(ig) = iwid

                    jm(1,ig) = imper(iper)
                    jm(2,ig) = imper(iper+1)
      if (iwid > 1) jm(3,ig) = imper(iper+2)
      if (iwid > 2) jm(4,ig) = imper(iper+3)
      if (iwid > 3) jm(5,ig) = imper(iper+4)
      if (iwid > 4) jm(6,ig) = imper(iper+5)

                    ju(1,ig) = iuper(iper)
      if (iwid > 1) ju(2,ig) = iuper(iper+1)
      if (iwid > 2) ju(3,ig) = iuper(iper+2)       
      if (iwid > 3) ju(4,ig) = iuper(iper+3)      
      if (iwid > 4) ju(5,ig) = iuper(iper+4)

      iper = iper + iwid

! Add required number of W,U,M points for this group

      mwa = mwa + iwid**2
      mua = mua + iwid * (3 * iwid - 1) / 2
      mma = mma + iwid * (iwid - 1) / 2

   enddo

enddo

ngrp = ig

return
end subroutine perim_add

!===========================================================================

subroutine perim_fill(ngr,mrloo,kma,kua,kwa,jm,ju,npts,ngrp,igsize)

use mem_grid,   only: nrows,mrows
use mem_ijtabs, only: itab_u,itab_w
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa
integer, intent(in) :: npts
integer, intent(inout) :: jm(nrows+1,npts)
integer, intent(inout) :: ju(nrows,npts)
integer, intent(in) :: ngrp
integer, intent(inout) :: igsize(npts)

integer :: jmo(nrows+1),juo(nrows)  ! temp storage for output row

integer :: irow, ig, iu, iw  ! last 2 special

! Determine width of transition zone at each point on perimeter

do irow = 1,mrows
   do ig = 1,ngrp
   
      if (igsize(ig) == 1) then
      
         call perim_fill_cent1(ngr,mrloo,kma,kua,kwa      &
                              ,jm(1,ig),ju(1,ig),jm(2,ig) )
      
      elseif (igsize(ig) == 2) then
      
         call perim_fill_cent2(ngr,mrloo,kma,kua,kwa      &
                              ,jm(1,ig),ju(1,ig),jm(2,ig) &
                              ,ju(2,ig),jm(3,ig)          &
                              ,jmo(1),juo(1),jmo(2)       )
         
      elseif (igsize(ig) == 3) then
      
         call perim_fill_right(ngr,mrloo,kma,kua,kwa      &
                              ,jm(1,ig),ju(1,ig),jm(2,ig) &
                              ,jmo(1),juo(1),jmo(2)       )

         call perim_fill_cent1(ngr,mrloo,kma,kua,kwa      &
                              ,jm(2,ig),ju(2,ig),jm(3,ig) )

         call perim_fill_left (ngr,mrloo,kma,kua,kwa      &
                              ,jm(3,ig),ju(3,ig),jm(4,ig) &
                              ,jmo(2),juo(2),jmo(3)       )
      
      elseif (igsize(ig) == 4) then
      
         call perim_fill_right(ngr,mrloo,kma,kua,kwa      &
                              ,jm(1,ig),ju(1,ig),jm(2,ig) &
                              ,jmo(1),juo(1),jmo(2)       )

         call perim_fill_cent2(ngr,mrloo,kma,kua,kwa      &
                              ,jm(2,ig),ju(2,ig),jm(3,ig) &
                              ,ju(3,ig),jm(4,ig)          &
                              ,jmo(2),juo(2),jmo(3)       )

         call perim_fill_left (ngr,mrloo,kma,kua,kwa      &
                              ,jm(4,ig),ju(4,ig),jm(5,ig) &
                              ,jmo(3),juo(3),jmo(4)       )
      
      elseif (igsize(ig) == 5) then
      
         call perim_fill_right(ngr,mrloo,kma,kua,kwa      &
                              ,jm(1,ig),ju(1,ig),jm(2,ig) &
                              ,jmo(1),juo(1),jmo(2)       )

         call perim_fill_right(ngr,mrloo,kma,kua,kwa      &
                              ,jm(2,ig),ju(2,ig),jm(3,ig) &
                              ,jmo(2),juo(2),jmo(3)       )

         call perim_fill_cent1(ngr,mrloo,kma,kua,kwa      &
                              ,jm(3,ig),ju(3,ig),jm(4,ig) )

         call perim_fill_left (ngr,mrloo,kma,kua,kwa      &
                              ,jm(4,ig),ju(4,ig),jm(5,ig) &
                              ,jmo(3),juo(3),jmo(4)       )

         call perim_fill_left (ngr,mrloo,kma,kua,kwa      &
                              ,jm(5,ig),ju(5,ig),jm(6,ig) &
                              ,jmo(4),juo(4),jmo(5)       )

      endif
      
      igsize(ig) = igsize(ig) - 1
      
      jm(1:igsize(ig)+1,ig) = jmo(1:igsize(ig)+1)
      ju(1:igsize(ig)  ,ig) = juo(1:igsize(ig)  )
      
   enddo
      
enddo

return
end subroutine perim_fill

!===========================================================================

subroutine perim_fill_cent1(ngr,mrloo,kma,kua,kwa,jm1,ju1,jm2)

use mem_ijtabs, only: itab_w, itab_u, ltab_w, ltab_u, nest_w, nest_u
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa  ! Index values for latest added points

integer, intent(in) :: jm1,jm2  ! Original border M indices
integer, intent(in) :: ju1      ! Original border U index

integer :: ju4,ju7

integer :: iw1,iw2
integer :: iu1,iu2,iu3,iu4,iu5
integer :: im1,im2,im3,im4

! Assign new I indices for all M, W points, and for U points that need no checks

im1 = jm1
im2 = nest_u(ju1)%im
im3 = jm2

iu4 = kua + 1  ! Newly added point
iw2 = kwa + 1  ! Newly added point

! Increment indices for newly added points

kua = kua + 1
kwa = kwa + 1

! Determine orientation of positive ju1 direction

if (jm1 == ltab_u(ju1)%im1) then  ! Positive ju1 points inward
   iw1 = ltab_u(ju1)%iw1
   iu1 = ju1
   iu2 = nest_u(ju1)%iu

   itab_u(iu1)%iw1 = iw1
   itab_u(iu2)%iw1 = iw2
else                              ! Positive ju1 points outward
   iw1 = ltab_u(ju1)%iw2
   iu1 = nest_u(ju1)%iu
   iu2 = ju1

   itab_u(iu1)%iw2 = iw1
   itab_u(iu2)%iw2 = iw2
endif

if (ju1 == ltab_w(iw1)%iu1) then
   iu3 = ltab_w(iw1)%iu2
   iu5 = ltab_w(iw1)%iu3
elseif (ju1 == ltab_w(iw1)%iu2) then
   iu3 = ltab_w(iw1)%iu3
   iu5 = ltab_w(iw1)%iu1
else
   iu3 = ltab_w(iw1)%iu1
   iu5 = ltab_w(iw1)%iu2
endif

! Fill remaining neighbor indices for U and W points

! nothing needed for iu3

if (iw1 == itab_u(iu5)%iw1) then
   im4 = itab_u(iu5)%im2
   itab_u(iu5)%iw1 = iw2
else
   im4 = itab_u(iu5)%im1
   itab_u(iu5)%iw2 = iw2
endif

itab_u(iu4)%im1 = im2
itab_u(iu4)%im2 = im4
itab_u(iu4)%iw1 = iw1
itab_u(iu4)%iw2 = iw2

itab_w(iw1)%iu1 = iu1
itab_w(iw1)%iu2 = iu3
itab_w(iw1)%iu3 = iu4

itab_w(iw2)%iu1 = iu2
itab_w(iw2)%iu2 = iu4
itab_w(iw2)%iu3 = iu5

! Fill mrl and inudp values for new W points

if (itab_w(iw1)%mrlw /= mrloo) go to 5

itab_w(iw2)%mrlw      = itab_w(iw1)%mrlw
itab_w(iw2)%mrlw_orig = itab_w(iw1)%mrlw_orig
itab_w(iw2)%inudp(1)  = ltab_w(iw1)%inudp(1)
         
! Fill loop indices for new U points

itab_u(iu4)%iup = iu4
call uloops('f',iu4, 1, 4, 5, 7, 8,11,12,13,14,15)
call uloops('n',iu4,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

return

5 continue

write(io6,*) ''
write(io6, '(1x,A)')    'Error in subroutine perim_fill_cent1.'
write(io6, '(1x,A,I0)') 'Current nested grid ', ngr
write(io6, '(1x,A)')    'crosses pre-existing grid boundary.'
write(io6, '(1x,A,I0)') 'iw1 = ',iw1
stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_cent1

!===========================================================================

subroutine perim_fill_cent2(ngr,mrloo,kma,kua,kwa  &
                           ,jm1,ju1,jm2,ju2,jm3,jmo1,juo1,jmo2)

use mem_grid, only:  xem, yem, zem
use mem_ijtabs, only: itab_w, itab_u, ltab_w, ltab_u, nest_w, nest_u
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa   ! Index values for latest added points

integer, intent(in) :: jm1,jm2,jm3   ! Original M indices

integer, intent(in) :: ju1,ju2       ! Original U indices

integer, intent(out) :: jmo1,juo1,jmo2  ! Original M,U indices for next row 

integer :: ju4,ju7

integer :: iw1,iw2,iw3,iw4,iw5,iw6          ! Temporary new W indices
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8  ! Temporary new U indices
integer :: iu9,iu10,iu11,iu12,iu13          ! Temporary new U indices
integer :: im1,im2,im3,im4,im5,im6,im7,im8  ! Temporary new M indices

! Assign new I indices for all M, W points, and for U points that need no checks

im1 = jm1
im2 = nest_u(ju1)%im
im3 = jm2
im4 = nest_u(ju2)%im
im5 = jm3

im7 = kma + 1  ! Newly added point

iu6  = kua + 1  ! Newly added point
iu8  = kua + 2  ! Newly added point
iu10 = kua + 3  ! Newly added point

iw2 = kwa + 1  ! Newly added point
iw4 = kwa + 2  ! Newly added point
iw6 = kwa + 3  ! Newly added point

! Determine orientation of positive ju1 direction to get iw1 (=jw1)

if (jm1 == ltab_u(ju1)%im1) then  ! Positive ju1 points inward
   iw1 = ltab_u(ju1)%iw1
   iu1 = ju1
   iu2 = nest_u(ju1)%iu

!   itab_u(iu1)%iw1 = iw1  ! not needed?
else                              ! Positive ju1 points outward
   iw1 = ltab_u(ju1)%iw2
   iu1 = nest_u(ju1)%iu
   iu2 = ju1

   itab_u(iu1)%iw2 = iw1
endif

! Determine neighbors of iw1 (=jw1)

if (ju1 == ltab_w(iw1)%iu1) then
   iu5 = ltab_w(iw1)%iu2
   ju4 = ltab_w(iw1)%iu3
elseif (ju1 == ltab_w(iw1)%iu2) then
   iu5 = ltab_w(iw1)%iu3
   ju4 = ltab_w(iw1)%iu1
else
   iu5 = ltab_w(iw1)%iu1
   ju4 = ltab_w(iw1)%iu2
endif

iu7 = ju4

! Determine iw3 (=jw2)

if (jm2 == ltab_u(ju4)%im1) then  ! Positive ju4 points to jw2
   im6 = ltab_u(ju4)%im2
   iw3 = ltab_u(ju4)%iw2
else                              ! Positive ju4 points to jw1
   im6 = ltab_u(ju4)%im1
   iw3 = ltab_u(ju4)%iw1
endif

! Revisit orientation of positive ju1 direction now that iw3 is known

if (jm1 == ltab_u(ju1)%im1) then  ! Positive ju1 points inward
   itab_u(iu2)%iw1 = iw3
else                              ! Positive ju1 points outward
   itab_u(iu2)%iw2 = iw3
endif

! Determine orientation of positive ju2 direction to get iw5 (=jw3)

if (jm2 == ltab_u(ju2)%im1) then  ! Positive ju2 points inward
   iw5 = ltab_u(ju2)%iw1

   iu3 = ju2
   iu4 = nest_u(ju2)%iu

   itab_u(iu3)%iw1 = iw4
   itab_u(iu4)%iw1 = iw6
else                              ! Positive ju2 points outward
   iw5 = ltab_u(ju2)%iw2

   iu3 = nest_u(ju2)%iu
   iu4 = ju2

   itab_u(iu3)%iw2 = iw4
   itab_u(iu4)%iw2 = iw6
endif

! Determine neighbors of iw5 (=jw3)

if (ju2 == ltab_w(iw5)%iu1) then
   iu9 = ltab_w(iw5)%iu2
   iu11 = ltab_w(iw5)%iu3
elseif (ju2 == ltab_w(iw5)%iu2) then
   iu9 = ltab_w(iw5)%iu3
   iu11 = ltab_w(iw5)%iu1
else
   iu9 = ltab_w(iw5)%iu1
   iu11 = ltab_w(iw5)%iu2
endif

! Determine ju7 as neighbor of jw2 (=iw3)

if (ju4 == ltab_w(iw3)%iu1) then
   ju7 = ltab_w(iw3)%iu2
elseif (ju4 == ltab_w(iw3)%iu2) then
   ju7 = ltab_w(iw3)%iu3
else
   ju7 = ltab_w(iw3)%iu1
endif

! Determine orientation of positive ju7 direction

if (im6 == ltab_u(ju7)%im1) then  ! Positive ju7 points inward
   im8 = ltab_u(ju7)%im2
   iu12 = ju7
   iu13 = kua + 4  ! Newly added point

   itab_u(iu12)%im1 = im6 ! not needed?
   itab_u(iu12)%im2 = im7
   itab_u(iu12)%iw2 = iw2

   itab_u(iu13)%im1 = im7
   itab_u(iu13)%im2 = im8
   itab_u(iu13)%iw2 = iw5
else
   im8 = ltab_u(ju7)%im1
   iu12 = kua + 4  ! Newly added point
   iu13 = ju7

   itab_u(iu12)%im1 = im7
   itab_u(iu12)%im2 = im6
   itab_u(iu12)%iw1 = iw2

   itab_u(iu13)%im1 = im8
   itab_u(iu13)%im2 = im7
   itab_u(iu13)%iw1 = iw5
endif

! Fill remaining neighbor indices for U and W points

! nothing needed for iu5

itab_u(iu6)%im1 = im2
itab_u(iu6)%im2 = im6
itab_u(iu6)%iw1 = iw1
itab_u(iu6)%iw2 = iw2

itab_u(iu7)%im1 = im2
itab_u(iu7)%im2 = im7
itab_u(iu7)%iw1 = iw2
itab_u(iu7)%iw2 = iw3

itab_u(iu8)%im1 = im3
itab_u(iu8)%im2 = im7
itab_u(iu8)%iw1 = iw3
itab_u(iu8)%iw2 = iw4

itab_u(iu9)%im1 = im4
itab_u(iu9)%im2 = im7
itab_u(iu9)%iw1 = iw4
itab_u(iu9)%iw2 = iw5

itab_u(iu10)%im1 = im4
itab_u(iu10)%im2 = im8
itab_u(iu10)%iw1 = iw5
itab_u(iu10)%iw2 = iw6

if (iw5 == ltab_u(iu11)%iw1) then
   itab_u(iu11)%iw1 = iw6
else
   itab_u(iu11)%iw2 = iw6
endif

itab_w(iw1)%iu1 = iu1
itab_w(iw1)%iu2 = iu5
itab_w(iw1)%iu3 = iu6

itab_w(iw2)%iu1 = iu6
itab_w(iw2)%iu2 = iu12
itab_w(iw2)%iu3 = iu7

itab_w(iw3)%iu1 = iu2
itab_w(iw3)%iu2 = iu7
itab_w(iw3)%iu3 = iu8

itab_w(iw4)%iu1 = iu3
itab_w(iw4)%iu2 = iu8
itab_w(iw4)%iu3 = iu9

itab_w(iw5)%iu1 = iu9
itab_w(iw5)%iu2 = iu13
itab_w(iw5)%iu3 = iu10

itab_w(iw6)%iu1 = iu4
itab_w(iw6)%iu2 = iu10
itab_w(iw6)%iu3 = iu11

! Fill earth coordinates for new M point

xem(im7) = .5 * (xem(im6) + xem(im8))
yem(im7) = .5 * (yem(im6) + yem(im8))
zem(im7) = .5 * (zem(im6) + zem(im8))

! Fill mrl and inudp values for new W points

if (itab_w(iw1)%mrlw /= mrloo) go to 5
if (itab_w(iw3)%mrlw /= mrloo) go to 5
if (itab_w(iw5)%mrlw /= mrloo) go to 5

itab_w(iw2)%mrlw      = itab_w(iw1)%mrlw
itab_w(iw2)%mrlw_orig = itab_w(iw1)%mrlw_orig
itab_w(iw2)%inudp(1)  = ltab_w(iw1)%inudp(1)
         
itab_w(iw4)%mrlw      = itab_w(iw3)%mrlw
itab_w(iw4)%mrlw_orig = itab_w(iw3)%mrlw_orig
itab_w(iw4)%inudp(1)  = ltab_w(iw3)%inudp(1)
         
itab_w(iw6)%mrlw      = itab_w(iw5)%mrlw
itab_w(iw6)%mrlw_orig = itab_w(iw5)%mrlw_orig
itab_w(iw6)%inudp(1)  = ltab_w(iw5)%inudp(1)
         
! Fill loop indices for new U points

itab_u(iu6)%iup = iu6
call uloops('f',iu6, 1, 4, 5, 7, 8,11,12,13,14,15)
call uloops('n',iu6,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

itab_u(iu8)%iup = iu8
call uloops('f',iu8, 1, 4, 5, 7, 8,11,12,13,14,15)
call uloops('n',iu8,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

itab_u(iu10)%iup = iu10
call uloops('f',iu10, 1, 4, 5, 7, 8,11,12,13,14,15)
call uloops('n',iu10,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

itab_u(kua+4)%iup = kua+4
call uloops('f',kua+4, 1, 4, 5, 7, 8,11,12,13,14,15)
call uloops('n',kua+4,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

! Set JM and JU group indices in preparation for next row

jmo1 = im6
jmo2 = im8
juo1 = ju7

nest_u(ju7)%im = im7
nest_u(ju7)%iu = kua + 4

! Increment indices for newly added points

kma = kma + 1
kua = kua + 4
kwa = kwa + 3

return

5 continue

write(io6,*) 'In subroutine perim_fill_cent2, current nested grid ',ngr
write(io6,*) 'crosses pre-existing grid boundary.'
write(io6,*) 'iw1 = ',iw1
write(io6,*) 'stopping model'
stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_cent2

!===========================================================================

subroutine perim_fill_right(ngr,mrloo,kma,kua,kwa,jm1,ju1,jm2,jmo1,juo1,jmo2)

use mem_grid,   only:  xem, yem, zem
use mem_ijtabs, only: itab_w, itab_u, ltab_w, ltab_u, nest_w, nest_u
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa  ! Index values for latest added points

integer, intent(in) :: jm1,jm2  ! Original M indices

integer, intent(in) :: ju1  ! Original U indices

integer, intent(out) :: jmo1,juo1,jmo2  ! Original M,U indices for next row 

integer :: ju3,ju5

integer :: iw1,iw2,iw3,iw4                      ! Temporary new W indices
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9  ! Temporary new U indices
integer :: im1,im2,im3,im4,im5,im6              ! Temporary new M indices

! Assign new I indices for all M, W points, and for U points that need no checks

im1 = jm1
im2 = nest_u(ju1)%im
im3 = jm2

im5 = kma + 1  ! Newly added point

iu4  = kua + 1  ! Newly added point
iu6  = kua + 2  ! Newly added point

iw2 = kwa + 1  ! Newly added point
iw4 = kwa + 2  ! Newly added point

! Determine orientation of positive ju1 direction to get iw1 (=jw1)

if (jm1 == ltab_u(ju1)%im1) then  ! Positive ju1 points inward
   iw1 = ltab_u(ju1)%iw1
   iu1 = ju1
   iu2 = nest_u(ju1)%iu

   itab_u(iu1)%iw1 = iw1 ! not needed?
else                              ! Positive ju1 points outward
   iw1 = ltab_u(ju1)%iw2
   iu1 = nest_u(ju1)%iu
   iu2 = ju1

   itab_u(iu1)%iw2 = iw1
endif

! Determine neighbors of iw1 (=jw1)

if (ju1 == ltab_w(iw1)%iu1) then
   iu3 = ltab_w(iw1)%iu2
   ju3 = ltab_w(iw1)%iu3
elseif (ju1 == ltab_w(iw1)%iu2) then
   iu3 = ltab_w(iw1)%iu3
   ju3 = ltab_w(iw1)%iu1
else
   iu3 = ltab_w(iw1)%iu1
   ju3 = ltab_w(iw1)%iu2
endif

iu5 = ju3

! Determine iw3 (=jw2)

if (jm2 == ltab_u(ju3)%im1) then  ! Positive ju3 points to jw2
   im4 = ltab_u(ju3)%im2
   iw3 = ltab_u(ju3)%iw2
else                              ! Positive ju3 points to jw1
   im4 = ltab_u(ju3)%im1
   iw3 = ltab_u(ju3)%iw1
endif

! Revisit orientation of positive ju1 direction now that iw3 is known

if (jm1 == ltab_u(ju1)%im1) then  ! Positive ju1 points inward
   itab_u(iu2)%iw1 = iw3
else                              ! Positive ju1 points outward
   itab_u(iu2)%iw2 = iw3
endif

! Determine neighbors of iw3 (=jw2)

if (ju3 == ltab_w(iw3)%iu1) then
   ju5 = ltab_w(iw3)%iu2
   iu7 = ltab_w(iw3)%iu3
elseif (ju3 == ltab_w(iw3)%iu2) then
   ju5 = ltab_w(iw3)%iu3
   iu7 = ltab_w(iw3)%iu1
else
   ju5 = ltab_w(iw3)%iu1
   iu7 = ltab_w(iw3)%iu2
endif

! Determine orientation of positive ju5 direction

if (im4 == ltab_u(ju5)%im1) then  ! Positive ju5 points inward
   im6 = ltab_u(ju5)%im2
   iu8 = ju5
   iu9 = kua + 3  ! Newly added point

   itab_u(iu8)%im1 = im4
   itab_u(iu8)%im2 = im5
   itab_u(iu8)%iw2 = iw2

   itab_u(iu9)%im1 = im5
   itab_u(iu9)%im2 = im6
   itab_u(iu9)%iw2 = iw4
else
   im6 = ltab_u(ju5)%im1
   iu8 = kua + 3  ! Newly added point
   iu9 = ju5

   itab_u(iu8)%im1 = im5
   itab_u(iu8)%im2 = im4
   itab_u(iu8)%iw1 = iw2

   itab_u(iu9)%im1 = im6
   itab_u(iu9)%im2 = im5
   itab_u(iu9)%iw1 = iw4
endif

! Fill remaining neighbor indices for U and W points

! nothing needed for iu3

itab_u(iu4)%im1 = im2
itab_u(iu4)%im2 = im4
itab_u(iu4)%iw1 = iw1
itab_u(iu4)%iw2 = iw2

itab_u(iu5)%im1 = im2
itab_u(iu5)%im2 = im5
itab_u(iu5)%iw1 = iw2
itab_u(iu5)%iw2 = iw3

itab_u(iu6)%im1 = im3
itab_u(iu6)%im2 = im5
itab_u(iu6)%iw1 = iw3
itab_u(iu6)%iw2 = iw4

if (iw3 == ltab_u(iu7)%iw1) then
   itab_u(iu7)%iw1 = iw4
else
   itab_u(iu7)%iw2 = iw4
endif

itab_w(iw1)%iu1 = iu1
itab_w(iw1)%iu2 = iu3
itab_w(iw1)%iu3 = iu4

itab_w(iw2)%iu1 = iu4
itab_w(iw2)%iu2 = iu8
itab_w(iw2)%iu3 = iu5

itab_w(iw3)%iu1 = iu2
itab_w(iw3)%iu2 = iu5
itab_w(iw3)%iu3 = iu6

itab_w(iw4)%iu1 = iu6
itab_w(iw4)%iu2 = iu9
itab_w(iw4)%iu3 = iu7

! Fill earth coordinates for new M point

xem(im5) = .5 * (xem(im4) + xem(im6))
yem(im5) = .5 * (yem(im4) + yem(im6))
zem(im5) = .5 * (zem(im4) + zem(im6))

! Fill mrl and inudp values for new W points

if (itab_w(iw1)%mrlw /= mrloo) go to 5
if (itab_w(iw3)%mrlw /= mrloo) go to 5

itab_w(iw2)%mrlw      = itab_w(iw1)%mrlw
itab_w(iw2)%mrlw_orig = itab_w(iw1)%mrlw_orig
itab_w(iw2)%inudp(1)  = ltab_w(iw1)%inudp(1)
         
itab_w(iw4)%mrlw      = itab_w(iw3)%mrlw
itab_w(iw4)%mrlw_orig = itab_w(iw3)%mrlw_orig
itab_w(iw4)%inudp(1)  = ltab_w(iw3)%inudp(1)
         
! Fill loop indices for new U points

itab_u(iu4)%iup = iu4
call uloops('f',iu4, 1, 4, 5, 7, 8,11,12,13,14,15)
call uloops('n',iu4,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

itab_u(iu6)%iup = iu6
call uloops('f',iu6, 1, 4, 5, 7, 8,11,12,13,14,15)
call uloops('n',iu6,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

itab_u(kua+3)%iup = kua+3
call uloops('f',kua+3, 1, 4, 5, 7, 8,11,12,13,14,15)
call uloops('n',kua+3,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

! Set JM and JU group indices in preparation for next row

jmo1 = im4
jmo2 = im6
juo1 = ju5

nest_u(ju5)%im = im5
nest_u(ju5)%iu = kua + 3

! Increment indices for newly added points

kma = kma + 1
kua = kua + 3
kwa = kwa + 2

return

5 continue

write(io6,*) 'In subroutine perim_fill_right, current nested grid ',ngr
write(io6,*) 'crosses pre-existing grid boundary.'
write(io6,*) 'iw1 = ',iw1
write(io6,*) 'stopping model'
stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_right

!===========================================================================

subroutine perim_fill_left(ngr,mrloo,kma,kua,kwa,jm1,ju1,jm2,jmo1,juo1,jmo2)

use mem_grid,   only:  xem, yem, zem
use mem_ijtabs, only: itab_w, itab_u, ltab_w, ltab_u, nest_w, nest_u
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa  ! Index values for latest added points

integer, intent(in) :: jm1,jm2  ! Original M indices
integer, intent(in) :: ju1      ! Original U indices

integer, intent(out) :: jmo1,juo1,jmo2  ! Original M,U indices for next row 


integer :: ju3,ju5

integer :: iw1,iw2,iw3,iw4                      ! Temporary new W indices
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9  ! Temporary new U indices
integer :: im1,im2,im3,im4,im5,im6              ! Temporary new M indices

! Assign new I indices for all M, W points, and for U points that need no checks

im1 = jm1
im2 = nest_u(ju1)%im
im3 = jm2

im5 = kma + 1  ! Newly added point

iu4  = kua + 1  ! Newly added point
iu6  = kua + 2  ! Newly added point

iw2 = kwa + 1  ! Newly added point
iw4 = kwa + 2  ! Newly added point

! Determine orientation of positive ju1 direction to get iw3 (=jw2)

if (jm1 == ltab_u(ju1)%im1) then  ! Positive ju1 points inward
   iw3 = ltab_u(ju1)%iw1
   iu1 = ju1
   iu2 = nest_u(ju1)%iu

   itab_u(iu1)%iw1 = iw2
   itab_u(iu2)%iw1 = iw4
else                              ! Positive ju1 points outward
   iw3 = ltab_u(ju1)%iw2
   iu1 = nest_u(ju1)%iu
   iu2 = ju1

   itab_u(iu1)%iw2 = iw2
   itab_u(iu2)%iw2 = iw4
endif

! Determine neighbors of iw3 (=jw2)

if (ju1 == ltab_w(iw3)%iu1) then
   ju3 = ltab_w(iw3)%iu2
   iu7 = ltab_w(iw3)%iu3
elseif (ju1 == ltab_w(iw3)%iu2) then
   ju3 = ltab_w(iw3)%iu3
   iu7 = ltab_w(iw3)%iu1
else
   ju3 = ltab_w(iw3)%iu1
   iu7 = ltab_w(iw3)%iu2
endif

iu5 = ju3

! Determine iw1 (=jw1)

if (jm1 == ltab_u(ju3)%im1) then  ! Positive ju3 points to jw2
   im6 = ltab_u(ju3)%im2
   iw1 = ltab_u(ju3)%iw1
else                              ! Positive ju3 points to jw1
   im6 = ltab_u(ju3)%im1
   iw1 = ltab_u(ju3)%iw2
endif

! Determine neighbors of iw1 (=jw1)

if (ju3 == ltab_w(iw1)%iu1) then
   iu3 = ltab_w(iw1)%iu2
   ju5 = ltab_w(iw1)%iu3
elseif (ju3 == ltab_w(iw1)%iu2) then
   iu3 = ltab_w(iw1)%iu3
   ju5 = ltab_w(iw1)%iu1
else
   iu3 = ltab_w(iw1)%iu1
   ju5 = ltab_w(iw1)%iu2
endif

! Determine orientation of positive ju5 direction

if (im6 == ltab_u(ju5)%im2) then  ! Positive ju5 points inward
   im4 = ltab_u(ju5)%im1
   iu8 = ju5
   iu9 = kua + 3  ! Newly added point

   itab_u(iu8)%im1 = im4
   itab_u(iu8)%im2 = im5
   itab_u(iu8)%iw2 = iw1

   itab_u(iu9)%im1 = im5
   itab_u(iu9)%im2 = im6
   itab_u(iu9)%iw2 = iw3
else
   im4 = ltab_u(ju5)%im2
   iu8 = kua + 3  ! Newly added point
   iu9 = ju5

   itab_u(iu8)%im1 = im5
   itab_u(iu8)%im2 = im4
   itab_u(iu8)%iw1 = iw1

   itab_u(iu9)%im1 = im6
   itab_u(iu9)%im2 = im5
   itab_u(iu9)%iw1 = iw3
endif

! Fill remaining neighbor indices for U and W points

! nothing needed for iu3

itab_u(iu4)%im1 = im1
itab_u(iu4)%im2 = im5
itab_u(iu4)%iw1 = iw1
itab_u(iu4)%iw2 = iw2

itab_u(iu5)%im1 = im2
itab_u(iu5)%im2 = im5
itab_u(iu5)%iw1 = iw2
itab_u(iu5)%iw2 = iw3

itab_u(iu6)%im1 = im2
itab_u(iu6)%im2 = im6
itab_u(iu6)%iw1 = iw3
itab_u(iu6)%iw2 = iw4

if (iw3 == ltab_u(iu7)%iw1) then
   itab_u(iu7)%iw1 = iw4
else
   itab_u(iu7)%iw2 = iw4
endif

itab_w(iw1)%iu1 = iu3
itab_w(iw1)%iu2 = iu8
itab_w(iw1)%iu3 = iu4

itab_w(iw2)%iu1 = iu1
itab_w(iw2)%iu2 = iu4
itab_w(iw2)%iu3 = iu5

itab_w(iw3)%iu1 = iu5
itab_w(iw3)%iu2 = iu9
itab_w(iw3)%iu3 = iu6

itab_w(iw4)%iu1 = iu2
itab_w(iw4)%iu2 = iu6
itab_w(iw4)%iu3 = iu7

! Fill earth coordinates for new M point

xem(im5) = .5 * (xem(im4) + xem(im6))
yem(im5) = .5 * (yem(im4) + yem(im6))
zem(im5) = .5 * (zem(im4) + zem(im6))

! Fill mrl and inudp values for new W points

if (itab_w(iw1)%mrlw /= mrloo) go to 5
if (itab_w(iw3)%mrlw /= mrloo) go to 5

itab_w(iw2)%mrlw      = itab_w(iw1)%mrlw
itab_w(iw2)%mrlw_orig = itab_w(iw1)%mrlw_orig
itab_w(iw2)%inudp(1)  = ltab_w(iw1)%inudp(1)
         
itab_w(iw4)%mrlw      = itab_w(iw3)%mrlw
itab_w(iw4)%mrlw_orig = itab_w(iw3)%mrlw_orig
itab_w(iw4)%inudp(1)  = ltab_w(iw3)%inudp(1)
         
! Fill loop indices for new U points

itab_u(iu4)%iup = iu4
call uloops('f',iu4, 1, 4, 5, 7, 8,11,12,13,14,15)
call uloops('n',iu4,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

itab_u(iu6)%iup = iu6
call uloops('f',iu6, 1, 4, 5, 7, 8,11,12,13,14,15)
call uloops('n',iu6,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

itab_u(kua+3)%iup = kua+3
call uloops('f',kua+3, 1, 4, 5, 7, 8,11,12,13,14,15)
call uloops('n',kua+3,16,20, 0, 0, 0, 0, 0, 0, 0, 0)

! Set JM and JU group indices in preparation for next row

jmo1 = im4
jmo2 = im6
juo1 = ju5

nest_u(ju5)%im = im5
nest_u(ju5)%iu = kua + 3

! Increment indices for newly added points

kma = kma + 1
kua = kua + 3
kwa = kwa + 2

return

5 continue

write(io6,*) 'In subroutine perim_fill_left, current nested grid ',ngr
write(io6,*) 'crosses pre-existing grid boundary.'
write(io6,*) 'iw1 = ',iw1
write(io6,*) 'stopping model'
stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_left

!===========================================================================

subroutine perim_mrow()

use mem_ijtabs, only: itab_w, itab_u
use mem_grid,   only: nwa
use misc_coms,  only: io6

implicit none

integer :: iper, iu, iw, iw1, iw2, iw3, iter, irow, jrow, mrow, mrowh

integer :: mrow_temp(nwa),mrowh_temp(nwa)

! Loop over all W points

do iw = 2,nwa

! Initialize all mrow values to zero.

   itab_w(iw)%mrow = 0
   mrow_temp(iw)   = 0
   mrowh_temp(iw)  = 0

! Set values on nested grid border to +/- 1.

   iw1 = itab_w(iw)%iw1
   iw2 = itab_w(iw)%iw2
   iw3 = itab_w(iw)%iw3
   
   if     (itab_w(iw)%mrlw < itab_w(iw1)%mrlw .or.  &
           itab_w(iw)%mrlw < itab_w(iw2)%mrlw .or.  &
           itab_w(iw)%mrlw < itab_w(iw3)%mrlw) then

      itab_w(iw)%mrow = 1
      itab_w(iw)%mrowh = 1
      mrow_temp(iw) = 1
      mrowh_temp(iw) = 1

   elseif (itab_w(iw)%mrlw > itab_w(iw1)%mrlw .or.  &
           itab_w(iw)%mrlw > itab_w(iw2)%mrlw .or.  &
           itab_w(iw)%mrlw > itab_w(iw3)%mrlw) then

      itab_w(iw)%mrow = -1
      itab_w(iw)%mrowh = -1
      mrow_temp(iw) = -1
      mrowh_temp(iw) = -1

   endif

enddo

do irow = 2,10  ! First row already done above
   jrow = mod(irow,2)

   do iw = 2,nwa

      if (itab_w(iw)%mrow == 0) then

! If IW is adjacent to any other IW cell with nonzero mrow, 
! set mrow for IW cell 

         iw1 = itab_w(iw)%iw1
         iw2 = itab_w(iw)%iw2
         iw3 = itab_w(iw)%iw3

! Check for positive mrow & mrowh values

         mrow = max(itab_w(iw1)%mrow  &
                   ,itab_w(iw2)%mrow  &
                   ,itab_w(iw3)%mrow)

         mrowh = max(itab_w(iw1)%mrowh  &
                   ,itab_w(iw2)%mrowh  &
                   ,itab_w(iw3)%mrowh)

         if (mrow > 0)  mrow_temp (iw) = mrow + jrow
         if (mrowh > 0) mrowh_temp(iw) = mrowh + 1

! Check for negative mrow & mrowh values

         mrow = min(itab_w(iw1)%mrow  &
                   ,itab_w(iw2)%mrow  &
                   ,itab_w(iw3)%mrow)

         mrowh = min(itab_w(iw1)%mrowh  &
                    ,itab_w(iw2)%mrowh  &
                    ,itab_w(iw3)%mrowh)

         if (mrow < 0)  mrow_temp (iw) = mrow - jrow
         if (mrowh < 0) mrowh_temp(iw) = mrowh - 1

      endif

   enddo

   do iw = 2,nwa
      itab_w(iw)%mrow  = mrow_temp(iw)
      itab_w(iw)%mrowh = mrowh_temp(iw)
   enddo

enddo

return
end subroutine perim_mrow







