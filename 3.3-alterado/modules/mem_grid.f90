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
Module mem_grid

   integer :: &  ! N values for full-domain reference on any process

      nza           &  ! Vertical number of all points
     ,nsw_max       &  ! Max # vert atm levels in IW column with sfc flux
     ,nma,nua,nwa      ! Horiz number of all M,U,W points in full domain

   integer :: &

      mza           &  ! Vertical number of all points
     ,mma,mua,mwa      ! Horiz number of all M,U,W points in sub-process

   integer, allocatable, dimension(:) ::  &

      lpu,lpw       &  ! Lowest prognosed U,W/T
     ,lcu           &  ! Lowest nonzero control volume for U
     ,lsw              ! number of W/T levels in contact with surface

   real, allocatable, dimension(:) ::  &

      zm    ,zt        &  ! Z coordinate at M,T point
     ,dzm   ,dzt       &  ! Delta Z at M,T point
     ,dzim  ,dzit      &  ! Delta Z inverse (1/dz) at M,T point

     ,xem  ,xeu  ,xew     &  ! XE coordinate at M,U,W point
     ,yem  ,yeu  ,yew     &  ! YE coordinate at M,U,W point
     ,zem  ,zeu  ,zew     &  ! ZE coordinate at M,U,W point

     ,unx, uny, unz     &  ! U point normal unit vector components
     ,wnx, wny, wnz     &  ! W point normal unit vector components

     ,utx, uty, utz     &  ! U point horizontal tangential unit vector components

     ,dtu              &  ! distance tangent to U edge (M1-M2 distance)
     ,dnu              &  ! dxy across U pt for PGF
     ,dniu             &  ! (1/dxy) across U pt for turb U XY flx
     ,arw0             &  ! Area of IW triangle at earth surface

     ,glatw ,glonw     &  ! Latitude/longitude at W point
     ,glatm ,glonm     &  ! Latitude/longitude at M point
     ,topm             &  ! Topography height at M point
     ,dts_w               ! pointwise acoustic timestep at W point

   real, allocatable, dimension(:,:) ::  &

      aru   ,arw           &  ! Aperture area of U,W face
     ,volui ,volwi            ! (1/Volume) of T,U,W cell

   real(kind=8), allocatable, dimension(:,:) ::  &

      volt                 &  ! Volume of T cell
     ,volti      

   integer :: impent(12) ! Scratch array for storing 12 pentagonal IM indices

   integer, parameter :: nrows = 5
   integer :: mrows

Contains

!===============================================================================

   subroutine alloc_gridz()

   implicit none

   allocate (zm   (mza) ,zt   (mza)  &
            ,dzm  (mza) ,dzt  (mza)  &
            ,dzim (mza) ,dzit (mza)  )

   return
   end subroutine alloc_gridz
   
!===============================================================================

   subroutine alloc_xyzem()

   implicit none

   allocate (xem(mma));  xem(1:mma) = 0.
   allocate (yem(mma));  yem(1:mma) = 0.
   allocate (zem(mma));  zem(1:mma) = 0.

   return
   end subroutine alloc_xyzem
   
!===============================================================================

   subroutine alloc_grid1()

   use misc_coms, only: io6

   implicit none
   
! Allocate and initialize arrays (xem, yem, zem are already allocated)

   allocate (lsw(mwa));  lsw(1:mwa) = 0

   allocate (xeu(mua));  xeu(1:mua) = 0.
   allocate (yeu(mua));  yeu(1:mua) = 0.
   allocate (zeu(mua));  zeu(1:mua) = 0.

   allocate (xew(mwa));  xew(1:mwa) = 0.
   allocate (yew(mwa));  yew(1:mwa) = 0.
   allocate (zew(mwa));  zew(1:mwa) = 0.

   allocate (unx(mua));  unx(1:mua) = 0.
   allocate (uny(mua));  uny(1:mua) = 0.
   allocate (unz(mua));  unz(1:mua) = 0.

   allocate (wnx(mwa));  wnx(1:mwa) = 0.
   allocate (wny(mwa));  wny(1:mwa) = 0.
   allocate (wnz(mwa));  wnz(1:mwa) = 0.

   allocate (utx(mua));  utx(1:mua) = 0.
   allocate (uty(mua));  uty(1:mua) = 0.
   allocate (utz(mua));  utz(1:mua) = 0.

   allocate (dnu  (mua));  dnu  (1:mua) = 0.
   allocate (dniu (mua));  dniu (1:mua) = 0.

   allocate (dtu (mua));  dtu (1:mua) = 0.
   allocate (arw0(mwa));  arw0(1:mwa) = 0.

   allocate (glatw(mwa));  glatw(1:mwa) = 0.
   allocate (glonw(mwa));  glonw(1:mwa) = 0.

   allocate  (topm(mma));   topm(1:mma) = 0.
   allocate (glatm(mma));  glatm(1:mma) = 0.
   allocate (glonm(mma));  glonm(1:mma) = 0.

   write(io6,*) 'finishing alloc_grid1'
            
   return
  
   end subroutine alloc_grid1

!===============================================================================

   subroutine alloc_grid2()

   use misc_coms, only: io6

   implicit none

! Allocate  and initialize arrays

   write(io6,*) 'alloc_grid2 ',mma,mua,mwa,mza

   allocate (lpu(mua));  lpu(1:mua) = 0
   allocate (lpw(mwa));  lpw(1:mwa) = 0  ! In vtables
   allocate (lcu(mua));  lcu(1:mua) = 0

   allocate (aru  (mza,mua));  aru  (1:mza,1:mua) = 0.
   allocate (arw  (mza,mwa));  arw  (1:mza,1:mwa) = 0.

   allocate (volt (mza,mwa));  volt (1:mza,1:mwa) = 0.
   allocate (volti(mza,mwa));  volti(1:mza,1:mwa) = 0.
   allocate (volui(mza,mua));  volui(1:mza,1:mua) = 0.
   allocate (volwi(mza,mwa));  volwi(1:mza,1:mwa) = 0.
            
   write(io6,*) 'finishing alloc_grid2'
            
   return
  
   end subroutine alloc_grid2

!===============================================================================

   subroutine dealloc_grid(nd)

   implicit none
   integer, intent(in) :: nd
   
   if (nd == 2) then

      if (allocated(lpu))  deallocate (lpu)
      if (allocated(lpw))  deallocate (lpw)
      if (allocated(lcu))  deallocate (lcu)
      if (allocated(lsw))  deallocate (lsw)

      if (allocated(zm)  ) deallocate (zm)
      if (allocated(zt)  ) deallocate (zt)
      if (allocated(dzm) ) deallocate (dzm)
      if (allocated(dzt) ) deallocate (dzt)
      if (allocated(dzim)) deallocate (dzim)
      if (allocated(dzit)) deallocate (dzit)

      if (allocated(dtu)  ) deallocate (dtu)
      if (allocated(dnu)  ) deallocate (dnu)
      if (allocated(dniu) ) deallocate (dniu)
      if (allocated(arw0) ) deallocate (arw0)

      if (allocated(glatw)) deallocate (glatw)
      if (allocated(glonw)) deallocate (glonw)

      if (allocated(topm) ) deallocate (topm)
      if (allocated(glatm)) deallocate (glatm)
      if (allocated(glonm)) deallocate (glonm)

      if (allocated(dts_w)) deallocate (dts_w)

      if (allocated(aru))   deallocate (aru) 
      if (allocated(arw))   deallocate (arw)
      
      if (allocated(volt))   deallocate (volt)
      if (allocated(volti))  deallocate (volti)
      if (allocated(volui))  deallocate (volui)
      if (allocated(volwi))  deallocate (volwi)
   
   endif

   return
   end subroutine dealloc_grid
   
End Module mem_grid

      
