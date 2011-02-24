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
subroutine topo_init()

use mem_grid, only: mma, topm, xem, yem, zem
use misc_coms, only: io6, deltax
use consts_coms, only: grav
use mem_ijtabs,  only: rotation_angle

!------------------------------------------
! Only for ncar test cases:
use leaf_coms, only: nslcon
use ncar_testcases_all, only: r8, ncar_testcase, mountain_rossby,  &
                              surface_geopotential
!------------------------------------------

implicit none

integer :: im

real :: hfwid
real :: hgt
real :: hfwid2

real :: glatm, glonm, raxis
real(r8) :: presss, uuu, vvv, tempn, phis, psfc

!------------------------------------------
! Only for ncar test cases:
ncar_testcase  = nslcon
!------------------------------------------

! Fill the TOPM array with a default value of 0 or modify it as desired.  
! If itopoflg is set to 1, these values will be overridden in the call to
! topo_database, which inputs a standard OLAM topography dataset.

hfwid = 10000.

! dudhia expts
! hfwid = 5. * .866 * deltax

! hgt = 405.
! hgt = 1012.
! end dudhia expts

hfwid2 = hfwid**2


print*, 'topo_init ',ncar_testcase,mma

do im = 2,mma
   topm(im) = 0.

!------------------------------------------
! Only for ncar test cases:

   if (ncar_testcase <= 2) then

! Find lat/lon of current M point

      raxis = sqrt(xem(im) ** 2 + yem(im) ** 2)  ! dist from earth axis

      glatm = atan2(zem(im),raxis)
      glonm = atan2(yem(im),xem(im))

      phis = surface_geopotential(real(glonm,r8),  &
                                  real(glatm,r8),  &
                                  rotation_angle   )

      topm(im) = phis / grav
      

      print*, 'topm5 ',im,topm(im),glatm,glonm,rotation_angle
      
   endif


   if (ncar_testcase == 5) then

! Find lat/lon of current M point

      raxis = sqrt(xem(im) ** 2 + yem(im) ** 2)  ! dist from earth axis

      glatm = atan2(zem(im),raxis)
      glonm = atan2(yem(im),xem(im))

      call mountain_Rossby (real(glonm,r8),   &
                            real(glatm,r8),   &
                            presss,           &
                            uuu, vvv, tempn,  &
                            phis,psfc         )

      topm(im) = phis / grav
      
!      print*, 'topm5 ',im,topm(im),glatm,glonm
      
   endif
!------------------------------------------

!   topm(im) = 200. * mod(im,4)

! SPECIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   topm(im) = max(0.,hgt * hfwid2 / (hfwid2 + xem(im)**2) - 1.)  
!   write(io6,*) 'topm ',im,xem(im),topm(im)
! TOPM = 0 AT LARGE DISTANCE FROM HILL
! END SPECIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

enddo

return
end subroutine topo_init
