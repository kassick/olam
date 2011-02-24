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
subroutine hs_hist(ihsflag)

use mem_basic
use mem_grid
use mem_ijtabs
use consts_coms

implicit none

integer, intent(in) :: ihsflag

integer, parameter :: modlevs = 25
integer :: ilat,k,iw,notavail
real :: wt1,wt2,fldval,rilat
integer, save :: jcall = 40
real, dimension(4) :: hcpn,vcpn,fldvals

real, save :: atime=0.,dummy(1)=0.

real, save, allocatable, dimension(:,:) :: uzonal_tavg,umerid_tavg  &
                                          ,temp_tavg,press_tavg     &
                                          ,temp2_tavg
real, save, dimension(90,modlevs) :: uzonal_tlon,umerid_tlon,temp_tlon  &
                                    ,press_tlon,temp2_tlon
real, save, dimension(90) :: wtlat
real, save, dimension(modlevs) :: pcol,pcolin,vctr1,vctr2,vctr3,vctr4  &
                                 ,vctr5,vctr6,vctr7,vctr8
character*1 :: title

real, save :: aspect = 0.7, scalelab = .016

! On first call, allocate averaging arrays and initialize to zero

if (jcall == 40) then
   jcall = 10
   allocate(uzonal_tavg(mza,mwa),umerid_tavg(mza,mwa)  &
           ,temp_tavg(mza,mwa),press_tavg(mza,mwa),temp2_tavg(mza,mwa))

   uzonal_tavg (1:mza,1:mwa) = 0.
   umerid_tavg (1:mza,1:mwa) = 0.
   temp_tavg   (1:mza,1:mwa) = 0.
   press_tavg  (1:mza,1:mwa) = 0.
   temp2_tavg  (1:mza,1:mwa) = 0.

! Fill pcol levels for vertical interpolation

   do k = 1,modlevs
!      pcol(k) = -25. + 26. * float(k)  ! range 1 to 989 for modlevs = 39
!      pcol(k) = -38. + 39. * float(k)  ! range 1 to 976 for modlevs = 26
!      pcol(k) = -31. + 32. * float(k)  ! range 1 to 993 for modlevs = 32
!      pcol(k) = -40. + 41. * float(k)  ! range 1 to 985 for modlevs = 25
   enddo
endif

do iw = 2,mwa
   do k = lpw(iw),mza-1
   
      call oplot_lib(1,k,iw,'VALUE','ZONAL_WINDW',fldval,notavail)    
      uzonal_tavg(k,iw) = uzonal_tavg(k,iw) * atime / (atime + 1.)  &
                        + fldval / (atime + 1.)  

      call oplot_lib(1,k,iw,'VALUE','MERID_WINDW',fldval,notavail)    
      umerid_tavg(k,iw) = umerid_tavg(k,iw) * atime / (atime + 1.)  &
                        + fldval / (atime + 1.) 

      call oplot_lib(1,k,iw,'VALUE','AIRTEMPK',fldval,notavail)    
      temp_tavg(k,iw)  = temp_tavg(k,iw) * atime / (atime + 1.)  &
                       + fldval / (atime + 1.) 
      temp2_tavg(k,iw) = temp2_tavg(k,iw) * atime / (atime + 1.)  & 
                       + fldval * fldval / (atime + 1.) 

      call oplot_lib(1,k,iw,'VALUE','PRESS',fldval,notavail)    
      press_tavg(k,iw) = press_tavg(k,iw) * atime / (atime + 1.)  & 
                       + fldval / (atime + 1.) 

   enddo
enddo

atime = atime + 1.

print*, 'atime ',atime

! Return if this is not the last history file read time

if (ihsflag == 0) return

! Reopen the current graphics output workstation

call o_reopnwk()

! Initialize lat/height arrays to zero prior to summation

uzonal_tlon (1:90,1:modlevs) = 0.
umerid_tlon (1:90,1:modlevs) = 0.
temp_tlon   (1:90,1:modlevs) = 0.
press_tlon  (1:90,1:modlevs) = 0.
temp2_tlon  (1:90,1:modlevs) = 0.

wtlat (1:90) = 0.

! Compute zonal averages of time-averaged fields

do iw = 2,mwa

   rilat = .5 * (glatw(iw) + 90.) + .5
   ilat = int(rilat)
   wt2 = rilat - real(ilat)

! Correct latitude index and weights for latitudes close to either pole

   if (ilat == 0) then
      ilat = 1
      wt2 = 0.
   elseif(ilat == 90) then
      ilat = 89
      wt2 = 1.
   endif

   wt1 = 1. - wt2

   wtlat(ilat)   = wtlat(ilat)   + wt1
   wtlat(ilat+1) = wtlat(ilat+1) + wt2

   do k = lpw(iw),mza-1

      uzonal_tlon(ilat  ,k-1) = uzonal_tlon(ilat  ,k-1) + wt1 * uzonal_tavg(k,iw)
      uzonal_tlon(ilat+1,k-1) = uzonal_tlon(ilat+1,k-1) + wt2 * uzonal_tavg(k,iw)

      umerid_tlon(ilat  ,k-1) = umerid_tlon(ilat  ,k-1) + wt1 * umerid_tavg(k,iw)
      umerid_tlon(ilat+1,k-1) = umerid_tlon(ilat+1,k-1) + wt2 * umerid_tavg(k,iw)

      temp_tlon  (ilat  ,k-1) = temp_tlon  (ilat  ,k-1) + wt1 * temp_tavg(k,iw)
      temp_tlon  (ilat+1,k-1) = temp_tlon  (ilat+1,k-1) + wt2 * temp_tavg(k,iw)

      press_tlon (ilat  ,k-1) = press_tlon (ilat  ,k-1) + wt1 * press_tavg(k,iw)
      press_tlon (ilat+1,k-1) = press_tlon (ilat+1,k-1) + wt2 * press_tavg(k,iw)

      temp2_tlon (ilat  ,k-1) = temp2_tlon (ilat  ,k-1) + wt1 * temp2_tavg(k,iw)
      temp2_tlon (ilat+1,k-1) = temp2_tlon (ilat+1,k-1) + wt2 * temp2_tavg(k,iw)

   enddo
enddo

! Normalize fields by wtlat array

do ilat = 1,90
   do k = 1,modlevs
      uzonal_tlon(ilat,k) = uzonal_tlon(ilat,k) / wtlat(ilat)
      umerid_tlon(ilat,k) = umerid_tlon(ilat,k) / wtlat(ilat)
      temp_tlon  (ilat,k) = temp_tlon  (ilat,k) / wtlat(ilat)
      press_tlon (ilat,k) = press_tlon (ilat,k) / wtlat(ilat)
      temp2_tlon (ilat,k) = temp2_tlon (ilat,k) / wtlat(ilat)

! Convert temp2_tlon to a variance

      temp2_tlon (ilat,k) = max(1.e-6,temp2_tlon (ilat,k)  &
                          - temp_tlon(ilat,k) * temp_tlon(ilat,k))

print*, 'hs_hist ',ilat,k,uzonal_tlon(ilat,k)


   enddo
enddo

! Vertically interpolate fields to constant pressure surfaces

!do ilat = 1,90
!   do k = 1,modlevs
!      vctr1(k)  = uzonal_tlon(ilat,modlevs+1-k)
!      vctr2(k)  = umerid_tlon(ilat,modlevs+1-k)
!      vctr3(k)  = temp_tlon  (ilat,modlevs+1-k)
!      vctr4(k)  = temp2_tlon (ilat,modlevs+1-k)
!      pcolin(k) = press_tlon (ilat,modlevs+1-k)
!   enddo

!   call hintrp_ee(modlevs,vctr1,pcolin,modlevs,vctr5,pcol)
!   call hintrp_ee(modlevs,vctr2,pcolin,modlevs,vctr6,pcol)
!   call hintrp_ee(modlevs,vctr3,pcolin,modlevs,vctr7,pcol)
!   call hintrp_ee(modlevs,vctr4,pcolin,modlevs,vctr8,pcol)

!   do k = 1,modlevs
!      uzonal_tlon(ilat,k) = vctr5(modlevs+1-k)
!      umerid_tlon(ilat,k) = vctr6(modlevs+1-k)
!      temp_tlon  (ilat,k) = vctr7(modlevs+1-k)
!      temp2_tlon (ilat,k) = vctr8(modlevs+1-k)
!   enddo
!enddo

! Set plot color (black)

call o_gsplci(10)
call o_gsfaci(10)
call o_sflush()

! Plot zonal/time averages

call plotback()

call oplot_xy2('1','a',aspect,scalelab    &
      ,1,dummy,dummy                      &
      ,'latitude (deg)','pressure (mb)'   &
      ,-90.,90.,10.,3   ,1000.,0.,-50.,5  )

! Fill horiz and vert coord values and field values

do ilat = 1,89
   do k = 1,modlevs-1  ! This loop centered at W pts

      hcpn(1) = ilat * 2 - 91.
      hcpn(2) = ilat * 2 - 89.
      hcpn(3) = hcpn(2)
      hcpn(4) = hcpn(1)

      vcpn(1) = press_tlon(ilat  ,k  )
      vcpn(2) = press_tlon(ilat+1,k  )
      vcpn(3) = press_tlon(ilat+1,k+1)
      vcpn(4) = press_tlon(ilat  ,k+1)

      fldvals(1) = uzonal_tlon(ilat  ,k  )
      fldvals(2) = uzonal_tlon(ilat+1,k  )
      fldvals(3) = uzonal_tlon(ilat+1,k+1)
      fldvals(4) = uzonal_tlon(ilat  ,k+1)

      call contpolyg(91,1,4,hcpn,vcpn,fldvals)
      call contpolyg(91,0,4,hcpn,vcpn,fldvals)

   enddo
enddo

call oplot_xy2('2','b',aspect,scalelab     &
      ,1,dummy,dummy                       &
      ,'latitude (deg)','pressure (mb)'    &
      ,-90.,90.,10.,3   ,1000.,0.,-50.,5  )

! Fill horiz and vert coord values and field values

do ilat = 1,89
   do k = 1,modlevs-1  ! Loop for this field is over W levels

      hcpn(1) = ilat * 2 - 91.
      hcpn(2) = ilat * 2 - 89.
      hcpn(3) = hcpn(2)
      hcpn(4) = hcpn(1)

      vcpn(1) = press_tlon(ilat  ,k  )
      vcpn(2) = press_tlon(ilat+1,k  )
      vcpn(3) = press_tlon(ilat+1,k+1)
      vcpn(4) = press_tlon(ilat  ,k+1)

      fldvals(1) = temp_tlon(ilat  ,k  )
      fldvals(2) = temp_tlon(ilat+1,k  )
      fldvals(3) = temp_tlon(ilat+1,k+1)
      fldvals(4) = temp_tlon(ilat  ,k+1)

      call contpolyg(92,1,4,hcpn,vcpn,fldvals)
      call contpolyg(92,0,4,hcpn,vcpn,fldvals)

   enddo
enddo

call oplot_xy2('3','c',aspect,scalelab    &
      ,1,dummy,dummy                      &
      ,'latitude (deg)','pressure (mb)'   &
      ,-90.,90.,10.,3   ,1000.,0.,-50.,5  )

! Fill horiz and vert coord values and field values

do ilat = 1,89
   do k = 1,modlevs-1  ! Loop for this field is over W levels

      hcpn(1) = ilat * 2 - 91.
      hcpn(2) = ilat * 2 - 89.
      hcpn(3) = hcpn(2)
      hcpn(4) = hcpn(1)

      vcpn(1) = press_tlon(ilat  ,k  )
      vcpn(2) = press_tlon(ilat+1,k  )
      vcpn(3) = press_tlon(ilat+1,k+1)
      vcpn(4) = press_tlon(ilat  ,k+1)

      fldvals(1) = umerid_tlon(ilat  ,k  )
      fldvals(2) = umerid_tlon(ilat+1,k  )
      fldvals(3) = umerid_tlon(ilat+1,k+1)
      fldvals(4) = umerid_tlon(ilat  ,k+1)

      call contpolyg(93,1,4,hcpn,vcpn,fldvals)
      call contpolyg(93,0,4,hcpn,vcpn,fldvals)

   enddo
enddo

call oplot_xy2('4','d',aspect,scalelab     &
      ,1,dummy,dummy                       &
      ,'latitude (deg)','pressure (mb)'    &
      ,-90.,90.,30.,1   ,1000.,0.,-50.,5  )

! Fill horiz and vert coord values and field values

do ilat = 1,89
   do k = 1,modlevs-1  ! Loop for this field is over W levels

      hcpn(1) = ilat * 2 - 91.
      hcpn(2) = ilat * 2 - 89.
      hcpn(3) = hcpn(2)
      hcpn(4) = hcpn(1)

      vcpn(1) = press_tlon(ilat  ,k  )
      vcpn(2) = press_tlon(ilat+1,k  )
      vcpn(3) = press_tlon(ilat+1,k+1)
      vcpn(4) = press_tlon(ilat  ,k+1)

      fldvals(1) = temp2_tlon(ilat  ,k  )
      fldvals(2) = temp2_tlon(ilat+1,k  )
      fldvals(3) = temp2_tlon(ilat+1,k+1)
      fldvals(4) = temp2_tlon(ilat  ,k+1)

      call contpolyg(94,1,4,hcpn,vcpn,fldvals)
      call contpolyg(94,0,4,hcpn,vcpn,fldvals)

   enddo
enddo

call o_frame()

!call plotback()
!call cpcnrc(uzonal_tlon, 90, 90, modlevs, 0, 0, 4.0, 0, -1, -682)
!call frame

!call plotback()
!call cpcnrc(umerid_tlon, 90, 90, modlevs, 0, 0, 0.5, 0, -1, -682)
!call frame

!call plotback()
!call cpcnrc(temp_tlon,   90, 90, modlevs, 0, 0, 5.0, 0, -1, -682)
!call frame

!call plotback()
!call cpcnrc(temp2_tlon,  90, 90, modlevs, 0, 0, 5.0, 0, -1, -682)
!call frame

!call plotback()
!call cpcnrc(press_tlon,  90, 90, modlevs, 0, 0,10.0, 0, -1, -682)
!call frame

! Close the current workstation. This allows viewing the complete
! meta file (including the last frame) during a run and in case the
! simulation crashes.

call o_clswk()

return
end subroutine hs_hist

