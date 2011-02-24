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
subroutine uwcomp(itask,iu,k,flux,wparwx2,uc2)

! Computes vertical flux of U momentum for 10-meter mountain experiment
! assuming background speed of 20 m/s and isothermal atmosphere at 250 K.

! THIS SUBROUTINE IS NOT DESIGNED FOR PARALLEL EXECUTION.  
! (This experiment is a small computational problem that does not need parallel execution.)

use mem_grid,   only: mza, zm
use misc_coms,  only: dtsm, time_istp8, deltax, dn01d, runtype, mdomain
use mem_ijtabs, only: mrls
use oplot_coms, only: op

implicit none

integer, intent(in) :: itask,iu,k
real, intent(in) :: flux,wparwx2,uc2

real, save, allocatable, dimension(:) :: vctr11,vctr12,vctr13,vctr14  &
                                        ,vctr16,vctr17,vctr18,zmkm

integer, save :: ncall_uwcomp=0

integer :: kk
real, parameter :: bvfreq = .01957, pio4 = 3.14159 / 4.,compfreq = 3600.
real :: tval, c1, chanwid2
real, save :: aspect = 1.0, scalelab = .016


if (mod(time_istp8,dble(compfreq)) < compfreq - 1.5 * dtsm(mrls)) return

! Reopen the current graphics output workstation if it is closed
   call o_reopnwk()


! Allocate local memory on first call only

if (ncall_uwcomp /= 1) then
   ncall_uwcomp = 1
   allocate (vctr11(mza),vctr12(mza),vctr13(mza),vctr14(mza)  &
            ,vctr16(mza),vctr17(mza),vctr18(mza),zmkm(mza))
            
   zmkm(1:mza) = 1.e-3 * zm(1:mza)
endif

c1 = 2. / sqrt(3.)  ! factor to get U_eastward from UC

if (itask == 0) then

   print*, ' '
   print*, 'UWCOMP:  K    ZM(K)     FLUX    FLUXBAR  FLUX-FLUXBAR'
   print*, ' '

   vctr11(1:mza) = 0.
   vctr12(1:mza) = 0.
   vctr13(1:mza) = 0.
   vctr14(1:mza) = 0.

   vctr16(1:mza) = 0.
   vctr17(1:mza) = 0.
   vctr18(1:mza) = 0.

elseif (itask == 1) then

! Double channel width

   if (mdomain == 3) chanwid2 = 3. * deltax
   if (mdomain == 4) chanwid2 = 6. * deltax

   vctr11(k) = vctr11(k) - flux    * c1 / chanwid2
   vctr12(k) = vctr12(k) - wparwx2 * c1 / chanwid2
   vctr13(k) = vctr13(k) + uc2     * c1 / 2.
   vctr14(k) = vctr14(k) + 1.  ! Counter over iu points

elseif (itask == 2) then

!  1m   tval = pio4 * 20. * bvfreq * (1.**2) * dn01d(2) ! dn01d still used; ok here
   tval = pio4 * 20. * bvfreq * (10.**2) * dn01d(2) ! dn01d still used; ok here
! 20m   tval = pio4 * 20. * bvfreq * (20.**2) * dn01d(2) ! dn01d still used; ok here
! 50m   tval = pio4 * 20. * bvfreq * (50.**2) * dn01d(2) ! dn01d still used; ok here

! Normalize flux quantities and print values

   do kk = 2,mza-2

      vctr16(kk) = vctr11(kk) * 100. / tval                           ! U*WM
      vctr17(kk) = vctr12(kk) * vctr13(kk) / vctr14(kk) * 100. / tval ! UB*WMB
      vctr18(kk) = vctr16(kk) - vctr17(kk)               ! U*WM - UB*WMB = U' * WM'

!!!      print 90, kk,zm(kk),vctr16(kk),vctr17(kk),vctr18(kk)
90    format(8x,i3,4f10.1)

   enddo

! Convert vctr18 from percentage to unit ratio, fill height array in km, 
! and plot profile 

   do kk = 2,mza-1

      vctr18(kk) = vctr18(kk) * .01
   enddo
   
   call plotback()

!----------------------------------------------------------------------
   call oplot_xy2('0','N',aspect,scalelab            &
                 ,mza-2,  vctr18(2),zmkm(2)          &
                 ,'NORMALIZED FLUX','Z [km]'         &
                 ,0.,1.1,.1,5  ,0.,zmkm(mza-1),1.,5  )
!-------------------------------------------------------------------

   call o_frame()

! Close the current workstation if output
! is to a NCAR graphics meta file. This allows viewing the complete
! meta file (including the last frame) during a run and in case the
! simulation crashes.

   if ((trim(runtype) .ne. 'PLOTONLY') .and. (op%plttype .eq. 0)) then
      call o_clswk()
   endif

endif

return
end subroutine uwcomp

