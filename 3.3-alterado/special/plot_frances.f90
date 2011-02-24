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
subroutine plot_frances(iplt)

use mem_ijtabs
use mem_basic
use mem_cuparm
use mem_grid
use mem_leaf
use mem_micro
use mem_radiate
use mem_scalar
use mem_tend
use mem_turb
use mem_nudge

use var_tables

use misc_coms
use oplot_coms
use consts_coms

implicit none

integer :: iplt

integer :: idat,num,nfhour,iwlp,iw
integer, save :: newcall=9
real, save, dimension(0:200) :: flat,flon
real :: zef,ref,xef,yef,xpt,ypt,press_lowest,bsize
character(len=2) :: title

integer, save :: iplottime=0,jplottime
real, save, dimension(1000) :: xmodel,ymodel,xobs,yobs
 
if (iplt /= 1) return

if (newcall /= 1) then
   newcall = 1
   open(32,file='frances_location00',status='old',form='formatted')
   do idat = 0,177
      read(32,21) num,flat(idat),flon(idat)
   enddo
   close(32)
endif
21 format(i4,2f8.2)

nfhour = nint(time8/3600.)

! Find "earth" coordinates of hurricane center

zef = erad * sin(flat(nfhour) * pio180)
ref = erad * cos(flat(nfhour) * pio180)  ! distance from earth center
xef = ref  * cos(flon(nfhour) * pio180)
yef = ref  * sin(flon(nfhour) * pio180)

! Transform hurricane earth coords to whatever projection is in use

call oplot_transform(iplt,xef,yef,zef,xpt,ypt)

! Set character line width

call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
call o_pcsetr('CL',2.)

! Set color

call o_gsplci(10)
call o_gstxci(10)
call o_gsfaci(10)
call o_sflush()

! Plot symbol for hurricane

!call o_pcsetr('CL',7.)
!write(title,'(a1)') '+'

!bsize = .04 * (op%hp2 - op%hp1)

!call o_plchhq (xpt,ypt,trim(adjustl(title)),bsize,0.,0.)
!call o_pcsetr('CL',2.)

!RETURN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Search for hurricane lowest pressure

press_lowest = 2.e5

!----------------------------------------------------------------------
do iw = 1,mwa
!----------------------------------------------------------------------

! Distance of this IW point from search center

   if (glatw(iw) < 15.  .or. glatw(iw) > 35. .or.  &
       glonw(iw) < -90. .or. glonw(iw) > -60.) cycle

   if (lpw(iw) == 2 .and. press(2,iw) < press_lowest) then
      press_lowest = press(2,iw)
      iwlp = iw
   endif
enddo


iplottime = iplottime + 1
xobs(iplottime) = xpt
yobs(iplottime) = ypt

! Find "earth" coordinates of modeled hurricane center

if (iplottime < 18) then
   zef = erad * sin(glatw(iwlp) * pio180)
   ref = erad * cos(glatw(iwlp) * pio180)  ! distance from earth center
   xef = ref  * cos(glonw(iwlp) * pio180)
   yef = ref  * sin(glonw(iwlp) * pio180)
elseif (iplottime == 18) then
   zef = erad * sin(34.4 * pio180)
   ref = erad * cos(34.4 * pio180)  ! distance from earth center
   xef = ref  * cos(-85.7 * pio180)
   yef = ref  * sin(-85.7 * pio180)
elseif (iplottime == 19) then
   zef = erad * sin(36.0 * pio180)
   ref = erad * cos(36.0 * pio180)  ! distance from earth center
   xef = ref  * cos(-86.1 * pio180)
   yef = ref  * sin(-86.1 * pio180)
elseif (iplottime == 20) then
   zef = erad * sin(37.4 * pio180)
   ref = erad * cos(37.4 * pio180)  ! distance from earth center
   xef = ref  * cos(-86.4 * pio180)
   yef = ref  * sin(-86.4 * pio180)
elseif (iplottime == 21) then
   zef = erad * sin(40.0 * pio180)
   ref = erad * cos(40.0 * pio180)  ! distance from earth center
   xef = ref  * cos(-85.2 * pio180)
   yef = ref  * sin(-85.2 * pio180)
endif

! Transform hurricane earth coords to whatever projection is in use

call oplot_transform(iplt,xef,yef,zef,xpt,ypt)

xmodel(iplottime) = xpt
ymodel(iplottime) = ypt

! Plot modeled and observed tracks, and connecting lines

!do jplottime = 1,iplottime
!   call o_frstpt(xmodel(jplottime),ymodel(jplottime))
!   call o_vector(xobs(jplottime),yobs(jplottime))
   
!   if (jplottime > 1) then
!      call o_frstpt(xmodel(jplottime),ymodel(jplottime))
!      call o_vector(xmodel(jplottime-1),ymodel(jplottime-1))
!      call o_frstpt(xobs(jplottime),yobs(jplottime))
!      call o_vector(xobs(jplottime-1),yobs(jplottime-1))
!   endif
!enddo

! Set character size and line width

bsize = .012 * (op%hp2 - op%hp1)
call o_pcsetr('CL',2.)

! Set color

call o_gsplci(10)
call o_gstxci(10)
call o_gsfaci(10)
call o_sflush()

do jplottime = 2,min(18,iplottime)
   write(title,'(i2)') jplottime-1
   call o_plchhq (xmodel(jplottime),ymodel(jplottime),trim(adjustl(title))  &
      ,bsize,0.,0.)
enddo

! Set color

call o_gsplci(215)
call o_gstxci(215)
call o_gsfaci(215)
call o_sflush()

write(title,'(a1)') '+'
call o_plchhq (xmodel(1),ymodel(1),trim(adjustl(title)),bsize,0.,0.)

do jplottime = 2,iplottime
   write(title,'(i2)') jplottime-1
call o_plchhq (xobs(jplottime),yobs(jplottime),trim(adjustl(title)),bsize,0.,0.)
enddo

return
end subroutine plot_frances

