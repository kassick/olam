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
real function walltime(wstart)
implicit none

real :: wstart
integer :: ii,ir

call system_clock(count=ii,count_rate=ir)
walltime=float(ii)/(float(ir)/1000.0) - wstart
end function walltime

!===============================================================================

subroutine rams_f_open(iunit, filenm, formt, stat, act, iclob)

! replaces old jclopen and jclget
! files are overwritten unless iclob (ICLOBBER) set to 1

use misc_coms, only: io6
implicit none

integer :: iunit, iclob
character(len=*) :: filenm, formt, stat, act
logical :: exans,opnd

!print*,'filenm,formt,stat1=',filenm,formt,stat

inquire(FILE=filenm,EXIST=exans)

!if(opnd) then
!   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!   print*,'!!!   Trying to open file name :'
!   print*,'!!!       ',filenm
!   print*,'!!!   but it is already opened. Run is ended.'
!   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!   stop 'rams_f_open - opened'
!endif

if(exans.and.iclob.eq.0.and.  &
     (act(1:4).eq.'WRIT'.or.act(1:4).eq.'writ')) then
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   print*,'!!!   trying to open file name :'
   print*,'!!!       ',filenm
   print*,'!!!   but it already exists. run is ended.'
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop 'rams_f_open - exists'
endif

!print*,'filenm,formt,stat2=',filenm(1:len_trim(filenm)),formt,stat
open(iunit,STATUS=stat,FILE=trim(filenm),FORM=formt)
write(io6,*) 'F_open - ',trim(filenm)

end subroutine rams_f_open
