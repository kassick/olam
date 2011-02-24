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
subroutine modsched()

use mem_ijtabs, only: nstp, mrls,  &
                      mrl_endl, mrl_ends, mrl_begl, mrl_begs, leafstep

use misc_coms,  only: io6, nacoust, ndtrat, dtlm, dtlong, dtsm, nqparm
use leaf_coms,  only: dt_leaf, mrl_leaf, isfcl
use sea_coms,   only: dt_sea

implicit none

integer :: ndts
integer :: ng
integer :: j
integer :: i
integer :: ii
integer :: iw
integer :: jstp
integer :: mrl

integer :: nshort(mrls)  ! automatic array
integer :: nlong(mrls)   ! automatic array

! Find number of short timesteps of most refined mesh done per long and short
!    timestep of every refinement level.

nshort(mrls) = 1
nlong(mrls) = nacoust(mrls)

do mrl = mrls-1,1,-1
   nlong(mrl) = nlong(mrl+1) * ndtrat(mrl+1)
   nshort(mrl) = nlong(mrl) / nacoust(mrl)
   ndts = nshort(mrl) / nshort(mrl+1)
   if (nshort(mrl) /= ndts * nshort(mrl+1)) then
      write(io6,'(/,a)')   'Short timesteps not integer ratio between consec MRLs.'
      write(io6,'(a,2i5)') 'mrl, nshort(mrl) = ',mrl,nshort(mrl)
      write(io6,'(a,2i5)') 'mrl+1, nshort(mrl+1) = ',mrl+1,nshort(mrl+1)
      stop 'stop modsched_1'
   endif
enddo

nstp = nlong(1)

do mrl = 1,mrls
   write(io6,*) 'modsched-0 ',mrl,nlong(mrl),nshort(mrl)
enddo

write(io6,'(/,a)') '=== Timestep Schedule ===='
write(io6,'(a,/)') '              jstp    mrl_l  mrl_s'

! Long timestep for each mesh refinement level

dtlm(1) = dtlong
dtsm(1) = dtlm(1) / nacoust(1)

do mrl = 2,mrls
   dtlm(mrl) = dtlm(mrl-1) / float(ndtrat(mrl))
   dtsm(mrl) = dtlm(mrl) / nacoust(mrl)
enddo

! Default:  Fill dt_leaf, mrl_leaf with values for MRL = 1

dt_leaf = dtlm(1)
mrl_leaf = 1
   
! Allocate mrl-schedule arrays

allocate (mrl_begl(nstp))
allocate (mrl_begs(nstp))
allocate (mrl_endl(nstp))
allocate (mrl_ends(nstp))

allocate (leafstep(nstp))

leafstep(1:nstp) = 0
   
! Fill mrl-table values

do jstp = 1,nstp

   mrl_begl(jstp) = 0
   mrl_begs(jstp) = 0
   mrl_endl(jstp) = 0
   mrl_ends(jstp) = 0

! Fill lowest MRL value for which loop processes are to be carried out when
! jstp has its current value.  (However, a zero MRL value indicates process
! is inactive on current jstp.)

   do mrl = mrls,1,-1
      if (mod(jstp-1, nlong(mrl)) == 0) mrl_begl(jstp) = mrl
      if (mod(jstp-1,nshort(mrl)) == 0) mrl_begs(jstp) = mrl

      if (mod(jstp  , nlong(mrl)) == 0) mrl_endl(jstp) = mrl
      if (mod(jstp  ,nshort(mrl)) == 0) mrl_ends(jstp) = mrl
   enddo

   write(io6,333) jstp,mrl_begl(jstp),mrl_begs(jstp),mrl_endl(jstp),mrl_ends(jstp)
   333 format('modsched0 ',5i7)

   if (mrl_ends(jstp) == 0 .or. mrl_begs(jstp) == 0) then
      write(io6,*) 'mrl_s value = 0 is not allowed'
      stop 'stop - modsched_2'
   endif

! Set LEAFSTEP = 1 to do leaf timestep for selected jstp value(s):
!    1. Always do leaf at end of long timestep for mrl = 1
!    2. Also do leaf at end of long timestep for any mrl > 1 whose dtlm > 30 s

   mrl = mrl_endl(jstp)
   if (isfcl == 1 .and. mrl > 0) then
      if (mrl == 1 .or. dtlm(mrl) > 30.) then
         leafstep(jstp) = 1

! Set leaf mrl and timestep according to highest selected mrl

         mrl_leaf = max(mrl_leaf,mrl)
         dt_leaf = min(dt_leaf,dtlm(mrl))
      endif
   endif

enddo

! Set dt_sea equal to dt_leaf

dt_sea = dt_leaf

return
end subroutine modsched
