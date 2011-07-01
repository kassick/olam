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
subroutine olam_mem_alloc()

use mem_basic,   only: alloc_basic, filltab_basic
use mem_cuparm,  only: alloc_cuparm, filltab_cuparm
use mem_micro,   only: alloc_micro, filltab_micro
use mem_radiate, only: alloc_radiate, filltab_radiate
use mem_addsc,   only: addsc, alloc_addsc, filltab_addsc
use mem_tend,    only: alloc_tend, filltab_tend
use mem_turb,    only: alloc_turb, filltab_turb
use mem_sflux,   only: mseaflux, mlandflux, filltab_sflux
use mem_grid,    only: mza, nsw_max, mma, mua, mwa, dts_w
use mem_nudge,   only: nudflag, alloc_nudge2, filltab_nudge
use mem_ijtabs,  only: mrls, filltab_itabs
use aero_coms,   only: alloc_aerosols

use var_tables,  only: vtab_r, vtab_par, scalar_tab, num_var, num_scalar

use misc_coms,   only: io6, naddsc, initial, idiffk, ilwrtyp, iswrtyp,  &
                       nqparm, nqparm_sh, dtsm

use micro_coms,  only: jnmb, ihail, igraup, iaggr, isnow, ipris, irain,  &
                       icloud, level
use max_dims,    only: maxsclr, maxvars

implicit none 

integer :: ng,nv,ntpts

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAM_O_MEM_ALLOC_IN,rst_buf)
#endif

! Allocate variable tables and initialize table counters to zero

allocate (vtab_r(maxvars))
allocate (vtab_par(maxvars))
allocate (scalar_tab(maxsclr))

num_var = 0
num_scalar = 0

! Allocate basic memory needed for 'INITIAL' and 'HISTORY' runs
! Fill variable tables

call filltab_itabs(mma,mua,mwa)  ! Already allocated

call alloc_basic(mza,mua,mwa)
call filltab_basic(mza,mua,mwa)

call alloc_cuparm(mza,mwa,maxval(nqparm(1:mrls)),maxval(nqparm_sh(1:mrls)))
call filltab_cuparm(mza,mwa) 

call alloc_micro(mza,mwa,level,icloud,irain,ipris,isnow,iaggr,igraup,ihail,jnmb) 
call filltab_micro(mza,mwa) 

call alloc_radiate(mza,mwa,ilwrtyp,iswrtyp) 
call filltab_radiate(mza,mwa,ilwrtyp,iswrtyp) 

call alloc_aerosols(mza,mwa,ilwrtyp,iswrtyp) 

call alloc_turb(mza,mwa,nsw_max,idiffk(1))
call filltab_turb(mza,mwa,nsw_max,idiffk(1)) 

call filltab_sflux(mseaflux,mlandflux)  ! Already allocated

if (initial == 2 .and. nudflag > 0) then
  call alloc_nudge2(mza)
  call filltab_nudge(mza)
endif

! Allocate any added Scalar types 

write(io6,*) 'start addsc alloc'

! Allocate length 1 of these datatypes by default

allocate(addsc(1))
if (naddsc > 0) then

! Deallocate datatypes, then re-alloc to correct length

   deallocate(addsc)
   allocate(addsc(naddsc))
   call alloc_addsc(mza,mwa,naddsc)
   call filltab_addsc(mza,mwa,naddsc)
endif

! Allocate tendency data type and fill scalar table.  These subroutine
! calls must occur after allocating all prognostic variables.

write(io6,*) 'start tendency alloc'

call alloc_tend(mza,mua,mwa,naddsc)
call filltab_tend(naddsc) 

! Default:  Fill dts_w with short timestep for MRL = 1

allocate (dts_w(mwa))

dts_w(1:mwa) = dtsm(1)

! Re-fill dts_w with correct value for each mesh refinement level

!do iw = 2,mwa
!   mrl = itab_w(iw)%mrl  !!! later???
!   dts_w(iw) = dtsm(mrl)   ! default

! Special section for CM process at FM timestep

!   if (itab_w(iw)%loop(41)) then  ! Use W52 as identifier
!      dts_w(iw) = dtsm(mrl+1)
!   endif

!enddo

write(io6,*) 'end alloc'
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_O_MEM_ALLOC_OUT,rst_buf)
#endif

return
end subroutine olam_mem_alloc

!===============================================================================

!subroutine dealloc_all()

!use mem_ijtabs
!use mem_basic
!use mem_cuparm
!use mem_grid
!use mem_leaf
!use mem_micro
!use mem_radiate
!use mem_addsc
!use mem_tend
!use mem_turb
!use mem_nudge

!use var_tables

!use misc_coms

! deallocate all model memory.  Used on dynamic balance

!deallocate(vtab_r,scalar_tab)

!call dealloc_basic
!call dealloc_cuparm
!call dealloc_micro
!call dealloc_radiate
!call dealloc_turb

!call dealloc_tend(naddsc)
!call dealloc_addsc(naddsc)

!return
!end
