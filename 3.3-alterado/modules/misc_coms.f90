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
Module misc_coms

use max_dims, only: maxsndg, maxgrds
use rastro_evts

type simtime
   integer year
   integer month
   integer date
   real time
   integer ifirst
end type simtime

type(simtime) :: current_time ! current simulation time

character(len=64) :: expnme
character(len=16) :: runtype
character(len=1)  :: timeunit
character(len=80) :: gridfile
character(len=80) :: hfilin
character(len=80) :: hfilepref
character(len=80) :: zonclim
character(len=80) :: topo_database
                    
integer :: io6
integer :: initial
integer :: lonrad
integer :: ngrids
integer :: nzp
integer :: icorflg
integer :: ilwrtyp
integer :: iswrtyp
integer :: mdomain
integer :: nsndg
integer :: iflag
integer :: iyear1
integer :: imonth1
integer :: idate1
integer :: ihour1
integer :: itime1
integer :: naddsc
integer :: ipsflg
integer :: itsflg
integer :: irtsflg
integer :: iusflg
integer :: ioutput
integer :: iclobber
integer :: itopoflg
integer :: ngrid
integer :: nzpp
integer :: nscl
integer :: nxp
integer :: iparallel

integer :: idiffk (maxgrds)
integer :: ndtrat (maxgrds)
integer :: nacoust(maxgrds)
integer :: nqparm (maxgrds)
integer :: nqparm_sh (maxgrds)

real :: wcldbs
real :: radfrq
real :: confrq
real :: ubmin
real :: dtlong
real :: topref
real :: polelat
real :: polelon
real :: cnum_vel
real :: cnum_sclr
real :: ztop
real :: dzrat
real :: dzmax
real :: frqstate
real :: deltax
real :: deltaz
real :: zbase
real :: p_sfc

real :: centlat(maxgrds)
real :: centlon(maxgrds)
real :: grdlen (maxgrds)
real :: grdwid (maxgrds)
real :: grdaxis(maxgrds)
real :: zkhkm  (maxgrds)
real :: xkhkm  (maxgrds)
real :: cflxy  (maxgrds)
real :: cflz   (maxgrds)
real :: csz    (maxgrds)
real :: csx    (maxgrds)
real :: akmin  (maxgrds)
real :: dtlm   (maxgrds)
real :: dtsm   (maxgrds)

real, save, allocatable :: u01d (:)
real, save, allocatable :: v01d (:)
real, save, allocatable :: pr01d(:)
real, save, allocatable :: th01d(:)
real, save, allocatable :: dn01d(:)
real, save, allocatable :: rt01d(:)

real :: us  (maxsndg)
real :: vs  (maxsndg)
real :: ts  (maxsndg)
real :: thds(maxsndg)
real :: ps  (maxsndg)
real :: hs  (maxsndg)
real :: rts (maxsndg)

real(kind=8) :: time8
real(kind=8) :: time_istp8
real(kind=8) :: timmax8
real(kind=8) :: s1900_init
real(kind=8) :: s1900_sim

Contains

!===============================================================================

   subroutine alloc_misc(mza)

   implicit none
   
   integer, intent(in) :: mza

#ifdef OLAM_RASTRO
call rst_event_i_f(OLAM_ALLOC_MISC_IN,mza)
#endif

   allocate (u01d (mza))
   allocate (v01d (mza))
   allocate (pr01d(mza))
   allocate (th01d(mza))
   allocate (dn01d(mza))
   allocate (rt01d(mza))
  
#ifdef OLAM_RASTRO 
call rst_event_i_f(OLAM_ALLOC_MISC_OUT,mza)
#endif
   return
   end subroutine alloc_misc
   
End Module


