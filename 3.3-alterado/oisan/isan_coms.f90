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
Module isan_coms

use max_dims, only: maxisfiles, maxpr, maxisdirs

!---------------------------------------------------------------------------
integer :: ioflgisz,ioflgvar,iszstage,ivrstage,iyear,imonth,idate  &
          ,ihour,ipoffset
!---------------------------------------------------------------------------
integer :: npd,lzon_bot,kzonoff,isdirs
real, dimension(maxpr+2) :: pcol_p,pcol_thet,pcol_pk,pcol_u,pcol_v,pcol_z &
   ,pcol_r,pcol_pi,pcol_temp,pcol_rt,pcol_thv

character(len=80) :: innpr
character(len=80), dimension(maxisdirs) :: iapr

integer :: nfgfiles, ifgfile
character(len=128) :: fnames_fg(maxisfiles)
character(len=14) :: ctotdate_fg(maxisfiles)
real(kind=8) :: s1900_fg(maxisfiles)

!---------------------------------------------------------------------------
!     Input pressure file header
!---------------------------------------------------------------------------
integer :: marker,isversion,iyy,imm,idd,ihh,itinc,inproj,ivertcoord
real    :: xnelat,xnelon,cntlat,cntlon,secondlat
!---------------------------------------------------------------------------
integer                   :: nprx,npry,nprz,nprz_rh
integer, dimension(maxpr) :: levpr
real                      :: xswlon,xswlat,gdatdx,gdatdy
real, dimension(maxpr)    :: pnpr
!---------------------------------------------------------------------------
character(len=8), dimension(50)       :: notid
character(len=128)         :: used_file
!---------------------------------------------------------------------------


End module isan_coms
