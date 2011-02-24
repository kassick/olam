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
subroutine sfluxes_offline()

  use ed_structure_defs
  use mem_leaf, only: land
  use consts_coms, only: p00i, rocp, cp
  use misc_coms, only: dtlm

  implicit none

  type(site), pointer :: cs
  type(patch), pointer :: ed_patch
  real :: vkmsfc
  real :: sfluxt
  real :: sfluxr
  real :: ustar0
  real :: thetaatm
  real :: thetacan
  real :: hgtmin

  cs => land%first_site
  do while(associated(cs))

     land%ustar(cs%iland) = 0.0

     ed_patch => cs%oldest_patch
     do while(associated(ed_patch))
        hgtmin = max(ed_patch%rough + 1.0e-3,  &
             max(50.0, cs%metinput%geoht))
!             max(ed_patch%veg_height, cs%metinput%geoht))
        call stars(hgtmin,         &
             ed_patch%rough  ,     &
             land%vels(cs%iland),  &
             land%rhos(cs%iland),  &
             cs%metinput%atm_tmp,  &
             cs%metinput%atm_shv,  &
             ed_patch%can_temp,    &
             ed_patch%can_shv ,    &
             vkmsfc,               &
             sfluxt,               &
             sfluxr,               &
             ustar0                )
        
        ed_patch%sxfer_t = dtlm(1) * sfluxt
        ed_patch%sxfer_r = dtlm(1) * sfluxr
        ed_patch%ustar = ustar0

        land%sxfer_t(cs%iland) = land%sxfer_t(cs%iland)  &
                               + dtlm(1) * sfluxt * ed_patch%area
        land%sxfer_r(cs%iland) = land%sxfer_r(cs%iland)  &
                               + dtlm(1) * sfluxr * ed_patch%area

        land%ustar(cs%iland) = land%ustar(cs%iland) + ustar0 * ed_patch%area

        ed_patch => ed_patch%younger
     enddo

     cs => cs%next_site
  enddo

  return
end subroutine sfluxes_offline
