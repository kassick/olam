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
real function mortality_rates(cc)
  
  use ed_structure_defs
  use pft_coms, only: mort1, mort2, mort3, plant_min_temp, frost_mort
  use disturbance_coms, only: treefall_disturbance_rate,   &
       treefall_hite_threshold

  implicit none

  type(cohort), target :: cc
  real :: mintemp
  real :: threshtemp
  real :: cold_mort

  mortality_rates = 0.0

  !--------Carbon Balance
  mortality_rates = mort1(cc%pft) / (1.0 + exp(mort2(cc%pft) * cc%cbr_bar))

  !--------Treefall
  if(cc%hite <= treefall_hite_threshold)  &
       mortality_rates = mortality_rates + treefall_disturbance_rate

  !--------Frost
  mintemp = cc%patchptr%avg_daily_temp - 273.15
  threshtemp = 5.0 + plant_min_temp(cc%pft)
  if(mintemp < threshtemp)then
     cold_mort = frost_mort
     if(mintemp > plant_min_temp(cc%pft))then
        cold_mort = cold_mort * (1.0 -   &
             (mintemp - plant_min_temp(cc%pft))  &
             /(threshtemp - plant_min_temp(cc%pft)))
     endif
     mortality_rates = mortality_rates + cold_mort
  endif

  !-------Density independent
  mortality_rates = mortality_rates + mort3(cc%pft)

  return
end function mortality_rates
