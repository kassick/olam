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
Module phenology_coms

real, parameter :: retained_carbon_fraction = 0.5  !  Before plants drop their leaves, they retain this fraction of their leaf carbon and nitrogen and put it into storage.

real, parameter :: theta_crit = 0.4  !  When soil porosity (relative to total soil porosity) drops below, this threshold, drought-deciduous plants drop their leaves.

! leaf offset parameters are from White et al. 1997 Global
! Biogeochemical Cycles 11(2) 217-234 
real, parameter :: dl_tr = 655.0     ! critical daylength in minutes
real, parameter :: st_tr1 = 284.3    !  critical soil temp
real, parameter :: st_tr2 = 275.15   ! second critical soil temp

! Phenology parameters for cold deciduous trees
! Botta et al. 2000, Global Change Biology, 6, 709--725
real, parameter :: phen_a = -68.0   
real, parameter :: phen_b = 638.0    
real, parameter :: phen_c = -0.01    
  
end Module phenology_coms
