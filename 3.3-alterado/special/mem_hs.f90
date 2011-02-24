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
Module mem_hs

   real, save, allocatable, dimension(:,:,:) :: theta_hs,rho_hs,temp_hs  &
      ,uzonal_hs,umerid_hs

   real, save, allocatable, dimension(:,:)   :: theta_avg_hs,rho_avg_hs  &
      ,temp_avg_hs  &
      ,uzonal_avg_hs,umerid_avg_hs,temp_var_hs,uspectra_hs,count_hs

   real, save, allocatable, dimension(:,:)   :: theta_tavg_hs,temp_tavg_hs  &
      ,uzonal_tavg_hs,umerid_tavg_hs,temp_tvar_hs,uspectrat_hs
                                         
   integer, parameter :: nzhs = 32
                                         
Contains

   subroutine alloc_hs()
   
   print*, 'allocating hs '

   allocate (theta_hs(146,90,nzhs),rho_hs(146,90,nzhs)  &
      ,temp_hs(146,90,nzhs)  &
      ,uzonal_hs(146,90,nzhs),umerid_hs(146,90,nzhs))

   allocate (theta_avg_hs(90,nzhs),rho_avg_hs(90,nzhs)  &
      ,temp_avg_hs(90,nzhs)  &
      ,uzonal_avg_hs(90,nzhs),umerid_avg_hs(90,nzhs),temp_var_hs(90,nzhs))
   allocate (count_hs(144,90))
   allocate (uspectra_hs(16,90))
      
   allocate (theta_tavg_hs(90,nzhs),temp_tavg_hs(90,nzhs)  &
      ,uzonal_tavg_hs(90,nzhs),umerid_tavg_hs(90,nzhs),temp_tvar_hs(90,nzhs))
   allocate (uspectrat_hs(16,90))

   theta_hs (1:146,1:90,1:nzhs) = 0.
   rho_hs   (1:146,1:90,1:nzhs) = 0.
   temp_hs  (1:146,1:90,1:nzhs) = 0.
   uzonal_hs(1:146,1:90,1:nzhs) = 0.
   umerid_hs(1:146,1:90,1:nzhs) = 0.
      
   theta_avg_hs (1:90,1:nzhs) = 0.
   rho_avg_hs   (1:90,1:nzhs) = 0.
   temp_avg_hs  (1:90,1:nzhs) = 0.
   uzonal_avg_hs(1:90,1:nzhs) = 0.
   umerid_avg_hs(1:90,1:nzhs) = 0.
   temp_var_hs  (1:90,1:nzhs) = 0.
      
   count_hs   (1:144,1:90) = 0.
   uspectra_hs(1:16,1:90) = 0.

   theta_tavg_hs (1:90,1:nzhs) = 0.
   temp_tavg_hs  (1:90,1:nzhs) = 0.
   uzonal_tavg_hs(1:90,1:nzhs) = 0.
   umerid_tavg_hs(1:90,1:nzhs) = 0.
   temp_tvar_hs  (1:90,1:nzhs) = 0.
      
   uspectrat_hs(1:16,1:90) = 0.

   return
   end subroutine alloc_hs

End Module mem_hs


   
