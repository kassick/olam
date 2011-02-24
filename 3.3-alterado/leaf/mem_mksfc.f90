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
Module mem_mksfc

   use leaf_coms, only: maxjml
   use mem_leaf,  only: itab_ul_vars, itab_wl_vars   

!----------------------------------------------------------------------------

   type (itab_ul_vars), allocatable :: mksfc_itab_uls(:)
   type (itab_ul_vars), allocatable :: mksfc_itab_ul(:)
   type (itab_ul_vars), allocatable :: mksfc_itab_us(:)

   type (itab_wl_vars), allocatable :: mksfc_itab_wls(:)
   type (itab_wl_vars), allocatable :: mksfc_itab_wl(:)
   type (itab_wl_vars), allocatable :: mksfc_itab_ws(:)

!----------------------------------------------------------------------------

   Type mksfc_vars
      integer, allocatable :: mpt_sea   (:) ! flag/counter of M pt for sea cell
      integer, allocatable :: mpt_land  (:) ! flag/counter of M pt for land cell
      integer, allocatable :: upt_sea   (:) ! flag/counter of U pt for sea cell
      integer, allocatable :: upt_land  (:) ! flag/counter of U pt for land cell
      integer, allocatable :: wpt_sea   (:) ! Counter of W pt for sea cell
      integer, allocatable :: wpt_land  (:) ! Counter of W pt for land cell

      integer, allocatable :: idatp     (:) ! integer array for storing database data
      integer, allocatable :: leaf_class(:) ! leaf (vegetation) class
      real, allocatable    :: area      (:) ! land/sea cell area [m^2]
      real, allocatable    :: xemls     (:) ! earth x coord of land/sea grid M pts [m]
      real, allocatable    :: yemls     (:) ! earth y coord of land/sea grid M pts [m]
      real, allocatable    :: zemls     (:) ! earth z coord of land/sea grid M pts [m]
      real, allocatable    :: xewls     (:) ! earth x coord of land/sea grid W pts [m]
      real, allocatable    :: yewls     (:) ! earth y coord of land/sea grid W pts [m]
      real, allocatable    :: zewls     (:) ! earth z coord of land/sea grid W pts [m]
   End Type

   Type mksfc_sea_vars
      integer, allocatable :: leaf_class(:) ! leaf (vegetation) class
      real, allocatable    :: area      (:) ! sea cell area [m^2]
      real, allocatable    :: xems      (:) ! earth x coord of sea grid M pts [m]
      real, allocatable    :: yems      (:) ! earth y coord of sea grid M pts [m]
      real, allocatable    :: zems      (:) ! earth z coord of sea grid M pts [m]
      real, allocatable    :: xews      (:) ! earth x coord of sea grid W pts [m]
      real, allocatable    :: yews      (:) ! earth y coord of sea grid W pts [m]
      real, allocatable    :: zews      (:) ! earth z coord of sea grid W pts [m]
   End Type

   Type mksfc_land_vars
      integer, allocatable :: leaf_class  (:) ! leaf (vegetation) class
      integer, allocatable :: ntext_soil(:,:) ! soil textural class
      real, allocatable    :: area        (:) ! land cell area [m^2]
      real, allocatable    :: xeml        (:) ! earth x coord of land grid M pts [m] 
      real, allocatable    :: yeml        (:) ! earth y coord of land grid M pts [m]
      real, allocatable    :: zeml        (:) ! earth z coord of land grid M pts [m]
      real, allocatable    :: xewl        (:) ! earth x coord of land grid W pts [m] 
      real, allocatable    :: yewl        (:) ! earth y coord of land grid W pts [m]
      real, allocatable    :: zewl        (:) ! earth z coord of land grid W pts [m]
   End Type

   type (mksfc_vars), save      :: mksfc
   type (mksfc_sea_vars), save  :: mksfc_sea
   type (mksfc_land_vars), save :: mksfc_land

End Module mem_mksfc
