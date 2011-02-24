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
Module max_dims

   integer, parameter ::    &

       maxgrds      = 20    & ! Max # of grids
      ,maxremote    = 30    & ! Max # of remote send/recv processes
      ,nzgmax       = 25    & ! Max # of soil levels
      ,maxsndg      = 200   & ! Max # of vertical levels for the input sounding
      ,maxvars      = 1000  & ! Max # of model variables 
      ,maxsclr      = 150   & ! Max # of additional scalars
      ,maxsstfiles  = 2000  & ! Max # of input SST files
      ,maxndvifiles = 2000  & ! Max # of input NDVI files
      ,maxisfiles   = 2000  & ! Max # of input data files for the isan stage
      ,maxisdirs    = 30    & ! Max # of directories that contain data files
      ,maxpr        = 100   & ! Max # of press levels in the input data files
      ,maxnplt      = 150   & ! Max # of fields to plot
      ,maxpltfiles  = 2000    ! Max # of input files for a plotonly run
End Module max_dims

