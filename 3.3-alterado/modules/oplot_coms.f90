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
Module oplot_coms

use max_dims, only: maxnplt, maxpltfiles

Type oplot_vars

! Variables copied or computed from $MODEL_PLOT namelist

   integer ::     &
       nplt       & ! Number of fields to plot
      ,nplt_files & ! Number of files to plot from
      ,plttype    & ! Plot output type (ncgm, ps, or pdf)
      ,pltorient  & ! Landscape or portrait
      ,vec_maxmrl   ! Highest mesh level to plot wind vectors

   real ::       &
       frqplt    & ! Time interval between plots
      ,dtvec     & ! Scaling time (seconds) for plotted velocity vectors
      ,headspeed & ! Scaling speed for length of arrow head in plotted velocity vectors
      ,xmin      & ! Frame window minimum x coord in x/y/z or lat/lon/z space
      ,xmax      & ! Frame window maximum x coord in x/y/z or lat/lon/z space
      ,ymin      & ! Frame window minimum y coord in x/y/z or lat/lon/z space
      ,ymax      & ! Frame window maximum y coord in x/y/z or lat/lon/z space
      ,plat3     & ! Plot projection pole latitude
      ,plon3     & ! Plot projection pole longitude
      ,conelat   & ! Plot cone center latitude
      ,conelon   & ! Plot cone center longitude
      ,coneang   & ! Plot cone angle
      ,viewazim    ! Plot view azimuth for vertical cross sections

   character(len=80) :: pltname

   character(len=80), dimension(maxpltfiles) :: plt_files

   character(len=20), dimension(maxnplt) ::  &
       fldname     ! Name of plotted field

   integer, dimension(maxnplt) ::  &
       icolortab   ! Color table number for plotting given field

   character(len=1), dimension(maxnplt) ::  &
       projectn        & ! Plot projection and cross section ['P','E','O','C','X','Y')
      ,contrtyp        & ! Contour type ['T', 'F', 'L']
      ,prtval          & ! Flag to print value ['P', 'N']

      ,pltindx         & ! Print ITABLE index ['I', 'J', 'N']
      ,vectbarb        & ! Plot vectors and/or windbarbs ['B', 'U', 'V']
      ,pltgrid         & ! Plot grid lines ['G', 'N']
      ,pltgrid_landsea & ! Plot land/sea grid lines ['g', 'N']
      ,pltdualgrid     & ! Plot grid lines ['D', 'N']
      ,pltborder       & ! Plot border, border + ticks + labels ['b', 't', 'N']
      ,maptyp          & ! Type of map plot ['M', 'm', 'N']
      ,pltcone         & ! Flag to plot cone circle ['C', 'N']
      ,windowin        & ! Flag to window in ['W', 'N']
      ,frameoff        & ! Flag to suppress frame call ['f', 'N']
      ,panel           & ! Panel number if reduced-size plot ['1', '2', '3', '4']
      ,colorbar        & ! Print colorbar ['c', 'N']
      ,labelbar        & ! Print title, title + info block ['t', 'i', 'N']
      ,pltlev            ! Plot on const press level or near sfc ['p', 's', 'N']
      
   real, dimension(maxnplt) ::  &
       slabloc     ! Z-coord of plot slab in x/y/z or lat/lon/z space

! Plot layout coordinates

   real ::     &
       h1      & ! Frame window minimum x coord in window space (range 0 to 1)
      ,h2      & ! Frame window maximum x coord in window space (range 0 to 1)
      ,v1      & ! Frame window minimum y coord in window space (range 0 to 1)
      ,v2      & ! Frame window maximum y coord in window space (range 0 to 1)
      ,hp1     & ! Panel window minimum x coord in window space (range 0 to 1)
      ,hp2     & ! Panel window maximum x coord in window space (range 0 to 1)
      ,vp1     & ! Panel window minimum y coord in window space (range 0 to 1)
      ,vp2     & ! Panel window maximum y coord in window space (range 0 to 1)
      ,fx1     & ! plot frame left side x coord
      ,fx2     & ! plot frame right side x coord
      ,fy1     & ! plot frame bottom y coord
      ,fy2     & ! plot frame top y coord

      ,cbx1    & ! colorbar left side x coord
      ,cbx2    & ! colorbar right side x coord
      ,cblx    & ! colorbar label left end x coord

      ,fnamey  & ! field name y coord

      ,xlaby   & ! x-axis label y coord
      ,xtlaby  & ! x-axis tick label y coord

      ,ylabx   & ! y-axis label x coord
      ,ytlabx  & ! y-axis tick label right end x coord

      ,timex   & ! elapsed time (sec) left end x coord
      ,timsy   & ! elapsed time (sec) y coord
      ,timdy   & ! elapsed time (day) y coord

      ,slabx   & ! slab left end x coord
      ,slaby   & ! slab y coord
      ,sminy   & ! slab min y coord
      ,smaxy     ! slab max y coord

! Variables set in default subroutine

   integer ::    &
       loopplot  &  ! Flag to plot DO loop indices [1 = YES, 0 = NO]
      ,icigrnd   &  ! "Underground" color index
      ,iplotback    ! Flag indicating whether plotback has been called since
                    !    last frame call [1 = YES, 0 = 'NO']

   character(len=16) ::  &
       dualpts   & ! Plot prognostic or xfer points ['PROG', 'XFER']
      ,relindx     ! Print ITABLE relative indices ['YES', 'NO']

   real ::       &
       psiz      & ! Plot letter/number size
      ,vsprd       ! Plot vertical distance to neighbor indices

! Variables set in oplot_lib subroutine

   integer ::      &
       ifill   ! Flag set to 1 if contrtyp = 'FILL'; 0 otherwise
      
   real ::          &
       fldval_min   & ! Minimum value in slab plot
      ,fldval_max   & ! Maximum value in slab plot
      ,fldvalv_min  & ! Minimum value in vector slab plot
      ,fldvalv_max    ! Maximum value in vector slab plot

   character(len=40) ::  &
       units       & ! Units of plotted field
      ,stagpt      & ! Location in staggered grid of value to be plotted
      ,dimens      & ! Dimensionality of array to be plotted
      ,label         ! Label of field to be plotted
      
! Variables used to save the current graphics attributes

   integer, dimension(14) :: i_att
   real,    dimension(7)  :: r_att

End Type

type (oplot_vars), save :: op

! Miscellaneous

   real, dimension(2) :: xepc,yepc,zepc  ! Earth coords of 2 pts of intersection of
                                         ! plot cone, earth surface, and sides of
                                         ! triangular grid column

End Module


