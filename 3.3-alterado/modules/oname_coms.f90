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
Module oname_coms

   use max_dims,  only: nzgmax, maxsndg, maxgrds, maxisdirs, &
                        maxnplt, maxpltfiles
   use mem_ed,    only: max_soi, max_ed_regions
        
   Type oname_vars

!!    RUNTYPE/NAME

      character(len=64) :: expnme   = 'OLAM test'
      character(len=16) :: runtype  = ''

!!    SIMULATION ENDING TIME

      character(len=1)  :: timeunit = ''
      real(kind=8)      :: timmax   = 0.d0
      
!!    START OF SIMULATION OR ISAN PROCESSING

      integer           :: itime1 = 0
      integer           :: idate1 = 0
      integer           :: imonth1 = 0
      integer           :: iyear1 = 0

!!    GRID SPECIFICATIONS

      integer           :: mdomain = 0
      integer           :: ngrids = 0
      integer           :: nzp = 0
      integer           :: nxp = 0
      real              :: dtlong = 0.0
      real              :: deltax = 0.0
      real              :: deltaz = 0.0
      real              :: dzrat = 1.0
      real              :: dzmax = 2000.0
      real              :: zbase = 0.0
      real, dimension(1003) :: zz = 0.0 ! hardwired dimension

!!    NESTED GRID DEFINITION

      real, dimension(maxgrds) :: centlat = 0.0
      real, dimension(maxgrds) :: centlon = 0.0
      real, dimension(maxgrds) :: grdlen  = 0.0
      real, dimension(maxgrds) :: grdwid  = 0.0
      real, dimension(maxgrds) :: grdaxis = 90.0
      
!!    TIMESTEP RATIOS

      integer, dimension(maxgrds) :: ndtrat  = 1
      integer, dimension(maxgrds) :: nacoust = 1
      
!!    VARIABLE INITIALIZATION INPUT

      integer :: initial = 0
      character(len=80) :: zonclim  = '../etc/ZONAVG_CLIMATE'

!!    NUDGING PARAMETERS

      integer :: nudflag  = 0
      integer :: nudrad   = 2
      real    :: tnudcent = 0.0

!!    GRID, HISTORY FILES

      integer :: ioutput  = 1
      integer :: iclobber = 0
      real    :: frqstate = 3600.0
      character(len=80) :: gridfile  = 'sfcfile/gridfile_0'
      character(len=80) :: hfilin    = ''
      character(len=80) :: hfilepref = 'hist/'

!!    TOPOGRAPHY INITIALIZATION

      integer           :: itopoflg      = 2
      character(len=80) :: topo_database = ''

!!    MODEL/NUMERICAL OPTIONS

      integer :: naddsc    = 0
      integer :: icorflg   = 1
      real    :: cnum_vel  = 0.2
      real    :: cnum_sclr = 0.0

!!    RALEIGH FRICTION PARAMETERS

      real :: rayf_zmin    = 30000.0
      real :: rayf_distim  = 0.0
      real :: rayf_expon   = 1.0
      real :: rayfw_zmin   = 30000.0
      real :: rayfw_distim = 0.0
      real :: rayfw_expon  = 1.0

!!    RADIATION PARAMETERIZATION PARAMETERS

      integer :: iswrtyp = 0
      integer :: ilwrtyp = 0
      real    :: radfrq  = 1800      
      integer :: lonrad  = 1
      
!!    CUMULUS PARAMETERIZATION PARAMETERS

      integer, dimension(maxgrds) :: nqparm = 0
      integer, dimension(maxgrds) :: nqparm_sh = 0
      real                        :: confrq = 1800.0
      real                        :: wcldbs = 0.01

!!    EDDY DIFFUSION PARAMETERS

      integer, dimension(maxgrds) :: idiffk = 2
      real,    dimension(maxgrds) :: zkhkm  = 3.0
      real,    dimension(maxgrds) :: xkhkm  = 3.0
      real,    dimension(maxgrds) :: csx    = 0.2
      real,    dimension(maxgrds) :: csz    = 0.2
      real,    dimension(maxgrds) :: akmin  = 0.0
                        
!!    MICROPHYSICS PARAMETERS

      integer :: level  = 1
      integer :: icloud = 4
      integer :: irain  = 2
      integer :: ipris  = 5
      integer :: isnow  = 2
      integer :: iaggr  = 2
      integer :: igraup = 2
      integer :: ihail  = 2
      real    :: cparm  = 0.3e+9
      real    :: rparm  = 1.0e-3
      real    :: pparm  = 0.0
      real    :: sparm  = 1.0e-2
      real    :: aparm  = 1.0e-2
      real    :: gparm  = 1.0e-2
      real    :: hparm  = 3.0e-2
      
!!    SOUNDING SPECIFICATION

      integer :: nsndg   = 0
      integer :: ipsflg  = 1
      integer :: itsflg  = 0
      integer :: irtsflg = 3
      integer :: iusflg  = 0
      real    :: hs      = 0.0
      real    :: p_sfc   = 1000.0

      real, dimension(5,maxsndg) :: sounding = 0.0

!!    LEAF VARIABLES

      integer :: isfcl    = 1
      integer :: nzg      = 11
      integer :: nzs      = 1
      integer :: ivegflg  = 2
      integer :: isoilflg = 2
      integer :: ndviflg  = 2
      integer :: isstflg  = 2
      integer :: iseaiceflg  = 2

      integer :: isoilstateinit = 0
      integer :: isoildepthflg = 0

      integer :: iupdndvi = 0
      integer :: iupdsst  = 0
      integer :: iupdseaice  = 0
      integer :: nvgcon   = 8
      integer :: nslcon   = 6
      real    :: seatmp   = 280.0

      real, dimension(nzgmax) :: slz = (/ -1.00, -.85, -.70, -.60, -.50, -.40, &
           -.30, -.20, -.15, -.10, -.05, (0.0, i=12,nzgmax) /)
      
      real, dimension(nzgmax) :: slmstr = (/ .35, .35, .35, .35, .35, .35, .35, &
           .35, .35, .35, .35, (0.0, i=12,nzgmax) /)
      
      character(len=80) :: landusefile = 'sfcfiles/landh'
      character(len=80) :: seafile     = 'sfcfiles/seah'

      character(len=80) :: veg_database  = ''
      character(len=80) :: soil_database = ''
      character(len=80) :: ndvi_database = ''
      character(len=80) :: sst_database  = ''
      character(len=80) :: seaice_database  = ''

      character(len=80) :: soilstate_db  = ''
      character(len=80) :: soildepth_db  = ''

!!    ED MODEL VARIABLES
      character(len=128) :: ed_hfilin = ''
      character(len=80) :: ed_inputs_dir = ''
      character(len=80) :: ed_offline_db = ''
      integer           :: n_soi         = 0
      integer           :: n_ed_region   = 0
      integer           :: ied_init_mode = 0
      integer           :: istoma_scheme = 0
      integer           :: iphen_scheme  = 0
      integer           :: n_plant_lim   = 0
      integer           :: n_decomp_lim  = 0
      integer           :: include_fire  = 0
      integer           :: ied_offline   = 0
      integer           :: metcyc1       = 0
      integer           :: metcyc2       = 0
      integer           :: ianth_disturb = 0
      real              :: treefall_disturbance_rate = 0.0
      real              :: runoff_time   = 0.0
      real, dimension(max_soi) :: soi_lat
      real, dimension(max_soi) :: soi_lon
      real, dimension(max_ed_regions) :: ed_reg_latmin
      real, dimension(max_ed_regions) :: ed_reg_latmax
      real, dimension(max_ed_regions) :: ed_reg_lonmin
      real, dimension(max_ed_regions) :: ed_reg_lonmax
      

!!    ISENTROPIC CONTROL

      integer           :: isdirs = 1
      character(len=80) :: iapr(maxisdirs) = ''

!!    MODEL_PLOT VARIABLES

      integer :: nplt       = 0
      integer :: nplt_files = 0
      real    :: frqplt     = 3600.0
      real    :: dtvec      = 1200.0
      real    :: headspeed  = 3.0
      integer :: plttype   = 0
      integer :: pltorient = 0
      character(len=80) :: pltname = 'gmeta'
      integer :: vec_maxmrl = maxgrds
      
!!    THE LIST OF FILES TO PLOT FROM

      character(len=80), dimension(maxpltfiles) :: plt_files = ''

!!    THE ARRAYS OF FIELDS TO PLOT

      character(len=20), dimension(maxnplt) :: fldname ,pltspec1,pltspec2
      integer,           dimension(maxnplt) :: icolortab
      real,              dimension(maxnplt) :: plotcoord1,plotcoord2,slabloc, &
                                               plotwid,viewazim
      character(len=1),  dimension(maxnplt) :: projectn

   End Type   

   type (oname_vars), save :: nl

End Module
