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
Module micro_coms

!-------------------------------------------------------------------------
! This module defines memory for parameters, tables, and other quantities
! that are initialized only once at the beginning of a model run for each
! compute-NODE process.  Afterward, during timesteps, these quantities are
! only read and not written, and may be safely treated as shared memory.
!--------------------------------------------------------------------------

!  The product [(nthz-1)  * dthz ] must equal 25.0
!  The product [(nrhhz-1) * drhhz] must equal 0.18
!  The product [(ntc-1)   * dtc  ] must equal 20.0
!  The product [(ndnc-1)  * ddnc ] must equal 20.e-6

integer, parameter :: maxgrds0 = 20  ! used for homog freezing table - change later

integer, parameter ::  &
    nthz   = 26    & ! # of temp values spanning haze nucleation table
   ,nrhhz  = 10    & ! # of R.H. values spanning haze nucleation table
   ,ngam   = 5000  & ! # of values in incomplete gamma function table
   ,ninc   = 201   & ! # of liquid fraction values spanning shed and melt tables
   ,ndns   = 15    & ! # of diameter values spanning shed table
   ,ntc    = 21    & ! # of temperature values spanning homogeneous freezing table
   ,ndnc   = 11    & ! # of diameter values spanning homogeneous freezing table
   ,nd1cc  = 30    & ! # of diameter values spanning cloud-cloud collection tables
   ,nd1cr  = 15    & ! # of cld diam values spanning cloud-rain collection tables
   ,nr2cr  = 10    & ! # of rain mixr values spanning cloud-rain collection tables
   ,nd2cr  = 30    & ! # of rain diam values spanning cloud-rain collection tables
   ,nr2rr  = 20    & ! # of rain mixr values spanning rain-rain collection tables
   ,nd2rr  = 20    & ! # of rain diam values spanning rain-rain collection tables
   ,ncat   = 7     & ! # of hydrometeor categories without distinguishing ice habits
   ,nhcat  = 15    & ! # of hydrometeor categories including ice habits
   ,npairc = 93    & ! # of pairs of species in number collection table
   ,npairr = 131   & ! # of pairs of species in mass collection table
   ,nembc  = 20    & ! # of diam values spanning main collection table (in 2 dims)
   ,neff   = 10      ! second array dimension for eff (# of types of collection efficiency) 
   
real, parameter ::  & 
    dtc   = 1.     & ! temperature increment in homogeneous nucleation table
   ,ddnc  = 2.e-6  & ! characteristic diam increment in homogeneous nucleation table
   ,dthz  = 1.     & ! temperature increment in haze nucleation table
   ,drhhz = .02      ! R.H. increment in haze nucleation table
   
!--------------------------------------------------------------------------
integer ::  &
    level   &  ! microphysics complexity level read from namelist
   ,icloud  &  ! cloud control flag read from namelist
   ,irain   &  ! rain control flag read from namelist
   ,ipris   &  ! pristine ice control flag read from namelist
   ,isnow   &  ! snow control flag read from namelist
   ,iaggr   &  ! aggregates control flag read from namelist
   ,igraup  &  ! graupel control flag read from namelist
   ,ihail      ! hail control flag read from namelist
   
integer :: mza0  ! micphys copy of number of model vertical levels (excludes top bnd)

integer, dimension(ncat) :: jnmb  ! category control parameter for each category

! table to map lhcat to lcat
integer, dimension(nhcat) :: lcat_lhcat = (/1,2,3,4,5,6,7,3,3,3,3,4,4,4,4/)

! tables to map collision pairs to index in mass and number collection tables
integer, dimension(nhcat,nhcat) :: ipairc,ipairr 

data ipairr/  &
     0,  0,  0,  1,  2,  3,  4,  0,  0,  0,  0,  5,  6,  7,  8,  &
     0,  0,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,  &
     0, 22, 23, 24,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
    25, 26, 27, 28,  0,  0,  0, 29, 30, 31, 32,  0,  0,  0,  0,  &
    33, 34, 35, 36,  0,  0,  0, 37, 38, 39, 40, 41, 42, 43, 44,  &
    45, 46, 47, 48, 49,  0,  0, 50, 51, 52, 53, 54, 55, 56, 57,  &
    58, 59, 60, 61, 62, 63,  0, 64, 65, 66, 67, 68, 69, 70, 71,  &
     0, 72,  0, 73,  0,  0,  0, 74,  0,  0,  0, 75, 76, 77, 78,  &
     0, 79,  0, 80,  0,  0,  0,  0, 81,  0,  0, 82, 83, 84, 85,  &
     0, 86,  0, 87,  0,  0,  0,  0,  0, 88,  0, 89, 90, 91, 92,  &
     0, 93,  0, 94,  0,  0,  0,  0,  0,  0, 95, 96, 97, 98, 99,  &
   100,101,102,  0,  0,  0,  0,103,104,105,106,107,  0,  0,  0,  &
   108,109,110,  0,  0,  0,  0,111,112,113,114,  0,115,  0,  0,  &
   116,117,118,  0,  0,  0,  0,119,120,121,122,  0,  0,123,  0,  &
   124,125,126,  0,  0,  0,  0,127,128,129,130,  0,  0,  0,131   /
   
data ipairc/  &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     0,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     4,  5,  6,  7,  0,  0,  0,  8,  9, 10, 11,  0,  0,  0,  0,  &
    12, 13, 14, 15, 16,  0,  0, 17, 18, 19, 20, 21, 22, 23, 24,  &
    25, 26, 27, 28, 29, 30,  0, 31, 32, 33, 34, 35, 36, 37, 38,  &
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,  &
     0, 54,  0,  0,  0,  0,  0, 55,  0,  0,  0,  0,  0,  0,  0,  &
     0, 56,  0,  0,  0,  0,  0,  0, 57,  0,  0,  0,  0,  0,  0,  &
     0, 58,  0,  0,  0,  0,  0,  0,  0, 59,  0,  0,  0,  0,  0,  &
     0, 60,  0,  0,  0,  0,  0,  0,  0,  0, 61,  0,  0,  0,  0,  &
    62, 63, 64,  0,  0,  0,  0, 65, 66, 67, 68, 69,  0,  0,  0,  &
    70, 71, 72,  0,  0,  0,  0, 73, 74, 75, 76,  0, 77,  0,  0,  &
    78, 79, 80,  0,  0,  0,  0, 81, 82, 83, 84,  0,  0, 85,  0,  &
    86, 87, 88,  0,  0,  0,  0, 89, 90, 91, 92,  0,  0,  0, 93   /

integer, dimension(31,100,2)    :: jhabtab
real ::      &
    cparm    & ! cloud parameter read from namelist
   ,rparm    & ! rain parameter read from namelist
   ,pparm    & ! pristine ice parameter read from namelist (obsolete)
   ,sparm    & ! snow parameter read from namelist
   ,aparm    & ! aggregates parameter read from namelist
   ,gparm    & ! graupel parameter read from namelist
   ,hparm    & ! hail parameter read from namelist
   ,rictmin  & ! minimum diameter-index for collection table lookup
   ,rictmax  & ! maximum diameter-index for collection table lookup
   ,dps      & ! cutoff diameter (microns) between pristine ice and snow
   ,dps2     & ! square of dps
   ,d1min    & ! min cloud droplet diam for cloud-cloud & cloud-rain collect tables
   ,r2min    & ! min rain mixrat for cloud-rain & rain-rain collect tables
   ,d2min    & ! min rain diam for cloud-rain & rain-rain collect tables
   ,d1max    & ! max cloud droplet diam for cloud-cloud & cloud-rain collect tables
   ,r2max    & ! max rain mixrat for cloud-rain & rain-rain collect tables
   ,d2max    & ! max rain diam for cloud-rain & rain-rain collect tables
   ,d1ecc    & ! cloud diameter increment for cloud-cloud collect tables
   ,d1ecr    & ! cloud diameter increment for cloud-rain collect tables
   ,r2ecr    & ! rain mixrat increment for cloud-rain collect tables
   ,r2err    & ! rain mixrat increment for rain-rain collect tables
   ,sedtime0 & ! min generalized timestep considered in sedimentation table
   ,sedtime1   ! max generalized timestep considered in sedimentation table
   
real, dimension(ncat) ::  &
    emb0     & ! minimum mean mass for each category
   ,emb1     & ! maximum mean mass for each category
   ,gnu      & ! gamma distribution width parameter for each category
   ,parm     & ! namelist hydrometeor parameter for each category
   ,emb0log  & ! log of emb0
   ,emb1log  & ! log of emb1
   ,dict     & ! log-mass increment of each category in collection tables
   ,rxmin      ! minimum bulk density for each category
   
real, dimension(nhcat) ::  &
    shapefac  & ! shape factor for each category
   ,cfmas     & ! mass power law coefficient for each category
   ,cfmasi    & ! inverse of cfmas
   ,pwmas     & ! mass power law exponent for each category
   ,pwmasi    & ! inverse of pwmas
   ,cfvt      & ! fall velocity power law coefficient for each category
   ,pwvt      & ! fall velocity power law exponent for each category
   ,emb2      & ! mean hydrometeor mass for each category for case when jnmb(lcat) = 2
   ,dpsmi     & ! mass of 125 micron diameter hydrometeor (used for pris ice & snow)
   ,cfemb0    & ! constant dependent on gamma function and mass power law
   ,cfen0     & ! constant dependent on gamma function and mass power law
   ,pwemb0    & ! constant dependent on gamma function and mass power law
   ,pwen0     & ! constant dependent on gamma function and mass power law
   ,vtfac     & ! constant dependent on gamma function and mass and fall power laws
   ,frefac1   & ! ventilation factor coefficient for constant term
   ,frefac2   & ! ventilation factor coefficient for linear term
   ,cfmasft   & ! constant dependent on gamma function and mass power law
   ,dnfac     & ! constant dependent on gamma function and mass power law 
   ,sipfac    & ! constant dependent on gamma function and mass power law
   ,ch2       & ! sedimentation table increment for each category
   ,ch3       & ! constant dependent on gamma function and mass power law
   ,cdp1      & ! constant dependent on mass and fall power laws
   ,pwvtmasi  & ! constant dependent on mass and fall power laws
   ,dispemb0  & ! minimum vertical displacement in sedim table for each category
   ,dispemb0i & ! inverse of dispemb0
   ,dispemb1    ! maximum vertical displacement in sedim table for each category

real, dimension(nembc,nembc,npairc) :: coltabc ! collection table for bulk number
real, dimension(nembc,nembc,npairr) :: coltabr ! collection table for bulk density

real, dimension(nrhhz,nthz)       :: frachz      ! Haze nucleation table
real, dimension(ndnc,ntc,maxgrds0) :: fracc       ! Homogeneous freezing table
real, dimension(4)                :: gamm,gamn1  ! complete gamma func of nu, nu+1
real, dimension(ngam,3)           :: gam         ! incomplete gamma func
real, dimension(ngam,2)           :: gaminc      ! incomplete gamma func
real, dimension(ngam)             :: gamsip13,gamsip24  ! incomplete gamma func
real, dimension(ninc)             :: rmlttab     ! melting table for bulk density
real, dimension(ninc,nhcat)       :: enmlttab    ! melting table for bulk number
real, dimension(ninc,ndns)        :: shedtab     ! shedding table
real, dimension(2)                :: sc,sk,sl    ! specific heats, latent heats
real, dimension(7)                :: sj          ! flag for rain, graupel, hail

real, dimension(nd1cc)             :: r1tabcc,c1tabcc,c2tabcc ! cloud-cloud coll tabs
real, dimension(nd1cr,nr2cr,nd2cr) :: r1tabcr,c1tabcr         ! cloud-rain coll tabs
real, dimension(nr2rr,nd2rr)       :: c2tabrr                 ! rain-rain coll tabs

! Sedimentation table section

integer, parameter ::  &
    nembfall = 20  & ! # of mass values spanning sedimentation table
   ,maxkfall = 4     ! # of grid levels of fall (in 1 timestep) spanning sedim table

real, allocatable, dimension(:) ::  &
    zmf    & ! micphys copy of z coord at M levels (expanded downward by maxkfall)
   ,dztf   & ! micphys copy of deltaz at T levels (expanded downward by maxkfall)
   ,dzitf    ! inverse of dztf

real, allocatable, dimension(:,:,:,:) ::  &
    pcpfillc  & ! sedim table for bulk number
   ,pcpfillr    ! sedim table for bulk density 


Contains

!===============================================================================

   subroutine alloc_sedimtab(mza)

   implicit none
   integer :: mza
   
   mza0 = mza - 1

   allocate (zmf  (2-maxkfall:mza0))
   allocate (dztf (2-maxkfall:mza0))
   allocate (dzitf(2-maxkfall:mza0))
   allocate (pcpfillc(mza0,maxkfall,nembfall,nhcat))
   allocate (pcpfillr(mza0,maxkfall,nembfall,nhcat))

   return
   end subroutine

End Module micro_coms








