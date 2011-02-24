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
subroutine oplot_lib(iplt,k,i,infotyp,fldname,fldval,notavail,ko)

use mem_ijtabs,  only: itab_w, itab_u, itab_m, itabg_w, itabg_u
use mem_basic,   only: umc, wmc, ump, uc, wc, rho, press, thil, theta, sh_w, sh_v
use mem_cuparm,  only: conprr, aconpr

use mem_grid,    only: dtu, arw0, aru, arw, volt, volui, volwi, xew, yew, zew, &
                       topm, glatw, glonw, lpw, mza, mua, mwa, zm, zt,         &
                       lpu, lcu, xem, yem, zem, xeu, yeu, zeu, dnu, volti,     &
                       unx, uny, unz, utx, uty, utz

use mem_leaf,    only: land

use mem_sea,     only: sea

use mem_micro,   only: sh_c, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h,  &
                       con_c, con_r, con_p, con_s, con_a, con_g, con_h,  &
                       con_ccn, con_ifn,  &
                       pcprr, pcprp, pcprs, pcpra, pcprg, pcprh,  &
                       accpr, accpp, accps, accpa, accpg, accph
                      
use mem_radiate, only: fthrd, rshort, rlong, rlongup, albedt
use mem_addsc,   only: addsc
use mem_tend,    only: umt, wmt
use mem_turb,    only: vkm, vkm_sfc, sflux_w, sflux_t, sflux_r
use mem_nudge,   only: rho_obs, theta_obs, shw_obs, uzonal_obs, umerid_obs,  &
                       rho_sim, theta_sim, shw_sim, uzonal_sim, umerid_sim

use misc_coms,   only: io6, pr01d, dn01d, th01d, time8, iparallel
use oplot_coms,  only: op
use consts_coms, only: p00, rocp, erad, piu180, cp, alvl, grav, omega2
use leaf_coms,   only: slcpd, nzg, slmsts, slz, mwl
use mem_sflux,   only: landflux, seaflux
use sea_coms,    only: mws

!------------------------------------------
! Only for ncar test cases:
!use ncar_testcases_all, only: ncar_testcase, ncar_choice
!------------------------------------------

implicit none

integer, intent(in) :: iplt,k,i
character(len=*), intent(in) :: infotyp,fldname
real, intent(out) :: fldval

integer, intent(out) :: notavail  ! 0 - variable is available
                                  ! 1 - variable is below ground
                                  ! 2 - variable is above model top
                                  ! 3 - variable is not available in this run
integer, intent(out) :: ko  


integer :: klev,kk,nls,ku,kt
real :: uavg,vavg,ue,ve,triu1,triu2,triu3,triv1,triv2,triv3,pressw
real :: vxu1,vyu1,vzu1,vxu2,vyu2,vzu2,vxu3,vyu3,vzu3
real :: vx,vy,vz,raxis,u,v
real :: tempk,fracliq
real :: tuu1,tuu2,tuu3,tuu4,vc,contrib

integer, parameter :: nfields = 291
integer :: nf,iu1,iu2,iu3,iu4,im1,im2,im3,iw1,iw2
integer :: itpn
integer :: iu,iw,kw
integer :: j,ka

real :: pressu1,pressu2,ucc,vcc,wtbotu,wttopu,wtbotw,wttopw,plev,topt
real :: ucc_init, vcc_init, vx_init, vy_init, vz_init, u_init, v_init

real :: theta_bar, theta_bar1

! Stored initial values for perturbation calculations

integer, save :: ncall = 0
real(kind=8), save, allocatable :: press_init(:,:)
real(kind=8), save, allocatable :: rho_init(:,:)
real,         save, allocatable :: theta_init(:,:)
real,         save, allocatable :: uc_init(:,:)

!!! end special

character(len=40), dimension(4,nfields) :: fldlib

!  fldname                        units                 stagpt  dimens
!-----------------------------------------------------------------------------
! ATMOSPHERE - 3D

data fldlib(1:4,  1:31)/  &
  'UMC'           ,'U3' ,'U-NORMAL MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'      & !p  1
 ,'WMC'           ,'W3' ,'W MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'             & !p  2
 ,'UMP'           ,'U3' ,'U-NORMAL MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'      & !p  3
 ,'UC'            ,'U3' ,'U-NORMAL VELOCITY',' (m s:S2:-1  )'                 & !p  4
 ,'WC'            ,'W3' ,'W VELOCITY',' (m s:S2:-1  )'                        & !p  5
 ,'RHO'           ,'T3' ,'AIR DENSITY',' (kg m:S2:-3  )'                      & !p  6
 ,'PRESS'         ,'T3' ,'PRESSURE',' (hPa)'                                  & !   7
 ,'THIL'          ,'T3' ,'ICE-LIQUID THETA',' (K)'                            & !p  8
 ,'THETA'         ,'T3' ,'THETA',' (K)'                                       & !p  9
 ,'AIRTEMPK'      ,'T3' ,'AIR TEMP',' (K)'                                    & !p 10
 ,'AIRTEMPC'      ,'T3' ,'AIR TEMP',' (C)'                                    & !p 11
 ,'SH_W'          ,'T3' ,'TOTAL WATER SPEC DENSITY',' (g kg:S2:-1  )'         & !p 12
 ,'SH_V'          ,'T3' ,'WATER VAPOR SPEC DENSITY',' (g kg:S2:-1  )'         & !p 13
 ,'SH_C'          ,'T3' ,'CLOUDWATER SPEC DENSITY',' (g kg:S2:-1  )'          & !p 14
 ,'SH_R'          ,'T3' ,'RAIN SPEC DENSITY',' (g kg:S2:-1  )'                & !p 15
 ,'SH_P'          ,'T3' ,'PRISTINE ICE SPEC DENSITY',' (g kg:S2:-1  )'        & !p 16
 ,'SH_S'          ,'T3' ,'SNOW SPEC DENSITY',' (g kg:S2:-1  )'                & !p 17
 ,'SH_A'          ,'T3' ,'AGGREGATES SPEC DENSITY',' (g kg:S2:-1  )'          & !p 18
 ,'SH_G'          ,'T3' ,'GRAUPEL SPEC DENSITY',' (g kg:S2:-1  )'             & !p 19
 ,'SH_H'          ,'T3' ,'HAIL SPEC DENSITY',' (g kg:S2:-1  )'                & !p 20
 ,'SH_C_P'        ,'T3' ,'CLOUD + PRIST ICE SPEC DENSITY',' (g kg:S2:-1  )'   & !p 21
 ,'SH_TOTCOND'    ,'T3' ,'CONDENSATE SPEC DENSITY',' (g kg:S2:-1  )'          & !p 22
 ,'CON_C'         ,'T3' ,'CLOUD DROPLET NUMBER CONCEN',' (# mg:S2:-1  )'      & !p 23
 ,'CON_R'         ,'T3' ,'RAIN NUMBER CONCEN',' (# kg:S2:-1  )'               & !p 24
 ,'CON_P'         ,'T3' ,'PRISTINE ICE NUMBER CONCEN',' (# kg:S2:-1  )'       & !p 25
 ,'CON_S'         ,'T3' ,'SNOW NUMBER CONCEN',' (# kg:S2:-1  )'               & !p 26
 ,'CON_A'         ,'T3' ,'AGGREGATES NUMBER CONCEN',' (# kg:S2:-1  )'         & !p 27
 ,'CON_G'         ,'T3' ,'GRAUPEL NUMBER CONCEN',' (# kg:S2:-1  )'            & !p 28
 ,'CON_H'         ,'T3' ,'HAIL NUMBER CONCEN',' (# kg:S2:-1  )'               & !p 29
 ,'CON_CCN'       ,'T3' ,'CCN NUMBER CONCEN',' (# mg:S2:-1  )'                & !p 30
 ,'CON_IFN'       ,'T3' ,'IFN NUMBER CONCEN',' (# kg:S2:-1  )'                / !p 31

data fldlib(1:4, 32:60)/  &
  'VKM'           ,'W3' ,'VERT TURB MOMENTUM K',' (N s m:S2:-2  )'            & !p 32
 ,'FTHRD'         ,'T3' ,'RADIATIVE THETA TENDENCY',' (K s:S2:-1  )'          & !p 33
 ,'UDIFF'         ,'U3' ,'U VELOCITY ERROR',' (m s:S2:-1  )'                  & !  34
 ,'ZONAL_WINDUP'  ,'U3' ,'ZONAL WIND PERT AT U',' (m s:S2:-1  )'              & !  35
 ,'MERID_WINDUP'  ,'U3' ,'MERIDIONAL WIND PERT U',' (m s:S2:-1  )'            & !  36
 ,'RVORTZM'       ,'P3' ,'REL VERT VORTICITY AT M',' (s:S2:-1  )'             & !p 37
 ,'TVORTZM'       ,'P3' ,'TOT VERT VORTICITY AT M',' (s:S2:-1  )'             & !p 38
 ,'DIVERG'        ,'T3' ,'HORIZONTAL DIVERGENCE',' (s:S2:-1  )'               & !p 39
 ,'VC'            ,'U3' ,'U-TANGENTIAL VELOCITY',' (m s:S2:-1  )'             & !p 40
 ,'SPEEDU'        ,'U3' ,'WIND SPEED AT U',' (m s:S2:-1  )'                   & !p 41
 ,'AZIMU'         ,'U3' ,'WIND AZIMUTH AT U',' (deg)'                         & !p 42
 ,'ZONAL_WINDU'   ,'U3' ,'ZONAL WIND AT U',' (m s:S2:-1  )'                   & !p 43
 ,'MERID_WINDU'   ,'U3' ,'MERIDIONAL WIND AT U',' (m s:S2:-1  )'              & !p 44
 ,'RVORTZW'       ,'T3' ,'REL VERT VORTICITY AT W',' (s:S2:-1  )'             & !  45
 ,'TVORTZW'       ,'T3' ,'TOT VERT VORTICITY AT W',' (s:S2:-1  )'             & !  46
 ,'SPEEDW'        ,'T3' ,'WIND SPEED AT W',' (m s:S2:-1  )'                   & !p 47
 ,'AZIMW'         ,'T3' ,'WIND AZIMUTH AT W',' (deg)'                         & !p 48
 ,'ZONAL_WINDW'   ,'T3' ,'ZONAL WIND AT W',' (m s:S2:-1  )'                   & !p 49
 ,'MERID_WINDW'   ,'T3' ,'MERIDIONAL WIND AT W',' (m s:S2:-1  )'              & !p 50
 ,'MASSFLUX'      ,'U3' ,'GRID CELL U-FACE MASS FLUX',' (kg s:S2:-1  )'       & !  51
 ,'PRESS_P'       ,'T3' ,'PRESSURE PERT',' (hPa)'                             & !  52
 ,'RHO_P'         ,'T3' ,'DENSITY PERT',' (kg m:S2:-3  )'                     & !  53
 ,'THETA_P'       ,'T3' ,'THETA PERT',' (K)'                                  & !  54
 ,'VMX'           ,'T3' ,'EARTH-X MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'       & !  55
 ,'VMY'           ,'T3' ,'EARTH-Y MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'       & !  56
 ,'UMT'           ,'U3' ,'U-NORM MOMENTUM TEND',' (kg m:S2:-2   s:S2:-2  )'   & !  57
 ,'WMT'           ,'W3' ,'W MOMENTUM TEND',' (kg m:S2:-2   s:S2:-2  )'        & !  58
 ,'WC100'         ,'W3' ,'W VELOCITY',' (cm s:S2:-1  )'                       & !p 59

! ATMOSPHERE - 2D

 ,'ZPLEV'         ,'T2' ,'HEIGHT OF CONST P SFC',' (m)'                       / !p 60

! ATMOSPHERE SURFACE - 2D

data fldlib(1:4, 61:87)/  &
  'RSHORT'        ,'T2' ,'SFC DOWNWARD SHORTWV FLX',' (W m:S2:-2  )'          & !  61
 ,'RLONG'         ,'T2' ,'SFC DOWNWARD LONGWV FLX',' (W m:S2:-2  )'           & !  62
 ,'RLONGUP'       ,'T2' ,'SFC UPWARD LONGWV FLX',' (W m:S2:-2  )'             & !  63
 ,'ALBEDT'        ,'T2' ,'NET GRID COLUMN SFC ALBEDO',' ( )'                  & !  64
 ,'VKM_SFC'       ,'T2' ,'SFC TURB K FOR MOMENTUM',' (N s m:S2:-2  )'         & !  65
 ,'SFLUX_W'       ,'T2' ,'SFC W MOMENTUM FLUX',' (N m:S2:-2  )'               & !  66
 ,'SENSFLUX'      ,'T2' ,'ATM SFC SENSIBLE HEAT FLUX',' (W m:S2:-2  )'        & !  67
 ,'VAPFLUX'       ,'T2' ,'ATM SFC VAPOR FLUX',' (kg m:S2:-2   s:S2:-1  )'     & !  68
 ,'LATFLUX'       ,'T2' ,'ATM SFC LATENT HEAT FLUX',' (W m:S2:-2  )'          & !  69
 ,'PCPRR'         ,'T2' ,'RAIN PRECIP RATE',' (kg m:S2:-2   h:S2:-1  )'       & !  70
 ,'PCPRP'         ,'T2' ,'PRIST ICE PCP RATE',' (kg m:S2:-2   h:S2:-1  )'     & !  71
 ,'PCPRS'         ,'T2' ,'SNOW PCP RATE',' (kg m:S2:-2   h:S2:-1  )'          & !  72
 ,'PCPRA'         ,'T2' ,'AGGREGATES PCP RATE',' (kg m:S2:-2   h:S2:-1  )'    & !  73
 ,'PCPRG'         ,'T2' ,'GRAUPEL PCP RATE',' (kg m:S2:-2   h:S2:-1  )'       & !  74
 ,'PCPRH'         ,'T2' ,'HAIL PCP RATE',' (kg m:S2:-2   h:S2:-1  )'          & !  75
 ,'PCPRTOT'       ,'T2' ,'MICROPHYS PCP RATE',' (kg m:S2:-2   h:S2:-1  )'     & !  76
 ,'CONPRR'        ,'T2' ,'CONV PCP RATE',' (kg m:S2:-2   h:S2:-1  )'          & !  77
 ,'TOTPRR'        ,'T2' ,'TOTAL PCP RATE',' (kg m:S2:-2   h:S2:-1  )'         & !  78
 ,'ACCPR'         ,'T2' ,'ACCUM RAIN',' (kg m:S2:-2  )'                       & !  79
 ,'ACCPP'         ,'T2' ,'ACCUM PRIST ICE',' (kg m:S2:-2  )'                  & !  80
 ,'ACCPS'         ,'T2' ,'ACCUM SNOW',' (kg m:S2:-2  )'                       & !  81
 ,'ACCPA'         ,'T2' ,'ACCUM AGGREGATES',' (kg m:S2:-2  )'                 & !  82
 ,'ACCPG'         ,'T2' ,'ACCUM GRAUPEL',' (kg m:S2:-2  )'                    & !  83
 ,'ACCPH'         ,'T2' ,'ACCUM HAIL',' (kg m:S2:-2  )'                       & !  84
 ,'ACCPTOT'       ,'T2' ,'ACCUM PRECIP',' (kg m:S2:-2  )'                     & !  85
 ,'ACONPR'        ,'T2' ,'ACCUM CONV PCP',' (kg m:S2:-2  )'                   & !  86
 ,'ATOTPR'        ,'T2' ,'TOTAL ACCUM PCP',' (kg m:S2:-2  )'                  / !  87

! GRID GEOMETRY - 3D

data fldlib(1:4, 88:92)/  &
  'ARU'           ,'U3' ,'AREA OF GRID CELL U-FACE',' (m:S2:2  )'             & !  88
 ,'ARW'           ,'W3' ,'AREA OF GRID CELL W-FACE',' (m:S2:2  )'             & !  89
 ,'VOLT'          ,'T3' ,'GRID T-CELL VOLUME',' (m:S2:3  )'                   & !  90
 ,'VOLU'          ,'U3' ,'GRID U-CELL VOLUME',' (m:S2:3  )'                   & !  91
 ,'VOLW'          ,'W3' ,'GRID W-CELL VOLUME',' (m:S2:3  )'                   / !  92

! GRID GEOMETRY - 2D

data fldlib(1:4, 93:123)/  &
  'TOPM'          ,'M2' ,'TOPOGRAPHY HEIGHT',' (m)'                           & !  93
 ,'GLATW'         ,'T2' ,'LATITUDE',' (deg)'                                  & !  94
 ,'GLONW'         ,'T2' ,'LONGITUDE',' (deg)'                                 & !  95
 ,'LPU'           ,'U2' ,'LOWEST PREDICTED U LEVEL',' ( )'                    & !  96
 ,'LPW'           ,'W2' ,'LOWEST PREDICTED W LEVEL',' ( )'                    & !  97
 ,'XEM'           ,'M2' ,'EARTH-X COORD OF M POINT',' ( )'                    & !  98
 ,'YEM'           ,'M2' ,'EARTH-Y COORD OF M POINT',' ( )'                    & !  99
 ,'ZEM'           ,'M2' ,'EARTH-Z COORD OF M POINT',' ( )'                    & ! 100
 ,'XEU'           ,'U2' ,'EARTH-X COORD OF U POINT',' ( )'                    & ! 101
 ,'YEU'           ,'U2' ,'EARTH-Y COORD OF U POINT',' ( )'                    & ! 102
 ,'ZEU'           ,'U2' ,'EARTH-Z COORD OF U POINT',' ( )'                    & ! 103
 ,'XEW'           ,'W2' ,'EARTH-X COORD OF W POINT',' ( )'                    & ! 104
 ,'YEW'           ,'W2' ,'EARTH-Y COORD OF W POINT',' ( )'                    & ! 105
 ,'ZEW'           ,'W2' ,'EARTH-Z COORD OF W POINT',' ( )'                    & ! 106
 ,'DNU'           ,'U2' ,'DNU',' (m)'                                         & ! 107
 ,'U_DIRU1'       ,'U2' ,'U_DIRU1',' ( )'                                     & ! 108
 ,'U_DIRU2'       ,'U2' ,'U_DIRU2',' ( )'                                     & ! 109
 ,'U_DIRU3'       ,'U2' ,'U_DIRU3',' ( )'                                     & ! 110
 ,'U_DIRU4'       ,'U2' ,'U_DIRU4',' ( )'                                     & ! 111

 ,'FUU5'          ,'U2' ,'FUU5',' ( )'                                        & ! 112
 ,'FUU6'          ,'U2' ,'FUU6',' ( )'                                        & ! 113
 ,'FUU7'          ,'U2' ,'FUU7',' ( )'                                        & ! 114
 ,'FUU8'          ,'U2' ,'FUU8',' ( )'                                        & ! 115
 ,'FUU9'          ,'U2' ,'FUU9',' ( )'                                        & ! 116
 ,'FUU10'         ,'U2' ,'FUU10',' ( )'                                       & ! 117
 ,'FUU11'         ,'U2' ,'FUU11',' ( )'                                       & ! 118
 ,'FUU12'         ,'U2' ,'FUU12',' ( )'                                       & ! 119

 ,'FUW3'          ,'U2' ,'FUW3',' ( )'                                        & ! 120
 ,'FUW4'          ,'U2' ,'FUW4',' ( )'                                        & ! 121
 ,'FUW5'          ,'U2' ,'FUW5',' ( )'                                        & ! 122
 ,'FUW6'          ,'U2' ,'FUW6',' ( )'                                        / ! 123

data fldlib(1:4, 124:152)/  &
  'PGC12'         ,'U2' ,'PGC12',' (m:S2:-1  )'                               & ! 124
 ,'PGC12b'        ,'U2' ,'PGC12b',' (m:S2:-1  )'                              & ! 125
 ,'PGC12c'        ,'U2' ,'PGC12c',' (m:S2:-1  )'                              & ! 126
 ,'PGC12d'        ,'U2' ,'PGC12d',' (m:S2:-1  )'                              & ! 127
 ,'PGC45'         ,'U2' ,'PGC45',' (m:S2:-1  )'                               & ! 128
 ,'PGC45b'        ,'U2' ,'PGC45b',' (m:S2:-1  )'                              & ! 129
 ,'PGC63'         ,'U2' ,'PGC63',' (m:S2:-1  )'                               & ! 130
 ,'PGC63c'        ,'U2' ,'PGC63c',' (m:S2:-1  )'                              & ! 131

 ,'VXU1_U'        ,'U2' ,'VXU1_U',' ( )'                                      & ! 132
 ,'VXU2_U'        ,'U2' ,'VXU2_U',' ( )'                                      & ! 133
 ,'VXU3_U'        ,'U2' ,'VXU3_U',' ( )'                                      & ! 134
 ,'VXU4_U'        ,'U2' ,'VXU4_U',' ( )'                                      & ! 135
 ,'VXW1_U'        ,'U2' ,'VXW1_U',' ( )'                                      & ! 136
 ,'VXW2_U'        ,'U2' ,'VXW2_U',' ( )'                                      & ! 137

 ,'VYU1_U'        ,'U2' ,'VYU1_U',' ( )'                                      & ! 138
 ,'VYU2_U'        ,'U2' ,'VYU2_U',' ( )'                                      & ! 139
 ,'VYU3_U'        ,'U2' ,'VYU3_U',' ( )'                                      & ! 140
 ,'VYU4_U'        ,'U2' ,'VYU4_U',' ( )'                                      & ! 141
 ,'VYW1_U'        ,'U2' ,'VYW1_U',' ( )'                                      & ! 142
 ,'VYW2_U'        ,'U2' ,'VYW2_U',' ( )'                                      & ! 143

 ,'W_DIRU1'       ,'T2' ,'W_DIRU1',' ( )'                                     & ! 144
 ,'W_DIRU2'       ,'T2' ,'W_DIRU2',' ( )'                                     & ! 145
 ,'W_DIRU3'       ,'T2' ,'W_DIRU3',' ( )'                                     & ! 146

 ,'FWU4'          ,'T2' ,'FWU4',' ( )'                                        & ! 147
 ,'FWU5'          ,'T2' ,'FWU5',' ( )'                                        & ! 148
 ,'FWU6'          ,'T2' ,'FWW6',' ( )'                                        & ! 149
 ,'FWU7'          ,'T2' ,'FWU7',' ( )'                                        & ! 150
 ,'FWU8'          ,'T2' ,'FWU8',' ( )'                                        & ! 151
 ,'FWU9'          ,'T2' ,'FWU9',' ( )'                                        / ! 152

data fldlib(1:4, 153:173)/  &
  'FWW1'          ,'T2' ,'FWW1',' ( )'                                        & ! 153
 ,'FWW2'          ,'T2' ,'FWW1',' ( )'                                        & ! 154
 ,'FWW3'          ,'T2' ,'FWW1',' ( )'                                        & ! 155

 ,'VXU1'          ,'T2' ,'VXU1',' ( )'                                        & ! 156
 ,'VXU2'          ,'T2' ,'VXU2',' ( )'                                        & ! 157
 ,'VXU3'          ,'T2' ,'VXU3',' ( )'                                        & ! 158
 ,'VXW'           ,'T2' ,'VXW',' ( )'                                         & ! 159

 ,'VYU1'          ,'T2' ,'VYU1',' ( )'                                        & ! 160
 ,'VYU2'          ,'T2' ,'VYU2',' ( )'                                        & ! 161
 ,'VYU3'          ,'T2' ,'VYU3',' ( )'                                        & ! 162
 ,'VYW'           ,'T2' ,'VYW',' ( )'                                         & ! 163

 ,'VZU1'          ,'T2' ,'VZU1',' ( )'                                        & ! 164
 ,'VZU2'          ,'T2' ,'VZU2',' ( )'                                        & ! 165
 ,'VZU3'          ,'T2' ,'VZU3',' ( )'                                        & ! 166
 ,'VZW'           ,'T2' ,'VZW',' ( )'                                         & ! 167

 ,'VXU1_W'        ,'T2' ,'VXU1_W',' ( )'                                      & ! 168
 ,'VXU2_W'        ,'T2' ,'VXU2_W',' ( )'                                      & ! 169
 ,'VXU3_W'        ,'T2' ,'VXU3_W',' ( )'                                      & ! 170

 ,'VYU1_W'        ,'T2' ,'VYU1_W',' ( )'                                      & ! 171
 ,'VYU2_W'        ,'T2' ,'VYU2_W',' ( )'                                      & ! 172
 ,'VYU3_W'        ,'T2' ,'VYU3_W',' ( )'                                      / ! 173

! GRID INDICES - 2D

data fldlib(1:4,174:195)/  &
  'MRLW'          ,'T2' ,'MRLW',' ( )'                                        & ! 174
 ,'MROW'          ,'T2' ,'MROW',' ( )'                                        & ! 175
 ,'IMGLOBE'       ,'M2' ,'IMGLOBE',' ( )'                                     & ! 176
 ,'IUGLOBE'       ,'U2' ,'IUGLOBE',' ( )'                                     & ! 177
 ,'IWGLOBE'       ,'W2' ,'IWGLOBE',' ( )'                                     & ! 178
 ,'U_IRANK'       ,'U2' ,'U_IRANK',' ( )'                                     & ! 179
 ,'W_IRANK'       ,'W2' ,'W_IRANK',' ( )'                                     & ! 180
 ,'IUP'           ,'U2' ,'IUP',' ( )'                                         & ! 181

 ,'U_IM1'         ,'U2' ,'U_IM1',' ( )'                                       & ! 182
 ,'U_IM2'         ,'U2' ,'U_IM2',' ( )'                                       & ! 183

 ,'U_IU1'         ,'U2' ,'U_IU1',' ( )'                                       & ! 184
 ,'U_IU2'         ,'U2' ,'U_IU2',' ( )'                                       & ! 185
 ,'U_IU3'         ,'U2' ,'U_IU3',' ( )'                                       & ! 186
 ,'U_IU4'         ,'U2' ,'U_IU4',' ( )'                                       & ! 187
 ,'U_IU5'         ,'U2' ,'U_IU5',' ( )'                                       & ! 188
 ,'U_IU6'         ,'U2' ,'U_IU6',' ( )'                                       & ! 189
 ,'U_IU7'         ,'U2' ,'U_IU7',' ( )'                                       & ! 190
 ,'U_IU8'         ,'U2' ,'U_IU8',' ( )'                                       & ! 191
 ,'U_IU9'         ,'U2' ,'U_IU9',' ( )'                                       & ! 192
 ,'U_IU10'        ,'U2' ,'U_IU10',' ( )'                                      & ! 193
 ,'U_IU11'        ,'U2' ,'U_IU11',' ( )'                                      & ! 194
 ,'U_IU12'        ,'U2' ,'U_IU12',' ( )'                                      / ! 195

data fldlib(1:4,196:218)/  &
  'U_IW1'         ,'U2' ,'U_IW1',' ( )'                                       & ! 196
 ,'U_IW2'         ,'U2' ,'U_IW2',' ( )'                                       & ! 197
 ,'U_IW3'         ,'U2' ,'U_IW3',' ( )'                                       & ! 198
 ,'U_IW4'         ,'U2' ,'U_IW4',' ( )'                                       & ! 199
 ,'U_IW5'         ,'U2' ,'U_IW5',' ( )'                                       & ! 200
 ,'U_IW6'         ,'U2' ,'U_IW6',' ( )'                                       & ! 201

 ,'IWP'           ,'W2' ,'IWP',' ( )'                                         & ! 202
 ,'INUDP1'        ,'W2' ,'INUDP1',' ( )'                                      & ! 203

 ,'W_IM1'         ,'W2' ,'W_IM1',' ( )'                                       & ! 204
 ,'W_IM2'         ,'W2' ,'W_IM2',' ( )'                                       & ! 205
 ,'W_IM3'         ,'W2' ,'W_IM3',' ( )'                                       & ! 206

 ,'W_IU1'         ,'W2' ,'W_IU1',' ( )'                                       & ! 207
 ,'W_IU2'         ,'W2' ,'W_IU2',' ( )'                                       & ! 208
 ,'W_IU3'         ,'W2' ,'W_IU3',' ( )'                                       & ! 209
 ,'W_IU4'         ,'W2' ,'W_IU4',' ( )'                                       & ! 210
 ,'W_IU5'         ,'W2' ,'W_IU5',' ( )'                                       & ! 211
 ,'W_IU6'         ,'W2' ,'W_IU6',' ( )'                                       & ! 212
 ,'W_IU7'         ,'W2' ,'W_IU7',' ( )'                                       & ! 213
 ,'W_IU8'         ,'W2' ,'W_IU8',' ( )'                                       & ! 214
 ,'W_IU9'         ,'W2' ,'W_IU9',' ( )'                                       & ! 215

 ,'W_IW1'         ,'W2' ,'W_IW1',' ( )'                                       & ! 216
 ,'W_IW2'         ,'W2' ,'W_IW2',' ( )'                                       & ! 217
 ,'W_IW3'         ,'W2' ,'W_IW3',' ( )'                                       / ! 218

! LAND_CELLS - 3D

data fldlib(1:4,219:242)/  &
  'SOIL_TEXT'     ,'L3G','SOIL TEXTURAL CLASS',' ( )'                         & ! 219
 ,'SOIL_ENERGY'   ,'L3G','SOIL ENERGY',' (J cm:S2:-3  )'                      & ! 220
 ,'SOIL_TEMPC'    ,'L3G','SOIL TEMP',' (C)'                                   & ! 221
 ,'SOIL_FRACLIQ'  ,'L3G','LIQUID FRACTION OF SOIL WATER',' ( )'               & ! 222
 ,'SOIL_WATER'    ,'L3G','SOIL WATER CONTENT',' ( )'                          & ! 223
 ,'SFWAT_MASS'    ,'L3S','SFCWATER MASS',' (kg m:S2:-2  )'                    & ! 224
 ,'SFWAT_ENERGY'  ,'L3S','SFCWATER ENERGY',' (J g:S2:-1  )'                   & ! 225
 ,'SFWAT_TEMPC'   ,'L3S','SFCWATER TEMP',' (C)'                               & ! 226
 ,'SFWAT_FRACLIQ' ,'L3S','SFCWATER LIQUID FRACTION',' ( )'                    & ! 227
 ,'SFWAT_DEPTH'   ,'L3S','SFCWATER DEPTH',' (m)'                              & ! 228

! LAND_CELLS - 2D

 ,'NLEV_SFWAT'    ,'L2' ,'NUMBER OF SFCWATER LAYERS',' ( )'                   & ! 229
 ,'LEAF_CLASS'    ,'L2' ,'LEAF CLASS',' ( )'                                  & ! 230
 ,'VEG_NDVIC'     ,'L2' ,'VEGETATION NDVI',' ( )'                             & ! 231
 ,'VEG_TEMPC'     ,'L2' ,'VEGETATION TEMP',' (C)'                             & ! 232
 ,'VEG_WATER'     ,'L2' ,'VEGETATION SFC WATER ',' (kg m:S2:-2  )'            & ! 233
 ,'STOM_RESIST'   ,'L2' ,'STOMATAL RESISTANCE',' (s m:S2:-1  )'               & ! 234
 ,'GROUND_SHV'    ,'L2' ,'EQUIL VAP SPEC DENSITY OF SOIL',' (g kg:S2:-1  )'   & ! 235
 ,'SOIL_DEPTH'    ,'L2' ,'SOIL DEPTH',' (m)'                                  & ! 236

! SEA_CELLS - 2D

 ,'SEATP'         ,'S2' ,'SEA SFC TEMP (PAST DATA)',' (K)'                    & ! 237
 ,'SEATF'         ,'S2' ,'SEA SFC TEMP (FUTURE DATA)',' (K)'                  & ! 238
 ,'SEATC'         ,'S2' ,'SEA SFC TEMP (CURRENT)',' (K)'                      & ! 239
 ,'SEAICEP'       ,'S2' ,'SEAICE FRACTION (PAST DATA)',' ( )'                 & ! 240
 ,'SEAICEF'       ,'S2' ,'SEAICE FRACTION (FUTURE DATA)',' ( )'               & ! 241
 ,'SEAICEC'       ,'S2' ,'SEAICE FRACTION (CURRENT)',' ( )'                   / ! 242

! LAND AND SEA CELLS - 2D

data fldlib(1:4,243:258)/  &
  'AREA'          ,'B2' ,'LAND/SEA CELL AREA',' (m:S2:2  )'                   & ! 243
 ,'ROUGH'         ,'B2' ,'NET ROUGHNESS HEIGHT',' (m)'                        & ! 244
 ,'CAN_TEMPC'     ,'B2' ,'CANOPY AIR TEMP',' (C)'                             & ! 245
 ,'CAN_SHV'       ,'B2' ,'CANOPY VAPOR SPEC DENSITY',' (g kg:S2:-1  )'        & ! 246
 ,'SFC_TEMPC'     ,'B2' ,'SOIL/SFCWATER/SEA SFC TEMP',' (C)'                  & ! 247
 ,'SFC_SSH'       ,'B2' ,'LAND/SEA SFC SAT VAP SPEC DENS',' (g kg:S2:-1  )'   & ! 248

 ,'SENSFLUX_LS'      ,'B2','L/S CAN TOP SENS HEAT FLX',' (W m:S2:-2  )'       & ! 249
 ,'VAPFLUX_LS'       ,'B2','L/S CAN TOP VAP FLX',' (kg m:S2:-2   s:S2:-1  )'  & ! 250
 ,'LATFLUX_LS'       ,'B2','L/S CAN TOP LAT HEAT FLX',' (W m:S2:-2  )'        & ! 251
 ,'RSHORT_LS'        ,'B2','L/S CAN TOP DOWN SW FLX',' (W m:S2:-2  )'         & ! 252
 ,'RSHORT_DIFFUSE_LS','B2','L/S CAN TOP DOWN DIFFUSE SW FLX',' (W m:S2:-2  )' & ! 253
 ,'RLONG_LS'         ,'B2','L/S CAN TOP DOWN LW FLX',' (W m:S2:-2  )'         & ! 254
 ,'RLONGUP_LS'       ,'B2','L/S CAN TOP UP LW FLX',' (W m:S2:-2  )'           & ! 255
 ,'RLONG_ALBEDO_LS'  ,'B2','L/S NET SFC LW ALBEDO',' ( )'                     & ! 256
 ,'ALBEDO_BEAM_LS'   ,'B2','L/S NET SFC BEAM ALBEDO',' ( )'                   & ! 257
 ,'ALBEDO_DIFFUSE_LS','B2','L/S NET SFC DIFFUSE ALBEDO',' ( )'                / ! 258

! Miscellaneous and new additions

data fldlib(1:4,259:291)/  &
  'RHO_OBS'       ,'T3' ,'NUDGING OBS AIR DENSITY',' (kg m:S2:-3  )'          & ! 259
 ,'THETA_OBS'     ,'T3' ,'NUDGING OBS THETA',' (K)'                           & ! 260
 ,'SHW_OBS'       ,'T3' ,'NUDGING OBS VAPOR SPEC DENSITY',' (g kg:S2:-1  )'   & ! 261
 ,'UZONAL_OBS'    ,'T3' ,'NUDGING OBS ZONAL WIND',' (m s:S2:-1  )'            & ! 262
 ,'UMERID_OBS'    ,'T3' ,'NUDGING OBS MERID WIND',' (m s:S2:-1  )'            & ! 263
 ,'RHO_SIM'       ,'T3' ,'NUDGING SIM AIR DENSITY',' (kg m:S2:-3  )'          & ! 264
 ,'THETA_SIM'     ,'T3' ,'NUDGING SIM THETA',' (K)'                           & ! 265
 ,'SHW_SIM'       ,'T3' ,'NUDGING SIM VAPOR SPEC DENSITY',' (g kg:S2:-1  )'   & ! 266
 ,'UZONAL_SIM'    ,'T3' ,'NUDGING SIM ZONAL WIND',' (m s:S2:-1  )'            & ! 267
 ,'UMERID_SIM'    ,'T3' ,'NUDGING SIM MERID WIND',' (m s:S2:-1  )'            & ! 268
 ,'PRESS_SFC'     ,'T2' ,'SURFACE PRESSURE',' (hPa)'                          & ! 269
 ,'AIRTEMPK_P'    ,'T3' ,'AIR TEMP PERT',' (K)'                               & !p270
 ,'PRESS_SFC_P'   ,'T2' ,'SURFACE PRESSURE PERT',' (hPa)'                     & ! 271
 ,'FCELL_ILSF'    ,'F2' ,'FLUXCELL INDEX',' '                                 & ! 272
 ,'FCELL_IWLS'    ,'F2' ,'FLUXCELL LAND/SEA INDEX',' '                        & ! 273
 ,'FCELL_IW'      ,'F2' ,'FLUXCELL ATM IW INDEX',' '                          & ! 274
 ,'FCELL_KW'      ,'F2' ,'FLUXCELL ATM KW INDEX',' '                          & ! 275
 ,'FCELL_AREA'    ,'F2' ,'FLUXCELL AREA',' (m:S2:2  )'                        & ! 276
 ,'FCELL_ARFATM'  ,'F2' ,'FLUXCELL ATM AREA FRACTION',' '                     & ! 277
 ,'FCELL_ARFLS'   ,'F2' ,'FLUXCELL LAND/SEA AREA FRACTION',' '                & ! 278
 ,'FCELL_SENS'    ,'F2' ,'FLUXCELL SENS HEAT FLUX',' (W m:S2:-2  )'           & ! 279
 ,'FCELL_VAP'     ,'F2' ,'FLUXCELL VAP FLUX',' (kg m:S2:-2   s:S2:-1  )'      & ! 280
 ,'FCELL_LAT'     ,'F2' ,'FLUXCELL LAT HEAT FLUX',' (W m:S2:-2  )'            & ! 281
 ,'FCELL_AIRTEMPC','F2' ,'FLUXCELL AIR TEMP',' (C)'                           & ! 282
 ,'FCELL_AIRTEMPK','F2' ,'FLUXCELL AIR TEMP',' (K)'                           & ! 283
 ,'ADDSC1'        ,'T3' ,'ADDED SCALAR #1 AMOUNT PER KG AIR',' '              & !p284
 ,'ADDSC2'        ,'T3' ,'ADDED SCALAR #2 AMOUNT PER KG AIR',' '              & !p285
 ,'ADDSC3'        ,'T3' ,'ADDED SCALAR #3 AMOUNT PER KG AIR',' '              & !p286
 ,'ADDSC4'        ,'T3' ,'ADDED SCALAR #4 AMOUNT PER KG AIR',' '              & !p287
 ,'ADDSC5'        ,'T3' ,'ADDED SCALAR #5 AMOUNT PER KG AIR',' '              & !p288
 ,'UCP'           ,'U3' ,'NORMAL WIND PERT AT U',' (m s:S2:-1  )'             & !p289
 ,'RVORTZMP'      ,'P3' ,'REL VERT VORTICITY PERT AT M',' (s:S2:-1  )'        & !p290
 ,'MROWH'         ,'T2' ,'MROWH',' ( )'                                       / ! 291

if (ncall /= 10) then
   ncall = 10
   allocate (press_init(mza,mwa))
   allocate (rho_init  (mza,mwa))
   allocate (theta_init(mza,mwa))
   allocate (uc_init   (mza,mua))

   press_init(:,:) = press(:,:)
   rho_init  (:,:) = rho  (:,:)
   theta_init(:,:) = theta(:,:)
   uc_init   (:,:) = uc   (:,:)

endif

if (infotyp == 'UNITS') then

write(io6,*) 'oplib ',iplt,trim(fldname)

   op%label = ' '  ! default null value
   op%units = ' '  ! default null value

   nf = 1
   do while (nf < nfields .and. trim(fldname) /= fldlib(1,nf))
      if (nf == nfields) then
         write(io6,*) 'Plot field ',trim(fldname),' not available in oplot_lib.'
         go to 1000
      endif
      nf = nf + 1
   enddo

   op%stagpt = fldlib(2,nf)(1:1)
   op%dimens = fldlib(2,nf)(2:3)
   op%label  = fldlib(3,nf)
   op%units  = fldlib(4,nf)

   op%fldval_min = 1.e30
   op%fldval_max = -1.e30
   op%fldvalv_min = 1.e30
   op%fldvalv_max = -1.e30

endif

notavail = 0

! Set default vertical index to input value

kt = k
ku = k
ko = k

wtbotw = 1.
wtbotu = 1.

wttopw = 0.
wttopu = 0.

! Reset vertical index for "horizontal" plot if not plotting at constant height

if ((infotyp == 'VALUE' .or. infotyp == 'VALUV') .and.  &
    (op%projectn(iplt) == 'L'  .or.  &
     op%projectn(iplt) == 'P'  .or.  &
     op%projectn(iplt) == 'O')) then

   if (op%pltlev(iplt) == 's') then

      if (op%stagpt == 'U') then
         ku = lpu(i)
      elseif (op%stagpt == 'W' .or. op%stagpt == 'T') then
         kt = lpw(i)
      endif
   
   elseif (op%pltlev(iplt) == 'p') then
     
      plev = op%slabloc(iplt) * 100.  ! convert from mb to Pa

      if (op%stagpt == 'U') then

         iw1 = itab_u(i)%iw1
         iw2 = itab_u(i)%iw2

         do ku = lpu(i)-1,mza-1
            if (.5 * (press(ku+1,iw1) + press(ku+1,iw2)) < plev) exit
         enddo

! Skip this (i,k) point and return with notavail = 1 or 2 if ku is 
! outside allowable range for vertical interpolation

         if (ku < lpu(i)) then
            notavail = 1
            return
         endif
      
         if (ku > mza - 2) then
            notavail = 2
            return
         endif
      
         pressu1 = .5 * (press(ku  ,iw1) + press(ku  ,iw2))
         pressu2 = .5 * (press(ku+1,iw1) + press(ku+1,iw2))

         wtbotu = (plev - pressu2) / (pressu1 - pressu2)    
         wttopu = 1. - wtbotu
   
      elseif (op%stagpt == 'W' .or. op%stagpt == 'T') then

         kt = lpw(i) - 1

         do while (press(kt+1,i) > plev .and. kt < mza-1)
            kt = kt + 1
         enddo

! Skip this (i,k) point and return with notavail = 1 if kt is outside allowable
! range for vertical interpolation

         if (kt < lpw(i)) then
            notavail = 1
            return
         endif

         if (kt > mza - 2) then
            notavail = 2
            return
         endif

         wtbotw = (plev - press(kt+1,i)) / (press(kt,i) - press(kt+1,i))
         wttopw = 1. - wtbotw

      endif
      
   endif
   
endif

! Make check of ku and kt in case not done above.  
! Here, compare ku with lcu instead of lpu.

if (op%stagpt == 'U') then
   ko = ku
   if (op%dimens == '3' .and. ku < lcu(i)) then
      notavail = 1
      return
   endif

   if (op%dimens == '3' .and. ku > mza - 1) then
      notavail = 2
      return
   endif
endif

if (op%stagpt == 'W' .or. op%stagpt == 'T') then
   ko = kt
   if (op%dimens == '3' .and. kt < lpw(i)) then
      notavail = 1
      return
   endif

   if (op%dimens == '3' .and. kt > mza - 1) then
      notavail = 2
      return
   endif
endif

if (op%stagpt == 'L') then
   if (mwl < 2) then
      notavail = 3
      return
   endif
endif

if (op%stagpt == 'S') then
   if (mws < 2) then
      notavail = 3
      return
   endif
endif

! Execute IF block below even when infotyp == 'UNITS'
! in order to check whether current plot field is available in this model run.
! For this type of call to this subroutine, (k,i) are passed in as (1,2).

if     (trim(fldname) == 'UMC'              ) then

   fldval = wtbotu * umc(ku  ,i) &
          + wttopu * umc(ku+1,i)

elseif (trim(fldname) == 'WMC'              ) then

! Need to re-examine use of kt-1 when kt = lpw

   fldval = wtbotw * .5 * (wmc(kt-1,i) + wmc(kt  ,i)) &
          + wttopw * .5 * (wmc(kt  ,i) + wmc(kt+1,i))

elseif (trim(fldname) == 'UMP'              ) then

   fldval = wtbotu * ump(ku  ,i) &
          + wttopu * ump(ku+1,i)

elseif (trim(fldname) == 'UC'               ) then

   fldval = wtbotu * uc(ku  ,i) &
          + wttopu * uc(ku+1,i)

elseif (trim(fldname) == 'UCP'               ) then

   fldval = uc(ku,i) - uc_init(ku,i)

elseif (trim(fldname) == 'VC'               ) then

   iu1 = itab_u(i)%iu1
   iu2 = itab_u(i)%iu2
   iu3 = itab_u(i)%iu3
   iu4 = itab_u(i)%iu4

   tuu1 = itab_u(i)%tuu1
   tuu2 = itab_u(i)%tuu2
   tuu3 = itab_u(i)%tuu3
   tuu4 = itab_u(i)%tuu4

   fldval = wtbotu * (uc(ku  ,iu1) * tuu1  &
                    + uc(ku  ,iu2) * tuu2  &
                    + uc(ku  ,iu3) * tuu3  &
                    + uc(ku  ,iu4) * tuu4) &
          + wttopu * (uc(ku+1,iu1) * tuu1  &
                    + uc(ku+1,iu2) * tuu2  &
                    + uc(ku+1,iu3) * tuu3  &
                    + uc(ku+1,iu4) * tuu4)

elseif (trim(fldname) == 'WC'               ) then

   fldval = wtbotw * .5 * (wc(kt-1,i) + wc(kt  ,i)) &
          + wttopw * .5 * (wc(kt  ,i) + wc(kt+1,i))

elseif (trim(fldname) == 'RHO'              ) then

   fldval = wtbotw * rho(kt  ,i) &
          + wttopw * rho(kt+1,i)

elseif (trim(fldname) == 'PRESS'            ) then

   fldval = press(kt,i) * .01

elseif (trim(fldname) == 'THIL'             ) then

   fldval = wtbotw * thil(kt  ,i) &
          + wttopw * thil(kt+1,i)

elseif (trim(fldname) == 'THETA'            ) then

   fldval = wtbotw * theta(kt  ,i) &
          + wttopw * theta(kt+1,i)

elseif (trim(fldname) == 'AIRTEMPK') then

   fldval = wtbotw * theta(kt  ,i) * (press(kt  ,i) / p00) ** rocp &
          + wttopw * theta(kt+1,i) * (press(kt+1,i) / p00) ** rocp

elseif (trim(fldname) == 'AIRTEMPK_P') then

! SPECIAL - NCAR TEST CASE 6

!   if (ncar_choice == 0) then
!      theta_bar  = 300. * exp(.0001 * zt(kt)   / grav)
!      theta_bar1 = 300. * exp(.0001 * zt(kt+1) / grav)
!   else
!      theta_bar  = 300. * exp(grav * zt(kt)   / (cp * 300.))
!      theta_bar1 = 300. * exp(grav * zt(kt+1) / (cp * 300.))
!   endif
      
!   fldval = wtbotw * (theta(kt  ,i) - theta_bar ) * (press(kt  ,i) / p00) ** rocp &
!          + wttopw * (theta(kt+1,i) - theta_bar1) * (press(kt+1,i) / p00) ** rocp

elseif (trim(fldname) == 'AIRTEMPC'         ) then

   fldval = wtbotw * theta(kt  ,i) * (press(kt  ,i) / p00) ** rocp  &
          + wttopw * theta(kt+1,i) * (press(kt+1,i) / p00) ** rocp - 273.15

elseif (trim(fldname) == 'SH_W'             ) then

   fldval = (wtbotw * sh_w(kt  ,i)  &
          +  wttopw * sh_w(kt+1,i)) * 1.e3

elseif (trim(fldname) == 'SH_V'             ) then

   fldval = (wtbotw * sh_v(kt  ,i)  &
          +  wttopw * sh_v(kt+1,i)) * 1.e3

elseif (trim(fldname) == 'SH_C'             ) then   

   if (.not. allocated(sh_c)) go to 1000

   fldval = (wtbotw * sh_c(kt  ,i)  &
          +  wttopw * sh_c(kt+1,i)) * 1.e3

elseif (trim(fldname) == 'SH_R'             ) then

   if (.not. allocated(sh_r)) go to 1000

   fldval = (wtbotw * sh_r(kt  ,i)  &
          +  wttopw * sh_r(kt+1,i)) * 1.e3

elseif (trim(fldname) == 'SH_P'             ) then

   if (.not. allocated(sh_p)) go to 1000

   fldval = (wtbotw * sh_p(kt  ,i)  &
          +  wttopw * sh_p(kt+1,i)) * 1.e3

elseif (trim(fldname) == 'SH_S'             ) then

   if (.not. allocated(sh_s)) go to 1000

   fldval = (wtbotw * sh_s(kt  ,i)  &
          +  wttopw * sh_s(kt+1,i)) * 1.e3

elseif (trim(fldname) == 'SH_A'             ) then

   if (.not. allocated(sh_a)) go to 1000

   fldval = (wtbotw * sh_a(kt  ,i)  &
          +  wttopw * sh_a(kt+1,i)) * 1.e3

elseif (trim(fldname) == 'SH_G'             ) then

   if (.not. allocated(sh_g)) go to 1000

   fldval = (wtbotw * sh_g(kt  ,i)  &
          +  wttopw * sh_g(kt+1,i)) * 1.e3

elseif (trim(fldname) == 'SH_H'             ) then

   if (.not. allocated(sh_h)) go to 1000

   fldval = (wtbotw * sh_h(kt  ,i)  &
          +  wttopw * sh_h(kt+1,i)) * 1.e3

elseif (trim(fldname) == 'SH_C_P'           ) then

   if (.not. allocated(sh_c)) go to 1000
   if (.not. allocated(sh_p)) go to 1000

   fldval = (wtbotw * (sh_c(kt  ,i) + sh_p(kt  ,i))  &
          +  wttopw * (sh_c(kt+1,i) + sh_p(kt+1,i))) * 1.e3

elseif (trim(fldname) == 'SH_TOTCOND'       ) then

   fldval = (wtbotw * (sh_w(kt  ,i) - sh_v(kt  ,i))  &
          +  wttopw * (sh_w(kt+1,i) - sh_v(kt+1,i))) * 1.e3

elseif (trim(fldname) == 'CON_C'            ) then

   if (.not. allocated(con_c)) go to 1000

   fldval = (wtbotw * con_c(kt  ,i)  &
          +  wttopw * con_c(kt+1,i)) * 1.e-6

elseif (trim(fldname) == 'CON_R'            ) then

   if (.not. allocated(con_r)) go to 1000

   fldval = wtbotw * con_r(kt  ,i) &
          + wttopw * con_r(kt+1,i)

elseif (trim(fldname) == 'CON_P'            ) then

   if (.not. allocated(con_p)) go to 1000

   fldval = wtbotw * con_p(kt  ,i) &
          + wttopw * con_p(kt+1,i)

elseif (trim(fldname) == 'CON_S'            ) then

   if (.not. allocated(con_s)) go to 1000

   fldval = wtbotw * con_s(kt  ,i) &
          + wttopw * con_s(kt+1,i)

elseif (trim(fldname) == 'CON_A'            ) then

   if (.not. allocated(con_a)) go to 1000

   fldval = wtbotw * con_a(kt  ,i) &
          + wttopw * con_a(kt+1,i)

elseif (trim(fldname) == 'CON_G'            ) then

   if (.not. allocated(con_g)) go to 1000

   fldval = wtbotw * con_g(kt  ,i) &
          + wttopw * con_g(kt+1,i)

elseif (trim(fldname) == 'CON_H'            ) then

   if (.not. allocated(con_h)) go to 1000

   fldval = wtbotw * con_h(kt  ,i) &
          + wttopw * con_h(kt+1,i)

elseif (trim(fldname) == 'CON_CCN'          ) then

   if (.not. allocated(con_ccn)) go to 1000

   fldval = (wtbotw * con_ccn(kt  ,i)  &
          +  wttopw * con_ccn(kt+1,i)) * 1.e-6

elseif (trim(fldname) == 'CON_IFN'          ) then

   if (.not. allocated(con_ifn)) go to 1000

   fldval = wtbotw * con_ifn(kt  ,i) &
          + wttopw * con_ifn(kt+1,i)

elseif (trim(fldname) == 'VKM'              ) then

   fldval = wtbotw * vkm(kt  ,i) &
          + wttopw * vkm(kt+1,i)

elseif (trim(fldname) == 'FTHRD'            ) then

   if (.not. allocated(fthrd)) go to 1000

   fldval = wtbotw * fthrd(kt  ,i) &
          + wttopw * fthrd(kt+1,i)

elseif (trim(fldname) == 'ARU'              ) then

   fldval = aru(ku,i)

elseif (trim(fldname) == 'ARW'              ) then

   fldval = arw(kt,i)

elseif (trim(fldname) == 'VOLT'             ) then

   fldval = volt(kt,i)

elseif (trim(fldname) == 'VOLU'             ) then

   fldval = 1. / volui(ku,i)

elseif (trim(fldname) == 'VOLW'             ) then

   fldval = 1. / volwi(kt,i)

elseif (trim(fldname) == 'RVORTZM'  .or.  &
        trim(fldname) == 'RVORTZMP' .or.  &
        trim(fldname) == 'TVORTZM') then

   fldval = 0.
   
   do itpn = 1,itab_m(i)%ntpn
      iu = itab_m(i)%iu(itpn)

      iu1 = itab_u(iu)%iu1
      iu2 = itab_u(iu)%iu2
      iu3 = itab_u(iu)%iu3
      iu4 = itab_u(iu)%iu4

      iw1 = itab_u(iu)%iw1
      iw2 = itab_u(iu)%iw2
      
      tuu1 = itab_u(iu)%tuu1
      tuu2 = itab_u(iu)%tuu2
      tuu3 = itab_u(iu)%tuu3
      tuu4 = itab_u(iu)%tuu4

! Set vertical index for possible interpolation vertically to pressure level

      ku = k

      wtbotu = 1.
      wttopu = 0.

      if ((op%projectn(iplt) == 'L'  .or.  &
           op%projectn(iplt) == 'P'  .or.  &
           op%projectn(iplt) == 'O') .and. &
           op%pltlev(iplt) == 'p') then

         ku = lpu(i) - 1

         do while (.5 * (press(ku+1,iw1) + press(ku+1,iw2)) > plev .and. ku < mza-1)
            ku = ku + 1
         enddo

! Skip this (i,k) point and return with notavail = 1 if ku is outside allowable
! range for vertical interpolation

         if (ku < lpu(iu)) then
            notavail = 1
            return
         endif

         if (ku > mza - 2) then
            notavail = 2
            return
         endif

         pressu1 = .5 * (press(ku  ,iw1) + press(ku  ,iw2))
         pressu2 = .5 * (press(ku+1,iw1) + press(ku+1,iw2))

         wtbotu = (plev - pressu2) / (pressu1 - pressu2)    
         wttopu = 1. - wtbotu

      endif

      if (aru(ku,iu) < 1.e-9) then
         notavail = 1
         return
      endif

      ucc = wtbotu * uc(ku  ,iu) &
          + wttopu * uc(ku+1,iu)

      vcc = wtbotu * (uc(ku  ,iu1) * tuu1  &
                    + uc(ku  ,iu2) * tuu2  &
                    + uc(ku  ,iu3) * tuu3  &
                    + uc(ku  ,iu4) * tuu4) &
          + wttopu * (uc(ku+1,iu1) * tuu1  &
                    + uc(ku+1,iu2) * tuu2  &
                    + uc(ku+1,iu3) * tuu3  &
                    + uc(ku+1,iu4) * tuu4)

      if (trim(fldname) == 'RVORTZMP') then

         ucc_init = wtbotu * uc_init(ku  ,iu) &
                  + wttopu * uc_init(ku+1,iu)

         vcc_init = wtbotu * (uc_init(ku  ,iu1) * tuu1  &
                            + uc_init(ku  ,iu2) * tuu2  &
                            + uc_init(ku  ,iu3) * tuu3  &
                            + uc_init(ku  ,iu4) * tuu4) &
                  + wttopu * (uc_init(ku+1,iu1) * tuu1  &
                            + uc_init(ku+1,iu2) * tuu2  &
                            + uc_init(ku+1,iu3) * tuu3  &
                            + uc_init(ku+1,iu4) * tuu4)

         ucc = ucc - ucc_init
         vcc = vcc - vcc_init

      endif

! Now reconstruct total wind vector projected onto vector from IW1 to IW2

      iw1 = itab_u(iu)%iw1
      iw2 = itab_u(iu)%iw2

      contrib = ucc * (unx(iu) * (xew(iw2) - xew(iw1))   &
                    +  uny(iu) * (yew(iw2) - yew(iw1))   &
                    +  unz(iu) * (zew(iw2) - zew(iw1)))  &
              + vcc * (utx(iu) * (xew(iw2) - xew(iw1))   &
                    +  uty(iu) * (yew(iw2) - yew(iw1))   &
                    +  utz(iu) * (zew(iw2) - zew(iw1)))

      if (i == itab_u(iu)%im1) then
         fldval = fldval + contrib
      else
         fldval = fldval - contrib
      endif

   enddo

   fldval = fldval / itab_m(i)%arm

   if (trim(fldname) == 'TVORTZM') then

      fldval = fldval + omega2 * zem(i) / erad  ! add earth vorticity at M point
   endif

elseif (trim(fldname) == 'RVORTZW' .or.  &
        trim(fldname) == 'TVORTZW'           ) then

   fldval = 0.

   do j = 1,3
      if (j == 1) iu = itab_w(i)%iu1
      if (j == 2) iu = itab_w(i)%iu2
      if (j == 3) iu = itab_w(i)%iu3

      if (k < lpu(iu)) then

         fldval = 0.
         exit
      endif

! Get horizontal tangential wind at U point

      iu1 = itab_u(iu)%iu1
      iu2 = itab_u(iu)%iu2
      iu3 = itab_u(iu)%iu3
      iu4 = itab_u(iu)%iu4

      tuu1 = itab_u(iu)%tuu1
      tuu2 = itab_u(iu)%tuu2
      tuu3 = itab_u(iu)%tuu3
      tuu4 = itab_u(iu)%tuu4

      vc = uc(k,iu1) * tuu1  &
         + uc(k,iu2) * tuu2  &
         + uc(k,iu3) * tuu3  &
         + uc(k,iu4) * tuu4

      if (i == itab_u(iu)%iw2) then
         fldval = fldval + vc * dtu(iu)
      else
         fldval = fldval - vc * dtu(iu)
      endif

   enddo

   fldval = fldval / arw0(i)

   if (trim(fldname) == 'TVORTZW') then
      fldval = fldval + omega2 * zew(i) / erad  ! add earth vorticity at W point
   endif

elseif (trim(fldname) == 'DIVERG'           ) then

   iu1 = itab_w(i)%iu1
   iu2 = itab_w(i)%iu2
   iu3 = itab_w(i)%iu3

   fldval = wtbotw * volti(kt  ,i)  &
                   * (umc(kt  ,iu1) * itab_w(i)%diru1 * aru(kt  ,iu1)  &
                   +  umc(kt  ,iu2) * itab_w(i)%diru2 * aru(kt  ,iu2)  &
                   +  umc(kt  ,iu3) * itab_w(i)%diru3 * aru(kt  ,iu3)) &
          + wttopw * volti(kt+1,i)  &
                   * (umc(kt+1,iu1) * itab_w(i)%diru1 * aru(kt+1,iu1)  &
                   +  umc(kt+1,iu2) * itab_w(i)%diru2 * aru(kt+1,iu2)  &
                   +  umc(kt+1,iu3) * itab_w(i)%diru3 * aru(kt+1,iu3)) 

elseif (trim(fldname) == 'SPEEDW'      .or.  &
        trim(fldname) == 'AZIMW'       .or.  &
        trim(fldname) == 'ZONAL_WINDW' .or.  &
        trim(fldname) == 'MERID_WINDW'       ) then

   iu1 = itab_w(i)%iu1
   iu2 = itab_w(i)%iu2
   iu3 = itab_w(i)%iu3

   vxu1 = itab_w(i)%vxu1; vyu1 = itab_w(i)%vyu1; vzu1 = itab_w(i)%vzu1
   vxu2 = itab_w(i)%vxu2; vyu2 = itab_w(i)%vyu2; vzu2 = itab_w(i)%vzu2
   vxu3 = itab_w(i)%vxu3; vyu3 = itab_w(i)%vyu3; vzu3 = itab_w(i)%vzu3

   vx = wtbotw * (vxu1 * uc(kt  ,iu1)  &
                + vxu2 * uc(kt  ,iu2)  &
                + vxu3 * uc(kt  ,iu3)) &
      + wttopw * (vxu1 * uc(kt+1,iu1)  &
                + vxu2 * uc(kt+1,iu2)  &
                + vxu3 * uc(kt+1,iu3))

   vy = wtbotw * (vyu1 * uc(kt  ,iu1)  &
                + vyu2 * uc(kt  ,iu2)  &
                + vyu3 * uc(kt  ,iu3)) &
      + wttopw * (vyu1 * uc(kt+1,iu1)  &
                + vyu2 * uc(kt+1,iu2)  &
                + vyu3 * uc(kt+1,iu3))

   vz = wtbotw * (vzu1 * uc(kt  ,iu1)  &
                + vzu2 * uc(kt  ,iu2)  &
                + vzu3 * uc(kt  ,iu3)) &
      + wttopw * (vzu1 * uc(kt+1,iu1)  &
                + vzu2 * uc(kt+1,iu2)  &
                + vzu3 * uc(kt+1,iu3))

   if (trim(fldname) == 'SPEEDW'             ) then
      fldval = sqrt(vx ** 2 + vy ** 2 + vz ** 2)
   else
      raxis = sqrt(xew(i) ** 2 + yew(i) ** 2)  ! dist from earth axis
      u = (vy * xew(i) - vx * yew(i)) / raxis
      v = vz * raxis / erad  &
        - (vx * xew(i) + vy * yew(i)) * zew(i) / (raxis * erad) 

      if (trim(fldname) == 'AZIMW'           ) then
         fldval = mod(450. - piu180 * atan2(v,u),360.)
      elseif (trim(fldname) == 'ZONAL_WINDW' ) then
         fldval = u
      elseif (trim(fldname) == 'MERID_WINDW' ) then
         fldval = v
      endif
   endif

elseif (trim(fldname) == 'SPEEDU'      .or.  &
        trim(fldname) == 'AZIMU'       .or.  &
        trim(fldname) == 'ZONAL_WINDU' .or.  &
        trim(fldname) == 'MERID_WINDU'  .or.  &
        trim(fldname) == 'ZONAL_WINDUP' .or.  &
        trim(fldname) == 'MERID_WINDUP'       ) then

   iu1 = itab_u(i)%iu1
   iu2 = itab_u(i)%iu2
   iu3 = itab_u(i)%iu3
   iu4 = itab_u(i)%iu4

   tuu1 = itab_u(i)%tuu1
   tuu2 = itab_u(i)%tuu2
   tuu3 = itab_u(i)%tuu3
   tuu4 = itab_u(i)%tuu4

! Tangential component at U point

   ucc = wtbotu * uc(ku  ,i) &
       + wttopu * uc(ku+1,i)

   vcc = wtbotu * (uc(ku  ,iu1) * tuu1  &
                 + uc(ku  ,iu2) * tuu2  &
                 + uc(ku  ,iu3) * tuu3  &
                 + uc(ku  ,iu4) * tuu4) &
       + wttopu * (uc(ku+1,iu1) * tuu1  &
                 + uc(ku+1,iu2) * tuu2  &
                 + uc(ku+1,iu3) * tuu3  &
                 + uc(ku+1,iu4) * tuu4)

   vx = unx(i) * ucc + utx(i) * vcc
   vy = uny(i) * ucc + uty(i) * vcc
   vz = unz(i) * ucc + utz(i) * vcc

   if (trim(fldname) == 'SPEEDU'             ) then
      fldval = sqrt(vx ** 2 + vy ** 2 + vz ** 2)
   else
      raxis = sqrt(xeu(i) ** 2 + yeu(i) ** 2)  ! dist from earth axis
      u = (vy * xeu(i) - vx * yeu(i)) / raxis
      v = vz * raxis / erad  &
        - (vx * xeu(i) + vy * yeu(i)) * zeu(i) / (raxis * erad) 

      if (trim(fldname) == 'AZIMU'           ) then
         fldval = mod(450. - piu180 * atan2(v,u),360.)
      elseif (trim(fldname) == 'ZONAL_WINDU' ) then
         fldval = u
      elseif (trim(fldname) == 'MERID_WINDU' ) then
         fldval = v

      else

! Tangential component of initial wind at U point

         ucc_init = wtbotu * uc_init(ku  ,i) &
                  + wttopu * uc_init(ku+1,i)

         vcc_init = wtbotu * (uc_init(ku  ,iu1) * tuu1  &
                            + uc_init(ku  ,iu2) * tuu2  &
                            + uc_init(ku  ,iu3) * tuu3  &
                            + uc_init(ku  ,iu4) * tuu4) &
                  + wttopu * (uc_init(ku+1,iu1) * tuu1  &
                            + uc_init(ku+1,iu2) * tuu2  &
                            + uc_init(ku+1,iu3) * tuu3  &
                            + uc_init(ku+1,iu4) * tuu4)

         vx_init = unx(i) * ucc_init + utx(i) * vcc_init
         vy_init = uny(i) * ucc_init + uty(i) * vcc_init
         vz_init = unz(i) * ucc_init + utz(i) * vcc_init

         raxis = sqrt(xeu(i) ** 2 + yeu(i) ** 2)  ! dist from earth axis
         u_init = (vy_init * xeu(i) - vx_init * yeu(i)) / raxis
         v_init = vz_init * raxis / erad  &
            - (vx_init * xeu(i) + vy_init * yeu(i)) * zeu(i) / (raxis * erad) 

         if     (trim(fldname) == 'ZONAL_WINDUP' ) then
            fldval = u - u_init
         elseif (trim(fldname) == 'MERID_WINDUP' ) then
            fldval = v - v_init
         endif

      endif

   endif

elseif (trim(fldname) == 'MASSFLUX'         ) then

   fldval = umc(ku,i) * aru(ku,i)

elseif (trim(fldname) == 'PRESS_P'          ) then

   fldval = press(kt,i) - press_init(kt,i)

elseif (trim(fldname) == 'RHO_P'            ) then

   fldval = wtbotw * (rho(kt  ,i) - rho_init(kt  ,i))  &
          + wttopw * (rho(kt+1,i) - rho_init(kt+1,i))

elseif (trim(fldname) == 'THETA_P'          ) then

   fldval = wtbotw * (theta(kt  ,i) - theta_init(kt  ,i))  &
          + wttopw * (theta(kt+1,i) - theta_init(kt+1,i))

elseif (trim(fldname) == 'VMX'              ) then

   iu1 = itab_w(i)%iu1
   iu2 = itab_w(i)%iu2
   iu3 = itab_w(i)%iu3

   fldval = itab_w(i)%vxu1 * umc(kt,iu1)  &
          + itab_w(i)%vxu2 * umc(kt,iu2)  &
          + itab_w(i)%vxu3 * umc(kt,iu3)  &
          + itab_w(i)%vxw * (wmc(kt-1,i) + wmc(kt,i))

elseif (trim(fldname) == 'VMY'              ) then

   iu1 = itab_w(i)%iu1
   iu2 = itab_w(i)%iu2
   iu3 = itab_w(i)%iu3

   fldval = itab_w(i)%vyu1 * umc(kt,iu1)  &
          + itab_w(i)%vyu2 * umc(kt,iu2)  &
          + itab_w(i)%vyu3 * umc(kt,iu3)  &
          + itab_w(i)%vyw * (wmc(kt-1,i) + wmc(kt,i))

elseif (trim(fldname) == 'UMT'              ) then

   fldval = umt(ku,i) * volui(ku,i)

elseif (trim(fldname) == 'WMT'              ) then

   fldval = wmt(kt,i) * volwi(kt,i)

elseif (trim(fldname) == 'WC100'            ) then

   fldval = (wtbotw * .5 * (wc(kt-1,i) + wc(kt  ,i)) &
          +  wttopw * .5 * (wc(kt  ,i) + wc(kt+1,i))) * 100.

elseif (trim(fldname) == 'TOPM'             ) then

   fldval = topm(i)

elseif (trim(fldname) == 'RSHORT'           ) then

   if (.not. allocated(rshort)) go to 1000

   fldval = rshort(i)

elseif (trim(fldname) == 'RLONG'            ) then

   if (.not. allocated(rlong)) go to 1000

   fldval = rlong(i)

elseif (trim(fldname) == 'RLONGUP'          ) then

   if (.not. allocated(rlongup)) go to 1000

   fldval = rlongup(i)

elseif (trim(fldname) == 'ALBEDT'           ) then

   if (.not. allocated(albedt)) go to 1000

   fldval = albedt(i)

elseif (trim(fldname) == 'GLATW'            ) then

   fldval = glatw(i)

elseif (trim(fldname) == 'GLONW'            ) then

   fldval = glonw(i)

elseif (trim(fldname) == 'VKM_SFC'          ) then

   fldval = vkm_sfc(i)

elseif (trim(fldname) == 'SFLUX_W'          ) then

   fldval = sflux_w(i)

elseif (trim(fldname) == 'SENSFLUX'         ) then

   fldval = sflux_t(i) * cp
   
elseif (trim(fldname) == 'VAPFLUX'          ) then

   fldval = sflux_r(i)

elseif (trim(fldname) == 'LATFLUX'          ) then

   fldval = sflux_r(i) * alvl

elseif (trim(fldname) == 'PCPRR'            ) then

   if (.not. allocated(pcprr)) go to 1000

   fldval = pcprr(i) * 3600.

elseif (trim(fldname) == 'PCPRP'            ) then

   if (.not. allocated(pcprp)) go to 1000

   fldval = pcprp(i)  * 3600.   

elseif (trim(fldname) == 'PCPRS'            ) then

   if (.not. allocated(pcprs)) go to 1000

   fldval = pcprs(i) * 3600.

elseif (trim(fldname) == 'PCPRA'            ) then

   if (.not. allocated(pcpra)) go to 1000

   fldval = pcpra(i) * 3600.

elseif (trim(fldname) == 'PCPRG'            ) then

   if (.not. allocated(pcprg)) go to 1000

   fldval = pcprg(i) * 3600.

elseif (trim(fldname) == 'PCPRH'            ) then

   if (.not. allocated(pcprh)) go to 1000

   fldval = pcprh(i) * 3600.

elseif (trim(fldname) == 'PCPRTOT'          ) then

   fldval = 0.

   if (allocated(pcprr)) fldval = fldval + pcprr(i)
   if (allocated(pcprp)) fldval = fldval + pcprp(i)
   if (allocated(pcprs)) fldval = fldval + pcprs(i)
   if (allocated(pcpra)) fldval = fldval + pcpra(i)
   if (allocated(pcprg)) fldval = fldval + pcprg(i)
   if (allocated(pcprh)) fldval = fldval + pcprh(i)

   fldval = fldval * 3600.

elseif (trim(fldname) == 'CONPRR'           ) then

   if (.not. allocated(conprr)) go to 1000

   fldval = conprr(i) * 3600.

elseif (trim(fldname) == 'TOTPRR'           ) then

   fldval = 0.

   if (allocated(pcprr)) fldval = fldval + pcprr(i)
   if (allocated(pcprp)) fldval = fldval + pcprp(i)
   if (allocated(pcprs)) fldval = fldval + pcprs(i)
   if (allocated(pcpra)) fldval = fldval + pcpra(i)
   if (allocated(pcprg)) fldval = fldval + pcprg(i)
   if (allocated(pcprh)) fldval = fldval + pcprh(i)
   if (allocated(conprr)) fldval = fldval + conprr(i)

   fldval = fldval * 3600.

elseif (trim(fldname) == 'ACCPR'            ) then

   if (.not. allocated(accpr)) go to 1000

   fldval = accpr(i)

elseif (trim(fldname) == 'ACCPP'            ) then

   if (.not. allocated(accpp)) go to 1000

   fldval = accpp(i)

elseif (trim(fldname) == 'ACCPS'            ) then

   if (.not. allocated(accps)) go to 1000

   fldval = accps(i)

elseif (trim(fldname) == 'ACCPA'            ) then

   if (.not. allocated(accpa)) go to 1000

   fldval = accpa(i)

elseif (trim(fldname) == 'ACCPG'            ) then

   if (.not. allocated(accpg)) go to 1000

   fldval = accpg(i)

elseif (trim(fldname) == 'ACCPH'            ) then

   if (.not. allocated(accph)) go to 1000

   fldval = accph(i)

elseif (trim(fldname) == 'ACCPTOT'          ) then

   fldval = 0.

   if (allocated(accpr)) fldval = fldval + accpr(i)
   if (allocated(accpp)) fldval = fldval + accpp(i)
   if (allocated(accps)) fldval = fldval + accps(i)
   if (allocated(accpa)) fldval = fldval + accpa(i)
   if (allocated(accpg)) fldval = fldval + accpg(i)
   if (allocated(accph)) fldval = fldval + accph(i)

elseif (trim(fldname) == 'ACONPR'           ) then

   if (.not. allocated(aconpr)) go to 1000

   fldval = aconpr(i)

elseif (trim(fldname) == 'ATOTPR'           ) then

   fldval = 0.

   if (allocated(accpr)) fldval = fldval + accpr(i)
   if (allocated(accpp)) fldval = fldval + accpp(i)
   if (allocated(accps)) fldval = fldval + accps(i)
   if (allocated(accpa)) fldval = fldval + accpa(i)
   if (allocated(accpg)) fldval = fldval + accpg(i)
   if (allocated(accph)) fldval = fldval + accph(i)
   if (allocated(aconpr)) fldval = fldval + aconpr(i)

elseif (trim(fldname) == 'ZPLEV'          ) then

   fldval = wtbotw * zt(kt  ) &
          + wttopw * zt(kt+1)

elseif (trim(fldname) == 'LPU'              ) then

   fldval = real(lpu(i))

elseif (trim(fldname) == 'LPW'              ) then

   fldval = real(lpw(i))

elseif (trim(fldname) == 'XEM'              ) then

   fldval = xem(i)

elseif (trim(fldname) == 'YEM'              ) then

   fldval = yem(i)

elseif (trim(fldname) == 'ZEM'              ) then

   fldval = zem(i)

elseif (trim(fldname) == 'XEU'              ) then

   fldval = xeu(i)

elseif (trim(fldname) == 'YEU'              ) then

   fldval = yeu(i)

elseif (trim(fldname) == 'ZEU'              ) then

   fldval = zeu(i)

elseif (trim(fldname) == 'XEW'              ) then

   fldval = xew(i)

elseif (trim(fldname) == 'YEW'              ) then

   fldval = yew(i)

elseif (trim(fldname) == 'ZEW'              ) then

   fldval = zew(i)

elseif (trim(fldname) == 'DNU'              ) then

   fldval = dnu(i)

elseif (trim(fldname) == 'U_DIRU1'          ) then

   fldval = itab_u(i)%diru1

elseif (trim(fldname) == 'U_DIRU2'          ) then

   fldval = itab_u(i)%diru2

elseif (trim(fldname) == 'U_DIRU3'          ) then

   fldval = itab_u(i)%diru3

elseif (trim(fldname) == 'U_DIRU4'          ) then

   fldval = itab_u(i)%diru4

elseif (trim(fldname) == 'FUU5'             ) then

   fldval = itab_u(i)%fuu5

elseif (trim(fldname) == 'FUU6'             ) then

   fldval = itab_u(i)%fuu6

elseif (trim(fldname) == 'FUU7'             ) then

   fldval = itab_u(i)%fuu7

elseif (trim(fldname) == 'FUU8'             ) then

   fldval = itab_u(i)%fuu8

elseif (trim(fldname) == 'FUU9'             ) then

   fldval = itab_u(i)%fuu9

elseif (trim(fldname) == 'FUU10'            ) then

   fldval = itab_u(i)%fuu10

elseif (trim(fldname) == 'FUU11'            ) then

   fldval = itab_u(i)%fuu11

elseif (trim(fldname) == 'FUU12'            ) then

   fldval = itab_u(i)%fuu12

elseif (trim(fldname) == 'FUW3'             ) then

   fldval = itab_u(i)%fuw3

elseif (trim(fldname) == 'FUW4'             ) then

   fldval = itab_u(i)%fuw4

elseif (trim(fldname) == 'FUW5'             ) then

   fldval = itab_u(i)%fuw5

elseif (trim(fldname) == 'FUW6'             ) then

   fldval = itab_u(i)%fuw6

elseif (trim(fldname) == 'PGC12'            ) then

   fldval = itab_u(i)%pgc12

elseif (trim(fldname) == 'PGC12b'           ) then

   fldval = itab_u(i)%pgc12b

elseif (trim(fldname) == 'PGC12c'           ) then

   fldval = itab_u(i)%pgc12c

elseif (trim(fldname) == 'PGC12d'           ) then

   fldval = itab_u(i)%pgc12d

elseif (trim(fldname) == 'PGC45'            ) then

   fldval = itab_u(i)%pgc45

elseif (trim(fldname) == 'PGC45b'           ) then

   fldval = itab_u(i)%pgc45b

elseif (trim(fldname) == 'PGC63'            ) then

   fldval = itab_u(i)%pgc63

elseif (trim(fldname) == 'PGC63c'           ) then

   fldval = itab_u(i)%pgc63c

elseif (trim(fldname) == 'VXU1_U'           ) then

   fldval = itab_u(i)%vxu1_u

elseif (trim(fldname) == 'VXU2_U'           ) then

   fldval = itab_u(i)%vxu2_u

elseif (trim(fldname) == 'VXU3_U'           ) then

   fldval = itab_u(i)%vxu3_u

elseif (trim(fldname) == 'VXU4_U'           ) then

   fldval = itab_u(i)%vxu4_u

elseif (trim(fldname) == 'VXW1_U'           ) then

   fldval = itab_u(i)%vxw1_u

elseif (trim(fldname) == 'VXW2_U'           ) then

   fldval = itab_u(i)%vxw2_u

elseif (trim(fldname) == 'VYU1_U'           ) then

   fldval = itab_u(i)%vyu1_u

elseif (trim(fldname) == 'VYU2_U'           ) then

   fldval = itab_u(i)%vyu2_u

elseif (trim(fldname) == 'VYU3_U'           ) then

   fldval = itab_u(i)%vyu3_u

elseif (trim(fldname) == 'VYU4_U'           ) then

   fldval = itab_u(i)%vyu4_u

elseif (trim(fldname) == 'VYW1_U'           ) then

   fldval = itab_u(i)%vyw1_u

elseif (trim(fldname) == 'VYW2_U'           ) then

   fldval = itab_u(i)%vyw2_u

elseif (trim(fldname) == 'W_DIRU1'          ) then

   fldval = itab_w(i)%diru1

elseif (trim(fldname) == 'W_DIRU2'          ) then

   fldval = itab_w(i)%diru2

elseif (trim(fldname) == 'W_DIRU3'          ) then

   fldval = itab_w(i)%diru3

elseif (trim(fldname) == 'FWU4'             ) then

   fldval = itab_w(i)%fwu4

elseif (trim(fldname) == 'FWU5'             ) then

   fldval = itab_w(i)%fwu5

elseif (trim(fldname) == 'FWU6'             ) then

   fldval = itab_w(i)%fwu6

elseif (trim(fldname) == 'FWU7'             ) then

   fldval = itab_w(i)%fwu7

elseif (trim(fldname) == 'FWU8'             ) then

   fldval = itab_w(i)%fwu8

elseif (trim(fldname) == 'FWU9'             ) then

   fldval = itab_w(i)%fwu9

elseif (trim(fldname) == 'FWW1'             ) then

   fldval = itab_w(i)%fww1

elseif (trim(fldname) == 'FWW2'             ) then

   fldval = itab_w(i)%fww2

elseif (trim(fldname) == 'FWW3'             ) then

   fldval = itab_w(i)%fww3

elseif (trim(fldname) == 'VXU1'             ) then

   fldval = itab_w(i)%vxu1

elseif (trim(fldname) == 'VXU2'             ) then

   fldval = itab_w(i)%vxu2

elseif (trim(fldname) == 'VXU3'             ) then

   fldval = itab_w(i)%vxu3

elseif (trim(fldname) == 'VXW'              ) then

   fldval = itab_w(i)%vxw

elseif (trim(fldname) == 'VYU1'             ) then

   fldval = itab_w(i)%vyu1

elseif (trim(fldname) == 'VYU2'             ) then

   fldval = itab_w(i)%vyu2

elseif (trim(fldname) == 'VYU3'             ) then

   fldval = itab_w(i)%vyu3

elseif (trim(fldname) == 'VYW'              ) then

   fldval = itab_w(i)%vyw

elseif (trim(fldname) == 'VZU1'             ) then

   fldval = itab_w(i)%vzu1

elseif (trim(fldname) == 'VZU2'             ) then

   fldval = itab_w(i)%vzu2

elseif (trim(fldname) == 'VZU3'             ) then

   fldval = itab_w(i)%vzu3

elseif (trim(fldname) == 'VZW'              ) then

   fldval = itab_w(i)%vzw

elseif (trim(fldname) == 'VXU1_W'           ) then

   fldval = itab_w(i)%vxu1_w

elseif (trim(fldname) == 'VXU2_W'           ) then

   fldval = itab_w(i)%vxu2_w

elseif (trim(fldname) == 'VXU3_W'           ) then

   fldval = itab_w(i)%vxu3_w

elseif (trim(fldname) == 'VYU1_W'           ) then

   fldval = itab_w(i)%vyu1_w

elseif (trim(fldname) == 'VYU2_W'           ) then

   fldval = itab_w(i)%vyu2_w

elseif (trim(fldname) == 'VYU3_W'           ) then

   fldval = itab_w(i)%vyu3_w

elseif (trim(fldname) == 'IUP'              ) then

   fldval = itab_u(i)%iup

elseif (trim(fldname) == 'U_IM1'            ) then

   fldval = itab_u(i)%im1

elseif (trim(fldname) == 'U_IM2'            ) then

   fldval = itab_u(i)%im2

elseif (trim(fldname) == 'U_IU1'            ) then

   fldval = itab_u(i)%iu1

elseif (trim(fldname) == 'U_IU2'            ) then

   fldval = itab_u(i)%iu2

elseif (trim(fldname) == 'U_IU3'            ) then

   fldval = itab_u(i)%iu3

elseif (trim(fldname) == 'U_IU4'            ) then

   fldval = itab_u(i)%iu4

elseif (trim(fldname) == 'U_IU5'            ) then

   fldval = itab_u(i)%iu5

elseif (trim(fldname) == 'U_IU6'            ) then

   fldval = itab_u(i)%iu6

elseif (trim(fldname) == 'U_IU7'            ) then

   fldval = itab_u(i)%iu7

elseif (trim(fldname) == 'U_IU8'            ) then

   fldval = itab_u(i)%iu8

elseif (trim(fldname) == 'U_IU9'            ) then

   fldval = itab_u(i)%iu9

elseif (trim(fldname) == 'U_IU10'           ) then

   fldval = itab_u(i)%iu10

elseif (trim(fldname) == 'U_IU11'           ) then

   fldval = itab_u(i)%iu11

elseif (trim(fldname) == 'U_IU12'           ) then

   fldval = itab_u(i)%iu12

elseif (trim(fldname) == 'U_IW1'            ) then

   fldval = itab_u(i)%iw1

elseif (trim(fldname) == 'U_IW2'            ) then

   fldval = itab_u(i)%iw2

elseif (trim(fldname) == 'U_IW3'            ) then

   fldval = itab_u(i)%iw3

elseif (trim(fldname) == 'U_IW4'            ) then

   fldval = itab_u(i)%iw4

elseif (trim(fldname) == 'U_IW5'            ) then

   fldval = itab_u(i)%iw5

elseif (trim(fldname) == 'U_IW6'            ) then

   fldval = itab_u(i)%iw6

elseif (trim(fldname) == 'IWP'              ) then

   fldval = itab_w(i)%iwp

elseif (trim(fldname) == 'INUDP1'           ) then

   fldval = itab_w(i)%inudp(1)

elseif (trim(fldname) == 'W_IM1'            ) then

   fldval = itab_w(i)%im1

elseif (trim(fldname) == 'W_IM2'            ) then

   fldval = itab_w(i)%im2

elseif (trim(fldname) == 'W_IM3'            ) then

   fldval = itab_w(i)%im3

elseif (trim(fldname) == 'W_IU1'            ) then

   fldval = itab_w(i)%iu1

elseif (trim(fldname) == 'W_IU2'            ) then

   fldval = itab_w(i)%iu2

elseif (trim(fldname) == 'W_IU3'            ) then

   fldval = itab_w(i)%iu3

elseif (trim(fldname) == 'W_IU4'            ) then

   fldval = itab_w(i)%iu4

elseif (trim(fldname) == 'W_IU5'            ) then

   fldval = itab_w(i)%iu5

elseif (trim(fldname) == 'W_IU6'            ) then

   fldval = itab_w(i)%iu6

elseif (trim(fldname) == 'W_IU7'            ) then

   fldval = itab_w(i)%iu7

elseif (trim(fldname) == 'W_IU8'            ) then

   fldval = itab_w(i)%iu8

elseif (trim(fldname) == 'W_IU9'            ) then

   fldval = itab_w(i)%iu9

elseif (trim(fldname) == 'W_IW1'            ) then

   fldval = itab_w(i)%iw1

elseif (trim(fldname) == 'W_IW2'            ) then

   fldval = itab_w(i)%iw2

elseif (trim(fldname) == 'W_IW3'            ) then

   fldval = itab_w(i)%iw3

elseif (trim(fldname) == 'SOIL_TEXT'        ) then

   fldval = real(land%ntext_soil(k,i))

elseif (trim(fldname) == 'SOIL_ENERGY'      ) then

   fldval = land%soil_energy(k,i) * 1.e-6

elseif (trim(fldname) == 'SOIL_TEMPC'       ) then

   call qwtk(land%soil_energy(k,i)        &
            ,land%soil_water(k,i)*1.e3    &
            ,slcpd(land%ntext_soil(k,i))  &
            ,tempk,fracliq)
   fldval = tempk - 273.15

elseif (trim(fldname) == 'SOIL_FRACLIQ'     ) then

   call qwtk(land%soil_energy(k,i)        &
            ,land%soil_water(k,i)*1.e3    &
            ,slcpd(land%ntext_soil(k,i))  &
            ,tempk,fracliq)
   fldval = fracliq

elseif (trim(fldname) == 'SOIL_WATER'       ) then

   fldval = land%soil_water(k,i)

elseif (trim(fldname) == 'SFWAT_MASS'       ) then

   fldval = land%sfcwater_mass(k,i)

elseif (trim(fldname) == 'SFWAT_ENERGY'     ) then

   fldval = land%sfcwater_energy(k,i) * 1.e-3

elseif (trim(fldname) == 'SFWAT_TEMPC'      ) then

   call qtk(land%sfcwater_energy(k,i),tempk,fracliq)
   fldval = tempk - 273.15

elseif (trim(fldname) == 'SFWAT_FRACLIQ'    ) then

   call qtk(land%sfcwater_energy(k,i),tempk,fracliq)
   fldval = fracliq

elseif (trim(fldname) == 'SFWAT_DEPTH'      ) then

   fldval = 0.
   do klev = 1,land%nlev_sfcwater(i)
      fldval = fldval + land%sfcwater_depth(klev,i)
   enddo

elseif (trim(fldname) == 'NLEV_SFWAT'       ) then

   fldval = real(land%nlev_sfcwater(i))

elseif (trim(fldname) == 'LEAF_CLASS'       ) then

   fldval = real(land%leaf_class(i))

elseif (trim(fldname) == 'VEG_NDVIC'        ) then

   fldval = land%veg_ndvic(i)

elseif (trim(fldname) == 'VEG_TEMPC'        ) then

   fldval = land%veg_temp(i) - 273.15

elseif (trim(fldname) == 'VEG_WATER'        ) then

   fldval = land%veg_water(i)

elseif (trim(fldname) == 'STOM_RESIST'      ) then

   fldval = land%stom_resist(i)

elseif (trim(fldname) == 'GROUND_SHV'       ) then

   fldval = land%ground_shv(i) * 1.e3

elseif (trim(fldname) == 'SOIL_DEPTH'       ) then

   fldval = -slz(land%lsl(i))

elseif (trim(fldname) == 'SEATP'            ) then

   fldval = sea%seatp(i)   

elseif (trim(fldname) == 'SEATF'            ) then

   fldval = sea%seatf(i)

elseif (trim(fldname) == 'SEATC'            ) then

   fldval = sea%seatc(i)

elseif (trim(fldname) == 'SEAICEP'          ) then

   fldval = sea%seaicep(i)   

elseif (trim(fldname) == 'SEAICEF'          ) then

   fldval = sea%seaicef(i)

elseif (trim(fldname) == 'SEAICEC'          ) then

   fldval = sea%seaicec(i)

elseif (trim(fldname) == 'AREA'             ) then

   if (op%stagpt == 'S') then
      fldval = sea%area(i)
   elseif (op%stagpt == 'L') then
      fldval = land%area(i)
   endif

elseif (trim(fldname) == 'ROUGH'            ) then

   if (op%stagpt == 'S') then
      fldval = sea%rough(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rough(i)
   endif

elseif (trim(fldname) == 'CAN_TEMPC'        ) then

   if (op%stagpt == 'S') then
      fldval = sea%can_temp(i) - 273.15
   elseif (op%stagpt == 'L') then
      fldval = land%can_temp(i) - 273.15
   endif

elseif (trim(fldname) == 'CAN_SHV'          ) then

   if (op%stagpt == 'S') then
      fldval = sea%can_shv(i) * 1.e3
   elseif (op%stagpt == 'L') then
      fldval = land%can_shv(i) * 1.e3
   endif

elseif (trim(fldname) == 'SFC_TEMPC'        ) then

   if (op%stagpt == 'S') then
      fldval = sea%seatc(i) - 273.15
   elseif (op%stagpt == 'L') then
      nls = land%nlev_sfcwater(i)

      if (nls > 0) then
         call qtk(land%sfcwater_energy(nls,i),tempk,fracliq)
         fldval = tempk - 273.15
      else
         call qwtk(land%soil_energy(nzg,i)        &
                  ,land%soil_water(nzg,i)*1.e3    &
                  ,slcpd(land%ntext_soil(nzg,i))  &
                  ,tempk,fracliq)
         fldval = tempk - 273.15
      endif
   endif

elseif (trim(fldname) == 'SFC_SSH'          ) then

   if (op%stagpt == 'S') then
      fldval = sea%surface_ssh(i) * 1.e3
   elseif (op%stagpt == 'L') then
      fldval = land%surface_ssh(i) * 1.e3
   endif

elseif (trim(fldname) == 'SENSFLUX_LS'      ) then

   if (op%stagpt == 'S') then
      fldval = sea%sxfer_t(i) * cp
   elseif (op%stagpt == 'L') then
      fldval = land%sxfer_t(i) * cp
   endif

elseif (trim(fldname) == 'VAPFLUX_LS'       ) then

   if (op%stagpt == 'S') then
      fldval = sea%sxfer_r(i)
   elseif (op%stagpt == 'L') then
      fldval = land%sxfer_r(i)
   endif

elseif (trim(fldname) == 'LATFLUX_LS'       ) then

   if (op%stagpt == 'S') then
      fldval = sea%sxfer_r(i) * alvl
   elseif (op%stagpt == 'L') then
      fldval = land%sxfer_r(i) * alvl
   endif

elseif (trim(fldname) == 'RSHORT_LS'        ) then

   if (op%stagpt == 'S') then
      fldval = sea%rshort(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rshort(i)
   endif

elseif (trim(fldname) == 'RSHORT_DIFFUSE_LS') then

   if (op%stagpt == 'S') then
      fldval = sea%rshort_diffuse(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rshort_diffuse(i)
   endif

elseif (trim(fldname) == 'RLONG_LS'         ) then

   if (op%stagpt == 'S') then
      fldval = sea%rlong(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rlong(i)
   endif

elseif (trim(fldname) == 'RLONGUP_LS'       ) then

   if (op%stagpt == 'S') then
      fldval = sea%rlongup(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rlongup(i)
   endif

elseif (trim(fldname) == 'RLONG_ALBEDO_LS'  ) then

   if (op%stagpt == 'S') then
      fldval = sea%rlong_albedo(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rlong_albedo(i)
   endif

elseif (trim(fldname) == 'ALBEDO_BEAM_LS'   ) then

   if (op%stagpt == 'S') then
      fldval = sea%albedo_beam(i)
   elseif (op%stagpt == 'L') then
      fldval = land%albedo_beam(i)
   endif

elseif (trim(fldname) == 'ALBEDO_DIFFUSE_LS') then

   if (op%stagpt == 'S') then
      fldval = sea%albedo_diffuse(i)
   elseif (op%stagpt == 'L') then
      fldval = land%albedo_diffuse(i)
   endif

! New additions

elseif (trim(fldname) == 'FCELL_ILSF'        ) then

   fldval = real(i)

elseif (trim(fldname) == 'FCELL_IWLS'        ) then

   if (op%stagpt == 'S') then
      fldval = real(seaflux(i)%iws)
   elseif (op%stagpt == 'L') then
      fldval = real(landflux(i)%iwl)
   endif

elseif (trim(fldname) == 'FCELL_IW'        ) then

   if (op%stagpt == 'S') then
      fldval = real(seaflux(i)%iw)
   elseif (op%stagpt == 'L') then
      fldval = real(landflux(i)%iw)
   endif

elseif (trim(fldname) == 'FCELL_KW'        ) then

   if (op%stagpt == 'S') then
      fldval = real(seaflux(i)%kw)
   elseif (op%stagpt == 'L') then
      fldval = real(landflux(i)%kw)
   endif

elseif (trim(fldname) == 'FCELL_AREA'        ) then

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%area
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%area
   endif

elseif (trim(fldname) == 'FCELL_ARFATM'        ) then

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%arf_atm
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%arf_atm
   endif

elseif (trim(fldname) == 'FCELL_ARFLS'        ) then

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%arf_sea
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%arf_land
   endif

elseif (trim(fldname) == 'FCELL_SENS'        ) then

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%sfluxt * cp
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%sfluxt * cp
   endif

elseif (trim(fldname) == 'FCELL_VAP'         ) then

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%sfluxr
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%sfluxr
   endif

elseif (trim(fldname) == 'FCELL_LAT'         ) then

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%sfluxr * alvl
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%sfluxr * alvl
   endif

elseif (trim(fldname) == 'FCELL_AIRTEMPC'        ) then

   if (op%stagpt == 'S') then
      iw = seaflux(i)%iw
      if (iparallel == 1) then
         iw  = itabg_w(iw)%iw_myrank
      endif
      kw = seaflux(i)%kw
      fldval = theta(kw,iw) * (press(kw,iw) / p00) ** rocp - 273.15
   elseif (op%stagpt == 'L') then
      iw = landflux(i)%iw
      if (iparallel == 1) then
         iw  = itabg_w(iw)%iw_myrank
      endif
      kw = landflux(i)%kw
      fldval = theta(kw,iw) * (press(kw,iw) / p00) ** rocp - 273.15
   endif

elseif (trim(fldname) == 'FCELL_AIRTEMPK'        ) then

   if (op%stagpt == 'S') then
      iw = seaflux(i)%iw
      kw = seaflux(i)%kw
      fldval = theta(kw,iw) * (press(kw,iw) / p00) ** rocp
   elseif (op%stagpt == 'L') then
      iw = landflux(i)%iw
      kw = landflux(i)%kw
      fldval = theta(kw,iw) * (press(kw,iw) / p00) ** rocp
   endif

elseif (trim(fldname) == 'RHO_OBS'       ) then

   if (.not. allocated(rho_obs)) go to 1000

   fldval = rho_obs(k,itab_w(i)%inudp(1))

elseif (trim(fldname) == 'THETA_OBS'       ) then

   if (.not. allocated(theta_obs)) go to 1000

   fldval = theta_obs(k,itab_w(i)%inudp(1))

elseif (trim(fldname) == 'SHW_OBS'       ) then

   if (.not. allocated(shw_obs)) go to 1000

   fldval = shw_obs(k,itab_w(i)%inudp(1)) * 1.e3

elseif (trim(fldname) == 'UZONAL_OBS'       ) then

   if (.not. allocated(uzonal_obs)) go to 1000

   fldval = uzonal_obs(k,itab_w(i)%inudp(1))

elseif (trim(fldname) == 'UMERID_OBS'       ) then

   if (.not. allocated(umerid_obs)) go to 1000

   fldval = umerid_obs(k,itab_w(i)%inudp(1))

elseif (trim(fldname) == 'RHO_SIM'       ) then

   if (.not. allocated(rho_sim)) go to 1000

   fldval = rho_sim(k,itab_w(i)%inudp(1))

elseif (trim(fldname) == 'THETA_SIM'       ) then

   if (.not. allocated(theta_sim)) go to 1000

   fldval = theta_sim(k,itab_w(i)%inudp(1))

elseif (trim(fldname) == 'SHW_SIM'       ) then

   if (.not. allocated(shw_sim)) go to 1000

   fldval = shw_sim(k,itab_w(i)%inudp(1)) * 1.e3

elseif (trim(fldname) == 'UZONAL_SIM'       ) then

   if (.not. allocated(uzonal_sim)) go to 1000

   fldval = uzonal_sim(k,itab_w(i)%inudp(1))

elseif (trim(fldname) == 'UMERID_SIM'       ) then

   if (.not. allocated(umerid_sim)) go to 1000

   fldval = umerid_sim(k,itab_w(i)%inudp(1))

elseif (trim(fldname) == 'PRESS_SFC' .or.  &
        trim(fldname) == 'PRESS_SFC_P') then

   im1 = itab_w(i)%im1
   im2 = itab_w(i)%im2
   im3 = itab_w(i)%im3
   
   ka = lpw(i)
   
   topt = (topm(im1) + topm(im2) + topm(im3)) / 3.

   fldval = (press(ka,i) + grav * rho(ka,i) * (zt(ka) - topt)) * .01

   if (trim(fldname) == 'PRESS_SFC_P') then
      
      fldval = fldval  &
         - (press_init(ka,i) + grav * rho_init(ka,i) * (zt(ka) - topt)) * .01

   endif

elseif (trim(fldname) == 'ADDSC1'          ) then

   if (.not. allocated(addsc(1)%sclp)) go to 1000

   fldval = wtbotw * addsc(1)%sclp(kt  ,i) &
          + wttopw * addsc(1)%sclp(kt+1,i)

elseif (trim(fldname) == 'ADDSC2'          ) then

   if (.not. allocated(addsc(2)%sclp)) go to 1000

   fldval = wtbotw * addsc(2)%sclp(kt  ,i) &
          + wttopw * addsc(2)%sclp(kt+1,i)

elseif (trim(fldname) == 'ADDSC3'          ) then

   if (.not. allocated(addsc(3)%sclp)) go to 1000

   fldval = wtbotw * addsc(3)%sclp(kt  ,i) &
          + wttopw * addsc(3)%sclp(kt+1,i)

elseif (trim(fldname) == 'ADDSC4'          ) then

   if (.not. allocated(addsc(4)%sclp)) go to 1000

   fldval = wtbotw * addsc(4)%sclp(kt  ,i) &
          + wttopw * addsc(4)%sclp(kt+1,i)

elseif (trim(fldname) == 'ADDSC5'          ) then

   if (.not. allocated(addsc(5)%sclp)) go to 1000

   fldval = wtbotw * addsc(5)%sclp(kt  ,i) &
          + wttopw * addsc(5)%sclp(kt+1,i)

elseif (trim(fldname) == 'MRLW'            ) then

   fldval = itab_w(i)%mrlw

elseif (trim(fldname) == 'MROW'            ) then

   fldval = itab_w(i)%mrow

elseif (trim(fldname) == 'MROWH'            ) then

   fldval = itab_w(i)%mrowh

elseif (trim(fldname) == 'IMGLOBE'            ) then

   fldval = itab_m(i)%imglobe

elseif (trim(fldname) == 'IUGLOBE'            ) then

   fldval = itab_u(i)%iuglobe

elseif (trim(fldname) == 'IWGLOBE'            ) then

   fldval = itab_w(i)%iwglobe

elseif (trim(fldname) == 'U_IRANK'            ) then

!   fldval = itabg_u(i)%irank
   fldval = itab_u(i)%irank

elseif (trim(fldname) == 'W_IRANK'            ) then

!   fldval = itabg_w(i)%irank
   fldval = itab_w(i)%irank

else

   go to 1000

endif

if (infotyp == 'VALUE') then
   if (op%fldval_min > fldval) op%fldval_min = fldval
   if (op%fldval_max < fldval) op%fldval_max = fldval
elseif (infotyp == 'VALUV') then
   if (op%fldvalv_min > fldval) op%fldvalv_min = fldval
   if (op%fldvalv_max < fldval) op%fldvalv_max = fldval
endif

! normal RETURN

return

! RETURN for "Plot field UNAVAILABLE"

1000 continue

notavail = 3

return
end subroutine oplot_lib

