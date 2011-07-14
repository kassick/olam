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

subroutine read_nl(file)

use oname_coms, only: nl
use misc_coms,  only: io6
use rastro_evts

implicit none

character(len=*), intent(in) :: file
character(len=9) :: line
integer :: iline, iplt
logical :: fexists
real, external :: walltime


namelist /OLAMIN/ nl

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_READ_NL_IN,file)
#endif

! OPEN THE NAMELIST FILE
inquire(file=file, exist=fexists)


if (.not. fexists) then
   write(io6,*) "The namelist file "//trim(file)//" is missing."
   stop "Stopping model run."
endif

open(10, status='OLD', file=file)
 
! READ GRID POINT AND OPTIONS INFORMATION FROM THE NAMELIST

read(10, nml=OLAMIN)

! READ THE MODEL PLOT FIELDS IN COLUMN FORM AT THE END OF THE NAMELIST

rewind(10)
do iline = 1,1000
   read(10,*,end=1999) line
   if (line(1:7) == 'FLDNAME') exit
enddo

read(10,*,end=1999) line

do iplt = 1,nl%nplt
   read(10,*) nl%fldname(iplt)     &
             ,nl%projectn(iplt)    &
             ,nl%icolortab(iplt)   &
             ,nl%pltspec2(iplt)    &
             ,nl%plotcoord1(iplt)  & 
             ,nl%plotcoord2(iplt)  &
             ,nl%slabloc(iplt)     &
             ,nl%plotwid(iplt)     &
             ,nl%viewazim(iplt) 
enddo

1999 continue

close(10)

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_READ_NL_OUT,file)
#endif

end subroutine read_nl

!===============================================================================

subroutine copy_nl(copy_type)

use max_dims,    only: maxgrds, nzgmax, maxisdirs
use oname_coms,  only: nl
use misc_coms,   only: io6, expnme, runtype, timeunit, timmax8, ndtrat,  &
                       nacoust, idiffk, zkhkm, xkhkm, csz, csx, akmin,   &
                       dtlong, initial, zonclim, topo_database,   &
                       gridfile, hfilin, ioutput, hfilepref, iclobber,   &
                       frqstate, naddsc, icorflg, ilwrtyp, iswrtyp, radfrq,  &
                       lonrad, nqparm, confrq, wcldbs, cnum_vel,  &
                       cnum_sclr, nsndg, ipsflg, itsflg, irtsflg, iusflg,  &
                       hs, p_sfc, us, vs, ts, ps, rts,  &
                       itime1, idate1, imonth1, iyear1, ngrids, nzp,  &
                       mdomain, itopoflg, nxp,  &
                       centlat, centlon, grdlen, grdwid, grdaxis,  &
                       dzrat, dzmax, deltax, deltaz, zbase,  &
                       current_time, nqparm_sh

use micro_coms,  only: level, icloud, irain, ipris, isnow, iaggr,  &
                       igraup, ihail, cparm, rparm, pparm, sparm, aparm,  &
                       gparm, hparm

use leaf_coms,   only: nvgcon, nslcon, slmstr, isoilflg, ndviflg,  &
                       isfcl, ivegflg, nzg, nzs, slz,  &
                       veg_database, soil_database,  &
                       ndvi_database, iupdndvi, landusefile, &
                       isoilstateinit, isoildepthflg,soilstate_db,soildepth_db

use sea_coms,    only: isstflg, sst_database, seatmp, seafile, iupdsst,  &
                       iseaiceflg, seaice_database, iupdseaice

use mem_ed,      only: n_soi, soi_lat, soi_lon, n_ed_region, ed_reg_latmin,  &
                       ed_reg_latmax, ed_reg_lonmin, ed_reg_lonmax
use oplot_coms,  only: op
use isan_coms,   only: iapr, isdirs
use mem_nudge,   only: tnudcent, nudflag, nudrad
use mem_rayf,    only: rayf_zmin, rayf_distim, rayf_expon,  &
                       rayfw_zmin, rayfw_distim, rayfw_expon
use ed_options,  only: ied_init_mode, istoma_scheme, iphen_scheme,  &
     ed_inputs_dir, n_plant_lim, n_decomp_lim, include_fire, ied_offline, &
     metcyc1, metcyc2, ed_offline_db,ianth_disturb, runoff_time, ed_hfilin

use disturbance_coms, only: treefall_disturbance_rate
use rastro_evts

implicit none

character(len=*) :: copy_type
integer :: i,j
real(kind=8) :: tfact

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_COPY_NL_IN,copy_type)
#endif

if (copy_type == 'ALL_CASES') then
	
! The namelist variables in this section will always be used from the
! namelist, regardless of which type of run this is (a history start or not).
! This allows model options to be changed if it is a history start.
	
   expnme   = nl%expnme
   runtype  = nl%runtype
   timeunit = nl%timeunit
   timmax8  = nl%timmax

!----------------------------------------------------------
! Convert timmax8 units to seconds if necessary

   if (timeunit == 'd' .or. timeunit == 'D') tfact = 86400.d0
   if (timeunit == 'h' .or. timeunit == 'H') tfact = 3600.d0
   if (timeunit == 'm' .or. timeunit == 'M') tfact = 60.d0
   if (timeunit == 's' .or. timeunit == 'S') tfact = 1.d0

   timmax8 = timmax8 * tfact
!----------------------------------------------------------

   do i = 1,maxgrds
      ndtrat(i)  = nl%ndtrat(i)
      nacoust(i) = nl%nacoust(i)
      idiffk(i)  = nl%idiffk(i)
      zkhkm(i)   = nl%zkhkm(i)
      xkhkm(i)   = nl%xkhkm(i)
      csz(i)     = nl%csz(i)
      csx(i)     = nl%csx(i)
      akmin(i)   = nl%akmin(i)
      nqparm(i)  = nl%nqparm(i)
      nqparm_sh(i)  = nl%nqparm_sh(i)
   enddo

   topo_database = nl%topo_database
   sst_database  = nl%sst_database
   seaice_database  = nl%seaice_database
   veg_database  = nl%veg_database
   soil_database = nl%soil_database
   ndvi_database = nl%ndvi_database

   soilstate_db = nl%soilstate_db
   soildepth_db = nl%soildepth_db

   n_soi         = nl%n_soi
   soi_lat       = nl%soi_lat
   soi_lon       = nl%soi_lon
   n_ed_region   = nl%n_ed_region
   ed_reg_latmin = nl%ed_reg_latmin
   ed_reg_latmax = nl%ed_reg_latmax
   ed_reg_lonmin = nl%ed_reg_lonmin
   ed_reg_lonmax = nl%ed_reg_lonmax
   ed_offline_db = trim(nl%ed_offline_db)
   ied_init_mode = nl%ied_init_mode
   metcyc1       = nl%metcyc1
   metcyc2       = nl%metcyc2
   ied_offline   = nl%ied_offline
   istoma_scheme = nl%istoma_scheme
   ianth_disturb = nl%ianth_disturb
   iphen_scheme  = nl%iphen_scheme
   ed_inputs_dir = trim(nl%ed_inputs_dir)
   ed_hfilin = trim(nl%ed_hfilin)
   n_plant_lim   = nl%n_plant_lim
   n_decomp_lim  = nl%n_decomp_lim
   include_fire  = nl%include_fire
   treefall_disturbance_rate = nl%treefall_disturbance_rate
   runoff_time   = nl%runoff_time

   dtlong        = nl%dtlong
   initial       = nl%initial
   zonclim       = nl%zonclim
   tnudcent      = nl%tnudcent
   nudflag       = nl%nudflag
   hfilin        = nl%hfilin
   ioutput       = nl%ioutput
   hfilepref     = nl%hfilepref
   iclobber      = nl%iclobber
   frqstate      = nl%frqstate
   gridfile      = nl%gridfile
   landusefile   = nl%landusefile
   seafile       = nl%seafile
   iupdsst       = nl%iupdsst
   iupdseaice    = nl%iupdseaice
   iupdndvi      = nl%iupdndvi
   naddsc        = nl%naddsc
   icorflg       = nl%icorflg
   rayf_zmin     = nl%rayf_zmin
   rayf_distim   = nl%rayf_distim
   rayf_expon    = nl%rayf_expon
   rayfw_zmin    = nl%rayfw_zmin
   rayfw_distim  = nl%rayfw_distim
   rayfw_expon   = nl%rayfw_expon
   ilwrtyp       = nl%ilwrtyp
   iswrtyp       = nl%iswrtyp
   radfrq        = nl%radfrq
   lonrad        = nl%lonrad
   confrq        = nl%confrq
   wcldbs        = nl%wcldbs
   seatmp        = nl%seatmp
   cnum_vel      = nl%cnum_vel
   cnum_sclr     = nl%cnum_sclr
   level         = nl%level
   icloud        = nl%icloud
   irain         = nl%irain
   ipris         = nl%ipris
   isnow         = nl%isnow
   iaggr         = nl%iaggr
   igraup        = nl%igraup
   ihail         = nl%ihail
   cparm         = nl%cparm
   rparm         = nl%rparm
   pparm         = nl%pparm
   sparm         = nl%sparm
   aparm         = nl%aparm
   gparm         = nl%gparm
   hparm         = nl%hparm
   nvgcon        = nl%nvgcon
   nslcon        = nl%nslcon
   nsndg         = nl%nsndg
   ipsflg        = nl%ipsflg
   itsflg        = nl%itsflg
   irtsflg       = nl%irtsflg
   iusflg        = nl%iusflg
   nudrad        = nl%nudrad
   isstflg       = nl%isstflg
   iseaiceflg    = nl%iseaiceflg
   isoilflg      = nl%isoilflg
   isoilstateinit = nl%isoilstateinit
   isoildepthflg = nl%isoildepthflg
   ndviflg       = nl%ndviflg
   isdirs        = nl%isdirs

   slmstr(1:nzgmax) = nl%slmstr(1:nzgmax)

   iapr(1:maxisdirs) = nl%iapr(1:maxisdirs)

   hs(1) = nl%hs
   p_sfc = nl%p_sfc
   do i = 1,nl%nsndg
      ps(i)  = nl%sounding(1,i)
      ts(i)  = nl%sounding(2,i)
      rts(i) = nl%sounding(3,i)
      us(i)  = nl%sounding(4,i)
      vs(i)  = nl%sounding(5,i)
   enddo

! Variables from $MODEL_PLOT namelist

   op%nplt       = nl%nplt
   op%nplt_files = nl%nplt_files
   op%frqplt     = nl%frqplt
   op%dtvec      = nl%dtvec
   op%headspeed  = nl%headspeed
   op%pltname    = nl%pltname
   op%plttype    = nl%plttype
   op%pltorient  = nl%pltorient
   op%vec_maxmrl = nl%vec_maxmrl

   do i = 1,op%nplt
      op%fldname(i)   = nl%fldname(i)
      op%projectn(i)  = nl%projectn(i)
      op%icolortab(i) = nl%icolortab(i)
      op%slabloc(i)   = nl%slabloc(i)

! Initialize to default null values...
      op%contrtyp(i)        = 'N'
      op%prtval(i)          = 'N'
      op%pltindx(i)         = 'N'         
      op%vectbarb(i)        = 'N'        
      op%pltgrid(i)         = 'N'       
      op%pltgrid_landsea(i) = 'N'       
      op%pltdualgrid(i)     = 'N'       
      op%pltborder(i)       = 'N'       
      op%colorbar(i)        = 'N'   
      op%labelbar(i)        = 'N'   
      op%maptyp(i)          = 'N'       
      op%pltcone(i)         = 'N'        
      op%pltlev(i)          = 'N'   
      op%windowin(i)        = 'N'   
      op%frameoff(i)        = 'N'  
      op%panel(i)           = 'N'   

      do j = 1,13
         if (nl%pltspec2(i)(j:j) == 'T') op%contrtyp(i) = 'T'         
         if (nl%pltspec2(i)(j:j) == 'F') op%contrtyp(i) = 'F'         
         if (nl%pltspec2(i)(j:j) == 'L') op%contrtyp(i) = 'L'         
         if (nl%pltspec2(i)(j:j) == 'O') op%contrtyp(i) = 'O'         

         if (nl%pltspec2(i)(j:j) == 'P') op%prtval(i) = 'P'         

         if (nl%pltspec2(i)(j:j) == 'I') op%pltindx(i) = 'I'         
         if (nl%pltspec2(i)(j:j) == 'J') op%pltindx(i) = 'J'         

         if (nl%pltspec2(i)(j:j) == 'B') op%vectbarb(i) = 'B'         
         if (nl%pltspec2(i)(j:j) == 'U') op%vectbarb(i) = 'U'         
         if (nl%pltspec2(i)(j:j) == 'u') op%vectbarb(i) = 'u'         
         if (nl%pltspec2(i)(j:j) == 'V') op%vectbarb(i) = 'V'         

         if (nl%pltspec2(i)(j:j) == 'G') op%pltgrid(i) = 'G'         
         if (nl%pltspec2(i)(j:j) == 'g') op%pltgrid_landsea(i) = 'g'         

         if (nl%pltspec2(i)(j:j) == 'D') op%pltdualgrid(i) = 'D'         
         if (nl%pltspec2(i)(j:j) == 'd') op%pltdualgrid(i) = 'd'         

         if (nl%pltspec2(i)(j:j) == 'b') op%pltborder(i) = 'b'         
         if (nl%pltspec2(i)(j:j) == 't') op%pltborder(i) = 't'         

         if (nl%pltspec2(i)(j:j) == 'n') op%labelbar(i) = 'n'         
         if (nl%pltspec2(i)(j:j) == 'i') op%labelbar(i) = 'i'         

         if (nl%pltspec2(i)(j:j) == 'c') op%colorbar(i) = 'c'         

         if (nl%pltspec2(i)(j:j) == 'M') op%maptyp(i) = 'M'         
         if (nl%pltspec2(i)(j:j) == 'm') op%maptyp(i) = 'm'         

         if (nl%pltspec2(i)(j:j) == 'C') op%pltcone(i) = 'C' 

         if (nl%pltspec2(i)(j:j) == 'p') op%pltlev(i) = 'p' 
         if (nl%pltspec2(i)(j:j) == 's') op%pltlev(i) = 's' 

         if (nl%pltspec2(i)(j:j) == 'W') op%windowin(i) = 'W' 

         if (nl%pltspec2(i)(j:j) == 'f') op%frameoff(i) = 'f' 

         if (nl%pltspec2(i)(j:j) == '1') op%panel(i) = '1' 
         if (nl%pltspec2(i)(j:j) == '2') op%panel(i) = '2' 
         if (nl%pltspec2(i)(j:j) == '3') op%panel(i) = '3' 
         if (nl%pltspec2(i)(j:j) == '4') op%panel(i) = '4' 
      enddo        

   enddo

   do i = 1,op%nplt_files
      op%plt_files(i) = nl%plt_files(i)
   enddo

elseif (copy_type == 'NOT_HISTORY') then
	
! The namelist variables in this section either must not be changed on a
! history restart or changing them would be irrelevant.  Thus, they are
! only copied to main model memory if this is not a history restart.

   itime1   = nl%itime1
   idate1   = nl%idate1
   imonth1  = nl%imonth1
   iyear1   = nl%iyear1
   ngrids   = nl%ngrids
   nzp      = nl%nzp
   nzg      = nl%nzg
   nzs      = nl%nzs
   mdomain  = nl%mdomain
   isfcl    = nl%isfcl
   itopoflg = nl%itopoflg
   ivegflg  = nl%ivegflg
   nxp      = nl%nxp

   do i = 1,maxgrds
      centlat(i) = nl%centlat(i)
      centlon(i) = nl%centlon(i)
      grdlen(i)  = nl%grdlen(i)
      grdwid(i)  = nl%grdwid(i)
      grdaxis(i) = nl%grdaxis(i)
   enddo
 
   dzrat     = nl%dzrat
   dzmax     = nl%dzmax
   deltax    = nl%deltax
   deltaz    = nl%deltaz
   zbase     = nl%zbase

   slz(1:nzgmax) = nl%slz(1:nzgmax)

   ! set current time to initial time here.  If this is a history run,
   ! reset current time in subroutine history_start.
   current_time%year = iyear1
   current_time%month = imonth1
   current_time%date = idate1
   current_time%time = int(itime1 * 0.01) * 3600.0   &
        + (itime1 * 0.01 - int(itime1*0.01))*100.0*60.0
endif

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_COPY_NL_OUT,copy_type)
#endif

end subroutine copy_nl

!===============================================================================

subroutine namelist_print()
  
use oname_coms, only: nl
use misc_coms,  only: io6

implicit none

namelist /OLAMIN/ nl

write(io6, nml=OLAMIN)

end subroutine namelist_print
