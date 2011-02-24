MODULE FUINPUT
implicit none
! PARAMETERS
integer, parameter :: nvx = 140!! MAX # of Model Levels
integer, parameter :: mccx= 26	!! MAX # of Cloud Conditions
integer, parameter :: mxovl=10  !! MAX # of Cloud Overlap cases!!!! FU20040901
integer, parameter :: mxat= 13  	!! Max # Aerosol Wavelengths to input
integer, parameter :: mxac= 8  	!! Max # Aerosol Constituents to input 


integer, parameter :: nv1x  = nvx + 1 
integer, parameter :: ndfsx = nvx, mdfsx = nvx + 1 
integer, parameter :: ndfs4x= 4*ndfsx, mbx = 18, mbsx = 6 
integer, parameter :: mbirx = 12 , ncx=8, mby = 10 

integer, parameter :: naer=25	!! # of aerosol property types
integer, parameter :: nrh=8 	!! # of aerosol Rel Humidity

integer, parameter :: iaform =3 !! Multi Constituent , Multi wavelength 
integer, parameter :: ivd=1 	!! USER PROVIDED AEROSOL VERTICAL PROFILE
integer, parameter :: ifg=0 	!! Aerosol property dependence on layer RH
integer, parameter :: nscm=4	!! # of solver modes(Clr, TOTAL,Pris,TOTALNoaot,Kato w/AOT ,KATO w/o AOT )

integer, parameter :: idkfr=2	!! Window K's scheme Kratz
integer, parameter :: iwtas=4	!! Continuum averaging scheme
real,parameter :: omit = -9999.0!! OMIT FLAG
!VLA parameters
integer, parameter ::  mfix  =  nvx ! Max # of FIXED/manditory Model levels
integer, parameter ::  mflo  =   10 ! Max # of Floating Model levels
integer, parameter ::  maxin =  nvx ! Max # of Input Profile levels
integer, parameter ::  maxrep=  nvx ! Max # of Output Levels to Report.
integer, parameter ::  mri = 1100
!----------------------------------------------------------

integer, save :: nv1,mb,mbs,mbir,nc,ndfs,mdfs,ndfs4,ndfs2 
logical ,save :: edding, quadra, hemisp, mquadr
logical ,save :: fourssl, foursir
logical, save ::  lscm(nscm)
logical, save ::  fl93i ,lband6a,lray

integer, save :: irobckd  ! H20 continuum index suggest set to 5 :ckd2.4
integer, save :: nhb  	! # if hidden bands (0-2) >2200cm-1
integer, save :: isolar_spectrum, isksw , isksolve , instrument

real ,save    ::   cfc_conc(3),umco2,umch4,umn2o 

!----
integer	 nv	! # of MODEL LAYERS
real ee(mbirx),ss ,pts,u0,ur
real dz(nvx) ! thickness of model layer in km

real ,dimension(nv1x):: pp,pt,ph,po
real ,dimension(nvx) :: pre,plwc,pde,piwc ,pdge

real, dimension(nvx,6,mccx):: optwat,optice

integer  	n_atau ,nac,itps(mxac)	
real 		a_wli(mxat)
real 		a_taus(mxat,mxac)
real	 	aprofs(nvx,mxac)

real ,dimension(nvx)::cf_sk,ccf_sk,nu_sk
real ,dimension(mccx,mxovl)::shapefactor
real,save :: taucorrection(mccx,mxovl)
!------------------------------------------------------------------------

! Direct Cloud property Insert
TYPE DPI_TYPE 
 sequence
logical ldpi
integer phase(nvx)
real plwc_iwc(nvx)
real   pre_de(nvx)
END TYPE DPI_TYPE 
!-------------------------------------------------------------------------


TYPE FUSUBCLD
 sequence
 real shapef
 real mn_lin_tau
END TYPE FUSUBCLD

TYPE FUCLDCND
sequence
real    cldfrac
integer novl    	!!!! FU20040901
integer iphase(mxovl)
integer icld_top(mxovl)
integer icld_bot(mxovl)
real    tau_vis(mxovl)
real    part_size(mxovl)

TYPE (FUSUBCLD) sc(mxovl)

TYPE (DPI_TYPE) dpi
END TYPE FUCLDCND

!===============================================================================
!!Variable  Level Algorithm (VLA)
!-------------------------------------------------------------------------------------

TYPE VLADEFINE
 sequence
 integer nfix		! # of FIXED levels Input
 real pfix(mfix)	! Array of Pressure(hPa) of fixed levels
 integer nflo		! # of FLOATING levels Input
 real pflo(mflo)	! Array of Pressure(hPa) of floating levels
 real cldpres(2,mccx,mxovl) ! Array of Pressure(hPa) Top(1) , Bot(2) for each cloud condition an overlap condition 
 integer nrep		! # of levels to find reporting level index upon output
 real report(maxrep)	! Array of Pressure levels to find reporting level indexes for
END TYPE VLADEFINE

TYPE (VLADEFINE) VD

TYPE VLAIN
 sequence
 integer nlev      ! # of Input levels in input profile
 real	pp(maxin)  !pressure profile (mb)
 real	pt(maxin)  !temperature profile (K)
 real	ph(maxin)  !H2o Mixing ratio profile (g/g)
 real	po(maxin)  !O3 Mixing ratio profile (g/g)
 real   hsfc       ! SURFACE GEOPOTENTIAL of Sounding.
END TYPE VLAIN

TYPE (VLAIN)  VI

TYPE VLAOUT
 sequence
 integer nlev	   ! # of Input levels in output profile
  real	pp(nv1x)  !pressure profile (mb)
 real	pt(nv1x)  !temperature profile (K)
 real	ph(nv1x)  !H2o Mixing ratio profile (g/g)
 real	po(nv1x)  !O3 Mixing ratio profile (g/g)
 integer cldindex(2,mccx,mxovl) ! Top(1)Bot(2) Cloud Layer Indexes for FuModel
 integer ireport(maxrep)	! indexes of reporting pressure levels.
 integer ri_pp(mri)		! reverse index of pressure Layer to nearest 1hpa
 integer ierr	!ERROR FLAG...  <0 warning   ,, >0 ERROR!!
END TYPE VLAOUT

TYPE (VLAOUT) VO





!------------------------------------------------------------------------
TYPE FUINPUTTYPE
sequence
integer	 nv		!# of model LAYERS
logical lscm(nscm)	!Select solver configurations modes
real	u0		!Cosine Solar Zenith angle ( 0-1)
real	ur 		!Cosine View Zenith angle ( 0-1) for LW Radiance
logical fl93i		! False New 98Ice cld (Dge) ::True=Old 1993 ice De
logical lband6a		! Treat Solar > 4um 
logical lray		! New Rayleigh Treatment 
integer irobckd		! Continuum Treatment 
integer nhb  		! # of hidden bands (0-2) >2200cm-1
real	sfcalb(15,2,0:mccx)	! New Surface albedo array(BAND,w w/oAOT, Clr:CldCnds)
real 	ee(mbirx)	!Surface Spectral Emissivity Longwave
real	ss		!Solar Constant (wm-2)
real	pts		!Skin Temperature (k)
real	pp(nv1x)	!pressure profile (mb)
real	pt(nv1x)	!temperature profile (K)
real	ph(nv1x)	!H2o Mixing ratio profile (g/g)
real	po(nv1x)	!O3 Mixing ratio profile (g/g)
real	umco2 		! Co2 concentration ( LW only)
real	umch4		! CH4 concentration ( LW only)
real	umn2o 		! N2O concentration ( LW only)
real	cfc_conc(3)	! CFC concentrations (LW only)
integer n_atau 		!# of Aerosol wavelength inputs 
integer nac		!# of Aerosol constituents
integer itps(mxac)	!# Aerosol types 
real 	a_wli(mxat)	!Wavelength of AOT input
real	a_taus(mxat,mxac) 	!AOT at Wavelength 
real	aprofs(nvx,mxac)	!fractional profiles of AOT for each constituent
integer isolar_spectrum ! Solar spectral shape assumption. : Suggest isolar_spectrum=4
integer isksw		! Correlated K option seiji k  : Suggest isksw=4
integer isksolve	! Solver option 0=Fu 1=Kato
integer instrument	! Ceres Window emmulation
logical fourssl		! Fu four-stream Solar
logical foursir		! Fu four-stream IR
integer wp_hgt_flag     ! 0 = CONSTANT MIX , 1=Fall off Linear w/height
 real   hsfc   		! OPTIONAL: Surface Geopotential Height (meters ) 
 real zz(nv1x)		! NOT AN INPUT COMPUTED BY thicks Geopotential height (meters)
 TYPE(FUCLDCND) 	fc(mccx)

 TYPE (VLADEFINE) VD
 TYPE (VLAIN)  VI
 TYPE (VLAOUT) VO
 
 TYPE (DPI_TYPE) DPI
integer ierror

END TYPE FUINPUTTYPE
	

TYPE(FUINPUTTYPE), save::	fi



CONTAINS 
!======================================================================================
subroutine set_default_options_fu
implicit none
!Trace Gas Concentrations

fi%umco2	     = 360. ! 360 ppmv CO2 Concetration (Affects LW only , SW Fixed)
fi%umch4	     = 1.75 ! 1.75 ppmv CH4 Concetration (Affects LW only , SW Fixed)
fi%umn2o	     = 0.31 ! 0.31 ppmv N2O Concetration (Affects LW only , SW Fixed)
fi%cfc_conc(1:3)     = (/0.268e-09 , 0.503e-09 ,0.105e-09/)

!Options 

fi%isolar_spectrum = 4		!4=Mod3.7_Newkur
fi%fl93i	   = .false.	!F=NEW ICE CLOUDS fu98 w dge
fi%lband6a	   = .true.	!T=Treat >4.0 micron Solar
fi%lray		   = .true.	!T=NEW Rayleigh Treatment
fi%irobckd	   = 5		!5=CKD2.4 Continuum
fi%nhb		   = 2		!2=Treat 2bands >2200cm-1
fi%isksw	   = 1		!1=Seiji Hitran2000 SW K's
fi%isksolve	   = 0		!0=Original Fu Solver 1=Kato Solver

!Full Suite of Calculations
fi%lscm	   = .true. ! USE all  SOLVER MODES ( CLEAR,TOTAL ,Pristine,NoaotCld )

fi%wp_hgt_flag = 0  !0 = cloud lwc constant w/height :: 1=cloudlwc linearly increases w/height, max @cldtop

fi%fc%dpi%ldpi =  .false. ! false do not directly insert properties from FI%DPI

fi%fourssl = .false.  ! False =Two Stream Solar ,,True = FOUR stream
fi%foursir = .false.  ! False = Two/four IR     ,,True = FOUR stream
fi%ierror = 0 ! initalize to no error condition
fi%vo%ierr= 0 ! initalize to no error condition vla algorithm

end subroutine set_default_options_fu
!=================================================================
subroutine assigninputs
implicit none


nv	= 	fi%nv
lscm	=	fi%lscm
u0	=	fi%u0
ur	=	fi%ur
fl93i	= 	fi%fl93i
lband6a =	fi%lband6a
lray 	=	fi%lray
irobckd = 	fi%irobckd
nhb  	= 	fi%nhb
! surface albedo directly input
ee	=	fi%ee
ss	=	fi%ss
pts	=	fi%pts
pp	=	fi%pp
pt	=	fi%pt
ph	=	fi%ph
po	=	fi%po
umco2	=	fi%umco2
umch4	=	fi%umch4
umn2o	=	fi%umn2o
cfc_conc=	fi%cfc_conc
n_atau	=	fi%n_atau
nac	=	fi%nac
itps	=	fi%itps
a_wli	=	fi%a_wli
a_taus	=	fi%a_taus
aprofs	=	fi%aprofs
isolar_spectrum	=fi%isolar_spectrum
isksw	= 	fi%isksw
isksolve=	fi%isksolve
instrument= 	fi%instrument
fourssl	=	fi%fourssl
foursir	=	fi%foursir

!Secondary...................
nv1	=	nv+1
ndfs 	= 	nv
mdfs 	= 	nv1
ndfs4 	= 	4*ndfs 
ndfs2	=  	2*ndfs
mb   	= 18
mbs  	= 6
mbir 	= 12 
nc   	= 8 

end subroutine assigninputs

!=================================================================
subroutine check_inputs_fu
  implicit none
  integer i,j
  
  fi%ierror = 0 !!! SET ERROR FLAG TO ZERO means OK..
  
  ! CHECK BASICS
  if ( fi%nv > nvx .or. fi%nv < 1 )	  fi%ierror=100 ! #of Model Levels out of range
  if ( fi%ierror .ne. 0 ) return
  if ( abs(fi%u0) > 1.0 )			  fi%ierror=101 ! Cos Sol out of range
  if ( abs(fi%ur) > 1.0 )			  fi%ierror=102 ! Cos View Zen of range
  if ( fi%ss < 1000 .or. fi%ss > 1500)	  fi%ierror=103 ! Solar Const out of range
  
  
  !CHECK ATMOSPHERE....
  do i=1,fi%nv+1
     if( fi%pp(i)< 1E-5  .or. fi%pp(i)>1100 )  fi%ierror=104 !Pressure(hPa) out of range
     if (fi%pt(i)< 100.  .or. fi%pt(i)> 450.)then
        fi%ierror=105 !Temperature (K) out of range
        print*,i,fi%pt(i)
     endif
     if (fi%ph(i)<0.0 .or. fi%ph(i)> 0.5 )  fi%ierror=106 !H20  Mix Ratio(g/g) out of range
     if (fi%po(i)<1.E-20 .or. fi%po(i)>1.E-2)then
        fi%ierror=107 !O3 Mix Ratio(g/g) out of range
        print*,'o3 mixing ratio',i,fi%po(i)
     endif
  enddo
  do i=1,fi%nv
     if (fi%pp(i)>fi%pp(i+1) )		  fi%ierror=108 !Pressure Not Monotonic increasing
  enddo
  
  if( fi%pts < 100. .or.  fi%pts>500.)	  fi%ierror=109 !Skin Temp out of range
  
  
  !CHECK AEROSOLS.....
  if(fi%nac < 0 .or. fi%nac > mxac )	  fi%ierror=110 !#of Aerosol Constituents out of range
  if ( fi%ierror .ne. 0 ) return
  if(fi%n_atau < 1 .or. fi%n_atau > mxat )  fi%ierror=111 !#of Aerosol Optical Depths out of range
  if ( fi%ierror .ne. 0 ) return
  
  do j = 1,fi%nac
     if (fi%itps(j) < 1 .or. fi%itps(j) > naer ) fi%ierror=112 !Aerosol Type out of range
     do i = 1,fi%n_atau 
        if ( fi%a_wli(i) < 0.2 .or. fi%a_wli(i)> 3.0 )   fi%ierror=113 !Aerosol wavelength (um) out of range
        if ( fi%a_taus(i,j) < 0.0 .or. fi%a_taus(i,j)> 50.0 ) then
           fi%ierror=114 !Aerosol Opt depth out of range
           print*,'Aerosol optical depth out of range'
           print*,fi%a_taus(i,j),i,j
           stop
        endif
     enddo
  enddo
     
  !!CHECK OPTION SETTINGS
  if ( fi%irobckd < 0 .or. fi%irobckd >5 ) fi%ierror=115 ! Continuum option out of range
  if ( fi%nhb < 0 .or. fi%nhb >2 ) fi%ierror=116 ! #of hidden thermal bands >2200cm-1 option out of range
  if (fi%isksw <0 .or. fi%isksw > 1 ) fi%ierror=117 ! SW K's :: 0=Fu    1=Seiji Kato HT200 K's in SW 
  if (fi%isksolve<0.or.fi%isksolve>1)  fi%ierror=118 ! Solver Mode 0=Fu 1=Seiji Kato overlap &cld tau shape factor
  if (fi%isolar_spectrum <0 .or. fi%isolar_spectrum >6) fi%ierror=128 ! Solar Spectrum shape source
  
  !CHECK CLOUDS
  
  
  CLOUDCONDITION :do i =1,mccx
     if ( fi%fc(i)%cldfrac < 0 .or. fi%fc(i)%cldfrac > 1.0001 ) fi%ierror=119 ! cloud fraction out of range
     
     if ( fi%fc(i)%cldfrac > 0 .and. fi%fc(i)%cldfrac <= 1.0001 .and. .not. fi%fc(i)%dpi%ldpi) then
        if ( fi%fc(i)%novl > mxovl .or. fi%fc(i)%novl < 1 ) fi%ierror=135 ! improper value of # of overlap conditions
        
        
        
        if ( fi%fc(i)%cldfrac == 0.0 ) cycle
        OVERLAP : do j=1,fi%fc(i)%novl
           if ( fi%fc(i)%iphase(j) < 1 .or. fi%fc(i)%iphase(j) > 2 ) fi%ierror=120 ! cloud phase out of range
           if ( fi%fc(i)%tau_vis(j) < 0 .or. fi%fc(i)%tau_vis(j) > 1000 ) fi%ierror=121 ! cloud tau_vis out of range
           if ( fi%fc(i)%iphase(j) == 1 .and. (fi%fc(i)%part_size(j) < 4 .or. fi%fc(i)%part_size(j) > 30) ) fi%ierror=122 ! WATER cloud part size out of range
           if ( fi%fc(i)%iphase(j) == 2 .and. (fi%fc(i)%part_size(j) < 20 .or. fi%fc(i)%part_size(j) > 300) ) fi%ierror=123 ! ICE cloud part size out of range
           
           if ( fi%fc(i)%icld_top(j) < 1 .or. fi%fc(i)%icld_bot(j) > fi%nv ) fi%ierror=124 ! Cloud layer placement error
           
           if ( fi%fc(i)%icld_top(j) > fi%fc(i)%icld_bot(j) ) fi%ierror=127 ! Cloud layer assignment order error
        enddo OVERLAP
        
     else
        if ( fi%fc(i)%cldfrac > 0 .and. fi%fc(i)%cldfrac <= 1.0001 .and. &
             fi%fc(i)%novl /= 0 .and. fi%fc(i)%dpi%ldpi  ) fi%ierror=136 ! improper value of # of overlap conditions
        
        
     endif
  enddo CLOUDCONDITION
  
  if ( sum(fi%fc(1:mccx)%cldfrac) > 1.0001)  fi%ierror=125 ! Total Cloud Fraction error 
  
  if ( fi%isksolve ==1 ) then
     do i =1,mccx
        if ( fi%fc(i)%cldfrac > 0 .and. .not. fi%fc(i)%dpi%ldpi ) then
           do j=1,fi%fc(i)%novl
              if ( fi%fc(i)%sc(j)%shapef <= 0 ) then
                 if ( fi%fc(i)%sc(j)%mn_lin_tau<= 0  ) fi%ierror=126 ! SUBSCALE Cloud Inputs ZERO
              endif
           enddo
        endif
     enddo
  endif

  if( fi%vo%ierr > 0 )then
     fi%ierror=137 !! fi%vo%ierr > 0 in VLA algorithm
     print*,fi%vo%ierr
     print*,'bad fi%vo%ierr, check_inputs_fu'
     stop
  endif

end subroutine check_inputs_fu
!=============================================================================================================
subroutine gwtsa_correction(icldcnd)
! "taucorrection" is the d(lntau) needed to correct the cloud optical depth
! for GWTSA algorithm based on biases computed using independent realizations
! of a gamma distribution of tau using the homogeneous solver for various
! Mean Tau, Nu, # of Cloud layers , Cos Sol 
! Corrections are largest for many layer clouds at low NU for mean tau in range 20:80.

USE FU_NU_TAU_NCL, only: taucorrect
implicit none
integer  icldcnd,ncl,ibx
real tau_vis_use,sfcalb_aprox
integer j 
real w_aprox(15) 
data w_aprox /0.000 , 0.000  ,0.000  ,0.000  ,0.002 ,0.027 , 0.099 , 0.114 , 0.169  ,0.145 &
            ,0.377 ,0.055 , 0.012 ,0.000 , 0.000 /


OVERLAP: do  j = 1,fi%fc(icldcnd)%novl
!----
sfcalb_aprox = sum( w_aprox(1:15)* fi%sfcalb ( 1:15 , 1 , icldcnd) ) 

tau_vis_use =  fi%fc(icldcnd)%tau_vis(j)
if ( fi%fc(icldcnd)%sc(j)%mn_lin_tau >0 ) tau_vis_use = fi%fc(icldcnd)%sc(j)%mn_lin_tau

ncl =  1+( fi%fc(icldcnd)%icld_bot(j) - fi%fc(icldcnd)%icld_top(j) )

taucorrection(icldcnd,j) = 0.0
if ( shapefactor(icldcnd,j) > 0 .and. shapefactor(icldcnd,j) < 100. ) &
taucorrection(icldcnd,j) = taucorrect(tau_vis_use,ncl,shapefactor(icldcnd,j),fi%u0,sfcalb_aprox )

! print'(a8,2i5,6f9.3)','TAUCORR',icldcnd,j,taucorrection(icldcnd,j), &
!  tau_vis_use,real(ncl),shapefactor(icldcnd,j),fi%u0,sfcalb_aprox
enddo OVERLAP
 
end subroutine gwtsa_correction
!=============================================================================================================
subroutine cloud_input_dpi(icldcnd,ib)
  !INITALIZE TO ZERO....
  integer icldcnd,ib
  plwc= 0.0  ; pre= 0.0  ;piwc= 0.0 ;pde= 0.0;pdge= 0.0 !! LAYERS
  
  where ( fi%fc(icldcnd)%dpi%phase(1:fi%nv) == 1 ) !! WATER PHASE
     plwc(1:fi%nv) = fi%fc(icldcnd)%dpi%plwc_iwc(1:fi%nv )
     pre(1:fi%nv) = fi%fc(icldcnd)%dpi%pre_de (1:fi%nv )
  endwhere
  
  where ( fi%fc(icldcnd)%dpi%phase(1:fi%nv) == 2 ) !!  !ICE PHASE
     piwc(1:fi%nv) = fi%fc(icldcnd)%dpi%plwc_iwc(1:fi%nv)
     pde(1:fi%nv) = fi%fc(icldcnd)%dpi%pre_de(1:fi%nv)
     pdge(1:fi%nv) =  -2.4+ 0.7566*fi%fc(icldcnd)%dpi%pre_de(1:fi%nv) + 9.22E-04* fi%fc(icldcnd)%dpi%pre_de(1:fi%nv)  ** 2
  endwhere



end subroutine cloud_input_dpi
!=============================================================================================================

subroutine cloud_input(icldcnd,ib)

  USE EXTRAS,only :tau_wp

  !USE FUINPUT ,only : &
  !fi ,nvx,dz, & !!! INPUT
  !plwc,pre,piwc,pde,pdge !!! OUTPUT
  
  implicit none
  integer, parameter:: idirection=1
  integer icldcnd,ib,j
  integer k,kb,ke
  
  real clc,dzcld,clwp
  
  real tau_vis_use
  real clc_check
  !INITALIZE TO ZERO....
  plwc= 0.0  ; pre= 0.0  ;piwc= 0.0 ;pde= 0.0;pdge= 0.0 !! LAYERS
  
  
  OVERLAP : do j = 1,fi%fc(icldcnd)%novl  !! OVERLAP CASES 
     
     kb= fi%fc(icldcnd)%icld_top(j)
     ke= fi%fc(icldcnd)%icld_bot(j)
     dzcld = sum( dz(kb:ke) ) *1000.
     
     
     if ( ib <= 6 .and. fi%isksolve==1 .and. fi%fc(icldcnd)%sc(j)%mn_lin_tau >0 ) then
        tau_vis_use = exp( log( fi%fc(icldcnd)%sc(j)%mn_lin_tau ) + taucorrection(icldcnd,j) )  
     else
        tau_vis_use = fi%fc(icldcnd)%tau_vis(j)
     endif
     
     
     
     clwp  = tau_wp(idirection,  &
          tau_vis_use,	  &
          fi%fc(icldcnd)%part_size(j),  &
          fi%fc(icldcnd)%iphase(j))
     
     
     
     
     
     !print*, 'Cloud Top   =',kb
     !print*, 'Cloud Bot   =',ke
     !print*, 'Cloud Pres  =',fi%pp(kb)),fi%pp(ke+1)
     !print*, 'Cloud DZ    =',dzcld
     !print*, 'Cloud TAU   =',fi%fc(icldcnd)%tau_vis(j)
     !print*, 'Cloud RE/De =',fi%fc(icldcnd)%part_size(j)
     !print*, 'Cloud Phase =',fi%fc(icldcnd)%iphase(j)
     !print*, 'Cloud LWP   =',clwp
     !print*, 'Cloud LWC   =',clc
     
     
     
     
     if ( clwp > 0) then 
        if ( fi%wp_hgt_flag == 0 )clc  = clwp/dzcld
        
        clc_check=0
        do k = kb,ke
           
           if(fi%wp_hgt_flag == 1 ) then 
              clc  = linear_clc(kb,ke,k) * clwp
           endif
           
           clc_check = clc_check + clc*dz(k)
           
           if ( fi%fc(icldcnd)%iphase(j) == 1 ) then ! WATER PHASE
              plwc ( k ) = clc
              pre  ( k ) = fi%fc(icldcnd)%part_size(j)
              
           elseif (fi%fc(icldcnd)%iphase(j) == 2 )then !ICE PHASE
              piwc ( k ) = clc
              pde  ( k ) = fi%fc(icldcnd)%part_size(j)
              pdge ( k ) = &
                   -2.4+ 0.7566*fi%fc(icldcnd)%part_size(j) + 9.22E-04*fi%fc(icldcnd)%part_size(j) ** 2
           else	 
              fi%ierror = 50 ! print*,   ' BAD Cloud Phase in cloud_input'
           endif !Phase
           
        enddo
        !print'(a10,3i5,f10.5)','clc_check',fi%wp_hgt_flag,icldcnd,j,clc_check
     endif !!clwp >0
     
  enddo OVERLAP



end    subroutine cloud_input
!===========================================================
real function linear_clc(kb,ke,k1)
implicit none
integer kb,ke,k,k1
real z0,z1,zbar(nvx),zsum

z0=0
zsum=0
do k = ke,kb,-1
z1=z0+dz(k)
zbar(k)=((z0+z1)*0.5)*dz(k)
zsum=zsum+zbar(k)
z0=z1
enddo

zbar(kb:ke) = zbar(kb:ke)/zsum
if ( sum(zbar(kb:ke)) > 1.001 ) stop ' Bad Sum linear_clc'

!do k = ke,kb,-1
!zbar(k)= zbar(k)*   dz(k)/sum( dz(kb:ke) )
!enddo
linear_clc = zbar(k1) / (dz(k1)*1000.)

!print*, 'linear_clc',sum(zbar(kb:ke))
!print*,k1,linear_clc
end function linear_clc

END MODULE FUINPUT

