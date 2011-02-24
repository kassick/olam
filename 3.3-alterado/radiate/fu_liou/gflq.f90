MODULE GENERATE_FULIOU_LEVELS
  USE FUINPUT,only :  fi,mfix
  USE EXTRAS, only : new_sfc_alt,z1_p1
  implicit none
  
  type gflq_type
     character(len=5) :: mode
     real :: hsfc   ! INPUT :: Surface Geopotential Height
     integer :: nlev  !! OUTPUT Number of levels
  end type gflq_type
  
  type ( gflq_type ) gflq
  
  
  real hh(mfix) !! OUTPUT Height levels
  real pp(mfix) !! OUTPUT Pressure levels
  


CONTAINS
  !==============================================
  subroutine generate_level_scheme
    real, parameter :: hmax = 65000 ! Max Height Meters 
    real, parameter :: sfcmindz = 1.0 ! Minimum thickness of surface layer.
    integer i
    real z0,p0,zx,dzz,px,sh, z0_moa
    real z0a
    real h0p,p0p,t0p,q0p,o0p
    
    ! Force levels to even Multiple of 120 meters as to not change vertical level alignment.
    h0p= FI%VI%hsfc
    p0p= FI%VI%pp(FI%VI%nlev)
    t0p= FI%VI%pt(FI%VI%nlev)
    q0p= FI%VI%ph(FI%VI%nlev)
    o0p= FI%VI%po(FI%VI%nlev)
    




    z0a =gflq%hsfc
    if ( z0a < -480. ) z0a =  -480.0
    fi%hsfc = z0a !!! ASSIGN FI%HSFC .. Sfc HEIGHT if generate_level_scheme called
    
    
    z0 = NINT(z0a/250) *250.0
    if ( z0a <  3000.) z0 = NINT(z0a/120) *120.0
    if ( z0a < -480. ) z0 =  -480.0
    
    
    
    
    if( z0a < h0p ) then  !; print*, ' SURFACE below Input sounding.Adjusting'
       
       call new_sfc_alt (  h0p, & !! INPUT PROFILE SURFACE  Geopotential
            p0p  , t0p  , q0p, &
            z0a , & ! I ! 1 meter below to make sure 
            FI%VI%pp(FI%VI%nlev+1), FI%VI%pt(FI%VI%nlev+1), FI%VI%ph(FI%VI%nlev+1) ) 
       
       FI%VI%po(FI%VI%nlev+1) = o0p 
       FI%VI%nlev = FI%VI%nlev +1 !! ADD A Level to FI%VI Profile 
       
       i = FI%VI%nlev ! use augmented profile 
       p0 = z1_p1 ( FI%VI%pp(i), FI%VI%pt(i) , z0a,    z0)
    else ! Find Surface pressure @ height Z0
       p0 =   z1_p1 ( p0p, t0p , h0p,    z0)
    endif
    
    
    zx = z0
    px = p0  !surface pressure
    MAKE_LEVELS : do  i=1,mfix
       
       if (gflq%mode == 'CALIP' ) then
          if (  zx <  3000 ) then
             dzz= 120 !!! 120Meter resolution  Below 3km
          elseif (  zx < 18000 .and. zx >=  3000) then
             dzz= 250 !!! 250Meter resolution 3-18km
          elseif (  zx < 30000 .and. zx >= 18000) then
             dzz= 500 !!! 500Meter resolution 18-30km
          elseif (  zx < 40000 .and. zx >= 30000) then 
             dzz= 1000 !!!1000Meter resolution 30-40km
          elseif (  zx < 50000 .and. zx >= 40000) then
             dzz= 2000 !!!2000Meter resolution 40-50km
          else
             dzz= 5000.		    !!!5000Meter resolution Above 50km
          endif
          
          !--------------------------------------------
       elseif (gflq%mode == 'CERES' ) then  
          if (  zx <  1000 ) then
             dzz= 250 !!! 250Meter resolution  Below 1km
          elseif (  zx < 3000 .and. zx >= 1000) then
             dzz= 500 !!! 500Meter resolution 1-3km
          elseif (  zx < 10000 .and. zx >= 3000) then
             dzz= 1000 !!! 1000Meter resolution 3-10km
             
          elseif (  zx < 30000 .and. zx >= 10000) then 
             dzz= 2500 !!!2500Meter resolution 10-30km
          else
             dzz= 5000.  !!!5000Meter resolution Above 30km
          endif
          
          !--------------------------------------------
       else
          stop ' gflq%mode'
       endif
       
       
       
       sh = get_scalehgt(px) ! Scale height in Meters
       !print'(i4,5f10.2)',i,zx,dzz,px,sh
       gflq%nlev = i
       pp(i) = px
       hh(i) = zx
       if ( zx >= hmax ) exit
       
       zx = zx+dzz
       px = pa(px,dzz,sh)
       
    enddo MAKE_LEVELS
    !------------------------------------------------------------------------
    pp(1:gflq%nlev) = pp(gflq%nlev:1:-1) !! Make Toa @level 1
    hh(1:gflq%nlev) = hh(gflq%nlev:1:-1)
    
    !do i = 1,gflq%nlev ; print'(i4,f10.0 ,f10.2,f10.0)',i,hh(i),pp(i), hh(i)-hh(i+1) ;enddo
    
    !Fixed Level Definition
    FI%VD%nfix = gflq%nlev
    FI%VD%pfix(1:FI%VD%nfix) = pp(1:gflq%nlev  )
    
    !Floating Level Definition 
    FI%VD%nflo = 0
    
    if ( abs( z0a - z0) > sfcmindz ) then  ! ADD SURFACE LAYER 
       FI%VD%nflo = 1
       FI%VD%pflo(FI%VD%nflo) = z1_p1 ( p0p, t0p , h0p,    z0a)
       !  print'(a10,4f8.2)','SFCFLO',z0a,FI%VD%pflo(FI%VD%nflo)
       if ( z0a > z0 ) FI%VD%nfix =  FI%VD%nfix -1 ! Remove a fixed layer
    endif
    
    !----------
    
    if ( gflq%mode == 'CALIP') then
       !Reporting level Definition
       FI%VD%nrep = FI%VD%nfix
       FI%VD%report(1:FI%VD%nfix ) = FI%VD%pfix(1:FI%VD%nfix )
       if ( FI%VD%nflo == 1 ) then ! Include Surface level as reporting level
          FI%VD%nrep = FI%VD%nrep+1
          FI%VD%report(FI%VD%nrep) =  FI%VD%pflo(FI%VD%nflo)
          hh(FI%VD%nrep ) = z0a
       endif
    endif
    
    if ( gflq%mode == 'CERES') then
       FI%VD%nflo = 4
       FI%VD%pflo(1: FI%VD%nflo)   = (/ z1_p1(p0p, t0p, h0p, z0a),500.,200.,70./)  ! decreasing in pressure
       FI%VD%nrep = 5
       FI%VD%report(1:FI%VD%nrep ) = (/ FI%VD%pfix(1), 70., 200.,500.,z1_p1(p0p, t0p, h0p, z0a)/)
       if ( z1_p1(p0p, t0p, h0p, z0a) ==  FI%VD%pfix(FI%VD%nfix ) ) FI%VD%nfix= FI%VD%nfix -1
    endif
    
    FI%VI%hsfc = z0a !! adjust definition of surface elevation in case this code is looped
    
    !print'(10f8.2)',FI%VD%pfix(1:FI%VD%nfix )
    !print'(10f8.2)',FI%VD%report(1:FI%VD%nfix )
    !print'(a10,4f8.2)','SFCFIX',z0,p0,FI%VD%pfix(FI%VD%nfix )
    
    !do i = 1,FI%VD%nrep; print'(i4,2f10.2)',i,FI%VD%report(i),hh(i) ; enddo
    
  end subroutine generate_level_scheme
  !==================================================================
  real function get_scalehgt(px)
    ! determine scale height of atm at input pressure level.
    real px,tmp(FI%VI%nlev),perr
    integer ix(1)
    
    tmp = abs( FI%VI%pp(1:FI%VI%nlev)-px )
    ix(1:1) = minloc(tmp)
    get_scalehgt = 29.2674 * FI%VI%pt(ix(1)) * &  !pt= temperature in (K)
         (1.0 +0.61* FI%VI%ph(ix(1)) )    !ph = H20 mixration in (g/g) 
    !print*, px,ix(1),FI%VI%pp(ix(1))
  end function get_scalehgt
  !-----------------------------------------------------------------
  real function  pa(pb,dzz,sh) 
    real pb,dzz,sh
    pa =pb*exp(-dzz/sh)  
  end function pa
  !-----------------------------------------------------------------
  real  function dz(p1,p2,sh)
    real p1,p2,sh
    dz = sh*log(p1/p2)
  end function dz
  
END  MODULE GENERATE_FULIOU_LEVELS
