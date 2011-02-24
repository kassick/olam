MODULE FULIOUMULTI
  use FUINPUT ,only : fi , set_default_options_fu, check_inputs_fu
  use FUOUTPUT,only : fo ,focc,fos,foscc ,ftoa,fsfc , fouv
  use  VLA_FU ,only : prepare_model_profile_fu,vla_interface_fu, &
       print_vla_in,print_vla_out,vo
  
contains
  !=========================================================================
  subroutine rad_multi_fu(iw)
    
    USE FUINPUT
    USE FUOUTPUT
    USE WSU_LUT_INTERFACE , only :  uvco,wsu_lut
    USE ICEDIRSFC ,only :dircorrect    
    use misc_coms, only: ilwrtyp, iswrtyp
    
    implicit none
    
    integer, intent(in) :: iw
    integer kg(mbx+2)
    integer i,j,ib,ig,mbn,iac,icp,icc,iccb,icce,idum
    real f0,fuq1,xx,asx
    real hk
    integer icmf
    integer ibt
    real eeband
    real p_mle_seiji2
    integer :: iband1, iband2
    
    real ffdr,ffdf
    common /dirdiff/ffdr(nv1x),ffdf(nv1x)

    real fiur
    common /radiance/ fiur(nv1x) 
    real hk1, fk1o3,sol_spect,fk1h2o
    common /band1/ hk1(10), fk1o3(10),sol_spect(0:7),fk1h2o(10)

    !  ******************************************************************
    !  kg(mb+2) is the number of intervals to perform the g-quadrature in
    !  each band to consider the nongray gaseous absorption.  In total,
    !  we need to perform "sum(kg)" spectral calculations in  the  scattering
    !  problem for each atmospheric profile.
    !  ******************************************************************
    data kg / 10, 8, 12, 7, 12, 5, &
         2, 3, 4, 4, 3, 5, 2, 10, 12, 7, 7, 8, &
         5,5 /

    call check_inputs_fu !! resets error flag!
    
    if (fi%ierror .ne. 0 )then
       print*,'ABORTING.  Bad Fu-Liou inputs.',fi%ierror
       stop
    endif
    
    call assigninputs !!! PULL GLOBAL VARIABLES FROM INPUT STRUCTURE

    !! Depending on K's option # of Solver loops needed changes...
    if ( idkfr == 0) kg(11:13) = (/3,5,2/)
    if ( idkfr == 1) kg(11:13) = (/80,40,10/)
    if ( idkfr == 2) kg(11:13) = (/5,5,2/)

    if ( isksw > 0)  kg(2:6) = (/7,8,7,8,7/) !ht02 sav

    !  *********************************************************************
    !  The following variables are declared here, instead of in the calling
    !  program, since any change in these would require modifications
    !  the code anyway.  A check is inserted here, to make sure the 
    !  number of layers given by the user, nv, is less than or equal
    !  to the number of levels nvx, given in the header file rad_0598.h.
    !  nvx is used for dimensioning arrays.  Uncomment to use.
    !  *********************************************************************
    nv1    = nv+1      ! Number of levels = LAYERS +1
    mb     = 18        ! Number of bands
    mbs    = 6         ! number of shortwave bands
    mbir   = 12        ! number of longwave bands
    nc     = 8         ! number of cloud types
    ndfs   = nv        ! number of layers
    mdfs   = nv1       ! number of levels
    ndfs4  = 4*ndfs    ! number of layers * 4
    ndfs2  = 2*ndfs    ! number of layers * 2
    f0     = 1.0 / 3.14159

    call flx_initalize

    call thicks

    call rayle2

    if ( u0 .le. 1.0e-4 ) then
       mbn = mbs + 1    ! Longwave only if sun is too low
    else
       mbn = 1          ! Shortwave and longwave
    endif

    call aerosol_init

    call ckd1_init( isolar_spectrum,lray )

    if ( lband6a ) then
       fuq1 = ss /  sol_spect(0)  
    else
       fuq1 = ss / ( sol_spect(0) - sol_spect(7) )
    endif



    if ( fi%isksolve == 1  ) then
       CLOUD_CONDITION00: do icc = 1,mccx
          if ( fi%fc(icc)%cldfrac > 0.0 ) then

             OVERLAP00: do j=1,fi%fc(icc)%novl
                if( fi%fc(icc)%sc(j)%shapef > 0 ) then 
                   shapefactor(icc,j) = fi%fc(icc)%sc(j)%shapef 
                else
                   shapefactor(icc,j)=	    &
                        p_mle_seiji2( fi%fc(icc)%sc(j)%mn_lin_tau ,&
                        log(fi%fc(icc)%tau_vis(j) ) )
                   if ( shapefactor(icc,j) > 1E+3)  shapefactor(icc,j) = 1E+3
                endif

             enddo OVERLAP00

             call gwtsa_correction(icc)

          endif
       enddo CLOUD_CONDITION00

    endif

    iband1 = mbn
    iband2 = mb + nhb
    if(iswrtyp /= 4)iband1 = mbs + 1
    if(ilwrtyp /= 4)iband2 = mbs

! SW only
!    BAND_LOOP : do  ibt = mbn,mbs
! All
!    BAND_LOOP : do  ibt = mbn,mb+nhb
! LW only
!    BAND_LOOP : do  ibt = mbs+1,mb+nhb
    BAND_LOOP : do  ibt = iband1, iband2

       !! WATCH CAREFULLY the use of ib and ibt
       ! ibt (1-20) band index over added hidden bands 
       ! ib (1-18) band index over reported bands
       ib=ibt
       if ( ibt == 19  ) ib =7 !!!! (19)2850:2500cm-1 
       if ( ibt == 20  ) ib =7 !!!! (20)2500:2200cm-1

       if (ib.ne.1) call aerosolxy (ib,'x')

       if (irobckd .eq. 0 ) then
          call gascon_mt_parm ( ib )  !! Ma and Tipping (2002)
       elseif (irobckd .eq.1) then !'Set Continuum Model 0= Ma_Tip 1=Roberts 2=CKD_2.1 3=NONE, 4=PARM CKD_2.1 5=PARM CKD_2.4'
          call gascon ( ib )  !! OLD ROBERTS CONTINUUM
       elseif (irobckd .eq.2) then
          call gascon_ckd ( ib ) !! EXACT CKD_2.1 (SLOW)
       elseif (irobckd .eq.3 ) then
          call gascon_off      !! Option to Turn off Continuum Absorbtion
       elseif (irobckd .gt. 3 ) then
          call gascon_ckd_parm(ib) !! Parameterized CKD_2.1 Cont. Abs. OR CKD_2.4
       endif

       if ( ibt >18) 	 call gascon_off 		
       if ( ibt .gt. mbs )  call planck ( ibt )


       !-------------------------------------------------------------
       K_LOOP: do  ig = 1, kg(ibt)

          !----------------
          !AEROSOLS
          if (ib.eq.1) call aerosolxy (ig,'y')

          !----------------
          !RAYLEIGH
          call rayle ( ib, ig )

          !----------------
          !GASES        
          if ( isksw > 0 .and. ( (ibt>=2 .and. ibt<=6) &
               .or. (ibt==1 .and. ig >=9) )  ) then
             call seijifu_ht02a_sav( ibt, ig, hk )
          else
             call gases ( ibt, ig, hk )
          endif

          !----------------
          !CLOUDS
          if (fi%lscm(2) .or. fi%lscm(4) ) then

             CLOUD_CONDITION1: do icc = 1,mccx
                if ( fi%fc(icc)%cldfrac > 0.0 ) then 

                   if(fi%fc(icc)%dpi%ldpi) then
                      call cloud_input_dpi(icc,ib)
                   else	
                      call cloud_input(icc,ib)
                   endif

                   !if ( ib == 1 .and. ig == 1 .or. ib > 1 ) then   !band
                   if (  ig == 1  ) then   !once per band
                      if ( fl93i ) then
                         call ice ( ib,icc )    !	fl93i = .true.  ! True = Old93 ice
                      else
                         call icenew ( ib,icc) ! 	fl93i = .false. !! False = New 98Ice cld 
                      endif

                      ! Use Yong Hu's Water cloud optical properties for SW bands
                      if ( ib <=6 ) then
                         call water_hu ( ib,icc ) ! CALL Yong Hu's WATER CLOUD OPTICS For SW bands ib=1:6
                      else
                         call water ( ib,icc )
                      endif

                   endif ! once per band
                endif ! CloudFrac >0
             enddo CLOUD_CONDITION1
          endif

          call comscp_0203(0,idum) ! 0=construct,
          ! 1=Clear,2=Total,3=Pristine,4=NoaerCLoud 
          
          call solver_configuration(ib,ibt,ig,f0,hk,fuq1,iw)
          
          
          !print*,'BAND',ib,ibt,ig
       enddo K_LOOP
    enddo BAND_LOOP


    call dircorrect 

    !--------------------------------------------------------------------
    ! Combine Cloud Condition Componets into Total Sky..	
    !--------------------------------------------------------------------

    call aux_flx_calc(0,1)	     ! Heat rate ,Net, Window emmulation ...
    fo(1) = focc(0,1)	     !!Clear
    fos(1) =foscc(0,1)
    
    call flx_cloudcombine(2,1,iw)  ! CLOUDY With Aerosols
    
    call aux_flx_calc(0,2)
    fo(3) = focc(0,2)	     !!Pristine 
    fos(3) =foscc(0,2)
    
    call flx_cloudcombine(4,2,iw)  ! CLOUDY NO Aerosols
    
!    call only_toa_sfc
    
!    call wsu_lut
    
    return

end subroutine rad_multi_fu
!================================================================================

!================================================================================
subroutine solver_configuration(ib,ibt,ig,f0,hk,fuq1,iw)

  USE FUINPUT
  USE FUOUTPUT
  implicit none
  integer ib,ibt,ig
  real f0,hk,fuq1
  real asx,eeband
  integer icmf,iccb,icce
  integer icc,icp 
  integer ibseq
  integer cc2_flag
  integer, intent(in) :: iw
  SOLVER_CONFIG : do icmf = 1,nscm

     if (icmf == 1 .or. icmf == 2 ) icp=1
     if (icmf == 3 .or. icmf == 4 ) icp=2

     if (fi%lscm(icmf) ) then 

        if (  icmf == 1 .or. icmf == 3 ) then
           iccb = 0   ;  icce = 0
        else 
           iccb = 1   ;  icce = mccx
        endif

        CLOUD_CONDITION2 : do icc=iccb,icce
           
           if(icc == 0)then
              cc2_flag = 1
           elseif(fi%fc(icc)%cldfrac > 0.0)then
              cc2_flag = 1
           else
              cc2_flag = 0
           endif
           if ( cc2_flag == 1 ) then 
              
              call comscp_0203(icmf,icc )  

              if ( ib .le. mbs ) then  !! SHORTWAVE CALLS

                 if ( ib == 1 ) then
                    ibseq=ig
                 elseif ( ib >= 2 .and. ib <=6) then
                    ibseq=ib+9
                 endif

                 asx = fi%sfcalb(ibseq,icp,icc)

                 if ( fourssl ) then
                    call qfts ( ib, asx, f0 )  !! FOUR STREAM SOLIR	
                 else
                    quadra = .false. ; hemisp = .false. ; edding = .true.  !! Original
                    !  quadra = .true. ; hemisp = .false. ; edding = .false.
                    if (isksolve == 0)   call qftsts ( ib, asx, f0 ) ! TWO STREAM SW FULIOU

                    if (isksolve == 1) call rapradgw_qftsts(asx,icc ) ! SEIJI_KATO SOLVER


                    if (fi%ierror .NE. 0) return
                 endif

                 call sw_spectral_integration(icc,icp,ibseq,hk,fuq1)

              else !!!THERMAL
                 !-------------------------------------------------------------------
                 eeband =   ee(ib-mbs)
                 !!if (ibt == 19 .or. ibt== 20) eeband = 1.0-fi%sfcalb(15,1)
                 if (ibt == 19 .or. ibt== 20 ) eeband = ee(1) 
                 if ( foursir ) then
                    call qfti ( ib, eeband, iw ) !FOUR STREAM IR
                 else
                    quadra=.false.;edding=.false.;hemisp=.true.;mquadr=.false.
                    call qftisf ( ib, eeband, iw) !TWO-FOUR STREAM COMB IR
                    ! call qftits ( ib, eeband ) !TWO-STREAM IR
                    if (fi%ierror .NE. 0) return
                 endif
                 
                 call lw_spectral_integration(icc,icp,ib,ibt,ig,hk,iw)
                 
              endif !SOLAR OR THERMAL

           endif
        enddo CLOUD_CONDITION2

     endif
  enddo SOLVER_CONFIG

end subroutine solver_configuration


END MODULE FULIOUMULTI
