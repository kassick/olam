MODULE testcases_jw

  !=======================================================================
  !
  !  Functions for setting up initial conditions for the Jablonowski-Williamson test case.
  !
  !  Given longitude (radians), latitude (radians), eta (pressure) and rotation_angle (degrees)
  !  the functions will return temperature, surface geopotential, zonal and meridional wind
  !  components, respectively.
  !
  !  lperturb=.FALSE. result in initial conditions for the steady-state test case.
  !  lperturb=.TRUE.  result in initial conditions for the baroclinic wave test case.
  !
  !     T   : FUNCTION temperature         (lon,lat,eta,rotation_angle)
  !     PHIS: FUNCTION surface_geopotential(lon,lat,rotation_angle)
  !     U   : FUNCTION u_wind              (lon,lat,eta,lperturb,rotation_angle)
  !     V   : FUNCTION v_wind              (lon,lat,eta,lperturb,rotation_angle)
  !     PS  : set to the constant p0
  !
  !  The non-rotated (rotation_angle=0) version of the test cases is described in: 
  !
  !                 Jablonowski, C., and D. L. Williamson, 2006: A baroclinic instability 
  !                 test case for atmospheric model dynamical cores. 
  !                 Quart. J. Roy. Meteor. Soc., 132, 2943-2975.
  !
  !                 Jablonowski, C., and D. L. Williamson, 2006: A Baroclinic Wave Test Case 
  !                 for Dynamical Cores of General Circulation Models: Model Intercomparisons, 
  !                 NCAR Technical Note, NCAR/TN-469+STR, 89 pp. 
  !  
  !  The rotated version simply rotates the initial conditions so that the spherical coordinate
  !  poles do not conicide with the earth's rotation axis. Thereby the Coriolis parameter is
  !  a function of latitude and longitude:
  !
  !      f = 2*Omega*(-cos(lon)*cos(lat)*sin(rotation_angle)+sin(lat)*cos(rotation_angle))
  !
  !  where Omega = 7.292 x 10E-5/s and rotation_angle is the angle between the flow direction
  !  and equator.
  !
  !  Author: Peter Hjort Lauritzen (NCAR, pel@ucar.edu)
  !          Christiane Jablonowski (University of Michigan, cjablono@umich.edu)
  !
  !=======================================================================

CONTAINS

  SUBROUTINE steady_state_baroclinic_wave
!=======================================================================
!     CAM3 grid and resolutions
!=======================================================================
      use cam3_grid
!=======================================================================
!     CAM3 variables
!=======================================================================
      use cam3_variables
      
      IMPLICIT NONE

      integer i, j, k
      real(r8) :: eta_c
!-----------------------------------------------------------------------
!    initialize the wind components
!-----------------------------------------------------------------------
      if (model_version.eq.1) then
!==========================================
!         EUL or SLD unstaggered grid
!==========================================
         do i=1,nlon
            do j = 1, nlat
               do k = 1, nlev
                  u(i,k,j) = u_wind(lon(i),lat(j),eta(k),choice_testcase == 2,rotation_angle)
                  v(i,k,j) = v_wind(lon(i),lat(j),eta(k),choice_testcase == 2,rotation_angle)
               enddo
            enddo
         enddo
      else
         
         !==========================================
         !         FV staggered grid
         !==========================================
         
         do k = 1, nlev
            do j = 1, nslat
               do i = 1, nlon
                  us(i,k,j) = u_wind(lon(i),slat(j),eta(k),choice_testcase == 2,rotation_angle)
               enddo
            enddo
         enddo
         do k = 1, nlev
            do j = 1, nlat
               do i = 1, nlon
                  vs(i,k,j) = v_wind(slon(i),lat(j),eta(k),choice_testcase == 2,rotation_angle)                  
               enddo
            enddo
         enddo
      endif

!-----------------------------------------------------------------------
!    initialize T (temperature)
!-----------------------------------------------------------------------
      do i=1,nlon
         do j = 1, nlat
            do k = 1, nlev
                t(i,k,j)   = temperature(lon(i),lat(j),eta(k),rotation_angle)
            enddo
         enddo
      enddo

!-----------------------------------------------------------------------
!    initialize PHIS (surface geopotential)
!-----------------------------------------------------------------------
      do j = 1, nlat
         do i=1,nlon
            phis(i,j) = surface_geopotential(lon(i),lat(j),rotation_angle) 
         enddo
      enddo

!-----------------------------------------------------------------------
!    initialize PS (surface pressure)
!-----------------------------------------------------------------------
     ps = p0  ! constant

!-----------------------------------------------------------------------
!      tracer q (specific humidity), in CAM automatically transported, therefore initialized
!-----------------------------------------------------------------------
     do i=1,nlon
       do j = 1, nlat 
         do k = 1, nlev
           q(i,k,j) = tracer_q(lon(i),lat(j),eta(k),rotation_angle)
         enddo
       enddo
     enddo

!-----------------------------------------------------------------------
!  Tracer variables in CAM3.5: 'TT_LW', 'TT_MD', 'TT_HI', 'TTRMD'
!-----------------------------------------------------------------------
     if (tracer_nr(1:1) == '0') then
       tc_variant = '0'
     else
!-----------------------------------------------------------------------
!     tracer q1
!-----------------------------------------------------------------------
      if (tracer_nr(1:1) == '1' .or. tracer_nr(2:2) == '1' .or. tracer_nr(3:3) == '1' .or. &
          tracer_nr(4:4) == '1' .or. tracer_nr(5:5) == '1') then
!      k = index(tracer_nr, '1', 'false')
        tc_variant = trim(tc_variant)//'1'             ! tracer q1 is test variant y=1
        eta_c = 0.6_r8                                 ! vertical position of the center of the tracer
        do i=1,nlon
           do j = 1, nlat
              do k = 1, nlev
                  tt_lw(i,k,j) = tracer_q1_q2(lon(i),lat(j),eta(k),rotation_angle, eta_c)
              enddo
           enddo
        enddo
      endif
!-----------------------------------------------------------------------
!     tracer q2
!-----------------------------------------------------------------------
      if (tracer_nr(1:1) == '2' .or. tracer_nr(2:2) == '2' .or. tracer_nr(3:3) == '2' .or. &
          tracer_nr(4:4) == '2' .or. tracer_nr(5:5) == '2') then
        tc_variant = trim(tc_variant)//'2'            ! tracer q2 is test variant y=2
        eta_c = 1._r8                                 ! centered at the surface
        do i=1,nlon
           do j = 1, nlat
              do k = 1, nlev
                  tt_md(i,k,j)   = tracer_q1_q2(lon(i),lat(j),eta(k),rotation_angle, eta_c)
              enddo
           enddo
        enddo
      endif
!-----------------------------------------------------------------------
!     tracer q3
!-----------------------------------------------------------------------
      if (tracer_nr(1:1) == '3' .or. tracer_nr(2:2) == '3' .or. tracer_nr(3:3) == '3' .or. &
          tracer_nr(4:4) == '3' .or. tracer_nr(5:5) == '3') then
        tc_variant = trim(tc_variant)//'3'            ! tracer q3 is test variant y=3
        do i=1,nlon
           do j = 1, nlat
              do k = 1, nlev
                  tt_hi(i,k,j)   = tracer_q3(lon(i),lat(j),eta(k),rotation_angle)
              enddo
           enddo
        enddo
      endif
!-----------------------------------------------------------------------
!     tracer q4
!-----------------------------------------------------------------------
      if (tracer_nr(1:1) == '4' .or. tracer_nr(2:2) == '4' .or. tracer_nr(3:3) == '4' .or. &
          tracer_nr(4:4) == '4' .or. tracer_nr(5:5) == '4') then
        tc_variant = trim(tc_variant)//'4'             ! tracer q4 is test variant y=4
        ttrmd = 1._r8                                  ! constant
      endif
    endif

  end subroutine steady_state_baroclinic_wave



!-----------------------------------------------------------------------
! Tracer q, absolute value of the relative vorticity of the unperturbed initial state
!           multiplied by 10^5
!-----------------------------------------------------------------------
  REAL(r8) FUNCTION tracer_q(lon,lat,eta,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta, lon, lat, rotation_angle
    REAL(r8) :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF
    tracer_q = abs(-4._r8 * u0/a * (cos((eta-eta0)*pi*0.5_r8))**1.5_r8 * sin(rot_lat) * &
               cos(rot_lat) * (2._r8-5._r8*(sin(rot_lat))**2)) * 1.e5_r8
    if (tracer_q < 1.e-9_r8) tracer_q = 0._r8  !  otherwise error in netcdf file

  END FUNCTION tracer_q

END MODULE testcases_jw
