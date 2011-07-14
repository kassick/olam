subroutine fuliou_radinit()

  use aero_coms, only: srbqex, rhinfo, sruqex, refu22, srhqex, sruqsc,  &
       sruqcb, truqex, truqsc, truqcb, q55u22, qxdust, refd25, refs25,   &
       srdqex, srsqex, tdust

  use rastro_evts
  implicit none

  character(len=256), parameter :: mie_file =   &
       '/cluster/home/dmm31/aerosols/miescatpar.abcdv2'
  character(len=256), parameter :: dust_file =   &
       '/cluster/home/dmm31/aerosols/dust8.tau9x8x13'
  integer :: n
  character(len=80) :: title
  real, dimension(6,11) :: rdum
  real, dimension(33,11) :: rdum2
  real, dimension(6,15) :: rdum3
  real, dimension(33,15) :: rdum4
  real, dimension(6,25) :: rdum5
  real, dimension(33,25) :: rdum6
  real, dimension(6,20,5) :: rdum7
  real, dimension(33,20,5) :: rdum8
  real, dimension(6,120) :: rdum9
  integer :: k
  integer :: na
  ! TROPOSPHERIC AEROSOL COMPOSITIONAL/TYPE PARAMETERS
  ! SO4    SEA    ANT    OCX    BCI    BCB    DST   VOL
  real, dimension(8), parameter :: refdry = (/0.200, 1.000, 0.300, 0.300,   &
       0.100, 0.100, 1.000,1.000/)
  real :: areff
  integer :: m
  integer :: n1
  integer :: kk
  integer :: n2
  real, dimension(8), parameter ::  &
       !MINERAL DUST PARAMETERS
       !                CLAY                  SILT
       redust=(/ 0.1, 0.2, 0.4, 0.8,   1.0, 2.0, 4.0, 8.0/)

#ifdef OLAM_RASTRO
character*1 :: rst_buf = '_'
call rst_event_s_f(OLAM_FULIOU_RADINIT_IN,rst_buf)
#endif
 
  open(12,file=trim(mie_file),form='formatted',status='old')
 3000 FORMAT(A80)
 3001 FORMAT(18X,6(F7.5,1X))
 3002 FORMAT(18X,6(F7.5,1X)/18X,6(F7.5,1X))
 3003 FORMAT(18X,6(F7.3,1X)/18X,5(F7.3,1X))
 3004 FORMAT(14X,7(F7.5,1X),4(/14X,7(F7.5,1X)))
 3005 FORMAT(/14X,7(F7.5,1X),4(/14X,7(F7.5,1X)))

  do n=1,11
     read (12,'(a80)') title
     read (12,'(18x,6(f7.5,1x))') (rdum(k,n),k=1,6)
     read (12,'(18x,6(f7.5,1x))') (rdum(k,n),k=1,6)
     read (12,'(18x,6(f7.5,1x))') (rdum(k,n),k=1,6)
  enddo
  read (12,'(18x,6(f7.5,1x)/18x,6(f7.5,1x))') (rdum(1,n),n=1,11)
  read (12,'(18x,6(f7.3,1x)/18x,5(f7.3,1x))') (rdum(1,n),n=1,11)
  read (12,'(18x,6(f7.3,1x)/18x,5(f7.3,1x))') (rdum(1,n),n=1,11)
  do n=1,11
     read (12,'(a80)') title
     read (12,'(14x,7(f7.5,1x),4(/14x,7(f7.5,1x)))') (rdum2(k,n),k=1,33)
     read (12,'(/14x,7(f7.5,1x),4(/14x,7(f7.5,1x)))') (rdum2(k,n),k=1,33)
     read (12,'(/14x,7(f7.5,1x),4(/14x,7(f7.5,1x)))') (rdum2(k,n),k=1,33)
  enddo

  do n=1,10
     if(n /= 6)then
        read (12,'(a80)') title
        read (12,'(18x,6(f7.5,1x))') (srbqex(k,n),k=1,6)
        read (12,'(18x,6(f7.5,1x))') (rdum(k,n),k=1,6)
        read (12,'(18x,6(f7.5,1x))') (rdum(k,n),k=1,6)
     endif
  enddo
  read (12,3002) (rdum2(1,n),n=1,5),(rdum2(1,n),n=7,10)
  read (12,3003) (rdum2(1,n),n=1,5),(rdum2(1,n),n=7,10)
  read (12,3003) (rdum2(1,n),n=1,5),(rdum2(1,n),n=7,10)
  do n=1,10
     if(n /= 6)then
        read (12,3000) title
        read (12,3004) (rdum2(k,n),k=1,33)
        read (12,3005) (rdum2(k,n),k=1,33)
        read (12,3005) (rdum2(k,n),k=1,33)
     endif
  enddo


  !                               cloud water, ice-non, ice-mie parameters
  !                               ----------------------------------------
  do n=1,15
     read (12,'(a80)') title
     read (12,'(18x,6(f7.5,1x))') (rdum3(k,n),k=1,6)
     read (12,'(18x,6(f7.5,1x))') (rdum3(k,n),k=1,6)
     read (12,'(18x,6(f7.5,1x))') (rdum3(k,n),k=1,6)
  enddo

  read (12,3006) (rdum3(1,n),n=1,15)
3006 format( 18x,6(f7.5,1x)/18x,6(f7.5,1x)/18x,6(f7.5,1x))
  read (12,3007) (rdum3(1,n),n=1,15)
3007 format( 18x,6(f7.3,1x)/18x,6(f7.3,1x)/18x,6(f7.3,1x))
  read (12,3007) (rdum3(1,n),n=1,15)
  do n=1,15
     read (12,3000) title
     read (12,3004) (rdum4(k,n),k=1,33)
     read (12,3005) (rdum4(k,n),k=1,33)
     read (12,3005) (rdum4(k,n),k=1,33)
     read (12,3005) (rdum4(k,n),k=1,33)
  enddo

  !                               desert dust 25 sizes, mie parameter data
  !                               ----------------------------------------
  do n=1,25
     read (12,3001) (srdqex(k,n),k=1,6)
     read (12,3001) (rdum5(k,n),k=1,6)
     read (12,3001) (rdum5(k,n),k=1,6)
  enddo
  read (12,3008) (rdum5(1,n),n=1,25)
3008 format( 18x,5(f7.5,1x),4(/18x,5(f7.5,1x)))
  read (12,3009) (refd25(n),n=1,25)
3009 format( 18x,12(f3.1,1x)/18x,12(f3.1,1x)/18x,f3.0)
  read (12,3010) (rdum5(1,n),n=1,25)
3010 format( 18x,12(f3.1,1x)/18x,12(f3.1,1x)/18x,f3.1)
  do n=1,25
     read (12,3000) title
     read (12,3004) (rdum6(k,n),k=1,33)
     read (12,3005) (rdum6(k,n),k=1,33)
     read (12,3005) (rdum6(k,n),k=1,33)
     read (12,3005) (rdum6(k,n),k=1,33)
  enddo
  
  !volcanic aerosol mie size, variance data
  !----------------------------------------
  do m=1,5
     if(m /= 4)then
        do n=1,20
           read (12,3001) (rdum7(k,n,m),k=1,6)
           read (12,3001) (rdum7(k,n,m),k=1,6)
           read (12,3001) (rdum7(k,n,m),k=1,6)
        enddo
        read (12,3011) (rdum7(1,n,m),n=1,20)
3011    format( 18x,5(f7.5,1x),3(/18x,5(f7.5,1x)))
        read (12,3012) (rdum7(1,n,m),n=1,20)
3012    format( 18x,12(f3.1,1x)/18x,8(f3.1,1x))
        read (12,3012) (rdum7(1,n,m),n=1,20)
        do n=1,20
           read (12,3000) title
           read (12,3004) (rdum8(k,n,m),k=1,33)
           read (12,3005) (rdum8(k,n,m),k=1,33)
           read (12,3005) (rdum8(k,n,m),k=1,33)
           read (12,3005) (rdum8(k,n,m),k=1,33)
        enddo
     endif
  enddo

  !sulfate aerosol, mie parameter 22-size data
  !-------------------------------------------
  do n=1,22
     read (12,3000) title
     read (12,3001) (sruqex(k,n),k=1,6)
     read (12,3001) (sruqsc(k,n),k=1,6)
     read (12,3001) (sruqcb(k,n),k=1,6)
  enddo
  read (12,3008) (q55u22(n),n=1,22)
  read (12,3013) (refu22(n),n=1,22)
3013 format( 18x,5(f7.3,1x),4(/18x,5(f7.3,1x)))
  read (12,3013) (rdum5(1,n),n=1,22)
  do n=1,22
     read (12,3000) title
     read (12,3004) (truqex(k,n),k=1,33)
     read (12,3005) (truqsc(k,n),k=1,33)
     read (12,3005) (truqcb(k,n),k=1,33)
     read (12,3005) (rdum6(k,n),k=1,33)
  enddo

  !soot aerosol, mie parameter 25-size data
  !----------------------------------------
  do n=1,25
     read (12,3000) title
     read (12,3001) (srsqex(k,n),k=1,6)
     read (12,3001) (rdum5(k,n),k=1,6)
     read (12,3001) (rdum5(k,n),k=1,6)
  enddo
  read (12,3008) (rdum5(1,n),n=1,25)
  read (12,3013) (refs25(n),n=1,25)
  read (12,3013) (rdum5(1,n),n=1,25)
  do n=1,25
     read (12,3000) title
     read (12,3004) (rdum6(k,n),k=1,33)
     read (12,3005) (rdum6(k,n),k=1,33)
     read (12,3005) (rdum6(k,n),k=1,33)
     read (12,3005) (rdum6(k,n),k=1,33)
  enddo
  
  !Seasalt aerosol, Mie parameter 22-size data
  !Nitrate aerosol, Mie parameter 22-size data
  !(Water) aerosol, Mie parameter 22-size data
  !Organic aerosol, Mie parameter 22-size data
  !-------------------------------------------
  n1=23
  do kk=1,4
     n2=n1+21
     do n=n1,n2
        read (12,3000) title
        read (12,3001) (sruqex(k,n),k=1,6)
        read (12,3001) (sruqsc(k,n),k=1,6)
        read (12,3001) (sruqcb(k,n),k=1,6)
     enddo
     read (12,3008) (q55u22(n),n=n1,n2)
     read (12,3013) (refu22(n),n=n1,n2)
     read (12,3013) (rdum9(1,n),n=n1,n2)
     n1=n2+1
  enddo
  n1=23
  do kk=1,4
     n2=n1+21
     do n=n1,n2
        read (12,3000) title
        read (12,3004) (truqex(k,n),k=1,33)
        read (12,3005) (truqsc(k,n),k=1,33)
        read (12,3005) (truqcb(k,n),k=1,33)
        read (12,3005) (rdum6(k,1),k=1,33)
     enddo
     n1=n2+1
  enddo
  
  !Sinyuk Desert Dust 25 sizes, Mie parameter data
  !-----------------------------------------------

  do n=1,25
     read (12,3001) (srdqex(k,n),k=1,6)
     read (12,3001) (rdum5(k,n),k=1,6)
     read (12,3001) (rdum5(k,n),k=1,6)
  enddo
  read (12,3008) (rdum5(1,n),n=1,25)
  read (12,3009) (refd25(n),n=1,25)
  read (12,3010) (rdum5(1,n),n=1,25)
  do n=1,25
     read (12,3000) title
     read (12,3004) (rdum6(k,n),k=1,33)
     read (12,3005) (rdum6(k,n),k=1,33)
     read (12,3005) (rdum6(k,n),k=1,33)
     read (12,3005) (rdum6(k,n),k=1,33)
  enddo

  close(12)

  do na = 1, 4
     areff=refdry(na)
     call getmie(na, areff, srhqex(1,1,na))
     rhinfo(1,2,na) = 1.0
  enddo

  do na = 5, 6
     areff=refdry(na)
     call getmie(na, areff, srbqex(1,na))
  enddo

  do n = 2, 190
     do na = 1, 4
        do k = 1, 6
           srhqex(k,n,na) = srhqex(k,1,na)
        enddo
        rhinfo(n,2,na) = rhinfo(1,2,na)
     enddo
  enddo

  do na = 1, 4
     call set_aerosol_relhum(refdry(na),na,sruqex,sruqsc,sruqcb, truqex,  &
          truqsc, truqcb, refu22, q55u22, srhqex(1,1,na), rhinfo(1,1,na))
  enddo

  do na = 1, 8
     call getmie(7, redust(na), qxdust(1,na))
  enddo

  open(12,file=trim(dust_file),form='unformatted',status='old',  &
       convert='big_endian')
  read(12)tdust
  close(12)

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_FULIOU_RADINIT_OUT,rst_buf)
#endif
  return
end subroutine fuliou_radinit

!==================================================================
subroutine getmie(na, areff, sqex)

  use aero_coms, only: sruqex, refu22, srdqex, refd25, srsqex, refs25

  implicit none

  real, intent(in) :: areff
  integer :: n0
  integer, intent(in) :: na
  integer :: k
  integer :: n
  integer :: nn
  real :: wts
  real :: wta
  real, dimension(8), parameter :: frsulf = (/ 0.0, 0.0, 0.0, 0.33, 0.0,  &
       0.0, 0.0, 1.0 /)
  real, dimension(25) :: qxaern
  real, dimension(6), intent(out) :: sqex

  if(na < 5)then
     n0=0
     if(na == 2) n0=22
     if(na == 3) n0=44
     if(na == 4) n0=88
     do k = 1, 6
        do n = 1, 22
           nn = n0 + n
           wts = frsulf(na)
           wta = 1.0 - wts
           qxaern(n) = sruqex(k,nn) * wta + sruqex(k,n) * wts
        enddo
        call spline_giss(refu22,qxaern,22,areff,sqex(k),1.0,1.0,1)
     enddo
  endif

  if(na.eq.5.or.na.eq.6) then  !   na : aerosol compositions bic,bcb
     do k=1,6
        do n=1,25
           qxaern(n)=srsqex(k,n)
        enddo
        call spline_giss(refs25,qxaern,25,areff,sqex(k),1.0,1.0,1)
     enddo
  endif

  if(na.eq.7) then                  !   na : aerosol composition dst
     do k=1,6
        do n=1,25
           qxaern(n)=srdqex(k,n)
        enddo
        call spline_giss(refd25,qxaern,25,areff,sqex(k),1.0,1.0,1)
     enddo
  endif


  return
end subroutine getmie

!=========================================================================

subroutine spline_giss(x,f,nxf,xx,ff,cuspwm,cuspwe,kxtrap)
  implicit none

  integer, intent(in)  :: nxf, kxtrap
  real,    intent(in)  :: x(nxf), f(nxf), xx, cuspwm,cuspwe
  real,    intent(out) :: ff
  
!---------------------------------------------------------------------
!
!    spline locates xx between points (f2,x2)(f3,x3) on 4-point spread
!       and returns 4-point cubic spline interpolated value ff = f(xx)
!
!    quadratic derivatives of spline are continuous at (f2,x2),(f3,x3)
!    (x-coordinate may be specified in increasing or decreasing order)
!
!---------------------------------------------------------------------
!
!    custom control parameters:  cuspwm,cuspwe,kxtrap
!------------------------------
!
!    in cases where data points are unevenly spaced and/or data points
!    exhibit abrupt changes in value, spline interpolation may produce
!    undesirable bulging of interpolated values. in more extreme cases
!    linear interpolation may be less problematic to use.
!
!    interpolation can be weighted between: cubic spline and linear by
!    adjusting weights cuspwm and cuspwe to values between 1.0 and 0.0
!
!    cuspwm = cubic spline weight at the (x2-x3) interval mid-point
!    cuspwe = cubic spline weight at the (x2-x3) interval end-points
!
!    for example, with:
!
!    cuspwm=1.0,cuspwe=1.0  ff returns cubic spline interpolated value
!    cuspwm=0.0,cuspwe=0.0  ff returns   linearly   interpolated value
!
!---------------------------------------------------------------------
!
!     extrapolation for xx outside of defined interval:  x(1)<->x(nxf)
!
!               kxtrap = 0    no extrapolation  (i.e., sets f(xx)=0.0)
!                        1    fixed extrapolation (f(xx) = edge value)
!                        2    linear extrapolation using 2 edge points
!
!---------------------------------------------------------------------

  real :: x1,x2,x3,x4, x21,x32,x43,x31,x42, betw,ffcusp,fflinr,cuspwt
  real :: f1,f2,f3,f4, f21,f32,f43,f3221,f4332,  a,b,c,d, xf,xe,xexm
  integer k
  
  k=2
  x2=x(k)
  x3=x(nxf-1)
  betw=(xx-x2)*(x3-xx)
  if(betw.le.0.0) go to 120
  
100 continue
  k=k+1
  x3=x(k)
  betw=(xx-x2)*(x3-xx)
  if(betw.ge.0.0) go to 110
  x2=x3
  go to 100
  
110 continue
  f3=f(k)
  f4=f(k+1)
  x4=x(k+1)
  f2=f(k-1)
  x2=x(k-1)
  f1=f(k-2)
  x1=x(k-2)
  x21=x2-x1
  x31=x3-x1
  x32=x3-x2
  x43=x4-x3
  x42=x4-x2
  f21=(f2-f1)/(x21*x21)
  f32=(f3-f2)/(x32*x32)
  f43=(f4-f3)/(x43*x43)
  f3221=(f32+f21)/x31*x21
  f4332=(f43+f32)/x42*x43
  a=f2
  b=x32*f3221
  c=3.0*f32-f3221-f3221-f4332
  d=(f3221+f4332-f32-f32)/x32
  xf=xx-x2
  
  !ffcusp= cubic spline interpolation result
  !-----------------------------------------
  
  ffcusp=a+xf*(b+xf*(c+xf*d))
  xe=(x3+x2-xx-xx)/x32
  if(xe.lt.0.0) xe=-xe
  xexm=xe**2
  cuspwt=(1.0-xexm)*cuspwm+xexm*cuspwe
  
  !fflinr= linear interpolation result
  !-----------------------------------
  fflinr=a+xf*f32*x32
  ff=ffcusp*cuspwt+fflinr*(1.0-cuspwt)
  go to 160
  
  !edge point interval interpolation and/or extrapolation
  !------------------------------------------------------
120 continue
  betw=(x2-xx)*(x3-x2)
  if(betw.lt.0.0) go to 140
  
  !x(1),x(2)  edge point interval interpolation
  !--------------------------------------------
  x1=x(1)
  f1=f(1)
  f2=f(2)
  x21=x2-x1
  f21=(f2-f1)/x21
  xf=xx-x1
  betw=(x2-xx)*xf
  if(betw.lt.0.0) go to 130
  f3=f(3)
  x3=x(3)
  x32=x3-x2
  x31=x3-x1
  c=((f3-f2)/x32-f21)/x31
  b=f21-x21*c
  a=f1
  ffcusp=a+xf*(b+xf*c)
  fflinr=a+xf*f21
  xe=1.0-2.0*xf/x21
  if(xe.lt.0.0) xe=-xe
  xexm=xe**2
  cuspwt=(1.0-xexm)*cuspwm+xexm*cuspwe
  ff=ffcusp*cuspwt+fflinr*(1.0-cuspwt)
  go to 160
  
130 continue
  !extrapolation for xx outside of interval x(1) - x(2)
  !----------------------------------------------------
  !if(kxtrap.eq.0)  (no extrapolation:  sets f(xx)=0.0)
  !if(kxtrap.eq.1)  (extrapolation at fixed edge value)
  !if(kxtrap.eq.2)  (2 edge point linear extrapolation)

  if(kxtrap.eq.0) ff=0.0
  if(kxtrap.eq.1) ff=f1
  if(kxtrap.eq.2) ff=f1+xf*f21
  go to 160
  
140 continue
  !x(nxf-1),x(nxf)  edge point interval interpolation
  !--------------------------------------------------
  f3=f(nxf)
  x3=x(nxf)
  f2=f(nxf-1)
  x2=x(nxf-1)
  x32=x3-x2
  f32=(f3-f2)/x32
  xf=xx-x3
  betw=(x2-xx)*(xx-x3)
  if(betw.lt.0.0) go to 150
  f1=f(nxf-2)
  x1=x(nxf-2)
  x21=x2-x1
  x31=x3-x1
  f21=(f2-f1)/x21
  xf=xx-x2
  
  !3 -point quadratic interpolation for edge intervals
  !--------------------------------------------------
  !
  !(edge option)     ----------------------------------------------
  !for linear interpolation within edge intervals
  !between x(1),x(2), and between x(nxf-1),x(nxf)
  !set the value of coefficient c below, to c=0.0
  !----------------------------------------------
  
  c=(f32-f21)/x31
  b=f21+x21*c
  a=f2
  ffcusp=a+xf*(b+xf*c)
  fflinr=a+xf*f32
  xe=1.0-2.0*xf/x32
  if(xe.lt.0.0) xe=-xe
  xexm=xe**2
  cuspwt=(1.0-xexm)*cuspwm+xexm*cuspwe
  ff=ffcusp*cuspwt+fflinr*(1.0-cuspwt)
  go to 160
  
150 continue
  !extrapolation for x outside of interval  x(nxf-1)-x(nxf)
  !--------------------------------------------------------
  !if(kxtrap.eq.0)  (no extrapolation:  sets f(xx)=0.0)
  !if(kxtrap.eq.1)  (extrapolation at fixed edge value)
  !if(kxtrap.eq.2)  (2 edge point linear extrapolation)
  
  if(kxtrap.eq.0) ff=0.0
  if(kxtrap.eq.1) ff=f3
  if(kxtrap.eq.2) ff=f3+xf*(f3-f2)/(x3-x2)
  
160 continue
  return
end subroutine spline_giss

!=================================================================
subroutine set_aerosol_relhum(reff0, naer, sruqex, sruqsc, sruqcb,   &
     truqex, truqsc, truqcb, refu22, q55u22, srhqex, rhdata)
  implicit none

  integer, intent(in) :: NAER
  integer, parameter :: KDREAD=13
  real, intent(in) :: REFF0
  real, dimension(6,110), intent(in) :: SRUQEX,SRUQSC,SRUQCB
  real, dimension(33,110), intent(in) :: TRUQEX,TRUQSC,TRUQCB
  real, dimension(110), intent(in) :: REFU22,Q55U22
  real SRHQSC(6,190),SRHQCB( 6,190),TRHQAB(33,190)
  real, dimension(6,190), intent(out) :: srhqex
  real, dimension(190,15), intent(out) :: rhdata

!     ------------------------------------------------------------------
!     REFF0  = Effective radius for dry aerosol seed size (in microns)
!     NAER   = Aerosol composition index
!     KDREAD = IO READ unit number for Q(m,r),g(m,r) data used by SETREL
!     ------------------------------------------------------------------
!     Aerosol index = NAER    Composition         Input data order = NNA
!                      1      SO4  Sulfate                            1
!                      2      SEA  Sea Salt                           2
!                      3      NO3  Nitrate                            3
!                                  Pure Water                         4
!                      4      ORG  Organic                            5
!     ------------------------------------------------------------------

  character*128, save :: dtfile='/cluster/home/dmm31/aerosols/oct2003.relhum.nr.Q633G633.table'
  logical, parameter :: qbinary=.false.  ; logical qexist
  integer i,j,j1,k,k1,in1,ir1,jdry,jwet,jhimax,khimax,maxdry,maxwet
  integer n,n0,n1,nn,np,nrhn1
  real x,xx,xi,xn0,xn1,xr1,ff,fi,gi,gd1,gd2,gw1,gw2,grh,qrh,rrh
  real rh,rhi,rr0,rd1,rd2,rw1,rw2,dwr,qd1,qd2,qw1,qw2,xdry,sdry
  real xwet,swet,qqdmax,qqwmax,rqdmax,rqwmax,q55dry,q63dry,dwn
  real aermas,ddry,dwet,reffi,rhrhi,sum,sumw,vd1,vd2,vw1,vw2
  real w1,w2,w3,w4,wd1,wd2,ww1,ww2,wtx,wty,wtz,wts,wta,xfdry
  real q55rh1,q55rh2,q55rh3,q55rh4,q550,q633,qgaerx,qscqcb

!     Output variables (RHDATA/RHINFO)

  real, dimension(8), parameter :: FRSULF=(/0.000, 0.000, 0.000, 0.330, 0.000, 0.000, 0.000, 1.00/)
  REAL    RHRHRH(190),RHTAUF(190),RHREFF(190) &
       ,RHWGM2(190),RHDGM2(190),RHTGM2(190) &
       ,RHXMFX(190),RHDENS(190),RHQ550(190) &
       ,TAUM2G(190),XNRRHX(190),ANBCM2(190) &
       ,COSBAR(190),PIZERO(190),ANGSTR(190),RHINFO(190,15)
  
  EQUIVALENCE (RHINFO(1, 1),RHRHRH(1)),(RHINFO(1, 2),RHTAUF(1))
  EQUIVALENCE (RHINFO(1, 3),RHREFF(1)),(RHINFO(1, 4),RHWGM2(1))
  EQUIVALENCE (RHINFO(1, 5),RHDGM2(1)),(RHINFO(1, 6),RHTGM2(1))
  EQUIVALENCE (RHINFO(1, 7),RHXMFX(1)),(RHINFO(1, 8),RHDENS(1))
  EQUIVALENCE (RHINFO(1, 9),RHQ550(1)),(RHINFO(1,10),TAUM2G(1))
  EQUIVALENCE (RHINFO(1,11),XNRRHX(1)),(RHINFO(1,12),ANBCM2(1))
  EQUIVALENCE (RHINFO(1,13),COSBAR(1)),(RHINFO(1,14),PIZERO(1))
  EQUIVALENCE (RHINFO(1,15),ANGSTR(1))
  
!     ------------------------------------------------------------------
!     RHDATA/    Local
!     RHINFO   Variable                   Description
!     ------   --------   ----------------------------------------------
!        1      RHRHRH    Relative humidity index RH (DO 110  0.0-0.999)
!        2      RHTAUF    Dry TAU multiplication factor due to RH effect
!        3      RHREFF    RH dependent effective radius
!        4      RHWGM2    Liquid water content (g/m2) per unit (dry) TAU
!        5      RHDGM2    Dry mass density     (g/m2) per unit (dry) TAU
!        6      RHTGM2    Total mass density   (g/m2) per unit (dry) TAU
!        7      RHXMFX    Dry mass fraction X of total aerosol mass
!        8      RHDENS    RH dependent density (g/cm3)
!        9      RHQ550    RH dependent Mie extinction efficiency (550nm)
!       10      TAUM2G    RH dependent TAU factor  (m2/g) of dry aerosol
!       11      XNRRHX    RH dependent real refractive index
!       12      ANBCM2    Aerosol Number density (Billion)/cm2
!       13      COSBAR    RH dependent Mie asymmetry parameter (visible)
!       14      PIZERO    RH dependent single scattering albedo(visible)
!       15      ANGSTR    Angstrom exponent = -(1-SRHQEX(5)/(0.55-0.815)
!     ------------------------------------------------------------------


!     Local variables

  REAL    R633NR(890),XNR(31),Q633NR(890,31),G633NR(890,31)
  REAL    Q880M1(890),G880M1(890),Q880M0(890),G880M0(890)
  REAL    Q880N1(890),Q880N0(890),R550NR(890),SMOOTH(890)
  REAL    RR0RHX(190),QRH633(190),GRH633(190),DNRX(190)


  REAL    QXAERN(33),QSAERN(33),QGAERN(33) &
       ,SR1QEX( 6),SR1QSC( 6),SR1QCB( 6) &
       ,SR2QEX( 6),SR2QSC( 6),SR2QCB( 6) &
       ,SR3QEX( 6),SR3QSC( 6),SR3QCB( 6) &
       ,SR4QEX( 6),SR4QSC( 6),SR4QCB( 6) &
       ,TR1QEX(33),TR1QSC(33),TR1QCB(33) &
       ,TR2QEX(33),TR2QSC(33),TR2QCB(33) &
       ,TR3QEX(33),TR3QSC(33),TR3QCB(33) &
       ,TR4QEX(33),TR4QSC(33),TR4QCB(33) &
       ,TRHQEX(33),TRHQSC(33),TRHQCB(33)
  
  integer, parameter, dimension(4) :: NRHCRY=(/38,47,28,38/)
  
  CHARACTER*8 AERTYP(4)
  DATA AERTYP/'Sulfate ','SeaSalt ','Nitrate ','Organic '/

!     ------------------------------------------------------------------
!     Hygroscopic aerosols (Sulfate,SeaSalt,Nitrate) physical properties
!     formulas from Tang and Munkelwitz (1994, 1996) in JGR 99, JGR 101.
!
!     AW=water activity RO=density  BX=growth factor RX=refractive index
!     SO4 = ammonium sulfate;   SEA = sea salt;   NO3 = ammonium nitrate
!     ------------------------------------------------------------------

!     functions

  REAL, external :: AWSO4,DWSO4,ROSO4,BXSO4,RXSO4,DRWSO4,DRDSO4
  REAL, external :: AWSEA,DWSEA,ROSEA,BXSEA,RXSEA,DRWSEA,DRDSEA
  REAL, external :: RRSEA,VVSEA,GXSEA
  REAL, external :: AWNO3,DWNO3,RONO3,BXNO3,R1NO3,R2NO3, DRXNO3
  REAL, external :: AWOCX,DWOCX,ROOCX,BXOCX,RXOCX,DRWOCX,DRDOCX


!     ------------------------------------------------------------------
!     Q,G Mie data (879x31) at 0.633 microns, use 31 points to cover the
!     refractive index from 1.30 to 1.60 with equal spacing of 0.01
!
!     Q,G data effective radius spans the range from 0.0 to 20.4 microns
!     in (3) segments of equally spaced data for optimized 4-point Cubic
!     Spline interpolation.  The equally spaced segments are as follows:
!
!     Index:    1 - 303   304 - 603   604 -  879   881 - 885   886 - 890
!     Reff:   0.00-3.02   3.04-9.02   9.04-20.04   2.98-3.04   8.96-9.08
!     Delta:     0.01        0.02        0.04         0.02        0.04
!
!     The last two intervals are constructed to accommodate transitions
!     between the (3) segments using 4-point Cubic Spline interpolation
!     ------------------------------------------------------------------


  inquire (file=trim(dtfile),exist=qexist)
  if(.not.qexist) dtfile='RH_QG_Mie  ' ! generic name used by GCM
  inquire (file=trim(dtfile),exist=qexist)
  if(.not.qexist)then
     print*,'setrel: no RH_QG files',255
     stop
  endif
  open(kdread, file=trim(dtfile), form='formatted', status='old')

  READ (KDREAD,7000) (XNR(J),J=1,31)
7000 FORMAT(12X,F5.3,30F8.3)
  DO I=1,880
     READ (KDREAD,7001) R633NR(I),(Q633NR(I,J),J=1,31)
7001 FORMAT(3X,F6.2,31F8.5)
  enddo
  READ (KDREAD,7000) (XNR(J),J=1,31)
  DO I=1,880
     READ (KDREAD,7001) R633NR(I),(G633NR(I,J),J=1,31)
  enddo
  close(kdread)
  
  J=880
  firstk: DO K=299,305
     IF(K.EQ.300)cycle firstk
     IF(K.EQ.302)cycle firstk
     J=J+1
     R633NR(J)=R633NR(K)
     DO I=1,31
        Q633NR(J,I)=Q633NR(K,I)
        G633NR(J,I)=G633NR(K,I)
     enddo
  enddo firstk
  secondk: DO K=600,606
     IF(K.EQ.601)cycle secondk
     IF(K.EQ.603)cycle secondk
     J=J+1
     R633NR(J)=R633NR(K)
     DO I=1,31
        Q633NR(J,I)=Q633NR(K,I)
        G633NR(J,I)=G633NR(K,I)
     enddo
  enddo secondk
        
!     Apply 13-point quadratic least-squares smoothing to large particle
!     portion of Mie Qx data to eliminate low-amplitude ripple in Q633NR
!     (Monotonic size dependence is needed for inverse Qx interpolation)
!     (Smoothing affects 4th decimal of Q633NR for large particle sizes)
!     ------------------------------------------------------------------
  DO I=1,31
     DO J=1,880
        SMOOTH(J)=Q633NR(J,I)
     END DO
     DO J=881,886
        SMOOTH(J)=SMOOTH(880)
     END DO
     firstj: DO J=250,880
        J1=J-2
        IF(SMOOTH(J).GE.SMOOTH(J-1))exit firstj
     END DO firstj
     DO J=J1,880
        SUM=4550.0/13.0*SMOOTH(J)
        DO K=1,6
           SUM=SUM+(4550.0/13.0-14*K*K)*(SMOOTH(J-K)+SMOOTH(J+K))
        END DO
        Q633NR(J,I)=SUM/2002.0
     enddo
  enddo

!                                     Set relative humidity RHRHRH scale
!                                     ----------------------------------
  DO I=1,190
     RHRHRH(I)=(I-1)/100.0
     IF(I.GT.91) RHRHRH(I)=0.90+(I-91)/1000.0
  enddo
  
!         Define RH (=AW), RO, BX, RX as functions of X for NAER aerosol
!         --------------------------------------------------------------
  NRHN1=NRHCRY(NAER)+1
  DO I=1,190
     RHI=RHRHRH(I)
     RR0RHX(I)=1.0
     RHXMFX(I)=1.0
     IF(NAER.EQ.1) THEN       !    Dry Sulfate refrac index and density
        XNRRHX(I)=1.526
        RHDENS(I)=1.760
        IF(I.LT.NRHN1) DNRX(I)=DRDSO4(RHI)
        IF(I.GE.NRHN1) DNRX(I)=DRWSO4(RHI)
     ENDIF
     IF(NAER.EQ.2) THEN       !    Dry SeaSalt refrac index and density
        XNRRHX(I)=1.490
        RHDENS(I)=2.165
        IF(I.LT.NRHN1) DNRX(I)=DRDSEA(RHI)
        IF(I.GE.NRHN1) DNRX(I)=DRWSEA(RHI)
     ENDIF
     IF(NAER.EQ.3) THEN       !    Dry Nitrate refrac index and density
        XNRRHX(I)=1.554
        RHDENS(I)=1.725
        DNRX(I)=DRXNO3(RHRHRH(I))
     ENDIF
     IF(NAER.EQ.4) THEN       !    Dry Organic refrac index and density
        XNRRHX(I)=1.526          !                  (representative value)
        RHDENS(I)=1.5            !                  (representative value)
        IF(I.LT.NRHN1) DNRX(I)=DRDOCX(RHI)
        IF(I.GE.NRHN1) DNRX(I)=DRWOCX(RHI)
     ENDIF
  enddo
     
!            Invert X, RO, BX, RX functions of (X) to be functions of RH
!            -----------------------------------------------------------
  I=191
  FF=1.0
  XX=0.0
  IF(NAER.EQ.1) GI=DWSO4(XX)
  IF(NAER.EQ.2) GI=DWSEA(XX)
  IF(NAER.EQ.3) GI=DWNO3(XX)
  IF(NAER.EQ.4) THEN
     FF=.9995
     XX=(1.0-FF)**.125
     GI=DWOCX(XX)
  END IF
112 I=I-1
  FI=RHRHRH(I)
  DO K=1,5
     XI=XX-(FF-FI)/GI
     IF(NAER.EQ.1) FF=AWSO4(XI)
     IF(NAER.EQ.2) FF=AWSEA(XI)
     IF(NAER.EQ.3) FF=AWNO3(XI)
     IF(NAER.EQ.4) FF=AWOCX(XI)
     IF(I.GT.0) THEN
     ENDIF
     XX=XI
     IF(NAER.EQ.1) GI=DWSO4(XX)
     IF(NAER.EQ.2) GI=DWSEA(XX)
     IF(NAER.EQ.3) GI=DWNO3(XX)
     IF(NAER.EQ.4) GI=DWOCX(XX)
  enddo
  RHXMFX(I)=XX
  IF(NAER.EQ.1) THEN          !       RH dependent Sulfate X,R,NR,RO
     RHDENS(I)=ROSO4(XX)
     RR0RHX(I)=BXSO4(XX)
     XNRRHX(I)=RXSO4(XX)
  ENDIF
  IF(NAER.EQ.2) THEN          !       RH dependent SeaSalt X,R,NR,RO
     RHDENS(I)=ROSEA(XX)
     RR0RHX(I)=BXSEA(XX)
     XNRRHX(I)=RXSEA(XX)
  ENDIF
  IF(NAER.EQ.3) THEN          !       RH dependent Nitrate X,R,NR,RO
     RHDENS(I)=RONO3(XX)
     RR0RHX(I)=BXNO3(XX)
     XNRRHX(I)=R1NO3(XX)
     IF(XX.GT.0.205) XNRRHX(I)=R2NO3(XX)
  ENDIF
  IF(NAER.EQ.4) THEN          !       RH dependent Organic X,R,NR,RO
     RHDENS(I)=ROOCX(XX)
     RR0RHX(I)=BXOCX(XX)
     XNRRHX(I)=RXOCX(XX)
  ENDIF
  IF(I.GT.NRHN1) GO TO 112
  
  !     ------------------------------------------------------------------
  !     Find Qdry(r),gdry(r) from Q(m,r),g(m,r) maps for each aerosol type
  !     Find Qwet(r),gwet(r) from Q(m,r),g(m,r) maps for each aerosol type
  !     also locate MAXDRY,MAXWET pts where Qdry(r),Qwet(r) are at maximum
  !          (M1 refers to mass fraction X of 1.0, i.e., "dry" aerosol)
  !          (M0 refers to mass fraction X of 0.0, i.e., "wet" aerosol)
  !     ------------------------------------------------------------------
  MAXDRY=1
  MAXWET=1
  QQDMAX=0.0
  QQWMAX=0.0
  XDRY=XNRRHX(1)
  !     IF(MCRYON.EQ.1) XDRY=XNRRHX(NRHN1) ! If "dry" = RHC reference line
  SDRY=XDRY*100.0-129
  JDRY=SDRY
  DDRY=SDRY-JDRY
  XWET=1.3330                     !  Pure water Nr = "wet" aerosol
  SWET=XWET*100.0-129
  JWET=SWET
  DWET=SWET-JWET
  DO I=1,880
     CALL SPLNI4(Q633NR,890,31,I,JDRY,DDRY,Q880M1(I))
     CALL SPLNI4(G633NR,890,31,I,JDRY,DDRY,G880M1(I))
     CALL SPLNI4(Q633NR,890,31,I,JWET,DWET,Q880M0(I))
     CALL SPLNI4(G633NR,890,31,I,JWET,DWET,G880M0(I))
     IF(Q880M1(I).GT.QQDMAX) THEN
        QQDMAX=Q880M1(I)
        MAXDRY=I
     ENDIF
     IF(Q880M0(I).GT.QQWMAX) THEN
        QQWMAX=Q880M0(I)
        MAXWET=I
     ENDIF
  enddo
  RQDMAX=R633NR(MAXDRY)
  RQWMAX=R633NR(MAXWET)

!     Define:  Qdry(r) and Qwet(r) at the reference wavelength of 550 nm
!              using refractive index off-set and size parameter scaling
!     ------------------------------------------------------------------
  XDRY=XNRRHX(1)*DNRX(1)             !      Dry aerosol Nr at 550 nm
!     IF(MCRYON.EQ.1) XDRY=XNRRHX(NRHN1) ! If "dry" = RHC reference line
  SDRY=XDRY*100.0-129
  JDRY=SDRY
  DDRY=SDRY-JDRY
  XWET=1.3330*1.001179     !       Pure water aerosol Nr at 550 nm
  SWET=XWET*100.0-129
  JWET=SWET
  DWET=SWET-JWET
  DO I=1,880
     CALL SPLNI4(Q633NR,890,31,I,JDRY,DDRY,Q880N1(I))
     CALL SPLNI4(Q633NR,890,31,I,JWET,DWET,Q880N0(I))
     R550NR(I)=R633NR(I)*(0.550/0.633) !  Size shift refers Q to 550 nm
  enddo
  CALL SPLINE_giss(R550NR,Q880N1,880,REFF0,Q55DRY,1.0,1.0,1)
  CALL SPLINE_giss(R633NR,Q880M1,880,REFF0,Q63DRY,1.0,1.0,1)
  
  !     Find Q(RH),g(RH) paths in Q(m,r),g(m,r) maps for seed size = REFF0
  !     2-coordinate paths defined via XN0=XNRRHX(I) & RR0=REFF0*RR0RHX(I)
  !     ------------------------------------------------------------------
  DO I=1,190
     XN0=XNRRHX(I)
     XN1=XN0*100.0-129
     IN1=XN1
     DWN=XN1-IN1
     RR0=REFF0*RR0RHX(I)
     IF(RR0.LT.0.01) RR0=0.01
     IF(RR0.LE.3.00) XR1=RR0*100.0+1
     IF(RR0.GT.3.00.AND.RR0.LT.3.04) XR1=RR0*50.0+732
     IF(RR0.GE.3.04.AND.RR0.LE.9.00) XR1=RR0*50.0+152
     IF(RR0.GT.9.00.AND.RR0.LT.9.08) XR1=RR0*25.0+662
     IF(RR0.GE.9.08) THEN
        XR1=RR0*25.0+378
        IF(XR1.GT.877.9999) XR1=877.9999
     ENDIF
     IR1=XR1
     DWR=XR1-IR1
     CALL SPLN44(Q633NR,890,31,IR1,DWR,IN1,DWN,QRH633(I))
     CALL SPLN44(G633NR,890,31,IR1,DWR,IN1,DWN,GRH633(I))
  enddo
     
  !     Define Q55(RH) by tracing path in Q(m,r) map for RH dependent size
  !     via 2-coordinate path XN0=XNRRHX(I)*DNRX(I), RR0=RRH(I)*(.633/.55)
  !     ------------------------------------------------------------------
  DO I=1,190
     XN0=XNRRHX(I)*DNRX(I)
     XN1=XN0*100.0-129
     IN1=XN1
     DWN=XN1-IN1
     RR0=REFF0*RR0RHX(I)*(0.633/0.550)
     IF(RR0.LT.0.01) RR0=0.01
     IF(RR0.LE.3.00) XR1=RR0*100.0+1
     IF(RR0.GT.3.00.AND.RR0.LT.3.04) XR1=RR0*50.0+732
     IF(RR0.GE.3.04.AND.RR0.LE.9.00) XR1=RR0*50.0+152
     IF(RR0.GT.9.00.AND.RR0.LT.9.08) XR1=RR0*25.0+662
     IF(RR0.GE.9.08) THEN
        XR1=RR0*25.0+378
        IF(XR1.GT.877.9999) XR1=877.9999
     ENDIF
     IR1=XR1
     DWR=XR1-IR1
     CALL SPLN44(Q633NR,890,31,IR1,DWR,IN1,DWN,RHQ550(I))
     RHREFF(I)=RR0RHX(I)*REFF0
  enddo
  
  !        Aerosol liquid water content is in kg/m2 per unit optical depth
  !     of dry aerosol with aerosol effective radius expressed in microns.
  !     ------------------------------------------------------------------
  
  DO I=1,190
     RHTAUF(I)=(RHQ550(I)/Q55DRY)*RR0RHX(I)**2
     AERMAS=1.33333333*RHREFF(I)*RHDENS(I)/RHQ550(I)*RHTAUF(I)
     RHTGM2(I)=AERMAS
     RHDGM2(I)=AERMAS*RHXMFX(I)
     RHWGM2(I)=RHTGM2(I)-RHDGM2(I)
     TAUM2G(I)=0.75/RHDENS(1)/RHREFF(1)*RHQ550(1)*RHTAUF(I)
     ANBCM2(I)=TAUM2G(I)/(1.5080*RHQ550(I)*RHREFF(I)**2)
  enddo
  
  !     Determination of RH dependent Mie scattering tables for GCM input.
  !     Find equivalent aersol dry sizes (RD1,RD2) and wet sizes (RW1,RW2)
  !     and corresponding weights to match the RH dependent Q(r) and g(r).
  !     Fits made to form: QRH=X*[Y*QD1+(1-Y)*QD2]+(1-X)*[Z*WD1+(1-Z)*WD2]
  !     ------------------------------------------------------------------
  J1=MAXWET
  JHIMAX=881-MAXWET
  K1=MAXDRY
  KHIMAX=881-MAXDRY
  NP=190-NRHN1+1
  DO I=1,190
     RHRHI=RHRHRH(I)
     XFDRY=RHXMFX(I)
     REFFI=RHREFF(I)
     RRH=RR0RHX(I)*REFF0
     GRH=GRH633(I)
     QRH=QRH633(I)
     QD1=QRH
     QD2=QRH
     QW1=QRH
     QW2=QRH
     IF(QW1.GT.QQWMAX) QW1=QQWMAX
     IF(QW2.GT.QQWMAX) QW2=QQWMAX
     CALL SPLINV(R633NR    ,Q880M0    ,MAXWET,RW1,QW1,1.0,1.0,1)
     CALL SPLINV(R633NR(J1),Q880M0(J1),JHIMAX,RW2,QW2,1.0,1.0,1)
     CALL SPLINE_giss(R633NR,G880M0,880,RW1,GW1,1.0,1.0,1)
     CALL SPLINE_giss(R633NR,G880M0,880,RW2,GW2,1.0,1.0,1)
     IF(I.GE.NRHN1.AND.QRH.GT.QQWMAX) THEN
        QD1=QQWMAX+(QRH-QQWMAX)/XFDRY ! QD1 such that  QRH=X*QD1+(1-X)*QW1
        QD2=2.3                     ! 2 dry sizes are used if QD1>QQWMAX
     ENDIF
     CALL SPLINV(R633NR    ,Q880M1    ,MAXDRY,RD1,QD1,1.0,1.0,1)
     CALL SPLINV(R633NR(K1),Q880M1(K1),KHIMAX,RD2,QD2,1.0,1.0,1)
     CALL SPLINE_giss(R633NR,G880M1,880,RD1,GD1,1.0,1.0,1)
     CALL SPLINE_giss(R633NR,G880M1,880,RD2,GD2,1.0,1.0,1)
     
     IF(I.LT.NRHN1) THEN         !          Pure dry aerosol region (1)
        WTX=1.0
        WTY=1.0
        IF(REFF0.GT.RQDMAX) WTY=0.0
        WTZ=1.0
     ELSE            !         Dry/wet weighted average regions (2)-(4)
        IF(QRH.LE.QQWMAX.AND.REFFI.LT.RW1) THEN  !   Small-size region (2)
           WTZ=1.0
           WTY=1.0
           WTX=(GRH-GW1)/(GD1-GW1)
        ENDIF
        IF(QRH.GT.QQWMAX) THEN                   !  Medium-size region (3)
           !     Fit form: QRH=X*(Y*QD1+(1-Y)*QD2)+(1-X)*QWmax & QRH=/QD1=/QD2=/QW1
           WTZ=1.0
           WTY=((GRH-GD2)*QRH*QD2+(GD2-GW1)*QD2*QW1+(GW1-GRH)*QRH*QW1) &
                /((GD1-GRH)*QRH*QD1+(GRH-GD2)*QRH*QD2+(GW1-GD1)*QD1*QW1 &
                +(GD2-GW1)*QD2*QW1)
           WTX=(QRH-QW1)/(WTY*(QD1-QD2)+(QD2-QW1))
        ENDIF
        IF(QRH.LE.QQWMAX.AND.REFFI.GT.RW1) THEN   !  Large size region (4)
           WTY=0.0
           WTZ=0.0
           WTX=(GRH-GW2)/(GD2-GW2)
        ENDIF
     ENDIF
     IF(REFFI.GT.RQWMAX.AND.RHRHI.GT.0.995) THEN   ! High RH region (5)
        WTY=0.0
        WTX=XFDRY
        WTZ=((GRH-GW2)-(GD2-GW2)*WTX)/((1.0-WTX)*(GW1-GW2))
     ENDIF
     
     VD1=WTX*WTY
     VD2=WTX*(1.0-WTY)
     VW1=WTZ*(1.0-WTX)
     VW2=(1.0-WTZ)*(1.0-WTX)
     RD1=MIN(RD1,10.0)
     RD2=MIN(RD2,10.0)
     RW1=MIN(RW1,10.0)
     RW2=MIN(RW2,10.0)
     
     !     Computed weight factors are for Lab reference wavelength of 633nm.
     !     Rescale spectral extinction to 550 nm & renormalize weight factors
     !     ------------------------------------------------------------------
     CALL SPLINE_GISS(R550NR,Q880N1,880,RD1,Q550,1.0,1.0,1)
     CALL SPLINE_GISS(R633NR,Q880M1,880,RD1,Q633,1.0,1.0,1)
     WD1=VD1*(Q550/Q633)
     CALL SPLINE_GISS(R550NR,Q880N1,880,RD2,Q550,1.0,1.0,1)
     CALL SPLINE_GISS(R633NR,Q880M1,880,RD2,Q633,1.0,1.0,1)
     WD2=VD2*(Q550/Q633)
     CALL SPLINE_GISS(R550NR,Q880N0,880,RW1,Q550,1.0,1.0,1)
     CALL SPLINE_GISS(R633NR,Q880M0,880,RW1,Q633,1.0,1.0,1)
     WW1=VW1*(Q550/Q633)
     CALL SPLINE_GISS(R550NR,Q880N0,880,RW2,Q550,1.0,1.0,1)
     CALL SPLINE_GISS(R633NR,Q880M0,880,RW2,Q633,1.0,1.0,1)
     WW2=VW2*(Q550/Q633)
     SUMW=WD1+WD2+WW1+WW2
     W1=WD1/SUMW
     W2=WD2/SUMW
     W3=WW1/SUMW
     W4=WW2/SUMW
     
     !     ------------------------------------------------------------------
     !     Tabulate relative humidity dependent solar, thermal Mie scattering
     !     parameters SRHQEX,SRHQSC,SRHQCS, TRHQAB for each aerosol type NAER
     !     These are mass weighted averages of equivalent dry and wet aerosol
     !     parameters for sizes matching the relative humidity dependent Q(r)
     !     ------------------------------------------------------------------
     
     N0=0                      !      Select Mie parameters for Sulfate
     IF(NAER.EQ.2) N0=22       !      Select Mie parameters for SeaSalt
     IF(NAER.EQ.3) N0=44       !      Select Mie parameters for Nitrate
     IF(NAER.EQ.4) N0=88       !      Select Mie parameters for Organic
     N1=N0+1
     DO K=1,6                            !   SW dry sizes RD1 & RD2
        DO N=1,22
           NN=N0+N
           WTS=FRSULF(NAER)
           WTA=1.0-WTS
           QXAERN(N)=SRUQEX(K,NN)*WTA+SRUQEX(K,N)*WTS
           QSAERN(N)=SRUQSC(K,NN)*WTA+SRUQSC(K,N)*WTS
           QGAERX=SRUQCB(K,NN)*SRUQSC(K,NN)*WTA+SRUQCB(K,N)*SRUQSC(K,N)*WTS
           QGAERN(N)=QGAERX/QSAERN(N)
        enddo
        CALL SPLINE_GISS(REFU22,QXAERN,22,RD1,SR1QEX(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QSAERN,22,RD1,SR1QSC(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QGAERN,22,RD1,SR1QCB(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QXAERN,22,RD2,SR2QEX(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QSAERN,22,RD2,SR2QSC(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QGAERN,22,RD2,SR2QCB(K),1.0,1.0,1)
     enddo
     
     DO K=1,33                           !   LW dry sizes RD1 & RD2
        DO N=1,22
           NN=N0+N
           WTS=FRSULF(NAER)
           WTA=1.0-WTS
           QXAERN(N)=TRUQEX(K,NN)*WTA+TRUQEX(K,N)*WTS
           QSAERN(N)=TRUQSC(K,NN)*WTA+TRUQSC(K,N)*WTS
           QGAERX=TRUQCB(K,NN)*TRUQSC(K,NN)*WTA+TRUQCB(K,N)*TRUQSC(K,N)*WTS
           QGAERN(N)=QGAERX/(QSAERN(N)+1d-10)
        enddo
        CALL SPLINE_GISS(REFU22,QXAERN,22,RD1,TR1QEX(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QSAERN,22,RD1,TR1QSC(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QGAERN,22,RD1,TR1QCB(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QXAERN,22,RD2,TR2QEX(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QSAERN,22,RD2,TR2QSC(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QGAERN,22,RD2,TR2QCB(K),1.0,1.0,1)
     enddo
     CALL SPLINE_GISS(REFU22,Q55U22(N1),22,RD1,Q55RH1,1.0,1.0,1)
     CALL SPLINE_GISS(REFU22,Q55U22(N1),22,RD2,Q55RH2,1.0,1.0,1)
     
     N0=66                    !   Select Mie parameters for pure water
     N1=N0+1
     DO K=1,6                           !   SW wet sizes RW1 & RW2
        DO N=1,22
           NN=N0+N
           QXAERN(N)=SRUQEX(K,NN)
           QSAERN(N)=SRUQSC(K,NN)
           QGAERN(N)=SRUQCB(K,NN)
        enddo
        CALL SPLINE_GISS(REFU22,QXAERN,22,RW1,SR3QEX(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QSAERN,22,RW1,SR3QSC(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QGAERN,22,RW1,SR3QCB(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QXAERN,22,RW2,SR4QEX(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QSAERN,22,RW2,SR4QSC(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QGAERN,22,RW2,SR4QCB(K),1.0,1.0,1)
     enddo
     
     DO K=1,33                          !   LW wet sizes RW1 & RW2
        DO N=1,22
           NN=N0+N
           QXAERN(N)=TRUQEX(K,NN)
           QSAERN(N)=TRUQSC(K,NN)
           QGAERN(N)=TRUQCB(K,NN)
        enddo
        CALL SPLINE_GISS(REFU22,QXAERN,22,RW1,TR3QEX(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QSAERN,22,RW1,TR3QSC(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QGAERN,22,RW1,TR3QCB(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QXAERN,22,RW2,TR4QEX(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QSAERN,22,RW2,TR4QSC(K),1.0,1.0,1)
        CALL SPLINE_GISS(REFU22,QGAERN,22,RW2,TR4QCB(K),1.0,1.0,1)
     enddo
     CALL SPLINE_GISS(REFU22,Q55U22(N1),22,RW1,Q55RH3,1.0,1.0,1)
     CALL SPLINE_GISS(REFU22,Q55U22(N1),22,RW2,Q55RH4,1.0,1.0,1)
     
     !      Weighted GCM SW Mie scattering parameters
     DO K=1,6
        SRHQEX(K,I)=W1*SR1QEX(K)+W2*SR2QEX(K)+W3*SR3QEX(K)+W4*SR4QEX(K)
        SRHQSC(K,I)=W1*SR1QSC(K)+W2*SR2QSC(K)+W3*SR3QSC(K)+W4*SR4QSC(K)
        QSCQCB     =W1*SR1QCB(K)*SR1QSC(K)+W2*SR2QCB(K)*SR2QSC(K)  &
             +W3*SR3QCB(K)*SR3QSC(K)+W4*SR4QCB(K)*SR4QSC(K)
        SRHQCB(K,I)=QSCQCB/SRHQSC(K,I)
     enddo
     !      Weighted GCM LW Mie scattering parameters
     DO K=1,33
        TRHQEX(K)=W1*TR1QEX(K)+W2*TR2QEX(K)+W3*TR3QEX(K)+W4*TR4QEX(K)
        TRHQSC(K)=W1*TR1QSC(K)+W2*TR2QSC(K)+W3*TR3QSC(K)+W4*TR4QSC(K)
        QSCQCB   =W1*TR1QCB(K)*TR1QSC(K)+W2*TR2QCB(K)*TR2QSC(K)  &
             +W3*TR3QCB(K)*TR3QSC(K)+W4*TR4QCB(K)*TR4QSC(K)
        TRHQCB(K)=QSCQCB/TRHQSC(K)
        TRHQAB(K,I)=TRHQEX(K)-TRHQSC(K)
     enddo
     
     COSBAR(I)=SRHQCB(6,I)
     PIZERO(I)=SRHQSC(6,I)/SRHQEX(6,I)
     ANGSTR(I)=-(1.0-SRHQEX(5,I))/(0.550-0.815)
     
     !              Transfer EQUIVALENCEd SETREL output information to RHDATA
     DO J=1,15
        RHDATA(I,J)=RHINFO(I,J)
     enddo
  enddo

  return
end subroutine set_aerosol_relhum

!======================================================================

real function awso4(x)
  implicit none
  real :: x

  !        Sulfate parametric formulas from Tang & Munkelwitz(94,96)
  AWSO4 = 1.0-0.2715*X+0.3113*X**2-2.336*X**3+1.412*X**4    ! TM94

  return 
end function awso4

real function dwso4(x)
  implicit none
  real :: x

  !        Sulfate parametric formulas from Tang & Munkelwitz(94,96)
  DWSO4=    -0.2715+0.6226*X   -7.008*X**2+5.648*X**3

  return 
end function dwso4

real function roso4(x)
  implicit none
  real :: x

  !        Sulfate parametric formulas from Tang & Munkelwitz(94,96)
  ROSO4=0.9971+0.592*X-0.05036*X**2+0.01024*X**3  ! TM94
  
  return 
end function roso4

real function bxso4(x)
  implicit none
  real :: x
  real, external :: roso4
  !        Sulfate parametric formulas from Tang & Munkelwitz(94,96)
  BXSO4=(1.0/X*1.760/ROSO4(X))**(1.0/3.0)             ! TM96

  return 
end function bxso4

real function rxso4(x)
  implicit none
  real :: x

  !        Sulfate parametric formulas from Tang & Munkelwitz(94,96)
  RXSO4=1.3330+0.16730*X-0.0395*X**2                       ! TM91

  return 
end function rxso4

real function drwso4(rh)
  implicit none
  real :: rh

  !        Sulfate parametric formulas from Tang & Munkelwitz(94,96)
  DRWSO4=1.002146-0.00149*RH+0.001*RH/(1.0+0.911*RH**10)

  return 
end function drwso4

real function drdso4(rh)
  implicit none
  real :: rh

  !        Sulfate parametric formulas from Tang & Munkelwitz(94,96)
  DRDSO4=1.002503     ! ratio of wet & dry nr(0.550) / nr(0.633)

  return 
end function drdso4

real function awsea(x)
  implicit none
  real :: x
  
  !        SeaSalt parametric formulas from Tang & Munkelwitz(94,96)
  AWSEA=1.0-0.6366*X+0.8624*X**2-11.58*X**3+15.18*X**4   ! TM96

  return
end function awsea

real function dwsea(x)
  implicit none
  real :: x
  DWSEA=     -0.6366+1.7248*X   -34.74*X**2+60.72*X**3
  
  

  return
end function dwsea

real function rosea(x)
  implicit none
  real :: x
  
  ROSEA=0.9971+0.741*X-0.3741*X**2+2.252*X**3-2.060*X**4   ! TM96

  return
end function rosea

real function bxsea(x)
  implicit none
  real :: x
  real, external :: rosea
  
  BXSEA=(1.0/X*2.165/ROSEA(X))**(1.0/3.0)


  return
end function bxsea

real function rrsea(x)
  implicit none
  real :: x
  
  RRSEA=3.70958+(8.95-3.70958)/(1.0+(1.0-X)/X*58.448/18.0)


  return
end function rrsea

real function vvsea(x)
  implicit none
  real :: x
  real, external :: rosea

  VVSEA=(18.0+(58.448-18.0)/(1.0+(1.0-X)/X*58.448/18.0))/ROSEA(X)


  return
end function vvsea

real function gxsea(x)
  implicit none
  real :: x
  real, external :: rrsea,vvsea

  GXSEA=SQRT((2.0*RRSEA(X)+VVSEA(X))/(VVSEA(X)-RRSEA(X))) ! TM96


  return
end function gxsea

real function rxsea(x)
  implicit none
  real :: x
  real, external :: gxsea

  RXSEA=1.333+(GXSEA(X)-1.333)*(1.490-1.333)/(1.544-1.333)


  return
end function rxsea

real function drwsea(rh)
  implicit none
  real :: rh
  
  DRWSEA=1.00212-0.001625*RH+0.00131*RH/(1.0+0.928*RH**3)


  return
end function drwsea

real function drdsea(rh)
  implicit none
  real :: rh
  
  DRDSEA=1.003007     ! ratio of wet & dry nr(0.550) / nr(0.633)


  return
end function drdsea

real function awno3(x)
  implicit none
  real :: x
  !        Nitrate parametric formulas from Tang & Munkelwitz(94,96)
  AWNO3=1.0-0.365*X-0.09155*X**2-0.2826*X**3      ! TM96
  return
end function awno3
real function dwno3(x)
  implicit none
  real :: x
  DWNO3=    -0.365  -0.1831*X   -0.8478*X**3
  return
end function dwno3
real function rono3(x)
  implicit none
  real :: x
  RONO3=0.9971+0.405*X+0.090*X**2                   ! TM96
  return
end function rono3
real function bxno3(x)
  implicit none
  real :: x
  real, external :: rono3
  BXNO3=(1.0/X*1.725/RONO3(X))**(1.0/3.0)             ! TM96
  return
end function bxno3
real function r1no3(x)
  implicit none
  real :: x
  R1NO3=1.3330+0.119*X          !  (X<0.205)              TWM81
  return
end function r1no3
real function r2no3(x)
  implicit none
  real :: x
  R2NO3=1.3285+0.145*X          !  (X>0.205)              TWM81
  return
end function r2no3
real function drxno3(rh)
  implicit none
  real :: rh
  DRXNO3=1.001179     ! ratio of wet & dry nr(0.550) / nr(0.633)
  return
end function drxno3
  
real function awocx(x)
  implicit none
  real :: x
  !        Organic Carbon - adapted from Sulfate parametric formulas
  !        yields growth factor G=1.1 at RH=0.84 Virkkula et al 1999
  AWOCX=1.0-X**8
  return
end function awocx
real function dwocx(x)
  implicit none
  real :: x
  DWOCX=-8.0*X**7
  return
end function dwocx
real function roocx(x)
  implicit none
  real :: x
  ROOCX=1.0+.5*X
  return
end function roocx
real function bxocx(x)
  implicit none
  real :: x
  real, external :: roocx
  BXOCX=(1.5/(X*ROOCX(X)))**(1.0/3.0)
  return
end function bxocx
real function rxocx(x)
  implicit none
  real :: x
  RXOCX=1.3330+.193*X
  return
end function rxocx
real function drwocx(rh)
  implicit none
  real :: rh
  DRWOCX=1.00253-0.00198*RH+0.00184*RH/(1.0+0.656*RH**1.1)
  return
end function drwocx
real function drdocx(rh)
  implicit none
  real :: rh
  DRDOCX=1.00253
  return
end function drdocx

!===========================================================

SUBROUTINE SPLNI4(Q,NI,NJ,IR,JN,DN,QQ)
  IMPLICIT NONE
  
  integer, intent(in)  ::   NI,NJ,  IR, JN
  real,  intent(in)  :: Q(NI,NJ),     DN
  real,  intent(out) :: QQ
  
  !nu   real*8,save :: CUSPWM=1., CUSPWE=1. ,CUSPWT,fflinr
  real ffcusp,a,b,c,d,xf,xe,xexm
  REAL f1,f2,f3,f4, f21,f32,f43,f3221,f4332
  
  F1=Q(IR,JN-1)
  F2=Q(IR,JN)
  F3=Q(IR,JN+1)
  F4=Q(IR,JN+2)
  F21=F2-F1
  F32=F3-F2
  F43=F4-F3
  F3221=0.5*(F32+F21)
  F4332=0.5*(F43+F32)
  A=F2
  B=F3221
  C=3.0*F32-F3221-F3221-F4332
  D=F3221+F4332-F32-F32
  XF=DN
  FFCUSP=A+XF*(B+XF*(C+XF*D))
  XE=1.0-XF-XF
  IF(XE.LT.0.0) XE=-XE
  XEXM=XE**2
  !=1   CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
  !nu   FFLINR=A+XF*F32
  QQ=FFCUSP ! *CUSPWT+FFLINR*(1.0-CUSPWT)
  RETURN
END SUBROUTINE SPLNI4

!=======================================================================
SUBROUTINE SPLN44(Q,NI,NJ,IR,DR,JN,DN,QQ)
  IMPLICIT NONE
  
  integer, intent(in)  ::   NI,NJ,  IR, JN
  real,  intent(in)  :: Q(NI,NJ), DR, DN
  real,  intent(out) :: QQ
  
  !nu   real*8,save :: CUSPWM=1., CUSPWE=1. ,CUSPWT,fflinr
  real QK(4),  ffcusp,a,b,c,d,xf,xe,xexm
  REAL f1,f2,f3,f4, f21,f32,f43,f3221,f4332
  integer k,kr,irm,irp

  K=0
  IRM=IR-1
  IRP=IR+2
  DO KR=IRM,IRP
     K=K+1
     F1=Q(KR,JN-1)
     F2=Q(KR,JN)
     F3=Q(KR,JN+1)
     F4=Q(KR,JN+2)
     F21=F2-F1
     F32=F3-F2
     F43=F4-F3
     F3221=0.5*(F32+F21)
     F4332=0.5*(F43+F32)
     A=F2
     B=F3221
     C=3.0*F32-F3221-F3221-F4332
     D=F3221+F4332-F32-F32
     XF=DN
     FFCUSP=A+XF*(B+XF*(C+XF*D))
     XE=1.0-XF-XF
     IF(XE.LT.0.0) XE=-XE
     XEXM=XE**2
     !=1   CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
     !nu   FFLINR=A+XF*F32
     QK(K)=FFCUSP ! *CUSPWT+FFLINR*(1.0-CUSPWT)
  enddo
  F1=QK(1)
  F2=QK(2)
  F3=QK(3)
  F4=QK(4)
  F21=F2-F1
  F32=F3-F2
  F43=F4-F3
  F3221=0.5*(F32+F21)
  F4332=0.5*(F43+F32)
  A=F2
  B=F3221
  C=3.0*F32-F3221-F3221-F4332
  D=F3221+F4332-F32-F32
  XF=DR
  FFCUSP=A+XF*(B+XF*(C+XF*D))
  XE=1.0-XF-XF
  IF(XE.LT.0.0) XE=-XE
  XEXM=XE**2
  !=1   CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
  !nu   FFLINR=A+XF*F32
  QQ=FFCUSP ! *CUSPWT+FFLINR*(1.0-CUSPWT)
  RETURN
END SUBROUTINE SPLN44

!==================================================================

SUBROUTINE SPLINV(X,F,NXF,XX,FF,CUSPWM,CUSPWE,KXTRAP)
  IMPLICIT NONE
  
  integer, intent(in)  :: NXF,              KXTRAP
  real , intent(in)  :: X(NXF),F(NXF),FF, CUSPWM,CUSPWE
  real , intent(out) :: XX

!---------------------------------------------------------------------
!    Inverse spline:
!    SPLINV locates FF between points (F2,X2)(F3,X3) on 4-point spread
!    and returns 4-point Cubic Spline value of XX such that FF = F(XX)
!
!    Quadratic Derivatives of Spline are continuous at (F2,X2),(F3,X3)
!    (X-Coordinate may be specified in increasing or decreasing order)
!
!---------------------------------------------------------------------
!
!    Custom Control Parameters:  CUSPWM,CUSPWE
!------------------------------
!
!    In cases where data points are unevenly spaced and/or data points
!    exhibit abrupt changes in value, Spline Interpolation may produce
!    undesirable bulging of interpolated values. In more extreme cases
!    Linear Interpolation may be less problematic to use.
!
!    Interpolation can be weighted between: Cubic Spline and Linear by
!    adjusting weights CUSPWM and CUSPWE to values between 1.0 and 0.0
!
!    CUSPWM = Cubic Spline Weight at the (X2-X3) Interval Mid-point
!    CUSPWE = Cubic Spline Weight at the (X2-X3) Interval End-points
!
!    For example, with:
!
!    CUSPWM=1.0,CUSPWE=1.0  FF returns Cubic Spline interpolated value
!    CUSPWM=0.0,CUSPWE=0.0  FF returns   Linearly   interpolated value
!
!---------------------------------------------------------------------
!
!     Extrapolation for XX outside of defined interval:  X(1)<->X(NXF)
!
!               KXTRAP = 0    No Extrapolation   (i.e., sets XX = 0.0)
!                        1    Fixed Extrapolation (sets XX=edge value)
!                        2    Linear Extrapolation using 2 edge points
!
!---------------------------------------------------------------------
!
!
!    NOTE:  F(X) is assumed to be monotonic between F(1) and F(NXF)
!
!------------------------------------------------------------------

  REAL x1,x2,x3,x4, x21,x32,x43,x31,x42, BETW,FFCUSP,FFLINR,CUSPWT
  REAL f1,f2,f3,f4, f21,f32,f43,f3221,f4332, a,b,c,d, xf,xe,xexm
  REAL DX,gg,xg,xy,deltx,slopec,slopel,slopes
  integer k,kk

  BETW=(F(2)-FF)*(F(NXF)-F(1))
  IF(BETW.GT.0.0) GO TO 120
  BETW=(FF-F(NXF-1))*(F(NXF)-F(1))
  IF(BETW.GT.0.0) GO TO 140
  
  DO K=3,NXF-1
     BETW=(FF-F(K-1))*(F(K)-FF)
     DX=(FF-F(K-1))/(F(K)-F(K-1))
     XX=X(K-1)+DX*(X(K)-X(K-1))
     IF(BETW.GE.0.0) GO TO 110
  enddo
     
110 CONTINUE
  DO KK=1,5
     X1=X(K-2)
     X2=X(K-1)
     X3=X(K)
     X4=X(K+1)
     F1=F(K-2)
     F2=F(K-1)
     F3=F(K)
     F4=F(K+1)
     X21=X2-X1
     X31=X3-X1
     X32=X3-X2
     X43=X4-X3
     X42=X4-X2
     F21=(F2-F1)/(X21*X21)
     F32=(F3-F2)/(X32*X32)
     F43=(F4-F3)/(X43*X43)
     F3221=(F32+F21)/X31*X21
     F4332=(F43+F32)/X42*X43
     A=F2
     B=X32*F3221
     C=3.0*F32-F3221-F3221-F4332
     D=(F3221+F4332-F32-F32)/X32
     XF=XX-X2
     
     !                             FFCUSP= Cubic Spline Interpolation Result
     !                             -----------------------------------------
     
     FFCUSP=A+XF*(B+XF*(C+XF*D))
     XE=(X3+X2-XX-XX)/X32
     IF(XE.LT.0.0) XE=-XE
     XEXM=XE**2
     CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
     
     !FFLINR= Linear Interpolation Result
     !                                   -----------------------------------
     FFLINR=A+XF*F32*X32
     GG=FFCUSP*CUSPWT+FFLINR*(1.0-CUSPWT)
     SLOPEC=B+2.0*C*XF+3.0*D*XF**2
     SLOPEL=F32*X32
     SLOPES=SLOPEC*CUSPWT+SLOPEL*(1.0-CUSPWT)
     XG=XF
     XY=XX
     DELTX=(GG-FF)/SLOPES
     XX=XF-(GG-FF)/SLOPES+X2
  enddo
  GO TO 160
     
  !Edge Point Interval Interpolation and/or Extrapolation
  !------------------------------------------------------
120 CONTINUE
  BETW=(F(1)-FF)*(F(NXF)-F(1))
  IF(BETW.GT.0.0) GO TO 130
  
  !                          F(1),F(2)  Edge Point Interval Interpolation
  !                          --------------------------------------------
  DO KK=2,6
     X1=X(1)
     X2=X(2)
     X3=X(3)
     F1=F(1)
     F2=F(2)
     F3=F(3)
     XX=X1+(FF-F(1))/(F(2)-F(1))*(X2-X1)
     XF=XX-X1
     X21=X2-X1
     F21=(F2-F1)/X21
     X32=X3-X2
     X31=X3-X1
     C=((F3-F2)/X32-F21)/X31
     B=F21-X21*C
     A=F1
     FFCUSP=A+XF*(B+XF*C)
     FFLINR=A+XF*F21
     XE=1.0-2.0*XF/X21
     XEXM=XE**2
     CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
     GG=FFCUSP*CUSPWT+FFLINR*(1.0-CUSPWT)
     SLOPEC=B+2.0*C*XF
     SLOPEL=F21
     SLOPES=SLOPEC*CUSPWT+SLOPEL*(1.0-CUSPWT)
     XG=XF
     DELTX=(GG-FF)/SLOPES
     XX=XF-(GG-FF)/SLOPES+X1
  enddo
  GO TO 160
  
130 CONTINUE
  !                  Extrapolation for FF Outside of Interval F(1) - F(2)
  !                  ----------------------------------------------------
  !                  IF(KXTRAP.EQ.0)  (No Extrapolation:   sets XX = 0.0)
  !                  IF(KXTRAP.EQ.1)  (Extrapolation at Fixed Edge Value)
  !                  IF(KXTRAP.EQ.2)  (2 Edge Point Linear Extrapolation)

  IF(KXTRAP.EQ.0) XX=0.0
  IF(KXTRAP.EQ.1) XX=X(1)
  IF(KXTRAP.EQ.2) XX=X(1)-(F(1)-FF)/(F(2)-F(1))*(X(2)-X(1))
  GO TO 160
  
140 CONTINUE
  BETW=(FF-F(NXF))*(F(NXF)-F(1))
  IF(BETW.GT.0.0) GO TO 150
  
  !                    F(NXF-1),F(NXF)  Edge Point Interval Interpolation
  !                    --------------------------------------------------
  DO KK=3,7
     X1=X(NXF-2)
     X2=X(NXF-1)
     X3=X(NXF)
     F1=F(NXF-2)
     F2=F(NXF-1)
     F3=F(NXF)
     XX=X2+(FF-F2)/(F3-F2)*(X3-X2)
     XF=XX-X2
     X32=X3-X2
     F32=(F3-F2)/X32
     X21=X2-X1
     X31=X3-X1
     F21=(F2-F1)/X21
     
     !                    3-Point Quadratic Interpolation for Edge Intervals
     !                    --------------------------------------------------
     !
     !      (Edge Option)     ----------------------------------------------
     !                        For Linear Interpolation within Edge Intervals
     !                        between F(1),F(2), and between F(NXF-1),F(NXF)
     !                        set the value of coefficient C below, to C=0.0
     !                        ----------------------------------------------
     
     C=(F32-F21)/X31
     B=F21+X21*C
     A=F2
     FFCUSP=A+XF*(B+XF*C)
     FFLINR=A+XF*F32
     XE=1.0-2.0*XF/X32
     IF(XE.LT.0.0) XE=-XE
     XEXM=XE**2
     CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
     GG=FFCUSP*CUSPWT+FFLINR*(1.0-CUSPWT)
     SLOPEC=B+2.0*C*XF
     SLOPEL=F21
     SLOPES=SLOPEC*CUSPWT+SLOPEL*(1.0-CUSPWT)
     XG=XF
     DELTX=(GG-FF)/SLOPES
     XX=XF-(GG-FF)/SLOPES+X2
  enddo
  GO TO 160
  
150 CONTINUE
  !              Extrapolation for F Outside of Interval  F(NXF-1)-F(NXF)
  !              --------------------------------------------------------
  !                  IF(KXTRAP.EQ.0)  (No Extrapolation:   sets XX = 0.0)
  !                  IF(KXTRAP.EQ.1)  (Extrapolation at Fixed Edge Value)
  !                  IF(KXTRAP.EQ.2)  (2 Edge Point Linear Extrapolation)
  
  IF(KXTRAP.EQ.0) XX=0.0
  IF(KXTRAP.EQ.1) XX=X(NXF)
  IF(KXTRAP.EQ.2) XX=X(NXF)  &
       -(F(NXF)-FF)/(F(NXF-1)-F(NXF))*(X(NXF-1)-X(NXF))
  
160 CONTINUE
  RETURN
END SUBROUTINE SPLINV
