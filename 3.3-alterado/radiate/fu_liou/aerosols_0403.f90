      subroutine aerosol_init
!
!                        8/14/95, 4/1/97 , 2/10/2000
!
!  **********************************************************************
!  Subroutine to create aerosol optical properties.  There are several
!  inputs and 6 outputs.
!
!    INPUTS FROM COMMON BLOCKS OR HEADER FILE:
!
!    a_tau(nwi) :  The input column aerosol optical depth
!    (real)           (common block "aer_tau" - see header file).
!
!    a_wli(nwi) :  Wavelength in microns corresponding to aerosol tau in "a_tau"
!
!    aprof(# layers): The input aerosol optical depth profile - LAYERS
!    (real)           (common block "aer_prof").
!
!    itp:       Aerosol type, given in header file rad_0598.h.
!
!    ifg:       The table will compute vertical distributions based on
!    (integer)  relative humidity (see explanation below).  If ifg is
!               set to 0, each layer will have properties calculated
!               based on the relative humidity of that layer.  If ifg
!               is set equal to another integer (1 through the number of
!               relative humidities given in the block data "aerosol")
!               the routine will calculate a vertical profile of optical
!               properties based on the relative humidity corresponding
!               to the index given.  The indices are: 1: 0%; 2: 50%;
!               3: 70%; 4: 80%; 5:90%; 6: 95%; 7: 98%; and 8: 99%.
!               If the number of relative humidities changes, these
!               numbers will have to be modified.
!
!    ivd:       Vertical tau distribution flag.  If set to zero, the
!               distribution is based on Jim Spinhirne's marine
!               distribution formulation, and no user input is required.
!               If set to one, the user's own vertical distribution is
!               used, and must be present in the array aprof(nlayers).
!               NOTE: This vertical distribution is used as a weighting
!               factor ONLY, to distribute input column optical depths!
!
!----------------------------------------------------------------------------
!    a_ssa, a_ext, a_asy:  Input single-scattering albedos, extinction
!           coefficients, and asymmetry parameters.  These variables
!           are dimensioned (# of bands, # of relative humidities,
!           # of aerosol types). An x or y is appended on these
!           variable names: if x, the numbers correspond to the 18
!           original bands.  If y, the numbers are for the 10
!           sub-intervals in the first shortwave band (.2-.7 microns).
!           All of these variables come from the block data statements
!           aerosol# (# corresponds to an integer, eg. aerosol1) and
!           are in common blocks aer_optx and aer_opty.
!
!    nv,mb,pp,pt,ph,dz: number of layers, number of bands, and the
!           pressure, temperature, humidity and thickness profiles.
!           These are shared by several subroutines.
!
!    OUTPUTS:
!
!    a_tau1,a_ssa1,a_asy1:  The optical depth, single-scattering albedo,
!       and asymmetry parameter vertical profiles for 18 bands.  These
!       are dimensioned (nvx, 18)  These are in the common block
!       aer_initx, which is shared by the subroutine "aerosolx".
!
!    a_tau2,a_ssa2,a_asy2:  Properties for SW band 1's 10 subintervals.
!       These are dimensioned (nvx, 10)  These are in the common block
!       aer_inity, which is shared by the subroutine "aerosoly".
!
!  **********************************************************************
      USE FUINPUT
      USE EXTRAS,only : sh_rh
      implicit none

      integer iq,mtop,n,m,ict,ix,iy,irh,krh,iac,itp
      real, dimension(mbx,nrh,naer) :: a_ssax,a_extx,a_asyx
      real, dimension(mby,nrh,naer) :: a_ssay,a_exty,a_asyy

      real, dimension(nvx) :: tauxxx
      real, dimension(nvx,mbx,mxac) :: a_tau1,a_ext1,a_ssa1,a_asy1
      real, dimension(nvx,mby,mxac) :: a_tau2,a_ext2,a_ssa2,a_asy2

      real ,dimension(nvx)  :: taux1,taux2,rh,ht,rhp
      real sumxxx

      real,dimension(mxat) :: a_tau
      real,dimension(nvx)  :: aprof
      real,dimension(nvx,mbx) :: wvd_x
      real,dimension(nvx,mby) :: wvd_y

      real p1,h1,z,sig,tp
      real rhx(nrh)
      real wts(4),tau3(2),tau3y(4)
      real aotf,wlf,sump,rirh
      real spinhirne_sig, spinhirne_tau

      common /aer_optx/ a_ssax,a_extx,a_asyx
      common /aer_opty/ a_ssay,a_exty,a_asyy
      common /aer_initx/ a_tau1,a_ssa1,a_asy1
      common /aer_inity/ a_tau2,a_ssa2,a_asy2
!      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)

      common /tau_spline_aot/ aotf(15),wlf(15)

      data rhx /0.,50.,70.,80.,90.,95.,98.,99./
      data wts /.23015,.28274,.25172,.23539/


!  Initialize.

       rh     = -9999.
       a_ssa1 = 0. ; a_ext1 = 0. ; a_asy1 = 0. ; a_tau1 = 0.
       a_ssa2 = 0. ; a_ext2 = 0. ; a_asy2 = 0. ; a_tau2 = 0.

!     if ( nac < 0 .or. nac > mxac ) fi%ierror= 11 ! print*, 'nac:# Aerosol Constituents'
!     if (n_atau<0 .or.n_atau>mxat) fi%ierror= 12 !print*, 'n_atau:# Aerosol Tau / Wavelengths'
!     if (ifg < 0 .or. ifg > 8)  fi%ierror= 13 ! print*, 'ifg: Aerosol RH% Flag'
      AEROSOL_CONSTITUENTS : do iac = 1,nac

!!!   a_wli(1:n_atau) = a_wlis(1:n_atau,iac)
      a_tau(1:n_atau) = a_taus(1:n_atau,iac)
      aprof(1:nv)  = aprofs(1:nv,iac)
      itp       = itps(iac)
!     print*, iac,itp
!     if ( itp < 1 .or. itp > naer ) fi%ierror= 14 ! print*, ' itp : Bad Aerosol Type'
!     print*,'CONSTITUENTS',iac,itp

! FOR Aerosol Optical Properties types that are constant with RH
      if (itp==1  .or. itp==2 .or. itp==3 .or. &
            itp==10 .or. itp==12 .or.itp==13 .or. itp==18 ) then
!!       Has already been filled in Block data
      else
      do krh=2,8
       a_extx(1:mbx,krh,itp)= a_extx(1:mbx,1,itp)
       a_ssax(1:mbx,krh,itp)= a_ssax(1:mbx,1,itp)
       a_asyx(1:mbx,krh,itp)= a_asyx(1:mbx,1,itp)

       a_exty(1:mby,krh,itp)= a_exty(1:mby,1,itp)
       a_ssay(1:mby,krh,itp)= a_ssay(1:mby,1,itp)
       a_asyy(1:mby,krh,itp)= a_asyy(1:mby,1,itp)

      enddo

      endif
!     if ( ifg .ne.0) print*,'CHECK',ifg,itp,a_ssax(1:mbx,ifg,itp)




!  ******************************************************************
!  Calculate heights at center of layer - find highest layer to place
!  aerosols (15 km) - calculate relative humidities of each layer as
!  needed.  Values of RH > 99% will be set equal to 99% to make table
!  lookup easier. "mtop" is the highest aerosol layer.
!  ******************************************************************
      z=0.
      m=nv
      iq=0
      do while (iq.eq.0)
       ht(m)=(z*2.+dz(m))/2.
       z=z+dz(m)
       if (z.gt.15.) then
         iq=1
         mtop=m
       endif
       p1=(pp(m)+pp(m+1))/2.
       tp=(pt(m)+pt(m+1))/2.
       h1=(ph(m)+ph(m+1))/2.
!       call ql_rh(rh(m),tp,p1,h1)
       rh(m)= sh_rh(tp,p1,h1)
       if (rh(m).gt.98.9) rh(m)=98.9
       if ((rh(m).lt..01).and.(rh(m).gt.-999.)) rh(m)=0.
       m=m-1
       end do

!  *************************************************************
!  Calculate vertical distribution of asymmetry, ss albedo and
!  extinction, based on aerosol type and relative humidity.
!  If ifg is not equal to 0, parameters  will corresponds to a
!  single RH, as described in header file. Loop 31 deals with
!  the 18 original bands, loop 32 with the 10 band 1 subintervals.
!  *************************************************************
      do 30 m=mtop,nv
       do 31 n=1,mbx
        if (rh(m).eq.-9999.) then
          a_ext1(m,n,iac)=-9999.
          a_ssa1(m,n,iac)=-9999.
          a_asy1(m,n,iac)=-9999.
        else
          if (ifg.eq.0) then          ! Dependence on layer RH.
            ict=2
            do while (rh(m).ge.rhx(ict))
             ict=ict+1
             end do
            a_ext1(m,n,iac)=a_extx(n,ict-1,itp)+(rh(m)-rhx(ict-1))/ &
           (rhx(ict)-rhx(ict-1))*(a_extx(n,ict,itp)-a_extx(n,ict-1,itp))
            a_ssa1(m,n,iac)=a_ssax(n,ict-1,itp)+(rh(m)-rhx(ict-1))/ &
           (rhx(ict)-rhx(ict-1))*(a_ssax(n,ict,itp)-a_ssax(n,ict-1,itp))
            a_asy1(m,n,iac)=a_asyx(n,ict-1,itp)+(rh(m)-rhx(ict-1))/ &
           (rhx(ict)-rhx(ict-1))*(a_asyx(n,ict,itp)-a_asyx(n,ict-1,itp))
        rhp(m) = rh(m)
          else                        ! Dependence on prescribed RH.
       ict=ifg
            a_ext1(m,n,iac)=a_extx(n,ict,itp)
            a_ssa1(m,n,iac)=a_ssax(n,ict,itp)
            a_asy1(m,n,iac)=a_asyx(n,ict,itp)
          endif
        endif
 31     continue
!-------------------------------------------
       do 32 n=1,mby
        if (rh(m).eq.-9999.) then
          a_ext2(m,n,iac)=-9999.
          a_ssa2(m,n,iac)=-9999.
          a_asy2(m,n,iac)=-9999.
        else
          if (ifg.eq.0) then          ! Dependence on layer RH.
            ict=2
            do while (rh(m).ge.rhx(ict))
             ict=ict+1
             end do
            a_ext2(m,n,iac)=a_exty(n,ict-1,itp)+(rh(m)-rhx(ict-1))/ &
           (rhx(ict)-rhx(ict-1))*(a_exty(n,ict,itp)-a_exty(n,ict-1,itp))
            a_ssa2(m,n,iac)=a_ssay(n,ict-1,itp)+(rh(m)-rhx(ict-1))/ &
           (rhx(ict)-rhx(ict-1))*(a_ssay(n,ict,itp)-a_ssay(n,ict-1,itp))
            a_asy2(m,n,iac)=a_asyy(n,ict-1,itp)+(rh(m)-rhx(ict-1))/ &
           (rhx(ict)-rhx(ict-1))*(a_asyy(n,ict,itp)-a_asyy(n,ict-1,itp))
          else                        ! Dependence on prescribed RH.
       ict=ifg
            a_ext2(m,n,iac)=a_exty(n,ict,itp)
            a_ssa2(m,n,iac)=a_ssay(n,ict,itp)
            a_asy2(m,n,iac)=a_asyy(n,ict,itp)
          endif
        endif
 32     continue

 30    continue

!  ******************************************************************
!  Vertical distribution of aerosol optical depths - CAGEX and CERES.
!       --------------------------------------------------------------
!       Use Spinhirne's vertical distribution of scattering properties
!       to calculate vertical distribution of optical depths.  The
!       distribution gives a scattering coefficient ("sig"). Use this,
!       along with the single-scattering albedo, to produce an
!       RH-dependent extinction coefficient (extx, exty, etc.), from
!       which optical depth is calculated (taux, tauy, etc.).  This
!       optical depth is summed (sum1, sumy2, sum, etc.) to give
!       column tau for weighting purposes.
!       --------------------------------------------------------------

      select case (ivd)
      case default
       fi%ierror = 15 ! print*, ' ivd : Aerosol Profile flag'
      case (0)  !! DEFAULT VERTICAL DISTRIBUTION Spinhirne

      sumxxx=0.0

        do  m=mtop,nv

       sig = spinhirne_sig( ht(m))
       tauxxx(m) = spinhirne_tau(sig,a_ssa2(m,9,iac),dz(m))
       sumxxx   = sumxxx + tauxxx(m)
!      print*,m,sig,a_ssa2(m,9,iac)
        enddo

      do m=mtop,nv
       tauxxx(m) = tauxxx(m)  / sumxxx
!!!    aprofs(m,iac) = tauxxx(m) !! See what the Sphinhirne profiles look like
      enddo

! ----------------------------------------------------------------
      case (1)   ! USER'S OWN VERTICAL DISTRIBUTION IVD=1
         sump =   sum( aprof(mtop:nv) )
       tauxxx(mtop:nv)= aprof(mtop:nv) / sump

         if(sump.eq.0.)fi%ierror = 16 !print*, 'No VERTICAL Profile OF A

      end select

!  ********************************************************************
!  IAFORM=2
!
!  Distribute optical depth spectrally into the first 2 Fu-Liou bands.
!  Band 1 will consist of the first 4 MFRSR bands, weighted with
!  respect to energy.  Band two will be the fifth MFRSR band.
!
!  Also, distribute optical depths into 4 of the 10 band 1 subintervals.
!  Subinterval 7 is directly inserted, since there is one MFRSR
!  measurement within the range of this band.  Subintervals 7 and 8
!  straddle the .497 micron MFRSR measurement, so interpolated values
!  are inserted into these, using .409 and .497 measurements for 7, and
!  .497 and .606 for 8.  Subinterval 10 contains two MFRSR measurements,
!  so it is filled using an energy-weighted average.  This is all
!  hardwired, so we need all of the MFRSR bands (.409, .497, .606, and
!  .661) for it to work. (The .855 micron band is also needed, but not
!  for this interval distribution.
!  ********************************************************************

      select case ( iaform )
       case default
       fi%ierror = 17 ! print*, ' iaform : Bad value of iaform '
      case(1)        ! CERES
!! No operations necessary

      case(2)        ! For CAGEX

        tau3(1)=a_tau(1)*wts(1)+a_tau(2)*wts(2)+ &
                a_tau(3)*wts(3)+a_tau(4)*wts(4)
        tau3(2)=a_tau(5)
        tau3y(1)=a_tau(1)      ! For subinterval 7 of 1st band (.409)
        tau3y(2)=a_tau(1)+.6705*(a_tau(2)-a_tau(1)) ! Subi 8 of band 1
        tau3y(3)=a_tau(2)+.4541*(a_tau(3)-a_tau(2)) ! Subi 9 of band 1
        tau3y(4)=a_tau(3)*.5175+a_tau(4)*.4825      ! Subi 10 of band 1

      case(3)        ! For AOT_SPLINEFIT

      if ( ifg == 0 ) then ! Find Aerosol weighted collumn mean RH index
       rirh=0
       do m =mtop,nv
       rirh = rirh + rhp(m)* tauxxx(m)  !! Aerosol Profile weighted mean
!      print*,m,rhp(m),tauxxx(m)
       enddo

         irh =1
         do ix= 1,7
         if( rirh >= rhx(ix) .and. rirh < rhx(ix+1) ) irh=ix
         enddo
         if( rirh >= rhx(8) )irh =8

      else  ! Use assigned RH index
       irh = ifg
      endif

! Can't handle ZERO in Log interpolation
      where ( a_tau(1:n_atau) .lt. 1.0E-20) a_tau(1:n_atau) = 1.0E-20

       call atau_spline_iaform3(mxat,n_atau,a_wli,a_tau,itp,irh)

!     write(22,'(a20,15f8.3)') 'AOT in Fu Bands',aotf(1:15)

!!! ACCOUNT FOR VERTICAL EXTINCTION VARIABILITY WITH HUMIDITY ABOUT THE MEAN RH "irh"
!!! ( IAFORM==3) only
      do iy = 1,mby
      wvd_y(mtop:nv,iy)=tauxxx(mtop:nv) &
                   *a_ext2(mtop:nv,iy,iac)/a_exty(iy,irh,itp)
      sump =   sum( wvd_y(mtop:nv,iy) )
      wvd_y(mtop:nv,iy) =  wvd_y(mtop:nv,iy) /sump
      enddo

      do ix = 1,mbx
      wvd_x(mtop:nv,ix)=tauxxx(mtop:nv) &
                   *a_ext1(mtop:nv,ix,iac)/a_extx(ix,irh,itp)
      sump =   sum( wvd_x(mtop:nv,ix) )
      wvd_x(mtop:nv,ix) =  wvd_x(mtop:nv,ix) /sump
      enddo

      end select


! ----------------------------------------------------------------
!       Use weighted optical depths  to distribute our input
!       column optical depths vertically and spectrally where needed.
!       For bands with "measured" input, we simply do the weighting.
!       For the remaining bands, we weight according to our vertically
!       distributed extinction coefficients (calculated in loop 30),
!       which carry all the spectral resolution we need.  a_tau1 is for
!       the 18 original bands, a_tau2 is for the 10 band 1 subintervals.
! ----------------------------------------------------------------
        VERTICAL : do  m=mtop,nv

      select case ( iaform )

      case(1)       ! For CERES

           a_tau1(m,1,iac)   = a_tau(1) * tauxxx(m)
           a_tau1(m,2:18,iac)= a_tau1(m,1,iac)* &
                       a_ext1(m,2:18,iac)/a_ext1(m,1,iac)

           a_tau2(m,9,iac)  = a_tau(1) * tauxxx(m)

           a_tau2(m,1:10,iac)=a_tau2(m,9,iac)* &
                       a_ext2(m,1:10,iac)/a_ext2(m,9,iac)

        case(2)        ! For CAGEX

          a_tau1(m,1:2,iac) = tau3(1:2) * tauxxx(m)
            a_tau1(m,3:18,iac)=a_tau1(m,2,iac)* &
                      a_ext1(m,3:18,iac)/a_ext1(m,2,iac)

          a_tau2(m,7:10,iac) = tau3y(1:4) * tauxxx(m)
            a_tau2(m,1:6,iac)  = a_tau2(m,7,iac)* &
                      a_ext2(m,1:6,iac)/a_ext2(m,7,iac)

       case(3)       ! For AOT_SPLINEFIT



!     a_tau2(m,1:10,iac) = aotf(1:10)  * tauxxx(m)
      a_tau2(m,1:10,iac) = aotf(1:10)  * wvd_y(m,1:10)
!     a_tau1(m,1,iac)    = aotf(9)     * tauxxx(m)
      a_tau1(m,1,iac)    = aotf(9)     * wvd_x(m,1)
!     a_tau1(m,2:6,iac)  = aotf(11:15) * tauxxx(m)
      a_tau1(m,2:6,iac)  = aotf(11:15) * wvd_x(m,2:6)
            a_tau1(m,7:18,iac) =a_tau1(m,2,iac)* &
                       a_ext1(m,7:18,iac)/a_ext1(m,2,iac)

         end select

!     print'(3I4,2f8.2,16f7.3)', m,iac,itp,dz(m),rh(m),
!     & (wvd_y(m,iy),iy=1,10),(wvd_x(m,ix),ix=1,6)

        enddo VERTICAL

!------------------------------------------------------------------------------
!!!--- Diagnostic Output of Atau
!     do ii=1,10
!     xxx=0
!      do jj=1,nv
!      xxx =xxx+ a_tau2(jj,ii,iac)
!      enddo
!     aotf(ii)=xxx
!     enddo

!     do ii=2,6
!     xxx=0
!      do jj=1,nv
!      xxx =xxx+ a_tau1(jj,ii,iac)
!      enddo
!     aotf(9+ii)=xxx
!     enddo

!     write(22,'(a20,15f8.3)') 'AOT in Fu Bands',aotf(1:15)

      enddo AEROSOL_CONSTITUENTS

      return
      end

!===========================================================================
      subroutine aerosolxy ( ib,cmode )
! *********************************************************************
!                      Modified 2/14/00
!
! tae, wae, and wwae are the optical depth, single scattering albedo,
! and expansion coefficients of the phase function ( 1, 2, 3, and 4 )
! due to the Mie scattering of aerosols for a given layer.
!
!  This subroutine is called for bands 2 - 18 (ib)
!  or vis subbands 1-10 (ig)
! *********************************************************************
      USE FUINPUT
      implicit none
      character*1 cmode
      integer i,ib,iac
      real x1,x2,x3,x4,y1,y2,y3,y4,tae,wae,wwae
      real ,dimension(nvx,18,mxac) :: a_tau1,a_ssa1,a_asy1
      real ,dimension(nvx,10,mxac) :: a_tau2,a_ssa2,a_asy2
      common /aer_initx/ a_tau1,a_ssa1,a_asy1
      common /aer_inity/ a_tau2,a_ssa2,a_asy2

      common /aer/ tae(nvx,mxac), wae(nvx,mxac), wwae(nvx,4,mxac)

      AEROSOL_CONSTITUENTS  : do iac=1,nac

      LEVELS : do  i = 1, nv
       select case (cmode)
      case ('x')
         tae(i,iac) = a_tau1(i,ib,iac)
         wae(i,iac) = a_ssa1(i,ib,iac)
         x1     = a_asy1(i,ib,iac)
      case ('y')
         tae(i,iac) = a_tau2(i,ib,iac)
         wae(i,iac) = a_ssa2(i,ib,iac)
         x1     = a_asy2(i,ib,iac)
       end select

       x2 = x1 * x1
       x3 = x2 * x1
       x4 = x3 * x1
       y1 = 3.0 * x1
       y2 = 5.0 * x2
       y3 = 7.0 * x3
       y4 = 9.0 * x4

       wwae(i,1,iac) = y1
       wwae(i,2,iac) = y2
       wwae(i,3,iac) = y3
       wwae(i,4,iac) = y4

      enddo LEVELS
      enddo AEROSOL_CONSTITUENTS

      return
      end
!----------------------------------------------------------------
      real function spinhirne_sig(ht)

      data sig0,a,ap,b,bp,f/0.025,0.4,2981.0,1.6,2.5,1.5e-7/

         t1=  sig0*(1+a)**2
         t4 = f*(1+ap)**2

         t2 = exp(ht/b)
         t3 = (a+exp(ht/b))**2
         t5 = exp(ht/bp)
         t6 = (a+exp(ht/bp))**2
         spinhirne_sig=t1*t2/t3+t4*t5/t6   ! scattering coefficient

      return
      end
!---------------------------------------------
      real function spinhirne_tau(sig,ssa,dz)
      ext = sig / ssa
      spinhirne_tau = ext / dz
      return
      end
!=========================================================================
      subroutine atau_spline_iaform3(mxat,nwi,wli,aoti,ityp,irh)
      parameter(nsub=5 ,nfuo=15,nwo=nsub*nfuo)
      common /aot_spect_5/  wlo2(5,15) , hkas(5,15) ,sflx (5,15)
!     parameter(nsub=25,nfuo=15 ,nwo=nsub*nfuo)
!     common /aot_spect_25/ wlo2(25,15) , hkas(25,15) ,sflx (25,15)  !!! Higer resolution Convolutio

      real ,dimension(mxat):: aoti,wli
      real ,dimension(nwo):: aoto,wlo
      real , dimension(nsub,nfuo) :: aoto2

      common /tau_spline_aot/ aotf(15),wlf(15)

!     wlo = reshape(wlo2,(/nwo/))
      kk=0
      do jj=1,nfuo
      do ii=1,nsub
      kk=kk+1
      wlo(kk) = wlo2(ii,jj)
      enddo;enddo


      call aot_ext &
         (mxat,nwi,aoti,wli,nwo,wlo,aoto,ityp,irh)

!     aoto2 = reshape(aoto,(/5,15/))
      kk=0
      do jj=1,nfuo
      do ii=1,nsub
      kk=kk+1
      aoto2(ii,jj)=aoto(kk)
      enddo;enddo


      wlf=0.0 ; aotf =0.0
      zord = 0.0
      do j=1,nfuo
      do i = 1,nsub
       wlf(j) =  wlf(j)+  wlo2(i,j) * hkas(i,j)
       aotf(j)= aotf(j)+ aoto2(i,j) * hkas(i,j)
       zord = zord + sflx (i,j)*exp(-aoto2(i,j))
      enddo
      enddo

!-  WRITE OUT interpolated AOTs

!     do i=1,nwo
!     write(11) d1,d2,wlo(i),aoto(i),float(irec),log(wlo(i)),log(aoto(i)),float(ityp)
!     enddo
!     print'(A6,f10.3, 3i4)','FLUX= ',zord ,nsub,ityp,irh
      return
      end
!----------------------------------------------------------------
      subroutine aot_ext &
       (mxat,nwin,aotin,wlin,nwo,wlo,aoto,ityp,irh)
        USE FUINPUT, only:fi
      parameter(ne=24)

      real ,dimension(mxat)::aotin,wlin
!     real ,allocatable,dimension(-100:100) :: aoti,wlix
      real ,dimension(-100:100) :: aoti,wlix
      real ,dimension(nwo) :: aoto,wlo
      real ,dimension(ne) :: wlp,extp

      common /dalm_ext/    wld(24) ,datd(24,8,3)
      common /mineral_ext/ wlt(24) ,datt(24,4:8)
      common /opac_ext/ wlopac(24) ,datopac(24,8,9:18)
      common /mineral_ext_0404/ datt2(24,19:25)
! Wavelength MICRONS
! wlix,aoti = Monotonicly increasing

      idtl=-1
      if ( ityp >= 1 .and. ityp <=3 )then ! d'Almedia
       wlp  = wld
       extp = datd(1:24,irh,ityp)
       idtl=1
      elseif ( ityp >=4 .and. ityp <= 8 ) then ! Tegen&Lacis
       wlp  = wlt
       extp = datt(1:24,ityp)
       idtl=2
      elseif (  ityp >=19 .and. ityp <= 25 ) then ! Tegen&Lacis_0404
       wlp  = wlt
       extp = datt2(1:24,ityp)

       idtl=2
      elseif ( ityp >=9 .and. ityp <= 18) then ! OPAC
       wlp  = wlopac
        if (ityp==10 .or. ityp==12 .or.ityp==13 .or.ityp==18 ) then
         extp = datopac(1:24,irh,ityp)
        else
         extp = datopac(1:24,  1,ityp)
        endif



       idtl=3
      else
       fi%ierror = 18 !print*, ' Bad Aerosol type'
      endif


      if ( nwin == 1) then
!      nes=-3 ; nel=19  ! 1 chan @ 500nm
       if( wlin(1) <= 0.325 .or. wlin(1) >= 0.675) &
       fi%ierror = 19 ! print*, ' OUT OF ALLOWABLE ARANGE'
       nes= -(wlin(1)-0.325)/0.05
       if( idtl == 3) nes=nes-1  ! OPAC starts at 0.25um instead of 0.30
       nel = 22+nes

!      print*,'NES NEL',nes,nel,wlin(1)

      else
       nes=0
      if(idtl==1) nel=8    !  >= 2um long d'Almedia
      if(idtl==2) nel=11   !  >= 2um long Tegin&Lacis
      if(idtl==3) nel=7    !  >= 2um long OPAC

      endif

      nb = ne+1-nel
      nwi= nwin+nel-nes+1
      iend=nwin+nel
!     print*, icall,'in AOTEXT',nes,iend

!     if ( allocated (aoti) ) deallocate ( aoti )
!     allocate( aoti(nes:iend) )


!     if ( allocated (wlix) ) deallocate ( wlix)
!     allocate(  wlix(nes:iend) )
!     if(icall == 2) print*,

      wlix(1:nwin) =wlin(1:nwin)
      aoti(1:nwin)=aotin(1:nwin)


      LONGSIDE : do i=1,ne
      if( wlix(nwin) >= wlp(i) .and.  wlix(nwin) <= wlp(i+1)) then
!      print*,i  ,wlp(i),extp(i) ,i+1,wlp(i+1),extp(i+1)

      ext_norm1= rlnintrp( wlp(i),wlp(i+1), &
              extp(i),extp(i+1),wlix(nwin))

!      print*,dy,dx,dx1,yy,ext_norm1
       exit
       endif
      enddo LONGSIDE
!      wlix(nwin+1:nwi) = wlp(nb:ne)
!      aoti(nwin+1:nwi) = aoti(nwin)*(extp(nb:ne)/ext_norm1)
       wlix(nwin+1:nwin+nel) = wlp(nb:ne)
       aoti(nwin+1:nwin+nel) = aoti(nwin)*(extp(nb:ne)/ext_norm1)

!         print*,'CHECK',nwin+1,nwi,nb,ne,nwin+nel
!         print*,wlp(nb:ne)
!------
      if (nwin == 1 ) then
!     print*,1,wlix(1),aoti(1)
      SHORTSIDE: do i=1,ne
      if( wlix(1) >= wlp(i) .and.  wlix(1) <= wlp(i+1)) then
!      print*,i  ,wlp(i),extp(i) , i+1,wlp(i+1),extp(i+1)

      ext_norm0= rlnintrp( wlp(i),wlp(i+1), &
              extp(i),extp(i+1),wlix(1))

!      print*,dy,dx,dx1,yy,ext_norm0
       exit
       endif
      enddo SHORTSIDE
         nq = -nes+1
         wlix(nes:0) = wlp(1:nq)
        aoti(nes:0) = aoti(1)*(extp(1:nq)/ext_norm0)
      else
          wlix(0)  = 0.001
         aoti(0) = 1
      endif
!------------------------------------------------------------------
!     print'(a18,40f7.3)','Wavelength input= ',wlix(nes:iend)
!     print'(a18,40f7.3)','       AOT input= ',aoti(nes:iend)

!     do i=nes,iend
!     write(10) d1,d2,wlix(i),aoti(i),float(irec),log(wlix(i)),log(aoti(i)),float(ityp)
!     enddo

      call aotspline(nwi,aoti(nes:iend),wlix(nes:iend),nwo,wlo,aoto)

!     print'(a18,500f7.3)','Wavelength Out= ',wlo
!     print'(a18,500f7.3)','       AOT Out= ',aoto


      return
      end
!===================================================================
!===================================================================
      real function rlnintrp(x1,x2,y1,y2, x)
       dx= log(x2) - log(x1)
       dy= log(y2) - log(y1)
       dx1=log(x)  - log(x1)
       yy= (dy/dx) * dx1
       rlnintrp = exp(log(y1)+yy)
      return
      end
!====================================================================
      subroutine aotspline(nwi,aoti,wli,nwo,wlo,aoto)

      real ,dimension(nwi)  :: aoti,wli

      real ,dimension(nwi+1):: xa,ya,y2a
      real ,dimension(nwo)  :: aoto,wlo,aa

      data yp1,ypn /1.0E+32,1.0E+32/

      nwi2=nwi+1

      xa(2:nwi+1)=log(wli(1:nwi))
      ya(2:nwi+1)=log(aoti(1:nwi))
      xa(1)=log(1.0E-6) !; xa(nwi2)=log(1.0E+6) !TENSION
      ya(1)=0      !; ya(nwi2)= ya(nwi+1)!TENSION

      call spline(xa,ya,nwi2,yp1,ypn,y2a)

      do iwo = 1,nwo
       x=log(wlo(iwo))
       call splint(xa,ya,y2a,nwi2,x,y)
       aoto(iwo)=exp(y)
      enddo

      return
      end
!------------------------------------------------------
      real function alphav(aot1,aot2,wl1,wl2)
      ar= aot1/aot2
      wr= wl1/wl2
      alphav = - log(ar)/ log(wr)
      return
      end
!---------------------------------------------------------------
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1) &
             -x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
             u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      USE FUINPUT, only:fi
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) fi%ierror = 20 ! print*, 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h** &
      2)/6.
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE polint(xa,ya,n,x,y,dy)
        USE FUINPUT, only:fi
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) fi%ierror = 21 !print*, 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .
!=============================================================
!===========================================================
      block data opac_extinctions
      common /opac_ext/ wl(24) ,edat(24,8,9:18)
       data wl / &
      0.2500E+00,0.3000E+00,0.3500E+00,0.4000E+00,0.4500E+00,0.5000E+00, &
      0.5500E+00,0.6000E+00,0.6500E+00,0.7000E+00,0.7500E+00,0.8000E+00, &
      0.9000E+00,0.1000E+01,0.1250E+01,0.1500E+01,0.1750E+01,0.2000E+01, &
      0.2500E+01,0.3000E+01,0.3200E+01,0.3390E+01,0.3500E+01,0.3750E+01/
 !--- -------------------------------------------------------
 !9)  inso Insoluble
       data (( edat(i,j, 9),i=1,24), j=1,1 ) / &
      0.9477E+00,0.9572E+00,0.9667E+00,0.9748E+00,0.9839E+00,0.9916E+00, &
      0.1000E+01,0.1008E+01,0.1016E+01,0.1024E+01,0.1031E+01,0.1038E+01, &
      0.1052E+01,0.1064E+01,0.1093E+01,0.1105E+01,0.1088E+01,0.1012E+01, &
      0.7983E+00,0.6625E+00,0.7897E+00,0.8403E+00,0.8668E+00,0.8205E+00/
 !--- -------------------------------------------------------
 !10) waso      Water Soluble    (8 RH%)
       data (( edat(i,j,10),i=1,24), j=1,8 ) / &
      0.2438E+01,0.2095E+01,0.1793E+01,0.1539E+01,0.1326E+01,0.1148E+01, &
      0.1000E+01,0.8739E+00,0.7689E+00,0.6782E+00,0.6032E+00,0.5249E+00, &
      0.4251E+00,0.3497E+00,0.2191E+00,0.1520E+00,0.9186E-01,0.4897E-01, &
      0.3146E-01,0.2746E-01,0.1577E-01,0.1314E-01,0.1167E-01,0.9260E-02, &
      0.2373E+01,0.2044E+01,0.1758E+01,0.1516E+01,0.1314E+01,0.1144E+01, &
      0.1000E+01,0.8785E+00,0.7759E+00,0.6882E+00,0.6135E+00,0.5395E+00, &
      0.4380E+00,0.3605E+00,0.2271E+00,0.1560E+00,0.9755E-01,0.5793E-01, &
      0.3256E-01,0.1309E+00,0.5573E-01,0.2602E-01,0.1914E-01,0.1332E-01, &
      0.2319E+01,0.2006E+01,0.1733E+01,0.1501E+01,0.1306E+01,0.1140E+01, &
      0.1000E+01,0.8812E+00,0.7803E+00,0.6942E+00,0.6199E+00,0.5477E+00, &
      0.4458E+00,0.3676E+00,0.2326E+00,0.1597E+00,0.1014E+00,0.6249E-01, &
      0.3361E-01,0.1628E+00,0.6956E-01,0.3125E-01,0.2246E-01,0.1529E-01, &
      0.2263E+01,0.1967E+01,0.1708E+01,0.1485E+01,0.1298E+01,0.1137E+01, &
      0.1000E+01,0.8840E+00,0.7847E+00,0.7002E+00,0.6264E+00,0.5557E+00, &
      0.4536E+00,0.3748E+00,0.2386E+00,0.1639E+00,0.1054E+00,0.6684E-01, &
      0.3483E-01,0.1857E+00,0.8036E-01,0.3571E-01,0.2545E-01,0.1714E-01, &
      0.2143E+01,0.1884E+01,0.1653E+01,0.1451E+01,0.1279E+01,0.1130E+01, &
      0.1000E+01,0.8897E+00,0.7942E+00,0.7128E+00,0.6405E+00,0.5724E+00, &
      0.4705E+00,0.3910E+00,0.2521E+00,0.1741E+00,0.1146E+00,0.7589E-01, &
      0.3791E-01,0.2186E+00,0.9814E-01,0.4402E-01,0.3134E-01,0.2095E-01, &
      0.2000E+01,0.1785E+01,0.1588E+01,0.1410E+01,0.1257E+01,0.1120E+01, &
      0.1000E+01,0.8968E+00,0.8062E+00,0.7286E+00,0.6585E+00,0.5932E+00, &
      0.4922E+00,0.4122E+00,0.2707E+00,0.1887E+00,0.1272E+00,0.8739E-01, &
      0.4252E-01,0.2441E+00,0.1154E+00,0.5345E-01,0.3848E-01,0.2582E-01, &
      0.1809E+01,0.1650E+01,0.1497E+01,0.1353E+01,0.1224E+01,0.1106E+01, &
      0.1000E+01,0.9072E+00,0.8240E+00,0.7520E+00,0.6856E+00,0.6241E+00, &
      0.5252E+00,0.4454E+00,0.3007E+00,0.2134E+00,0.1481E+00,0.1055E+00, &
      0.5081E-01,0.2683E+00,0.1372E+00,0.6728E-01,0.4955E-01,0.3371E-01, &
      0.1686E+01,0.1563E+01,0.1437E+01,0.1314E+01,0.1202E+01,0.1097E+01, &
      0.1000E+01,0.9142E+00,0.8366E+00,0.7687E+00,0.7051E+00,0.6462E+00, &
      0.5496E+00,0.4702E+00,0.3238E+00,0.2331E+00,0.1648E+00,0.1196E+00, &
      0.5785E-01,0.2815E+00,0.1521E+00,0.7776E-01,0.5825E-01,0.4012E-01/
 !--- -------------------------------------------------------
 !11) soot Soot
       data (( edat(i,j,11),i=1,24), j=1,1 ) / &
      0.2447E+01,0.2188E+01,0.1837E+01,0.1555E+01,0.1331E+01,0.1153E+01, &
      0.1000E+01,0.8818E+00,0.7906E+00,0.7082E+00,0.6445E+00,0.5904E+00, &
      0.5087E+00,0.4453E+00,0.3412E+00,0.2767E+00,0.2367E+00,0.2055E+00, &
      0.1639E+00,0.1398E+00,0.1284E+00,0.1218E+00,0.1189E+00,0.1105E+00/
 !--- -------------------------------------------------------
 !12) ssam Sea Salt (Accumulation Mode)  (8 RH%)
       data (( edat(i,j,12),i=1,24), j=1,8 ) / &
      0.8801E+00,0.9114E+00,0.9354E+00,0.9580E+00,0.9733E+00,0.9887E+00, &
      0.1000E+01,0.1002E+01,0.1005E+01,0.1003E+01,0.9963E+00,0.9846E+00, &
      0.9618E+00,0.9232E+00,0.8315E+00,0.7230E+00,0.6190E+00,0.5382E+00, &
      0.3864E+00,0.4589E+00,0.3172E+00,0.2808E+00,0.2656E+00,0.2284E+00, &
      0.8883E+00,0.9080E+00,0.9321E+00,0.9508E+00,0.9715E+00,0.9866E+00, &
      0.1000E+01,0.1013E+01,0.1021E+01,0.1026E+01,0.1030E+01,0.1030E+01, &
      0.1023E+01,0.1009E+01,0.9498E+00,0.8710E+00,0.7792E+00,0.6931E+00, &
      0.4881E+00,0.6444E+00,0.6217E+00,0.5226E+00,0.4772E+00,0.4012E+00, &
      0.8884E+00,0.9093E+00,0.9300E+00,0.9493E+00,0.9672E+00,0.9853E+00, &
      0.1000E+01,0.1012E+01,0.1023E+01,0.1031E+01,0.1036E+01,0.1040E+01, &
      0.1041E+01,0.1032E+01,0.9874E+00,0.9212E+00,0.8368E+00,0.7534E+00, &
      0.5404E+00,0.6945E+00,0.6969E+00,0.5975E+00,0.5486E+00,0.4654E+00, &
      0.8932E+00,0.9113E+00,0.9305E+00,0.9508E+00,0.9671E+00,0.9850E+00, &
      0.1000E+01,0.1016E+01,0.1025E+01,0.1035E+01,0.1044E+01,0.1050E+01, &
      0.1054E+01,0.1052E+01,0.1020E+01,0.9637E+00,0.8874E+00,0.8075E+00, &
      0.5909E+00,0.7366E+00,0.7591E+00,0.6629E+00,0.6126E+00,0.5246E+00, &
      0.9051E+00,0.9182E+00,0.9355E+00,0.9528E+00,0.9682E+00,0.9843E+00, &
      0.1000E+01,0.1014E+01,0.1028E+01,0.1040E+01,0.1050E+01,0.1059E+01, &
      0.1072E+01,0.1077E+01,0.1070E+01,0.1036E+01,0.9774E+00,0.9089E+00, &
      0.6961E+00,0.8117E+00,0.8691E+00,0.7870E+00,0.7377E+00,0.6452E+00, &
      0.9158E+00,0.9311E+00,0.9419E+00,0.9586E+00,0.9717E+00,0.9857E+00, &
      0.1000E+01,0.1015E+01,0.1025E+01,0.1037E+01,0.1052E+01,0.1062E+01, &
      0.1080E+01,0.1094E+01,0.1109E+01,0.1098E+01,0.1064E+01,0.1013E+01, &
      0.8218E+00,0.8883E+00,0.9782E+00,0.9207E+00,0.8777E+00,0.7878E+00, &
      0.9361E+00,0.9476E+00,0.9567E+00,0.9668E+00,0.9761E+00,0.9905E+00, &
      0.1000E+01,0.1010E+01,0.1019E+01,0.1030E+01,0.1041E+01,0.1051E+01, &
      0.1075E+01,0.1090E+01,0.1126E+01,0.1144E+01,0.1141E+01,0.1122E+01, &
      0.9894E+00,0.9762E+00,0.1094E+01,0.1079E+01,0.1052E+01,0.9811E+00, &
      0.9487E+00,0.9578E+00,0.9655E+00,0.9729E+00,0.9804E+00,0.9933E+00, &
      0.1000E+01,0.1008E+01,0.1014E+01,0.1023E+01,0.1033E+01,0.1042E+01, &
      0.1064E+01,0.1079E+01,0.1119E+01,0.1148E+01,0.1164E+01,0.1166E+01, &
      0.1091E+01,0.1025E+01,0.1148E+01,0.1162E+01,0.1150E+01,0.1102E+01/
 !--- -------------------------------------------------------
 !13) sscm Sea Salt (Coarse Mode)   (8 RH%)
       data (( edat(i,j,13),i=1,24), j=1,8 ) / &
      0.9667E+00,0.9741E+00,0.9797E+00,0.9810E+00,0.9885E+00,0.9930E+00, &
      0.1000E+01,0.9993E+00,0.1010E+01,0.1015E+01,0.1013E+01,0.1015E+01, &
      0.1029E+01,0.1039E+01,0.1062E+01,0.1075E+01,0.1107E+01,0.1123E+01, &
      0.1170E+01,0.1175E+01,0.1208E+01,0.1222E+01,0.1228E+01,0.1241E+01, &
      0.9759E+00,0.9817E+00,0.9864E+00,0.9907E+00,0.9941E+00,0.9969E+00, &
      0.1000E+01,0.1001E+01,0.1008E+01,0.1011E+01,0.1014E+01,0.1016E+01, &
      0.1023E+01,0.1027E+01,0.1041E+01,0.1060E+01,0.1073E+01,0.1093E+01, &
      0.1134E+01,0.1094E+01,0.1129E+01,0.1156E+01,0.1167E+01,0.1188E+01, &
      0.9786E+00,0.9834E+00,0.9872E+00,0.9894E+00,0.9929E+00,0.9969E+00, &
      0.1000E+01,0.1003E+01,0.1006E+01,0.1010E+01,0.1013E+01,0.1014E+01, &
      0.1019E+01,0.1028E+01,0.1040E+01,0.1051E+01,0.1070E+01,0.1082E+01, &
      0.1121E+01,0.1083E+01,0.1116E+01,0.1140E+01,0.1151E+01,0.1172E+01, &
      0.9798E+00,0.9855E+00,0.9882E+00,0.9925E+00,0.9953E+00,0.9978E+00, &
      0.1000E+01,0.1005E+01,0.1007E+01,0.1010E+01,0.1015E+01,0.1014E+01, &
      0.1021E+01,0.1023E+01,0.1040E+01,0.1051E+01,0.1064E+01,0.1078E+01, &
      0.1113E+01,0.1079E+01,0.1109E+01,0.1130E+01,0.1141E+01,0.1159E+01, &
      0.9830E+00,0.9866E+00,0.9899E+00,0.9929E+00,0.9965E+00,0.9983E+00, &
      0.1000E+01,0.1003E+01,0.1006E+01,0.1009E+01,0.1012E+01,0.1014E+01, &
      0.1018E+01,0.1024E+01,0.1035E+01,0.1042E+01,0.1054E+01,0.1065E+01, &
      0.1097E+01,0.1069E+01,0.1095E+01,0.1111E+01,0.1120E+01,0.1136E+01, &
      0.9843E+00,0.9864E+00,0.9885E+00,0.9925E+00,0.9954E+00,0.9977E+00, &
      0.1000E+01,0.1002E+01,0.1005E+01,0.1006E+01,0.1010E+01,0.1010E+01, &
      0.1015E+01,0.1016E+01,0.1029E+01,0.1038E+01,0.1044E+01,0.1056E+01, &
      0.1078E+01,0.1060E+01,0.1080E+01,0.1093E+01,0.1100E+01,0.1112E+01, &
      0.9872E+00,0.9898E+00,0.9923E+00,0.9948E+00,0.9966E+00,0.9977E+00, &
      0.1000E+01,0.1002E+01,0.1004E+01,0.1002E+01,0.1007E+01,0.1009E+01, &
      0.1013E+01,0.1016E+01,0.1025E+01,0.1032E+01,0.1040E+01,0.1046E+01, &
      0.1060E+01,0.1052E+01,0.1066E+01,0.1075E+01,0.1079E+01,0.1087E+01, &
      0.9882E+00,0.9875E+00,0.9913E+00,0.9877E+00,0.9953E+00,0.9983E+00, &
      0.1000E+01,0.1001E+01,0.1003E+01,0.1004E+01,0.1006E+01,0.1007E+01, &
      0.1010E+01,0.1013E+01,0.1020E+01,0.1026E+01,0.1032E+01,0.1038E+01, &
      0.1049E+01,0.1044E+01,0.1056E+01,0.1063E+01,0.1066E+01,0.1071E+01/
 !--- -------------------------------------------------------
 !14) minm Mineral Dust (Nucleation Mode)
       data (( edat(i,j,14),i=1,24), j=1,1 ) / &
      0.1000E+01,0.9711E+00,0.9279E+00,0.8725E+00,0.8129E+00,0.7512E+00, &
      0.6916E+00,0.6347E+00,0.5817E+00,0.5327E+00,0.4879E+00,0.4470E+00, &
      0.3760E+00,0.3175E+00,0.2129E+00,0.1474E+00,0.1051E+00,0.7735E-01, &
      0.4556E-01,0.3771E-01,0.2699E-01,0.2199E-01,0.1993E-01,0.1462E-01/
 !--- -------------------------------------------------------
 !15) miam Mineral Dust (Accumulation Mode)
       data (( edat(i,j,15),i=1,24), j=1,1 ) / &
      0.9037E+00,0.9193E+00,0.9354E+00,0.9513E+00,0.9682E+00,0.9837E+00, &
      0.1000E+01,0.1015E+01,0.1031E+01,0.1045E+01,0.1056E+01,0.1069E+01, &
      0.1088E+01,0.1103E+01,0.1117E+01,0.1105E+01,0.1074E+01,0.1030E+01, &
      0.9115E+00,0.7974E+00,0.7466E+00,0.7089E+00,0.6875E+00,0.6277E+00/
 !--- -------------------------------------------------------
 !16) micm Mineral Dust (Coarse Mode)
       data (( edat(i,j,16),i=1,24), j=1,1 ) / &
      0.9742E+00,0.9790E+00,0.9836E+00,0.9878E+00,0.9922E+00,0.9959E+00, &
      0.1000E+01,0.1004E+01,0.1008E+01,0.1012E+01,0.1014E+01,0.1018E+01, &
      0.1026E+01,0.1032E+01,0.1048E+01,0.1065E+01,0.1082E+01,0.1096E+01, &
      0.1128E+01,0.1151E+01,0.1167E+01,0.1178E+01,0.1184E+01,0.1199E+01/
 !--- -------------------------------------------------------
 !17) mitr Mineral Dust (Transported Mode)
       data (( edat(i,j,17),i=1,24), j=1,1 ) / &
      0.9303E+00,0.9420E+00,0.9537E+00,0.9652E+00,0.9773E+00,0.9879E+00, &
      0.1000E+01,0.1010E+01,0.1021E+01,0.1032E+01,0.1042E+01,0.1052E+01, &
      0.1073E+01,0.1089E+01,0.1124E+01,0.1145E+01,0.1158E+01,0.1155E+01, &
      0.1123E+01,0.1056E+01,0.1031E+01,0.1010E+01,0.9970E+00,0.9592E+00/
 !--- -------------------------------------------------------
 !18) suso  Sulfate Droplets  (8 RH%)
       data (( edat(i,j,18),i=1,24), j=1,8 ) / &
      0.1589E+01,0.1522E+01,0.1418E+01,0.1303E+01,0.1190E+01,0.1092E+01, &
      0.1000E+01,0.9143E+00,0.8376E+00,0.7654E+00,0.6996E+00,0.6399E+00, &
      0.5385E+00,0.4523E+00,0.2965E+00,0.1989E+00,0.1377E+00,0.9764E-01, &
      0.4830E-01,0.9007E-01,0.1102E+00,0.1189E+00,0.1158E+00,0.9204E-01, &
      0.1348E+01,0.1324E+01,0.1275E+01,0.1210E+01,0.1140E+01,0.1070E+01, &
      0.1000E+01,0.9327E+00,0.8690E+00,0.8095E+00,0.7523E+00,0.6989E+00, &
      0.6053E+00,0.5239E+00,0.3671E+00,0.2634E+00,0.1894E+00,0.1396E+00, &
      0.6786E-01,0.2235E+00,0.1501E+00,0.1073E+00,0.9444E-01,0.7255E-01, &
      0.1270E+01,0.1260E+01,0.1228E+01,0.1178E+01,0.1121E+01,0.1062E+01, &
      0.1000E+01,0.9397E+00,0.8817E+00,0.8269E+00,0.7733E+00,0.7227E+00, &
      0.6326E+00,0.5533E+00,0.3966E+00,0.2903E+00,0.2118E+00,0.1581E+00, &
      0.7746E-01,0.2542E+00,0.1694E+00,0.1142E+00,0.9794E-01,0.7432E-01, &
      0.1213E+01,0.1214E+01,0.1193E+01,0.1154E+01,0.1107E+01,0.1055E+01, &
      0.1000E+01,0.9457E+00,0.8920E+00,0.8410E+00,0.7904E+00,0.7420E+00, &
      0.6550E+00,0.5774E+00,0.4212E+00,0.3129E+00,0.2309E+00,0.1740E+00, &
      0.8606E-01,0.2758E+00,0.1859E+00,0.1223E+00,0.1035E+00,0.7800E-01, &
      0.1130E+01,0.1144E+01,0.1137E+01,0.1114E+01,0.1082E+01,0.1044E+01, &
      0.1000E+01,0.9549E+00,0.9094E+00,0.8652E+00,0.8201E+00,0.7764E+00, &
      0.6957E+00,0.6218E+00,0.4674E+00,0.3562E+00,0.2682E+00,0.2055E+00, &
      0.1039E+00,0.3110E+00,0.2181E+00,0.1419E+00,0.1185E+00,0.8881E-01, &
      0.1047E+01,0.1070E+01,0.1078E+01,0.1072E+01,0.1055E+01,0.1030E+01, &
      0.1000E+01,0.9671E+00,0.9313E+00,0.8955E+00,0.8577E+00,0.8197E+00, &
      0.7476E+00,0.6795E+00,0.5295E+00,0.4157E+00,0.3210E+00,0.2510E+00, &
      0.1311E+00,0.3533E+00,0.2636E+00,0.1744E+00,0.1455E+00,0.1094E+00, &
      0.9649E+00,0.9946E+00,0.1015E+01,0.1023E+01,0.1023E+01,0.1014E+01, &
      0.1000E+01,0.9820E+00,0.9595E+00,0.9353E+00,0.9080E+00,0.8792E+00, &
      0.8213E+00,0.7633E+00,0.6247E+00,0.5105E+00,0.4084E+00,0.3287E+00, &
      0.1810E+00,0.4156E+00,0.3391E+00,0.2356E+00,0.1991E+00,0.1520E+00, &
      0.9232E+00,0.9524E+00,0.9764E+00,0.9920E+00,0.1001E+01,0.1004E+01, &
      0.1000E+01,0.9924E+00,0.9792E+00,0.9644E+00,0.9463E+00,0.9258E+00, &
      0.8810E+00,0.8331E+00,0.7091E+00,0.5986E+00,0.4932E+00,0.4068E+00, &
      0.2352E+00,0.4729E+00,0.4135E+00,0.3016E+00,0.2589E+00,0.2015E+00/
            end
!=============================================================
      block data tegen_lacis_ext

      common /mineral_ext/ wl(24) ,dat(24,4:8)

      data wl / &
         0.30,   0.35,   0.40,   0.45,   0.50,   0.55, &
         0.60,   0.65,   0.70,   0.80,   1.00,   1.25, &
         1.50,   2.00,   2.51,   2.61,   2.83,   2.96, &
         3.04,   3.26,   3.47,   3.69,   3.90,   4.11/

      data dat/ &
!  4
       8.78E-01,9.06E-01,9.36E-01,9.67E-01,1.00E+00,1.03E+00, &
       1.06E+00,1.08E+00,1.09E+00,1.09E+00,1.02E+00,8.70E-01, &
       7.04E-01,4.36E-01,2.08E-01,1.95E-01,1.91E-01,1.73E-01, &
       1.51E-01,1.35E-01,9.41E-02,7.95E-02,6.45E-02,5.36E-02, &
!  5
       9.41E-01,9.56E-01,9.70E-01,9.85E-01,1.00E+00,1.02E+00, &
       1.03E+00,1.05E+00,1.07E+00,1.10E+00,1.18E+00,1.25E+00, &
       1.27E+00,1.19E+00,8.64E-01,8.36E-01,7.68E-01,7.25E-01, &
       6.93E-01,6.92E-01,5.44E-01,4.96E-01,4.41E-01,3.93E-01, &
!  6
       9.65E-01,9.75E-01,9.83E-01,9.92E-01,1.00E+00,1.01E+00, &
       1.02E+00,1.02E+00,1.03E+00,1.05E+00,1.08E+00,1.12E+00, &
       1.17E+00,1.27E+00,1.34E+00,1.34E+00,1.29E+00,1.29E+00, &
       1.30E+00,1.33E+00,1.27E+00,1.24E+00,1.21E+00,1.17E+00, &
!  7
       9.78E-01,9.84E-01,9.89E-01,9.95E-01,1.00E+00,1.00E+00, &
       1.01E+00,1.01E+00,1.02E+00,1.03E+00,1.05E+00,1.07E+00, &
       1.09E+00,1.13E+00,1.19E+00,1.20E+00,1.21E+00,1.22E+00, &
       1.24E+00,1.25E+00,1.29E+00,1.31E+00,1.34E+00,1.36E+00, &
!  8
       9.86E-01,9.90E-01,9.93E-01,9.97E-01,1.00E+00,1.00E+00, &
       1.01E+00,1.01E+00,1.01E+00,1.02E+00,1.03E+00,1.04E+00, &
       1.05E+00,1.08E+00,1.10E+00,1.11E+00,1.11E+00,1.12E+00, &
       1.12E+00,1.13E+00,1.14E+00,1.15E+00,1.17E+00,1.18E+00  &
       /

      end
!============================================================
        block data tegen_lacis_ext_0404
        common /mineral_ext_0404/ dat(24,19:25)
        data dat/ &
 !  19  
         0.877,   0.905,   0.935,   0.967,   1.000,   1.031, &
         1.057,   1.077,   1.090,   1.092,   1.013,   0.837, &
         0.651,   0.384,   0.208,   0.194,   0.190,   0.172, &
         0.150,   0.134,   0.094,   0.079,   0.064,   0.053, &
 !  20  
         0.941,   0.953,   0.967,   0.982,   1.000,   1.016, &
         1.032,   1.047,   1.065,   1.103,   1.179,   1.249, &
         1.261,   1.131,   0.862,   0.834,   0.766,   0.723, &
         0.692,   0.690,   0.543,   0.495,   0.440,   0.392, &
 !  21  
         0.966,   0.975,   0.984,   0.992,   1.000,   1.012, &
         1.014,   1.021,   1.033,   1.047,   1.084,   1.127, &
         1.179,   1.280,   1.339,   1.341,   1.294,   1.288, &
         1.301,   1.334,   1.268,   1.245,   1.210,   1.168, &
 !  22  
         0.979,   0.985,   0.990,   0.996,   1.000,   1.006, &
         1.008,   1.014,   1.020,   1.031,   1.047,   1.078, &
         1.091,   1.135,   1.192,   1.201,   1.211,   1.223, &
         1.238,   1.253,   1.294,   1.315,   1.339,   1.361, &
 !  23  
         0.986,   0.990,   0.993,   0.997,   1.000,   1.003, &
         1.007,   1.009,   1.013,   1.019,   1.029,   1.042, &
         1.054,   1.078,   1.106,   1.105,   1.112,   1.118, &
         1.123,   1.131,   1.144,   1.155,   1.165,   1.176, &
 !  24  
         0.900,   0.931,   0.960,   0.983,   1.000,   1.010, &
         1.014,   1.012,   1.002,   0.969,   0.854,   0.682, &
         0.522,   0.304,   0.164,   0.153,   0.149,   0.135, &
         0.118,   0.106,   0.074,   0.062,   0.051,   0.042, &
 !  25  
         0.927,   0.944,   0.963,   0.981,   1.000,   1.020, &
         1.040,   1.059,   1.077,   1.105,   1.134,   1.125, &
         1.089,   1.012,   0.927,   0.919,   0.875,   0.862, &
         0.864,   0.885,   0.833,   0.821,   0.807,   0.790  &
       /

        end block data tegen_lacis_ext_0404
!=============================================================
      block data dalmedia_ext

      common /dalm_ext/ wl(24) ,dat(24,8,3)

      data wl / &
         0.30,   0.35,   0.40,   0.45,   0.50,   0.55, &
         0.60,   0.65,   0.70,   0.75,   0.80,   0.90, &
         1.00,   1.25,   1.50,   1.75,   2.00,   2.50, &
         3.00,   3.20,   3.39,   3.50,   3.78,   4.00/
      data dat/ &
!  1
       2.07E-04,2.08E-04,2.08E-04,2.07E-04,2.07E-04,2.08E-04, &
       2.09E-04,2.11E-04,2.12E-04,2.12E-04,2.12E-04,2.09E-04, &
       2.04E-04,1.90E-04,1.79E-04,1.71E-04,1.67E-04,1.63E-04, &
       1.72E-04,1.69E-04,1.62E-04,1.66E-04,1.67E-04,1.67E-04, &
       2.45E-04,2.46E-04,2.45E-04,2.44E-04,2.43E-04,2.43E-04, &
       2.44E-04,2.45E-04,2.46E-04,2.44E-04,2.43E-04,2.38E-04, &
       2.31E-04,2.14E-04,2.01E-04,1.91E-04,1.85E-04,1.79E-04, &
       1.90E-04,1.87E-04,1.83E-04,1.83E-04,1.84E-04,1.84E-04, &
       3.52E-04,3.50E-04,3.50E-04,3.51E-04,3.50E-04,3.47E-04, &
       3.47E-04,3.47E-04,3.49E-04,3.50E-04,3.51E-04,3.51E-04, &
       3.48E-04,3.28E-04,3.10E-04,2.92E-04,2.80E-04,2.63E-04, &
       2.85E-04,2.82E-04,2.71E-04,2.68E-04,2.65E-04,2.65E-04, &
       7.98E-04,7.93E-04,7.87E-04,7.86E-04,7.82E-04,7.80E-04, &
       7.83E-04,7.86E-04,7.86E-04,7.84E-04,7.82E-04,7.83E-04, &
       7.90E-04,8.10E-04,8.05E-04,7.79E-04,7.47E-04,6.70E-04, &
       7.23E-04,7.47E-04,7.11E-04,6.94E-04,6.68E-04,6.55E-04, &
       1.14E-03,1.12E-03,1.12E-03,1.11E-03,1.11E-03,1.11E-03, &
       1.11E-03,1.11E-03,1.11E-03,1.11E-03,1.11E-03,1.11E-03, &
       1.11E-03,1.13E-03,1.15E-03,1.14E-03,1.11E-03,9.98E-04, &
       1.05E-03,1.10E-03,1.06E-03,1.04E-03,9.94E-04,9.67E-04, &
       1.69E-03,1.67E-03,1.66E-03,1.64E-03,1.64E-03,1.63E-03, &
       1.63E-03,1.62E-03,1.62E-03,1.61E-03,1.62E-03,1.63E-03, &
       1.63E-03,1.62E-03,1.65E-03,1.68E-03,1.67E-03,1.54E-03, &
       1.56E-03,1.65E-03,1.64E-03,1.61E-03,1.53E-03,1.49E-03, &
       2.88E-03,2.87E-03,2.86E-03,2.83E-03,2.82E-03,2.80E-03, &
       2.78E-03,2.76E-03,2.74E-03,2.76E-03,2.75E-03,2.74E-03, &
       2.74E-03,2.76E-03,2.75E-03,2.78E-03,2.83E-03,2.79E-03, &
       2.69E-03,2.83E-03,2.89E-03,2.88E-03,2.81E-03,2.73E-03, &
       4.24E-03,4.27E-03,4.26E-03,4.26E-03,4.24E-03,4.21E-03, &
       4.18E-03,4.16E-03,4.14E-03,4.13E-03,4.11E-03,4.09E-03, &
       4.09E-03,4.06E-03,4.09E-03,4.08E-03,4.10E-03,4.22E-03, &
       4.02E-03,4.16E-03,4.26E-03,4.30E-03,4.31E-03,4.25E-03, &
!  2
       1.76E-05,1.57E-05,1.40E-05,1.25E-05,1.11E-05,9.95E-06, &
       8.91E-06,8.01E-06,7.21E-06,6.53E-06,5.78E-06,4.79E-06, &
       4.05E-06,2.66E-06,1.46E-06,1.02E-06,1.33E-06,4.35E-07, &
       3.49E-07,2.36E-07,2.02E-07,1.90E-07,1.54E-07,1.39E-07, &
       1.76E-05,1.57E-05,1.40E-05,1.25E-05,1.11E-05,9.95E-06, &
       8.91E-06,8.01E-06,7.21E-06,6.53E-06,5.79E-06,4.80E-06, &
       4.05E-06,2.66E-06,1.46E-06,1.02E-06,1.33E-06,4.37E-07, &
       3.50E-07,2.38E-07,2.03E-07,1.91E-07,1.55E-07,1.40E-07, &
       1.89E-05,1.69E-05,1.50E-05,1.34E-05,1.19E-05,1.07E-05, &
       9.57E-06,8.61E-06,7.75E-06,7.02E-06,6.23E-06,5.17E-06, &
       4.36E-06,2.87E-06,1.59E-06,1.12E-06,1.40E-06,4.73E-07, &
       5.67E-07,3.34E-07,2.44E-07,2.23E-07,1.77E-07,1.59E-07, &
       2.54E-05,2.27E-05,2.03E-05,1.81E-05,1.62E-05,1.45E-05, &
       1.31E-05,1.18E-05,1.06E-05,9.64E-06,8.61E-06,7.17E-06, &
       6.06E-06,4.01E-06,2.35E-06,1.65E-06,1.81E-06,6.73E-07, &
       1.72E-06,8.68E-07,4.78E-07,4.04E-07,3.02E-07,2.66E-07, &
       3.71E-05,3.35E-05,3.01E-05,2.71E-05,2.44E-05,2.20E-05, &
       1.99E-05,1.80E-05,1.64E-05,1.49E-05,1.35E-05,1.13E-05, &
       9.58E-06,6.42E-06,4.02E-06,2.84E-06,2.68E-06,1.11E-06, &
       4.01E-06,2.02E-06,1.01E-06,8.25E-07,5.90E-07,5.07E-07, &
       4.64E-05,4.22E-05,3.82E-05,3.46E-05,3.13E-05,2.83E-05, &
       2.57E-05,2.34E-05,2.14E-05,1.95E-05,1.77E-05,1.49E-05, &
       1.27E-05,8.61E-06,5.57E-06,3.96E-06,3.51E-06,1.52E-06, &
       5.97E-06,3.08E-06,1.54E-06,1.24E-06,8.75E-07,7.43E-07, &
       5.89E-05,5.40E-05,4.93E-05,4.50E-05,4.10E-05,3.73E-05, &
       3.41E-05,3.12E-05,2.86E-05,2.62E-05,2.39E-05,2.03E-05, &
       1.74E-05,1.20E-05,7.99E-06,5.74E-06,4.85E-06,2.21E-06, &
       8.84E-06,4.73E-06,2.40E-06,1.94E-06,1.38E-06,1.16E-06, &
       7.31E-05,6.77E-05,6.22E-05,5.72E-05,5.24E-05,4.80E-05, &
       4.41E-05,4.05E-05,3.73E-05,3.44E-05,3.15E-05,2.70E-05, &
       2.32E-05,1.62E-05,1.11E-05,8.06E-06,6.61E-06,3.11E-06, &
       1.23E-05,6.85E-06,3.56E-06,2.88E-06,2.05E-06,1.72E-06, &
!  3
       1.16E-05,1.03E-05,9.18E-06,8.17E-06,7.28E-06,6.49E-06, &
       5.81E-06,5.22E-06,4.70E-06,4.25E-06,3.77E-06,3.13E-06, &
       2.64E-06,1.74E-06,9.67E-07,6.83E-07,8.71E-07,3.00E-07, &
       2.42E-07,1.67E-07,1.44E-07,1.35E-07,1.11E-07,1.00E-07, &
       1.16E-05,1.03E-05,9.20E-06,8.18E-06,7.28E-06,6.50E-06, &
       5.82E-06,5.23E-06,4.70E-06,4.26E-06,3.78E-06,3.13E-06, &
       2.64E-06,1.74E-06,9.69E-07,6.85E-07,8.72E-07,3.02E-07, &
       2.44E-07,1.68E-07,1.45E-07,1.36E-07,1.12E-07,1.01E-07, &
       1.25E-05,1.11E-05,9.88E-06,8.79E-06,7.83E-06,6.98E-06, &
       6.25E-06,5.62E-06,5.06E-06,4.58E-06,4.07E-06,3.38E-06, &
       2.85E-06,1.88E-06,1.06E-06,7.49E-07,9.23E-07,3.28E-07, &
       3.90E-07,2.32E-07,1.72E-07,1.58E-07,1.28E-07,1.15E-07, &
       1.68E-05,1.49E-05,1.33E-05,1.19E-05,1.06E-05,9.48E-06, &
       8.51E-06,7.66E-06,6.92E-06,6.27E-06,5.61E-06,4.67E-06, &
       3.94E-06,2.61E-06,1.55E-06,1.09E-06,1.18E-06,4.55E-07, &
       1.14E-06,5.75E-07,3.19E-07,2.71E-07,2.05E-07,1.81E-07, &
       2.45E-05,2.20E-05,1.97E-05,1.77E-05,1.59E-05,1.43E-05, &
       1.29E-05,1.17E-05,1.06E-05,9.67E-06,8.71E-06,7.30E-06, &
       6.19E-06,4.15E-06,2.60E-06,1.85E-06,1.73E-06,7.27E-07, &
       2.63E-06,1.31E-06,6.55E-07,5.33E-07,3.82E-07,3.28E-07, &
       3.08E-05,2.78E-05,2.51E-05,2.26E-05,2.05E-05,1.85E-05, &
       1.67E-05,1.52E-05,1.39E-05,1.26E-05,1.15E-05,9.65E-06, &
       8.21E-06,5.55E-06,3.60E-06,2.56E-06,2.26E-06,9.87E-07, &
       3.93E-06,2.00E-06,9.86E-07,7.92E-07,5.59E-07,4.75E-07, &
       3.98E-05,3.62E-05,3.28E-05,2.98E-05,2.70E-05,2.45E-05, &
       2.23E-05,2.04E-05,1.86E-05,1.71E-05,1.55E-05,1.32E-05, &
       1.12E-05,7.70E-06,5.14E-06,3.68E-06,3.10E-06,1.40E-06, &
       5.86E-06,3.06E-06,1.52E-06,1.22E-06,8.51E-07,7.16E-07, &
       4.99E-05,4.58E-05,4.18E-05,3.81E-05,3.48E-05,3.17E-05, &
       2.90E-05,2.66E-05,2.44E-05,2.24E-05,2.06E-05,1.75E-05, &
       1.51E-05,1.04E-05,7.15E-06,5.15E-06,4.20E-06,1.95E-06, &
       8.20E-06,4.44E-06,2.24E-06,1.79E-06,1.25E-06,1.05E-06/

      end
!------------------------------------------------------------
      block data aerosol_convolve5

      common /aot_spect_5/  wlo(5,15) , hkas(5,15) ,sflx(5,15)
      data wlo / &
       0.1794,0.1878,0.1970,0.2073,0.2186, &
       0.2265,0.2301,0.2339,0.2378,0.2418, &
       0.2475,0.2551,0.2632,0.2717,0.2809, &
       0.2869,0.2894,0.2920,0.2946,0.2972, &
       0.3008,0.3053,0.3100,0.3149,0.3199, &
       0.3257,0.3323,0.3391,0.3462,0.3537, &
       0.3642,0.3783,0.3935,0.4100,0.4279, &
       0.4429,0.4539,0.4656,0.4778,0.4908, &
       0.5059,0.5233,0.5419,0.5620,0.5836, &
       0.6033,0.6206,0.6389,0.6583,0.6789, &
       0.7236,0.8026,0.9009,1.0267,1.1933, &
       1.3414,1.4358,1.5444,1.6708,1.8198, &
       1.9512,2.0513,2.1622,2.2857,2.4242, &
       2.5740,2.7360,2.9197,3.1299,3.3727, &
       3.5524,3.6430,3.7383,3.8388,3.9448 &
       /
      data hkas / &
       1.9732E-02,8.2461E-03,1.9433E-02,2.5239E-01,7.0020E-01, &
       1.7698E-01,1.9552E-01,1.8982E-01,2.0086E-01,2.3681E-01, &
       6.1144E-02,9.1518E-02,2.0913E-01,3.0821E-01,3.3000E-01, &
       1.1853E-01,1.6684E-01,2.5436E-01,2.2410E-01,2.3617E-01, &
       1.3174E-01,1.8637E-01,2.1818E-01,2.1883E-01,2.4487E-01, &
       1.6818E-01,1.9811E-01,1.7183E-01,2.1417E-01,2.4771E-01, &
       1.3629E-01,1.4266E-01,1.6518E-01,2.7083E-01,2.8504E-01, &
       1.6885E-01,1.9578E-01,2.0319E-01,2.1598E-01,2.1620E-01, &
       1.7878E-01,1.8362E-01,1.9993E-01,2.1178E-01,2.2589E-01, &
       1.9345E-01,1.9568E-01,2.0132E-01,2.0166E-01,2.0789E-01, &
       1.9354E-01,2.0025E-01,2.0400E-01,2.0302E-01,1.9919E-01, &
       2.1754E-01,2.1221E-01,2.0669E-01,1.9305E-01,1.7051E-01, &
       2.3672E-01,2.1723E-01,1.9925E-01,1.8215E-01,1.6466E-01, &
       2.4692E-01,2.2312E-01,1.9962E-01,1.7636E-01,1.5398E-01, &
       2.2023E-01,2.0999E-01,1.9844E-01,1.9036E-01,1.8098E-01 &
       /
      data sflx / &
       1.4056E-02,5.8741E-03,1.3843E-02,1.7979E-01,4.9878E-01, &
       1.5997E-01,1.7673E-01,1.7157E-01,1.8156E-01,2.1404E-01, &
       4.0480E-01,6.0590E-01,1.3845E+00,2.0405E+00,2.1848E+00, &
       7.6075E-01,1.0708E+00,1.6325E+00,1.4383E+00,1.5157E+00, &
       2.1362E+00,3.0220E+00,3.5377E+00,3.5483E+00,3.9706E+00, &
       5.6607E+00,6.6683E+00,5.7836E+00,7.2089E+00,8.3378E+00, &
       1.4768E+01,1.5458E+01,1.7899E+01,2.9348E+01,3.0887E+01, &
       1.9903E+01,2.3078E+01,2.3950E+01,2.5459E+01,2.5484E+01, &
       3.2273E+01,3.3147E+01,3.6093E+01,3.8232E+01,4.0778E+01, &
       2.9661E+01,3.0003E+01,3.0869E+01,3.0920E+01,3.1876E+01, &
       9.5588E+01,9.8903E+01,1.0076E+02,1.0027E+02,9.8377E+01, &
       3.3902E+01,3.3072E+01,3.2211E+01,3.0086E+01,2.6573E+01, &
       1.2208E+01,1.1203E+01,1.0275E+01,9.3939E+00,8.4915E+00, &
       7.0905E+00,6.4071E+00,5.7325E+00,5.0645E+00,4.4219E+00, &
       1.2380E+00,1.1805E+00,1.1155E+00,1.0701E+00,1.0173E+00 &
       /


      end
!----------------------------------------------------------------
      block data aerosol_convolve_25

      common /aot_spect_25/  wlo(25,15) , hkas(25,15) ,sflx(25,15)
      data wlo / &
       0.1762,0.1778,0.1794,0.1810,0.1826, &
       0.1843,0.1860,0.1878,0.1896,0.1914, &
       0.1932,0.1951,0.1970,0.1990,0.2010, &
       0.2030,0.2051,0.2073,0.2094,0.2116, &
       0.2139,0.2162,0.2186,0.2210,0.2235, &
       0.2251,0.2258,0.2265,0.2272,0.2279, &
       0.2287,0.2294,0.2301,0.2309,0.2316, &
       0.2324,0.2332,0.2339,0.2347,0.2355, &
       0.2362,0.2370,0.2378,0.2386,0.2394, &
       0.2402,0.2410,0.2418,0.2427,0.2435, &
       0.2446,0.2461,0.2475,0.2490,0.2505, &
       0.2520,0.2535,0.2551,0.2567,0.2583, &
       0.2599,0.2615,0.2632,0.2648,0.2665, &
       0.2682,0.2700,0.2717,0.2735,0.2753, &
       0.2772,0.2790,0.2809,0.2828,0.2847, &
       0.2860,0.2865,0.2869,0.2874,0.2879, &
       0.2884,0.2889,0.2894,0.2899,0.2904, &
       0.2910,0.2915,0.2920,0.2925,0.2930, &
       0.2935,0.2940,0.2946,0.2951,0.2956, &
       0.2961,0.2966,0.2972,0.2977,0.2982, &
       0.2991,0.3000,0.3009,0.3018,0.3027, &
       0.3036,0.3045,0.3054,0.3064,0.3073, &
       0.3082,0.3092,0.3101,0.3111,0.3120, &
       0.3130,0.3140,0.3150,0.3159,0.3169, &
       0.3179,0.3189,0.3199,0.3210,0.3220, &
       0.3232,0.3245,0.3258,0.3271,0.3284, &
       0.3297,0.3310,0.3323,0.3337,0.3350, &
       0.3364,0.3378,0.3392,0.3406,0.3420, &
       0.3434,0.3448,0.3463,0.3477,0.3492, &
       0.3507,0.3522,0.3537,0.3552,0.3567, &
       0.3590,0.3617,0.3643,0.3671,0.3698, &
       0.3726,0.3755,0.3784,0.3813,0.3843, &
       0.3874,0.3905,0.3936,0.3968,0.4000, &
       0.4033,0.4067,0.4101,0.4135,0.4170, &
       0.4206,0.4243,0.4280,0.4317,0.4356, &
       0.4387,0.4408,0.4429,0.4451,0.4473, &
       0.4495,0.4518,0.4540,0.4563,0.4586, &
       0.4609,0.4633,0.4656,0.4680,0.4705, &
       0.4729,0.4754,0.4779,0.4804,0.4830, &
       0.4855,0.4881,0.4908,0.4934,0.4961, &
       0.4996,0.5029,0.5062,0.5096,0.5130, &
       0.5165,0.5200,0.5236,0.5272,0.5309, &
       0.5346,0.5383,0.5422,0.5460,0.5500, &
       0.5540,0.5580,0.5621,0.5663,0.5705, &
       0.5748,0.5792,0.5836,0.5881,0.5927, &
       0.5969,0.6002,0.6035,0.6069,0.6103, &
       0.6137,0.6172,0.6207,0.6243,0.6279, &
       0.6316,0.6352,0.6390,0.6428,0.6466, &
       0.6504,0.6544,0.6583,0.6623,0.6664, &
       0.6705,0.6747,0.6789,0.6832,0.6875, &
       0.6962,0.7096,0.7236,0.7381,0.7532, &
       0.7690,0.7854,0.8026,0.8205,0.8392, &
       0.8588,0.8794,0.9009,0.9235,0.9473, &
       0.9724,0.9988,1.0267,1.0562,1.0874, &
       1.1206,1.1558,1.1933,1.2333,1.2762, &
       1.3070,1.3240,1.3414,1.3592,1.3776, &
       1.3965,1.4158,1.4358,1.4562,1.4773, &
       1.4990,1.5214,1.5444,1.5681,1.5926, &
       1.6179,1.6439,1.6708,1.6987,1.7274, &
       1.7572,1.7879,1.8198,1.8529,1.8871, &
       1.9139,1.9324,1.9512,1.9704,1.9900, &
       2.0101,2.0305,2.0513,2.0725,2.0942, &
       2.1164,2.1390,2.1622,2.1858,2.2099, &
       2.2346,2.2599,2.2857,2.3121,2.3392, &
       2.3669,2.3952,2.4242,2.4540,2.4845, &
       2.5145,2.5439,2.5740,2.6048,2.6364, &
       2.6688,2.7020,2.7360,2.7709,2.8066, &
       2.8433,2.8810,2.9197,2.9595,3.0003, &
       3.0423,3.0855,3.1299,3.1756,3.2227, &
       3.2712,3.3212,3.3727,3.4258,3.4807, &
       3.5174,3.5348,3.5524,3.5702,3.5881, &
       3.6062,3.6245,3.6430,3.6617,3.6805, &
       3.6996,3.7189,3.7383,3.7580,3.7779, &
       3.7979,3.8183,3.8388,3.8595,3.8805, &
       3.9017,3.9231,3.9448,3.9667,3.9888 &
       /
      data hkas / &
       5.2799E-03,5.0829E-03,4.6490E-03,4.0241E-03,3.4829E-03, &
       2.8976E-03,2.3480E-03,1.8516E-03,1.3749E-03,9.3841E-04, &
       5.9685E-04,3.2885E-04,1.3047E-04,2.2220E-05,5.5991E-03, &
       2.7673E-02,3.1899E-02,3.7975E-02,4.7693E-02,8.3620E-02, &
       1.1867E-01,1.4080E-01,1.2000E-01,1.8402E-01,1.6904E-01, &
       5.0121E-02,4.9452E-02,4.6279E-02,4.0956E-02,3.7858E-02, &
       2.6053E-02,3.0795E-02,4.6226E-02,3.7153E-02,3.5560E-02, &
       4.7180E-02,3.7850E-02,4.0573E-02,4.5254E-02,3.5590E-02, &
       3.0742E-02,4.1064E-02,4.1940E-02,4.0586E-02,4.7516E-02, &
       3.2378E-02,4.4837E-02,3.5456E-02,3.5664E-02,4.2917E-02, &
       1.6466E-02,1.5248E-02,1.2129E-02,1.2876E-02,1.0915E-02, &
       1.5288E-02,1.2174E-02,1.1518E-02,1.5566E-02,2.3291E-02, &
       3.5784E-02,2.7197E-02,2.5602E-02,2.8189E-02,6.8271E-02, &
       7.4253E-02,7.6509E-02,7.4379E-02,7.6454E-02,6.3299E-02, &
       4.4999E-02,7.9458E-02,5.3948E-02,3.1778E-02,9.4411E-02, &
       2.9705E-02,2.9759E-02,2.6074E-02,1.4967E-02,7.4578E-03, &
       2.3107E-02,3.1405E-02,3.1772E-02,3.8404E-02,2.7981E-02, &
       2.5950E-02,3.7908E-02,4.2120E-02,5.6514E-02,6.0159E-02, &
       5.4489E-02,5.7136E-02,5.3536E-02,4.5756E-02,5.3907E-02, &
       5.1124E-02,4.7270E-02,5.0278E-02,5.2757E-02,5.0463E-02, &
       3.0474E-02,3.1978E-02,2.5582E-02,3.1001E-02,2.6341E-02, &
       2.9467E-02,2.5155E-02,3.6947E-02,4.1472E-02,3.9135E-02, &
       3.7515E-02,3.9334E-02,4.3622E-02,4.5644E-02,3.7261E-02, &
       3.9516E-02,6.1067E-02,4.4004E-02,4.7490E-02,4.7719E-02, &
       5.4569E-02,3.5377E-02,5.9282E-02,4.1212E-02,4.8834E-02, &
       3.3110E-02,2.7957E-02,2.4541E-02,3.1640E-02,3.8452E-02, &
       4.3699E-02,3.7692E-02,4.5887E-02,3.6735E-02,4.2471E-02, &
       4.0508E-02,4.0549E-02,3.1263E-02,3.2911E-02,3.6471E-02, &
       4.0021E-02,3.6207E-02,4.0685E-02,3.9095E-02,4.8464E-02, &
       4.7468E-02,4.8137E-02,5.3969E-02,4.7248E-02,5.4817E-02, &
       2.4834E-02,2.1951E-02,2.6050E-02,2.7710E-02,3.4553E-02, &
       3.1715E-02,2.9451E-02,2.7336E-02,3.6347E-02,3.1034E-02, &
       2.3427E-02,2.7237E-02,3.7576E-02,2.6206E-02,3.2126E-02, &
       5.2125E-02,5.2807E-02,5.1698E-02,5.6664E-02,5.8765E-02, &
       6.0530E-02,5.9962E-02,5.9826E-02,5.6824E-02,5.3245E-02, &
       2.9942E-02,3.4314E-02,3.1576E-02,3.4169E-02,3.6263E-02, &
       3.5425E-02,3.8424E-02,4.0153E-02,3.9218E-02,3.9656E-02, &
       4.1712E-02,4.0169E-02,4.1568E-02,4.1185E-02,4.1190E-02, &
       4.2675E-02,4.1942E-02,4.2627E-02,4.4241E-02,4.5120E-02, &
       4.5103E-02,4.4734E-02,3.9963E-02,4.4448E-02,4.4182E-02, &
       3.6238E-02,3.5264E-02,3.4319E-02,3.6506E-02,3.6858E-02, &
       3.7418E-02,3.5769E-02,3.4340E-02,3.7866E-02,3.7059E-02, &
       4.0566E-02,3.8579E-02,4.0826E-02,3.9780E-02,4.1088E-02, &
       4.1249E-02,4.2418E-02,4.2017E-02,4.2620E-02,4.3139E-02, &
       4.3749E-02,4.5394E-02,4.5299E-02,4.6414E-02,4.5227E-02, &
       3.7723E-02,3.8411E-02,3.8671E-02,3.8195E-02,3.9511E-02, &
       3.9418E-02,3.9055E-02,3.8240E-02,3.9599E-02,3.9985E-02, &
       3.9300E-02,4.0631E-02,3.9929E-02,4.0581E-02,4.0475E-02, &
       4.0624E-02,4.0457E-02,4.1359E-02,3.7988E-02,4.1377E-02, &
       4.1472E-02,4.1668E-02,4.1660E-02,4.1838E-02,4.1834E-02, &
       3.8155E-02,3.8333E-02,3.8502E-02,3.9057E-02,3.9013E-02, &
       3.9461E-02,3.9758E-02,4.0111E-02,4.0243E-02,4.0392E-02, &
       3.9720E-02,4.0458E-02,4.1065E-02,4.1312E-02,4.1326E-02, &
       4.1240E-02,4.1154E-02,4.0773E-02,4.0305E-02,4.0042E-02, &
       4.0018E-02,4.0038E-02,3.9956E-02,3.9901E-02,3.9668E-02, &
       4.4043E-02,4.3956E-02,4.3412E-02,4.2920E-02,4.2591E-02, &
       4.2487E-02,4.2451E-02,4.2278E-02,4.1980E-02,4.1849E-02, &
       4.1556E-02,4.1621E-02,4.1312E-02,4.0978E-02,4.0676E-02, &
       4.0101E-02,3.9593E-02,3.8711E-02,3.8058E-02,3.7376E-02, &
       3.6255E-02,3.5384E-02,3.4555E-02,3.3467E-02,3.2390E-02, &
       4.9085E-02,4.8014E-02,4.7452E-02,4.6084E-02,4.6143E-02, &
       4.5126E-02,4.4490E-02,4.3771E-02,4.2825E-02,4.1722E-02, &
       4.0817E-02,4.0346E-02,3.9877E-02,3.8934E-02,3.8461E-02, &
       3.8034E-02,3.7265E-02,3.6323E-02,3.5726E-02,3.5049E-02, &
       3.4210E-02,3.3481E-02,3.2891E-02,3.2285E-02,3.1588E-02, &
       5.1217E-02,5.0595E-02,4.9284E-02,4.8233E-02,4.7849E-02, &
       4.6488E-02,4.5516E-02,4.4685E-02,4.3738E-02,4.2612E-02, &
       4.1889E-02,4.0913E-02,3.9819E-02,3.8998E-02,3.8111E-02, &
       3.7139E-02,3.6048E-02,3.5235E-02,3.4379E-02,3.3489E-02, &
       3.2579E-02,3.1537E-02,3.0736E-02,2.9824E-02,2.9085E-02, &
       4.5112E-02,4.4731E-02,4.4278E-02,4.3820E-02,4.3350E-02, &
       4.2925E-02,4.2577E-02,4.2229E-02,4.1713E-02,4.1188E-02, &
       4.0646E-02,4.0021E-02,3.9857E-02,3.8890E-02,3.8999E-02, &
       3.8626E-02,3.8315E-02,3.7974E-02,3.7671E-02,3.7171E-02, &
       3.6752E-02,3.6471E-02,3.6046E-02,3.5548E-02,3.5090E-02 &
       /
      data sflx/ &
       3.3009E-03,3.1778E-03,2.9065E-03,2.5158E-03,2.1775E-03, &
       1.8116E-03,1.4679E-03,1.1576E-03,8.5959E-04,5.8668E-04, &
       3.7314E-04,2.0560E-04,8.1568E-05,1.3892E-05,3.5005E-03, &
       1.7301E-02,1.9943E-02,2.3742E-02,2.9817E-02,5.2278E-02, &
       7.4193E-02,8.8026E-02,7.5023E-02,1.1505E-01,1.0568E-01, &
       4.4282E-02,4.3691E-02,4.0888E-02,3.6184E-02,3.3448E-02, &
       2.3018E-02,2.7208E-02,4.0841E-02,3.2825E-02,3.1417E-02, &
       4.1683E-02,3.3441E-02,3.5846E-02,3.9982E-02,3.1444E-02, &
       2.7161E-02,3.6280E-02,3.7054E-02,3.5858E-02,4.1981E-02, &
       2.8606E-02,3.9614E-02,3.1325E-02,3.1509E-02,3.7918E-02, &
       1.0101E-01,9.3537E-02,7.4402E-02,7.8986E-02,6.6956E-02, &
       9.3786E-02,7.4679E-02,7.0655E-02,9.5489E-02,1.4288E-01, &
       2.1951E-01,1.6683E-01,1.5705E-01,1.7292E-01,4.1880E-01, &
       4.5550E-01,4.6934E-01,4.5627E-01,4.6900E-01,3.8830E-01, &
       2.7604E-01,4.8743E-01,3.3094E-01,1.9494E-01,5.7915E-01, &
       1.7599E-01,1.7631E-01,1.5448E-01,8.8671E-02,4.4184E-02, &
       1.3690E-01,1.8606E-01,1.8824E-01,2.2753E-01,1.6577E-01, &
       1.5374E-01,2.2459E-01,2.4955E-01,3.3482E-01,3.5642E-01, &
       3.2283E-01,3.3851E-01,3.1718E-01,2.7109E-01,3.1938E-01, &
       3.0289E-01,2.8005E-01,2.9787E-01,3.1256E-01,2.9897E-01, &
       4.7450E-01,4.9792E-01,3.9834E-01,4.8271E-01,4.1015E-01, &
       4.5883E-01,3.9168E-01,5.7529E-01,6.4574E-01,6.0936E-01, &
       5.8414E-01,6.1246E-01,6.7922E-01,7.1071E-01,5.8018E-01, &
       6.1529E-01,9.5085E-01,6.8518E-01,7.3945E-01,7.4301E-01, &
       8.4968E-01,5.5085E-01,9.2306E-01,6.4171E-01,7.6037E-01, &
       1.0836E+00,9.1500E-01,8.0321E-01,1.0356E+00,1.2585E+00, &
       1.4302E+00,1.2336E+00,1.5018E+00,1.2023E+00,1.3900E+00, &
       1.3258E+00,1.3271E+00,1.0232E+00,1.0772E+00,1.1937E+00, &
       1.3099E+00,1.1850E+00,1.3316E+00,1.2795E+00,1.5862E+00, &
       1.5536E+00,1.5755E+00,1.7664E+00,1.5464E+00,1.7941E+00, &
       2.6150E+00,2.3115E+00,2.7431E+00,2.9179E+00,3.6384E+00, &
       3.3396E+00,3.1012E+00,2.8785E+00,3.8273E+00,3.2679E+00, &
       2.4668E+00,2.8681E+00,3.9567E+00,2.7595E+00,3.3828E+00, &
       5.4888E+00,5.5606E+00,5.4437E+00,5.9667E+00,6.1879E+00, &
       6.3738E+00,6.3140E+00,6.2996E+00,5.9836E+00,5.6067E+00, &
       3.4766E+00,3.9842E+00,3.6664E+00,3.9674E+00,4.2105E+00, &
       4.1133E+00,4.4615E+00,4.6622E+00,4.5536E+00,4.6045E+00, &
       4.8432E+00,4.6640E+00,4.8265E+00,4.7820E+00,4.7826E+00, &
       4.9550E+00,4.8700E+00,4.9494E+00,5.1368E+00,5.2390E+00, &
       5.2370E+00,5.1941E+00,4.6402E+00,5.1609E+00,5.1301E+00, &
       6.4948E+00,6.3201E+00,6.1509E+00,6.5429E+00,6.6059E+00, &
       6.7062E+00,6.4108E+00,6.1545E+00,6.7865E+00,6.6419E+00, &
       7.2704E+00,6.9142E+00,7.3170E+00,7.1295E+00,7.3641E+00, &
       7.3928E+00,7.6024E+00,7.5305E+00,7.6385E+00,7.7317E+00, &
       7.8410E+00,8.1357E+00,8.1187E+00,8.3186E+00,8.1057E+00, &
       5.7948E+00,5.9004E+00,5.9404E+00,5.8673E+00,6.0694E+00, &
       6.0551E+00,5.9994E+00,5.8742E+00,6.0830E+00,6.1422E+00, &
       6.0369E+00,6.2414E+00,6.1336E+00,6.2338E+00,6.2175E+00, &
       6.2403E+00,6.2148E+00,6.3533E+00,5.8355E+00,6.3561E+00, &
       6.3706E+00,6.4007E+00,6.3995E+00,6.4268E+00,6.4263E+00, &
       1.8893E+01,1.8981E+01,1.9065E+01,1.9339E+01,1.9317E+01, &
       1.9539E+01,1.9686E+01,1.9861E+01,1.9927E+01,2.0000E+01, &
       1.9668E+01,2.0033E+01,2.0333E+01,2.0456E+01,2.0463E+01, &
       2.0420E+01,2.0378E+01,2.0189E+01,1.9958E+01,1.9827E+01, &
       1.9815E+01,1.9825E+01,1.9784E+01,1.9757E+01,1.9642E+01, &
       6.9768E+00,6.9631E+00,6.8768E+00,6.7990E+00,6.7469E+00, &
       6.7303E+00,6.7247E+00,6.6973E+00,6.6500E+00,6.6293E+00, &
       6.5828E+00,6.5933E+00,6.5442E+00,6.4914E+00,6.4434E+00, &
       6.3523E+00,6.2719E+00,6.1322E+00,6.0288E+00,5.9207E+00, &
       5.7431E+00,5.6052E+00,5.4739E+00,5.3015E+00,5.1310E+00, &
       2.6183E+00,2.5612E+00,2.5312E+00,2.4582E+00,2.4613E+00, &
       2.4071E+00,2.3732E+00,2.3348E+00,2.2844E+00,2.2255E+00, &
       2.1772E+00,2.1521E+00,2.1271E+00,2.0768E+00,2.0516E+00, &
       2.0288E+00,1.9878E+00,1.9375E+00,1.9057E+00,1.8696E+00, &
       1.8248E+00,1.7859E+00,1.7544E+00,1.7221E+00,1.6850E+00, &
       1.5186E+00,1.5002E+00,1.4613E+00,1.4301E+00,1.4187E+00, &
       1.3784E+00,1.3496E+00,1.3249E+00,1.2968E+00,1.2635E+00, &
       1.2420E+00,1.2131E+00,1.1806E+00,1.1563E+00,1.1300E+00, &
       1.1012E+00,1.0688E+00,1.0447E+00,1.0193E+00,9.9296E-01, &
       9.6597E-01,9.3507E-01,9.1134E-01,8.8430E-01,8.6237E-01, &
       2.6980E-01,2.6752E-01,2.6481E-01,2.6207E-01,2.5926E-01, &
       2.5672E-01,2.5464E-01,2.5256E-01,2.4947E-01,2.4633E-01, &
       2.4309E-01,2.3935E-01,2.3837E-01,2.3259E-01,2.3324E-01, &
       2.3101E-01,2.2915E-01,2.2711E-01,2.2530E-01,2.2231E-01, &
       2.1980E-01,2.1812E-01,2.1558E-01,2.1260E-01,2.0986E-01 &
       /

      end

