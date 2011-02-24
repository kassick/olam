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
subroutine interp_hs()

use mem_hs
use mem_basic
use mem_grid
use mem_ijtabs
use consts_coms

implicit none

integer :: ilat,ilon,k,iw,i,j,kup,im1,im2,im3,iwave
integer :: npts
real :: wt1,wtr,wt2,wt3,xep,yep,zep,alon,alat,fldval,store1,bandmass

real, save :: atime=0.

!!!!!!!!!! FFT & CONTOUR stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer, parameter :: n=144,np1=145,np2=146,lot=2880,ntrigs=217,nwork=417600  
! lot = 90 * nzhs
! ntrigs = 3*n/2+1
! nwork = (np1)*lot
! Initialize all fields on HS grid to zero prior to summing
real, dimension(np2,1) :: a
real, dimension(ntrigs) :: trigs
real, dimension(nwork) :: work
integer, dimension(13) :: ifax

integer :: isign,jump,inc,notavail

integer, save :: kzdt=90, mzdt=90, nset=0, nhgh=0, ndsh=0
real, save    :: flow=0., fhgh=0., finc=0.

isign = -1
jump = np2
inc = 1
call fftfax(n,ifax,trigs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Zero out both 3D lat/lon grid arrays plus count array prior to interpolation 
! from triangle mesh

do ilat = 1,90
   do ilon = 1,144
      count_hs(ilon,ilat) = 0.
      do k = 1,nzhs
         uzonal_hs(ilon,ilat,k) = 0.
         theta_hs (ilon,ilat,k) = 0.
         rho_hs   (ilon,ilat,k) = 0.
         temp_hs  (ilon,ilat,k) = 0.
      enddo
   enddo
enddo

! Loop over all W columns of triangle mesh

do iw = 1,mwa

   npts = int (5. / cos(glatw(iw)*pio180))

   im1 = itab_w(iw)%im1
   im2 = itab_w(iw)%im2
   im3 = itab_w(iw)%im3

! Assume that j spans from im1 vertex to im2-im3 side

   do j = 1,npts
      do i = 1,npts
         do kup = 1,2
            if (kup == 1) then
            	
! Upper point weights

               wt1 = real(3*(npts-j)+1) / real(3*npts)
               wtr = real(2*i-1) / real(2*npts)

            else

               if (i == npts) cycle

! Lower point weights

               wt1 = real(3*(npts-j)+2) / real(3*npts)
               wtr = real(2*i) / real(2*npts)
               
            endif

            wt2 = (1. - wt1) * wtr
            wt3 = (1. - wt1) * (1. - wtr)

! W point (xep,yep,zep) coordinates
         
            xep = wt1 * xem(im1)  &
                + wt2 * xem(im2)  &
                + wt3 * xem(im3)

            yep = wt1 * yem(im1)  &
                + wt2 * yem(im2)  &
                + wt3 * yem(im3)

            zep = wt1 * zem(im1)  &
                + wt2 * zem(im2)  &
                + wt3 * zem(im3)

! Transform from Earth to EC (lat/lon) coordinates

            call e_ec(xep,yep,zep,alon,alat)  ! alon/alat in degrees

! Find lat/lon grid cell index (cells are 2.0 deg lat and 2.5 deg lon)

            ilon = int (.4 * (alon + 180.)) + 1
            ilat = int (.5 * (alat + 90.)) + 1
         
! Add contribution of both value and number to current ilon/ilat cell

            count_hs(ilon,ilat) = count_hs(ilon,ilat) + 1.
         
            do k = 1,nzhs

               call oplot_lib(1,k+1,iw,'VALUE','ZONAL_WIND',fldval,notavail)    
               uzonal_hs(ilon,ilat,k) = uzonal_hs(ilon,ilat,k) + fldval
               call oplot_lib(1,k+1,iw,'VALUE','MERID_WIND',fldval,notavail)    
               umerid_hs(ilon,ilat,k) = umerid_hs(ilon,ilat,k) + fldval
               theta_hs(ilon,ilat,k) = theta_hs(ilon,ilat,k) + theta(k+1,iw)
               rho_hs(ilon,ilat,k) = rho_hs(ilon,ilat,k) + rho(k+1,iw)
               call oplot_lib(1,k+1,iw,'VALUE','AIRTEMPK',fldval,notavail)    
               temp_hs(ilon,ilat,k) = temp_hs(ilon,ilat,k) + fldval
                     
            enddo   ! end k loop

         enddo  ! end kup loop
      enddo     ! end i loop
   enddo        ! end j loop

enddo           ! end iw loop

! Normalize both fields on HS grid

do k = 1,nzhs
   do ilat = 1,90
      do ilon = 1,144
         uzonal_hs(ilon,ilat,k) = uzonal_hs(ilon,ilat,k) / count_hs(ilon,ilat)
         umerid_hs(ilon,ilat,k) = umerid_hs(ilon,ilat,k) / count_hs(ilon,ilat)
         theta_hs (ilon,ilat,k) = theta_hs (ilon,ilat,k) / count_hs(ilon,ilat)
         rho_hs   (ilon,ilat,k) = rho_hs   (ilon,ilat,k) / count_hs(ilon,ilat)
         temp_hs  (ilon,ilat,k) = temp_hs  (ilon,ilat,k) / count_hs(ilon,ilat)
      enddo

! Set extra 2 points with wrap-around values

      uzonal_hs(145,ilat,k) = uzonal_hs(1,ilat,k)
      umerid_hs(145,ilat,k) = umerid_hs(1,ilat,k)
      theta_hs (145,ilat,k) = theta_hs (1,ilat,k)
      rho_hs   (145,ilat,k) = rho_hs   (1,ilat,k)
      temp_hs  (145,ilat,k) = temp_hs  (1,ilat,k)

      uzonal_hs(146,ilat,k) = uzonal_hs(2,ilat,k)
      umerid_hs(146,ilat,k) = umerid_hs(2,ilat,k)
      theta_hs (146,ilat,k) = theta_hs (2,ilat,k)
      rho_hs   (146,ilat,k) = rho_hs   (2,ilat,k)
      temp_hs  (146,ilat,k) = temp_hs  (2,ilat,k)

   enddo
enddo

! Compute and plot averages

do k = 1,nzhs
   do ilat = 1,90
      uzonal_avg_hs(ilat,k) = 0.
      umerid_avg_hs(ilat,k) = 0.
      theta_avg_hs (ilat,k) = 0.
      rho_avg_hs   (ilat,k) = 0.
      temp_avg_hs  (ilat,k) = 0.
      do ilon = 1,144
         uzonal_avg_hs(ilat,k) = uzonal_avg_hs(ilat,k)  &
                               + uzonal_hs(ilon,ilat,k) / 144.
         umerid_avg_hs(ilat,k) = umerid_avg_hs(ilat,k)  &
                               + umerid_hs(ilon,ilat,k) / 144.
         theta_avg_hs (ilat,k) = theta_avg_hs (ilat,k)  &
                               + theta_hs (ilon,ilat,k) / 144.
         rho_avg_hs   (ilat,k) = rho_avg_hs (ilat,k)  &
                               + rho_hs   (ilon,ilat,k) / 144.
         temp_avg_hs (ilat,k)  = temp_avg_hs  (ilat,k)  &
                               + temp_hs  (ilon,ilat,k) / 144.
      enddo
!!!!!!!!!!!!!!!!! TIME AVERAGES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      uzonal_tavg_hs(ilat,k) = uzonal_tavg_hs(ilat,k) * atime / (atime + 1.)  &
                             + uzonal_avg_hs(ilat,k) / (atime + 1.)  
      umerid_tavg_hs(ilat,k) = umerid_tavg_hs(ilat,k) * atime / (atime + 1.)  &
                             + umerid_avg_hs(ilat,k) / (atime + 1.)  
      theta_tavg_hs (ilat,k) = theta_tavg_hs (ilat,k) * atime / (atime + 1.)  &
                             + theta_avg_hs (ilat,k) / (atime + 1.)  
      temp_tavg_hs  (ilat,k) = temp_tavg_hs  (ilat,k) * atime / (atime + 1.)  &
                             + temp_avg_hs  (ilat,k) / (atime + 1.)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   enddo
enddo

! Plot latitude-height dependence of amplitude
   

call plotback()
call cpcnrc(uzonal_avg_hs, kzdt, mzdt, nzhs  &
   , flow, fhgh, 4.0, nset, -1, -682)
call frame

call plotback()
call cpcnrc(uzonal_tavg_hs, kzdt, mzdt, nzhs  &
   , flow, fhgh, 4.0, nset, -1, -682)
call frame

call plotback()
call cpcnrc(umerid_avg_hs, kzdt, mzdt, nzhs  &
   , flow, fhgh, 0.5, nset, -1, -682)
call frame

call plotback()
call cpcnrc(umerid_tavg_hs, kzdt, mzdt, nzhs  &
   , flow, fhgh, 0.5, nset, -1, -682)
call frame

call plotback()
call cpcnrc(theta_avg_hs,  kzdt, mzdt, nzhs  &
   , flow, fhgh, 5.0, nset, -1, ndsh)
call frame
   
call plotback()
call cpcnrc(theta_tavg_hs,  kzdt, mzdt, nzhs  &
   , flow, fhgh, 5.0, nset, -1, ndsh)
call frame
   
call plotback()
call cpcnrc(temp_avg_hs,  kzdt, mzdt, nzhs  &
   , flow, fhgh, 5.0, nset, -1, ndsh)
call frame
   
call plotback()
call cpcnrc(temp_tavg_hs,  kzdt, mzdt, nzhs  &
   , flow, fhgh, 5.0, nset, -1, ndsh)
call frame
   
! Compute eddy variance of temperature

do k = 1,nzhs
   do ilat = 1,90
      store1 = 0.
      do ilon = 1,144
         store1 = store1 + (temp_hs(ilon,ilat,k) - temp_avg_hs(ilat,k))**2
      enddo
      temp_var_hs(ilat,k) = store1 / 143.
!!!!!!!!!!!!!!!!! TIME AVERAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      temp_tvar_hs (ilat,k) = temp_tvar_hs(ilat,k) * atime / (atime + 1.)  &
                            + temp_var_hs (ilat,k) / (atime + 1.)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   enddo
enddo

! Plot eddy variance of temperature

call plotback()
call cpcnrc(temp_var_hs,  kzdt, mzdt, nzhs  &
   , flow, fhgh, 5.0, nset, -1, ndsh)
call frame

call plotback()
call cpcnrc(temp_tvar_hs,  kzdt, mzdt, nzhs  &
   , flow, fhgh, 5.0, nset, -1, ndsh)
call frame

! Compute eddy variance of zonal wind

do k = 1,nzhs
   do ilat = 1,90
      do ilon = 1,144
         uzonal_hs(ilon,ilat,k)  &
            = (uzonal_hs(ilon,ilat,k) - uzonal_avg_hs(ilat,k)) ** 2
      enddo
   enddo
enddo

! Apply FFT in zonal direction for eddy variance of zonal wind

call fft99(uzonal_hs,work,trigs,ifax,inc,jump,n,lot,isign)

! Compute amplitude of each wavenumber and average in height

uspectra_hs(1:16,1:90) = 0.

do ilat = 1,90
! sum avg mass in latitude band
   bandmass = 0.
   do k = 1,nzhs
      bandmass = bandmass + rho_avg_hs(ilat,k) * dzt(k+1)
   enddo
! reconstruct and sum spectra
   do iwave = 2,16
      do k = 1,nzhs
         uspectra_hs(iwave,ilat) = uspectra_hs(iwave,ilat)  &
             + sqrt(uzonal_hs(2*iwave-1,ilat,k)**2    &
              + uzonal_hs(2*iwave  ,ilat,k)**2)   &
! Weighting by mass at k level
             * rho_avg_hs(ilat,k) * dzt(k+1)
      enddo
! Normalize by mass in vertical      
      uspectra_hs(iwave,ilat) = uspectra_hs(iwave,ilat) / bandmass
!!!!!!!!!!!!!!!!! TIME AVERAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      uspectrat_hs(iwave,ilat) = uspectrat_hs(iwave,ilat) * atime / (atime + 1.)  &
                               + uspectra_hs(iwave,ilat) / (atime + 1.)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   enddo
enddo

! Plot latitude-height dependence of amplitude
   
call plotback()
call cpcnrc(uspectra_hs, 16, 16, 90  &
   , flow, fhgh, 2.0, nset, -1, ndsh)
call frame

call plotback()
call cpcnrc(uspectrat_hs, 16, 16, 90  &
   , flow, fhgh, 2.0, nset, -1, ndsh)
call frame

! Compute and plot all time averages

atime = atime + 1.

return
end

!**************************************************************************

subroutine fft99(a,work,trigs,ifax,inc,jump,n,lot,isign)


!
! PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE
!              WILL PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
!              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
!              TRANSFORMS, I.E.  GIVEN A SET OF REAL DATA VECTORS, THE
!              PACKAGE RETURNS A SET OF 'HALF-COMPLEX' FOURIER
!              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE
!              TRANSFORMS MUST BE AN EVEN NUMBER GREATER THAN 4 THAT HAS
!              NO OTHER FACTORS EXCEPT POSSIBLY POWERS OF 2, 3, AND 5.
!              THIS IS AN ALL FORTRAN VERSION OF THE CRAYLIB PACKAGE
!              THAT IS MOSTLY WRITTEN IN CAL.
!
!              THE PACKAGE FFT99F CONTAINS SEVERAL USER-LEVEL ROUTINES:

!            SUBROUTINE FFTFAX
!                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE
!                BEFORE A SEQUENCE OF CALLS TO THE FFT ROUTINES
!                (PROVIDED THAT N IS NOT CHANGED).

!            SUBROUTINES FFT99 AND FFT991
!                TWO FFT ROUTINES THAT RETURN SLIGHTLY DIFFERENT
!                ARRANGEMENTS OF THE DATA IN GRIDPOINT SPACE.


! ACCESS       THIS FORTRAN VERSION MAY BE ACCESSED WITH

!                   *FORTRAN,P=XLIB,SN=FFT99F

!              TO ACCESS THE CRAY OBJECT CODE, CALLING THE USER ENTRY
!              POINTS FROM A CRAY PROGRAM IS SUFFICIENT.  THE SOURCE
!              FORTRAN AND CAL CODE FOR THE CRAYLIB VERSION MAY BE
!              ACCESSED USING

!                   FETCH P=CRAYLIB,SN=FFT99
!                   FETCH P=CRAYLIB,SN=CAL99

! USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 1,
!              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF
!              CALLS TO TRANSFORM A GIVEN SET OF REAL VECTORS OF LENGTH
!              N TO A SET OF 'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS
!              OF LENGTH N IS

!                   DIMENSION IFAX(13),TRIGS(3*N/2+1),A(M*(N+2)),
!                  +          WORK(M*(N+1))

!                   CALL FFTFAX (N, IFAX, TRIGS)
!                   CALL FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)

!              SEE THE INDIVIDUAL WRITE-UPS FOR FFTFAX, FFT99, AND
!              FFT991 BELOW, FOR A DETAILED DESCRIPTION OF THE
!              ARGUMENTS.

! HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN
!              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED
!              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.

!-----------------------------------------------------------------------

! SUBROUTINE FFTFAX (N,IFAX,TRIGS)

! PURPOSE      A SET-UP ROUTINE FOR FFT99 AND FFT991.  IT NEED ONLY BE
!              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO THE FFT
!              ROUTINES (PROVIDED THAT N IS NOT CHANGED).

! ARGUMENT     IFAX(13),TRIGS(3*N/2+1)
! DIMENSIONS

! ARGUMENTS

! ON INPUT     N
!               AN EVEN NUMBER GREATER THAN 4 THAT HAS NO PRIME FACTOR
!               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE
!               THE DOCUMENTATION FOR FFT99 AND FFT991 FOR THE
!               DEFINITIONS OF THE TRANSFORMS).

!              IFAX
!               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED
!               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING
!               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN A MILLION.

!              TRIGS
!               A FLOATING POINT ARRAY OF DIMENSION 3*N/2 IF N/2 IS
!               EVEN, OR 3*N/2+1 IF N/2 IS ODD.

! ON OUTPUT    IFAX
!               CONTAINS THE FACTORIZATION OF N/2.  IFAX(1) IS THE
!               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
!               IN IFAX(2),IFAX(3),...  IF FFTFAX IS CALLED WITH N ODD,
!               OR IF N HAS ANY PRIME FACTORS GREATER THAN 5, IFAX(1)
!               IS SET TO -99.

!              TRIGS
!               AN ARRAY OF TRIGNOMENTRIC FUNCTION VALUES SUBSEQUENTLY
!               USED BY THE FFT ROUTINES.

!-----------------------------------------------------------------------

! SUBROUTINE FFT991 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
!                       AND
! SUBROUTINE FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)

! PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
!              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
!              TRANSFORMS, USING ORDINARY SPATIAL ORDER OF GRIDPOINT
!              VALUES (FFT991) OR EXPLICIT CYCLIC CONTINUITY IN THE
!              OF REAL DATA VECTORS, THE PACKAGE RETURNS A SET OF
!              'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS, OR VICE
!              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE AN EVEN
!              NUMBER THAT HAS NO OTHER FACTORS EXCEPT POSSIBLY POWERS
!              OF 2, 3, AND 5.  THESE VERSION OF FFT991 AND FFT99 ARE
!              OPTIMIZED FOR USE ON THE CRAY-1.
!
! ARGUMENT     A(M*(N+2)), WORK(M*(N+1)), TRIGS(3*N/2+1), IFAX(13)
! DIMENSIONS
!
! ARGUMENTS
!
! ON INPUT     A
!               AN ARRAY OF LENGTH M*(N+2) CONTAINING THE INPUT DATA
!               OR COEFFICIENT VECTORS.  THIS ARRAY IS OVERWRITTEN BY
!               THE RESULTS.
!
!              WORK
!               A WORK ARRAY OF DIMENSION M*(N+1)
!
!              TRIGS
!               AN ARRAY SET UP BY FFTFAX, WHICH MUST BE CALLED FIRST.
!
!               AN ARRAY SET UP BY FFTFAX, WHICH MUST BE CALLED FIRST.
!
!              INC
!               THE INCREMENT (IN WORDS) BETWEEN SUCCESSIVE ELEMENTS OF
!               EACH DATA OR COEFFICIENT VECTOR (E.G.  INC=1 FOR
!               CONSECUTIVELY STORED DATA).
!
!              JUMP
!               THE INCREMENT (IN WORDS) BETWEEN THE FIRST ELEMENTS OF
!               SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-1,
!               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8
!               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF
!               INC AND JUMP, SEE THE EXAMPLES BELOW.
!
!              N
!               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF
!               TRANSFORMS, BELOW).
!
!              M
!               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.
!
!              ISIGN
!               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO
!                    GRIDPOINT VALUES.
!               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER
!                    COEFFICIENTS.
!
! ON OUTPUT    A
!               IF ISIGN = +1, AND M COEFFICIENT VECTORS ARE SUPPLIED
!               EACH CONTAINING THE SEQUENCE:
!
!               A(0),B(0),A(1),B(1),...,A(N/2),B(N/2)  (N+2 VALUES)
!
!               THEN THE RESULT CONSISTS OF M DATA VECTORS EACH
!               CONTAINING THE CORRESPONDING N+2 GRIDPOINT VALUES:
!
!               FOR FFT991, X(0), X(1), X(2),...,X(N-1),0,0.
!               FOR FFT99, X(N-1),X(0),X(1),X(2),...,X(N-1),X(0).
!                   (EXPLICIT CYCLIC CONTINUITY)
!
!               WHEN ISIGN = +1, THE TRANSFORM IS DEFINED BY:
!                 X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!                 WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!                 AND I=SQRT (-1)
!
!               IF ISIGN = -1, AND M DATA VECTORS ARE SUPPLIED EACH
!               CONTAINING A SEQUENCE OF GRIDPOINT VALUES X(J) AS
!               DEFINED ABOVE, THEN THE RESULT CONSISTS OF M VECTORS
!               EACH CONTAINING THE CORRESPONDING FOURIER COFFICIENTS
!               A(K), B(K), 0 .LE. K .LE N/2.
!
!               WHEN ISIGN = -1, THE INVERSE TRANSFORM IS DEFINED BY:
!                 C(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*EXP(-2*I*J*K*PI/N))
!                 WHERE C(K)=A(K)+I*B(K) AND I=SQRT(-1)
!
!               A CALL WITH ISIGN=+1 FOLLOWED BY A CALL WITH ISIGN=-1
!               (OR VICE VERSA) RETURNS THE ORIGINAL DATA.
!
!               NOTE: THE FACT THAT THE GRIDPOINT VALUES X(J) ARE REAL
!               IMPLIES THAT B(0)=B(N/2)=0.  FOR A CALL WITH ISIGN=+1,
!               IT IS NOT ACTUALLY NECESSARY TO SUPPLY THESE ZEROS.
!
! EXAMPLES      GIVEN 19 DATA VECTORS EACH OF LENGTH 64 (+2 FOR EXPLICIT
!               CYCLIC CONTINUITY), COMPUTE THE CORRESPONDING VECTORS OF
!               FOURIER COEFFICIENTS.  THE DATA MAY, FOR EXAMPLE, BE
!               ARRANGED LIKE THIS:
!
! FIRST DATA   A(1)=    . . .                A(66)=             A(70)
! VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
!
! SECOND DATA  A(71)=   . . .                                  A(140)
! VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
!
!               AND SO ON.  HERE INC=1, JUMP=70, N=64, M=19, ISIGN=-1,
!               AND FFT99 SHOULD BE USED (BECAUSE OF THE EXPLICIT CYCLIC
!               CONTINUITY).
!
!               ALTERNATIVELY THE DATA MAY BE ARRANGED LIKE THIS:
!
!                FIRST         SECOND                          LAST
!                DATA          DATA                            DATA
!                VECTOR        VECTOR                          VECTOR
!
!                 A(1)=         A(2)=                           A(19)=
!
!                 X(63)         X(63)       . . .               X(63)
!        A(20)=   X(0)          X(0)        . . .               X(0)
!        A(39)=   X(1)          X(1)        . . .               X(1)
!                  .             .                               .
!                  .             .                               .
!                  .             .                               .
!
!               IN WHICH CASE WE HAVE INC=19, JUMP=1, AND THE REMAINING
!               PARAMETERS ARE THE SAME AS BEFORE.  IN EITHER CASE, EACH
!               COEFFICIENT VECTOR OVERWRITES THE CORRESPONDING INPUT
!               DATA VECTOR.
!
!-----------------------------------------------------------------------
      real, dimension(n) :: A,WORK,TRIGS
      integer, dimension(*) :: IFAX

!     SUBROUTINE "FFT99" - MULTIPLE FAST REAL PERIODIC TRANSFORM
!     CORRESPONDING TO OLD SCALAR ROUTINE FFT9
!     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
!     (1970), 315-337)
!
!     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*LOT
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(N-1),X(0),X(1),X(2),...,X(N),X(0)
!         I.E. EXPLICIT CYCLIC CONTINUITY; (N+2) LOCATIONS REQUIRED
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!     PARALLEL
!
!     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!

! THE FOLLOWING CALL IS FOR MONITORING LIBRARY USE AT NCAR
      NFAX=IFAX(1)
      NX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30

!     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=INC+1
      JBASE=1
      DO L=1,LOT
         I=IBASE
         J=JBASE
         DO M=1,N
            WORK(J)=A(I)
            I=I+INC
            J=J+1
         enddo
         IBASE=IBASE+JUMP
         JBASE=JBASE+NX
      enddo

      IGO=60
      GO TO 40

!     PREPROCESSING (ISIGN=+1)
!     ------------------------

   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60

!     COMPLEX TRANSFORM
!     -----------------

   40 CONTINUE
      IA=INC+1
      LA=1
      DO K=1,NFAX
         IF (IGO.EQ.60) GO TO 60
   50    CONTINUE
         CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,  &
           INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
         IGO=60
         GO TO 70
   60    CONTINUE
         CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,  &
            2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
         IGO=50
   70    CONTINUE
         LA=LA*IFAX(K+1)
      enddo

      IF (ISIGN.EQ.-1) GO TO 130

!     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=IA
      DO L=1,LOT
         I=IBASE
         J=JBASE

         DO M=1,N
            A(J)=WORK(I)
            I=I+1
            J=J+INC
         enddo
         IBASE=IBASE+NX
         JBASE=JBASE+JUMP
      enddo

!     FILL IN CYCLIC BOUNDARY POINTS
  110 CONTINUE
      IA=1
      IB=N*INC+1

      DO L=1,LOT
         A(IA)=A(IB)
         A(IB+INC)=A(IA+INC)
         IA=IA+JUMP
         IB=IB+JUMP
      enddo
      GO TO 140

!     POSTPROCESSING (ISIGN=-1):
!     --------------------------

  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)

  140 CONTINUE
      RETURN
      END

!**************************************************************************

      SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      real, dimension(n) :: A,WORK,TRIGS

!     SUBROUTINE FFT99A - PREPROCESSING STEP FOR FFT99, ISIGN=+1
!     (SPECTRAL TO GRIDPOINT TRANSFORM)

      NH=N/2
      NX=N+1
      INK=INC+INC

!     A(0) AND A(N/2)
      IA=1
      IB=N*INC+1
      JA=1
      JB=2

      DO L=1,LOT
         WORK(JA)=A(IA)+A(IB)
         WORK(JB)=A(IA)-A(IB)
         IA=IA+JUMP
         IB=IB+JUMP
         JA=JA+NX
         JB=JB+NX
      enddo

!     REMAINING WAVENUMBERS
      IABASE=2*INC+1
      IBBASE=(N-2)*INC+1
      JABASE=3
      JBBASE=N-1

      DO K=3,NH,2
         IA=IABASE
         IB=IBBASE
         JA=JABASE
         JB=JBBASE
         C=TRIGS(N+K)
         S=TRIGS(N+K+1)

         DO L=1,LOT
            WORK(JA)=(A(IA)+A(IB))-  &
               (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
            WORK(JB)=(A(IA)+A(IB))+  &
               (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
            WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+  &
               (A(IA+INC)-A(IB+INC))
            WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))-  &
               (A(IA+INC)-A(IB+INC))
            IA=IA+JUMP
            IB=IB+JUMP
            JA=JA+NX
            JB=JB+NX
         enddo
         IABASE=IABASE+INK
         IBBASE=IBBASE-INK
         JABASE=JABASE+2
         JBBASE=JBBASE-2
      enddo

      IF (IABASE.NE.IBBASE) GO TO 50
!     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE

      DO L=1,LOT
         WORK(JA)=2.0*A(IA)
         WORK(JA+1)=-2.0*A(IA+INC)
         IA=IA+JUMP
         JA=JA+NX
      enddo

   50 CONTINUE
      RETURN
      END

!**************************************************************************

      SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
      DIMENSION WORK(N),A(N),TRIGS(N)

!     SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN=-1
!     (GRIDPOINT TO SPECTRAL TRANSFORM)

      NH=N/2
      NX=N+1
      INK=INC+INC

!     A(0) AND A(N/2)
      SCALE=1.0/FLOAT(N)
      IA=1
      IB=2
      JA=1
      JB=N*INC+1

      DO L=1,LOT
         A(JA)=SCALE*(WORK(IA)+WORK(IB))
         A(JB)=SCALE*(WORK(IA)-WORK(IB))
         A(JA+INC)=0.0
         A(JB+INC)=0.0
         IA=IA+NX
         IB=IB+NX
         JA=JA+JUMP
         JB=JB+JUMP
      enddo

!     REMAINING WAVENUMBERS
      SCALE=0.5*SCALE
      IABASE=3
      IBBASE=N-1
      JABASE=2*INC+1
      JBBASE=(N-2)*INC+1
!
      DO K=3,NH,2
         IA=IABASE
         IB=IBBASE
         JA=JABASE
         JB=JBBASE
         C=TRIGS(N+K)
         S=TRIGS(N+K+1)

         DO L=1,LOT
            A(JA)=SCALE*((WORK(IA)+WORK(IB))  &
              +(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
            A(JB)=SCALE*((WORK(IA)+WORK(IB))  &
              -(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
            A(JA+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))  &
               +(WORK(IB+1)-WORK(IA+1)))
            A(JB+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))  &
               -(WORK(IB+1)-WORK(IA+1)))
            IA=IA+NX
            IB=IB+NX
            JA=JA+JUMP
            JB=JB+JUMP
         enddo
         IABASE=IABASE+2
         IBBASE=IBBASE-2
         JABASE=JABASE+INK
         JBBASE=JBBASE-INK
      enddo

      IF (IABASE.NE.IBBASE) GO TO 50
!     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      SCALE=2.0*SCALE

      DO L=1,LOT
         A(JA)=SCALE*WORK(IA)
         A(JA+INC)=-SCALE*WORK(IA+1)
         IA=IA+NX
         JA=JA+JUMP
      enddo

   50 CONTINUE
      RETURN
      END

!**************************************************************************

      SUBROUTINE FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
      DIMENSION A(N),WORK(N),TRIGS(N),IFAX(1)

!     SUBROUTINE "FFT991" - MULTIPLE REAL/HALF-COMPLEX PERIODIC
!     FAST FOURIER TRANSFORM

!     SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO
!     THAT IN MRFFT2

!     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
!     (1970), 315-337)

!     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*LOT
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL

!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED

!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1)

!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!     PARALLEL

!     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER

!     DEFINITION OF TRANSFORMS:
!     -------------------------

!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)

!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))

! THE FOLLOWING CALL IS FOR MONITORING LIBRARY USE AT NCAR
      NFAX=IFAX(1)
      NX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30

!     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=1
      JBASE=1
      DO L=1,LOT
         I=IBASE
         J=JBASE

         DO M=1,N
            WORK(J)=A(I)
            I=I+INC
           J=J+1
         enddo
         IBASE=IBASE+JUMP
         JBASE=JBASE+NX
      enddo

      IGO=60
      GO TO 40

!     PREPROCESSING (ISIGN=+1)
!     ------------------------

   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60

!     COMPLEX TRANSFORM
!     -----------------

   40 CONTINUE
      IA=1
      LA=1
      DO K=1,NFAX
         IF (IGO.EQ.60) GO TO 60
   50    CONTINUE
         CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,  &
           INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
         IGO=60
         GO TO 70
   60    CONTINUE
         CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,  &
            2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
         IGO=50
   70    CONTINUE
         LA=LA*IFAX(K+1)
      enddo

      IF (ISIGN.EQ.-1) GO TO 130

!     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=1
      DO L=1,LOT
         I=IBASE
         J=JBASE

         DO M=1,N
            A(J)=WORK(I)
            I=I+1
            J=J+INC
         enddo
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
     enddo

!     FILL IN ZEROS AT END
  110 CONTINUE
      IB=N*INC+1

      DO L=1,LOT
         A(IB)=0.0
         A(IB+INC)=0.0
         IB=IB+JUMP
      enddo
      GO TO 140

!     POSTPROCESSING (ISIGN=-1):
!     --------------------------

  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)

  140 CONTINUE
      RETURN
      END

!**************************************************************************

      SUBROUTINE FFTFAX(N,IFAX,TRIGS)
      DIMENSION IFAX(13),TRIGS(1)

! MODE 3 IS USED FOR REAL/HALF-COMPLEX TRANSFORMS.  IT IS POSSIBLE
! TO DO COMPLEX/COMPLEX TRANSFORMS WITH OTHER VALUES OF MODE, BUT
! DOCUMENTATION OF THE DETAILS WERE NOT AVAILABLE WHEN THIS ROUTINE
! WAS WRITTEN.

      DATA MODE /3/
      CALL FAX (IFAX, N, MODE)
      I = IFAX(1)
      IF (IFAX(I+1) .GT. 5 .OR. N .LE. 4) IFAX(1) = -99
      IF (IFAX(1) .LE. 0 )PRINT*,' FFTFAX - INVALID N ',N
      CALL FFTRIG (TRIGS, N, MODE)
      RETURN
      END

!**************************************************************************

      SUBROUTINE FAX(IFAX,N,MODE)
      DIMENSION IFAX(10)
      NN=N
      IF (IABS(MODE).EQ.1) GO TO 10
      IF (IABS(MODE).EQ.8) GO TO 10
      NN=N/2
      IF ((NN+NN).EQ.N) GO TO 10
      IFAX(1)=-99
      RETURN
   10 K=1
!     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
!     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
!     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
!     NOW FIND REMAINING FACTORS
   50 L=5
      INC=2
!     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 L=L+INC
      INC=6-INC
      GO TO 60
   80 IFAX(1)=K-1
!     IFAX(1) CONTAINS NUMBER OF FACTORS
      NFAX=IFAX(1)
!     SORT FACTORS INTO ASCENDING ORDER
      IF (NFAX.EQ.1) GO TO 110
      DO II=2,NFAX
         ISTOP=NFAX+2-II
         DO I=2,ISTOP
            IF (IFAX(I+1).GE.IFAX(I)) GO TO 90
            ITEM=IFAX(I)
            IFAX(I)=IFAX(I+1)
            IFAX(I+1)=ITEM
   90       continue
         enddo
      enddo
  110 CONTINUE
      RETURN
      END

!**************************************************************************

      SUBROUTINE FFTRIG(TRIGS,N,MODE)
      DIMENSION TRIGS(1)
      PI=2.0*ASIN(1.0)
      IMODE=IABS(MODE)
      NN=N
      IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2
      DEL=(PI+PI)/FLOAT(NN)
      L=NN+NN
      DO I=1,L,2
         ANGLE=0.5*FLOAT(I-1)*DEL
         TRIGS(I)=COS(ANGLE)
         TRIGS(I+1)=SIN(ANGLE)
      enddo
      IF (IMODE.EQ.1) RETURN
      IF (IMODE.EQ.8) RETURN
      DEL=0.5*DEL
      NH=(NN+1)/2
      L=NH+NH
      LA=NN+NN
      DO I=1,L,2
         ANGLE=0.5*FLOAT(I-1)*DEL
         TRIGS(LA+I)=COS(ANGLE)
         TRIGS(LA+I+1)=SIN(ANGLE)
      enddo
      IF (IMODE.LE.3) RETURN
      DEL=0.5*DEL
      LA=LA+NN
      IF (MODE.EQ.5) GO TO 40
      DO I=2,NN
         ANGLE=FLOAT(I-1)*DEL
         TRIGS(LA+I)=2.0*SIN(ANGLE)
      enddo
      RETURN
   40 CONTINUE
      DEL=0.5*DEL
      DO I=2,N
         ANGLE=FLOAT(I-1)*DEL
         TRIGS(LA+I)=SIN(ANGLE)
      enddo
      RETURN
      END

!**************************************************************************

      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)
      DIMENSION A(N),B(N),C(N),D(N),TRIGS(N)

!     SUBROUTINE "VPASSM" - MULTIPLE VERSION OF "VPASSA"
!     PERFORMS ONE PASS THROUGH DATA
!     AS PART OF MULTIPLE COMPLEX FFT ROUTINE
!     A IS FIRST REAL INPUT VECTOR
!     B IS FIRST IMAGINARY INPUT VECTOR
!     C IS FIRST REAL OUTPUT VECTOR
!     D IS FIRST IMAGINARY OUTPUT VECTOR
!     TRIGS IS PRECALCULATED TABLE OF SINES " COSINES
!     INC1 IS ADDRESSING INCREMENT FOR A AND B
!     INC2 IS ADDRESSING INCREMENT FOR C AND D
!     INC3 IS ADDRESSING INCREMENT BETWEEN A"S & B"S
!     INC4 IS ADDRESSING INCREMENT BETWEEN C"S & D"S
!     LOT IS THE NUMBER OF VECTORS
!     N IS LENGTH OF VECTORS
!     IFAC IS CURRENT FACTOR OF N
!     LA IS PRODUCT OF PREVIOUS FACTORS

      DATA SIN36/0.587785252292473/,COS36/0.809016994374947/,  &
          SIN72/0.951056516295154/,COS72/0.309016994374947/,  &
          SIN60/0.866025403784437/

      M=N/IFAC
      IINK=M*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.GT.4) RETURN
      GO TO (10,50,90,130),IGO

!     CODING FOR FACTOR 2

   10 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      DO L=1,LA
         I=IBASE
         J=JBASE

         DO IJK=1,LOT
            C(JA+J)=A(IA+I)+A(IB+I)
            D(JA+J)=B(IA+I)+B(IB+I)
            C(JB+J)=A(IA+I)-A(IB+I)
            D(JB+J)=B(IA+I)-B(IB+I)
            I=I+INC3
            J=J+INC4
         enddo
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      enddo
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO K=LA1,M,LA
         KB=K+K-2
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         DO L=1,LA
            I=IBASE
            J=JBASE

            DO IJK=1,LOT
               C(JA+J)=A(IA+I)+A(IB+I)
               D(JA+J)=B(IA+I)+B(IB+I)
               C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
               D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
               I=I+INC3
               J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         enddo
         JBASE=JBASE+JUMP
      enddo
      RETURN

!     CODING FOR FACTOR 3

   50 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      DO L=1,LA
         I=IBASE
         J=JBASE

         DO IJK=1,LOT
            C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
            D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
            C(JB+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))
            C(JC+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))
            D(JB+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I)))
            D(JC+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I)))
            I=I+INC3
            J=J+INC4
         enddo
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      enddo
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO K=LA1,M,LA
         KB=K+K-2
         KC=KB+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         DO L=1,LA
            I=IBASE
            J=JBASE

            DO IJK=1,LOT
               C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
               D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
               C(JB+J)=  &
                  C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))  &
                 -S1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
               D(JB+J)=  &
                  S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))  &
                 +C1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
               C(JC+J)=  &
                  C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))  &
                 -S2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
               D(JC+J)=  &
                  S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))  &
                 +C2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
               I=I+INC3
               J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         enddo
         JBASE=JBASE+JUMP
      enddo
      RETURN

!     CODING FOR FACTOR 4

   90 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      DO L=1,LA
         I=IBASE
         J=JBASE

         DO IJK=1,LOT
            C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
            C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
            D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
            D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
            C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
            C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
            D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
            D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
            I=I+INC3
            J=J+INC4
         enddo
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      enddo
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO K=LA1,M,LA
         KB=K+K-2
         KC=KB+KB
         KD=KC+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         C3=TRIGS(KD+1)
         S3=TRIGS(KD+2)
         DO L=1,LA
            I=IBASE
            J=JBASE

            DO IJK=1,LOT
               C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
               D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
               C(JC+J)=  &
                  C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))  &
                 -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
               D(JC+J)=  &
                  S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))  &
                 +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
               C(JB+J)=  &
                  C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))  &
                 -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
               D(JB+J)=  &
                  S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))  &
                 +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
               C(JD+J)=  &
                  C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))  &
                 -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
               D(JD+J)=  &
                  S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))  &
                 +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
               I=I+INC3
               J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         enddo
         JBASE=JBASE+JUMP
      enddo
      RETURN

!     CODING FOR FACTOR 5

  130 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      IE=ID+IINK
      JE=JD+JINK
      DO L=1,LA
         I=IBASE
         J=JBASE

         DO IJK=1,LOT
            C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
            D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
            C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))  &
             -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
            C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))  &
             +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
            D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))  &
             +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
            D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))  &
             -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
            C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))  &
             -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
            C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))  &
             +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
            D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))  &
             +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
            D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))  &
             -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
            I=I+INC3
            J=J+INC4
         enddo
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      enddo
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO K=LA1,M,LA
         KB=K+K-2
         KC=KB+KB
         KD=KC+KB
         KE=KD+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         C3=TRIGS(KD+1)
         S3=TRIGS(KD+2)
         C4=TRIGS(KE+1)
         S4=TRIGS(KE+2)
         DO L=1,LA
            I=IBASE
            J=JBASE

            DO IJK=1,LOT
               C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
               D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
               C(JB+J)=  &
                  C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))  &
                    -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))  &
                 -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))  &
                    +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
               D(JB+J)=  &
                  S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))  &
                    -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))  &
                 +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))  &
                    +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
               C(JE+J)=  &
                  C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))  &
                    +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))  &
                 -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))  &
                    -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
               D(JE+J)=  &
                  S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))  &
                   +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))  &
                +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))  &
                   -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
               C(JC+J)=  &
                  C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))  &
                    -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))  &
                 -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))  &
                    +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
               D(JC+J)=  &
                  S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))  &
                    -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))  &
                 +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))  &
                    +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
               C(JD+J)=  &
                  C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))  &
                    +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))  &
                 -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))  &
                    -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
               D(JD+J)=  &
                  S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))  &
                    +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))  &
                 +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))  &
                    -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
               I=I+INC3
               J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         enddo
         JBASE=JBASE+JUMP
      enddo
      RETURN
      END

