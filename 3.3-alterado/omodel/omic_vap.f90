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
subroutine thrmstr(iw0,lpw0,k1,k2  &
   ,press0,thil0,rhow,rhoi,exner0,tair,theta0,rhov,rhovstr,tairstrc  &
   ,rx,qx,sa)

use micro_coms,  only: mza0, ncat
use consts_coms, only: p00i, rocp, alvl, alvi, cpi4, cpi, cp253i
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iw0,lpw0

integer, intent(in) :: k1(10)
integer, intent(in) :: k2(10)

real, intent(in)  :: press0   (mza0)
real, intent(in ) :: thil0    (mza0)
real, intent(in ) :: rhoi     (mza0)
real, intent(out) :: exner0   (mza0)
real, intent(out) :: tair     (mza0)
real, intent(out) :: theta0   (mza0)
real, intent(out) :: rhov     (mza0)
real, intent(out) :: rhovstr  (mza0)
real, intent(out) :: tairstrc (mza0)

real(kind=8), intent(in) :: rhow(mza0)

real, intent(in) :: rx        (mza0,ncat)
real, intent(in) :: qx        (mza0,ncat)

real, intent(out) :: sa       (mza0,9)

integer :: k,lcat
real :: fracliq,tcoal,tairstr

real :: rholiq (mza0) ! automatic array
real :: rhoice (mza0) ! automatic array
real :: qhydm  (mza0) ! automatic array
real :: til    (mza0) ! automatic array

do k = lpw0,mza0
   exner0(k) = (press0(k) * p00i) ** rocp  ! defined WITHOUT CP factor
   theta0(k) = thil0(k)
   tair(k) = theta0(k) * exner0(k)
   rhov(k) = rhow(k)
enddo

do k = k1(10),k2(10)
   til(k) = thil0(k) * exner0(k)
   rholiq(k) = 0.
   rhoice(k) = 0.
enddo

do lcat = 1,2
   do k = k1(lcat),k2(lcat)
      rholiq(k) = rholiq(k) + rx(k,lcat)
   enddo
enddo

do lcat = 3,5
   do k = k1(lcat),k2(lcat)
      rhoice(k) = rhoice(k) + rx(k,lcat)
   enddo
enddo

do lcat = 6,7
   do k = k1(lcat),k2(lcat)
      call qtc(qx(k,lcat),tcoal,fracliq)
      rholiq(k) = rholiq(k) + rx(k,lcat) * fracliq
      rhoice(k) = rhoice(k) + rx(k,lcat) * (1. - fracliq)
   enddo
enddo

do k = k1(10),k2(10)
   qhydm(k) = alvl * rholiq(k) + alvi * rhoice(k)
   rhovstr(k) = rhow(k) - rholiq(k) - rhoice(k)
   sa(k,1) = til(k) * qhydm(k) / (1.e-12 + rholiq(k) + rhoice(k))  ! stays the same
enddo

do k = k1(10),k2(10)
   if (tair(k) > 253.) then
!     tairstr = .5 * (til(k) + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
! Change in tairstr computation since qhydm is now J/m^3 instead of J/kg:    
      tairstr = .5  &
         * (til(k) + sqrt(til(k) * (til(k) + cpi4 * qhydm(k) * rhoi(k))))
      sa(k,1) = sa(k,1) * cpi / (2. * tairstr - til(k))  ! stays the same
   else
!     tairstr = til(k) * (1. + qhydm(k) * cp253i)
! Change in tairstr computation since qhydm is now J/m^3 instead of J/kg:    
      tairstr = til(k) * (1. + qhydm(k) * rhoi(k) * cp253i)
      sa(k,1) = sa(k,1) * cp253i  ! stays the same
  endif
  tairstrc(k) = tairstr - 273.15
enddo

return
end subroutine thrmstr

!===============================================================================

subroutine diffprep(lcat,k1,k2                    &
   ,pi4dt,jhcat,sa,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz  &
   ,rx,cx,qr,emb,rhoa,rhov,rhovsrefp,rdynvsci,vapdif,thrmcon,sumuy,sumuz)

use micro_coms, only: rxmin, frefac1, pwmasi, frefac2, cdp1, sl, sj, sc, sk,  &
                      mza0, ncat
use misc_coms,  only: io6

implicit none

integer, intent(in) :: lcat

integer, intent(in) :: k1(10)
integer, intent(in) :: k2(10)

real, intent(in) :: pi4dt

integer, intent(in) :: jhcat (mza0,ncat)

real, intent(in   ) :: sa  (mza0,9)

real, intent(out  ) :: sb  (mza0,ncat)
real, intent(out  ) :: sd  (mza0,ncat)
real, intent(out  ) :: se  (mza0,ncat)
real, intent(out  ) :: sf  (mza0,ncat)
real, intent(out  ) :: sg  (mza0,ncat)
real, intent(inout) :: sh  (mza0,ncat)
real, intent(inout) :: sm  (mza0,ncat)
real, intent(out  ) :: ss  (mza0,ncat)
real, intent(out  ) :: su  (mza0,ncat)
real, intent(out  ) :: sw  (mza0,ncat)
real, intent(out  ) :: sy  (mza0,ncat)
real, intent(out  ) :: sz  (mza0,ncat)

real, intent(in   ) :: rx  (mza0,ncat)
real, intent(in   ) :: cx  (mza0,ncat)
real, intent(in   ) :: qr  (mza0,ncat)
real, intent(in   ) :: emb (mza0,ncat)

real, intent(in   ) :: rhov      (mza0)
real, intent(in   ) :: rhovsrefp (mza0,2)
real, intent(in   ) :: rdynvsci  (mza0)
real, intent(in   ) :: vapdif    (mza0)
real, intent(in   ) :: thrmcon   (mza0)

real(kind=8), intent(in) :: rhoa(mza0)

real, intent(inout) :: sumuy     (mza0)
real, intent(inout) :: sumuz     (mza0)

integer :: k,mynum,if1,if4,if6,if8,lhcat
real :: fre,scdei

real, dimension(mza0) :: ttest ! automatic array

if (lcat <= 2) then
   if1 = 1
   if4 = 4
   if6 = 6
   if8 = 8
else
   if1 = 2
   if4 = 5
   if6 = 7
   if8 = 9
endif

do k = k1(lcat),k2(lcat)
   lhcat = jhcat(k,lcat)

   if (rx(k,lcat) < rxmin(lcat)) cycle

   fre = frefac1(lhcat) * emb(k,lcat) ** pwmasi(lhcat)  &
      + rdynvsci(k) * frefac2(lhcat) * emb(k,lcat) ** cdp1(lhcat)

   sb(k,lcat) = cx(k,lcat) * fre * pi4dt  ! stays the same (rhoa factor removed)
   su(k,lcat) = vapdif(k) * sb(k,lcat)    ! stays the same
   sd(k,lcat) = sh(k,lcat) * rx(k,lcat)   ! went up by rhoa factor
   se(k,lcat) = su(k,lcat) * sa(k,if6) + sb(k,lcat) * thrmcon(k) * rhoa(K)
    ! se went up by rhoa factor (rhoa factor had to be inserted in second term)
   sf(k,lcat) = su(k,lcat) * sl(if1) - sb(k,lcat) * sa(k,2)  ! stays the same
   sg(k,lcat) = su(k,lcat) * sa(k,if8) + sb(k,lcat) * sa(k,3)  &
              + sj(lcat) * qr(k,lcat) ! went up by rhoa factor
!     + lambda_j [Joules/m^3 added by radiative heating this timestep]
   scdei = 1. / (sc(if1) * sd(k,lcat) + se(k,lcat)) ! went down by rhoa factor
   ss(k,lcat) = sf(k,lcat) * scdei   ! went down by rhoa factor
   sw(k,lcat) = (sg(k,lcat) - sk(if1) * sd(k,lcat)) * scdei  ! stays the same
   ttest(k) = ss(k,lcat) * rhov(k) + sw(k,lcat)  ! stays the same

enddo

if (lcat >= 3 .and. lcat <= 5) then
   do k = k1(lcat),k2(lcat)
      if (rx(k,lcat) < rxmin(lcat)) cycle
      if (ttest(k) >= 0.) then
         sm(k,lcat) = 0.
         sh(k,lcat) = 1.
         sd(k,lcat) = sh(k,lcat) * rx(k,lcat)
         scdei = 1. / (sc(if1) * sd(k,lcat) + se(k,lcat))
         ss(k,lcat) = sf(k,lcat) * scdei
         sw(k,lcat) = (sg(k,lcat) - sk(if1) * sd(k,lcat)) * scdei
      else
         sm(k,lcat) = 1.
      endif
   enddo
endif

if (lcat >= 6) then
   do k = k1(lcat),k2(lcat)
      if (rx(k,lcat) < rxmin(lcat))cycle
      if (ttest(k) >= 0.) then
         sm(k,lcat) = 0.
      else
         sm(k,lcat) = 1.
      endif
   enddo
endif

do k = k1(lcat),k2(lcat)
   if (rx(k,lcat) < rxmin(lcat))cycle
   sy(k,lcat) = rhovsrefp(k,if1) * sm(k,lcat) * sw(k,lcat) - sa(k,if4)
   sz(k,lcat) = 1. - rhovsrefp(k,if1) * ss(k,lcat) * sm(k,lcat)
   sumuy(k) = sumuy(k) + su(k,lcat) * sy(k,lcat)
   sumuz(k) = sumuz(k) + su(k,lcat) * sz(k,lcat)
enddo

return
end subroutine diffprep

!===============================================================================

subroutine vapdiff (j1,j2,rhov,rhovstr,sumuy,sumuz)

use micro_coms, only: mza0
use misc_coms,  only: io6

implicit none

integer, intent(in) :: j1
integer, intent(in) :: j2

real, intent(out) :: rhov    (mza0)
real, intent(in ) :: rhovstr (mza0)
real, intent(in ) :: sumuy   (mza0)
real, intent(in ) :: sumuz   (mza0)

integer :: k

do k = j1,j2
   rhov(k) = (rhovstr(k) + sumuy(k)) / (1.0 + sumuz(k))
enddo

return
end subroutine vapdiff

!===============================================================================

subroutine vapflux(iw0,lcat,k1,k2  &
   ,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap  &
   ,rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr)

use micro_coms, only: mza0, ncat, rxmin, sc, sk
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iw0,lcat

integer, intent(in) :: k1(10)
integer, intent(in) :: k2(10)

real, intent(in   ) :: sa  (mza0,9)

real, intent(in   ) :: sd  (mza0,ncat)
real, intent(in   ) :: se  (mza0,ncat)
real, intent(in   ) :: sf  (mza0,ncat)
real, intent(in   ) :: sg  (mza0,ncat)
real, intent(in   ) :: sm  (mza0,ncat)
real, intent(in   ) :: ss  (mza0,ncat)
real, intent(in   ) :: su  (mza0,ncat)
real, intent(in   ) :: sw  (mza0,ncat)
real, intent(in   ) :: sy  (mza0,ncat)
real, intent(in   ) :: sz  (mza0,ncat)

real, intent(inout) :: rx  (mza0,ncat)
real, intent(inout) :: cx  (mza0,ncat)
real, intent(out  ) :: qx  (mza0,ncat)
real, intent(out  ) :: qr  (mza0,ncat)
real, intent(out  ) :: tx  (mza0,ncat)
real, intent(out  ) :: vap (mza0,ncat)

real, intent(in   ) :: rhovsrefp (mza0,2)
real, intent(inout) :: rhov      (mza0)
real, intent(in   ) :: rhovstr   (mza0)
real, intent(inout) :: sumuy     (mza0)
real, intent(inout) :: sumuz     (mza0)
real, intent(inout) :: sumvr     (mza0)

integer :: k,if1,if4
real :: rxx

if (lcat <= 2) then
   if1 = 1
   if4 = 4
else
   if1 = 2
   if4 = 5
endif

do k = k1(lcat),k2(lcat)

   if (rx(k,lcat) < rxmin(lcat)) cycle

   tx(k,lcat) = (ss(k,lcat) * rhov(k) + sw(k,lcat)) * sm(k,lcat)
   vap(k,lcat) = su(k,lcat) * (rhov(k) + sa(k,if4) - rhovsrefp(k,if1) * tx(k,lcat))

   if (vap(k,lcat) > -rx(k,lcat)) then

      rxx = rx(k,lcat) + vap(k,lcat)

      if (sm(k,lcat) > .5) then
         qx(k,lcat) = sc(if1) * tx(k,lcat) + sk(if1)
         qr(k,lcat) = qx(k,lcat) * rxx
      else
         qx(k,lcat) = (rhov(k) * sf(k,lcat) + sg(k,lcat)  &
                    - tx(k,lcat) * se(k,lcat)) / sd(k,lcat)
         qx(k,lcat) = min(350000.,max(-100000.,qx(k,lcat)))
         qr(k,lcat) = qx(k,lcat) * rxx
      endif

   endif

!bob Now also do the following section if pristine ice totally melts:
! evaporate it too.

   if ((lcat == 3 .and. qx(k,lcat) > 330000.) .or.  &
      vap(k,lcat) <= -rx(k,lcat)) then

      sumuy(k) = sumuy(k) - su(k,lcat) * sy(k,lcat)
      sumuz(k) = sumuz(k) - su(k,lcat) * sz(k,lcat)
      sumvr(k) = sumvr(k) + rx(k,lcat)
      rhov(k) = (rhovstr(k) + sumuy(k) + sumvr(k)) / (1.0 + sumuz(k))

      vap(k,lcat) = - rx(k,lcat)
      tx(k,lcat) = 0.
      rx(k,lcat) = 0.
      cx(k,lcat) = 0.
      qx(k,lcat) = 0.
      qr(k,lcat) = 0.
   else
      rx(k,lcat) = rxx
   endif

enddo

return
end subroutine vapflux

!===============================================================================

subroutine psxfer(iw0,j1,j2,jhcat,vap,rx,cx,qx,qr)

use micro_coms, only: rxmin, dnfac, pwmasi, gam, dps2, dps, gnu, gamn1,  &
                      pwmas, dpsmi, mza0, ncat, emb0, emb1
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iw0,j1,j2

integer, intent(in   ) :: jhcat (mza0,ncat)

real, intent(in   ) :: vap (mza0,ncat)
real, intent(inout) :: rx  (mza0,ncat)
real, intent(inout) :: cx  (mza0,ncat)
real, intent(in   ) :: qx  (mza0,ncat)
real, intent(inout) :: qr  (mza0,ncat)

integer :: k,lhcat,it
real :: embx,dn,xlim,dvap,dqr,dnum

do k = j1,j2

   if (vap(k,3) > 0.) then

      lhcat = jhcat(k,3)
      embx = max(emb0(3),min(emb1(3),rx(k,3) / max(1.e-9,cx(k,3))))

      dn = dnfac(lhcat) * embx ** pwmasi(lhcat)
      it = nint(dn * 1.e6)

      xlim = gam(it,3) * dps2 * (dps / dn) ** (gnu(3) - 1.)  &
         / (gamn1(3) * pwmas(lhcat) * dn ** 2)

      dvap = min(rx(k,3),vap(k,3) * (xlim + gam(it,1) / gamn1(3))) ! up by rhoa
      dqr = dvap * qx(k,3)                                         ! up by rhoa
      dnum = dvap * min(dpsmi(lhcat),1./embx)                      ! up by rhoa

      rx(k,3) = rx(k,3) - dvap
      cx(k,3) = cx(k,3) - dnum
      qr(k,3) = qr(k,3) - dqr
      rx(k,4) = rx(k,4) + dvap
      cx(k,4) = cx(k,4) + dnum
      qr(k,4) = qr(k,4) + dqr

   elseif (vap(k,4) < 0. .and. rx(k,4) >= 0.) then

      lhcat = jhcat(k,4)
      embx = max(emb0(4),min(emb1(4),rx(k,4) / max(1.e-3,cx(k,4))))

      dn = dnfac(lhcat) * embx ** pwmasi(lhcat)
      it = nint(dn * 1.e6)

      xlim = gam(it,3) * dps2 * (dps / dn) ** (gnu(4) - 1.)  &
         / (gamn1(4) * pwmas(lhcat) * dn ** 2)

      dvap = max(-rx(k,4),vap(k,4) * xlim)     ! up by rhoa
      dqr = dvap * qx(k,4)                     ! up by rhoa
      dnum = dvap * max(dpsmi(lhcat),1./embx)  ! up by rhoa

      rx(k,3) = rx(k,3) - dvap
      cx(k,3) = cx(k,3) - dnum
      qr(k,3) = qr(k,3) - dqr
      rx(k,4) = rx(k,4) + dvap
      cx(k,4) = cx(k,4) + dnum
      qr(k,4) = qr(k,4) + dqr

   endif

enddo

return
end subroutine psxfer

!===============================================================================

subroutine newtemp(j1,j2  &
   ,tairstrc,rhoi,rhovstr,rhov,exner0,tairc,tair,theta0,rhovslair,rhovsiair,sa)

use micro_coms, only: mza0
use misc_coms,  only: io6

implicit none

integer, intent(in) :: j1
integer, intent(in) :: j2

real, intent(in   ) :: tairstrc  (mza0)
real, intent(in   ) :: rhoi      (mza0)
real, intent(in   ) :: rhovstr   (mza0)
real, intent(in   ) :: rhov      (mza0)
real, intent(in   ) :: exner0    (mza0)
real, intent(out  ) :: tairc     (mza0)
real, intent(out  ) :: tair      (mza0)
real, intent(out  ) :: theta0    (mza0)
real, intent(out  ) :: rhovslair (mza0)
real, intent(out  ) :: rhovsiair (mza0)

real, intent(in   ) :: sa (mza0,9)

real, external :: rhovsl,rhovsi

integer :: k

do k = j1,j2
   tairc(k) = tairstrc(k) + sa(k,1) * rhoi(k) * (rhovstr(k) - rhov(k)) ! rhoi inserted
   tair(k)  = tairc(k) + 273.15
   theta0(k) = tair(k) / exner0(k)
   rhovslair(k) = rhovsl(tairc(k))
   rhovsiair(k) = rhovsi(tairc(k))
enddo

return
end subroutine newtemp
