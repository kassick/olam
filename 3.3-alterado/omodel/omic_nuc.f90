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
subroutine cldnuc(lpw0  &
   ,rx,cx,rhov,rhoa,tairc,cccnx,wc0,rhovslair)

use micro_coms, only: jnmb, parm, emb0, cparm, mza0, ncat, icloud
use misc_coms,  only: io6

implicit none

integer, intent(in) :: lpw0

real, intent(inout) :: rx   (mza0,ncat)
real, intent(inout) :: cx   (mza0,ncat)
real, intent(inout) :: rhov (mza0)

real, intent(in) :: tairc    (mza0)
real, intent(in) :: cccnx    (mza0)
real, intent(in) :: wc0      (mza0)
real, intent(in) :: rhovslair(mza0)

real(kind=8), intent(in) :: rhoa(mza0)

integer :: k,jtemp,jw,jconcen,iw,iconc
real :: rnuc,excessrhov,rhocnew,tab,concen_tab,cxadd  &
   ,tairc_nuc,w_nuc,rjw,wtw2,wtw1,concen_nuc,rjconcen,wtconcen2,wtconcen1

real, dimension(9,7,7) :: cldnuctab
data ((cldnuctab(iw,iconc,1),iw=1,9),iconc=1,7)/                 &
  .307,  .520,  .753,  .919,  .990,  .990,  .990,  .990,  .990,  &
  .230,  .426,  .643,  .860,  .969,  .990,  .990,  .990,  .990,  &
  .164,  .336,  .552,  .777,  .940,  .990,  .990,  .990,  .990,  &
  .098,  .254,  .457,  .701,  .892,  .979,  .990,  .990,  .990,  &
  .045,  .145,  .336,  .614,  .822,  .957,  .990,  .990,  .990,  &
  .018,  .073,  .206,  .426,  .672,  .877,  .969,  .990,  .990,  &
  .008,  .027,  .085,  .206,  .280,  .336,  .336,  .336,  .906/
  
data ((cldnuctab(iw,iconc,2),iw=1,9),iconc=1,7)/                 &
  .230,  .426,  .643,  .860,  .969,  .990,  .990,  .990,  .990,  &
  .164,  .336,  .552,  .777,  .930,  .990,  .990,  .990,  .990,  &
  .112,  .254,  .457,  .701,  .877,  .974,  .990,  .990,  .990,  &
  .073,  .184,  .365,  .583,  .822,  .949,  .990,  .990,  .990,  &
  .038,  .112,  .254,  .489,  .727,  .906,  .982,  .990,  .990,  &
  .015,  .054,  .145,  .365,  .614,  .841,  .957,  .990,  .990,  &
  .005,  .018,  .073,  .184,  .395,  .614,  .800,  .940,  .990/
  
data ((cldnuctab(iw,iconc,3),iw=1,9),iconc=1,7)/                 &
  .164,  .336,  .552,  .800,  .949,  .990,  .990,  .990,  .990,  &
  .128,  .254,  .457,  .701,  .892,  .979,  .990,  .990,  .990,  &
  .085,  .184,  .365,  .583,  .822,  .949,  .990,  .990,  .990,  &
  .054,  .128,  .280,  .489,  .727,  .906,  .982,  .990,  .990,  &
  .027,  .085,  .206,  .395,  .643,  .841,  .963,  .990,  .990,  &
  .012,  .038,  .112,  .280,  .520,  .777,  .930,  .990,  .990,  &
  .004,  .015,  .054,  .145,  .365,  .614,  .822,  .949,  .990/
  
data ((cldnuctab(iw,iconc,4),iw=1,9),iconc=1,7)/                 &
  .145,  .280,  .489,  .727,  .919,  .990,  .990,  .990,  .990,  &
  .098,  .206,  .395,  .614,  .841,  .963,  .990,  .990,  .990,  &
  .063,  .145,  .307,  .520,  .753,  .919,  .990,  .990,  .990,  &
  .038,  .098,  .230,  .426,  .643,  .860,  .963,  .990,  .990,  &
  .022,  .063,  .164,  .336,  .552,  .777,  .930,  .990,  .990,  &
  .010,  .027,  .085,  .230,  .426,  .701,  .877,  .974,  .990,  &
  .003,  .012,  .038,  .112,  .280,  .552,  .777,  .940,  .990/
  
data ((cldnuctab(iw,iconc,5),iw=1,9),iconc=1,7)/                 &
  .112,  .230,  .457,  .701,  .892,  .982,  .990,  .990,  .990,  &
  .073,  .164,  .336,  .552,  .800,  .940,  .990,  .990,  .990,  &
  .054,  .112,  .254,  .457,  .672,  .877,  .979,  .990,  .990,  &
  .032,  .085,  .184,  .365,  .583,  .800,  .940,  .990,  .990,  &
  .018,  .045,  .128,  .254,  .457,  .701,  .892,  .979,  .990,  &
  .008,  .022,  .073,  .184,  .365,  .614,  .822,  .949,  .990,  &
  .003,  .010,  .032,  .098,  .230,  .489,  .727,  .906,  .979/
  
data ((cldnuctab(iw,iconc,6),iw=1,9),iconc=1,7)/                 &
  .098,  .206,  .395,  .643,  .860,  .974,  .990,  .990,  .990,  &
  .063,  .145,  .307,  .520,  .753,  .930,  .990,  .990,  .990,  &
  .045,  .098,  .206,  .395,  .643,  .841,  .963,  .990,  .990,  &
  .027,  .063,  .145,  .307,  .520,  .753,  .919,  .990,  .990,  &
  .015,  .038,  .098,  .230,  .426,  .643,  .841,  .963,  .990,  &
  .007,  .018,  .054,  .145,  .307,  .552,  .777,  .919,  .990,  &
  .003,  .008,  .027,  .073,  .206,  .395,  .672,  .860,  .969/
  
data ((cldnuctab(iw,iconc,7),iw=1,9),iconc=1,7)/                 &
  .098,  .206,  .365,  .614,  .841,  .969,  .990,  .990,  .990,  &
  .054,  .128,  .280,  .489,  .727,  .906,  .990,  .990,  .990,  &
  .038,  .085,  .184,  .365,  .583,  .822,  .957,  .990,  .990,  &
  .022,  .063,  .128,  .280,  .457,  .701,  .892,  .982,  .990,  &
  .012,  .038,  .085,  .184,  .365,  .583,  .822,  .949,  .990,  &
  .005,  .018,  .045,  .128,  .280,  .489,  .727,  .892,  .979,  &
  .002,  .007,  .022,  .063,  .164,  .365,  .614,  .822,  .949/

if (jnmb(1) == 1 .or. jnmb(1) == 4) then

! cloud number per kg_air specified in parm(1)   

   rnuc = parm(1) * emb0(1)   ! kg_nuc / kg_air

   do k = lpw0,mza0
      excessrhov = rhov(k) - 1.0001 * rhovslair(k)  ! up by rhoa factor
      rhocnew = 0.                                  ! up by rhoa factor
      if (excessrhov > 0.) then                     ! up by rhoa factor
         rhocnew = min(rnuc*real(rhoa(k)),.5*excessrhov)  ! up by rhoa factor
         rx(k,1) = rx(k,1) + rhocnew                ! up by rhoa factor
         rhov(k) = rhov(k) - rhocnew                ! up by rhoa factor
         cx(k,1) = min(parm(1),rx(k,1) / emb0(1))   ! up by rhoa factor
      endif
   enddo

elseif (jnmb(1) >= 5) then

! cloud number predicted from ccn field

   do k = lpw0,mza0
      excessrhov = rhov(k) - 1.0001 * rhovslair(k)      ! up by rhoa factor

      if (excessrhov > 0.) then             ! up by rhoa factor

          tairc_nuc = tairc(k)
          if (tairc_nuc < -30.) then
             tairc_nuc = -30.
          elseif (tairc_nuc > 30.) then
             tairc_nuc = 30.
          endif
          jtemp = nint(.1 * (tairc_nuc + 30.)) + 1

          w_nuc = wc0(k)
          if (w_nuc < .010001) then 
             w_nuc = .010001
          elseif (w_nuc > 99.99) then
             w_nuc = 99.99
          endif
          rjw = 2. * log10(100. * w_nuc) + 1.
          jw = int(rjw)
          wtw2 = rjw - float(jw)
          wtw1 = 1. - wtw2

          if (jnmb(1) == 5) then
             concen_nuc = cparm       ! #/kg specified in cparm 
          elseif (jnmb(1) == 6) then
           ! concen_nuc = prof(k)     ! #/kg specified in prof(k)
          elseif (jnmb(1) == 7) then
             concen_nuc = cccnx(k)    ! #/kg copied from con_ccn(k,iw0)
          else
             write(io6,*) 'icloud set to value greater than 7: ',icloud
             write(io6,*) 'stopping model '
             stop 'icloud'
          endif  

          if (concen_nuc < 10.001e6) then
             concen_nuc = 10.001e6
          elseif (concen_nuc > 9999.e6) then
             concen_nuc = 9999.e6
          endif
          
          rjconcen = 2. * log10(1.e-7 * concen_nuc) + 1.
          jconcen = int(rjconcen)
          wtconcen2 = rjconcen - float(jconcen)
          wtconcen1 = 1. - wtconcen2

          tab = wtconcen1 * (wtw1 * cldnuctab(jw  ,jconcen  ,jtemp)   &
                          +  wtw2 * cldnuctab(jw+1,jconcen  ,jtemp))  &
              + wtconcen2 * (wtw1 * cldnuctab(jw  ,jconcen+1,jtemp)   &
                          +  wtw2 * cldnuctab(jw+1,jconcen+1,jtemp))

          concen_tab = concen_nuc * tab * rhoa(k)   ! #/m^3 (up by rhoa)
          
! Nucleate cloud droplets only if concen_tab > existing cloud concentration

          if (concen_tab > cx(k,1)) then
             cxadd = concen_tab - cx(k,1)           ! up by rhoa factor
             if (cxadd > excessrhov / emb0(1)) cxadd = excessrhov / emb0(1) ! up by rhoa
             cx(k,1) = cx(k,1) + cxadd              ! up by rhoa factor
             rx(k,1) = rx(k,1) + excessrhov         ! up by rhoa factor
          endif

      endif

   enddo

else
   write(io6,*) 'icloud not allowed to be 2 or 3'
   write(io6,*) 'stopping model '
   stop 'icloud'
endif

return
end subroutine cldnuc

!===============================================================================

subroutine icenuc(k1,k2,lpw0,mrl0  &
   ,jhcat,rx,cx,qr,emb,vap,tx,rhov,rhoa,press0,dynvisc,thrmcon  &
   ,tair,tairc,rhovslair,rhovsiair,cifnx,dtl0)

use micro_coms,  only: mza0, ncat, dnfac, pwmasi, rxmin, ndnc, ddnc, dtc,  &
                       fracc, drhhz, dthz, frachz, ipris, emb0
                      
use consts_coms, only: cice
use misc_coms,   only: io6

implicit none

integer, intent(inout) :: k1(10)
integer, intent(inout) :: k2(10)

integer, intent(in) :: lpw0
integer, intent(in) :: mrl0

integer, intent(inout) :: jhcat(mza0,ncat)

real, intent(inout) :: rx  (mza0,ncat)
real, intent(inout) :: cx  (mza0,ncat)
real, intent(inout) :: qr  (mza0,ncat)
real, intent(in   ) :: emb (mza0,ncat)
real, intent(in   ) :: vap (mza0,ncat)
real, intent(in   ) :: tx  (mza0,ncat)

real, intent(in) :: rhov     (mza0)
real, intent(in) :: press0   (mza0)
real, intent(in) :: dynvisc  (mza0)
real, intent(in) :: thrmcon  (mza0)
real, intent(in) :: tair     (mza0)
real, intent(in) :: tairc    (mza0)
real, intent(in) :: rhovslair(mza0)
real, intent(in) :: rhovsiair(mza0)
real, intent(in) :: cifnx    (mza0)

real(kind=8), intent(in) :: rhoa(mza0)

real, intent(in) :: dtl0

integer :: k,idnc,itc,irhhz,ithz,nt,ns

real :: dn1,fraccld,ridnc,ssi0,wdnc2,tc,ritc,wtc2  &
       ,pbvi,ptvi,pdvi,ptotvi,fracifn,cldnuc,cldnucr,rhhz,haznuc  &
       ,rirhhz,wrhhz2,thz,rithz,wthz2,frachaz,ssi,diagni  &
       ,vapnuc,vapnucr,availvap,relhum

! Define ssi0 to be maximum supersaturation with respect to ice for
! determining total number of IFN that can nucleate in Meyers' formula
data ssi0/0.40/
save

! implement paul's immersion freezing of rain here.  This would
! replace mike's homogeneous freezing of rain which was in h03.

do k = k1(1),k2(1)

!  Homogeneous ice nucleation of cloud droplets

! define dn locally from emb

   dn1 = dnfac(1) * emb(k,1) ** pwmasi(1)

   fraccld = 0.

   if (rx(k,1) > rxmin(1) .and. tairc(k) <= -30.01) then

      ridnc = max(1.,min(real(ndnc-1),dn1 / ddnc))
      idnc = int(ridnc)
      wdnc2 = ridnc - float(idnc)

      tc = max(-49.99,tairc(k))
      ritc = (tc + 50.00) / dtc + 1.0
      itc = int(ritc)
      wtc2 = ritc - float(itc)
      fraccld = (1.-wdnc2) * (1.-wtc2) * fracc(idnc  ,itc  ,mrl0)  &
              +     wdnc2  * (1.-wtc2) * fracc(idnc+1,itc  ,mrl0)  &
              + (1.-wdnc2) *     wtc2  * fracc(idnc  ,itc+1,mrl0)  &
              +     wdnc2  *     wtc2  * fracc(idnc+1,itc+1,mrl0)

   endif

!  Heterogeneous contact ice nucleation of cloud droplets by diffusio-
!  phoresis, thermophoresis, and Brownian motion (transport of IN)

   call contnuc (rx(k,1),cx(k,1),tx(k,1),vap(k,1),press0(k)  &
      ,dynvisc(k),thrmcon(k),tair(k),tairc(k)  &
      ,pbvi,ptvi,pdvi,ptotvi,dn1,dtl0,rxmin(1))

! progIFN: Scale ptotvi returned from contnuc by prognosed IFN fraction

!::later   ptotvi = ptotvi * fracifn

! MIKE ADDED THIS COMMENTED ccinp(k)=ccinp(k)-ptotvi, but
! probably do not want sink of ccinp here.

   cldnuc = ptotvi + max(0.,fraccld * cx(k,1) - cx(k,3))  ! up by rhoa factor
!         cldnucr = cldnuc * emb(k,1)
   cldnucr = min(rx(k,1),ptotvi * emb(k,1) + fraccld * rx(k,1))  ! up by rhoa factor

   rx(k,3) = rx(k,3) + cldnucr  ! up by rhoa factor
   rx(k,1) = rx(k,1) - cldnucr  ! up by rhoa factor
   cx(k,3) = cx(k,3) + cldnuc   ! up by rhoa factor
   cx(k,1) = cx(k,1) - cldnuc   ! up by rhoa factor

enddo

! DEMOTT'S NEW SCHEME: In 4.3 and beyond, assume that it gives #/KG

!  Homogeneous nucleation of haze

do k = lpw0,mza0
   rhhz = rhov(k) / rhovslair(k)  ! stays the same
   haznuc = 0.
   if (rhhz > 0.82 .and. tairc(k) <= -35.01) then
      rirhhz = min(0.1799,rhhz-0.82) / drhhz + 1.0
      irhhz = int(rirhhz)
      wrhhz2 = rirhhz - float(irhhz)
      thz = max(-59.99,tairc(k))
      rithz = (thz + 60.00) / dthz + 1.0
      ithz = int(rithz)
      wthz2 = rithz - float(ithz)
      frachaz = (1.-wrhhz2) * (1.-wthz2) * frachz(irhhz  ,ithz  )  &
              +     wrhhz2  * (1.-wthz2) * frachz(irhhz+1,ithz  )  &
              + (1.-wrhhz2) *     wthz2  * frachz(irhhz  ,ithz+1)  &
              +     wrhhz2  *     wthz2  * frachz(irhhz+1,ithz+1)
      frachaz = 1. - exp(-frachaz * dtl0)
! OPTION 1
      haznuc = frachaz * 300.e6
! OPTION 2
!           haznuc = frachaz * caero(k)
   endif

! meyers -  no cld aerosol source or sink here

!  Heterogeneous nucleation by deposition condensation freezing
!  with deposition nuclei.  In 4.3 and beyond, assume that it gives #/kg.

   ssi = min(ssi0,rhov(k) / rhovsiair(k) - 1.)

   if (ssi > 0. .and. tairc(k) <= -5.) then
      fracifn = exp(12.96 * (ssi - ssi0))
   else
      fracifn = 0.
   endif

! Diagnose maximum number of IFN to activate based on ipris

   if (ipris == 5) then
      diagni = fracifn * 1.e5
   elseif (ipris == 6) then
      diagni = fracifn * rhoa(k) ** 5.4 * 1.e5
   elseif (ipris == 7) then
      diagni = fracifn * cifnx(k)
   endif

! orig Meyers formula:     +      diagni = exp(6.269 + 12.96 * ssi)

!  Combine nucleation types, and limit amounts
! vapnuc is #/kg_air and vapnucr is kg/kg_air

! BEGIN MIKE'S SECTION FOR LIMITING NUMBER OF CRYSTALS NUCLEATED
! BY NUMBER OF ICE CRYSTALS PRESENT ALREADY

   vapnuc = max(0.,(haznuc + diagni) * real(rhoa(k)) - cx(k,3))   ! incld rhoa
   vapnucr = vapnuc * emb0(3)
   if (vapnucr > 0.) then
      availvap = .5 * (rhov(k) - rhovsiair(k))
      if (vapnucr > availvap) then
         vapnucr = min(vapnucr, max(0.,availvap))
      endif
   endif
   vapnuc = vapnucr / emb0(3)

   rx(k,3) = rx(k,3) + vapnucr
   qr(k,3) = qr(k,3) + vapnucr * cice * tairc(k)
   cx(k,3) = cx(k,3) + vapnuc

enddo

! here mike has the habit diagnosis. option 1 is to use habit
! at cloud top, option 2 is to use new habit at each level.
! need to consider other options.  how about method of formation?
! my question about how much of habit is due to existing ice
! structure, and how much is due to current growth environment
! (temp and supsat). relevant supsat is wrt liquid?

return
end subroutine icenuc

!===============================================================================

subroutine contnuc (rx,cx,tx,vap,press  &
   ,dynvisc,thrmcon,tair,tairc,pbvi,ptvi,pdvi,ptotvi,dn1,dtl,rxmin1)

implicit none

real, intent(in) :: rx,cx,tx,vap,press,dynvisc,thrmcon,tair,tairc
real, intent(in) :: dn1,dtl,rxmin1
real, intent(out) :: pbvi,ptvi,pdvi,ptotvi

real :: aka,raros,ana,akn,dfar,f1,f2,ft
data aka,raros/5.39e-3,3.e-7/

!  Heterogeneous contact ice nucleation of cloud droplets by diffusio-
!  phoresis, thermophoresis, and Brownian motion (transport of IN)
!
!  ana   = # IN per kg available for contact freezing (from Meyers et al. 1992
!          where ana was interpreted as # per m^3)
!  akn   = Knudsen number (Walko et al. 1995, Eq. 58)
!          [2.28e-5 = mfp * p00 / 293.15]
!  raros = aerosol radius = 3.e-7 m from Cotton et al. (1986)
!  dfar  = aerosol diffusivity (Pruppacher and Klett Eq. 12-15)
!          [7.32e-25 = Boltzmann constant / (6 pi)]
!  f1    = "function 1" (Walko et al. 1995 Eq. 55) multiplied by delta t
!          [#/m^3]
!  f2    = "function 2" (Walko et al. 1995 Eq. 56)
!  ft    = "function ft" (Walko et al. 1995 Eq. 57)
!  pbvi  = Brownian motion nucleation amount this timestep [#/m^3]
!  ptvi  = Thermophoretic nucleation amount this timestep [#/m^3]
!  pdvi  = Diffusiophoretic nucleation amount this timestep [#/m^3],
!          reformulated to use vapor diffusion directly.  Factor of 1.2
!          is (1+sigma_va x_a) from Pruppacher and Klett Eq. 12-102
!          divided by .622, the molecular weight ratio between water and air.

ptotvi = 0.

if (tx <= -2. .and. rx > rxmin1) then

   ana = exp(4.11 - .262 * tx)
   akn = 2.28e-5 * tair / (press * raros)
   dfar = 7.32e-25 * tair * (1.+ akn) / (raros * dynvisc)
   f1 = 6.28318 * dn1 * cx * ana * dtl  ! up by rhoa
   f2 = thrmcon * (tairc - tx) / press
   ft = .4 * (1. + 1.45 * akn + .4 * akn * exp(-1. / akn))  &
      * (thrmcon + 2.5 * akn * aka)  &
      / ((1. + 3. * akn)  &
      * (2. * thrmcon + 5. * aka * akn + aka))
   pbvi = f1 * dfar                     ! up by rhoa
   ptvi = f1 * f2 * ft                  ! up by rhoa
   pdvi = 1.2 * ana * vap               ! up by rhoa
   ptotvi = max(0.,pbvi + ptvi + pdvi)  ! up by rhoa

endif
return
end subroutine contnuc
