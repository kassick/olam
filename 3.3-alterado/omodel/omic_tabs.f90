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
subroutine haznuc()

use micro_coms, only: nthz, dthz, nrhhz, drhhz, frachz
use misc_coms,  only: io6

implicit none

integer :: ithz,irhhz,k
real :: denccn,gnuccn,dnccn,ddccn,rhhz,c1hz,c2hz,c3hz,bhz,dm,sum  &
       ,dccn,y,dum,thz
real :: gammln

!  Haze nucleation table

denccn = 1.769
gnuccn = 1.
dnccn =   .075E-4
ddccn = .005e-4
do ithz = 1,nthz
   thz = -60. + dthz * float(ithz - 1)
   do irhhz = 1,nrhhz
      rhhz = 0.82 + drhhz * float(irhhz - 1)
      c1hz = (3.14159 * denccn / 6.) ** (-.333333)
      c2hz = -14.65 - 1.045 * thz
      c3hz = -492.35 - 8.34 * thz - 0.0608 * thz ** 2
      bhz = min(38., max(-38., c2hz + c3hz * (1. - rhhz)))
      dm = c1hz * 10 ** (-bhz/6.)

      sum = 0.
      dccn = 0.
      do k=1,200
         dccn = dccn + ddccn
         y=dccn / dnccn
         dum=min(50., (dccn / dm) ** 6)
         sum = sum + y ** (gnuccn-1.) * exp(-y) * (1. - exp(-dum))
      enddo
      frachz(irhhz,ithz) = sum*ddccn/(exp(gammln(gnuccn))*dnccn)
   enddo
enddo

return
end subroutine haznuc

!===============================================================================

subroutine homfrzcl(dtl,mrl)

use micro_coms, only: ntc, dtc, ndnc, ddnc, fracc
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mrl
real, intent(in) :: dtl

integer :: itc,k,idnc
real :: gnuc,ddc,ajlso,dnc,sum,dc,v1,tc,y
real :: gammln

!  Make table for homogeneous freezing of cloud droplets
!  Need gnuc = gnu(1) ???

gnuc = 1.
ddc = 0.5e-6
do itc = 1,ntc
   tc = -50. + dtc * float(itc-1)
   y = -(606.3952+tc*(52.6611+tc*(1.7439+tc*(.0265+tc*1.536e-4))))
   ajlso = 1.e6 * 10. ** y
   do idnc = 1,ndnc
      dnc = ddnc * float(idnc)
      sum = 0.
      dc = 0.
      do k = 1,2000
         dc = dc + ddc
         v1 = 0.523599 * dc ** 3
         sum = sum + (dc / dnc) ** (gnuc - 1.) * exp(-dc / dnc)  &
            * (1. - exp(-ajlso * v1 * dtl))
      enddo
      fracc(idnc,itc,mrl) = sum * ddc / (exp(gammln(gnuc)) * dnc)
   enddo
enddo
return
end subroutine homfrzcl

!===============================================================================

subroutine mksedim_tab(m1,zm,dzt,dzit)

use micro_coms, only: sedtime0, sedtime1, mza0, zmf, dztf, dzitf, maxkfall,  &
                      nhcat, lcat_lhcat, dispemb0, cfvt, emb0, cfmasi,  &
                      pwvt, pwmasi, dispemb0i, dispemb1, ch2, nembfall, gnu,  &
                      pwmas, pcpfillc, pcpfillr, emb1
use misc_coms,  only: io6

implicit none

integer, intent(in) :: m1
real, intent(in) :: zm(m1),dzt(m1),dzit(m1)

integer, parameter :: nbin=50
integer :: iembs,lcat,lhcat,k,kkf,ibin,kk,jbin

real :: dmbodn,diam0,diam1,fac1,fac3,sumc,sumr,diam,fac2,fac4  &
       ,disp,ztopnew,zbotnew,fallin,delzsfc,dispemb,dispmax,dispmx
real :: gammln,gammp
real, dimension(nbin) :: cbin,rbin,reldisp

! To cover the possible range over all timesteps, define sedtime0 and sedtime1
! here as 0.1 seconds and 3000 seconds.  The former is supposed to be
! less than 0.7 of the shortest timestep on any grid (sqrt(rhoi) never exceeds
! 0.7) and the latter is the longest timestep expected to ever be used (600
! seconds) times a factor of 5 for the largest value of sqrt(rhoi).

sedtime0 = .1
sedtime1 = 3000.
dispmax = 500.

! Fill zmf and ztf arrays to easily handle "underground" points

do k = 1,mza0
   zmf(k) = zm(k)
   dztf(k) = dzt(k)
   dzitf(k) = dzit(k)
enddo

do k = 0,2-maxkfall,-1
   zmf(k) = 2. * zmf(k+1) - zmf(k+2)
   dztf(k+1) = zmf(k+1) - zmf(k)
   dzitf(k+1) = 1. / dztf(k+1)
enddo

! Loop over hydrometeor categories

do lhcat = 1,nhcat
   lcat = lcat_lhcat(lhcat)

   dispemb0(lhcat) = sedtime0 * cfvt(lhcat)  &
      * (emb0(lcat) * cfmasi(lhcat)) ** (pwvt(lhcat) * pwmasi(lhcat))
      
   dispemb0i(lhcat) = 1. / dispemb0(lhcat)

   dispemb1(lhcat) = sedtime1 * cfvt(lhcat)  &
      * (emb1(lcat) * cfmasi(lhcat)) ** (pwvt(lhcat) * pwmasi(lhcat))

!Bob (10/24/00):  Limit dispemb1 to a maximum of dispmax

   if (dispemb1(lhcat) > dispmax) dispemb1(lhcat) = dispmax

   ch2(lhcat) = float(nembfall-1) / log10(dispemb1(lhcat) * dispemb0i(lhcat))

! Loop over bins, filling them with fractional number, fractional mass,
! and displacement quotient relative to emb.

   dmbodn = (exp(gammln(gnu(lcat) + pwmas(lhcat))  &
      - gammln(gnu(lcat)))) ** pwmasi(lhcat)
   diam0 = 0.06 * dmbodn
   diam1 = 1.0 * dmbodn
   fac1 = gammp(gnu(lcat),diam0)
   fac3 = gammp(gnu(lcat) + pwmas(lhcat),diam0)
   sumc = 0.
   sumr = 0.

   do jbin = 1,nbin

      diam = diam0 * (diam1 / diam0) ** (float(jbin)/float(nbin))
      fac2 = gammp(gnu(lcat),diam)
      fac4 = gammp(gnu(lcat) + pwmas(lhcat),diam)
      cbin(jbin) = fac2 - fac1
      rbin(jbin) = fac4 - fac3
      fac1 = fac2
      fac3 = fac4
      sumc = sumc + cbin(jbin)
      sumr = sumr + rbin(jbin)
      reldisp(jbin) = diam ** pwvt(lhcat)

   enddo

   do jbin = 1,nbin
      cbin(jbin) = cbin(jbin) / sumc
      rbin(jbin) = rbin(jbin) / sumr
   enddo

! Loop over displacement distance for size emb.

   do iembs = 1,nembfall
      dispemb = dispemb0(lhcat)  &
         * (dispemb1(lhcat) * dispemb0i(lhcat))  &
         ** (float(iembs-1) / float(nembfall-1))

! Zero out concentration and mass fill arrays and surface precip array
! before accumulation.

      do k = 1,mza0
         do kkf = 1,maxkfall
            pcpfillc(k,kkf,iembs,lhcat) = 0.
            pcpfillr(k,kkf,iembs,lhcat) = 0.
         enddo
      enddo

! Loop over vertical grid index.

      do k = 2,mza0

!Bob (10/24/00):  Limit disp distance to (maxkfall-1) levels

         dispmx = min(dispmax,zmf(k-1) - zmf(k-maxkfall))

! Loop over bins

         do ibin = 1,nbin
            disp = dispemb * reldisp(ibin)
            if (disp > dispmx) disp = dispmx

            ztopnew = zmf(k) - disp
            zbotnew = zmf(k-1) - disp

! Loop over grid cells that a parcel falls into, including the one it starts from.

            do kkf = 1,maxkfall

               kk = k + 1 - kkf
               if (zbotnew > zmf(kk)) go to 50

               if (ztopnew <= zmf(kk-1)) then
                  fallin = 0.
               else
                  fallin = dzitf(kk) *  &
                     (min(zmf(kk),ztopnew) - max(zmf(kk-1),zbotnew))
               endif

               pcpfillc(k,kkf,iembs,lhcat) = pcpfillc(k,kkf,iembs,lhcat)  &
                  + fallin * cbin(ibin)

               pcpfillr(k,kkf,iembs,lhcat) = pcpfillr(k,kkf,iembs,lhcat)  &
                  + fallin * rbin(ibin)

            enddo

50          continue

         enddo

      enddo
   enddo
enddo

return
end subroutine mksedim_tab

!===============================================================================

subroutine tabmelt()

use micro_coms, only: nhcat, lcat_lhcat, gnu, rmlttab, enmlttab, ndns,  &
                      shedtab, cfmas, pwmas, cfvt, pwvt, shapefac, ninc, ncat
use misc_coms,  only: io6

implicit none

integer, parameter :: nbins=500

integer :: lhcat,lcat,ndns1,ibin,inc,iter,idns
real :: dn,gammaa,totfmg,totmass,vtx,fre,totqm,qmgoal,qmnow,totmdqdt,deltat  &
       ,pliqmass,picemass,critmass,vk

real, dimension(nbins) :: db,fmg,pmass,binmass,dqdt,q
real, dimension(ncat) :: dmean
real :: gammln
data dmean/20.e-6,500.e-6,30.e-6,500.e-6,500.e-6,500.e-6,8000.e-6/
data vk/0.2123e-04/

do lhcat = 1,nhcat
   lcat = lcat_lhcat(lhcat)

   dn = dmean(lcat) / gnu(lcat)
   gammaa = exp(gammln(gnu(lcat)))

   rmlttab(1) = 0.0
   rmlttab(ninc) = 1.0
   enmlttab(1,lhcat) = 0.0
   enmlttab(ninc,lhcat) = 1.0

   ndns1 = 1
   if (lcat == 7) ndns1 = ndns

   do idns = 1,ndns1
      shedtab(1,idns) = 0.0
      shedtab(ninc,idns) = 0.0

      if (ndns1 > 1) dn = 1.e-3 * float(idns) / gnu(lcat)

      totfmg = 0.
      totmass = 0.
      do ibin = 1,nbins
         db(ibin) = 0.02 * dn * (float(ibin) - 0.5)
         fmg(ibin) = (db(ibin) / dn) ** (gnu(lcat) - 1.)  &
            / (dn * gammaa) * exp(-db(ibin) / dn)
         totfmg = totfmg + fmg(ibin)
         q(ibin) = 0.
         pmass(ibin) = cfmas(lhcat) * db(ibin) ** pwmas(lhcat)
         binmass(ibin) = pmass(ibin) * fmg(ibin)
         totmass = totmass + binmass(ibin)
         vtx = cfvt(lhcat) * db(ibin) ** pwvt(lhcat)
         fre = (1.0 + 0.229 * sqrt(vtx * db(ibin) / vk))  &
            * shapefac(lhcat)
         dqdt(ibin) = db(ibin) ** (1. - pwmas(lhcat)) * fre
      enddo
      totqm = totmass * 80.

      do inc = 2,ninc-1
         qmgoal = totqm * float(inc-1) / float(ninc-1)

         do iter = 1,2
            qmnow = 0.
            totmdqdt = 0.
            do ibin = 1,nbins
               if(q(ibin) < 79.9999)then
                  totmdqdt = totmdqdt + binmass(ibin) * dqdt(ibin)
               endif
               qmnow = qmnow + q(ibin) * binmass(ibin)
            enddo
            deltat = max(0.,(qmgoal - qmnow) / totmdqdt)
            do ibin = 1,nbins
               q(ibin) = min(80.,q(ibin) + dqdt(ibin) * deltat)
            enddo
         enddo

!  For the current inc value (representing total liquid fraction), compute
!  melted mixing ratio (rmlttab) and number (enmlttab) from totally-melted
!  bins and compute shedded mixing ratio (shedtab) from partially-melted bins.

         if(idns == 7)then
            rmlttab(inc) = 0.0
            do ibin = 1,nbins
               if(q(ibin) > 79.9)then
                  rmlttab(inc) = rmlttab(inc) + binmass(ibin)
               endif
            enddo
            rmlttab(inc) = rmlttab(inc) / totmass
         endif

         if(idns == 7 .or. ndns1 == 1)then
            enmlttab(inc,lhcat) = 0.0
            do ibin = 1,nbins
               if(q(ibin) > 79.9)then
                  enmlttab(inc,lhcat) = enmlttab(inc,lhcat) + fmg(ibin)
               endif
            enddo
            enmlttab(inc,lhcat) = enmlttab(inc,lhcat) / totfmg
         endif

         if(lcat == 7)then
            shedtab(inc,idns) = 0.0
!                  do ibin = kbin,nbins
            do ibin = 1,nbins
               if(q(ibin) <= 79.9)then
                  pliqmass = pmass(ibin) * q(ibin) / 80.
                  picemass = pmass(ibin) - pliqmass
                  critmass = .268e-3 + .1389 * picemass
                  shedtab(inc,idns) = shedtab(inc,idns)  &
                     + max(0.0, pliqmass - critmass) * fmg(ibin)
               endif
            enddo
            shedtab(inc,idns) = shedtab(inc,idns) / totmass
         endif

      enddo
   enddo
enddo
return
end subroutine tabmelt

!===============================================================================

subroutine mkcoltb()

use micro_coms, only: nhcat, lcat_lhcat, gnu, pwmas, emb0, cfmasi, pwmasi,  &
                      emb1, ipairc, ipairr, pwvt, nembc, cfvt, cfmas,  &
                      coltabc, coltabr
use misc_coms,  only: io6

implicit none

integer, parameter :: ndx=20
integer :: ihx,ix,ihy,iy,iemby,iembx,idx
real :: gxm,dnminx,dnmaxx,dxlo,dxhi,gyn,gyn1,gyn2,gynp,gynp1,gynp2,gym  &
       ,dnminy,dnmaxy,dny,vny,dnx,ans
real :: gammln,xj,bans,ratmin,ratmax

real, dimension(ndx) :: dx,fx,gx

ratmin = 1.
ratmax = 1.

do ihx = 1,nhcat
   ix = lcat_lhcat(ihx)

   gxm = exp(gammln(gnu(ix)) - gammln(gnu(ix) + pwmas(ihx)))
   dnminx = ((emb0(ix) * cfmasi(ihx)) * gxm) ** pwmasi(ihx)
   dnmaxx = ((emb1(ix) * cfmasi(ihx)) * gxm) ** pwmasi(ihx)
   dxlo = .01 * dnminx
   dxhi = 10. * dnmaxx

do ihy = 1,nhcat

   write(io6,*) 'ihx,ihy,ipairc,ipairr',ihx,ihy,ipairc(ihx,ihy)  &
      ,ipairr(ihx,ihy)

   iy = lcat_lhcat(ihy)

   if (ipairc(ihx,ihy) > 0 .or. ipairr(ihx,ihy) > 0) then
      gyn = exp(gammln(gnu(iy)))
      gyn1 = exp(gammln(gnu(iy) + 1.)) / gyn
      gyn2 = exp(gammln(gnu(iy) + 2.)) / gyn
      gynp = exp(gammln(gnu(iy) + pwvt(ihy))) / gyn
      gynp1 = exp(gammln(gnu(iy) + pwvt(ihy) + 1.)) / gyn
      gynp2 = exp(gammln(gnu(iy) + pwvt(ihy) + 2.)) / gyn

      gym = exp(gammln(gnu(iy)) - gammln(gnu(iy) + pwmas(ihy)))
      dnminy = ((emb0(iy) * cfmasi(ihy)) * gym) ** pwmasi(ihy)
      dnmaxy = ((emb1(iy) * cfmasi(ihy)) * gym) ** pwmasi(ihy)

      do iemby = 1,nembc
         dny = dnminy * (dnmaxy / dnminy) ** (float(iemby-1)  &
            / float(nembc-1))
         vny = cfvt(ihy) * dny ** pwvt(ihy)
         do iembx = 1,nembc

            dnx = dnminx * (dnmaxx / dnminx) ** (float(iembx-1)  &
               / float(nembc-1))
            do idx = 1,ndx
               dx(idx) = dxlo * (dxhi / dxlo)  &
                  ** (float(idx-1) / float(ndx-1))
                  
               fx(idx) = xj(dx(idx),cfvt(ihx),pwvt(ihx),cfvt(ihy)  &
                  ,pwvt(ihy),vny,dnx,dny,gnu(ix),gnu(iy)  &
                  ,gyn1,gyn2,gynp,gynp1,gynp2)
               gx(idx) = fx(idx) * cfmas(ihx) * dx(idx) ** pwmas(ihx)

            enddo
            if (ipairc(ihx,ihy) > 0) then
               call avint(dx,fx,ndx,dxlo,dxhi,ans)
!nonlog10                     coltabc(iembx,iemby,ipairc(ihx,ihy))=max(0.,ans)
               coltabc(iembx,iemby,ipairc(ihx,ihy))=  &
                  -log10(max(1.e-30,ans))
            endif
            if (ipairr(ihx,ihy) > 0) then
               call avint(dx,gx,ndx,dxlo,dxhi,ans)
!nonlog10                     coltabr(iembx,iemby,ipairr(ihx,ihy))=max(0.,ans)
               coltabr(iembx,iemby,ipairr(ihx,ihy))=  &
                  -log10(max(1.e-30,ans))
            endif

         enddo
      enddo
   endif
enddo
enddo

return
end subroutine mkcoltb

!===============================================================================

subroutine mkcoltb_brute()

use micro_coms, only: nhcat, lcat_lhcat, gnu, pwmas, emb0, cfmasi, pwmasi,  &
                      emb1, ipairc, ipairr, nembc, cfvt, pwvt, cfmas,  &
                      coltabc, coltabr
use misc_coms,  only: io6

implicit none

integer, parameter :: ndx=100,ndy=100
integer :: ihx,ix,ihy,iy,iemby,iembx,idx,idy
real :: gxm,dnminx,dnmaxx,dxlo,dxhi,gxn,gyn,gym  &
       ,dnminy,dnmaxy,dny,dnx,bint,sum_num,sum_xmass,vx,vy,dx,dy  &
       ,fgamx,emx,dx1,dx2,fgamy,emy,dy1,dy2,dyhi,dylo
real, external :: gammln

! Loop over colliding category X

do ihx = 1,nhcat
   ix = lcat_lhcat(ihx)

   gxm = exp(gammln(gnu(ix)) - gammln(gnu(ix) + pwmas(ihx)))
   dnminx = ((emb0(ix) * cfmasi(ihx)) * gxm) ** pwmasi(ihx)
   dnmaxx = ((emb1(ix) * cfmasi(ihx)) * gxm) ** pwmasi(ihx)
   dxlo = .01 * dnminx
   dxhi = 10. * dnmaxx

! Loop over colliding category Y

   do ihy = 1,nhcat

! Check if colliding pair (X,Y) is to be considered

      if (ipairc(ihx,ihy) == 0 .and. ipairr(ihx,ihy) == 0) cycle

      iy = lcat_lhcat(ihy)

      gym = exp(gammln(gnu(iy)) - gammln(gnu(iy) + pwmas(ihy)))
      dnminy = ((emb0(iy) * cfmasi(ihy)) * gym) ** pwmasi(ihy)
      dnmaxy = ((emb1(iy) * cfmasi(ihy)) * gym) ** pwmasi(ihy)
      dylo = .01 * dnminy
      dyhi = 100. * dnmaxy

! Loop over all mean-mass table values for category Y

      do iemby = 1,nembc
         dny = dnminy * (dnmaxy / dnminy) ** (float(iemby-1)  &
            / float(nembc-1))

! Loop over all mean-mass table values for category X

         do iembx = 1,nembc
            dnx = dnminx * (dnmaxx / dnminx) ** (float(iembx-1)  &
               / float(nembc-1))

            sum_num = 0.    ! Initialize integral number sum to zero
            sum_xmass = 0.  ! Initialize integral xmass sum to zero
            
! Loop over spectrum of Y diameters for current mean-mass Y value

            do idy = 1,ndy
               dy1 = dylo * (dyhi / dylo) ** (float(idy-1) / float(ndy))
               dy2 = dylo * (dyhi / dylo) ** (float(idy) / float(ndy))
               dy = .5 * (dy1 + dy2)
               vy = cfvt(ihy) * dy ** pwvt(ihy)
               emy = cfmas(ihy) * dy ** pwmas(ihy)  ! not used

               gyn = exp(gammln(gnu(iy)))
               fgamy = (dy / dny) ** (gnu(iy) - 1) * exp(-dy / dny)  &
                  / (gyn * dny)

! Loop over spectrum of X diameters for current mean-mass X value

               do idx = 1,ndx
                  dx1 = dxlo * (dxhi / dxlo) ** (float(idx-1) / float(ndx))
                  dx2 = dxlo * (dxhi / dxlo) ** (float(idx) / float(ndx))
                  dx = .5 * (dx1 + dx2)
                  vx = cfvt(ihx) * dx ** pwvt(ihx)
                  emx = cfmas(ihx) * dx ** pwmas(ihx)

                  gxn = exp(gammln(gnu(ix)))
                  fgamx = (dx / dnx) ** (gnu(ix) - 1) * exp(-dx / dnx)  &
                     / (gxn * dnx)

! BINT is (integrand * del_dx * del_dy)

                  bint = (dx + dy) ** 2 * abs(vx - vy) * fgamx * fgamy  &
                     * (dy2 - dy1) * (dx2 - dx1)

                  sum_num = sum_num + bint
                  sum_xmass = sum_xmass + bint * emx
               enddo
            enddo

! sum_num and sum_xmass are the definite integral sums of number and mass
! for the current X and Y mean-mass diameter pair.  Enter these in tables.

            if (ipairc(ihx,ihy) > 0) then
!nonlog10      coltabc(iembx,iemby,ipairc(ihx,ihy)) = sum_num
               coltabc(iembx,iemby,ipairc(ihx,ihy)) =  &
                  -log10(max(1.e-30,sum_num))
            endif

            if (ipairr(ihx,ihy) > 0) then
!nonlog10      coltabr(iembx,iemby,ipairr(ihx,ihy)) = sum_xmass
               coltabr(iembx,iemby,ipairr(ihx,ihy)) =  &
                  -log10(max(1.e-30,sum_xmass))
            endif

         enddo
      enddo
   enddo
enddo

return
end subroutine mkcoltb_brute

!===============================================================================

real function xj(dx,cvx,pvx,cvy,pvy,vny,dnx,dny,xnu,ynu  &
                ,gyn1,gyn2,gynp,gynp1,gynp2)
implicit none
real :: dx,cvx,pvx,cvy,pvy,vny,dnx,dny,xnu,ynu,gyn1,gyn2,gynp,gynp1,gynp2  &
       ,dnxi,rdx,vx,dxy,ynup
real :: gammln,gammp,gammq
dnxi = 1. / dnx
rdx = dx * dnxi
vx = cvx * dx ** pvx
dxy = (vx / cvy) ** (1. / pvy) / dny
ynup = ynu + pvy

if (rdx < 38.) then
   xj=exp(-rdx-gammln(xnu)-gammln(ynu))*rdx**(xnu-1.)*dnxi*(  &
       vx*(dx*dx*(gammp(ynu,dxy)-gammq(ynu,dxy))  &
         +2.*dx*dny*gyn1*(gammp(ynu+1.,dxy)-gammq(ynu+1.,dxy))  &
         +dny*dny*gyn2*(gammp(ynu+2.,dxy)-gammq(ynu+2.,dxy)))  &
     -vny*(dx*dx*gynp*(gammp(ynup,dxy)-gammq(ynup,dxy))  &
         +2.*dx*dny*gynp1*(gammp(ynup+1.,dxy)-gammq(ynup+1.,dxy))  &
         +dny*dny*gynp2*(gammp(ynup+2.,dxy)-gammq(ynup+2.,dxy))))
else
   xj = 0.
endif
return
end

!===============================================================================

subroutine make_autotab()

use micro_coms, only: d1min, d1max, d1ecc, nd1cc, d1ecr, nd1cr, r2min,  &
                      r2max, r2ecr, nr2cr, r2err, nr2rr, d2min, d2max,  &
                      gnu, r1tabcc, c1tabcc, c2tabcc, r1tabcr, c1tabcr,  &
                      nd2cr, nd2rr, c2tabrr
use misc_coms,  only: io6

implicit none

integer, parameter :: ibins=36,ithresh=15
integer :: i,k,id1cc,id1cr,ir2cr,id2cr,ir2rr,id2rr
real :: r2,en1,en2,en1i,en1i2,d1,r1,d2minx,d2ecr,d2,sum1,sum10,sun10,sun1  &
       ,sun20,sum20,sun2,sum2,d2err

real, dimension(ibins+1) :: x,diam
real, dimension(ibins) :: ank0,amk0,ank,amk,ank1,amk1,ank2,amk2
real, dimension(ibins,ibins,3) :: akbarx
real, dimension(ibins,ibins) :: akbar

! This subroutine works in cgs units.

! read in mass grid x(k+1)=2*x(k), diameters (diam) and collection kernel kbar

! GNF: kernels for ice with cloud?

!     open(53,file='kbarn',status='old')
!      open(53,file='kbarf',status='old')
call data(x,diam,akbar,ibins)
!      close(53)

akbarx(:,:,:) = 0.0
do i=1,ibins
   do k=1,ibins
      if(i > ithresh .and. k > ithresh) then
         akbarx(i,k,3) = akbar(i,k)
      elseif(i <= ithresh .and. k <= ithresh) then
         akbarx(i,k,1) = akbar(i,k)
      else
         akbarx(i,k,2) = akbar(i,k)
      endif
   enddo
enddo

! d1min and d1max are equivalent to dmb0 and dmb1, but may have different
! values.

d1min = 4.e-4
d1max = 50.e-4
d1ecc = log10 (d1max / d1min) / float(nd1cc-1)
d1ecr = log10 (d1max / d1min) / float(nd1cr-1)

r2min = .01e-6
r2max = 20.e-6
r2ecr = log10 (r2max / r2min) / float(nr2cr-1)
r2err = log10 (r2max / r2min) / float(nr2rr-1)

d2min = 1.e-2
d2max = 1.

! Start 1 cc loop for dm1, dn1, and dn2.

r2 = .01e-06
en1 = 100.
en2 = 1.e-6
en1i = 1. / en1
en1i2 = en1i ** 2

do id1cc = 1,nd1cc

   d1 = d1min + (d1max - d1min) * float(id1cc-1) / float(nd1cc-1)
   r1 = en1 * .5236 * d1 ** 3

   call initg(r1,r2,en1,en2,gnu(1),gnu(2),diam,x,amk0,ank0  &
      ,ank1,amk1,ank2,amk2,ibins,ithresh)
   call sumn(ank0,amk0,1,ithresh,ibins,sun10,sum10)
   call sumn(ank0,amk0,ithresh+1,ibins,ibins,sun20,sum20)
   call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,1))
   call sumn(ank,amk,1,ithresh,ibins,sun1,sum1)
   call sumn(ank,amk,ithresh+1,ibins,ibins,sun2,sum2)

   r1tabcc(id1cc) = max(0.,(sum10-sum1) * en1i2)  ! loss of cloud mass (> 0)
   c1tabcc(id1cc) = max(0.,(sun10-sun1) * en1i2)  ! loss of cloud # (> 0)
   c2tabcc(id1cc) = max(0.,(sun2-sun20) * en1i2)  ! gain of rain # (> 0)

enddo

! Start 3 cr loops for dm1 and dn1.

do id1cr = 1,nd1cr
   d1 = d1min * 10. ** (d1ecr * float(id1cr-1))
   r1 = en1 * .5236 * d1 ** 3

   do ir2cr = 1,nr2cr
      r2 = r2min * 10. ** (r2ecr * float(ir2cr-1))
      d2minx = max(d2min,(r2 / (.1 * .5236)) ** .333333)
      d2ecr = alog10(d2max / d2minx) / float(nd2cr-1)

      do id2cr = 1,nd2cr
         d2 = d2minx * 10. ** (d2ecr * float(id2cr-1))
         en2 = r2 / (.5236 * d2 ** 3)

         call initg(r1,r2,en1,en2,gnu(1),gnu(2),diam,x,amk0,ank0  &
            ,ank1,amk1,ank2,amk2,ibins,ithresh)
         call sumn(ank0,amk0,1,ithresh,ibins,sun10,sum10)
         call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,2))
         call sumn(ank,amk,1,ithresh,ibins,sun1,sum1)

         r1tabcr(id1cr,ir2cr,id2cr)  &
            = alog10(max(1.e-20,(sum10-sum1) * en1i))  ! loss of cloud mass (> 0)
         c1tabcr(id1cr,ir2cr,id2cr)  &
            = alog10(max(1.e-20,(sun10-sun1) * en1i))  ! loss of cloud # (> 0)

      enddo
   enddo
enddo

! Start 2 rr loops for dn2.

d1 = 4.e-4
r1 = en1 * .5236 * d1 ** 3

do ir2rr = 1,nr2rr
   r2 = r2min * 10. ** (r2err * float(ir2rr-1))
   d2minx = max(d2min,(r2 / (.1 * .5236)) ** .333333)
   d2err = alog10(d2max / d2minx) / float(nd2rr-1)

   do id2rr = 1,nd2rr
      d2 = d2minx * 10. ** (d2err * float(id2rr-1))
      en2 = r2 / (.5236 * d2 ** 3)

      call initg(r1,r2,en1,en2,gnu(1),gnu(2),diam,x,amk0,ank0  &
         ,ank1,amk1,ank2,amk2,ibins,ithresh)
      call sumn(ank0,amk0,ithresh+1,ibins,ibins,sun20,sum20)
      call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,3))
      call sumn(ank,amk,ithresh+1,ibins,ibins,sun2,sum2)

      c2tabrr(ir2rr,id2rr) = alog10(max(1.e-25,sun20-sun2))  ! loss of rain # (> 0)

   enddo
enddo
return
end subroutine make_autotab

!===============================================================================

subroutine sxy(x,amkd,ankd,amk,ank,akbar)

implicit none

integer, parameter :: ibins=36
integer :: i,ik,k,l
real, dimension(ibins+1) :: x
real, dimension(ibins,ibins) :: akbar
real, dimension(ibins) :: xave,ankd,ank,amkd,amk,am2,am3,am4,psi,f
real :: ap,pi,dm,dn,sm1,sm2,sm3,sm4,sm5,sn1,sn2,sn3,sn4,dm4,dm2

data ap/1.062500000/
data pi/3.141592654/

do l=1,ibins
   if (ankd(l) > 0) then
      xave(l)=amkd(l)/ankd(l)
   else
      xave(l)=0.
   endif
enddo

do k=1,ibins

! calculation of the 2nd, 3rd, and 4th moments of the mass distribution
! based on equation 8 in reference.

   am2(k)=ap*xave(k)*amkd(k)
   am3(k)=ap*ap*xave(k)*am2(k)
   am4(k)=ap*ap*ap*xave(k)*am3(k)

! these functions come out of the linear approximations used to integrate
! over partial bins.  they are defined:
!      psi(k) = nk(k+1)
!        f(k) = nk(k)
! where nk is the distribution function.  see equation 13 in reference.

   psi(k)=2./x(k)*(amkd(k)/x(k)-ankd(k))
   f(k)=2./x(k)*(2.*ankd(k)-amkd(k)/x(k))

! zeroing the tendencies on the moments.

   sm1=0.
   sm2=0.
   sm3=0.
   sm4=0.
   sm5=0.
   sn1=0.
   sn2=0.
   sn3=0.
   sn4=0.

! calculation of tendencies on moments

   do i=k,ibins
      dm=akbar(i,k)*(am2(k)*ankd(i)+amkd(k)*amkd(i))
      dn=akbar(i,k)*(ankd(k)*amkd(i)+amkd(k)*ankd(i))

      sm5=sm5+dm
      sn4=sn4+dn
   enddo

   if (k > 1) then
      sm3=akbar(k-1,k-1)*(am2(k-1)*ankd(k-1)+amkd(k-1)**2)
      sn2=akbar(k-1,k-1)*ankd(k-1)*amkd(k-1)
      dn=sn2
      dm=sm3
   endif

   do i=1,k-1
      dm4=akbar(k,i)*(ankd(k)*am2(i)+amkd(k)*amkd(i))
      sm4=sm4+dm4


      if (xave(k) >= x(k)) then
         dm2=akbar(k,i)*(4.*x(k)**2*psi(k)*amkd(i)  &
            +0.5*x(k)*(4.*psi(k)+f(k))*am2(i)  &
            -(psi(k)-f(k))*am3(i)  &
            -0.5/(x(k))*(psi(k)-f(k))*am4(i))
         sm2=sm2+dm2
         dn=akbar(k,i)*(2.*x(k)*psi(k)*amkd(i)  &
            +0.5*f(k)*am2(i)  &
            -0.5/(x(k))*(psi(k)-f(k))*am3(i))
         sn3=sn3+dn

      endif
   enddo

   do i=1,k-2
      ik=k-1
      if (xave(ik) > x(ik)) then

         dm=akbar(ik,i)*(4.*x(ik)**2*psi(ik)*amkd(i)  &
            +x(ik)/2.*(4.*psi(ik)+f(ik))*am2(i)  &
            -(psi(ik)-f(ik))*am3(i)  &
            -0.5/(x(ik))*(psi(ik)-f(ik))*am4(i))
         sm1=sm1+dm
         dn=akbar(ik,i)*(2.*x(ik)*psi(ik)*amkd(i)  &
            +0.5*f(ik)*am2(i)  &
            -0.5/(x(ik))*(psi(ik)-f(ik))*am3(i))
         sn1=sn1+dn

      endif
   enddo

   amk(k)=amkd(k)+sm1-sm2+sm3+sm4-sm5
   ank(k)=ankd(k)+sn1+sn2-sn3-sn4

enddo

return
end subroutine sxy

!===============================================================================

subroutine data(x,diam,akbar,ibins)

implicit none

integer :: l,i,j,kount,ibins,n
real :: pi,vpi,p,ap
real, dimension(ibins+1) :: x,diam
real, dimension(ibins,ibins) :: akbar
real, dimension(36,36) :: aabar

data (aabar( 1,n),n=1, 1) /-.47757E-01  /
data (aabar( 2,n),n=1, 2) /-.26460E+00,-.47965E-01 /
data (aabar( 3,n),n=1, 3) /-.82258E+00,-.26760E+00,-.20453E-01 /
data (aabar( 4,n),n=1, 4) /-.19050E+01,-.82072E+00,-.11992E+00, .78909E-01 /
data (aabar( 5,n),n=1, 5) /-.39171E+01,-.18915E+01,-.33270E+00, .41936E+00  &
   ,.34801E+00 /
data (aabar( 6,n),n=1, 6) /-.76415E+01,-.38808E+01,-.73737E+00, .14121E+01  &
   ,.18851E+01,.99793E+00 /
data (aabar( 7,n),n=1, 7) /-.14595E+02,-.75638E+01,-.14861E+01, .33598E+01  &
   ,.61219E+01, .54314E+01, .24751E+01 /
data (aabar( 8,n),n=1, 8) /-.27720E+02,-.14442E+02,-.28741E+01, .69895E+01  &
   ,.14394E+02, .17479E+02, .13500E+02, .57110E+01 /
data (aabar( 9,n),n=1, 9) /-.52737E+02,-.27428E+02,-.54729E+01, .13703E+02  &
   ,.29792E+02, .40971E+02, .43267E+02, .31185E+02, .12630E+02 /
data (aabar(10,n),n=1,10) /-.10083E+03,-.52188E+02,-.10391E+02, .26218E+02  &
   ,.58283E+02, .84686E+02, .10128E+03, .99726E+02, .69014E+02, .27176E+02 /
data (aabar(11,n),n=1,11) /-.19396E+03,-.99799E+02,-.19790E+02, .49801E+02  &
   ,.11143E+03, .16558E+03, .20922E+03, .23326E+03, .22039E+03, .14858E+03  &
   ,.57396E+02 /
data (aabar(12,n),n=1,12) /-.37536E+03,-.19200E+03,-.37896E+02, .94692E+02  &
   ,.21165E+03, .31650E+03, .40896E+03, .48169E+03, .51524E+03, .47402E+03  &
   ,.31389E+03, .11962E+03 /
data (aabar(13,n),n=1,13) /-.73047E+03,-.37164E+03,-.73015E+02, .18089E+03  &
   ,.40253E+03, .60115E+03, .78166E+03, .94143E+03, .10638E+04, .11078E+04  &
   ,.10008E+04, .65436E+03, .24691E+03 /
data (aabar(14,n),n=1,14) /-.14285E+04,-.72333E+03,-.14152E+03, .34764E+03  &
   ,.76925E+03, .11434E+04, .14846E+04, .17993E+04, .20789E+04, .22870E+04  &
   ,.23385E+04, .20854E+04, .13509E+04, .50600E+03 /
data (aabar(15,n),n=1,15) /-.41365E+04,-.20869E+04,-.40697E+03, .99310E+03  &
   ,.21878E+04, .32394E+04, .41995E+04, .51084E+04, .59888E+04, .68297E+04  &
   ,.75528E+04, .79583E+04, .76785E+04, .62489E+04, .76776E+03 /
data (aabar(16,n),n=1,16) / .63760E+04, .64739E+04, .65970E+04, .67516E+04  &
   ,.69451E+04, .71861E+04, .74835E+04, .78448E+04, .82709E+04, .87453E+04  &
   ,.92111E+04, .95276E+04, .94079E+04, .83797E+04, .26045E+04, .89777E+03 /
data (aabar(17,n),n=1,17) / .62974E+04, .63746E+04, .64717E+04, .65934E+04  &
   ,.67457E+04, .69355E+04, .71702E+04, .74571E+04, .78005E+04, .81957E+04  &
   ,.86163E+04, .89879E+04, .91399E+04, .87394E+04, .46530E+04, .26045E+04  &
   ,.89777E+03 /
data (aabar(18,n),n=1,18) / .62353E+04, .62963E+04, .63729E+04, .64689E+04  &
   ,.65889E+04, .67383E+04, .69233E+04, .71502E+04, .74238E+04, .77446E+04  &
   ,.81009E+04, .84538E+04, .87067E+04, .86514E+04, .59471E+04, .46530E+04  &
   ,.26045E+04, .89777E+03 /
data (aabar(19,n),n=1,19) / .61862E+04, .62344E+04, .62949E+04, .63707E+04  &
   ,.64653E+04, .65831E+04, .67290E+04, .69080E+04, .71250E+04, .73819E+04  &
   ,.76742E+04, .79815E+04, .82491E+04, .83524E+04, .66125E+04, .59471E+04  &
   ,.46530E+04, .26045E+04, .89777E+03 /
data (aabar(20,n),n=1,20) / .61474E+04, .61855E+04, .62334E+04, .62932E+04  &
   ,.63679E+04, .64608E+04, .65759E+04, .67172E+04, .68887E+04, .70932E+04  &
   ,.73291E+04, .75856E+04, .78311E+04, .79911E+04, .68735E+04, .66125E+04  &
   ,.59471E+04, .46530E+04, .26045E+04, .89777E+03 /
data (aabar(21,n),n=1,21) / .61166E+04, .61468E+04, .61847E+04, .62320E+04  &
   ,.62910E+04, .63644E+04, .64552E+04, .65668E+04, .67023E+04, .68644E+04  &
   ,.70531E+04, .72625E+04, .74738E+04, .76415E+04, .69140E+04, .68735E+04  &
   ,.66125E+04, .59471E+04, .46530E+04, .26045E+04, .89777E+03 /
data (aabar(22,n),n=1,22) / .60923E+04, .61162E+04, .61462E+04, .61836E+04  &
   ,.62303E+04, .62883E+04, .63600E+04, .64481E+04, .65553E+04, .66836E+04  &
   ,.68338E+04, .70027E+04, .71786E+04, .73330E+04, .68498E+04, .69140E+04  &
   ,.68735E+04, .66125E+04, .59471E+04, .46530E+04, .26045E+04, .89777E+03 /
data (aabar(23,n),n=1,23) / .60730E+04, .60919E+04, .61157E+04, .61453E+04  &
   ,.61823E+04, .62281E+04, .62848E+04, .63545E+04, .64392E+04, .65408E+04  &
   ,.66601E+04, .67953E+04, .69391E+04, .70729E+04, .67447E+04, .68498E+04  &
   ,.69140E+04, .68735E+04, .66125E+04, .59471E+04, .46530E+04, .26045E+04  &
   ,.89777E+03 /
data (aabar(24,n),n=1,24) / .60577E+04, .60727E+04, .60915E+04, .61150E+04  &
   ,.61443E+04, .61806E+04, .62254E+04, .62805E+04, .63475E+04, .64279E+04  &
   ,.65225E+04, .66304E+04, .67467E+04, .68590E+04, .66311E+04, .67447E+04  &
   ,.68498E+04, .69140E+04, .68735E+04, .66125E+04, .59471E+04, .46530E+04  &
   ,.26045E+04, .89777E+03 /
data (aabar(25,n),n=1,25) / .77967E+04, .78122E+04, .78316E+04, .78560E+04  &
   ,.78863E+04, .79242E+04, .79713E+04, .80294E+04, .81008E+04, .81878E+04  &
   ,.82924E+04, .84158E+04, .85571E+04, .87104E+04, .86265E+04, .88325E+04  &
   ,.90719E+04, .93363E+04, .95996E+04, .98007E+04, .98157E+04, .94274E+04  &
   ,.83361E+04, .63023E+04, .57988E+03 /
data (aabar(26,n),n=1,26) / .69349E+04, .69458E+04, .69595E+04, .69766E+04  &
   ,.69979E+04, .70244E+04, .70573E+04, .70978E+04, .71473E+04, .72072E+04  &
   ,.72788E+04, .73623E+04, .74565E+04, .75566E+04, .74715E+04, .76064E+04  &
   ,.77647E+04, .79435E+04, .81311E+04, .82983E+04, .83827E+04, .82640E+04  &
   ,.77406E+04, .65488E+04, .15807E+04, .51662E+03 /
data (aabar(27,n),n=1,27) / .61704E+04, .61781E+04, .61877E+04, .61997E+04  &
   ,.62147E+04, .62333E+04, .62562E+04, .62843E+04, .63186E+04, .63598E+04  &
   ,.64086E+04, .64648E+04, .65271E+04, .65912E+04, .65100E+04, .65961E+04  &
   ,.66976E+04, .68135E+04, .69382E+04, .70574E+04, .71390E+04, .71194E+04  &
   ,.68816E+04, .62379E+04, .26526E+04, .14083E+04, .46025E+03 /
data (aabar(28,n),n=1,28) / .54916E+04, .54971E+04, .55038E+04, .55123E+04  &
   ,.55228E+04, .55357E+04, .55517E+04, .55712E+04, .55949E+04, .56232E+04  &
   ,.56562E+04, .56938E+04, .57345E+04, .57747E+04, .57001E+04, .57533E+04  &
   ,.58161E+04, .58880E+04, .59661E+04, .60426E+04, .61008E+04, .61062E+04  &
   ,.59940E+04, .56500E+04, .31742E+04, .23632E+04, .12546E+04, .41004E+03 /
data (aabar(29,n),n=1,29) / .48886E+04, .48924E+04, .48971E+04, .49031E+04  &
   ,.49104E+04, .49195E+04, .49306E+04, .49441E+04, .49604E+04, .49797E+04  &
   ,.50020E+04, .50269E+04, .50530E+04, .50774E+04, .50108E+04, .50422E+04  &
   ,.50792E+04, .51212E+04, .51667E+04, .52111E+04, .52447E+04, .52486E+04  &
   ,.51861E+04, .49913E+04, .32935E+04, .28279E+04, .21054E+04, .11177E+04  &
   ,.36530E+03 /
data (aabar(30,n),n=1,30) / .43524E+04, .43551E+04, .43585E+04, .43626E+04  &
   ,.43678E+04, .43741E+04, .43818E+04, .43912E+04, .44024E+04, .44155E+04  &
   ,.44304E+04, .44467E+04, .44631E+04, .44771E+04, .44188E+04, .44361E+04  &
   ,.44561E+04, .44786E+04, .45022E+04, .45241E+04, .45384E+04, .45339E+04  &
   ,.44893E+04, .43663E+04, .31847E+04, .29342E+04, .25193E+04, .18757E+04  &
   ,.99579E+03, .32545E+03 /
data (aabar(31,n),n=1,31) / .38756E+04, .38775E+04, .38799E+04, .38828E+04  &
   ,.38864E+04, .38908E+04, .38961E+04, .39026E+04, .39102E+04, .39191E+04  &
   ,.39290E+04, .39395E+04, .39494E+04, .39568E+04, .39066E+04, .39149E+04  &
   ,.39241E+04, .39340E+04, .39435E+04, .39507E+04, .39516E+04, .39392E+04  &
   ,.39006E+04, .38129E+04, .29707E+04, .28372E+04, .26141E+04, .22445E+04  &
   ,.16710E+04, .88715E+03, .28994E+03 /
data (aabar(32,n),n=1,32) / .30106E+04, .30118E+04, .30132E+04, .30149E+04  &
   ,.30171E+04, .30197E+04, .30229E+04, .30266E+04, .30309E+04, .30357E+04  &
   ,.30408E+04, .30456E+04, .30491E+04, .30494E+04, .30032E+04, .30013E+04  &
   ,.29981E+04, .29929E+04, .29844E+04, .29706E+04, .29480E+04, .29111E+04  &
   ,.28504E+04, .27503E+04, .20956E+04, .19717E+04, .17926E+04, .15245E+04  &
   ,.11243E+04, .55892E+03, .22135E+03, .00000E+00 /
data (aabar(33,n),n=1,33) / .23888E+04, .23895E+04, .23903E+04, .23914E+04  &
   ,.23927E+04, .23943E+04, .23962E+04, .23983E+04, .24007E+04, .24033E+04  &
   ,.24057E+04, .24075E+04, .24077E+04, .24049E+04, .23645E+04, .23582E+04  &
   ,.23497E+04, .23382E+04, .23225E+04, .23007E+04, .22699E+04, .22258E+04  &
   ,.21613E+04, .20655E+04, .15572E+04, .14495E+04, .13057E+04, .11057E+04  &
   ,.82055E+03, .41732E+03, .17636E+03, .00000E+00, .00000E+00 /
data (aabar(34,n),n=1,34) / .18955E+04, .18959E+04, .18964E+04, .18971E+04  &
   ,.18979E+04, .18988E+04, .18999E+04, .19011E+04, .19024E+04, .19036E+04  &
   ,.19045E+04, .19047E+04, .19033E+04, .18990E+04, .18647E+04, .18567E+04  &
   ,.18462E+04, .18326E+04, .18145E+04, .17905E+04, .17581E+04, .17138E+04  &
   ,.16525E+04, .15662E+04, .11695E+04, .10771E+04, .95995E+03, .80547E+03  &
   ,.59516E+03, .30433E+03, .13287E+03, .00000E+00, .00000E+00, .00000E+00 /
data (aabar(35,n),n=1,35) / .15041E+04, .15044E+04, .15047E+04, .15051E+04  &
   ,.15056E+04, .15061E+04, .15067E+04, .15074E+04, .15080E+04, .15084E+04  &
   ,.15085E+04, .15079E+04, .15058E+04, .15011E+04, .14725E+04, .14642E+04  &
   ,.14536E+04, .14399E+04, .14221E+04, .13988E+04, .13682E+04, .13274E+04  &
   ,.12725E+04, .11976E+04, .88682E+03, .80897E+03, .71340E+03, .59221E+03  &
   ,.43360E+03, .22074E+03, .97241E+02, .00000E+00, .00000E+00, .00000E+00  &
   ,.00000E+00 /
data (aabar(36,n),n=1,36) / .11936E+04, .11938E+04, .11940E+04, .11942E+04  &
   ,.11945E+04, .11948E+04, .11951E+04, .11954E+04, .11957E+04, .11957E+04  &
   ,.11954E+04, .11943E+04, .11921E+04, .11876E+04, .11640E+04, .11563E+04  &
   ,.11464E+04, .11337E+04, .11174E+04, .10963E+04, .10689E+04, .10330E+04  &
   ,.98554E+03, .92214E+03, .67808E+03, .61344E+03, .53582E+03, .44015E+03  &
   ,.31885E+03, .16089E+03, .70536E+02, .00000E+00, .00000E+00, .00000E+00  &
   ,.00000E+00, .00000E+00 /

data pi,vpi/3.141592654,0.52359877/

! calculating the mass categories (c.g.s) with the lowest diameter
! of 3.125 microns and mass doubling every bin
!
diam(1)=1.5625*2.e-04
x(1)=pi/6.*diam(1)**3.
do l=2,ibins+1
   x(l)=2.*x(l-1)
   diam(l)=(6./pi*x(l))**(0.333333333333)
enddo

l=1
p=2.0**(1/l)
ap=0.5+(p+1.0)*(p+1.0)/(8.0*p)

!
! long's collection kernel as calculated by myself (1-36) with (1-36)
! weighted for (x+y)
!
kount=0
do i=1,ibins
   do j=1,i
      kount=kount+1
      akbar(i,j)=aabar(i,j)
      if (akbar(i,j) < 0.) akbar(i,j)=0.
      akbar(36,i)=0.
   enddo
enddo

do j=1,ibins
   do i=1,j
      akbar(i,j)=akbar(j,i)
   enddo
enddo

return
end subroutine data

!===============================================================================

subroutine initg(r1,r2,n1,n2,gnu1,gnu2,diam,x,amk,ank  &
       ,ank1,amk1,ank2,amk2,ibins,ithresh)
       
implicit none

integer :: i,ibins,ithresh
real :: r1,r2,n1,n2,gnu1,gnu2,dn1,dn2,trunc,trunc1,fac1,fac2,pi  &
       ,gamp,dmean,sum,sumn,xntot,xr3,ex1
real, dimension(ibins+1) :: x,diam
real, dimension(ibins) :: amk,ank,amk1,ank1,amk2,ank2

real :: gammln,gammp

! * *
! * Initial double gamma distribution: n(D) = n1(D) + n2(D)
! * *

data pi/3.141592654/
data ex1/0.333333333/

! * *
! * gamma spectrum
! * *

gamp = exp(gammln(gnu1))
dmean = (6.*r1/(pi*n1))**ex1   !mass mean diam
dn1 = dmean * (exp(gammln(gnu1) - gammln(gnu1+3.))) ** ex1

do i=1,ibins
  ank1(i)=0.
  amk1(i)=0.
  ank2(i)=0.
  amk2(i)=0.
enddo

sum=0.
sumn=0.
do i=1,ithresh
 fac1=gammp(gnu1,diam(i)/dn1)
 fac2=gammp(gnu1,diam(i+1)/dn1)
 trunc1=fac2-fac1
 ank1(i)=n1*trunc1

 fac1=gammp(gnu1+3.,(diam(i)/dn1))
 fac2=gammp(gnu1+3.,(diam(i+1)/dn1))
 trunc=fac2-fac1
 amk1(i)=r1*trunc

 sum=sum+amk1(i)
 sumn=sumn+ank1(i)
enddo

! * Scale to exactly r1,n1

!      do i=1,ithresh
!        ank1(i)=ank1(i)*n1/sumn
!        amk1(i)=amk1(i)*r1/sum
!      enddo

gamp = exp(gammln(gnu2))
dmean = (6. * r2 / (pi * n2)) ** ex1   !mass mean diam
dn2 = dmean * (exp(gammln(gnu2) - gammln(gnu2+3.))) ** ex1

sum=0.
sumn=0.

do i=ithresh+1,ibins
 fac1=gammp(gnu2,diam(i)/dn2)
 fac2=gammp(gnu2,diam(i+1)/dn2)
 trunc=fac2-fac1
 ank2(i)=n2*trunc

 fac1=gammp(gnu2+3.,(diam(i)/dn2))
 fac2=gammp(gnu2+3.,(diam(i+1)/dn2))
 trunc=fac2-fac1
 amk2(i)=r2*trunc

 sum=sum+amk2(i)
 sumn=sumn+ank2(i)

enddo

! * Scale to exactly r2,n2

!      do i=ithresh+1,ibins
!        ank2(i)=ank2(i)*n2/sumn
!        amk2(i)=amk2(i)*r2/sum
!      enddo

do i=1,ibins
   ank(i)=ank1(i)+ank2(i)
   amk(i)=amk1(i)+amk2(i)
enddo

    xntot  = 0.0
    xr3=0.

    do i=1,ibins
       xr3       = xr3 + amk(i)
       xntot     = xntot + ank(i)
    enddo
!          write(io6,*) 'LWC,N,gnus,ibins ',xr3,xntot,gnu1,gnu2,ibins

return
end subroutine initg

!===============================================================================

subroutine sumn(ank,amk,imin,imax,ibins,sun,sum)

implicit none

integer :: imin,imax,ibins
real, dimension(ibins) :: ank,amk
real :: sun,sum

integer :: i

  sum=0.
  sun=0.

  do i=imin,imax
   sun=sun+ank(i)
   sum=sum+amk(i)
  enddo

return
end subroutine sumn

!===============================================================================

subroutine tabhab()

use micro_coms, only: jhabtab
use misc_coms,  only: io6

implicit none

integer, parameter :: nhab=0
integer :: it,is

if (nhab ==  0) write(io6,*) 'VARIABLE HABIT PREDICTION'
if (nhab ==  3) write(io6,*) 'ASSUMED HABIT IS COLUMNS'
if (nhab ==  8) write(io6,*) 'ASSUMED HABIT IS HEX PLATES'
if (nhab ==  9) write(io6,*) 'ASSUMED HABIT IS DENDRITES'
if (nhab == 10) write(io6,*) 'ASSUMED HABIT IS NEEDLES'
if (nhab == 11) write(io6,*) 'ASSUMED HABIT IS ROSETTES'
!c    if (nhab .eq.  x) write(io6,*) 'ASSUMED HABIT IS SPHERES'

! nt is temp, ns = satur (liq)

do it = 1,31
   do is = 1,100
      if (nhab == 0) then
         if (it >= 0 .and. it <= 2) then
            if (is <= 95) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 8
               jhabtab(it,is,2) = 12
            endif
         else if(it > 2 .and. it <= 4) then
            if (is < 90) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 8
               jhabtab(it,is,2) = 12
            endif
         else if(it > 4 .and. it <= 6) then
            if (is < 85) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 10
               jhabtab(it,is,2) = 14
            endif
         else if(it > 6 .and. it <= 9) then
            if (is < 90) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 10
               jhabtab(it,is,2) = 14
            endif
         else if(it > 9 .and. it <= 22) then
            if (is < 90) then
               jhabtab(it,is,1) = 8
               jhabtab(it,is,2) = 12
            else
               jhabtab(it,is,1) = 9
               jhabtab(it,is,2) = 13
            endif
         elseif(it > 22 .and. it <= 30) then
            if (is < 80) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 10
               jhabtab(it,is,2) = 14
            endif
         elseif(it > 30) then
            if (is < 90) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 11
               jhabtab(it,is,2) = 15
            endif
         endif
      else
         jhabtab(it,is,1) = nhab
         jhabtab(it,is,2) = nhab + 4
         if (nhab == 3) jhabtab(it,is,2) = 4
      endif
   enddo
enddo
return
end subroutine tabhab
