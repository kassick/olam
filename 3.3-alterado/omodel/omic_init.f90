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
subroutine micinit()

use mem_basic,   only: theta, press, rho

use mem_micro,   only: sh_r, sh_p, sh_s, sh_a, sh_g, sh_h,  &
                       accpr, accpp, accps, accpa, accpg, accph,  &
                       pcprr, pcprp, pcprs, pcpra, pcprg, pcprh,  &
                       con_c, con_r, con_p, con_s, con_a, con_g, con_h,  &
                       con_ccn, con_ifn, q2, q6, q7, pcpgr, qpcpgr, dpcpgr

use misc_coms,   only: io6, dtlm
use consts_coms, only: p00, rocp
use mem_ijtabs,  only: jtab_w, mrls
use mem_grid,    only: mza, zm, dzt, dzit

use micro_coms,  only: level, irain, ipris, isnow, iaggr, igraup, ihail, jnmb,  &
                       icloud, ncat, gnu, emb0, emb1, nhcat,  &
                       cfmas, pwmas, cfvt, pwvt, npairc, coltabc, nembc,  &
                       npairr, coltabr, alloc_sedimtab

implicit none

real :: tair(mza)  ! automatic array

integer :: i,j,iw,k,ilhcat,ilcat,lcat,lhcat,nip,idum,nd1,nd2,mrl
logical :: l1,l2
character*80 dataline,cname

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAM_MICINIT_IN,rst_buf)
#endif

call micinit_gam()

if (level < 3) return

call make_autotab()
call haznuc()
call tabmelt()
call tabhab()

call alloc_sedimtab(mza)
call mksedim_tab(mza,zm,dzt,dzit)

do mrl = 1,mrls
   call homfrzcl(dtlm(mrl),mrl)
enddo

! Initialize 3D and 2D microphysics fields

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   tair(1:mza) = theta(1:mza,iw) * (press(1:mza,iw) / p00) ** rocp

   if (irain >= 1)  then
      sh_r(1:mza,iw) = 0.
      accpr(iw) = 0.
      pcprr(iw) = 0.
      q2(1:mza,iw) = tair(1:mza) - 193.15
   endif
   
   if (ipris >= 1)  then
      sh_p(1:mza,iw) = 0.
      accpp(iw) = 0.
      pcprp(iw) = 0.
   endif

   if (isnow >= 1)  then
      sh_s(1:mza,iw) = 0.
      accps(iw) = 0.
      pcprs(iw) = 0.
   endif

   if (iaggr >= 1)  then
      sh_a(1:mza,iw) = 0.
      accpa(iw) = 0.
      pcpra(iw) = 0.
   endif

   if (igraup >= 1) then
      sh_g(1:mza,iw) = 0.
      accpg(iw) = 0.
      pcprg(iw) = 0.
      q6(1:mza,iw) = .5 * min(0.,tair(1:mza) - 273.15)
   endif

   if (ihail >= 1)  then
      sh_h(1:mza,iw) = 0.
      accph(iw) = 0.
      pcprh(iw) = 0.
      q7(1:mza,iw) = .5 * min(0.,tair(1:mza) - 273.15)
   endif

   if (jnmb(1) == 5) con_c(1:mza,iw) = 0.
   if (jnmb(2) == 5) con_r(1:mza,iw) = 0.
   if (jnmb(3) == 5) con_p(1:mza,iw) = 0.
   if (jnmb(4) == 5) con_s(1:mza,iw) = 0.
   if (jnmb(5) == 5) con_a(1:mza,iw) = 0.
   if (jnmb(6) == 5) con_g(1:mza,iw) = 0.
   if (jnmb(7) == 5) con_h(1:mza,iw) = 0.
   if (icloud  == 7) con_ccn(1:mza,iw) = 0.  !?? NEEDS FIELD SPECIFICATION HERE
   if (ipris   == 7) con_ifn(1:mza,iw) =  &
           1.e5 * rho(1:mza,iw) ** 5.4  ! This is default, could change init values

   pcpgr(iw)  = 0.
   qpcpgr(iw) = 0.
   dpcpgr(iw) = 0.

enddo
call rsub('Wc',7)

! Make collection table

call mkcoltb()

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_MICINIT_OUT,rst_buf)
#endif

return

end subroutine micinit

!===============================================================================

subroutine micinit_gam()

use micro_coms,  only: icloud, irain, ipris, isnow, iaggr, igraup, ihail,  &
                       parm, nhcat, shapefac, cfmas, pwmas, cfvt, pwvt,  &
                       ncat, emb0, emb1, gnu, rxmin, level, sl, sc, sj,  &
                       cparm, rparm, sparm, aparm, gparm, hparm, sk,  &
                       dps, dps2, rictmin, rictmax, nembc, lcat_lhcat,  &
                       emb0log, emb1log, emb2, cfmasi, pwmasi, pwen0,  &
                       pwemb0, ch3, cdp1, pwvtmasi, jnmb, cfemb0, cfen0,  &
                       dnfac, vtfac, frefac1, frefac2, sipfac, cfmasft,  &
                       dict, dpsmi, gamm, gamn1, ngam, gam, gaminc,  &
                       gamsip13, gamsip24
                      
use consts_coms, only: alvl, alvi, alli
use misc_coms,   only: io6

implicit none

integer :: lhcat,khcat,lcat,igam,nd1,nd2,nip,ilcat,ilhcat,idum
real :: c1,glg,glg1,glg2,glgm,glgc,glgmv,gym,flngi,dpsi,embsip,dnsip
real :: gammln,gammp,gammq

real, dimension(9,15) :: dstprms

data dstprms/ &
!-------------------------------------------------------------------------
!shapefac  cfmas  pwmas     cfvt   pwvt    dmb0     dmb1  gnu  rxmin
!------------------------------------------------------------------------- 
    .5,     524.,    3.,   3173.,    2.,  2.e-6,  40.e-6,  9., 1.e-12, & !cloud
    .5,     524.,    3.,    149.,    .5,  .1e-3,   5.e-3,  2.,  1.e-9, & !rain
  .179,    110.8,  2.91, 5.769e5,  1.88, 15.e-6, 125.e-6,  2., 1.e-12, & !pris col
  .179, 2.739e-3,  1.74, 188.146,  .933,  .1e-3,  10.e-3,  2.,  1.e-9, & !snow col
    .5,     .496,   2.4,   3.084,    .2,  .1e-3,  10.e-3,  2.,  1.e-9, & !aggreg
    .5,     157.,    3.,    93.3,    .5,  .1e-3,   5.e-3,  2.,  1.e-9, & !graup
    .5,     471.,    3.,    161.,    .5,  .8e-3,  10.e-3,  2.,  1.e-9, & !hail 
 .0429,    .8854,   2.5,    316.,  1.01,    .00,     .00, 00.,    00., & !pris hex
 .3183,  .377e-2,    2.,    316.,  1.01,    .00,     .00, 00.,    00., & !pris den
 .1803,  1.23e-3,   1.8, 5.769e5,  1.88,    .00,     .00, 00.,    00., & !pris ndl
    .5,    .1001, 2.256,  3.19e4,  1.66,    .00,     .00, 00.,    00., & !pris ros
 .0429,    .8854,   2.5,   4.836,   .25,    .00,     .00, 00.,    00., & !snow hex
 .3183,  .377e-2,    2.,   4.836,   .25,    .00,     .00, 00.,    00., & !snow den
 .1803,  1.23e-3,   1.8, 188.146,  .933,    .00,     .00, 00.,    00., & !snow ndl
    .5,    .1001, 2.256, 1348.38, 1.241,    .00,     .00, 00.,    00.  / !snow ros

! Initialize arrays based on microphysics namelist parameters

parm(1) = cparm
parm(2) = rparm
parm(3) = 100.e-6  ! Obsolete
parm(4) = sparm
parm(5) = aparm
parm(6) = gparm
parm(7) = hparm

if (icloud <= 1) parm(1) = .3e9
if (irain  == 1) parm(2) = .1e-2
if (ipris  == 1) parm(3) = 100.e-6  ! Obsolete
if (isnow  == 1) parm(4) = .1e-2
if (iaggr  == 1) parm(5) = .1e-2
if (igraup == 1) parm(6) = .1e-2
if (ihail  == 1) parm(7) = .3e-2

!  Copy individual arrays from above data table

do lhcat = 1,nhcat
   shapefac(lhcat) = dstprms(1,lhcat)
   cfmas   (lhcat) = dstprms(2,lhcat)
   pwmas   (lhcat) = dstprms(3,lhcat)
   cfvt    (lhcat) = dstprms(4,lhcat)
   pwvt    (lhcat) = dstprms(5,lhcat)
enddo

do lcat = 1,ncat
   emb0  (lcat) = cfmas(lcat) * dstprms(6,lcat) ** pwmas(lcat)
   emb1  (lcat) = cfmas(lcat) * dstprms(7,lcat) ** pwmas(lcat)
   gnu   (lcat) = dstprms(8,lcat)
   rxmin (lcat) = dstprms(9,lcat)
enddo

if (level < 3) RETURN

! Initialize constants for vapor diffusion

sl(1) = alvl
sl(2) = alvi
sc(1) = 4186.
sc(2) = 2093.  ! 2106 is correct value
sj(1) = 0
sj(2) = 1
sj(3) = 0
sj(4) = 0
sj(5) = 0
sj(6) = 1
sj(7) = 1
sk(1) = alli
sk(2) = 0.

dps = 125.e-6
dps2 = dps ** 2
rictmin = 1.0001
rictmax = 0.9999 * float(nembc)

do lhcat = 1,nhcat
   lcat = lcat_lhcat(lhcat)

   emb0log(lcat) = log(emb0(lcat))
   emb1log(lcat) = log(emb1(lcat))

   emb2  (lhcat)   = cfmas(lhcat) * parm(lcat) ** pwmas(lhcat)
   cfmasi(lhcat)   = 1. / cfmas(lhcat)
   pwmasi(lhcat)   = 1. / pwmas(lhcat)
   pwen0(lhcat)    = 1. / (pwmas(lhcat) + 1.)
   pwemb0(lhcat)   = pwmas(lhcat) / (pwmas(lhcat) + 1.)
   ch3(lhcat)      = pwvt(lhcat) * pwmasi(lhcat)
   cdp1(lhcat)     = pwmasi(lhcat) * (1.5 + .5 * pwvt(lhcat))
   pwvtmasi(lhcat) = pwvt(lhcat) * pwmasi(lhcat)

   c1 = 1.5 + .5 * pwvt(lhcat)
   glg = gammln(gnu(lcat))
   glg1 = gammln(gnu(lcat) + 1.)
   glg2 = gammln(gnu(lcat) + 2.)
   glgm = gammln(gnu(lcat) + pwmas(lhcat))
   glgc = gammln(gnu(lcat) + c1)
   glgmv = gammln(gnu(lcat) + pwmas(lhcat) + pwvt(lhcat))

   if (jnmb(lcat) == 3) then
      cfemb0(lhcat) = cfmas(lhcat) * exp(glgm - glg)  &
         ** pwen0(lhcat) * (1. / parm(lcat)) ** pwemb0(lhcat)
      cfen0(lhcat) = parm(lcat) * (exp(glg - glgm) / parm(lcat))  &
         ** pwen0(lhcat)
   endif

   dnfac(lhcat) = (cfmasi(lhcat) * exp(glg - glgm)) ** pwmasi(lhcat)

   vtfac(lhcat) = cfvt(lhcat) * exp(glgmv - glgm)  &
      * (cfmasi(lhcat) * exp(glg - glgm)) ** (pwvt(lhcat) *pwmasi(lhcat))

   frefac1(lhcat) = shapefac(lhcat) * exp(glg1 - glg)  &
      * (cfmasi(lhcat) * exp(glg - glgm)) ** pwmasi(lhcat)

   frefac2(lhcat) = shapefac(lhcat) * 0.229 * sqrt(cfvt(lcat))  &
      * (cfmasi(lhcat) * exp(glg - glgm)) ** (pwmasi(lhcat) * c1)  &
      * exp(glgc - glg)

   sipfac(lhcat) = .785 * exp(glg2 - glg)  &
      * (cfmasi(lhcat) * exp(glg - glgm)) ** (2. * pwmasi(lhcat))

   cfmasft(lhcat) = cfmas(lhcat) * exp(gammln  &
      (gnu(lcat) + pwmas(lhcat)) - gammln(gnu(lcat)))

   dict(lcat) = float(nembc-1) / (emb1log(lcat) - emb0log(lcat))

   dpsmi(lhcat) = 1. / (cfmas(lhcat) * dps ** pwmas(lhcat))
   if (lhcat <= 4) gamm(lhcat) = exp(glg)
   if (lhcat <= 4) gamn1(lhcat) = exp(glg1)

! gam1   :  the integral of the pristine distribution from dps to infty
! gam2   :  the integral of the snow dist. from 0 to dps
! gam3   :  values of the exponential exp(-dps/dn)

enddo

flngi = 1. / float(ngam)
do igam = 1,ngam
   dpsi = dps * 1.e6 / float(igam)

   gam(igam,1) = gammq(gnu(3) + 1., dpsi)
   gam(igam,2) = gammp(gnu(4) + 1., dpsi)
   gam(igam,3) = exp(-dpsi)

   GAMINC(igam,1) = GAMMQ(GNU(3),dpsi)
   GAMINC(igam,2) = GAMMP(GNU(4),dpsi)

   embsip = emb1(1) * float(igam) * flngi
   dnsip = dnfac(1) * embsip ** pwmasi(1)
   gamsip13(igam) = gammp(gnu(1),13.e-6/dnsip)
   gamsip24(igam) = gammq(gnu(1),24.e-6/dnsip)
enddo

return
end subroutine micinit_gam

!===============================================================================

subroutine jnmbinit()

use micro_coms, only: level, jnmb, icloud, irain, ipris, isnow, iaggr,  &
                      igraup, ihail
use misc_coms,  only: io6

implicit none

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAN_JNMBINIT_IN,rst_buf)
#endif

if (level /= 3) then

   if (level <= 1) then
      jnmb(1) = 0
   else
      jnmb(1) = 4
   endif

   jnmb(2) = 0
   jnmb(3) = 0
   jnmb(4) = 0
   jnmb(5) = 0
   jnmb(6) = 0
   jnmb(7) = 0

else

   jnmb(1) = icloud
   jnmb(2) = irain
   jnmb(3) = ipris
   jnmb(4) = isnow
   jnmb(5) = iaggr
   jnmb(6) = igraup
   jnmb(7) = ihail

   if (icloud == 1) jnmb(1) = 4
   if (irain  == 1) jnmb(2) = 2
   if (ipris  == 1) jnmb(3) = 5
   if (isnow  == 1) jnmb(4) = 2
   if (iaggr  == 1) jnmb(5) = 2
   if (igraup == 1) jnmb(6) = 2
   if (ihail  == 1) jnmb(7) = 2

   if (irain == 5 .or. isnow == 5 .or. iaggr == 5 .or.  &
      igraup == 5 .or. ihail == 5) then

      if (irain  >= 1) jnmb(2) = 5
      if (isnow  >= 1) jnmb(4) = 5
      if (iaggr  >= 1) jnmb(5) = 5
      if (igraup >= 1) jnmb(6) = 5
      if (ihail  >= 1) jnmb(7) = 5

   endif

endif
#ifdef OLAM_RASTRO
call rst_event_s_f(OLAN_JNMBINIT_OUT,rst_buf)
#endif
return
end subroutine jnmbinit

