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
subroutine auto_accret(j1,j2  &
   ,dtl0,rx,cx,qx,emb,denfac,rxfer,qrxfer,enxfer)

use micro_coms, only: mza0, ncat, rxmin, cfmasi, pwmasi, d1min, d1max,  &
                      r2min, r2max, d2min, d2max, nd1cc, d1ecr, r2ecr,  &
                      nd2cr, r2ecr, nd2rr, r1tabcc, c1tabcc, c2tabcc,  &
                      r1tabcr, c1tabcr, c2tabrr, r2err
use misc_coms,  only: io6

implicit none

integer, intent(in) :: j1
integer, intent(in) :: j2

real, intent(in) :: dtl0

real, intent(in) :: rx  (mza0,ncat)
real, intent(in) :: cx  (mza0,ncat)
real, intent(in) :: qx  (mza0,ncat)
real, intent(in) :: emb (mza0,ncat)

real, intent(in) :: denfac(mza0)

real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: k,id1cc,id1cr,id1crn,ir2cr,id2cr,ir2rr,id2rr
real :: dtl3,dtl6,dmb1cgs,dmb2cgs,r2cgs,en1cgs,ad1,ar2,d2minx,ad2  &
       ,bd1,br2,bd2,d2e,bd1cc,bd1cr,br2cr,bd2cr,br2rr,bd2rr,wd1cr,wr2dr  &
       ,wr2rr,wd2rr,tm1cc,tn1cc,tn2cc,tm1cr,tn1cr,tn2rr,en1cgs_2  &
       ,um1cc,un1cc,un2cc,um1cr,un1cr,un2rr,um2,un1  &
       ,wr2cr,df_dtl3,df_dtl6

dtl3 = 1.e3 * dtl0
dtl6 = 1.e6 * dtl0

do k = j1,j2

   if (rx(k,1) < rxmin(1)) cycle

! This subroutine works in cgs units, so convert inputs from mks

   dmb1cgs = 100. * (emb(k,1) * cfmasi(1)) ** pwmasi(1)  ! cm
   dmb2cgs = 100. * (emb(k,2) * cfmasi(2)) ** pwmasi(2)  ! cm
   r2cgs = 1.e-3 * rx(k,2)  ! removed rhoa factor        ! g/cm^3
   en1cgs = 1.e-6 * cx(k,1) ! removed rhoa factor        ! #/cm^3

   ad1 = max(d1min,min(d1max,dmb1cgs))
   ar2 = max(r2min,min(r2max,r2cgs))
    d2minx = max(d2min,(r2cgs / (.1 * .5236)) ** pwmasi(2))
   ad2 = max(d2minx,min(d2max,dmb2cgs))

   bd1 = alog10(ad1/d1min)
   br2 = alog10(ar2/r2min)
   bd2 = alog10(ad2/d2minx)
    d2e =  alog10(d2max / d2minx)

   bd1cc = float(nd1cc-1) * (ad1 - d1min) / (d1max - d1min) + 1.
   bd1cr = bd1 / d1ecr + 1.
   br2cr = br2 / r2ecr + 1.
   bd2cr = bd2 / d2e * float(nd2cr-1) + 1.
   br2rr = br2 / r2err + 1.
   bd2rr = bd2 / d2e * float(nd2rr-1) + 1.

!      id1cc  =  int(bd1cc)
   id1cc  =  nint(bd1cc)
   id1cr  =  int(bd1cr)
   id1crn = nint(bd1cr)
   ir2cr  =  int(br2cr)
   id2cr  = nint(bd2cr)
   ir2rr  =  int(br2rr)
   id2rr  =  int(bd2rr)

   wd1cr = bd1cr - float(id1cr)
   wr2cr = br2cr - float(ir2cr)
   wr2rr = br2rr - float(ir2rr)
   wd2rr = bd2rr - float(id2rr)

   tm1cc =                            r1tabcc(id1cc)

   tn1cc =                            c1tabcc(id1cc)

   tn2cc =                            c2tabcc(id1cc)

   tm1cr = (1.-wd1cr) * ((1.-wr2cr) * r1tabcr(id1cr  ,ir2cr  ,id2cr)   &
         +                   wr2cr  * r1tabcr(id1cr  ,ir2cr+1,id2cr))  &
         +     wd1cr  * ((1.-wr2cr) * r1tabcr(id1cr+1,ir2cr  ,id2cr)   &
         +                   wr2cr  * r1tabcr(id1cr+1,ir2cr+1,id2cr))

   tn1cr =               (1.-wr2cr) * c1tabcr(id1crn,ir2cr  ,id2cr)  &
         +                   wr2cr  * c1tabcr(id1crn,ir2cr+1,id2cr)

   tn2rr = (1.-wd2rr) * ((1.-wr2rr) * c2tabrr(ir2rr  ,id2rr  )   &
         +                   wr2rr  * c2tabrr(ir2rr+1,id2rr  ))  &
         +     wd2rr  * ((1.-wr2rr) * c2tabrr(ir2rr  ,id2rr+1)   &
         +                   wr2rr  * c2tabrr(ir2rr+1,id2rr+1))

   en1cgs_2 = en1cgs ** 2

   df_dtl3 = denfac(k) * dtl3  ! includes density factor on fall velocity
   df_dtl6 = denfac(k) * dtl6  ! includes density factor on fall velocity

   um1cc = tm1cc * en1cgs_2      * df_dtl3  ! kg/m^3 this timestep
   un1cc = tn1cc * en1cgs_2      * df_dtl6  ! #/m^3  this timestep
   un2cc = tn2cc * en1cgs_2      * df_dtl6  ! #/m^3  this timestep
   um1cr = 10. ** tm1cr * en1cgs * df_dtl3  ! kg/m^3 this timestep
   un1cr = 10. ** tn1cr * en1cgs * df_dtl6  ! #/m^3  this timestep
   un2rr = 10. ** tn2rr          * df_dtl6  ! #/m^3  this timestep

   um2 = min(rx(k,1),(um1cc + um1cr))  ! limit mass xfer to available cloud mass
   un1 = min(cx(k,1),(un1cc + un1cr))  ! limit cloud number loss to available number

   rxfer(k,1,2)  =  rxfer(k,1,2) + um2
   qrxfer(k,1,2) = qrxfer(k,1,2) + um2 * qx(k,1)
   enxfer(k,1,1) = enxfer(k,1,1) + un1 - un2cc
   enxfer(k,1,2) = enxfer(k,1,2) + un2cc

! no collis breakup yet - do not use next line but use col(2,2) in 3d micro instead

!cc         enxfer(k,2,2) = enxfer(k,2,2) + un2rr

! aerosol loss here?

enddo
return
end subroutine auto_accret

!===============================================================================

subroutine effxy(lpw0,k1,k2,rx,qr,emb,tx,eff)

use micro_coms, only: mza0, ncat, jnmb, rxmin, neff
use misc_coms,  only: io6

implicit none

integer, intent(in) :: lpw0

integer, intent(in) :: k1(10)
integer, intent(in) :: k2(10)

real, intent(in)  :: rx  (mza0,ncat)
real, intent(in)  :: qr  (mza0,ncat)
real, intent(in)  :: emb (mza0,ncat)
real, intent(in)  :: tx  (mza0,ncat)

real, intent(out) :: eff (mza0,neff)

integer :: k
real :: dmr
save

!     1 = rp,rs,ra,rg,rh

if (jnmb(2) >= 1 .and. jnmb(3) >= 1) then
   eff(lpw0:mza0,1) = 1.0
endif

!     2 = cs,ca

if (jnmb(2) >= 1 .or. jnmb(3) >= 1) then
   do k = k1(1),k2(1)

! Rough fit from Pruppacher and Klett Fig. 14-14 p. 496:
!  close to curve for 404 microns.  Replace with auto_accret eventually.

      if (emb(k,1) > 9.e-13) then
         eff(k,2) = min(1.,30. * (emb(k,1) - 9.e-13) ** .15)
      else
         eff(k,2) = 0.
      endif
   enddo
endif

!     3 = rr

if (jnmb(2) >= 1) then
   do k = k1(2),k2(2)

      if (rx(k,2) < rxmin(2)) cycle

! rain breakup (old)

!            dmr = dn(k,2) * gnu2
!            if (dmr .lt. .0006) then
!               eff(k,3) = 1.0
!            elseif (dmr .gt. .001446) then
!               eff(k,3) = -5.0
!            else
!               eff(k,3) = exp(2300. * (dmr - .0006))
!            endif

! rain breakup (new - temporary; eventually combine with autoconv/accret

      if (emb(k,2) < .113e-6) then
         eff(k,3) = 1.0
      elseif (emb(k,2) > .158e-5) then
         eff(k,3) = -5.0
      else
         eff(k,3) = 2. - exp(.1326e7 * (emb(k,2) - .113e-6))
      endif

   enddo

endif

!     4 = pp,ps,pa

if (jnmb(5) >= 1) then
   do k = k1(3),k2(3)
      if (abs(tx(k,3)+14.) <= 2.) then
         eff(k,4) = 1.4
      else
         eff(k,4) = min(.2,10. ** (.035 * tx(k,3) - .7))
      endif

   enddo

!     5 = ss,sa

   do k = k1(4),k2(4)
      if (abs(tx(k,4)+14.) <= 2.) then
         eff(k,5) = 1.4
      else
         eff(k,5) = min(.2,10. ** (.035 * tx(k,4) - .7))
      endif
   enddo

!     6 = aa

   do k = k1(5),k2(5)

      if (rx(k,5) < rxmin(5)) cycle

      if (abs(tx(k,5)+14.) <= 2.) then
         eff(k,6) = 1.4
      elseif (tx(k,5) >= -1.) then
         eff(k,6) = 1.
      else
         eff(k,6) = min(.2,10. ** (.035 * tx(k,5) - .7))
      endif

   enddo
endif

!     7 = pg,sg,ag,gg,gh

if (jnmb(6) >= 1) then
   do k = k1(6),k2(6)
      if (qr(k,6) > 0.) then
         eff(k,7) = 1.0
      else
         eff(k,7) = min(.2,10. ** (.035 * tx(k,6) - .7))
      endif
   enddo
endif

!     8 = ph,sh,ah,gh

if (jnmb(7) >= 1) then
   do k = k1(7),k2(7)

      if (rx(k,7) < rxmin(7)) cycle

      if (qr(k,7) > 0.) then
         eff(k,8) = 1.0
      else
         eff(k,8) = min(.2,10. ** (.035 * tx(k,7) - .7))
      endif

   enddo
endif

!     9 = cg,ch

if (jnmb(2) >= 1 .or. jnmb(3) >= 1) then
   do k = k1(1),k2(1)


! Rough fit from Pruppacher and Klett Fig. 14-11 p. 485:
!  close to curves for 142 and 305 microns.  Replace with auto_accret eventually.

      if (emb(k,1) > 3.4e-14) then
         eff(k,9) = min(1.,1426. * (emb(k,1) - 3.4e-14) ** .28)
      else
         eff(k,9) = 0.
      endif
   enddo
endif

!     10 = hh (trial)

if (jnmb(7) >= 1) then
   do k = k1(7),k2(7)
      eff(k,10) = max(0.,.1 + .005 * tx(k,7))
   enddo
endif

return
end subroutine effxy

!===============================================================================

!c      SUBROUTINE EFFAB(m1,MX,MY,EFF,DIAX,DIAY,TMPX,TMPY,TDEW)
!c      DIMENSION EFF(*),PT(*),DIAX(*),DIAY(*),tmpx(*),tmpy(*),tdew(*)
!
!c      IF (MX.EQ.6.AND.MY.EQ.6) THEN
!c        DO K=2,M1
!c          EFF(K)=MIN(.2,10.**(.035*(MAX(TMPX(K),TMPY(K))
!c     +          -273.15)-.7))
!c        ENDDO
!c      ELSE
!c        DO K=2,M1
!c          TMP=MAX(TMPX(K),TMPY(K))
!c          EFF(K)=MAX(MIN(1.,(TMP-273.06)*1.E10),
!c     +           MIN(10.**(.035*TMP-273.15)-.7),.2))
!c         ENDDO
!c         IF (MX.LE.5.AND.MY.LE.5) THEN
!c           DO K=2,M1
!c            EFF(K)=MAX(EFF(K),MIN(1.4,
!c     +        (MAX(0.,(2.-ABS(TMPX(k)-259.15)))*(TDEW(k)-TMPX(k))
!c     +        +MAX(0.,(2.-ABS(TMPY(k)-259.15)))*(TDEW(k)-TMPY(k)))
!c     +        *1.E10))
!c            ENDDO
!c         ENDIF
!c      ENDIF
!c      RETURN
!c      END

!===============================================================================

subroutine cols(mx,mc1,j1,j2  &
   ,jhcat,ict1,ict2,wct1,wct2,rx,cx,eff,colfac,enxfer)

use micro_coms, only: mza0, ncat, rxmin, ipairc, coltabc, neff
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mx
integer, intent(in) :: mc1
integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat (mza0,ncat)
integer, intent(in) :: ict1  (mza0,ncat)
integer, intent(in) :: ict2  (mza0,ncat)

real, intent(in) :: wct1 (mza0,ncat)
real, intent(in) :: wct2 (mza0,ncat)
real, intent(in) :: rx   (mza0,ncat)
real, intent(in) :: cx   (mza0,ncat)
real, intent(in) :: eff  (mza0,neff)

real, intent(in) :: colfac(mza0)

real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: ipc,k
real :: colnum,tabval

do k = j1,j2

   if (rx(k,mx) < rxmin(mx)) cycle

   ipc = ipairc(jhcat(k,mx),jhcat(k,mx))

   tabval  &
   = wct1(k,mx) ** 2              * coltabc(ict1(k,mx),ict1(k,mx),ipc)  &
   + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc)  &
   + wct2(k,mx) ** 2              * coltabc(ict2(k,mx),ict2(k,mx),ipc)

   colnum = colfac(k) * eff(k,mc1) * cx(k,mx) ** 2 * 10. ** (-tabval)
   enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(0.5 * cx(k,mx),colnum)

enddo
return
end subroutine cols

!===============================================================================

subroutine col3344(mx,mz,mc1,j1,j2  &
   ,jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac2,rxfer,qrxfer,enxfer)

use micro_coms, only: mza0, ncat, rxmin, ipairr, ipairc, jnmb, neff,  &
                      coltabc, coltabr
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mx
integer, intent(in) :: mz
integer, intent(in) :: mc1
integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat (mza0,ncat)
integer, intent(in) :: ict1  (mza0,ncat)
integer, intent(in) :: ict2  (mza0,ncat)

real, intent(in) :: wct1 (mza0,ncat)
real, intent(in) :: wct2 (mza0,ncat)
real, intent(in) :: rx   (mza0,ncat)
real, intent(in) :: cx   (mza0,ncat)
real, intent(in) :: qx   (mza0,ncat)
real, intent(in) :: eff  (mza0,neff)

real, intent(in) :: colfac2(mza0)

real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: k,ip,ipc
real :: c1,tabvalx,colamt,tabvaln,colnum

do k = j1,j2

   if (rx(k,mx) < rxmin(mx)) cycle
   
   ip = ipairr(jhcat(k,mx),jhcat(k,mx))
   ipc = ipairc(jhcat(k,mx),jhcat(k,mx))

   c1 = eff(k,mc1) * colfac2(k) * cx(k,mx) ** 2

   tabvalx  &
    = wct1(k,mx) ** 2              * coltabr(ict1(k,mx),ict1(k,mx),ip)  &
    + 2. * wct1(k,mx) * wct2(k,mx) * coltabr(ict1(k,mx),ict2(k,mx),ip)  &
    + wct2(k,mx) ** 2              * coltabr(ict2(k,mx),ict2(k,mx),ip)

   colamt = min(rx(k,mx),c1 * 10. ** (-tabvalx))
   rxfer(k,mx,mz) = rxfer(k,mx,mz) + colamt
   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + colamt * qx(k,mx)

   if (jnmb(mz) < 5) cycle

   tabvaln  &
      = wct1(k,mx) ** 2              * coltabc(ict1(k,mx),ict1(k,mx),ipc)  &
      + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc)  &
      + wct2(k,mx) ** 2              * coltabc(ict2(k,mx),ict2(k,mx),ipc)

      colnum = min(0.5 * cx(k,mx),c1 * 10. ** (-tabvaln))
      enxfer(k,mx,mz) = enxfer(k,mx,mz) + colnum
      enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnum

enddo
return
end subroutine col3344

!===============================================================================

subroutine col3443(j1,j2  &
   ,jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

use micro_coms, only: mza0, ncat, rxmin, neff, ipairc, ipairr, coltabc, coltabr
use misc_coms,  only: io6

implicit none

integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat (mza0,ncat)
integer, intent(in) :: ict1  (mza0,ncat)
integer, intent(in) :: ict2  (mza0,ncat)

real, intent(in) :: wct1  (mza0,ncat)
real, intent(in) :: wct2  (mza0,ncat)
real, intent(in) :: rx    (mza0,ncat)
real, intent(in) :: cx    (mza0,ncat)
real, intent(in) :: qx    (mza0,ncat)
real, intent(in) :: eff   (mza0,neff)

real, intent(in) :: colfac(mza0)

real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: k,jhcatx,jhcaty,ipxy,ipyx,ipc
real :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum

do k = j1,j2

   if (rx(k,3) < rxmin(3) .or. rx(k,4) < rxmin(4)) cycle

   jhcatx = jhcat(k,3)
   jhcaty = jhcat(k,4)
   ipxy = ipairr(jhcatx,jhcaty)
   ipyx = ipairr(jhcaty,jhcatx)
   ipc  = ipairc(jhcatx,jhcaty)

   c1 = eff(k,4) * colfac(k) * cx(k,3) * cx(k,4)

   tabvalx  &
     = wct1(k,3) * wct1(k,4) * coltabr (ict1(k,3),ict1(k,4),ipxy)  &
     + wct2(k,3) * wct1(k,4) * coltabr (ict2(k,3),ict1(k,4),ipxy)  &
     + wct1(k,3) * wct2(k,4) * coltabr (ict1(k,3),ict2(k,4),ipxy)  &
     + wct2(k,3) * wct2(k,4) * coltabr (ict2(k,3),ict2(k,4),ipxy)
   rcx = min(rx(k,3),c1 * 10. ** (-tabvalx))

   tabvaly  &
     = wct1(k,4) * wct1(k,3) * coltabr (ict1(k,4),ict1(k,3),ipyx)  &
     + wct2(k,4) * wct1(k,3) * coltabr (ict2(k,4),ict1(k,3),ipyx)  &
     + wct1(k,4) * wct2(k,3) * coltabr (ict1(k,4),ict2(k,3),ipyx)  &
     + wct2(k,4) * wct2(k,3) * coltabr (ict2(k,4),ict2(k,3),ipyx)
   rcy = min(rx(k,4),c1 * 10. ** (-tabvaly))

   rxfer(k,3,5) = rxfer(k,3,5) + rcx
   qrxfer(k,3,5) = qrxfer(k,3,5) + rcx * qx(k,3)

   rxfer(k,4,5) = rxfer(k,4,5) + rcy
   qrxfer(k,4,5) = qrxfer(k,4,5) + rcy * qx(k,4)

   tabvaln  &
       = wct1(k,3) * wct1(k,4) * coltabc (ict1(k,3),ict1(k,4),ipc)  &
       + wct2(k,3) * wct1(k,4) * coltabc (ict2(k,3),ict1(k,4),ipc)  &
       + wct1(k,3) * wct2(k,4) * coltabc (ict1(k,3),ict2(k,4),ipc)  &
       + wct2(k,3) * wct2(k,4) * coltabc (ict2(k,3),ict2(k,4),ipc)
   colnum = c1 * 10. ** (-tabvaln)

   if (cx(k,3) > cx(k,4)) then
      enxfer(k,4,5) = enxfer(k,4,5) + min(cx(k,4),colnum)
      enxfer(k,3,3) = enxfer(k,3,3) + min(cx(k,3),colnum)
   else
      enxfer(k,3,5) = enxfer(k,3,5) + min(cx(k,3),colnum)
      enxfer(k,4,4) = enxfer(k,4,4) + min(cx(k,4),colnum)
   endif

! also loss for aerosol

enddo
return
end subroutine col3443

!===============================================================================

subroutine col1(mx,my,mz,mc4,j1,j2  &
   ,jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

use micro_coms, only: mza0, ncat, neff, rxmin, ipairr, ipairc, jnmb,  &
                      coltabc, coltabr
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mx
integer, intent(in) :: my
integer, intent(in) :: mz
integer, intent(in) :: mc4
integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat (mza0,ncat)
integer, intent(in) :: ict1  (mza0,ncat)
integer, intent(in) :: ict2  (mza0,ncat)

real, intent(in) :: wct1  (mza0,ncat)
real, intent(in) :: wct2  (mza0,ncat)
real, intent(in) :: rx    (mza0,ncat)
real, intent(in) :: cx    (mza0,ncat)
real, intent(in) :: qx    (mza0,ncat)
real, intent(in) :: eff   (mza0,neff)

real, intent(in) :: colfac(mza0)

real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: k,ipxy,ipc
real :: c1,tabvalx,rcx,tabvaln,colnum

do k = j1,j2

   if (rx(k,mx) < rxmin(mx) .or. rx(k,my) < rxmin(my)) cycle

   ipxy = ipairr(jhcat(k,mx),jhcat(k,my))
   ipc  = ipairc(jhcat(k,mx),jhcat(k,my))

   c1 = eff(k,mc4) * colfac(k) * cx(k,mx) * cx(k,my)

   tabvalx  &
     = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
     + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
     + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
     + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)

   rcx = min(rx(k,mx),c1 * 10. ** (-tabvalx))
   rxfer(k,mx,mz) = rxfer(k,mx,mz) + rcx
   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + rcx * qx(k,mx)

   if (jnmb(mx) < 5) cycle

   tabvaln  &
     = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
     + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
     + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
     + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

   colnum = c1 * 10. ** (-tabvaln)
   enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(colnum,cx(k,mx))

! also loss for aerosol

enddo
return
end subroutine col1

!===============================================================================

subroutine col2(my,mz,mc2,j1,j2  &
   ,jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac  &
   ,rxfer,qrxfer,enxfer,dtl0)

use micro_coms, only: mza0, ncat, rxmin, ipairr, ipairc, coltabr, coltabc,  &
                      sipfac, pwmasi, emb1, gamsip13, gamsip24, emb0, neff
use misc_coms,  only: io6

implicit none

integer, intent(in) :: my
integer, intent(in) :: mz
integer, intent(in) :: mc2
integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat (mza0,ncat)
integer, intent(in) :: ict1  (mza0,ncat)
integer, intent(in) :: ict2  (mza0,ncat)

real, intent(in) :: wct1  (mza0,ncat)
real, intent(in) :: wct2  (mza0,ncat)
real, intent(in) :: rx    (mza0,ncat)
real, intent(in) :: cx    (mza0,ncat)
real, intent(in) :: qx    (mza0,ncat)
real, intent(in) :: emb   (mza0,ncat)
real, intent(in) :: eff   (mza0,neff)

real, intent(in) :: colfac(mza0)

real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

real, intent(in) :: dtl0

integer :: k,jhcatx,jhcaty,ipxy,ipyx,ipc,it
real :: c1,c2,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum0,colnum,rcoal  &
       ,qrcx,qrcy,qrcoal,qcoal,fracliq,tcoal,coalliq,coalice,area,cn13,cn24  &
       ,sip,rsip,qrsip,rfinlz,xtoz

real, dimension(15) ::  alpha,beta
!            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
data alpha /00.,00.,00., 1., 1., 1., 1.,00.,00.,00.,00., 1., 1., 1., 1./
data beta  /00.,00.,00.,1.5,1.1,0.0,0.0,00.,00.,00.,00.,1.2,1.1,1.1,1.3/

do k = j1,j2

   if (rx(k,1) < rxmin(1) .or. rx(k,my) < rxmin(my)) cycle

   jhcatx = jhcat(k,1)
   jhcaty = jhcat(k,my)
   ipxy = ipairr(jhcatx,jhcaty)
   ipyx = ipairr(jhcaty,jhcatx)
   ipc  = ipairc(jhcatx,jhcaty)

   c2 = colfac(k) * cx(k,1) * cx(k,my)
   c1 = eff(k,mc2) * c2

   tabvalx = wct1(k,1) * wct1(k,my) * coltabr (ict1(k,1),ict1(k,my),ipxy)  &
           + wct2(k,1) * wct1(k,my) * coltabr (ict2(k,1),ict1(k,my),ipxy)  &
           + wct1(k,1) * wct2(k,my) * coltabr (ict1(k,1),ict2(k,my),ipxy)  &
           + wct2(k,1) * wct2(k,my) * coltabr (ict2(k,1),ict2(k,my),ipxy)

   rcx = min(rx(k,1),c1 * 10. ** (-tabvalx))

   tabvaly = wct1(k,my) * wct1(k,1) * coltabr (ict1(k,my),ict1(k,1),ipyx)  &
           + wct2(k,my) * wct1(k,1) * coltabr (ict2(k,my),ict1(k,1),ipyx)  &
           + wct1(k,my) * wct2(k,1) * coltabr (ict1(k,my),ict2(k,1),ipyx)  &
           + wct2(k,my) * wct2(k,1) * coltabr (ict2(k,my),ict2(k,1),ipyx)

   rcy = min(rx(k,my),c1 * 10. ** (-tabvaly))

   tabvaln = wct1(k,1) * wct1(k,my) * coltabc (ict1(k,1),ict1(k,my),ipc)  &
           + wct2(k,1) * wct1(k,my) * coltabc (ict2(k,1),ict1(k,my),ipc)  &
           + wct1(k,1) * wct2(k,my) * coltabc (ict1(k,1),ict2(k,my),ipc)  &
           + wct2(k,1) * wct2(k,my) * coltabc (ict2(k,1),ict2(k,my),ipc)

   colnum0 = c2 * 10. ** (-tabvaln)
   colnum = colnum0 * eff(k,mc2)

   rcoal = rcx + rcy
   qrcx = rcx * qx(k,1)
   qrcy = rcy * qx(k,my)
   qrcoal = qrcx + qrcy
   qcoal = qrcoal / (1.e-13 + rcoal)

   call qtc(qcoal,tcoal,fracliq)

   coalliq = rcoal * fracliq
   coalice = rcoal - coalliq

! secondary ice production: cn24 is the number fraction of collected cloud
! droplets larger than 24 microns and is obtained from an incomplete gamma
! function table.  cn13 is the fraction of collected cloud droplets
! smaller than 13 microns.  area is cross section area of collecting ice
! per m^3 of atmospheric volume.

   if (tcoal > -8. .and. tcoal < -3.) then

      area = cx(k,my) * sipfac(jhcaty) * emb(k,my) ** (2.*pwmasi(jhcaty)) ! remd rhoa
      it = nint(emb(k,1) / emb1(1) * 5000.)
      cn13 = colnum * gamsip13(it) / (area * dtl0)
      cn24 = min(cx(k,1),colnum0) * gamsip24(it)
      sip = 9.1e-10 * cn24 * cn13 ** .93           ! has units of #/m^3
      if (tcoal < -5.) then
         sip = 0.33333 * (tcoal + 8.) * sip
      else
         sip = -0.5 * (tcoal + 3.) * sip
      endif

      rsip = sip * emb0(3)           ! has units of kg/m^3
      qrsip = qcoal * rsip           ! has units of J/m^3

      rcoal = rcoal - rsip
      qrcoal = qrcoal - qrsip

      enxfer(k,1,3) = enxfer(k,1,3) + sip
      rxfer(k,1,3)  = rxfer(k,1,3)  + rsip
      qrxfer(k,1,3) = qrxfer(k,1,3) + qrsip

   endif

! ALWAYS NEED (ALPHA + BETA) .GE. 1 but in the (rare) case that
! fracliq may be a little larger than fracx due to collected
! liquid being above 0C, need (ALPHA + BETA) to be at least 1.1
! or 1.2, or need ALPHA itself to be at least 1.0.

   rfinlz = min(rcoal,alpha(jhcaty) * coalliq + beta(jhcaty) * rcx)

   xtoz = min(rcx,rfinlz)

   rxfer(k,1,mz) = rxfer(k,1,mz) + xtoz
   rxfer(k,1,my) = rxfer(k,1,my) + rcx - xtoz
   if (my /= mz) rxfer(k,my,mz) = rxfer(k,my,mz) + rfinlz - xtoz

   qrxfer(k,1,mz) = qrxfer(k,1,mz) + qx(k,1) * xtoz
   qrxfer(k,1,my) = qrxfer(k,1,my) + qx(k,1) * (rcx - xtoz)
   if (my /= mz) qrxfer(k,my,mz) = qrxfer(k,my,mz) + qx(k,my) * (rfinlz - xtoz)

   enxfer(k,1,1) = enxfer(k,1,1) + min(colnum,cx(k,1))
   if (my /= mz) enxfer(k,my,mz) = enxfer(k,my,mz)  &
      + (rfinlz - xtoz) * min(colnum,cx(k,my)) / max(1.e-20,rcy)
      
! BUT NEED TO CHANGE THE ABOVE FOR 177 COLLECTION BECAUSE X = Y

! also include loss of aerosol

enddo
return
end subroutine col2

!===============================================================================

subroutine col3(my,mz,j1,j2  &
   ,jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

use micro_coms, only: mza0, ncat, rxmin, ipairr, ipairc, jnmb, neff,  &
                      coltabr, coltabc
use misc_coms,  only: io6

implicit none

integer, intent(in) :: my
integer, intent(in) :: mz
integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat (mza0,ncat)
integer, intent(in) :: ict1  (mza0,ncat)
integer, intent(in) :: ict2  (mza0,ncat)

real, intent(in) :: wct1  (mza0,ncat)
real, intent(in) :: wct2  (mza0,ncat)
real, intent(in) :: rx    (mza0,ncat)
real, intent(in) :: cx    (mza0,ncat)
real, intent(in) :: qx    (mza0,ncat)
real, intent(in) :: eff   (mza0,neff)

real, intent(in) :: colfac(mza0)

real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: k,ipxy,ipyx,ipc,jhcaty
real :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum,colnumx,colnumy,coalnum  &
       ,rcoal,qrcx,qrcy,qrcoal,qcoal,fracliq,coalliq,coalice,xtoz  &
       ,rfinlz,tcoal,cfinlz
real, dimension(15) :: alpha,beta

!            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
data alpha /00.,00., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1./
data beta  /00.,00., 2., 2., 2., 1., 0., 2., 2., 2., 2., 2., 2., 2., 2./

do k = j1,j2

   if (rx(k,2) < rxmin(2) .or. rx(k,my) < rxmin(my)) cycle

   jhcaty = jhcat(k,my)
   ipxy = ipairr(jhcat(k,2),jhcaty)
   ipyx = ipairr(jhcaty,jhcat(k,2))
   ipc  = ipairc(jhcat(k,2),jhcaty)
   c1 = eff(k,1) * colfac(k) * cx(k,2) * cx(k,my)

   tabvalx = wct1(k,2) * wct1(k,my) * coltabr (ict1(k,2),ict1(k,my),ipxy)  &
           + wct2(k,2) * wct1(k,my) * coltabr (ict2(k,2),ict1(k,my),ipxy)  &
           + wct1(k,2) * wct2(k,my) * coltabr (ict1(k,2),ict2(k,my),ipxy)  &
           + wct2(k,2) * wct2(k,my) * coltabr (ict2(k,2),ict2(k,my),ipxy)

   rcx = min(rx(k,2),c1 * 10. ** (-tabvalx))

   tabvaly = wct1(k,my) * wct1(k,2) * coltabr (ict1(k,my),ict1(k,2),ipyx)  &
           + wct2(k,my) * wct1(k,2) * coltabr (ict2(k,my),ict1(k,2),ipyx)  &
           + wct1(k,my) * wct2(k,2) * coltabr (ict1(k,my),ict2(k,2),ipyx)  &
           + wct2(k,my) * wct2(k,2) * coltabr (ict2(k,my),ict2(k,2),ipyx)

   rcy = min(rx(k,my),c1 * 10. ** (-tabvaly))

   if (jnmb(2) >= 5) then
      tabvaln = wct1(k,2) * wct1(k,my) * coltabc (ict1(k,2),ict1(k,my),ipc)  &
              + wct2(k,2) * wct1(k,my) * coltabc (ict2(k,2),ict1(k,my),ipc)  &
              + wct1(k,2) * wct2(k,my) * coltabc (ict1(k,2),ict2(k,my),ipc)  &
              + wct2(k,2) * wct2(k,my) * coltabc (ict2(k,2),ict2(k,my),ipc)

      colnum = c1 * 10. ** (-tabvaln)
      colnumx = min(cx(k,2),colnum)
      colnumy = min(cx(k,my),colnum)
      coalnum = min(colnumx,colnumy)
   endif

   rcoal = rcx + rcy
   qrcx = rcx * qx(k,2)
   qrcy = rcy * qx(k,my)
   qrcoal = qrcx + qrcy
   qcoal = qrcoal / (1.e-20 + rcoal)

   call qtc(qcoal,tcoal,fracliq)

   coalliq = rcoal * fracliq
   coalice = rcoal - coalliq

   if (fracliq >= .99) then

      rxfer(k,my,2) = rxfer(k,my,2) + rcy
      qrxfer(k,my,2) = qrxfer(k,my,2) + qrcy
      if (jnmb(2) >= 5) enxfer(k,my,my) = enxfer(k,my,my) + colnumy

   else

      rfinlz = min(rcoal, alpha(jhcaty) * coalliq + beta(jhcaty) * rcx)

      xtoz = min(rcx,rfinlz)

      rxfer(k,2,mz) = rxfer(k,2,mz) + xtoz
      rxfer(k,2,my) = rxfer(k,2,my) + rcx - xtoz
      if (my /= mz) rxfer(k,my,mz) = rxfer(k,my,mz) + rfinlz - xtoz

! NEED TO USE QCOAL TO TRANSFER Q?

      qrxfer(k,2,mz) = qrxfer(k,2,mz) + qx(k,2) * xtoz
      qrxfer(k,2,my) = qrxfer(k,2,my) + qx(k,2) * (rcx - xtoz)
      if (my /= mz) qrxfer(k,my,mz) = qrxfer(k,my,mz)  &
         + qx(k,my) * (rfinlz - xtoz)

      if (jnmb(2) >= 5) then
         if (my == mz) then
            enxfer(k,2,2) = enxfer(k,2,2) + colnumx
         elseif (colnumy >= colnumx) then
            cfinlz = coalnum * rfinlz / (rcoal + 1.e-20)
            enxfer(k,2,mz) = enxfer(k,2,mz) + cfinlz
            enxfer(k,2,2) = enxfer(k,2,2) + colnumx - cfinlz
            enxfer(k,my,my) = enxfer(k,my,my) + colnumy
         else
            cfinlz = coalnum * rfinlz / (rcoal + 1.e-20)
            enxfer(k,my,mz) = enxfer(k,my,mz) + cfinlz
            enxfer(k,2,2) = enxfer(k,2,2) + colnumx
            enxfer(k,my,my) = enxfer(k,my,my) + colnumy - cfinlz
         endif
      endif

   endif

enddo

! also include loss of aerosol

return
end subroutine col3

!===============================================================================

subroutine colxfers(k1,k2,rx,cx,qr,rxfer,qrxfer,enxfer)

use micro_coms, only: mza0, ncat, jnmb
use misc_coms,  only: io6

implicit none

integer, intent(in) :: k1(10)
integer, intent(in) :: k2(10)

real, intent(inout) :: rx    (mza0,ncat)
real, intent(inout) :: cx    (mza0,ncat)
real, intent(inout) :: qr    (mza0,ncat)
real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: lcat,kd1,kd2,jcat,k

real :: rloss (mza0) ! automatic array
real :: enloss(mza0) ! automatic array

!  All rxfer values are nonnegative.

do lcat = 1,7
   if (jnmb(lcat) >= 1) then
      kd1 = k1(lcat)
      kd2 = k2(lcat)

      do k = kd1,kd2
         rloss(k) = 0.
         enloss(k) = 0.
      enddo

      do jcat = 1,7
! change this to include enxfer of the same categories
         if (jnmb(jcat) >= 1) then
            if (lcat /= jcat) then
               do k = kd1,kd2
                  rloss(k) = rloss(k) + rxfer(k,lcat,jcat)
               enddo
            endif
            do k = kd1,kd2
               enloss(k) = enloss(k) + enxfer(k,lcat,jcat)  ! #/m^3
            enddo
         endif
      enddo

      do k = kd1,kd2
         rloss(k) = min(1.,rx(k,lcat) / max(1.e-20,rloss(k)))
         enloss(k) = min(1.,cx(k,lcat) / max(1.e-10,enloss(k)))   ! dimensionless
      enddo

      do jcat = 1,7
         if (jnmb(jcat) >= 1) then
            if (lcat /= jcat) then
               do k = kd1,kd2
                  rxfer(k,lcat,jcat) = rxfer(k,lcat,jcat)*rloss(k)
                  qrxfer(k,lcat,jcat)=qrxfer(k,lcat,jcat)*rloss(k)
               enddo
            endif
            do k = kd1,kd2
               enxfer(k,lcat,jcat) = enxfer(k,lcat,jcat)*enloss(k)
            enddo
         endif
      enddo
   endif
enddo

do lcat = 1,7

   if (jnmb(lcat) >= 1) then

      kd1 = k1(lcat)
      kd2 = k2(lcat)

      do jcat = 1,7
         if (jnmb(jcat) >= 1 .and. lcat /= jcat) then
            do k = kd1,kd2
               rx(k,lcat) = rx(k,lcat) - rxfer(k,lcat,jcat)
               rx(k,jcat) = rx(k,jcat) + rxfer(k,lcat,jcat)
               qr(k,lcat) = qr(k,lcat) - qrxfer(k,lcat,jcat)
               qr(k,jcat) = qr(k,jcat) + qrxfer(k,lcat,jcat)
               cx(k,lcat) = cx(k,lcat) - enxfer(k,lcat,jcat)
               cx(k,jcat) = cx(k,jcat) + enxfer(k,lcat,jcat)
            enddo
         endif
      enddo

      if (jnmb(lcat) >= 5) then
         do k = kd1,kd2
            cx(k,lcat) = cx(k,lcat) - enxfer(k,lcat,lcat)
         enddo
      endif

   endif
enddo
return
end subroutine colxfers


