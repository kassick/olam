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
subroutine scalar_transport(umarusc,wmarwsc,rho_old)

use mem_ijtabs, only: istp, jtab_w, mrl_endl
use mem_grid,   only: mza, mua, mwa, zt, zm, dzim
use misc_coms,  only: io6

implicit none

real, intent(in) :: umarusc(mza,mua)
real, intent(in) :: wmarwsc(mza,mwa)
real, intent(in) :: rho_old(mza,mwa)

integer :: j,iw,k,mrl

! Automatic arrays:

real :: zwt1(mza)
real :: zwt2(mza)

! Compute vertical advective weights for a stretched grid

do k = 1,mza-1
   zwt1(k) = (zt(k+1) - zm(k)) * dzim(k)
   zwt2(k) =  (zm(k) - zt(k))  * dzim(k)
enddo

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
do j = 1,jtab_w(26)%jend(mrl); iw = jtab_w(26)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call scalar_transport0(iw,umarusc,wmarwsc,zwt1,zwt2,rho_old)

enddo     ! end horizontal loop over W/T points
endif
call rsub('W',26)

return
end subroutine scalar_transport

!===============================================================================

subroutine scalar_transport0(iw,umarusc,wmarwsc,zwt1,zwt2,rho_old)

use var_tables, only: num_scalar, scalar_tab
use mem_ijtabs, only: itab_w
use misc_coms,  only: io6, dtlm, cnum_sclr
use mem_turb,   only: vkh, hkm, sxfer_rk
use mem_grid,   only: mza, mua, mwa, lpw, lsw,  &
                      dniu, volt, aru, arw, volwi, dzim, volti

implicit none

integer, intent(in) :: iw

real, intent(in) :: umarusc(mza,mua)
real, intent(in) :: wmarwsc(mza,mwa)
real, intent(in) :: rho_old(mza,mwa)

real, intent(in) :: zwt1(mza)
real, intent(in) :: zwt2(mza)

real, pointer :: scp(:,:)
real, pointer :: sct(:,:)

integer :: n,k,ka,ks,iu1,iu2,iu3,iw1,iw2,iw3
real :: dtl,dtl2,dtli,hdtli,hdniu1,hdniu2,hdniu3,diru1,diru2,diru3

! Automatic arrays:

real, dimension(mza) :: akodz,dtomass,vctr5,vctr6,vctr7,vctr8,vctr9
real, dimension(mza) :: fadvz
real, dimension(mza) :: flxu1,flxu2,flxu3
real, dimension(mza) :: tmass,tmass1,tmass2,tmass3
real, dimension(mza) :: umass1,umass2,umass3
real, dimension(mza) :: evac,evac1,evac2,evac3
real, dimension(mza) :: cnum_u1,cnum_u2,cnum_u3
real, dimension(mza) :: akodn1,akodn2,akodn3,cnum_w
real, dimension(mza) :: flx1,flx2,flx3
real, dimension(mza) :: del_scp

iu1 = itab_w(iw)%iu1; iu2 = itab_w(iw)%iu2; iu3 = itab_w(iw)%iu3
iw1 = itab_w(iw)%iw1; iw2 = itab_w(iw)%iw2; iw3 = itab_w(iw)%iw3
diru1 = itab_w(iw)%diru1; diru2 = itab_w(iw)%diru2; diru3 = itab_w(iw)%diru3

! Initial computations for this column that need not be repeated for
! each transported scalar field

dtl = dtlm(itab_w(iw)%mrlw)
dtl2 = dtl * 2.
dtli = 1. / dtl
hdtli = .5 * dtli
ka = lpw(iw)
   
hdniu1 = .5 * dniu(iu1)
hdniu2 = .5 * dniu(iu2)
hdniu3 = .5 * dniu(iu3)

! Vertical loop over T levels 

do k = ka,mza-1

! Horizontal advective mass flux from 3 UM components

   flxu1(k) = diru1 * umarusc(k,iu1)
   flxu2(k) = diru2 * umarusc(k,iu2)
   flxu3(k) = diru3 * umarusc(k,iu3)

! mass in control volumes of 4 T cells

   tmass(k)  = rho_old(k,iw)  * volt(k,iw)
   tmass1(k) = rho_old(k,iw1) * volt(k,iw1)
   tmass2(k) = rho_old(k,iw2) * volt(k,iw2)
   tmass3(k) = rho_old(k,iw3) * volt(k,iw3)

! mass in control volumes of 3 UM cells

   umass1(k) = tmass(k) + tmass1(k)
   umass2(k) = tmass(k) + tmass2(k)
   umass3(k) = tmass(k) + tmass3(k)

! limiting mass evacuation rate in 4 T cells for advective flux limiters
 
   evac(k)  = hdtli * tmass(k)
   evac1(k) = hdtli * tmass1(k)
   evac2(k) = hdtli * tmass2(k)
   evac3(k) = hdtli * tmass3(k)

! Scalar horizontal advective courant numbers from 3 UM components
! (double time in numerator compensates for double mass in denominator
! since flux used for T cell) 

   cnum_u1(k) = flxu1(k) * dtl2 / umass1(k)
   cnum_u2(k) = flxu2(k) * dtl2 / umass2(k)
   cnum_u3(k) = flxu3(k) * dtl2 / umass3(k)

   if (cnum_sclr > abs(cnum_u1(k))) cnum_u1(k) = sign(cnum_sclr,cnum_u1(k)) 
   if (cnum_sclr > abs(cnum_u2(k))) cnum_u2(k) = sign(cnum_sclr,cnum_u2(k)) 
   if (cnum_sclr > abs(cnum_u3(k))) cnum_u3(k) = sign(cnum_sclr,cnum_u3(k)) 

! Scalar horizontal turbulent flux coefficients

   akodn1(k) = aru(k,iu1) * (hkm(k,iw1) + hkm(k,iw)) * hdniu1
   akodn2(k) = aru(k,iu2) * (hkm(k,iw2) + hkm(k,iw)) * hdniu2
   akodn3(k) = aru(k,iu3) * (hkm(k,iw3) + hkm(k,iw)) * hdniu3

   dtomass(k) = dtl / tmass(k)

enddo

! Vertical loop over W levels

do k = ka,mza-2

! Prepare for vertical advection:  half vertical courant number for scalars

   cnum_w(k) = wmarwsc(k,iw) * dtl * volwi(k,iw)  &
      / (rho_old(k,iw) + rho_old(k+1,iw))

   if (.5*cnum_sclr > abs(cnum_w(k))) then
      cnum_w(k) = sign(.5*cnum_sclr,cnum_w(k))
   endif

! Prepare for vertical diffusion - Fill tri-diagonal matrix coefficients

   akodz(k) = arw(k,iw) * vkh(k,iw) * dzim(k)
   vctr5(k) = - akodz(k) * dtomass(k)
   vctr7(k) = - akodz(k) * dtomass(k+1)
   vctr6(k) = 1. - vctr5(k) - vctr7(k)
enddo
   
! Loop over scalars here

do n = 1,num_scalar

! Point SCP and SCT to scalar table arrays

   scp => scalar_tab(n)%var_p
   sct => scalar_tab(n)%var_t

! Vertical loop over W levels: 
! Fill r.h.s. for tri-diagonal matrix eqn for vertical diffusion

   do k = ka,mza-2
      vctr8(k) = akodz(k) * (scp(k,iw) - scp(k+1,iw))
   enddo

! Special case for total water sh_w:  apply surface vapor flux

   if (scalar_tab(n)%name == 'SH_W') then

! Vertical loop over T levels that are adjacent to surface

      do ks = 1,lsw(iw)
         k = ka + ks - 1

! Apply surface vapor xfer [kg_vap] directly to SCT [kg_vap / (m^3 s)]

         sct(k,iw) = sct(k,iw) + dtli * volti(k,iw) * sxfer_rk(ks,iw)

! Change in SCP from surface xfer

         del_scp(k) = sxfer_rk(ks,iw) / (rho_old(k,iw) * volt(k,iw))

! Zero out sxfer_rk(ks,iw) now that it has been transferred to the atm

         sxfer_rk(ks,iw) = 0.  

      enddo

! Lowest T level that is not adjacent to surface

      del_scp(ka+lsw(iw)) = 0.

! Vertical loop over W levels that are adjacent to surface

      do ks = 1,lsw(iw)
         k = ka + ks - 1

! Change in vctr8 from surface vapor xfer

         vctr8(k) = vctr8(k) + akodz(k) * (del_scp(k) - del_scp(k+1))
      enddo

   endif

! Solve tri-diagonal matrix equation

   call tridiffo(mza,ka,mza-2,vctr5,vctr6,vctr7,vctr8,vctr9)

! Set bottom and top vertical internal turbulent fluxes to zero

   vctr9(ka-1) = 0.
   vctr9(mza-1) = 0.

! Set bottom and top vertical advective fluxes to zero

   fadvz(ka-1) = 0.
   fadvz(mza-1) = 0.

! Vertical loop over W levels

   do k = ka,mza-2

! Compute scalar vertical advective flux

      fadvz(k) = wmarwsc(k,iw)      &
         * (zwt1(k) * scp(k,iw) + zwt2(k) * scp(k+1,iw)  &
         +  cnum_w(k) * (scp(k,iw) - scp(k+1,iw)))

! Modify scalar vertical advective flux if necessary to retain positive-definiteness of advected scalar field

      if (wmarwsc(k,iw) > 0.) then
         if (fadvz(k) > evac(k) * scp(k,iw))  &
             fadvz(k) = wmarwsc(k,iw) * scp(k,iw)
      else
         if (-fadvz(k) > evac(k+1) * scp(k+1,iw))  &
              fadvz(k) = wmarwsc(k,iw) * scp(k+1,iw)
      endif

   enddo

! Vertical loop over T levels

   do k = ka,mza-1

! Compute scalar horizontal advective fluxes

      flx1(k) = flxu1(k) * .5 *  &
         (scp(k,iw1) + scp(k,iw) + cnum_u1(k) * (scp(k,iw1) - scp(k,iw)))

      flx2(k) = flxu2(k) * .5 *  &
         (scp(k,iw2) + scp(k,iw) + cnum_u2(k) * (scp(k,iw2) - scp(k,iw)))

      flx3(k) = flxu3(k) * .5 *  &
         (scp(k,iw3) + scp(k,iw) + cnum_u3(k) * (scp(k,iw3) - scp(k,iw)))

! Modify scalar horizontal advective fluxes if necessary to retain
! positive-definiteness of advected scalar field

      if (flxu1(k) > 0.) then
         if (flx1(k) > evac1(k) * scp(k,iw1))  &
             flx1(k) = flxu1(k) * scp(k,iw1)
      else
         if (-flx1(k) > evac(k) * scp(k,iw))   &
              flx1(k) = flxu1(k) * scp(k,iw)
      endif

      if (flxu2(k) > 0.) then
         if (flx2(k) > evac2(k) * scp(k,iw2))  &
             flx2(k) = flxu2(k) * scp(k,iw2)
      else
         if (-flx2(k) > evac(k) * scp(k,iw))   &
              flx2(k) = flxu2(k) * scp(k,iw)
      endif

      if (flxu3(k) > 0.) then
         if (flx3(k) > evac3(k) * scp(k,iw3))  &
             flx3(k) = flxu3(k) * scp(k,iw3)
      else
         if (-flx3(k) > evac(k) * scp(k,iw))   &
              flx3(k) = flxu3(k) * scp(k,iw)
      endif

! Add contributions to scalar tendency from horizontal and vertical 
! advection and diffusion

      sct(k,iw) = sct(k,iw) + volti(k,iw) * (     &

            flx1(k) + flx2(k) + flx3(k)           & ! horizontal advection

         +  akodn1(k) * (scp(k,iw1) - scp(k,iw))  & ! horiz diff across ARU1
         +  akodn2(k) * (scp(k,iw2) - scp(k,iw))  & ! horiz diff across ARU2
         +  akodn3(k) * (scp(k,iw3) - scp(k,iw))  & ! horiz diff across ARU3

         +  fadvz(k-1) - fadvz(k)                 & ! vertical advection
         +  vctr9(k-1) - vctr9(k)                 ) ! vertical diffusion

   enddo
      
enddo  ! end loop over scalars

return
end subroutine scalar_transport0
