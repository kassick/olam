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
subroutine fldsisan(iaction, o_rho, o_theta, o_shv, o_uzonal, o_umerid)

use mem_nudge,   only: nudflag, nudrad, mnudp,  &
                       rho_obsp, theta_obsp, shw_obsp, uzonal_obsp,   &
                       umerid_obsp, rho_obsf, theta_obsf, shw_obsf,   &
                       uzonal_obsf, umerid_obsf
use mem_basic,   only: umc, ump, uc, thil, sh_w, sh_v, wmc, wc,  &
                       theta, rho, press
use mem_grid,   only:  mwa, nza, mza, mwa, mua, zm, zt, lpw, aru, volt,  &
                       unx, uny, unz, xeu, yeu, zeu
use misc_coms,  only: io6, deltax, deltaz, dzrat, dzmax, iparallel, runtype
use mem_micro,   only: sh_c
use micro_coms,  only: level
use mem_ijtabs,  only: jtab_u, jtab_w, itab_u, itab_w
use consts_coms, only: pc1, rdry, rvap,cpocv,erad

implicit none

integer, intent(in)  :: iaction

integer :: j,iw,k,iu,iw1,iw2,iup,ng,iuv,ngrd,inudp,inudp1

real :: rcloud,qlatu,qlonu,dummy,qlatv,qlonv,temp,rvls  &
       ,cosuv,sinuv,qlatuv,qlonuv,uu,vv,angle
real :: alph_p
real :: volnudpi

!real :: ug,vg,raxis,uvgr,uvgx,uvgy,uvgz

! Automatic arrays

real(kind=8), intent(in) :: o_rho   (mza,mwa)
real,         intent(in) :: o_theta (mza,mwa)
real,         intent(in) :: o_shv   (mza,mwa)
real,         intent(in) :: o_uzonal(mza,mwa)
real,         intent(in) :: o_umerid(mza,mwa)

real :: volnudp(mza,mnudp)

real :: ug(mza),vg(mza),raxis(mza),uvgr(mza),uvgx(mza),uvgy(mza),uvgz(mza)

! Check whether initializing model

if (iaction == 0 .and. runtype == 'INITIAL') then

! If initializing the model, fill the main model arrays
! and initialize related arrays

   call psub()
!----------------------------------------------------------------------
   do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)
!---------------------------------------------------------------------
   call qsub('W',iw)
      do k = 1,mza

         rho  (k,iw) = o_rho(k,iw)
         theta(k,iw) = o_theta(k,iw)
         sh_v (k,iw) = o_shv(k,iw)

         thil(k,iw) = theta(k,iw)
         sh_w(k,iw) = sh_v(k,iw)   ! no condensate considered here
         if (level > 1) then
            sh_c(k,iw) = 0.        ! no condensate considered here
         endif

! alph_p is like alpha_press in thiltend_long except that the (theta/thil) 
! factor is excluded here (it's equal to 1)

         alph_p = pc1 * (((1. - sh_w(k,iw)) * rdry + sh_v(k,iw) * rvap))  &
                ** cpocv

         press(k,iw) = alph_p * (rho(k,iw) * thil(k,iw)) ** cpocv

         wc(k,iw) = 0.
         wmc(k,iw) = 0.
         
      enddo
   enddo
   call rsub('Wb',7)

   if (iparallel == 1) then
      call mpi_send_w('I',wc)  ! Send W group; wc used as dummy array here
      call mpi_recv_w('I',wc)  ! Recv W group; wc used as dummy array here
   endif

! UMC, UC

   call psub()
!----------------------------------------------------------------------
   do j = 1,jtab_u(7)%jend(1); iu = jtab_u(7)%iu(j)
      iw1 = itab_u(iu)%iw1; iw2 = itab_u(iu)%iw2
!----------------------------------------------------------------------
   call qsub('U',iu)

      if (iw1 < 2) write(io6,*) 'iu,iw1,iw2',iu,iw1,iw2
      if (iw2 < 2) write(io6,*) 'iu,iw2,iw1',iu,iw2,iw1

      if (iw1 < 2) iw1 = iw2
      if (iw2 < 2) iw2 = iw1

! Average winds to U point and rotate at U point

      do k = 1,mza
         ug(k) = .5 * (o_uzonal(k,iw1) + o_uzonal(k,iw2))
         vg(k) = .5 * (o_umerid(k,iw1) + o_umerid(k,iw2))

         raxis(k) = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)  ! dist from earth axis
         uvgr(k) = -vg(k) * zeu(iu) / erad  ! radially outward from axis

         uvgx(k) = (-ug(k) * yeu(iu) + uvgr(k) * xeu(iu)) / raxis(k) 
         uvgy(k) = ( ug(k) * xeu(iu) + uvgr(k) * yeu(iu)) / raxis(k) 
         uvgz(k) =   vg(k) * raxis(k) / erad 

         uc(k,iu) = uvgx(k) * unx(iu) + uvgy(k) * uny(iu) + uvgz(k) * unz(iu)
         umc(k,iu) = uc(k,iu) * .5 * (rho(k,iw1) + rho(k,iw2))
      enddo

   enddo
   
   call rsub('Ub',7)

! Set UMC, UC to zero wherever ARU = 0.

   call psub()
!----------------------------------------------------------------------
   do iu = 2,mua
!----------------------------------------------------------------------
   call qsub('U',iu) !QQQQQ
      do k = 1,mza-1
         if (aru(k,iu) < 1.e-9) then
            umc(k,iu) = 0.
            uc(k,iu)  = 0.
         endif
      enddo
   enddo
   call rsub('Uc',0)

   if (iparallel == 1) then
      call mpi_send_u('I')  ! Send U group
      call mpi_recv_u('I')  ! Recv U group
   endif

! Set UMP to UMC

   ump(:,:) = umc(:,:)

! Print out initial state column from column 2

   iw = 2

write(io6,*)' '
write(io6,*)'========================================================================'
write(io6,*)'                    OLAM INITIAL STATE COLUMN (vari)'
write(io6,*)'========================================================================'
write(io6,*)'   zm(m)    k     zt(m)   press(Pa)   rho(kg/m3)   theta(K)   sh_w(g/kg)'
write(io6,*)'========================================================================'
write(io6,*)' '

   do k = mza-1,2,-1
      write(io6, '(f10.2,1x,9(''-------''))') zm(k)
      write(io6, '(10x,i5,f10.2,f11.2,f11.4,f12.2,f11.4)')  &
          k,zt(k),press(k,iw),rho(k,iw),theta(k,iw),sh_w(k,iw)*1.e3
   enddo

   write(io6, '(f10.2,1x,9(''-------''))') zm(1)
   write(io6,*)' '

endif

! Return if nudging flag is not activated

if (nudflag < 1 .or. nudrad < 1) return

! If we got here, nudging will be done in this model run, so fill nudging arrays.

! Swap future data time into past data time if necessary.  (Only polygon nudging
! has been implemented.)

if (iaction == 1) then

   do inudp = 2,mnudp
      do k = 2,mza-1
            rho_obsp(k,inudp) =    rho_obsf(k,inudp)
          theta_obsp(k,inudp) =  theta_obsf(k,inudp)
            shw_obsp(k,inudp) =    shw_obsf(k,inudp)
         uzonal_obsp(k,inudp) = uzonal_obsf(k,inudp)
         umerid_obsp(k,inudp) = umerid_obsf(k,inudp)
      enddo
   enddo

endif

! Zero out nudging polygon arrays and volume counter prior to summing

do inudp = 2,mnudp
   do k = 2,mza-1
          volnudp(k,inudp) = 0.

         rho_obsf(k,inudp) = 0.
       theta_obsf(k,inudp) = 0.
         shw_obsf(k,inudp) = 0.
      uzonal_obsf(k,inudp) = 0.
      umerid_obsf(k,inudp) = 0.
   enddo
enddo

! Sum data to nudging polygon arrays

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)
   inudp1 = itab_w(iw)%inudp(1)
!---------------------------------------------------------------------
call qsub('W',iw)

   do k = lpw(iw),mza-1
          volnudp(k,inudp1) =     volnudp(k,inudp1) + volt(k,iw) 

         rho_obsf(k,inudp1) =    rho_obsf(k,inudp1) + o_rho   (k,iw)*volt(k,iw)
       theta_obsf(k,inudp1) =  theta_obsf(k,inudp1) + o_theta (k,iw)*volt(k,iw)
         shw_obsf(k,inudp1) =    shw_obsf(k,inudp1) + o_shv   (k,iw)*volt(k,iw)
      uzonal_obsf(k,inudp1) = uzonal_obsf(k,inudp1) + o_uzonal(k,iw)*volt(k,iw)
      umerid_obsf(k,inudp1) = umerid_obsf(k,inudp1) + o_umerid(k,iw)*volt(k,iw)

   enddo

enddo
call rsub('Wbb',7)

! Compute average for each polygon nudge point

do inudp = 2,mnudp
   do k = 2,mza-1

                  volnudpi = 1. / max(1.,volnudp(k,inudp))

         rho_obsf(k,inudp) =    rho_obsf(k,inudp) * volnudpi
       theta_obsf(k,inudp) =  theta_obsf(k,inudp) * volnudpi
         shw_obsf(k,inudp) =    shw_obsf(k,inudp) * volnudpi
      uzonal_obsf(k,inudp) = uzonal_obsf(k,inudp) * volnudpi
      umerid_obsf(k,inudp) = umerid_obsf(k,inudp) * volnudpi

   enddo
enddo

return
end subroutine fldsisan

!===========================================================================

subroutine obs_nudge(rhot)

use mem_nudge,   only: tnudcent, mnudp,                                    &
                       rho_sim, rho_obs, rho_obsp, rho_obsf,               &
                       theta_sim, theta_obs, theta_obsp, theta_obsf,       &
                       uzonal_sim, uzonal_obs, uzonal_obsp, uzonal_obsf,   &
                       umerid_sim, umerid_obs, umerid_obsp, umerid_obsf,   &
                       shw_sim, shw_obs, shw_obsp, shw_obsf
use mem_basic,   only: uc, rho, theta, sh_w
use mem_grid,    only: mza, mwa, xeu, yeu, zeu, xew, yew, zew,  &
                       unx, uny, unz, lpu, lpw, volt
use misc_coms,   only: io6, time8, s1900_sim
use mem_ijtabs,  only: istp, jtab_u, jtab_w, itab_u, itab_w, mrl_begl
use consts_coms, only: erad, eradi
use mem_tend,    only: umt, thilt, sh_wt
use isan_coms,   only: ifgfile, s1900_fg
use rastro_evts

implicit none

! Nudge selected model fields (rho, thil, sh_w, umc) to observed data
! using polygon filtering

real, intent(inout) :: rhot(mza,mwa)

integer :: inudp,k,j,iu1,iu2,iu3,iw,inudp1,inudp2,inudp3,iw1,iw2,iu,mrl

real :: volnudp (mza,mnudp) ! automatic array
real :: uzonalt (mza,mwa)   ! automatic array
real :: umeridt (mza,mwa)   ! automatic array

real :: volnudpi,tp,tf,tnudi
real :: vxu1,vyu1,vzu1
real :: vxu2,vyu2,vzu2
real :: vxu3,vyu3,vzu3
real :: vx,vy,vz,raxis,uzonal,umerid,ugt,vgt,uvgrt,uvgxt,uvgyt,uvgzt
real :: fnudp1,fnudp2,fnudp3

real :: umass,umassoraxis

!sp
real :: umtbef

#ifdef OLAM_RASTRO
character(len=*) :: rst_buf = '_'
call rst_event_s_f(OLAM_OBS_NUDGE_IN,rst_buf)
#endif


! Check whether it is time to nudge

mrl = mrl_begl(istp)
if (mrl < 1) then

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_OBS_NUDGE_OUT,rst_buf)
#endif
	return
endif

! Time interpolation coefficients

tf = real ( (s1900_sim         - s1900_fg(ifgfile-1)) &
          / (s1900_fg(ifgfile) - s1900_fg(ifgfile-1)) )

tp = 1. - tf
tnudi = 1. / tnudcent

! Zero out polygon arrays for model fields and volume counter prior to summing

do inudp = 2,mnudp
   do k = 2,mza-1
         rho_sim(k,inudp) = 0.
       theta_sim(k,inudp) = 0.
         shw_sim(k,inudp) = 0.
      uzonal_sim(k,inudp) = 0.
      umerid_sim(k,inudp) = 0.
      
      volnudp(k,inudp) = 0.
   enddo
enddo

! Sum model values to polygon arrays

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(23)%jend(mrl); iw = jtab_w(23)%iw(j)

    iu1 = itab_w(iw)%iu1;   iu2 = itab_w(iw)%iu2;   iu3 = itab_w(iw)%iu3      	

   vxu1 = itab_w(iw)%vxu1; vyu1 = itab_w(iw)%vyu1; vzu1 = itab_w(iw)%vzu1
   vxu2 = itab_w(iw)%vxu2; vyu2 = itab_w(iw)%vyu2; vzu2 = itab_w(iw)%vzu2
   vxu3 = itab_w(iw)%vxu3; vyu3 = itab_w(iw)%vyu3; vzu3 = itab_w(iw)%vzu3

   inudp1 = itab_w(iw)%inudp(1)
!---------------------------------------------------------------------
call qsub('W',iw)

   do k = lpw(iw),mza-1

! Evaluate zonal and meridional wind components from model

      vx = vxu1 * uc(k,iu1) + vxu2 * uc(k,iu2) + vxu3 * uc(k,iu3)
      vy = vyu1 * uc(k,iu1) + vyu2 * uc(k,iu2) + vyu3 * uc(k,iu3)
      vz = vzu1 * uc(k,iu1) + vzu2 * uc(k,iu2) + vzu3 * uc(k,iu3)

      raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis
      uzonal = (vy * xew(iw) - vx * yew(iw)) / raxis
      umerid = vz * raxis / erad  &
             - (vx * xew(iw) + vy * yew(iw)) * zew(iw) / (raxis * erad) 

! Sum model fields and volume to polygon arrays

         rho_sim(k,inudp1) =    rho_sim(k,inudp1) + rho  (k,iw) * volt(k,iw)
       theta_sim(k,inudp1) =  theta_sim(k,inudp1) + theta(k,iw) * volt(k,iw)
         shw_sim(k,inudp1) =    shw_sim(k,inudp1) + sh_w (k,iw) * volt(k,iw)
      uzonal_sim(k,inudp1) = uzonal_sim(k,inudp1) + uzonal      * volt(k,iw)
      umerid_sim(k,inudp1) = umerid_sim(k,inudp1) + umerid      * volt(k,iw)

         volnudp(k,inudp1) =    volnudp(k,inudp1) + volt(k,iw) 

   enddo

enddo
call rsub('Wa',23)

! Horizontal loop over nudging polygons

do inudp = 2,mnudp

! Vertical loop over nudging polygons

   do k = 2,mza-1

! Inverse volume of nudging cell

                 volnudpi = 1. / max(1.,volnudp(k,inudp))

! Compute model average for each polygon nudge point

         rho_sim(k,inudp) =    rho_sim(k,inudp) * volnudpi
       theta_sim(k,inudp) =  theta_sim(k,inudp) * volnudpi
         shw_sim(k,inudp) =    shw_sim(k,inudp) * volnudpi
      uzonal_sim(k,inudp) = uzonal_sim(k,inudp) * volnudpi
      umerid_sim(k,inudp) = umerid_sim(k,inudp) * volnudpi

! Interpolate observational fields in time

         rho_obs(k,inudp) = tp *    rho_obsp(k,inudp) + tf *    rho_obsf(k,inudp)
       theta_obs(k,inudp) = tp *  theta_obsp(k,inudp) + tf *  theta_obsf(k,inudp)
         shw_obs(k,inudp) = tp *    shw_obsp(k,inudp) + tf *    shw_obsf(k,inudp)
      uzonal_obs(k,inudp) = tp * uzonal_obsp(k,inudp) + tf * uzonal_obsf(k,inudp)
      umerid_obs(k,inudp) = tp * umerid_obsp(k,inudp) + tf * umerid_obsf(k,inudp)

   enddo
enddo

! Loop over all W columns, find 3 neighboring polygon points for each,
! and interpolate (obs - model) differences at each polygon point to the W point

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(23)%jend(mrl); iw = jtab_w(23)%iw(j)
   inudp1 = itab_w(iw)%inudp(1);  fnudp1 = itab_w(iw)%fnudp(1)
   inudp2 = itab_w(iw)%inudp(2);  fnudp2 = itab_w(iw)%fnudp(2)
   inudp3 = itab_w(iw)%inudp(3);  fnudp3 = itab_w(iw)%fnudp(3)
!---------------------------------------------------------------------
call qsub('W',iw)

   do k = lpw(iw),mza-1

      rhot(k,iw) = rhot(k,iw)  + tnudi * volt(k,iw) * (  &
                    fnudp1 * (rho_obs(k,inudp1) - rho_sim(k,inudp1))  &
                 +  fnudp2 * (rho_obs(k,inudp2) - rho_sim(k,inudp2))  &
                 +  fnudp3 * (rho_obs(k,inudp3) - rho_sim(k,inudp3))  )

     thilt(k,iw) = thilt(k,iw) + tnudi * volt(k,iw) * rho(k,iw) * (  &
                    fnudp1 * (theta_obs(k,inudp1) - theta_sim(k,inudp1))  &
                 +  fnudp2 * (theta_obs(k,inudp2) - theta_sim(k,inudp2))  &
                 +  fnudp3 * (theta_obs(k,inudp3) - theta_sim(k,inudp3))  )
                 
     sh_wt(k,iw) = sh_wt(k,iw) + tnudi * rho(k,iw) * (  &
                    fnudp1 * (shw_obs(k,inudp1) - shw_sim(k,inudp1))  &
                 +  fnudp2 * (shw_obs(k,inudp2) - shw_sim(k,inudp2))  &
                 +  fnudp3 * (shw_obs(k,inudp3) - shw_sim(k,inudp3))  )

   uzonalt(k,iw) = tnudi * (  &
                    fnudp1 * (uzonal_obs(k,inudp1) - uzonal_sim(k,inudp1))  &
                 +  fnudp2 * (uzonal_obs(k,inudp2) - uzonal_sim(k,inudp2))  &
                 +  fnudp3 * (uzonal_obs(k,inudp3) - uzonal_sim(k,inudp3))  )

   umeridt(k,iw) = tnudi * (  &
                    fnudp1 * (umerid_obs(k,inudp1) - umerid_sim(k,inudp1))  &
                 +  fnudp2 * (umerid_obs(k,inudp2) - umerid_sim(k,inudp2))  &
                 +  fnudp3 * (umerid_obs(k,inudp3) - umerid_sim(k,inudp3))  )

                
!if (itab_w(iu)%iwglobe == 804 .or. itab_w(iu)%iwglobe == 805) then
!   write(io6,*) 'nuda0 '
!   write(io6,*) 'nuda1 ',k,iw
!   write(io6,*) 'nuda2 ',inudp1,inudp2,inudp3
!   write(io6,*) 'nuda3 ',tnudi,fnudp1,fnudp2,fnudp3
!   write(io6,*) 'nuda4 ',umerid_obs(k,inudp1),umerid_obs(k,inudp2),  & 
!                          umerid_obs(k,inudp3) 
!   write(io6,*) 'nuda5 ',umerid_sim(k,inudp1),umerid_sim(k,inudp2),  &
!                          umerid_sim(k,inudp3) 
!endif

   enddo

enddo
call rsub('Wb',23)

! UMC, UC

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_u(13)%jend(mrl); iu = jtab_u(13)%iu(j)
   iw1 = itab_u(iu)%iw1; iw2 = itab_u(iu)%iw2
!----------------------------------------------------------------------
call qsub('U',iu)

   if (iw1 < 2) iw1 = iw2
   if (iw2 < 2) iw2 = iw1

! Average winds to U point and rotate at U point

   do k = lpu(iu),mza-1
   
      ugt = .5 * (uzonalt(k,iw1) + uzonalt(k,iw2))
      vgt = .5 * (umeridt(k,iw1) + umeridt(k,iw2))

      raxis = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)  ! dist from earth axis

      umass = rho(k,iw1) * volt(k,iw1) + rho(k,iw2) * volt(k,iw2)
      umassoraxis = umass / raxis
   
      uvgrt = -vgt * zeu(iu) * eradi  ! radially outward from axis

      uvgxt = (-ugt * yeu(iu) + uvgrt * xeu(iu)) * umassoraxis
      uvgyt = ( ugt * xeu(iu) + uvgrt * yeu(iu)) * umassoraxis
      uvgzt =   vgt * raxis * umass * eradi 

      umtbef = umt(k,iu)

      umt(k,iu) = umt(k,iu)  &
                + uvgxt * unx(iu) + uvgyt * uny(iu) + uvgzt * unz(iu)
                
!if (itab_u(iu)%iuglobe == 1206) then
!   write(io6,*) 'nud0 '
!   write(io6,*) 'nud1 ',k,iw1,iw2
!   write(io6,*) 'nud2 ',ugt,vgt
!   write(io6,*) 'nud3 ',uzonalt(k,iw1),uzonalt(k,iw2),  &
!                        umeridt(k,iw1),umeridt(k,iw2)
!   write(io6,*) 'nud4 ',umt(k,iu),umtbef
!endif



   enddo

enddo
call rsub('U',13)

#ifdef OLAM_RASTRO
call rst_event_s_f(OLAM_OBS_NUDGE_OUT,rst_buf)
#endif

return
end subroutine obs_nudge
