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
Module mem_micro

   real, allocatable :: sh_c   (:,:) ! cloud water spec hum [kg_cld/kg_air]
   real, allocatable :: sh_r   (:,:) ! rain spec dens [kg_rain/kg_air]
   real, allocatable :: sh_p   (:,:) ! pristine ice spec dens[kg_pris/kg_air]
   real, allocatable :: sh_s   (:,:) ! snow spec dens [kg_snow/kg_air]
   real, allocatable :: sh_a   (:,:) ! aggregates spec dens [kg_agg/kg_air]
   real, allocatable :: sh_g   (:,:) ! graupel spec dens [kg_graup/kg_air]
   real, allocatable :: sh_h   (:,:) ! hail spec dens [kg_hail/kg_air]
   real, allocatable :: con_c  (:,:) ! cloud drop number conc [#_cld/kg_air]
   real, allocatable :: con_r  (:,:) ! rain number conc [#_rain/kg_air]
   real, allocatable :: con_p  (:,:) ! pristine ice number conc [#_pris/kg_air]
   real, allocatable :: con_s  (:,:) ! snow number conc [#_snow/kg_air]
   real, allocatable :: con_a  (:,:) ! aggregates number conc [#_aggr/kg_air]
   real, allocatable :: con_g  (:,:) ! graupel number conc [#_graup/kg_air]
   real, allocatable :: con_h  (:,:) ! hail number conc [#_hail/kg_air]
   real, allocatable :: con_ccn(:,:) ! CCN number conc [#_ccn/kg_air]
   real, allocatable :: con_ifn(:,:) ! IFN number conc [#_ifn/kg_air]
   real, allocatable :: q2     (:,:) ! rain internal energy [J/kg]
   real, allocatable :: q6     (:,:) ! graupel internal energy [J/kg]
   real, allocatable :: q7     (:,:) ! hail internal energy [J/kg]

   real, allocatable :: accpr (:) ! sfc rain accum [kg/m^2]
   real, allocatable :: accpp (:) ! sfc pristine ice accum [kg/m^2]
   real, allocatable :: accps (:) ! sfc snow accum [kg/m^2]
   real, allocatable :: accpa (:) ! sfc aggregates  accum [kg/m^2]
   real, allocatable :: accpg (:) ! sfc graupel accum [kg/m^2]
   real, allocatable :: accph (:) ! sfc hail accum [kg/m^2]
   real, allocatable :: pcprr (:) ! sfc rain pcp rate [kg/(m^2 s)]
   real, allocatable :: pcprp (:) ! sfc pristine ice pcp rate [kg/(m^2 s)]
   real, allocatable :: pcprs (:) ! sfc snow pcp rate [kg/(m^2 s)]
   real, allocatable :: pcpra (:) ! sfc aggregates pcp rate [kg/(m^2 s)]
   real, allocatable :: pcprg (:) ! sfc graupel pcp rate [kg/(m^2 s)]
   real, allocatable :: pcprh (:) ! sfc hail pcp rate [kg/(m^2 s)]
   real, allocatable :: pcpgr (:) ! sfc total pcp rate [kg/(m^2 s)]
   real, allocatable :: qpcpgr(:) ! sfc total pcp energy flux [J/(m^2 s)]
   real, allocatable :: dpcpgr(:) ! sfc total pcp depth accum rate [m/s]

Contains

!===============================================================================

   subroutine alloc_micro(mza,mwa,level  &
      ,icloud,irain,ipris,isnow,iaggr,igraup,ihail,jnmb)

   implicit none

   integer, intent(in) :: mza,mwa,level
   integer, intent(in) :: icloud,irain,ipris,isnow,iaggr,igraup,ihail
   integer, intent(in) :: jnmb(7)

! Allocate arrays based on options (if necessary)
! Initialize arrays to zero

   if (level >= 2 ) then
      allocate (sh_c(mza,mwa)) ; sh_c(1:mza,1:mwa) = 0.
   endif

   if (level >= 3) then

      if (irain >= 1)  then
         allocate (sh_r(mza,mwa)) ; sh_r(1:mza,1:mwa) = 0.
         allocate (q2(mza,mwa))   ; q2  (1:mza,1:mwa) = 0.
         allocate (accpr(mwa))    ; accpr(1:mwa) = 0.
         allocate (pcprr(mwa))    ; pcprr(1:mwa) = 0.
      endif

      if (ipris >= 1)  then
         allocate (sh_p(mza,mwa)) ; sh_p(1:mza,1:mwa) = 0.
         allocate (accpp(mwa))    ; accpp(1:mwa) = 0.
         allocate (pcprp(mwa))    ; pcprp(1:mwa) = 0.
      endif

      if (isnow >= 1)  then
         allocate (sh_s(mza,mwa)) ; sh_s (1:mza,1:mwa) = 0.
         allocate (accps(mwa))    ; accps(1:mwa) = 0.
         allocate (pcprs(mwa))    ; pcprs(1:mwa) = 0.
      endif

      if (iaggr >= 1)  then
         allocate (sh_a(mza,mwa)) ; sh_a (1:mza,1:mwa) = 0.
         allocate (accpa(mwa))    ; accpa(1:mwa) = 0.
         allocate (pcpra(mwa))    ; pcpra(1:mwa) = 0.
      endif

      if (igraup >= 1) then
         allocate (sh_g(mza,mwa)) ; sh_g(1:mza,1:mwa) = 0.
         allocate (q6(mza,mwa))   ; q6  (1:mza,1:mwa) = 0.
         allocate (accpg(mwa))    ; accpg(1:mwa) = 0.
         allocate (pcprg(mwa))    ; pcprg(1:mwa) = 0.
      endif

      if (ihail >= 1)  then
         allocate (sh_h(mza,mwa)) ; sh_h (1:mza,1:mwa) = 0.
         allocate (q7(mza,mwa))   ; q7   (1:mza,1:mwa) = 0.
         allocate (accph(mwa))    ; accph(1:mwa) = 0.
         allocate (pcprh(mwa))    ; pcprh(1:mwa) = 0.
      endif

      if (jnmb(1) == 5) then
          allocate (con_c(mza,mwa)) ; con_c(1:mza,1:mwa) = 0.
      endif

      if (jnmb(2) == 5) then
         allocate (con_r(mza,mwa)) ; con_r(1:mza,1:mwa) = 0.
      endif

      if (jnmb(3) == 5) then
         allocate (con_p(mza,mwa)) ; con_p(1:mza,1:mwa) = 0.
      endif

      if (jnmb(4) == 5) then
         allocate (con_s(mza,mwa)) ; con_s(1:mza,1:mwa) = 0.
      endif

      if (jnmb(5) == 5) then
         allocate (con_a(mza,mwa)) ; con_a(1:mza,1:mwa) = 0.
      endif

      if (jnmb(6) == 5) then
         allocate (con_g(mza,mwa)) ; con_g(1:mza,1:mwa) = 0.
      endif

      if (jnmb(7) == 5) then
         allocate (con_h(mza,mwa)) ; con_h(1:mza,1:mwa) = 0.
      endif

      if (icloud == 7)  then
         allocate (con_ccn(mza,mwa)) ; con_ccn(1:mza,1:mwa) = 0.
      endif

      if (ipris == 7)  then
         allocate (con_ifn(mza,mwa)) ; con_ifn(1:mza,1:mwa) = 0.
      endif

      allocate (pcpgr(mwa))  ; pcpgr (1:mwa) = 0.
      allocate (qpcpgr(mwa)) ; qpcpgr(1:mwa) = 0.
      allocate (dpcpgr(mwa)) ; dpcpgr(1:mwa) = 0.

   endif

   return
   end subroutine alloc_micro

!===============================================================================

   subroutine dealloc_micro()

   implicit none

   if (allocated(sh_c))    deallocate (sh_c)
   if (allocated(sh_r))    deallocate (sh_r)
   if (allocated(sh_p))    deallocate (sh_p)
   if (allocated(sh_s))    deallocate (sh_s)
   if (allocated(sh_a))    deallocate (sh_a)
   if (allocated(sh_g))    deallocate (sh_g)
   if (allocated(sh_h))    deallocate (sh_h)
   if (allocated(con_c))   deallocate (con_c)
   if (allocated(con_r))   deallocate (con_r)
   if (allocated(con_p))   deallocate (con_p)
   if (allocated(con_s))   deallocate (con_s)
   if (allocated(con_a))   deallocate (con_a)
   if (allocated(con_g))   deallocate (con_g)
   if (allocated(con_h))   deallocate (con_h)
   if (allocated(con_ccn)) deallocate (con_ccn)
   if (allocated(con_ifn)) deallocate (con_ifn)
   if (allocated(q2))      deallocate (q2)
   if (allocated(q6))      deallocate (q6)
   if (allocated(q7))      deallocate (q7)

   if (allocated(accpr))   deallocate (accpr)
   if (allocated(accpp))   deallocate (accpp)
   if (allocated(accps))   deallocate (accps)
   if (allocated(accpa))   deallocate (accpa)
   if (allocated(accpg))   deallocate (accpg)
   if (allocated(accph))   deallocate (accph)
   if (allocated(pcprr))   deallocate (pcprr)
   if (allocated(pcprp))   deallocate (pcprp)
   if (allocated(pcprs))   deallocate (pcprs)
   if (allocated(pcpra))   deallocate (pcpra)
   if (allocated(pcprg))   deallocate (pcprg)
   if (allocated(pcprh))   deallocate (pcprh)
   if (allocated(pcpgr))   deallocate (pcpgr)
   if (allocated(qpcpgr))  deallocate (qpcpgr)
   if (allocated(dpcpgr))  deallocate (dpcpgr)

   return
   end subroutine dealloc_micro

!===============================================================================

   subroutine filltab_micro(mza,mwa)

   use var_tables, only: vtables

   implicit none

   integer, intent(in) :: mza,mwa

   integer :: ndims
   integer, dimension(2) :: idims

   ndims = 2
   idims(1) = mza
   idims(2) = mwa

   if (allocated(sh_c))    call vtables(ndims,idims,'SH_C    :hist:mpt1'  &
         ,rvara2=sh_c)
   if (allocated(sh_r))    call vtables(ndims,idims,'SH_R    :hist:mpt1'  &
         ,rvara2=sh_r)
   if (allocated(sh_p))    call vtables(ndims,idims,'SH_P    :hist:mpt1'  &
         ,rvara2=sh_p)
   if (allocated(sh_s))    call vtables(ndims,idims,'SH_S    :hist:mpt1'  &
         ,rvara2=sh_s)
   if (allocated(sh_a))    call vtables(ndims,idims,'SH_A    :hist:mpt1'  &
         ,rvara2=sh_a)
   if (allocated(sh_g))    call vtables(ndims,idims,'SH_G    :hist:mpt1'  &
         ,rvara2=sh_g)
   if (allocated(sh_h))    call vtables(ndims,idims,'SH_H    :hist:mpt1'  &
         ,rvara2=sh_h)
   if (allocated(con_c))   call vtables(ndims,idims,'CON_C   :hist:mpt1'  &
         ,rvara2=con_c)
   if (allocated(con_r))   call vtables(ndims,idims,'CON_R   :hist:mpt1'  &
         ,rvara2=con_r)
   if (allocated(con_p))   call vtables(ndims,idims,'CON_P   :hist:mpt1'  &
         ,rvara2=con_p)
   if (allocated(con_s))   call vtables(ndims,idims,'CON_S   :hist:mpt1'  &
         ,rvara2=con_s)
   if (allocated(con_a))   call vtables(ndims,idims,'CON_A   :hist:mpt1'  &
         ,rvara2=con_a)
   if (allocated(con_g))   call vtables(ndims,idims,'CON_G   :hist:mpt1'  &
         ,rvara2=con_g)
   if (allocated(con_h))   call vtables(ndims,idims,'CON_H   :hist:mpt1'  &
         ,rvara2=con_h)
   if (allocated(con_ccn)) call vtables(ndims,idims,'CON_CCN :hist:mpt1'  &
         ,rvara2=con_ccn)
   if (allocated(con_ifn)) call vtables(ndims,idims,'CON_IFN :hist:mpt1'  &
         ,rvara2=con_ifn)
   if (allocated(q2))      call vtables(ndims,idims,'Q2      :hist:mpt1'  &
         ,rvara2=q2)
   if (allocated(q6))      call vtables(ndims,idims,'Q6      :hist:mpt1'  &
         ,rvara2=q6)
   if (allocated(q7))      call vtables(ndims,idims,'Q7      :hist:mpt1'  &
         ,rvara2=q7)

   ndims = 1
   idims(1) = mwa

   if (allocated(accpr))  call vtables(ndims,idims,'ACCPR  :hist',rvara1=accpr)
   if (allocated(accpp))  call vtables(ndims,idims,'ACCPP  :hist',rvara1=accpp)
   if (allocated(accps))  call vtables(ndims,idims,'ACCPS  :hist',rvara1=accps)
   if (allocated(accpa))  call vtables(ndims,idims,'ACCPA  :hist',rvara1=accpa)
   if (allocated(accpg))  call vtables(ndims,idims,'ACCPG  :hist',rvara1=accpg)
   if (allocated(accph))  call vtables(ndims,idims,'ACCPH  :hist',rvara1=accph)
   if (allocated(pcprr))  call vtables(ndims,idims,'PCPRR  :hist',rvara1=pcprr)
   if (allocated(pcprp))  call vtables(ndims,idims,'PCPRP  :hist',rvara1=pcprp)
   if (allocated(pcprs))  call vtables(ndims,idims,'PCPRS  :hist',rvara1=pcprs)
   if (allocated(pcpra))  call vtables(ndims,idims,'PCPRA  :hist',rvara1=pcpra)
   if (allocated(pcprg))  call vtables(ndims,idims,'PCPRG  :hist',rvara1=pcprg)
   if (allocated(pcprh))  call vtables(ndims,idims,'PCPRH  :hist',rvara1=pcprh)
   if (allocated(pcpgr))  call vtables(ndims,idims,'PCPGR  :hist',rvara1=pcpgr)
   if (allocated(qpcpgr)) call vtables(ndims,idims,'QPCPGR :hist',rvara1=qpcpgr)
   if (allocated(dpcpgr)) call vtables(ndims,idims,'DPCPGR :hist',rvara1=dpcpgr)

   return
   end subroutine filltab_micro

End Module mem_micro
