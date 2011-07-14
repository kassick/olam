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
Module mem_ijtabs

   use max_dims, only: maxgrds, maxremote

   integer, parameter :: nloops_m =  4           ! max # of DO loops for M pts
   integer, parameter :: nloops_u = 25 + maxremote ! max # of DO loops for U pts
   integer, parameter :: nloops_w = 35 + maxremote ! max # of DO loops for W pts
   integer, parameter :: maxtpn   =  8 ! max # (+1) of U/W pts around M pt 

   integer :: nstp  ! # of finest grid acoustic timesteps in coarse grid dtlong
   integer :: istp  ! Current timestep counter from 1 to nstp
   integer :: mrls  ! Number of active mesh refinement levels (MRLs)

   integer, allocatable :: mrl_begl(:)  ! MRL at beginning of long timestep
   integer, allocatable :: mrl_begs(:)  ! MRL at beginning of short timestep
   integer, allocatable :: mrl_endl(:)  ! MRL at end of long timestep
   integer, allocatable :: mrl_ends(:)  ! MRL at end of short timestep

   integer, allocatable :: leafstep(:)  ! flag to run leaf on any sub-timestep

   Type itab_m_vars             ! data structure for M pts (individual rank)
      logical :: loop(nloops_m) = .false. ! flag to perform each DO loop at this M pt

      integer :: ntpn  = 0      ! number of U/W neighbors of this M pt
      integer :: itopm = 1      ! M point from which to copy this M pt's topo
      integer :: iw(maxtpn) = 1 ! array of W neighbors of this M pt
      integer :: iu(maxtpn) = 1 ! array of U neighbors of this M pt
      integer :: imglobe = 1    ! global index of this M pt (in parallel case)

      real    :: arm = 0.0      ! polygon area bounded by W pts around this M pt
   End Type itab_m_vars

   Type itab_u_vars             ! data structure for U pts (individual rank)
      logical :: loop(nloops_u) = .false. ! flag to perform each DO loop at this M pt

      integer :: iup=1          ! U pt from which to copy this U pt's values
      integer :: im1=1, im2=1   ! neighbor M pts of this U pt
      integer :: iu1=1, iu2=1, iu3=1,  iu4=1,  iu5=1,  iu6=1 ! neighbor U pts
      integer :: iu7=1, iu8=1, iu9=1, iu10=1, iu11=1, iu12=1 ! neighbor U pts
      integer :: iw1=1, iw2=1, iw3=1,  iw4=1,  iw5=1,  iw6=1 ! neighbor W pts
      integer :: irank = -1     ! rank of parallel process at this U pt
      integer :: iuglobe = 1    ! global index of this U pt (in parallel case)
      integer :: mrlu=0         ! mesh refinement level of this U pt

      real :: diru1=0., diru2=0., diru3=0., diru4=0. ! pos direction of U neighbors
      real :: fuu5=0., fuu6=0., fuu7=0., fuu8=0.     ! proj coefs of U neighbors
      real :: fuu9=0., fuu10=0., fuu11=0., fuu12=0.  ! proj coefs of U neighbors
      real :: fuw3=0., fuw4=0., fuw5=0., fuw6=0.     ! proj coefs of W neighbors
      real :: tuu1=0., tuu2=0., tuu3=0., tuu4=0.     ! proj coefs of U neighbors
      real :: pgc12=0., pgc45=0., pgc63=0., pgc12b=0.    ! PGF proj coefs
      real :: pgc45b=0., pgc12c=0., pgc63c=0., pgc12d=0. ! PGF proj coefs
      real :: vxu1_u=0., vxu2_u=0., vxu3_u=0., vxu4_u=0., vxw1_u=0., vxw2_u=0. ! proj coefs
      real :: vyu1_u=0., vyu2_u=0., vyu3_u=0., vyu4_u=0., vyw1_u=0., vyw2_u=0. ! proj coefs
   End Type itab_u_vars

   Type itab_w_vars                   ! data structure for W pts (individual rank)
      logical :: loop(nloops_w) = .false. ! flag to perform each DO loop at this W pt

      integer :: iwp=1                ! W pt from which to copy this W pt's values
      integer :: im1=1, im2=1, im3=1  ! neighbor M pts of this W pt
      integer :: iu1=1, iu2=1, iu3=1, iu4=1, iu5=1, iu6=1, iu7=1, iu8=1, iu9=1 ! neighbor U pts
      integer :: iw1=1, iw2=1, iw3=1    ! neighbor W pts
      integer :: irank = -1             ! rank of parallel process at this W pt
      integer :: iwglobe = 1            ! global index of this W pt (in parallel case)
      integer :: mrlw=0, mrlw_orig=0    ! mesh refinement level of this W pt
      integer :: mrow=0, mrowh=0        ! Full and half row number outside nest

      real :: diru1=0.,  diru2=0., diru3=0.      ! pos direction of U neighbors
      real :: fwu4=0., fwu5=0., fwu6=0., fwu7=0., fwu8=0., fwu9=0.  ! proj coefs of U neighbors
      real :: fww1=0., fww2=0., fww3=0.         ! proj coefs of W neighbors
      real :: vxu1=0., vxu2=0., vxu3=0., vxw=0. ! proj coefs of U/W neighbors
      real :: vyu1=0., vyu2=0., vyu3=0., vyw=0. ! proj coefs of U/W neighbors
      real :: vzu1=0., vzu2=0., vzu3=0., vzw=0. ! proj coefs of U/W neighbors
      real :: vxu1_w=0., vxu2_w=0., vxu3_w=0.   ! proj coefs of U/W neighbors
      real :: vyu1_w=0., vyu2_w=0., vyu3_w=0.   ! proj coefs of U/W neighbors

      integer :: inudp(3) = 1  ! local nudpoly pts
      real    :: fnudp(3) = 0. ! local nudpoly coeffs
   End Type itab_w_vars

   Type itabg_m_vars            ! data structure for M pts (global)
      integer :: im_myrank = -1 ! local (parallel subdomain) index of this M pt
      integer :: irank = -1     ! rank of parallel process at this M pt
   End Type itabg_m_vars

   Type itabg_u_vars            ! data structure for U pts (global)
      integer :: iu_myrank = -1 ! local (parallel subdomain) index of this U pt
      integer :: irank = -1     ! rank of parallel process at this U pt
   End Type itabg_u_vars

   Type itabg_w_vars            ! data structure for W pts (global)
      integer :: iw_myrank = -1 ! local (parallel subdomain) index of this W pt
      integer :: irank = -1     ! rank of parallel process at this W pt
   End Type itabg_w_vars

   Type nest_u_vars         ! temporary U-pt data structure for spawning nested grids
      integer :: im=0, iu=0 ! new M/U pts attached to this U pt 
   End Type nest_u_vars
   
   Type nest_w_vars         ! temporary W-pt data structure for spawning nested grids
      integer :: iw1=0, iw2=0, iw3=0, iu1=0, iu2=0, iu3=0  ! new U/W pts attached to this W pt
   End Type nest_w_vars

   type (itab_m_vars), allocatable :: itab_m(:)
   type (itab_u_vars), allocatable :: itab_u(:)
   type (itab_w_vars), allocatable :: itab_w(:)

   type (itab_m_vars), allocatable :: ltab_m(:)
   type (itab_u_vars), allocatable :: ltab_u(:)
   type (itab_w_vars), allocatable :: ltab_w(:)

   type (itabg_m_vars), allocatable, target :: itabg_m(:)
   type (itabg_u_vars), allocatable, target :: itabg_u(:)
   type (itabg_w_vars), allocatable, target :: itabg_w(:)

   type (nest_u_vars), allocatable :: nest_u(:)
   type (nest_w_vars), allocatable :: nest_w(:)
   
   Type jtab_m_vars
      integer, allocatable :: im(:)
      integer, allocatable :: jend(:)
   End Type jtab_m_vars

   Type jtab_u_vars
      integer, allocatable :: iu(:)
      integer, allocatable :: jend(:)
   End Type jtab_u_vars

   Type jtab_w_vars
      integer, allocatable :: iw(:)
      integer, allocatable :: jend(:)
   End Type jtab_w_vars

   type (jtab_m_vars) :: jtab_m(nloops_m)
   type (jtab_u_vars) :: jtab_u(nloops_u)
   type (jtab_w_vars) :: jtab_w(nloops_w)

Contains

!===============================================================================

   subroutine alloc_itabs(mma,mua,mwa)

   implicit none

   integer, intent(in) :: mma,mua,mwa

   allocate (itab_m(mma))
   allocate (itab_u(mua))
   allocate (itab_w(mwa))

   return
   end subroutine alloc_itabs

!===============================================================================

   subroutine filltab_itabs(mma,mua,mwa)

   use var_tables, only: vtables
   use misc_coms,  only: iparallel, runtype

   implicit none

   integer, intent(in) :: mma,mua,mwa

   integer :: ndims
   integer :: idims(2)
   character(20) :: action

!  THESE ONLY NEED TO BE WRITTEN TO HISTORY FILE FOR PARALLEL RUNS,
!  AND READ FOR PLOTONLY OR PARCOMBINE RUNS.

   if (iparallel==1 .or. runtype=='PLOTONLY' .or. runtype=='PARCOMBINE') then

      if (iparallel == 1) then
         action = ':hist:nohist'
      else
         action = ':hist'
      endif

      ndims = 1
      idims(2) = 1

      idims(1) = mma

!------------------------------------------------------------------------
      if (allocated(itab_m)) call vtables(ndims,idims  &
           ,'IMGLOBE '//trim(action) ,ivara1=itab_m(:)%imglobe)
!------------------------------------------------------------------------

      idims(1) = mua

!------------------------------------------------------------------------
      if (allocated(itab_u)) call vtables(ndims,idims  &
           ,'IUGLOBE '//trim(action) ,ivara1=itab_u(:)%iuglobe)
!------------------------------------------------------------------------
      if (allocated(itab_u)) call vtables(ndims,idims  &
           ,'IRANKU  '//trim(action)  ,ivara1=itab_u(:)%irank)
!------------------------------------------------------------------------

      idims(1) = mwa

!------------------------------------------------------------------------
      if (allocated(itab_w)) call vtables(ndims,idims  &
           ,'IWGLOBE '//trim(action) ,ivara1=itab_w(:)%iwglobe)
!------------------------------------------------------------------------
      if (allocated(itab_w)) call vtables(ndims,idims  &
           ,'IRANKW '//trim(action)  ,ivara1=itab_w(:)%irank)
!------------------------------------------------------------------------

   endif

   return
   end subroutine filltab_itabs

!===============================================================================

   subroutine fill_jtabs(mma,mua,mwa)

   use misc_coms,  only: io6, nqparm
   use rastro_evts

   implicit none

   integer, intent(in) :: mma,mua,mwa

   integer :: iw,iu,im,k,nl,mrl
   integer :: iloop,iw1,iw2,mrl0,jend

#ifdef OLAM_RASTRO
   call rst_event_iii_f(OLAM_FILL_JTABS_IN,mma,mua,mwa)
#endif

! Allocate and zero-fill jtab_ws_mpi%jend

   do iloop = 1,nloops_m
      allocate (jtab_m(iloop)%jend(mrls))
      jtab_m(iloop)%jend(1:mrls) = 0
   enddo
   
   do iloop = 1,nloops_u
      allocate (jtab_u(iloop)%jend(mrls))
      jtab_u(iloop)%jend(1:mrls) = 0
   enddo
   
   do iloop = 1,nloops_w
      allocate (jtab_w(iloop)%jend(mrls))
      jtab_w(iloop)%jend(1:mrls) = 0
   enddo

! Compute and store jtab%jend(1)

   do iloop = 1,nloops_m
      jtab_m(iloop)%jend(1) = 0
      do im = 2,mma
         if (itab_m(im)%loop(iloop)) then
            jtab_m(iloop)%jend(1) = jtab_m(iloop)%jend(1) + 1
         endif
      enddo
      jtab_m(iloop)%jend(1) = max(1,jtab_m(iloop)%jend(1))
   enddo

   do iloop = 1,nloops_u
      jtab_u(iloop)%jend(1) = 0
      do iu = 2,mua
         if (itab_u(iu)%loop(iloop)) then
            jtab_u(iloop)%jend(1) = jtab_u(iloop)%jend(1) + 1
         endif
      enddo
      jtab_u(iloop)%jend(1) = max(1,jtab_u(iloop)%jend(1))
   enddo

   do iloop = 1,nloops_w
      jtab_w(iloop)%jend(1) = 0
      do iw = 2,mwa
         if (itab_w(iw)%loop(iloop)) then
            jtab_w(iloop)%jend(1) = jtab_w(iloop)%jend(1) + 1
         endif
      enddo
      jtab_w(iloop)%jend(1) = max(1,jtab_w(iloop)%jend(1))
   enddo

! Allocate and zero-fill JTAB_M%IM, JTAB_U%IU, JTAB_W%IW

   do iloop = 1,nloops_m
      jend = jtab_m(iloop)%jend(1)
      allocate (jtab_m(iloop)%im(jend))
      jtab_m(iloop)%im(1:jend) = 0
   enddo

   do iloop = 1,nloops_u
      jend = jtab_u(iloop)%jend(1)
      allocate (jtab_u(iloop)%iu(jend))
      jtab_u(iloop)%iu(1:jend) = 0
   enddo

   do iloop = 1,nloops_w
      jend = jtab_w(iloop)%jend(1)
      allocate (jtab_w(iloop)%iw(jend))
      jtab_w(iloop)%iw(1:jend) = 0
   enddo

! Initialize JTAB%JEND counters to zero

   do iloop = 1,nloops_m
      jtab_m(iloop)%jend(1:mrls) = 0
   enddo

   do iloop = 1,nloops_u
      jtab_u(iloop)%jend(1:mrls) = 0
   enddo

   do iloop = 1,nloops_w
      jtab_w(iloop)%jend(1:mrls) = 0
   enddo

! ///////////////////////////////////////////////////////////////////////////
! Compute JTAB_M%IM
! ///////////////////////////////////////////////////////////////////////////

! Grid level 1 only

   do im = 2,mma
      do iloop = 1,4
         if (itab_m(im)%loop(iloop)) then
            jtab_m(iloop)%jend(1) = jtab_m(iloop)%jend(1) + 1
            jtab_m(iloop)%im(jtab_m(iloop)%jend(1)) = im
         endif
      enddo
   enddo
   
! ///////////////////////////////////////////////////////////////////////////
! Compute JTAB_U%IU
! ///////////////////////////////////////////////////////////////////////////

! MRL-independent loops

   do iu = 2,mua
      do iloop = 1,10
         if (itab_u(iu)%loop(iloop)) then
            jtab_u(iloop)%jend(1) = jtab_u(iloop)%jend(1) + 1
            jtab_u(iloop)%iu(jtab_u(iloop)%jend(1)) = iu
         endif
      enddo      
   enddo

! MRL-dependent loops

   do mrl = mrls,1,-1
      do iu = 2,mua
         do iloop = 11,nloops_u
            if (itab_u(iu)%loop(iloop) .and. itab_u(iu)%mrlu == mrl) then
               jtab_u(iloop)%jend(1:mrl) = jtab_u(iloop)%jend(1:mrl) + 1
               jtab_u(iloop)%iu(jtab_u(iloop)%jend(1)) = iu
            endif
         enddo
      enddo
   enddo

! ///////////////////////////////////////////////////////////////////////////
! Compute JTAB_W%IW
! ///////////////////////////////////////////////////////////////////////////

! MRL-independent loops

   do iw = 2,mwa
      do iloop = 1,10
         if (itab_w(iw)%loop(iloop)) then
            jtab_w(iloop)%jend(1) = jtab_w(iloop)%jend(1) + 1
            jtab_w(iloop)%iw(jtab_w(iloop)%jend(1)) = iw
         endif
      enddo      
   enddo

! MRL-dependent loops

   do mrl = mrls,1,-1
      do iw = 2,mwa
         do iloop = 11,nloops_w

            if (itab_w(iw)%loop(iloop) .and. itab_w(iw)%mrlw == mrl) then
               jtab_w(iloop)%jend(1:mrl) = jtab_w(iloop)%jend(1:mrl) + 1
               jtab_w(iloop)%iw(jtab_w(iloop)%jend(1)) = iw
            endif

         enddo
      enddo
   enddo

! NOTE:  For cumulus parameterization (W loop 15), MRL dependence of the
! parameterization computations is selected in the parameterization itself.
! Here, loop 15 MRL dependence is as for other loops in order to copy
! THSRC and RTSRC to tendency arrays at the proper time for each MRL.

#ifdef OLAM_RASTRO
   call rst_event_iii_f(OLAM_FILL_JTABS_OUT,mma,mua,mwa)
#endif
   return
   end subroutine fill_jtabs

End Module mem_ijtabs

