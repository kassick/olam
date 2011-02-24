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
Module var_tables

!-------------------------------------------------------------------

type var_tables_r
   
   integer,      pointer :: ivar1_p(:)
   integer,      pointer :: ivar2_p(:,:)
   integer,      pointer :: ivar3_p(:,:,:)
   real,         pointer :: rvar1_p(:)
   real,         pointer :: rvar2_p(:,:)
   real,         pointer :: rvar3_p(:,:,:)
   real(kind=8), pointer :: dvar1_p(:)
   real(kind=8), pointer :: dvar2_p(:,:)
   real(kind=8), pointer :: dvar3_p(:,:,:)
   
   character (len=32) :: name

   integer :: ndims
   integer :: idims(3)
   integer :: ihist
   integer :: nohist
   integer :: impt1
   integer :: impt2
   integer :: impt3
   
end type

type(var_tables_r), allocatable :: vtab_r(:)  ! Allocated to maxvars

integer :: num_var

!-------------------------------------------------------------------

type scalar_table
   
   real, pointer :: var_p(:,:)
   real, pointer :: var_t(:,:)

   character (len=32) :: name
   
end type

type(scalar_table), allocatable :: scalar_tab(:)

integer :: num_scalar

!-------------------------------------------------------------------

type ED_table

   integer,      pointer :: ivar1_p(:)
   integer,      pointer :: ivar2_p(:,:)
   real,         pointer :: rvar1_p(:)
   real,         pointer :: rvar2_p(:,:)
   real(kind=8), pointer :: dvar2_p(:,:)
   
   integer :: mavg ! write-out on the monthly average?
   integer :: yavg ! write-out on the yearly average?
   integer :: ndims
   integer :: idims(3) 
   integer :: ihist
   integer :: impt1
   integer :: impt2
   integer :: impt3

   character (len=32) :: name

end type

type(ED_table), allocatable :: vtab_ED(:)

integer :: num_ED

!-------------------------------------------------------------------

type var_table_par

   real, pointer :: rvar2_p(:,:)

end type

type(var_table_par), allocatable :: vtab_par(:)

integer :: nvar_par = 0

Contains

!===============================================================================

   subroutine vtables(ndims, idims, tabstr, ivara1, ivara2, ivara3,  &
                                            rvara1, rvara2, rvara3,  &
                                            dvara1, dvara2, dvara3)

   use max_dims,  only: maxvars
   use misc_coms, only: io6

   implicit none

   integer, intent(in) :: ndims
   integer, intent(in) :: idims(ndims)
   character (len=*), intent(in) :: tabstr

   integer,      optional, target :: ivara1(:)
   integer,      optional, target :: ivara2(:,:)
   integer,      optional, target :: ivara3(:,:,:)
   real,         optional, target :: rvara1(:)
   real,         optional, target :: rvara2(:,:)
   real,         optional, target :: rvara3(:,:,:)
   real(kind=8), optional, target :: dvara1(:)
   real(kind=8), optional, target :: dvara2(:,:)
   real(kind=8), optional, target :: dvara3(:,:,:)

   character (len=80) :: line
   character (len=1)  :: toksep=':', cdimen,ctype
   character (len=32) :: tokens(10)
   character (len=8)  :: cname,ctab
   
   integer :: ntok,nt,nv
   
   call tokenize1(tabstr,tokens,ntok,toksep)
   
   num_var = num_var + 1
   if (num_var > maxvars) then
      write(io6,*) 'There are more variables than the variable tables are'
      write(io6,*) 'allocated to handle. Please increase maxvars in max_dims.f90'
      write(io6,*) 'Maxvars = ', maxvars
      stop    'vtables: too many variables'
   endif
   nv = num_var

   vtab_r(nv)%name = tokens(1)
   vtab_r(nv)%ndims = ndims
   vtab_r(nv)%idims(1:ndims) = idims(1:ndims)

   vtab_r(nv)%ihist = 0
   vtab_r(nv)%nohist = 0
   vtab_r(nv)%impt1 = 0
   vtab_r(nv)%impt2 = 0
   vtab_r(nv)%impt3 = 0

   do nt = 2,ntok

      ctab = tokens(nt)         

      if     (ctab == 'hist') then
         vtab_r(nv)%ihist = 1
      elseif (ctab == 'nohist') then
         vtab_r(nv)%nohist = 1
      elseif (ctab == 'mpt1') then
         vtab_r(nv)%impt1 = 1
      elseif (ctab == 'mpt2') then
         vtab_r(nv)%impt2 = 1
      elseif (ctab == 'mpt3') then
         vtab_r(nv)%impt3 = 1
      else
         write(io6,*) 'Illegal table specification for var:', tokens(1),ctab
         stop 'bad var table'
      endif

   enddo
  
   nullify ( vtab_r(nv)%ivar1_p )
   nullify ( vtab_r(nv)%ivar2_p )
   nullify ( vtab_r(nv)%ivar3_p )
   nullify ( vtab_r(nv)%rvar1_p )
   nullify ( vtab_r(nv)%rvar2_p )
   nullify ( vtab_r(nv)%rvar3_p )
   nullify ( vtab_r(nv)%dvar1_p )
   nullify ( vtab_r(nv)%dvar2_p )
   nullify ( vtab_r(nv)%dvar3_p )

! Assign pointer according to type of data input

! integers
   if     (present(ivara1)) then
        vtab_r(nv)%ivar1_p  => ivara1
   elseif (present(ivara2)) then
        vtab_r(nv)%ivar2_p  => ivara2
   elseif (present(ivara3)) then
        vtab_r(nv)%ivar3_p  => ivara3

! reals

   elseif (present(rvara1)) then
        vtab_r(nv)%rvar1_p  => rvara1
   elseif (present(rvara2)) then
        vtab_r(nv)%rvar2_p  => rvara2

!--------------------------------------------
! Parallel communication table for scalars

      if (vtab_r(nv)%impt1 == 1) then
         nvar_par = nvar_par + 1
         vtab_par(nvar_par)%rvar2_p => rvara2
      endif
!--------------------------------------------

   elseif (present(rvara3)) then
        vtab_r(nv)%rvar3_p  => rvara3

! real(8)

   elseif (present(dvara1)) then
        vtab_r(nv)%dvar1_p  => dvara1
   elseif (present(dvara2)) then
        vtab_r(nv)%dvar2_p  => dvara2
   elseif (present(dvara3)) then
        vtab_r(nv)%dvar3_p  => dvara3

   else
      write(io6,*) 'Incorrect or missing data field argument in vtables'
      stop 'vtables: bad data field'
   endif

   return
   end subroutine vtables

!===============================================================================

   subroutine vtables_ED(ndims,idims,tabstr,ivara1,rvara1,ivara2,rvara2,dvara2)

   use misc_coms, only: io6

   implicit none

   integer, intent(in) :: ndims
   integer, intent(in) :: idims(ndims)
   character (len=*), intent(in) :: tabstr

   integer,      optional, target :: ivara1(:)
   integer,      optional, target :: ivara2(:,:)
   real,         optional, target :: rvara1(:)
   real,         optional, target :: rvara2(:,:)
   real(kind=8), optional, target :: dvara2(:,:)

   integer, target :: ivar
   real, target :: rvar

   character (len=80) :: line
   character (len=1)  :: toksep=':', cdimen,ctype
   character (len=32) :: tokens(10)
   character (len=8)  :: cname,ctab
   
   integer :: ntok,nt,nv
   
   call tokenize1(tabstr,tokens,ntok,toksep)
   
   num_ED = num_ED + 1
   nv = num_ED

   vtab_ED(nv)%name = tokens(1)
   vtab_ED(nv)%ndims = ndims
   vtab_ED(nv)%idims(1:ndims) = idims(1:ndims)

   vtab_ED(nv)%ihist = 0
   vtab_ED(nv)%impt1 = 0
   vtab_ED(nv)%impt2 = 0
   vtab_ED(nv)%impt3 = 0
   vtab_ED(nv)%mavg  = 0
   vtab_ED(nv)%yavg  = 0

  do nt = 2,ntok

      ctab = tokens(nt)         

      if (ctab == 'hist') then
         vtab_ED(nv)%ihist = 1
      elseif (ctab == 'mpt1') then
         vtab_ED(nv)%impt1 = 1
      elseif (ctab == 'mpt2') then
         vtab_ED(nv)%impt2 = 1
      elseif (ctab == 'mpt3') then
         vtab_ED(nv)%impt3 = 1
      elseif (ctab == 'mavg') then
         vtab_ED(nv)%mavg  = 1
      elseif (ctab == 'yavg') then
         vtab_ED(nv)%yavg  = 1
      else
         write(io6,*) 'Illegal table specification for var:', tokens(1),ctab
         stop 'bad ED var table'
      endif

   enddo
  
   nullify ( vtab_ED(nv)%ivar1_p )
   nullify ( vtab_ED(nv)%rvar1_p )
   nullify ( vtab_ED(nv)%ivar2_p )
   nullify ( vtab_ED(nv)%rvar2_p )
   nullify ( vtab_ED(nv)%dvar2_p )

! Assign pointer according to type of data input

   if     (present(ivara1)) then
        vtab_ED(nv)%ivar1_p => ivara1
   elseif (present(rvara1)) then
        vtab_ED(nv)%rvar1_p => rvara1
   elseif (present(ivara2)) then
        vtab_ED(nv)%ivar2_p => ivara2
   elseif (present(rvara2)) then
        vtab_ED(nv)%rvar2_p => rvara2
   elseif (present(dvara2)) then
        vtab_ED(nv)%dvar2_p => dvara2
   else
      write(io6,*) 'Incorrect or missing data field argument in vtables_ED'
      stop 'vtables_ED: bad data field'
   endif

   return
   end subroutine vtables_ED

!===============================================================================
   
   subroutine vtables_scalar(varp,vart,tabstr)

   use max_dims,  only: maxsclr
   use misc_coms, only: io6

   implicit none

   real, target :: varp(:,:)
   real, target :: vart(:,:)
   
   character (len=*), intent(in) :: tabstr

   character (len=80) :: line
   character (len=1)  :: toksep=':'
   character (len=32) :: tokens(10)
   character (len=32) :: cname
   
   integer :: ntok,nv,ns
     
   call tokenize1(tabstr,tokens,ntok,toksep)

   cname = tokens(1)

!    Fill in existing table slot or make new scalar slot

   num_scalar = num_scalar + 1

   if (num_scalar > maxsclr) then
      write(io6,*) 'There are more scalars than the scalar tables are'
      write(io6,*) 'allocated to handle. Please increase maxsclr in max_dims.f90'
      write(io6,*) 'Maxsclr = ', maxsclr
      stop    'vtables_scalar: too many scalars'
   endif

   nv = num_scalar
   scalar_tab(nv)%name = cname

   scalar_tab(nv)%var_p => varp
   scalar_tab(nv)%var_t => vart
  
   return
   end subroutine vtables_scalar

End Module var_tables

