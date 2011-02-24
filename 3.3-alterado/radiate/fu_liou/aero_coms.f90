Module aero_coms

  integer, allocatable, dimension(:,:,:) :: kdeliq
  real, dimension(9,6,72,46) :: a6jday
  real, dimension(9,8,72,46) :: ddjday
  real, dimension(6,10) :: srbqex
  real, dimension(6,8) :: qxdust
  real, dimension(190,2,4) :: rhinfo
  real, dimension(6,190,4) :: srhqex
  real, dimension(6,120) :: sruqex
  real, dimension(6,120) :: sruqsc
  real, dimension(6,120) :: sruqcb
  real, dimension(33,120) :: truqex
  real, dimension(33,120) :: truqsc
  real, dimension(33,120) :: truqcb
  real, dimension(120) :: refu22
  real, dimension(120) :: q55u22
  real, dimension(6,25) :: srdqex
  real, dimension(25) :: refd25
  real, dimension(6,25) :: srsqex
  real, dimension(25) :: refs25
  real, dimension(72,46,9,8,12) :: tdust

  integer, parameter :: nac = 8  ! number of aerosol constituents 
  integer, parameter :: n_atau = 1 ! number of wavelengths at which aerosol optical depth is specified
  integer, allocatable, dimension(:) :: itps
  real, allocatable, dimension(:) :: a_wli
  real, allocatable, dimension(:,:) :: a_taus

Contains

  subroutine alloc_aerosols(mza,mwa,ilwrtyp,iswrtyp)

    implicit none

    integer, intent(in) :: mza
    integer, intent(in) :: mwa
    integer, intent(in) :: ilwrtyp
    integer, intent(in) :: iswrtyp

    if(ilwrtyp == 4 .or. iswrtyp == 4)then

       ! Check to see if we have already allocated memory
       if(.not.allocated(kdeliq))then

          allocate(kdeliq(mza,mwa,4))

          if(nac > 0)then
             allocate(itps(nac))
             allocate(a_wli(n_atau))
             allocate(a_taus(n_atau,nac))
             !  4: 0.5 um mineral dust
             !  5: 1 um mineral dust
             !  6: 2 um mineral dust
             !  7: 4 um mineral dust
             !  8: 8 um mineral dust
             ! 10: water-soluble; sulfates, ammoniums, organic carbon
             ! 11: soot; black carbon (insoluble)
             ! 12: sea salt (accumulation mode)
             ! 13: sea salt (coarse mode)
             ! 14: mineral dust, nucleation mode
             ! 15: mineral dust, accumulation mode
             ! 16: mineral dust, coarse mode
             ! 17: mineral dust, transported mode
             if(nac == 8)then
                itps(1) = 10
                itps(2) = 11 
                itps(3) = 12
                itps(4) = 4 
                itps(5) = 5
                itps(6) = 6 
                itps(7) = 7
                itps(8) = 8 
             else
                print*,'bad nac in aero_coms.'
             endif
             if(n_atau == 1)then
                !  a_wli(1) = 3.1 ! (um) 
                !  a_wli(2) = 1.85
                !  a_wli(3) = 1.375
                !  a_wli(4) = 1.055
                !  a_wli(5) = 0.815
                a_wli(1) = 0.535
             else
                print*,'bad n_atau in aero_coms.'
             endif
          endif

       endif
       kdeliq(1:mza,1:mwa,1:4) = 1

    endif

    return
  end subroutine alloc_aerosols

end Module aero_coms
