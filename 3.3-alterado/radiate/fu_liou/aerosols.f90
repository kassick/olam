subroutine updaer(jyeara, jjdaya)

  use aero_coms, only: a6jday, ddjday, tdust

  implicit none

  integer, intent(in) :: jyeara
  integer, intent(in) :: jjdaya

  character(len=256), parameter :: aero_dir = '/cluster/home/dmm31/aerosols/'
  character(len=256), dimension(4), parameter :: rdfile = (/  &
       'sep2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850', &
       'sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990', &
       'sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990', &
       'sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990'  &
       /)
  character(len=256) :: aero_file

  integer :: ia
  integer :: ndd
  integer :: idd
  integer :: m
  character(len=80) :: xtitle
  real, dimension(72,46,9) :: amon
  real, dimension(72,46,9,12,10) :: predd
  real, dimension(72,46,9,12,8) :: suidd
  real, dimension(72,46,9,12,8) :: ocidd
  real, dimension(72,46,9,12,8) :: bcidd
  integer :: l
  integer :: j
  integer :: i
  real, dimension(72,46,9,0:12,6) :: a6year
  !     TROP AEROSOL 1850 BACKGROUND, INDUSTRIAL & BIO-BURNING PARAMETERS
  real, dimension(13) :: aermix=(/   &
!      Pre-Industrial+Natural 1850 Level  Industrial Process  BioMBurn
!      ---------------------------------  ------------------  --------
!       1    2    3    4    5    6    7    8    9   10   11   12   13
!      SNP  SBP  SSP  ANP  ONP  OBP  BBP  SUI  ANI  OCI  BCI  OCB  BCB
       1.0, 1.0, .26, 1.0, 2.5, 2.5, 1.9, 1.0, 1.0, 2.5, 1.9, 2.5, 1.9/)
  real, dimension(8) ::  &
!                TROPOSPHERIC AEROSOL COMPOSITIONAL/TYPE PARAMETERS
!                  SO4    SEA    ANT    OCX    BCI    BCB    DST   VOL
       drym2g = (/4.667, 0.866, 4.448, 5.017, 9.000, 9.000, 1.000,1.000/)

  !=================================================================
  real :: wtani
  real :: wtocb
  real :: wtbcb
  real, external :: glopop
  real :: wt75
  integer :: iys
  integer :: iyc
  integer :: jys
  integer :: jyc
  real :: swti
  real :: swtj
  real :: cwti
  real :: cwtj
  real :: xmi
  integer :: mi
  real :: wtmi
  real :: wtmj
  integer :: n
  integer :: mj
  real, dimension(9,6) :: global_mean

  !       Aerosol data pressure levels
  !  1010.,934.,854.,720.,550.,390.,255.,150.,70.,10.

  do ia = 1, 4
     aero_file = trim(aero_dir)//trim(rdfile(ia))
     open(12,file=trim(aero_file),form='unformatted',status='old',  &
          convert='big_endian')
     ndd = 8
     if(ia == 1)ndd = 10
     do idd = 1, ndd
        do m = 1, 12
           read(12)xtitle,amon
           do l = 1, 9
              do j = 1, 46
                 do i = 1, 72
                    if(ia == 1)predd(i,j,l,m,idd) = amon(i,j,l)
                    if(ia == 2)suidd(i,j,l,m,idd) = amon(i,j,l)
                    if(ia == 3)ocidd(i,j,l,m,idd) = amon(i,j,l)
                    if(ia == 4)bcidd(i,j,l,m,idd) = amon(i,j,l)
                 enddo
              enddo
           enddo
        enddo
     enddo
     close(12)
  enddo

  do m = 1, 12
     do l = 1, 9
        do j = 1, 46
           do i = 1, 72
              a6year(i,j,l,m,1) = (aermix(1) * predd(i,j,l,m,1) +   &
                   aermix(2) * predd(i,j,l,m,2)) * 1000.0 * drym2g(1)
              a6year(i,j,l,m,2) = aermix(3) * predd(i,j,l,m,3) * 1000.0 *   &
                   drym2g(2)
              a6year(i,j,l,m,3) = aermix(4) * predd(i,j,l,m,4) * 1000.0 *   &
                   drym2g(3)
              a6year(i,j,l,m,4) = (aermix(5) * predd(i,j,l,m,5) +  &
                   aermix(6) * predd(i,j,l,m,6)) * 1000.0 * drym2g(4)
              a6year(i,j,l,m,5) = 0.0
              a6year(i,j,l,m,6) = aermix(7) * predd(i,j,l,m,7) * 1000.0 *   &
                   drym2g(6)
           enddo
        enddo
     enddo
  enddo

  do l = 1, 9
     do j = 1, 46
        do i = 1, 72
           a6year(i,j,l,0,1) = a6year(i,j,l,12,1)
           a6year(i,j,l,0,2) = a6year(i,j,l,12,2)
           a6year(i,j,l,0,3) = a6year(i,j,l,12,3)
           a6year(i,j,l,0,4) = a6year(i,j,l,12,4)
           a6year(i,j,l,0,5) = a6year(i,j,l,12,5)
           a6year(i,j,l,0,6) = a6year(i,j,l,12,6)
        enddo
     enddo
  enddo

  if(jyeara > 1850)then
     wtani = glopop(jyeara)
     wtocb = min(1.0, (jyeara-1850)/140.0)
     wtbcb = wtocb
     do m = 1, 12
        do l = 1, 9
           do j = 1, 46
              do i = 1, 72
                 a6year(i,j,l,m,3) = a6year(i,j,l,m,3) +   &
                      aermix(9) * wtani * predd(i,j,l,m,8) * 1000.0 * drym2g(3)
                 a6year(i,j,l,m,4) = a6year(i,j,l,m,4) +   &
                      aermix(12) * wtocb * predd(i,j,l,m,9) * 1000.0 *   &
                      drym2g(4)
                 a6year(i,j,l,m,6) = a6year(i,j,l,m,6) +   &
                      aermix(13) * wtbcb * predd(i,j,l,m,10) * 1000.0 *   &
                      drym2g(6)
              enddo
           enddo
        enddo
     enddo
     wtani = glopop(jyeara-1)
     wtocb = min(139.0/140.0, (jyeara-1851)/140.0)
     wtbcb = wtocb
     m = 12
     do l = 1, 9 
        do j = 1, 46
           do i = 1, 72
              a6year(i,j,l,0,3) = a6year(i,j,l,0,3) +   &
                   aermix(9) * wtani * predd(i,j,l,m,8) * 1000.0 * drym2g(3)
              a6year(i,j,l,0,4) = a6year(i,j,l,0,4) +   &
                   aermix(12) * wtocb * predd(i,j,l,m,9) * 1000.0 * drym2g(4)
              a6year(i,j,l,0,6) = a6year(i,j,l,0,6) +   &
                   aermix(13) * wtbcb * predd(i,j,l,m,10) * 1000.0 * drym2g(6)
           enddo
        enddo
     enddo
  endif

  if(jyeara > 1850 .and. jyeara < 1876)then
     wt75 = (jyeara - 1850) / 25.0
     do m = 1, 12
        do l = 1, 9
           do j = 1, 46
              do i = 1, 72
                 a6year(i,j,l,m,1) = a6year(i,j,l,m,1) + wt75 *   &
                      suidd(i,j,l,m,1) * aermix(8) * 1000.0 * drym2g(1)
                 a6year(i,j,l,m,4) = a6year(i,j,l,m,4) + wt75 *   &
                      ocidd(i,j,l,m,1) * aermix(10) * 1000.0 * drym2g(4)
                 a6year(i,j,l,m,5) = a6year(i,j,l,m,5) + wt75 *   &
                      bcidd(i,j,l,m,1) * aermix(11) * 1000.0 * drym2g(5)
              enddo
           enddo
        enddo
     enddo
     wt75 = (jyeara - 1851) / 25.0
     m = 12
     do l = 1, 9
        do j = 1, 46
           do i = 1, 72
              a6year(i,j,l,0,1) = a6year(i,j,l,0,1) + wt75 *   &
                   suidd(i,j,l,m,1) * aermix(8) * 1000.0 * drym2g(1)
              a6year(i,j,l,0,4) = a6year(i,j,l,0,4) + wt75 *   &
                   ocidd(i,j,l,m,1) * aermix(10) * 1000.0 * drym2g(4)
              a6year(i,j,l,0,5) = a6year(i,j,l,0,5) + wt75 *   &
                   bcidd(i,j,l,m,1) * aermix(11) * 1000.0 * drym2g(5)
           enddo
        enddo
     enddo
  endif

  if(jyeara > 1875)then
     call strend(jyeara, iys, jys, swti, swtj)
     call ctrend(jyeara, iyc, jyc, cwti, cwtj)
     do m = 1, 12
        do l = 1, 9
           do j = 1, 46
              do i = 1, 72
                 a6year(i,j,l,m,1) = a6year(i,j,l,m,1) + aermix(8) *   &
                      1000.0 * drym2g(1) * (swti * suidd(i,j,l,m,iys) +  &
                      swtj * suidd(i,j,l,m,jys))
                 a6year(i,j,l,m,4) = a6year(i,j,l,m,4) + aermix(10) *   &
                      1000.0 * drym2g(4) * (cwti * ocidd(i,j,l,m,iyc) +  &
                      cwtj * ocidd(i,j,l,m,jyc))
                 a6year(i,j,l,m,5) = a6year(i,j,l,m,5) + aermix(11) *   &
                      1000.0 * drym2g(5) * (cwti * bcidd(i,j,l,m,iyc) +  &
                      cwtj * bcidd(i,j,l,m,jyc))
              enddo
           enddo
        enddo
     enddo
     call strend(jyeara-1, iys, jys, swti, swtj)
     call ctrend(jyeara-1, iyc, jyc, cwti, cwtj)
     m = 12
     do l = 1, 9
        do j = 1, 46
           do i = 1, 72
              a6year(i,j,l,0,1) = a6year(i,j,l,0,1) + aermix(8) *   &
                   1000.0 * drym2g(1) * (swti * suidd(i,j,l,m,iys) +  &
                   swtj * suidd(i,j,l,m,jys))
              a6year(i,j,l,0,4) = a6year(i,j,l,0,4) + aermix(10) *   &
                   1000.0 * drym2g(4) * (cwti * ocidd(i,j,l,m,iyc) +  &
                   cwtj * ocidd(i,j,l,m,jyc))
              a6year(i,j,l,0,5) = a6year(i,j,l,0,5) + aermix(11) *   &
                   1000.0 * drym2g(5) * (cwti * bcidd(i,j,l,m,iyc) +  &
                   cwtj * bcidd(i,j,l,m,jyc))
           enddo
        enddo
     enddo
  endif
     
  xmi = (jjdaya + jjdaya + 31 - (jjdaya + 15) / 61 + (jjdaya + 14) / 61) / 61.0
  mi = xmi
  wtmj = xmi - mi
  wtmi = 1.0 - wtmj
  if(mi > 11)mi = 0
  mj = mi + 1
  do j = 1, 46
     do i = 1, 72
        do n = 1, 6
           do l = 1, 9
              a6jday(l,n,i,j) = wtmi * a6year(i,j,l,mi,n) +   &
                   wtmj * a6year(i,j,l,mj,n)
           enddo
        enddo
     enddo
  enddo

!  global_mean(:,:) = 0.0
!  do j = 1, 46
!     do i = 1, 72
!        do n = 1, 6
!           do l = 1, 9
!              global_mean(l,n) = global_mean(l,n) + a6jday(l,n,i,j)
!           enddo
!        enddo
!     enddo
!  enddo
!  do n = 1, 6
!     do l = 1, 9
!        print*,n,l,global_mean(l,n)
!     enddo
!  enddo

!  call dust(jjdaya)
!  call volcano(jyeara,jjdaya)

  xmi = (jjdaya + jjdaya + 31 - (jjdaya + 15) / 61 + (jjdaya + 14) / 61) / 61.0
  mi = xmi
  wtmj = xmi - mi
  wtmi = 1.0 - wtmj
  if(mi < 1)mi=12
  if(mi > 12) mi = 1
  mj = mi + 1
  if(mj > 12) mj = 1
  do j = 1, 46
     do i = 1, 72
        do n = 1, 8
           do l = 1, 9
              ddjday(l,n,i,j) = wtmi * tdust(i,j,l,n,mi) + wtmj *   &
                   tdust(i,j,l,n,mj)
           enddo
        enddo
     enddo
  enddo

  return
end subroutine updaer

!==============================================================

function glopop(jyear)

  implicit none

  real :: glopop

!     ----------------------------------------------------------------
!     GLOPOP = normalized global population trend set to unity in 1990
!              based on UN statistics & population projections to 2050
!
!     GLOPOP = 0.000 for 1850 and earlier
!            = 1.000 for 1990
!            = 1.658 for 2050 and later
!     ----------------------------------------------------------------

  integer, intent(in) :: jyear
  real :: GPNORM = 5.27-1.26 ,DNGPOP(21), GPOP(21) = (/  &
!               1850                     1900                     1950
       1.26,1.33,1.41,1.49,1.57,1.65,1.75,1.86,2.07,2.30,2.52  &
!                                        2000                     2050
       ,3.02,3.70,4.44,5.27,6.06,6.79,7.50,8.11,8.58,8.91/)
  integer i
  integer :: iy
  real :: xy
  real :: dy

  do i = 1, 21
     dngpop(i)=(gpop(i)-gpop(1))/gpnorm
  enddo
  xy=(jyear-1840)/10.0
  IY=XY
  DY=XY-IY
  IF(IY < 1) THEN
     IY=1
     DY=0.0
  ENDIF
  IF(IY > 20) THEN
     IY=20
     DY=1.0
  ENDIF
  GLOPOP=DNGPOP(IY)+DY*(DNGPOP(IY+1)-DNGPOP(IY))
  RETURN
END FUNCTION GLOPOP

!========================================================================

subroutine ctrend(jyear,idec,jdec,cwti,cwtj)
  implicit none

!-------------------------------------------------------------------
!     black carbon interdecadal tau interpolation is based on linear
!     tau trend (between decadal global taumaps) with a superimposed
!     intra-decadal time dependence scaled to the black carbon total
!     emission rate.
!
!        input: jyear   (julian year)
!
!                 ctrend coefficients refer to sep2003_oci_koch maps
!                 ctrend coefficients refer to sep2003_bci_koch maps
!                 --------------------------------------------------
!
!                 map=  1850 1875 1900 1925 1950 1960 1970 1980 1990
!       output:  idec=   (0)   1    2    3    4    5    6    7    8
!                jdec=  idec + 1    (returned idec,jdec are (1 to 8)
!
!                cwti=   (multiplicative weight for bc datamap idec)
!                cwtj=   (multiplicative weight for bc datamap jdec)
!
!        note:  time dependence is linear before 1950. industrial bc
!               is assumed 0 in 1850 so cwti=0, and idec is set to 1
!-------------------------------------------------------------------

  integer, intent(in)  :: jyear
  integer, intent(out) :: idec,jdec
  real,  intent(out) :: cwti,cwtj

!    global annual emissions of bc u   emission (mt/yr)

  real, parameter, dimension(5,45) :: bce = reshape( (/  &
!      year    hard_coal    brown_coal      diesel        total
       50.0, 2.280581713, 0.4449132979, 0.1599090248, 2.885536671, &
       51.0, 2.443193913, 0.4855868816, 0.1884280443, 3.117194653, &
       52.0, 2.473641872, 0.5115299225, 0.2027695477, 3.187930107, &
       53.0, 2.481340885, 0.5448409319, 0.2149295360, 3.241089582, &
       54.0, 2.505670071, 0.5780177116, 0.2343477309, 3.317960978, &
       55.0, 2.698692560, 0.6238067150, 0.2733324766, 3.595800638, &
       56.0, 2.855226278, 0.6531309485, 0.3043369055, 3.812692404, &
       57.0, 2.975781679, 0.6821750998, 0.3207367063, 3.978575468, &
       58.0, 3.341105223, 0.7035279870, 0.3370627165, 4.381746292, &
       59.0, 3.638528824, 0.7075053453, 0.3695519567, 4.715488434, &
       60.0, 3.770926714, 0.7416650057, 0.3832504749, 4.896034241, &
       61.0, 3.392980337, 0.7805693150, 0.4217525721, 4.595387459, &
       62.0, 3.288835049, 0.8179932237, 0.4603823125, 4.567360401, &
       63.0, 3.359177589, 0.8604368567, 0.5090782642, 4.728550911, &
       64.0, 3.432664871, 0.8952696323, 0.5388473868, 4.866865158, &
       65.0, 3.529418945, 0.8819132447, 0.5785927773, 4.989773750, &
       66.0, 3.577459812, 0.8817394972, 0.6323299408, 5.091631413, &
       67.0, 3.418204546, 0.8635972142, 0.6592246890, 4.941041946, &
       68.0, 3.452457905, 0.8943673372, 0.7338049412, 5.080585003, &
       69.0, 3.626069546, 0.9298774004, 0.7889106274, 5.344810009, &
       70.0, 3.264039755, 0.9229136109, 0.8880128860, 5.074741840, &
       71.0, 3.437611580, 0.9374827743, 0.9531223178, 5.328329086, &
       72.0, 3.473345757, 0.7836616039, 1.0180075170, 5.274850368, &
       73.0, 3.495583296, 0.8056778908, 1.1174367670, 5.418928623, &
       74.0, 3.506143808, 0.8251076341, 1.0828053950, 5.413989067, &
       75.0, 3.906814098, 0.8527192473, 1.0454736950, 5.804963112, &
       76.0, 4.005736828, 0.8900613785, 1.1400985720, 6.035901546, &
       77.0, 4.236912251, 0.9103702307, 1.2190728190, 6.366260529, &
       78.0, 4.459666252, 0.9303293228, 1.2408012150, 6.630728722, &
       79.0, 4.697422504, 0.9856286645, 1.3019220830, 6.984815121, &
       80.0, 4.796229839, 0.9959300756, 1.2336660620, 7.026207924, &
       81.0, 4.789204121, 1.0459070210, 1.1664049630, 7.001126766, &
       82.0, 4.872739315, 1.0975246430, 1.1601715090, 7.130136490, &
       83.0, 4.983223438, 1.1424025300, 1.1732926370, 7.298912525, &
       84.0, 5.265352249, 1.2178678510, 1.2251536850, 7.708741188, &
       85.0, 5.763637543, 1.2965050940, 1.2428865430, 8.303324699, &
       86.0, 5.924767494, 1.3386499880, 1.2930148840, 8.556744576, &
       87.0, 6.155550480, 1.3738890890, 1.3162037130, 8.845513344, &
       88.0, 6.379704475, 1.3670797350, 1.3813229800, 9.127896309, &
       89.0, 6.594299316, 1.4169263840, 1.4029121400, 9.414231300, &
       90.0, 6.566919804, 1.4685817960, 1.4224120380, 9.458042145, &
       91.0, 6.661097050, 1.2067918780, 1.4163945910, 9.284657478, &
       92.0, 7.737902641, 1.3509917260, 1.4471185210, 10.53625107, &
       93.0, 7.393332005, 1.2448183300, 1.4543261530, 10.09271908, &
       94.0, 7.515841007, 1.2333894970, 1.4780857560, 10.22745800 &
       /), (/5,45/) )

  real xdec
  integer ibcdec,jbcdec,ijyear

  if(jyear.lt.1876) then
     cwtj=(jyear-1850)/25.0
     if(cwtj.lt.0.0) cwtj=0.0
     cwti=0.0
     idec=1
     jdec=1
     return
  endif

  if(jyear.lt.1950) then
     xdec=(jyear-1850)/25.0
     idec=xdec
     jdec=idec+1
     cwtj=xdec-idec
     cwti=1.0-cwtj
     return
  endif

  if(jyear.lt.1990) then
     idec=(jyear-1910)/10
     jdec=idec+1
     ibcdec=1+(idec-4)*10
     jbcdec=ibcdec+10
     ijyear=jyear-1949
     cwtj=(bce(5,ijyear)-bce(5,ibcdec))  &
          /(bce(5,jbcdec)-bce(5,ibcdec))
     cwti=1.0-cwtj
     return
  endif
  
  if(jyear.gt.1989) then
     idec=7
     jdec=8
     ijyear=jyear-1949
     if(ijyear.gt.45) ijyear=45
     cwtj=bce(5,ijyear)/bce(5,41)
     cwti=0.0
  endif

  return
end subroutine ctrend

subroutine strend(jyear,idec,jdec,swti,swtj)
  implicit none

!-------------------------------------------------------------------
!     anthropogenic sulfate inter-decadal tau interpolation is based
!     on a linear tau trend (between decadal global tau-maps) with a
!     superimposed intradecadal time dependence scaled in proportion
!     to the anthropogenic sulfate global emission rate.
!
!        input: jyear   (julian year)
!
!                 ctrend coefficients refer to sep2003_sui_koch maps
!                 --------------------------------------------------
!
!                 map=  1850 1875 1900 1925 1950 1960 1970 1980 1990
!       output:  idec=   (0)   1    2    3    4    5    6    7    8
!                jdec=  idec + 1    (returned idec,jdec are (1 to 8)
!
!                swti=  (multiplicative weight for sui datamap idec)
!                swtj=  (multiplicative weight for sui datamap jdec)
!
!        note:  time dependence linear before 1950.   industrial sui
!               is assumed 0 in 1850 so swti=0, and idec is set to 1
!-------------------------------------------------------------------

  integer, intent(in)  :: jyear
  integer, intent(out) :: idec,jdec
  real,  intent(out) :: swti,swtj

!     global emission of sulfate

!     emission (mt/yr)
!               year      anthropogenic_sulfate natural_sulfate
  real, parameter, dimension(3,41) :: sue = reshape( (/  &
       1950.0,     30.46669769,           14.4, &
       1951.0,     32.38347244,           14.4, &
       1952.0,     32.18632889,           14.4, &
       1953.0,     32.83379745,           14.4, &
       1954.0,     32.79270935,           14.4, &
       1955.0,     35.79611969,           14.4, &
       1956.0,     39.93603897,           14.4, &
       1957.0,     38.68806839,           14.4, &
       1958.0,     39.35904312,           14.4, &
       1959.0,     41.06065369,           14.4, &
       1960.0,     42.67050934,           14.4, &
       1961.0,     41.32410431,           14.4, &
       1962.0,     41.80470276,           14.4, &
       1963.0,     43.26312637,           14.4, &
       1964.0,     44.68368530,           14.4, &
       1965.0,     45.81701660,           14.4, &
       1966.0,     46.61584091,           14.4, &
       1967.0,     46.42276001,           14.4, &
       1968.0,     47.77438354,           14.4, &
       1969.0,     49.30817032,           14.4, &
       1970.0,     52.81050873,           14.4, &
       1971.0,     52.95043945,           14.4, &
       1972.0,     54.10167694,           14.4, &
       1973.0,     55.93037415,           14.4, &
       1974.0,     57.31056213,           14.4, &
       1975.0,     58.52788162,           14.4, &
       1976.0,     59.71361542,           14.4, &
       1977.0,     62.59599304,           14.4, &
       1978.0,     61.98198318,           14.4, &
       1979.0,     64.71042633,           14.4, &
       1980.0,     65.28986359,           14.4, &
       1981.0,     63.23768234,           14.4, &
       1982.0,     62.88000488,           14.4, &
       1983.0,     61.45023346,           14.4, &
       1984.0,     63.85008621,           14.4, &
       1985.0,     66.47412872,           14.4, &
       1986.0,     68.00902557,           14.4, &
       1987.0,     69.87956238,           14.4, &
       1988.0,     70.52937317,           14.4, &
       1989.0,     72.06355286,           14.4, &
       1990.0,     71.29174805,           14.4 &
       /), (/3,41/) )

  real :: xdec
  integer isudec,jsudec,ijyear

  if(jyear.lt.1876) then
     swtj=(jyear-1850)/25.0
     if(swtj.lt.0.0) swtj=0.0
     swti=0.0
     idec=1
     jdec=1
     return
  endif

  if(jyear.lt.1950) then
     xdec=(jyear-1850)/25.0
     idec=xdec
     jdec=idec+1
     swtj=xdec-idec
     swti=1.0-swtj
     return
  endif

  if(jyear.lt.1990) then
     idec=(jyear-1910)/10
     jdec=idec+1
     isudec=1+(idec-4)*10
     jsudec=isudec+10
     ijyear=jyear-1949
     swtj=(sue(2,ijyear)-sue(2,isudec))  &
          /(sue(2,jsudec)-sue(2,isudec))
     swti=1.0-swtj
     return
  endif
  
  if(jyear.gt.1989) then
     idec=7
     jdec=8
     ijyear=jyear-1949
     if(ijyear.gt.41) ijyear=41
     swtj=sue(2,ijyear)/sue(2,41)
     swti=0.0
  endif

  return
end subroutine strend

!================================================================
subroutine getaer(nrad, iw, rhl, tlm, pl, waso, soot, sea_salt)

  use aero_coms, only: kdeliq, a6jday, srbqex, rhinfo, srhqex
  use mem_grid, only: glatw, glonw

  implicit none

  integer, intent(in) :: nrad
  integer, intent(in) :: iw
  real, dimension(nrad), intent(in) :: pl ! pressure (mb)
  real, dimension(nrad), intent(in) :: tlm ! mean layer temperature [K]
  real, dimension(nrad), intent(in) :: rhl ! relative humidity [0--1]

  real, dimension(nrad,6), intent(out) :: waso
  real, dimension(nrad,6), intent(out) :: soot
  real, dimension(nrad,6), intent(out) :: sea_salt

  integer, dimension(nrad,8) :: nrhnan
  integer :: k
  real :: xrh
  integer :: nrh
  integer :: na
  real :: rhdna
  real, external :: rhdtna
  ! Crystallization RH
  real, dimension(4), parameter :: rhc = (/ 0.38, 0.47, 0.28, 0.38 /)
  integer :: ilat
  integer :: ilon
  !r PLBA09 Vert. Layering for tropospheric aerosols/dust (reference)
  real, parameter, dimension(10) :: plba09=(/  &
       1010.0, 934.0, 854.0, 720.0, 550.0, 390.0, 255.0, 150.0, 70.0, 10.0/)
  !   Layer   1      2      3      4      5      6      7      8     9
  integer :: my_layer
  integer :: l
  real, dimension(nrad,6) :: ataulx
  real, dimension(nrad,6) :: sraext
  integer :: iband
  real :: rhftau

  ! First ascertain the humidity effect
  nrhnan(:,:) = 1
  do k = 1, nrad
     if(rhl(k) > 0.9005)then
        xrh = (rhl(k) - 0.899499) * 1000.0
        nrh = xrh + 90
        if(nrh > 189)nrh=189
     else
        xrh = rhl(k) * 100.5
        nrh = xrh
        if(nrh < 0)nrh = 0
     endif
     do na = 1, 4
        if(kdeliq(k,iw,na) == 0)then
           rhdna = rhdtna(tlm(k),na)
           if(rhl(k) > rhdna)kdeliq(k,iw,na) = 1
        else
           if(rhl(k) < rhc(na))kdeliq(k,iw,na) = 0
        endif
        nrhnan(k,na) = nrh * kdeliq(k,iw,na) + 1
     enddo
  enddo

  ! Now map from the aerosol grid to the model grid
  ! Aerosol latitude bins go from 90S to 90N, 46 points.
  if(glatw(iw) > -88.0)then
     ilat = int(0.25 * (glatw(iw) + 88.0)) + 2
  else
     ilat = 1
  endif
  ilat = max(1,min(46,ilat))
  ! Aerosol longitude bins go from -180 to 180, 72 points.
  ilon = int(0.2 * (glonw(iw) + 180.0)) + 1
  ilon = max(1,min(72,ilon))

  ! Do the vertical mapping and fill the tau array
  do k = 1, nrad
     ! Find layer
     my_layer = 0
     find_layer: do l = 1, 9
        if(pl(k) < plba09(l) .and. pl(k) > plba09(l+1))then
           my_layer = l
           exit find_layer
        endif
     enddo find_layer
     if(my_layer == 0 .and. pl(k) > plba09(1))my_layer = 1
     if(my_layer /= 0)then
        ataulx(k,1:6) = a6jday(my_layer,1:6,ilon,ilat)
     else
        ataulx(k,1:6) = 0.0
     endif
  enddo

  ! (solar bci,bcb components) : 6 wavelength bands:
  !  The nominal Mie scattering spectral band subdivisions are:
  !                       -------------NIR------------      VIS
  !                   L=    1     2     3     4     5        6
  !           WavA (nm)=  2200  1500  1250   860   770      300
  !           WavB (nm)=  4000  2200  1500  1250   860      770

  do k = 1, nrad
     do iband = 1, 6
        sraext(k,iband) = 2.0 * (srbqex(iband,5) * ataulx(k,5) +   &
             srbqex(iband,6) * ataulx(k,6))
     enddo
  enddo

  soot = sraext
  waso(:,:) = 0.0

  do k = 1, nrad
     do iband = 1, 6
        do na = 1, 4
           rhftau = rhinfo(nrhnan(k,na),2,na) * ataulx(k,na)
           sraext(k,iband) = sraext(k,iband) +   &
                srhqex(iband,nrhnan(k,na),na) * rhftau
           if(na /= 2)then
              waso(k,iband) = waso(k,iband) +   &
                   srhqex(iband,nrhnan(k,na),na) * rhftau
           else
              sea_salt(k,iband) = srhqex(iband,nrhnan(k,na),na) * rhftau
           endif

        enddo
     enddo
  enddo

  return
end subroutine getaer

!=============================================================================

real function rhdtna(tk,na)
  implicit none

  integer, intent(in) :: na
  real,  intent(in) :: tk

  real, dimension(4), parameter :: acoef = (/ 0.8, 0.75, 0.62, 0.8 /)
  real, dimension(4), parameter :: bcoef = (/ 25.0, 80.0, 852.0, 25.0 /)

  rhdtna = min(1.0, acoef(na) * exp(bcoef(na) * (298.0 - tk) / (298.0 * tk)))

  return
end function rhdtna

!================================================================
subroutine getdst(iw, nrad, pl, dust)

  use mem_grid, only: glatw, glonw
  use aero_coms, only: ddjday, qxdust

  implicit none

  integer, intent(in) :: iw
  integer :: ilat
  integer :: ilon
  integer :: k
  integer, intent(in) :: nrad
  integer :: my_layer
  integer :: l
  real, dimension(nrad), intent(in) :: pl ! pressure (mb)
  !r PLBA09 Vert. Layering for tropospheric aerosols/dust (reference)
  real, parameter, dimension(10) :: plba09=(/  &
       1010.0, 934.0, 854.0, 720.0, 550.0, 390.0, 255.0, 150.0, 70.0, 10.0/)
  !   Layer   1      2      3      4      5      6      7      8     9
  real, dimension(nrad,8) :: dtaulx
  real, dimension(nrad,6,8), intent(out) :: dust
  integer :: iband

  ! Map the dust grid to the model grid
  ! Dust latitude bins go from 90S to 90N, 46 points.
  if(glatw(iw) > -88.0)then
     ilat = int(0.25 * (glatw(iw) + 88.0)) + 2
  else
     ilat = 1
  endif
  ilat = max(1,min(46,ilat))
  ! Aerosol longitude bins go from -180 to 180, 72 points.
  ilon = int(0.2 * (glonw(iw) + 180.0)) + 1
  ilon = max(1,min(72,ilon))

  ! Do the vertical mapping and fill the tau array
  do k = 1, nrad
     ! Find layer
     my_layer = 0
     find_layer: do l = 1, 9
        if(pl(k) < plba09(l) .and. pl(k) > plba09(l+1))then
           my_layer = l
           exit find_layer
        endif
     enddo find_layer
     if(my_layer == 0 .and. pl(k) > plba09(1))my_layer = 1
     if(my_layer /= 0)then
        dtaulx(k,1:8) = ddjday(my_layer,1:8,ilon,ilat)
     else
        dtaulx(k,1:8) = 0.0
     endif
  enddo

  do k = 1, nrad
     do iband = 1, 6
        do l = 1, 8
           dust(k,iband,l) = qxdust(iband,l) * dtaulx(k,l)
        enddo
     enddo
  enddo

  return
end subroutine getdst
