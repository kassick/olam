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
subroutine makesfc()

! This subroutine generates LAND and SEA files for OLAM runs. 

use mem_grid,    only: nza, nma, nua, nwa, xem, yem, zem, xew, yew, zew, arw0

use misc_coms,   only: io6, runtype

use sea_coms,    only: nms, nus, nws, maxjms, seafile

use leaf_coms,   only: nml, nul, nwl, nzg, nzs, nslcon, nvgcon, maxjml,  &
                       isfcl, ivegflg, landusefile,  &
                       isoilflg, soil_database, veg_database

use mem_mksfc,   only: mksfc, mksfc_sea, mksfc_land,  &
                       mksfc_itab_uls, mksfc_itab_ul, mksfc_itab_us,  &
                       mksfc_itab_wls, mksfc_itab_wl, mksfc_itab_ws

use mem_ijtabs,  only: itab_u, itab_w

use consts_coms, only: erad, piu180

use leaf_db,     only: leaf_database_read

implicit none

integer :: k
integer :: ndims, idims(1)
integer :: j
integer :: im, iu, iw
integer :: ius, iul
integer :: im1ls, im2ls, iw1ls, iw2ls
integer :: nmls, nuls, nwls
integer :: iuls, iwls, imls, jmls
integer :: iws, ims
integer :: iwl, iml
integer :: idatq, datsoil
integer :: nzg0
integer :: nws0
integer :: nwl0
integer :: nmls0
integer :: ivegflg0, isoilflg0

character(len=1) :: dummy
character(len=128) :: flnm

! Define array dimensions for land/sea cells and cell nodes to be those of
! OLAM atmospheric grid as defined in current OLAMIN file.

nmls = nma
nuls = nua
nwls = nwa

! Allocate mksfc arrays for combined land/sea grid

allocate (mksfc_itab_uls(nuls))
allocate (mksfc_itab_wls(nwls))

allocate (mksfc%mpt_sea   (nmls))
allocate (mksfc%mpt_land  (nmls))
allocate (mksfc%upt_sea   (nuls))
allocate (mksfc%upt_land  (nuls))
allocate (mksfc%wpt_sea   (nwls))
allocate (mksfc%wpt_land  (nwls))

allocate (mksfc%idatp     (nwls))
allocate (mksfc%leaf_class(nwls))
allocate (mksfc%area      (nwls))

allocate (mksfc%xemls(nmls))
allocate (mksfc%yemls(nmls))
allocate (mksfc%zemls(nmls))

allocate (mksfc%xewls(nwls))
allocate (mksfc%yewls(nwls))
allocate (mksfc%zewls(nwls))

! Initialize mksfc arrays

do imls = 1,nmls
   im = imls

   mksfc%mpt_sea (imls) = 1
   mksfc%mpt_land(imls) = 1

   mksfc%xemls(imls) = xem(im)
   mksfc%yemls(imls) = yem(im)
   mksfc%zemls(imls) = zem(im)
enddo

do iuls = 2,nuls
   iu = iuls

   mksfc_itab_uls(iuls)%im1 = itab_u(iu)%im1
   mksfc_itab_uls(iuls)%im2 = itab_u(iu)%im2
   mksfc_itab_uls(iuls)%iw1 = itab_u(iu)%iw1
   mksfc_itab_uls(iuls)%iw2 = itab_u(iu)%iw2

   mksfc%upt_sea (iuls) = 1
   mksfc%upt_land(iuls) = 1
enddo

do iwls = 2,nwls

   iw = iwls

   mksfc_itab_wls(iwls)%jm = 3

   mksfc_itab_wls(iwls)%im(1) = itab_w(iw)%im1
   mksfc_itab_wls(iwls)%im(2) = itab_w(iw)%im2
   mksfc_itab_wls(iwls)%im(3) = itab_w(iw)%im3

   mksfc_itab_wls(iwls)%iu(1) = itab_w(iw)%iu1
   mksfc_itab_wls(iwls)%iu(2) = itab_w(iw)%iu2
   mksfc_itab_wls(iwls)%iu(3) = itab_w(iw)%iu3

   mksfc%wpt_sea (iwls) = 1
   mksfc%wpt_land(iwls) = 1

   mksfc%area(iwls) = arw0(iw)

   mksfc%xewls(iwls) = xew(iw)
   mksfc%yewls(iwls) = yew(iw)
   mksfc%zewls(iwls) = zew(iw)
enddo

if (ivegflg == 2) then

! If ivegflg == 2, fill sea/land cell areas and IW values, plus
! leaf_class for land cells, from default value defined in OLAMIN.
! User customization can be done here.

   mksfc%leaf_class(2:nwls) = nvgcon

else

! If ivegflg == 1, fill sea/land cells with leaf_class from leaf database

   call leaf_database_read(nwls,               &
                           nmls,               &
                           mksfc_itab_wls,     &
                           mksfc%xemls,        &
                           mksfc%yemls,        &
                           mksfc%zemls,        &
                           trim(veg_database), &
                           trim(veg_database), &
                           'leaf_class',       &
                           idatp=mksfc%idatp   )

   do iwls = 2,nwls
      call datp_datq(mksfc%idatp(iwls),idatq)
      mksfc%leaf_class(iwls) = idatq
   enddo

endif

! Initialize global sea/land cell counters to 1

nws = 1
nwl = 1

! Loop over all land/sea cells

do iwls = 2,nwls

! Count up individual land & sea cells.
! Flag M and U points for association with land or sea cells.

   if (mksfc%leaf_class(iwls) <= 1) then
      nws = nws + 1

      do j = 1,mksfc_itab_wls(iwls)%jm
         imls = mksfc_itab_wls(iwls)%im(j)
         iuls = mksfc_itab_wls(iwls)%iu(j)
         mksfc%mpt_sea(imls) = 1
         mksfc%upt_sea(iuls) = 1
         mksfc%wpt_sea(iwls) = nws
      enddo

   else  
      nwl = nwl + 1

      do j = 1,mksfc_itab_wls(iwls)%jm
         imls = mksfc_itab_wls(iwls)%im(j)
         iuls = mksfc_itab_wls(iwls)%iu(j)
         mksfc%mpt_land(imls) = 1
         mksfc%upt_land(iuls) = 1
         mksfc%wpt_land(iwls) = nwl
      enddo
   endif

enddo

! Count up land & sea M and U pts

nms = 1
nml = 1

do imls = 2,nmls
   if (mksfc%mpt_sea (imls) == 1) nms = nms + 1
   if (mksfc%mpt_land(imls) == 1) nml = nml + 1
enddo

nus = 1
nul = 1

do iuls = 2,nuls
   if (mksfc%upt_sea (iuls) == 1) nus = nus + 1
   if (mksfc%upt_land(iuls) == 1) nul = nul + 1
enddo

! Now that final nws, nwl, nus, nul, nms, and nml values are determined, 
! allocate separate sea and land cell arrays.  Copy data from temporary 
! arrays to these.

allocate (mksfc_itab_us(nus))
allocate (mksfc_itab_ws(nws))

allocate (mksfc_sea%leaf_class(nws))
allocate (mksfc_sea%area      (nws))

allocate (mksfc_sea%xems(nms))
allocate (mksfc_sea%yems(nms))
allocate (mksfc_sea%zems(nms))

allocate (mksfc_sea%xews(nws))
allocate (mksfc_sea%yews(nws))
allocate (mksfc_sea%zews(nws))

allocate (mksfc_itab_ul(nul))
allocate (mksfc_itab_wl(nwl))

allocate (mksfc_land%leaf_class    (nwl))
allocate (mksfc_land%ntext_soil(nzg,nwl))
allocate (mksfc_land%area          (nwl))

allocate (mksfc_land%xeml(nml))
allocate (mksfc_land%yeml(nml))
allocate (mksfc_land%zeml(nml))

allocate (mksfc_land%xewl(nwl))
allocate (mksfc_land%yewl(nwl))
allocate (mksfc_land%zewl(nwl))

! Convert land & sea M pt flags to M pt indices.  Copy earth coordinate values 
! from combined land + sea arrays to individual land & sea arrays.

ims = 1
iml = 1

do imls = 2,nmls
   if (mksfc%mpt_sea(imls) == 1) then
      ims = ims + 1

      mksfc%mpt_sea(imls) = ims
      
      mksfc_sea%xems(ims) = mksfc%xemls(imls)
      mksfc_sea%yems(ims) = mksfc%yemls(imls)
      mksfc_sea%zems(ims) = mksfc%zemls(imls)
   endif

   if (mksfc%mpt_land(imls) == 1) then
      iml = iml + 1

      mksfc%mpt_land(imls) = iml

      mksfc_land%xeml(iml) = mksfc%xemls(imls)
      mksfc_land%yeml(iml) = mksfc%yemls(imls)
      mksfc_land%zeml(iml) = mksfc%zemls(imls)
   endif
enddo

! Convert land & sea U pt flags to U pt indices.

ius = 1
iul = 1

do iuls = 2,nuls
   if (mksfc%upt_sea(iuls) == 1) then
      ius = ius + 1

      mksfc%upt_sea(iuls) = ius
   endif

   if (mksfc%upt_land(iuls) == 1) then
      iul = iul + 1

      mksfc%upt_land(iuls) = iul
   endif
enddo

! Copy data from combined land + sea arrays to individual land & sea arrays.

do iwls = 2,nwls

   if (mksfc%leaf_class(iwls) <= 1) then

      iws = mksfc%wpt_sea(iwls)

      mksfc_sea%leaf_class(iws) = mksfc%leaf_class(iwls)
      mksfc_sea%area      (iws) = mksfc%area      (iwls)

      mksfc_itab_ws(iws)%jm = mksfc_itab_wls(iwls)%jm

      do j = 1,mksfc_itab_ws(iws)%jm
         imls = mksfc_itab_wls(iwls)%im(j)
         iuls = mksfc_itab_wls(iwls)%iu(j)
         mksfc_itab_ws(iws)%im(j) = mksfc%mpt_sea(imls)
         mksfc_itab_ws(iws)%iu(j) = mksfc%upt_sea(iuls)
      enddo

      mksfc_sea%xews(iws) = mksfc%xewls(iwls)
      mksfc_sea%yews(iws) = mksfc%yewls(iwls)
      mksfc_sea%zews(iws) = mksfc%zewls(iwls)

   else

      iwl = mksfc%wpt_land(iwls)

      mksfc_land%leaf_class(iwl) = mksfc%leaf_class(iwls)
      mksfc_land%area      (iwl) = mksfc%area      (iwls)

      mksfc_itab_wl(iwl)%jm = mksfc_itab_wls(iwls)%jm

      do j = 1,mksfc_itab_wl(iwl)%jm
         imls = mksfc_itab_wls(iwls)%im(j)
         iuls = mksfc_itab_wls(iwls)%iu(j)
         mksfc_itab_wl(iwl)%im(j) = mksfc%mpt_land(imls)
         mksfc_itab_wl(iwl)%iu(j) = mksfc%upt_land(iuls)
      enddo

      mksfc_land%xewl(iwl) = mksfc%xewls(iwls)
      mksfc_land%yewl(iwl) = mksfc%yewls(iwls)
      mksfc_land%zewl(iwl) = mksfc%zewls(iwls)

   endif

enddo

! Copy U pt neighbor data from combined land + sea arrays to individual 
! land & sea arrays.

do iuls = 2,nuls

   im1ls = mksfc_itab_uls(iuls)%im1
   im2ls = mksfc_itab_uls(iuls)%im2
   iw1ls = mksfc_itab_uls(iuls)%iw1
   iw2ls = mksfc_itab_uls(iuls)%iw2

   ius = mksfc%upt_sea(iuls)
   iul = mksfc%upt_land(iuls)
   
   if (ius > 1) then
      mksfc_itab_us(ius)%im1 = mksfc%mpt_sea(im1ls)
      mksfc_itab_us(ius)%im2 = mksfc%mpt_sea(im2ls)
      mksfc_itab_us(ius)%iw1 = mksfc%wpt_sea(iw1ls)
      mksfc_itab_us(ius)%iw2 = mksfc%wpt_sea(iw2ls)
   endif

   if (iul > 1) then
      mksfc_itab_ul(iul)%im1 = mksfc%mpt_land(im1ls)
      mksfc_itab_ul(iul)%im2 = mksfc%mpt_land(im2ls)
      mksfc_itab_ul(iul)%iw1 = mksfc%wpt_land(iw1ls)
      mksfc_itab_ul(iul)%iw2 = mksfc%wpt_land(iw2ls)
   endif

enddo

! Initialize soil textural class

if (isoilflg == 2) then

! If soilflg == 2, fill land cells with default horizontally homogeneous
! soil textural class value defined in OLAMIN.
! User customization can be done here.

   mksfc_land%ntext_soil(1:nzg,1:nwl) = nslcon

else
   
! If soilflg == 2, read soil textural class from database

   call leaf_database_read(nwl,                 &
                           nml,                 &
                           mksfc_itab_wl,       &
                           mksfc_land%xeml,     &
                           mksfc_land%yeml,     &
                           mksfc_land%zeml,     &
                           trim(soil_database), &
                           trim(soil_database), &
                           'soil_text',         &
                           idatp=mksfc%idatp    )

! Loop over all land cells (already defined and filled with leaf_class)

   do iwl = 2,nwl

      call datp_datsoil(mksfc%idatp(iwl),datsoil)

! For now, assign single-level FAO textural class to all soil layers.

      do k = 1,nzg
         mksfc_land%ntext_soil(k,iwl) = datsoil
      enddo

   enddo

endif

! Write land file (contains LEAF CLASS, SOIL TEXTURAL CLASS, and grid info)   
write(io6,'(/,a)') 'calling landfile_write'

call landfile_write()

! Write sea file (contains LEAF CLASS and grid info)   

write(io6,'(/,a)') 'calling seafile_write'

call seafile_write()

write(io6,'(/,a)') 'called seafile_write'

! Deallocate mksfc arrays

deallocate (mksfc_itab_uls)
deallocate (mksfc_itab_wls)

deallocate (mksfc%mpt_sea)
deallocate (mksfc%mpt_land)
deallocate (mksfc%upt_sea)
deallocate (mksfc%upt_land)
deallocate (mksfc%wpt_sea)
deallocate (mksfc%wpt_land)

deallocate (mksfc%idatp)
deallocate (mksfc%leaf_class)
deallocate (mksfc%area)

deallocate (mksfc%xemls)
deallocate (mksfc%yemls)
deallocate (mksfc%zemls)

deallocate (mksfc%xewls)
deallocate (mksfc%yewls)
deallocate (mksfc%zewls)

deallocate (mksfc_itab_us)
deallocate (mksfc_itab_ws)

deallocate (mksfc_sea%leaf_class)
deallocate (mksfc_sea%area)

deallocate (mksfc_sea%xems)
deallocate (mksfc_sea%yems)
deallocate (mksfc_sea%zems)

deallocate (mksfc_sea%xews)
deallocate (mksfc_sea%yews)
deallocate (mksfc_sea%zews)

deallocate (mksfc_itab_ul)
deallocate (mksfc_itab_wl)

deallocate (mksfc_land%leaf_class)
deallocate (mksfc_land%ntext_soil)
deallocate (mksfc_land%area)

deallocate (mksfc_land%xeml)
deallocate (mksfc_land%yeml)
deallocate (mksfc_land%zeml)

deallocate (mksfc_land%xewl)
deallocate (mksfc_land%yewl)
deallocate (mksfc_land%zewl)

return
end subroutine makesfc

!==========================================================================

subroutine datp_datq(idatp,idatq)

! This subroutine maps the input idatp classes to a smaller set idatq
! which represents the full set of LEAF-2 or LEAF-3 classes for which 
! LSP values are
! defined.

implicit none

integer, intent(in) :: idatp
integer, intent(out) :: idatq

integer :: catb(0:100)
integer :: catb_leaf3(0:100)

!  Olson Global Ecosystems dataset OGE_2 (96 classes) mapped to LEAF-3 classes
!  (see leaf3_document).

!-------------------------------------------!
data catb/ 0,                             & !
          19, 8, 4, 5, 6, 7, 9, 3,11,16,  & !  0
          10, 2,17, 1, 0,12,13,14,18, 4,  & ! 10
           4, 4,14,14, 6, 6, 4, 7, 7,15,  & ! 20
          15, 6, 7, 7,15,16,16,16,16, 8,  & ! 30
           8, 8,18,17,17,12,12, 7,10, 3,  & ! 40
          10,10,11,14,18,18,18,18,13, 6,  & ! 50
           5, 4,11,12, 0, 0, 0, 0, 3, 2,  & ! 60
           3,20, 0,17,17,17, 4,14, 7, 3,  & ! 70
           3, 3, 3, 3, 3, 3, 8,12, 7, 6,  & ! 80
          18,15,15,15, 4, 5, 0, 0, 0, 0   / ! 90  ! 97 & 98 not used
!-------------------------------------------!     ! 99 is Goode Homolosine 
                                                  !    empty space
!          1  2  3  4  5  6  7  8  9 10           ! 100 is missing data
                                                  ! Map all of these to ocean (idatq=0)
idatq = catb(idatp)

return
end subroutine datp_datq

!==========================================================================

subroutine datp_datsoil(idatp,idatsoil)

! This subroutine maps the input idatp soil classes to a smaller set idatsoil
! which represents the full set of LEAF-2 classes for which soil parameter
! values are defined.

implicit none

integer, intent(in) :: idatp
integer, intent(out) :: idatsoil

integer :: catb(0:133)

! (Bob 9/14/2000) This table maps FAO soil units numbered 0 to 132, plus our 
! own missing value designated 133, to the USDA soil textural classes.  FAO 
! classes [0] (ocean), [1, 7, 27, 42, 66, 69, 77, 78, 84, 98, 113] (soils not 
! used in original FAO dataset), [132] (water), and [133] (our missing value) 
! are all mapped to a default class of sandy clay loam in case they happen to
! correspond to a land surface area in the landuse dataset that RAMS uses to
! define land area.  We wrote missing value class 133 to the RAMS FAO files
! whenever a negative value (which is out of range of defined values) was
! found in the original FAO dataset, which occurred in about 2.6% of the
! pixels.  For the remaining FAO classes, a cross reference table to Zobler 
! soil texture classes that was provided, plus our own cross referencing table
! from Zobler to USDA classes listed below, provides the mapping from FAO to 
! USDA.  In this mapping, we use only organic USDA classes and omit nonorganic
! classes since most soils do contain organic matter, and organic content 
! information is not provided in the Zobler classes.

!  Zobler Class              USDA Class

!  1  Coarse                 2  Loamy sand
!  2  Medium                 4  Silt loam
!  3  Fine                   8  Clay loam
!  4  Coarse-medium          3  Sandy loam
!  5  Coarse-fine            6  Sandy clay loam
!  6  Medium-fine            7  Silty clay loam
!  7  Coarse-medium-fine     6  Sandy clay loam
!  8  Organic matter         5  Loam

!                            1  Sand (not used)
!                            9  Sandy clay (not used)
!                           10  Silty clay (not used)
!                           11  Clay (not used)
!                           12  Peat (not used)

!-------------------------------------------!
data catb/ 6,                             & !
           6, 4, 4, 7, 7, 8, 6, 4, 4, 4,  & !   0
           7, 4, 4, 4, 8, 4, 8, 4, 4, 8,  & !  10
           4, 2, 4, 4, 4, 4, 6, 8, 8, 8,  & !  20
           4, 8, 8, 2, 6, 4, 7, 4, 4, 3,  & !  30
           4, 6, 7, 4, 4, 4, 4, 4, 4, 4,  & !  40
           4, 4, 4, 4, 4, 4, 2, 4, 4, 2,  & !  50
           4, 3, 4, 2, 7, 6, 4, 4, 6, 8,  & !  60
           8, 7, 2, 5, 4, 5, 6, 6, 4, 2,  & !  70
           2, 2, 4, 6, 2, 2, 2, 2, 2, 4,  & !  80
           2, 2, 2, 4, 2, 4, 3, 6, 2, 7,  & !  90
           4, 4, 4, 8, 8, 8, 3, 7, 4, 4,  & ! 100
           4, 3, 6, 4, 2, 4, 4, 4, 2, 2,  & ! 110
           2, 4, 6, 4, 4, 7, 7, 6, 3, 2,  & ! 120
           2, 6, 6                        / ! 130
!-------------------------------------------!
!          1  2  3  4  5  6  7  8  9 10

idatsoil = catb(idatp)

return
end subroutine datp_datsoil

!==========================================================================

subroutine landfile_write()

use leaf_coms, only: nzg, nml, nul, nwl, maxjml, landusefile, slz

use mem_mksfc,  only: mksfc_land, mksfc_itab_ul, mksfc_itab_wl

use hdf5_utils, only: shdf5_open, shdf5_orec, shdf5_close
use misc_coms,  only: io6, ngrids, nxp,  &
                      centlat, centlon, grdlen, grdwid, grdaxis

implicit none

integer :: ndims
integer :: idims(2)
integer, save :: iclobber1 = 1

character(len=128) :: flnm

integer :: iwl
integer :: iscr(maxjml,nwl)

! Write land file header information

flnm=trim(landusefile)//'-S'//'.h5'

call shdf5_open(flnm,'W',iclobber1)

ndims = 1
idims(1) = 1

print*, 'landfile_write1 ',nml,nul,nwl,nzg


call shdf5_orec(ndims,idims,'nzg'   ,ivars=nzg)
call shdf5_orec(ndims,idims,'nml'   ,ivars=nml)
call shdf5_orec(ndims,idims,'nul'   ,ivars=nul)
call shdf5_orec(ndims,idims,'nwl'   ,ivars=nwl)
call shdf5_orec(ndims,idims,'ngrids',ivars=ngrids)
call shdf5_orec(ndims,idims,'nxp'   ,ivars=nxp)

idims(1) = ngrids

call shdf5_orec(ndims,idims,'centlat',rvara=centlat)

call shdf5_orec(ndims,idims,'centlon',rvara=centlon)
call shdf5_orec(ndims,idims,'grdlen',rvara=grdlen)
call shdf5_orec(ndims,idims,'grdwid',rvara=grdwid)
call shdf5_orec(ndims,idims,'grdaxis',rvara=grdaxis)
idims(1) = nzg

call shdf5_orec(ndims,idims,'slz',rvara=slz)

! Write arrays to land file

ndims = 1
idims(1) = nwl

call shdf5_orec(ndims,idims,'land_area' ,rvara=mksfc_land%area)
call shdf5_orec(ndims,idims,'leaf_class',ivara=mksfc_land%leaf_class)
call shdf5_orec(ndims,idims,'jml'       ,ivara=mksfc_itab_wl(:)%jm)
call shdf5_orec(ndims,idims,'xewl'      ,rvara=mksfc_land%xewl)
call shdf5_orec(ndims,idims,'yewl'      ,rvara=mksfc_land%yewl)
call shdf5_orec(ndims,idims,'zewl'      ,rvara=mksfc_land%zewl)

idims(1) = nul

call shdf5_orec(ndims,idims,'itab_ul%im1',ivara=mksfc_itab_ul(:)%im1)
call shdf5_orec(ndims,idims,'itab_ul%im2',ivara=mksfc_itab_ul(:)%im2)
call shdf5_orec(ndims,idims,'itab_ul%iw1',ivara=mksfc_itab_ul(:)%iw1)
call shdf5_orec(ndims,idims,'itab_ul%iw2',ivara=mksfc_itab_ul(:)%iw2)

idims(1) = nml

call shdf5_orec(ndims,idims,'xeml',rvara=mksfc_land%xeml)
call shdf5_orec(ndims,idims,'yeml',rvara=mksfc_land%yeml)
call shdf5_orec(ndims,idims,'zeml',rvara=mksfc_land%zeml)

ndims = 2
idims(1) = nzg
idims(2) = nwl

call shdf5_orec(ndims,idims,'ntext_soil',ivara=mksfc_land%ntext_soil)

ndims = 2
idims(1) = maxjml
idims(2) = nwl

! Copy itab_wl%im to scratch array for output

do iwl = 1,nwl
   iscr(1:maxjml,iwl) = mksfc_itab_wl(iwl)%im(1:maxjml)
enddo

call shdf5_orec(ndims,idims,'itab_wl%im',ivara=iscr)

! Copy itab_wl%iu to scratch array for output

do iwl = 1,nwl
   iscr(1:maxjml,iwl) = mksfc_itab_wl(iwl)%iu(1:maxjml)
enddo

call shdf5_orec(ndims,idims,'itab_wl%iu',ivara=iscr)

call shdf5_close()

return
end subroutine landfile_write

!==========================================================================

subroutine seafile_write()

use sea_coms,  only: nms, nus, nws, maxjms, seafile

use mem_mksfc,  only: mksfc_sea, mksfc_itab_us, mksfc_itab_ws

use hdf5_utils, only: shdf5_open, shdf5_orec, shdf5_close
use misc_coms,  only: io6, ngrids, nxp,  &
                      centlat, centlon, grdlen, grdwid, grdaxis

implicit none

integer :: ndims, idims(2)
integer, save :: iclobber1 = 1

character(len=128) :: flnm

integer :: iws
integer :: iscr(maxjms,nws)

! Write sea file header information

flnm=trim(seafile)//'-S'//'.h5'

call shdf5_open(flnm,'W',iclobber1)

ndims = 1
idims(1) = 1

call shdf5_orec(ndims,idims,'nms'   ,ivars=nms)
call shdf5_orec(ndims,idims,'nus'   ,ivars=nus)
call shdf5_orec(ndims,idims,'nws'   ,ivars=nws)
call shdf5_orec(ndims,idims,'ngrids',ivars=ngrids)
call shdf5_orec(ndims,idims,'nxp'   ,ivars=nxp)

idims(1) = ngrids

call shdf5_orec(ndims,idims,'centlat',rvara=centlat)
call shdf5_orec(ndims,idims,'centlon',rvara=centlon)
call shdf5_orec(ndims,idims,'grdlen' ,rvara=grdlen)
call shdf5_orec(ndims,idims,'grdwid' ,rvara=grdwid)
call shdf5_orec(ndims,idims,'grdaxis',rvara=grdaxis)

! Write arrays to sea file

ndims = 1
idims(1) = nws

call shdf5_orec(ndims,idims,'sea_area',rvara=mksfc_sea%area)
call shdf5_orec(ndims,idims,'jms'     ,ivara=mksfc_itab_ws(:)%jm)
call shdf5_orec(ndims,idims,'xews'    ,rvara=mksfc_sea%xews)
call shdf5_orec(ndims,idims,'yews'    ,rvara=mksfc_sea%yews)
call shdf5_orec(ndims,idims,'zews'    ,rvara=mksfc_sea%zews)

idims(1) = nus

call shdf5_orec(ndims,idims,'itab_us%im1',ivara=mksfc_itab_us(:)%im1)
call shdf5_orec(ndims,idims,'itab_us%im2',ivara=mksfc_itab_us(:)%im2)
call shdf5_orec(ndims,idims,'itab_us%iw1',ivara=mksfc_itab_us(:)%iw1)
call shdf5_orec(ndims,idims,'itab_us%iw2',ivara=mksfc_itab_us(:)%iw2)

idims(1) = nms

call shdf5_orec(ndims,idims,'xems',rvara=mksfc_sea%xems)
call shdf5_orec(ndims,idims,'yems',rvara=mksfc_sea%yems)
call shdf5_orec(ndims,idims,'zems',rvara=mksfc_sea%zems)

ndims = 2
idims(1) = maxjms
idims(2) = nws

! Copy itab_ws%im to temporary array for output

do iws = 1,nws
   iscr(1:maxjms,iws) = mksfc_itab_ws(iws)%im(1:maxjms)
enddo

call shdf5_orec(ndims,idims,'itab_ws%im',ivara=iscr)

! Copy itab_ws%iu to scratch array for output

do iws = 1,nws
   iscr(1:maxjms,iws) = mksfc_itab_ws(iws)%iu(1:maxjms)
enddo

call shdf5_orec(ndims,idims,'itab_ws%iu',ivara=iscr)

call shdf5_close()

return
end subroutine seafile_write

!==========================================================================

subroutine landfile_read()

use leaf_coms,  only: nzg, nml, nul, nwl, mml, mul, mwl, maxjml, landusefile, slz

use mem_leaf,   only: land, itab_ul, itab_wl, alloc_land_grid

use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
use misc_coms,  only: io6, ngrids, nxp,  &
                      centlat, centlon, grdlen, grdwid, grdaxis

implicit none

integer :: ndims
integer :: idims(2)
integer, save :: iclobber1 = 1

character(len=128) :: flnm

logical there

integer :: iwl
integer, allocatable :: iscr(:,:)
real, external :: walltime

!-------------------------------------------------------------------------------
! STEP 1: Open land file and read 4 array dimensions
!-------------------------------------------------------------------------------

flnm = trim(landusefile)//'-S'//'.h5'

write(io6,*) 'Checking leaf file ',trim(flnm)

inquire(file=flnm,exist=there)

if (.not. there) then
   write(io6,*) 'Land file was not found - stopping run'
   stop 'stop: no land file'
endif

call shdf5_open(flnm,'R')

ndims = 1
idims(1) = 1

call shdf5_irec(ndims,idims,'nml',ivars=nml)
call shdf5_irec(ndims,idims,'nul',ivars=nul)
call shdf5_irec(ndims,idims,'nwl',ivars=nwl)
call shdf5_irec(ndims,idims,'nzg',ivars=nzg)

write(io6, '(/,a)')   '==============================================='
write(io6, '(a)')     'Reading from land file:'
write(io6, '(a,4i8)') '  nml, nul, nwl, nzg = ', nml, nul, nwl, nzg
write(io6, '(a,/)')   '==============================================='

!-------------------------------------------------------------------------------
! STEP 2: Allocate land grid arrays
!-------------------------------------------------------------------------------

call alloc_land_grid(nml,nul,nwl,nzg,1)
allocate (iscr(maxjml,nwl))

!-------------------------------------------------------------------------------
! STEP 3: Read arrays from land file
!-------------------------------------------------------------------------------

idims(1) = nzg

call shdf5_irec(ndims,idims,'slz',rvara=slz)

! Write arrays to land file

ndims = 1
idims(1) = nwl

call shdf5_irec(ndims,idims,'land_area' ,rvara=land%area)
call shdf5_irec(ndims,idims,'leaf_class',ivara=land%leaf_class)
call shdf5_irec(ndims,idims,'jml'       ,ivara=itab_wl(:)%jm)
call shdf5_irec(ndims,idims,'xewl'      ,rvara=land%xewl)
call shdf5_irec(ndims,idims,'yewl'      ,rvara=land%yewl)
call shdf5_irec(ndims,idims,'zewl'      ,rvara=land%zewl)

idims(1) = nul

call shdf5_irec(ndims,idims,'itab_ul%im1',ivara=itab_ul(:)%im1)
call shdf5_irec(ndims,idims,'itab_ul%im2',ivara=itab_ul(:)%im2)
call shdf5_irec(ndims,idims,'itab_ul%iw1',ivara=itab_ul(:)%iw1)
call shdf5_irec(ndims,idims,'itab_ul%iw2',ivara=itab_ul(:)%iw2)

idims(1) = nml

call shdf5_irec(ndims,idims,'xeml',rvara=land%xeml)
call shdf5_irec(ndims,idims,'yeml',rvara=land%yeml)
call shdf5_irec(ndims,idims,'zeml',rvara=land%zeml)

ndims = 2
idims(1) = nzg
idims(2) = nwl

call shdf5_irec(ndims,idims,'ntext_soil',ivara=land%ntext_soil)

ndims = 2
idims(1) = maxjml
idims(2) = nwl

call shdf5_irec(ndims,idims,'itab_wl%im',ivara=iscr)

! Copy input scratch array to itab_wl%im

do iwl = 1,nwl
   itab_wl(iwl)%im(1:maxjml) = iscr(1:maxjml,iwl)
enddo

call shdf5_irec(ndims,idims,'itab_wl%iu',ivara=iscr)

! Copy input scratch array to itab_wl%iu

do iwl = 1,nwl
   itab_wl(iwl)%iu(1:maxjml) = iscr(1:maxjml,iwl)
enddo

call shdf5_close()

mml = nml
mul = nul
mwl = nwl

return
end subroutine landfile_read

!==========================================================================

subroutine seafile_read()

use sea_coms,   only: nms, nus, nws, mms, mus, mws, maxjms, seafile
use mem_sea,    only: sea, itab_us, itab_ws, alloc_sea_grid
use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_close
use misc_coms,  only: io6

implicit none

integer :: ndims, idims(2)

character(len=128) :: flnm

logical :: there

integer :: iws
integer, allocatable :: iscr(:,:)
real, external :: walltime


!-------------------------------------------------------------------------------
! STEP 1: Open SEAFILE and read 3 array dimensions
!-------------------------------------------------------------------------------

flnm = trim(seafile)//'-S'//'.h5'

write(io6,*)  'Checking sea file ',trim(flnm)

inquire(file=flnm,exist=there)

if (.not. there) then
   write(io6,*) 'SEA file was not found - stopping run'
   stop 'stop: no sea file'
endif

call shdf5_open(flnm,'R')

ndims = 1
idims(1) = 1

call shdf5_irec(ndims,idims,'nms',ivars=nms)
call shdf5_irec(ndims,idims,'nus',ivars=nus)
call shdf5_irec(ndims,idims,'nws',ivars=nws)

write(io6, '(/,a)')   '====================================='
write(io6, '(a)')     'Reading from sea file:'
write(io6, '(a,3i8)') '  nms, nus, nws = ', nms, nus, nws
write(io6, '(a,/)')   '====================================='

!-------------------------------------------------------------------------------
! STEP 2: Allocate SEA grid arrays
!-------------------------------------------------------------------------------

call alloc_sea_grid(nms,nus,nws,1)
allocate (iscr(maxjms,nws))

!-------------------------------------------------------------------------------
! STEP 3: Read arrays from sea file
!-------------------------------------------------------------------------------

ndims = 1
idims(1) = nws

call shdf5_irec(ndims,idims,'sea_area',rvara=sea%area)
call shdf5_irec(ndims,idims,'jms'     ,ivara=itab_ws(:)%jm)
call shdf5_irec(ndims,idims,'xews'    ,rvara=sea%xews)
call shdf5_irec(ndims,idims,'yews'    ,rvara=sea%yews)
call shdf5_irec(ndims,idims,'zews'    ,rvara=sea%zews)

idims(1) = nus

call shdf5_irec(ndims,idims,'itab_us%im1',ivara=itab_us(:)%im1)
call shdf5_irec(ndims,idims,'itab_us%im2',ivara=itab_us(:)%im2)
call shdf5_irec(ndims,idims,'itab_us%iw1',ivara=itab_us(:)%iw1)
call shdf5_irec(ndims,idims,'itab_us%iw2',ivara=itab_us(:)%iw2)

idims(1) = nms

call shdf5_irec(ndims,idims,'xems',rvara=sea%xems)
call shdf5_irec(ndims,idims,'yems',rvara=sea%yems)
call shdf5_irec(ndims,idims,'zems',rvara=sea%zems)

ndims = 2
idims(1) = maxjms
idims(2) = nws

iscr(:,:) = 0

call shdf5_irec(ndims,idims,'itab_ws%im',ivara=iscr)

! Copy temporary scratch array to itab_ws%im

do iws = 1,nws
   itab_ws(iws)%im(1:maxjms) = iscr(1:maxjms,iws)
enddo

call shdf5_irec(ndims,idims,'itab_ws%iu',ivara=iscr)

! Copy temporary scratch array to itab_ws%iu

do iws = 1,nws
   itab_ws(iws)%iu(1:maxjms) = iscr(1:maxjms,iws)
enddo

call shdf5_close()

mms = nms
mus = nus
mws = nws

return
end subroutine seafile_read
