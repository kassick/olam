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
Module leaf_db

Contains

subroutine leaf_database_read(mwl,mml,itab_wl,xeml,yeml,zeml,  &
                              ofn,ofn2,iaction,idatp,datp)

use consts_coms, only: erad, piu180
use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec
use misc_coms,   only: io6
use mem_leaf,    only: itab_wl_vars

implicit none

integer, intent(in) :: mwl
integer, intent(in) :: mml

type (itab_wl_vars), intent(in) :: itab_wl(mwl)

real, intent(in) :: xeml(mml), yeml(mml), zeml(mml)

character(len=*), intent(in) :: ofn,ofn2
character(len=*), intent(in) :: iaction

integer, optional, intent(out) :: idatp(mwl)
real,    optional, intent(out) :: datp(mwl)

! automatic arrays

integer :: ptable(mwl)
real    :: glat(mwl), glon(mwl)

integer :: nio, njo
integer :: niosh, njosh
integer :: isoc, iwoc
integer :: isocpt, isocpo
integer :: iwocph, iwocpt, iwocpo
integer :: io, jo
integer :: ind
integer :: ind1, ind2
integer :: io1, io2, jo1, jo2
integer :: ifiles, jfiles
integer :: ifile, jfile
integer :: io_full, jo_full
integer :: jtime
integer :: iwl
integer :: im, jm
integer :: nperdeg
integer :: ndims, idims(2)

real :: offpix 
real :: offlat, offlon
real :: wlat1, wlon1
real :: wio1, wio2, wjo1, wjo2
real :: rio_full, rjo_full
real :: xewl, yewl, zewl

integer(kind=4), allocatable :: idato (:,:)
real(kind=4), allocatable :: dato (:,:)

integer, allocatable :: nump    (:,:)
integer, allocatable :: numpind1(:,:)
integer, allocatable :: numpind2(:,:) ! (ifiles,jfiles)

character :: title1*3
character :: title2*4
character :: title3*80

character(len=1) :: dummy

character(len=80) :: flnm
character(len=80) :: line
character(len=80) :: filemonth
character(len=80) :: line2

logical l1,l2

! Open, read, and close OGE, FAO, or NDVI dataset header file

call rams_f_open(29,trim(ofn)//'HEADER','FORMATTED','OLD','READ',0)
read(29,*) nio,njo,nperdeg
close(29)

! Compute number of pixels in a shift to adjacent file (niosh and njosh).
! Compute number of files in database that span all latitudes and
! longitudes on earth. [This algorithm will change when multiple resolutions
! of the SRTM data become available.]

if (mod(nio,nperdeg) == 2) then
   offpix = .5
   niosh = nio - 2
   njosh = njo - 2
else
   offpix = 0.
   niosh = nio - 1
   njosh = njo - 1
endif

! Compute number of files in input dataset that span all latitudes and
! longitudes on earth.  

ifiles = 360 * nperdeg / niosh
jfiles = 180 * nperdeg / njosh
   
! Allocate 5 arrays.

allocate (nump    (ifiles,jfiles))
allocate (numpind1(ifiles,jfiles))
allocate (numpind2(ifiles,jfiles))
allocate (idato   (nio,njo))
allocate (dato    (nio,njo))

do jfile = 1,jfiles
   do ifile = 1,ifiles
      nump(ifile,jfile) = 0
   enddo
enddo

! Get file index (ifile,jfile) within full dataset and count number of p 
! points (nump) that occur in each file

do iwl = 2,mwl

! Determine latitude & longitude of land/sea cells

   xewl = 0.
   yewl = 0.
   zewl = 0.

   do jm = 1,itab_wl(iwl)%jm
      im = itab_wl(iwl)%im(jm)

      xewl = xewl + xeml(im)
      yewl = yewl + yeml(im)
      zewl = zewl + zeml(im)
   enddo

   xewl = xewl / itab_wl(iwl)%jm
   yewl = yewl / itab_wl(iwl)%jm
   zewl = zewl / itab_wl(iwl)%jm

   glat(iwl) = asin(zewl / erad) * piu180
   glon(iwl) = atan2(yewl,xewl) * piu180

   wlat1 = max(-89.9999,min(89.9999,glat(iwl)))
   wlon1 = glon(iwl)

   if (wlon1 >=  180.) wlon1 = wlon1 - 360.
   if (wlon1 <  -180.) wlon1 = wlon1 + 360.
   wlon1 = max(-179.9999,min(179.9999,wlon1))

   rio_full = (wlon1 + 180.) * nperdeg ! must ignore pixel offset here
   rjo_full = (wlat1 +  90.) * nperdeg ! must ignore pixel offset here

   io_full = int(rio_full)
   jo_full = int(rjo_full)

   ifile = io_full / niosh + 1
   jfile = jo_full / njosh + 1

   nump(ifile,jfile) = nump(ifile,jfile) + 1  ! Summation to # pts 
                               ! filled by file (ifile,jfile) in dataset

enddo

! Set up array index values for ptable array

ind = 1
do jfile = 1,jfiles
   do ifile = 1,ifiles
      numpind1(ifile,jfile) = ind
      numpind2(ifile,jfile) = ind
      ind = ind + nump(ifile,jfile)
   enddo
enddo

! Fill ptable array

do iwl = 2,mwl

   wlat1 = max(-89.9999,min(89.9999,glat(iwl)))
   wlon1 = glon(iwl)

   if (wlon1 >=  180.) wlon1 = wlon1 - 360.
   if (wlon1 <  -180.) wlon1 = wlon1 + 360.
   wlon1 = max(-179.9999,min(179.9999,wlon1))

   rio_full = (wlon1 + 180.) * nperdeg ! must ignore pixel offset here
   rjo_full = (wlat1 +  90.) * nperdeg ! must ignore pixel offset here

   io_full = int(rio_full)
   jo_full = int(rjo_full)

   ifile = io_full / niosh + 1
   jfile = jo_full / njosh + 1

   ind = numpind2(ifile,jfile)
   ptable(ind) = iwl
   numpind2(ifile,jfile) = numpind2(ifile,jfile) + 1 ! Sums to ending index
         
enddo

! Read files and extract data

do jfile = 1,jfiles
   do ifile = 1,ifiles
   
      ind1 = numpind1(ifile,jfile)
      ind2 = numpind2(ifile,jfile)
   
      if (ind2 > ind1) then
         iwoc = (ifile - 1) * niosh / nperdeg - 180.  ! SW longitude of current file
         isoc = (jfile - 1) * njosh / nperdeg -  90.  ! SW latitude of current file

! Construct filename

         isocpt = abs(isoc) / 10
         isocpo = abs(isoc) - isocpt*10
         iwocph = abs(iwoc) / 100
         iwocpt = (abs(iwoc) - iwocph * 100) / 10
         iwocpo = abs(iwoc) - iwocph * 100 - iwocpt * 10
         
         if (isoc >= 0) then
            write(title1,'(2i1,a1)') isocpt,isocpo,'N'
         else
            write(title1,'(2i1,a1)') isocpt,isocpo,'S'
         endif
         
         if (iwoc >= 0) then
            write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'E'
         else
            write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'W'
         endif

         title3 = trim(ofn2)//title1//title2//'.h5'
             
         inquire(file=trim(title3),exist=l1,opened=l2)

         write(io6,*) 'gf ',ifile,ifiles,jfile,jfiles,trim(title3)    


! Read file

         if (l1) then
            write(io6,*) 'getting file ',trim(title3)    

            call shdf5_open(trim(title3),'R')

            ndims = 2
            idims(1) = nio
            idims(2) = njo

            if (trim(iaction) == 'leaf_class') then
               call shdf5_irec(ndims,idims,'oge2',ivara=idato)
            elseif (trim(iaction) == 'soil_text') then
               call shdf5_irec(ndims,idims,'fao',ivara=idato)
            elseif (trim(iaction) == 'ndvi') then
               call shdf5_irec(ndims,idims,'ndvi',rvara=dato)
            else
               write(io6,*) 'incorrect action specified in leaf_database'
               write(io6,*) 'stopping run'
               stop 'stop landuse_input1'
            endif

            call shdf5_close()
         else
            write(io6,*) 'In landuse_input, ',iaction,' file is missing'
            write(io6,*) 'Filename = ',trim(title3)
            write(io6,*) 'Stopping model run'
            stop 'stop_landuse_input2'
         endif

         do ind = ind1,ind2-1

            iwl = ptable(ind)         

            wlat1 = max(-89.9999,min(89.9999,glat(iwl)))
            wlon1 = glon(iwl)
            
            if (wlon1 >=  180.) wlon1 = wlon1 - 360.
            if (wlon1 <= -180.) wlon1 = wlon1 + 360.

            rio_full = (wlon1 + 180.) * nperdeg ! must ignore pixel offset here
            rjo_full = (wlat1 +  90.) * nperdeg ! must ignore pixel offset here

            io_full = int(rio_full)
            jo_full = int(rjo_full)

            io1 = mod(io_full,niosh) + 1
            jo1 = mod(jo_full,njosh) + 1

            wio2 = rio_full - float(io_full) + offpix
            wjo2 = rjo_full - float(jo_full) + offpix
           
! At this point, io1, jo1, wio2, and wjo2 are correct if offpix = 0,
! but need correction if offpix = .5

            if (wio2 > 1.) then
               wio2 = wio2 - 1.
               io1 = io1 + 1
            endif
            
            if (wjo2 > 1.) then
               wjo2 = wjo2 - 1.
               jo1 = jo1 + 1
            endif
            
! This is end of correction
            
            io2 = io1 + 1
            jo2 = jo1 + 1

            wio1 = 1. - wio2
            wjo1 = 1. - wjo2

            if (trim(iaction) == 'ndvi') then

! Interpolate from 4 surrounding values

               datp(iwl)  & 
                  = wio1 * (wjo1 * dato(io1,jo1) + wjo2 * dato(io1,jo2))  &
                  + wio2 * (wjo1 * dato(io2,jo1) + wjo2 * dato(io2,jo2))

            elseif (trim(iaction) == 'leaf_class' .or.  &
                    trim(iaction) == 'soil_text') then

! Use nearest data point - do not interpolate

               io = io2
               jo = jo2
               if (wio2 < .5) io = io1
               if (wjo2 < .5) jo = jo1
               
               idatp(iwl) = idato(io,jo)

            endif

         enddo

     endif
   enddo
enddo

deallocate(nump,numpind1,numpind2,idato,dato)

return
end subroutine leaf_database_read

End module
