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
subroutine topo_database_read()

use misc_coms,   only: io6, topo_database
use mem_grid,    only: mma, xem, yem, zem, topm, glatm, glonm
use mem_ijtabs,  only: jtab_m, itab_m
use consts_coms, only: erad, piu180, pi1, pi2
use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec

implicit none

integer, parameter :: nradp = 2, np = 3*nradp*(nradp-1)+1 ! np = # data values avgd per topm pt

integer, allocatable :: ptable(:)
integer, allocatable :: nump(:,:)
integer, allocatable :: numpind1(:,:)
integer, allocatable :: numpind2(:,:)

real, allocatable :: dato(:,:)
real, allocatable :: datp(:,:)
real, allocatable :: pmlat(:,:)
real, allocatable :: pmlon(:,:)

integer :: im,ip,jp,nio,njo,niosh,njosh,nperdeg  &
   ,isoc,iwoc  &
   ,isocpt,isocpo,iwocph,iwocpt,iwocpo,lb  &
   ,ind,j2d,j1d,ind1,ind2,io1,jo1,io2,jo2  &
   ,ifiles,jfiles,ifile,jfile,ptab,ptab0,j,io_full,jo_full  &
   ,iw,iradp,nazimp,iazimp

integer :: ndims,idims(2)

real :: offlat,offlon  &
       ,pmlat1,pmlon1,dlato,dlono,wio2,wjo2,wio1,wjo1,yp,xp,c1,c2  &
       ,rio_full,rjo_full,radm,radp,azimpoff,azimp,offpix

character(len=3) :: title1
character(len=4) :: title2
character(len=80) :: title3,hfn,ofn,title

logical l1,l2

write(io6,*) 'Starting topography database read'

! Check topography file prefix
   
lb = len_trim(topo_database)
if (lb <= 0) then
   write(io6,*) '==================================================='
   write(io6,*) '|  Problem in topo_database: Input data prefix incorrect !'
   write(io6,*) '|  File prefix:',topo_database
   write(io6,*) '==================================================='
   stop 'topo_database_file'
endif

! Open, read, and close topography dataset header file

call rams_f_open(29,trim(topo_database)//'HEADER','FORMATTED','OLD','READ',0)
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

ifiles = 360 * nperdeg / niosh
jfiles = 180 * nperdeg / njosh
   
! Allocate 8 arrays

allocate (nump    (ifiles,jfiles)) ! # pts filled by file(ifile,jfile)
allocate (numpind1(ifiles,jfiles))
allocate (numpind2(ifiles,jfiles))
allocate (ptable  (np*mma))

allocate (dato(nio,njo))
allocate (datp(np,mma))
allocate (pmlat(np,mma))
allocate (pmlon(np,mma))

! Fill pmlat and pmlon arrays with offset latitudes and longitudes of all p points

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_m(1)%jend(1); im = jtab_m(1)%im(j)
!----------------------------------------------------------------------
call qsub('M',im)

! Compute radius from current M point for averaging topo data values

   radm = sqrt(itab_m(im)%arm / pi1)

! Loop over array (array size options = 1,7,19,37) of pts within radius 
! RADM - get lat/lon for each pt
   
   ip = 0
   do iradp = 1,nradp
      radp = real(2 * iradp - 2) / real(2 * nradp - 1) * radm
      nazimp = 6 * (iradp - 1)
      azimpoff = .5 * mod(iradp,2)  ! offset to stagger azimp of alt. rows
      do iazimp = 1,nazimp
         azimp = (real(iazimp) + azimpoff) / real(nazimp) * pi2  
         xp = radp * cos(azimp)
         yp = radp * sin(azimp)
         ip = ip + 1
         call xy_ll(pmlat(ip,im),pmlon(ip,im),glatm(im),glonm(im),xp,yp)
      enddo
   enddo
enddo
call rsub('M-topodatabase_a',1)

do jfile = 1,jfiles
   do ifile = 1,ifiles
      nump(ifile,jfile) = 0
   enddo
enddo

! Get file index (ifile,jfile) within full dataset and count number of p 
! points (nump) that occur in each file

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_m(1)%jend(1); im = jtab_m(1)%im(j)
!----------------------------------------------------------------------
call qsub('M',im)
   do ip = 1,np

      pmlat1 = max(-89.9999,min(89.9999,pmlat(ip,im)))
      pmlon1 = pmlon(ip,im)
            
      if (pmlon1 >=  180.) pmlon1 = pmlon1 - 360.
      if (pmlon1 <  -180.) pmlon1 = pmlon1 + 360.

      rio_full = (pmlon1 + 180.) * nperdeg ! must ignore pixel offset here
      rjo_full = (pmlat1 +  90.) * nperdeg ! must ignore pixel offset here

      io_full = int(rio_full)
      jo_full = int(rjo_full)

      ifile = io_full / niosh + 1
      jfile = jo_full / njosh + 1

      nump(ifile,jfile) = nump(ifile,jfile) + 1  ! Summation to # pts 
                                     ! filled by file (ifile,jfile) in dataset

   enddo
enddo
call rsub('M-topodatabase_b',1)

! Set up array index values for ptable array

ind = 1
do jfile = 1,jfiles
   do ifile = 1,ifiles
      numpind1(ifile,jfile) = ind   ! Beginning pt index for file (ifile,jfile)
      numpind2(ifile,jfile) = ind   ! For now, same as numpind1
      ind = ind + nump(ifile,jfile)
   enddo
enddo

! Fill ptable array

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_m(1)%jend(1); im = jtab_m(1)%im(j)
!----------------------------------------------------------------------
call qsub('M',im)

   j1d = (im - 1) * np
   do ip = 1,np

      pmlat1 = max(-89.9999,min(89.9999,pmlat(ip,im)))
      pmlon1 = pmlon(ip,im)

      if (pmlon1 >=  180.) pmlon1 = pmlon1 - 360.
      if (pmlon1 <= -180.) pmlon1 = pmlon1 + 360.

      rio_full = (pmlon1 + 180.) * nperdeg ! must ignore pixel offset here
      rjo_full = (pmlat1 +  90.) * nperdeg ! must ignore pixel offset here

      io_full = int(rio_full)
      jo_full = int(rjo_full)

      ifile = io_full / niosh + 1
      jfile = jo_full / njosh + 1

      ind = numpind2(ifile,jfile)

      ptable(ind) = j1d + ip  ! point index in 1-d
      numpind2(ifile,jfile) = numpind2(ifile,jfile) + 1 ! Sums to ending index
            
   enddo

enddo
call rsub('M-topodatabase_c',1)

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

         title3 = trim(topo_database)//title1//title2//'.h5'
         
         inquire(file=trim(title3),exist=l1,opened=l2)

! Read file

         if (l1) then
            write(io6,*) 'topo_database_read2 ',trim(title3)

            call shdf5_open(title3,'R')
            ndims=2 ; idims(1)=nio ; idims(2)=njo
            call shdf5_irec(ndims,idims,'topo',rvara=dato)
            call shdf5_close()
         else
            write(io6,*) 'Topography file is missing'
            write(io6,*) 'Topography filename = ',title3(1:lb)
            write(io6,*) 'Stopping model run'
            stop 'missing topo file'
         endif
         
         do ind = ind1,ind2-1

            ptab = ptable(ind)         
            ptab0 = ptab - 1

            im = ptab0 / np + 1
            j1d = (im - 1) * np

            ip = ptab - j1d

            pmlat1 = max(-89.9999,min(89.9999,pmlat(ip,im)))
            pmlon1 = pmlon(ip,im)
            
            if (pmlon1 >=  180.) pmlon1 = pmlon1 - 360.
            if (pmlon1 <= -180.) pmlon1 = pmlon1 + 360.

            rio_full = (pmlon1 + 180.) * nperdeg ! must ignore pixel offset here
            rjo_full = (pmlat1 +  90.) * nperdeg ! must ignore pixel offset here

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

            datp(ip,im)                                                 &
               = wio1 * (wjo1 * dato(io1,jo1) +  wjo2 * dato(io1,jo2))  &
               + wio2 * (wjo1 * dato(io2,jo1) +  wjo2 * dato(io2,jo2))
               
         enddo
         
      endif
   enddo
enddo

! Average datp points together to get each topm value

c1 = 1. / float(np)

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_m(1)%jend(1); im = jtab_m(1)%im(j)
!----------------------------------------------------------------------
call qsub('M',im)

   topm(im) = 0.
   do ip = 1,np
      topm(im) = topm(im) + datp(ip,im)
   enddo
   topm(im) = topm(im) * c1

! Special for israel:  add 400 m
!!!   topm(im) = topm(im) + 400.
! end special


enddo
call rsub('M-topodatabase_d',1)

deallocate(nump,numpind1,numpind2,ptable,dato,datp,pmlat,pmlon)

return
end subroutine topo_database_read
