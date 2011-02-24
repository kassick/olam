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
integer function lastchar(str)
implicit none
character(len=*) :: str
integer :: n,ln
! returns last non-blank character position from a string

ln=len(str)
do n=ln,1,-1
   if(str(n:n).ne.' ') then
      lastchar=n
      return
   endif
enddo
lastchar=0

return
end

!===============================================================================

subroutine deblank(str1,str2,nch)
implicit none
character(len=*) :: str1,str2
integer :: n,ln,nch

! strips blanks from a string and returns number of chars

str2=' '
ln=len(str1)
nch=0
do n=1,ln
   if(str1(n:n).ne.' ') then
      nch=nch+1
      str2(nch:nch)=str1(n:n)
   endif
enddo

return
end

!===============================================================================

subroutine char_strip_var(line,var,line2)
implicit none
character(len=*) :: line,var,line2
integer :: nn,ncl,nb

! removes instances of a substring from a string

ncl=len(line)
do nn=1,ncl
   if(line(nn:nn).ne.' ') then
      nb=index(line(nn:),' ')
      var=line(nn:nn+nb-1)
      goto 25
   endif
enddo
25 continue
line2=line(nn+nb-1:)

return
end

!===============================================================================

subroutine tokenize(str1,tokens,ntok,toksep,nsep)
implicit none
integer :: nsep,ntok
character(len=*) :: str1,tokens(*)
character(len=1) :: toksep(nsep)

character(len=256) :: str
integer :: npt,nch,nc,ns

! this routine "parses" character string str into different pieces
! or tokens by looking for  possible token separators (toks
! str contains nch characters.  the number of tokens identified is nto
! the character string tokens are stored in tokens.

ntok=0
npt=1
call deblank(str1,str,nch)
do nc=1,nch
   do ns=1,nsep
      if(str(nc:nc).eq.toksep(ns).or.nc.eq.nch) then
         if(nc-npt.ge.1)then
            ntok=ntok+1
            tokens(ntok)=str(npt:nc-1)
            if(nc.eq.nch.and.str(nc:nc).ne.toksep(ns)) then
               tokens(ntok)=str(npt:nc)
               goto 10
            endif
         endif
         ntok=ntok+1
         tokens(ntok)=str(nc:nc)
         npt=nc+1
         goto 10
      endif
   enddo
10      continue
enddo
return
end

!===============================================================================

subroutine tokenize1(str1,tokens,ntok,toksep)
implicit none
integer :: ntok
character(len=*) :: str1,tokens(*)
character(len=1) :: toksep

character(len=256) :: str
integer :: nch,ist,npt,nc

! this routine "parses" character string str into different pieces
! or tokens by looking for  possible token separators (toks
! str contains nch characters.  the number of tokens identified is nto
! the character string tokens are stored in tokens.

call deblank(str1,str,nch)
!print*,'ttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt'

ist=1
if(str(1:1).eq.toksep) ist=2
npt=ist
ntok=0
do nc=ist,nch
         !print*,'ttttttttttt0:',ist,nch,nc,npt
         !print*,'ttttttttttt0:',str(nc:nc),toksep
   if(str(nc:nc) == toksep .or. nc == nch) then
      if(nc-npt >= 1) then
         ntok=ntok+1
         tokens(ntok)=str(npt:nc-1)
         !print*,'ttttttttttt1:',ntok,npt,nc,tokens(ntok)
         if(nc == nch .and. str(nc:nc) /= toksep) then
            tokens(ntok)=str(npt:nc)
            !print*,'ttttttttttt2:',ntok,npt,nc,tokens(ntok)
            exit
         endif
         npt=nc+1
      elseif(nc == nch) then
         ntok=ntok+1
         tokens(ntok)=str(npt:nc)
         !print*,'ttttttttttt3:',ntok,npt,nc,tokens(ntok)
         exit
      endif
   endif
enddo

return
end

!===============================================================================

integer function letter (str)
implicit none
character*(*) str

! First character alpha check - test to see if the first character of
! the string STR is alphabetic: LETTER = 0 if 'no', = 1 if 'yes'.

letter=0
if((str(1:1).ge.'A'.and.str(1:1).le.'Z').or.  &
   (str(1:1).ge.'a'.and.str(1:1).le.'z')) letter=1
   
return
end

!===============================================================================

integer function number (str)
implicit none
character(len=*) :: str

! First character number check - test to see if the first character of
! the string STR is numeric:  NUMBER = 0 if 'no', = 1 if 'yes' (includ
! a decimal point or minus sign).

number=0
if(str(1:1).ge.'0'.and.str(1:1).le.'9') number=1
if(str(1:1).eq.'.'.or.str(1:1).eq.'-') number=1

return
end

!===============================================================================

integer function letint (str)
implicit none
character(len=*) :: str

! First character integer variable check - test to see if the first
! character of STR is an I, J, K, L, M, or N, or i, j, k, l ,m or n:
! LETINT = 0 if 'no', = 1 if 'yes'

letint=0
if((str(1:1).ge.'I'.and.str(1:1).le.'N').or.  &
   (str(1:1).ge.'i'.and.str(1:1).le.'n')) letint=1
   
return
end

!===============================================================================

integer function letquo (str)
implicit none
character(len=*) :: str

! First character quote check - test to see if the first character
! of STR is a quote:  LETQUO = 0 if 'no', = 1 if 'yes'.

letquo=0
if(str(1:1).eq.'''') letquo=1

return
end
