!=======================================================================
	MODULE SEIJI_KLUT_SAV
	parameter ( nv = 3  ,mxn =17)
! nv = # of Variables in look up table.
! mxn = Max # of nodes for any single variable
! nva = # of nodes for each variables
! anodes = The values at each node in the table; for all variables

	logical ichoose(nv,2**nv)
	integer,save:: nva(nv)
	real,save:: anodes(mxn,nv)
	CONTAINS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------------------------------------------------------------------
	subroutine tf_table(n1,n2,ib)
	implicit none
	integer n1,n2,i,ii
	logical ib(n1,n2)

	do i=1,n2
	do ii=1,n1
	ib(ii,i)=btest(i-1,ii-1)
	enddo
!	print*,i,ib(1:n1,i)
	enddo

	return
	end subroutine tf_table

!-----------------------------------------------------------------------
	integer function ilui(iva)
	implicit none
	integer i,k
	integer  iva(nv)
	integer  mva(nv)

	mva(1)=nva(1)
	do i=2,nv
	 mva(i)=mva(i-1)*nva(i)
	enddo

	k= iva(1)
	do i=2,nv
	k= k+ (iva(i)-1)*mva(i-1)
	enddo

	ilui=k
	return
	end function ilui

!----------------------------------------------------------------
	real function getsk(val,table)
	implicit none
	logical,save:: notcalled
	data notcalled /.true./
	real table(13*11*17) ! For a small table 
	real value_table, swx,wx,rmin,rmax,rlui
	integer ,save:: i,j,k,iv,in
	integer,save::  iva(nv,2**nv)
	logical,save :: reinit
	integer,save :: ksav(2**nv)
	real ,save ::   wxsav(2**nv) 
	real ,save ::   rvasav(nv) !!RAINER

	real  wb(nv),wt(nv),rva(nv)
	real  val(nv)
	integer  ibase(nv),itop(nv)
! VAL(3) Pressure (Hpa)
! VAL(2) Temperature(k)	
! VAL(1) Natural log Number H2O concentration

!---------------------------------------------------------Begin init
!----------------------------------------------------------------BEGIN CALL ONCE
if ( notcalled ) then
	nva =(/13,11,17/)

!PRESSURE
	anodes(1:17,3) = (/1050.0, 834.0, 662.0, 526.0, 418.0, 332.0, &
                            264.0, 168.0, 106.0,  66.9,  42.2, &
                             26.6,  16.8,  10.6,   1.0,   0.1,  0.01/)
	anodes(1:17,3)=log(anodes(17:1:-1,3)) !!!! LOG PRESSURE & FLIP SMALL to LARGE (see mksktbl.f too)


!TEMPERATURE
	anodes(1:11,2) = (/120.0, 140.0, 160.0, 180.0, 200.0, 220.0, &
                         240.0, 260.0, 280.0, 300.0, 320.0/)

!H2O CONCENTRATION
	anodes(1:13,1) = (/0.15e26, 0.47e25, &
                         0.15e25, 0.47e24, 0.15e24, 0.47e23, &
                         0.15e23, 0.47e22, 0.15e22, 0.47e21, &
                         0.15e21, 0.47e20, 0.15e20/)
! FLIP SMALL to LARGE and H20Conc to log(H20MMR) (see mksktbl.f too)
	anodes(1:13,1)=log(anodes(13:1:-1,1))


	call tf_table(nv,2**nv,ichoose)
	notcalled = .false.
	
	endif
!-------------------------------------------------------------------END CALL ONCE
!---------------------------------------------------------End init
	call limit(val)
	val(3)=log(val(3)) !!! LOG PRESSURE
	reinit= .false. !! ASSUME OLD VAL
!!!!		reinit= .true. !! reinit each entry
!--------------------
	 VSET : do iv=1,nv 
	rmin = anodes(1,iv)
	rmax = anodes(nva(iv),iv)
	if (val(iv) <= rmin )then
	  val(iv) = rmin
	  rva(iv) = 1
	elseif (val(iv) >= rmax ) then
	  val(iv) = rmax
	  rva(iv) = nva(iv)
	else

	FNODES : do in = 1,nva(iv)-1
	if( val(iv) >= anodes(in,iv) .and. val(iv) < anodes(in+1,iv) ) then
	rva(iv) = in+( val(iv)-anodes(in,iv)) / (anodes(in+1,iv)-anodes(in,iv))
	endif
	enddo FNODES
	endif

	if( rva(iv) .ne. rvasav(iv) ) reinit= .true.
!	print*,iv,rva(iv),reinit
	enddo VSET
!---------------------

	 ibase = floor(rva)
	 where( ibase < 1) ibase=1

	 itop= ceiling(rva)
	 where( itop > nva) itop=nva

	wb=1.0-(rva-ibase)
	wt=1.0-wb
	
!	print'(A10,10f5.1)','Exact',rva(1:nv)
!	print'(A10,10I5)','BASE',ibase
!	print'(A10,10I5)','TOP',itop
!	print'(A10,10f5.1)','Wb',wb	
!	print'(A10,10f5.1)','Wt',wt

	swx=0
	rlui =0
!	print*,'---'

	NODE : do j=1,2**nv
	wx=1.0


	if (reinit) then !!!!!!
	rvasav(1:nv) =rva(1:nv) !!RAINER
	 DIMENSION: do i=1,nv

	 if( ichoose(i,j) ) then
	  iva(i,j) = itop(i)
	  wx=wx*wt(i)
	 else
	  iva(i,j) = ibase(i)
	  wx=wx*wb(i)
	 endif

	 enddo DIMENSION
	  k = ilui(iva(1,j))
	    ksav(j) = k
	   wxsav(j) = wx
!!	  rvasav(j) =rva(j) !!RAINER
	endif    !!!!!!!!!!!!!!!!!!!

	wx = wxsav(j)
	k  = ksav(j)

	 if( wx > 1.0E-10 ) then

!	 k = ilui(iva(1,j))

	 value_table = table(k)

	 swx       = swx  +               wx
	 rlui      = rlui + value_table * wx

!	print'(2I4,2f10.5,ES10.2,10I4)',j,k,wx,swx,rlui,iva(1:nv,j)
	if (swx >= .9999)exit
	endif

	enddo NODE
	getsk = rlui/swx
	return
	end function getsk
!--------------------------------------------------------------
!--------------------------------------------------------------
subroutine limit(val)
implicit none
real val(3),tlim,pbar,tbar

pbar=val(3)
tbar=val(2)
tlim=val(2)
if(pbar <=0.01 		      .and.tbar>=280) 	tlim=280

if(pbar >=0.10 .and.pbar<=1   .and.tbar <180) 	tlim=180
if(pbar  >1.00 .and.pbar<=332 .and.tbar <160) 	tlim=160

if(pbar >=11.0 .and.pbar<=168 .and.tbar >280) 	tlim=280
if(pbar >168 .and.pbar <= 270 .and.tbar>=300) 	tlim=300

if(pbar >= 330 .and. pbar < 662 .and.tbar<=180) tlim=180
if(pbar >= 662 		.and.tbar<=200) tlim=200

val(2)= tlim
return
end subroutine limit

	END MODULE SEIJI_KLUT_SAV


!!==============================================================
!! INTERFACE WITH FULIOU CODE

	subroutine seijifu_ht02a_sav(ib0,ig0,skhk1)
	
	USE FUINPUT
	use SEIJI_KLUT_SAV,ONLY :getsk 
	USE SKTBL_HT02a

	implicit none
	real skhk1

	real hk1,fk1o3,sol_spect,fk1h2o
	real tg

	common /band1/ hk1(10), fk1o3(10),sol_spect(0:7),fk1h2o(10)
!	common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
	common /gas/ tg(nvx)

	real,save:: val(3) ,tglevel(nv1x,10,6),skhksfr(10,6) 
	
	integer i,ib,ig,ib0,ig0
	real qbar,density,conc_to_mmr,fq
	integer kg0(6),kg1(6)
	logical,save :: reinitk
	
		if (ib0 ==1 .and. ig0 ==9 ) reinitk =.true.
	if ( reinitk ) then

	kg0(1:6) = (/ 9,1,1,1,1,1/)
	kg1(1:6) = (/10,7,8,7,8,7/) !ht02
!--
	LEVEL: do i=1,nv1

	BANDINIT : do ib=1,6
	KINIT : do ig= kg0(ib),kg1(ib)

	val(3)= pp(i)
	val(2)= pt(i)
	qbar  = ph(i)
	if (qbar < 1E-08) qbar=1E-08
	density= 1.275*(VAL(3)/1000.0)*(273.15/VAL(2))
	conc_to_mmr = (1.0/6.02214199E+26)*18./density

	val(1)= log(qbar/conc_to_mmr) !!! MMR to log(#Conc)

	  if (ib == 1 .and. ig == 9 ) then

	    tglevel(i,ig,ib) = getsk(val,sktbl00(1,1,1, 1) )
	   skhksfr(ig,ib)=skhks00(1)
	elseif (ib == 1 .and. ig == 10 ) then

	   tglevel(i,ig,ib) = getsk(val,sktbl01(1,1,1, 1) )
	   skhksfr(ig,ib)=skhks01(1)
!-------
	elseif ( ib == 2  ) then

	 tglevel(i,ig,ib) = getsk(val,sktbl02(1,1,1,ig) )
	 skhksfr(ig,ib)=skhks02(ig)
	elseif ( ib == 3 ) then 

	 tglevel(i,ig,ib) = getsk(val,sktbl03(1,1,1,ig) )
	 skhksfr(ig,ib)=skhks03(ig)
	elseif ( ib == 4 ) then 

	 tglevel(i,ig,ib) = getsk(val,sktbl04(1,1,1,ig) )
	 skhksfr(ig,ib)=skhks04(ig)
	elseif ( ib == 5 ) then

	 tglevel(i,ig,ib) = getsk(val,sktbl05(1,1,1,ig) )
	 skhksfr(ig,ib)=skhks05(ig)
	elseif ( ib == 6  ) then

	 tglevel(i,ig,ib) = getsk(val,sktbl06(1,1,1,ig) )
	 skhksfr(ig,ib)=skhks06(ig)
	else
!	 print*,ib,ig
!	 print*,'BAD BAND # seiji_k.f90'
!	 stop
	endif
	
!	print'(I4,3f9.3,f8.4,E12.2)',ig,val(2:4),skhksfr(ig,ib),tg(i)
	enddo KINIT
	enddo BANDINIT

	enddo LEVEL
	reinitk=.false.
	endif ! reinitk
!----------------------------------------------------------------------
!----------------------------------------------------------------------



	LAYER : do i=1,nv
	tg(i)=(tglevel(i,ig0,ib0)+tglevel(i+1,ig0,ib0))*0.5 * dz(i)*1000. ! TAU PROFILE
	enddo LAYER


	skhk1 = skhksfr(ig0,ib0)*sol_spect(ib0)
	if ( ib0 ==5 .and. lband6a )then
	 skhk1 = skhksfr(ig0,ib0)*(sol_spect(5)+sol_spect(7))
	endif

!ADD FU O3 K
	if ( ib0 == 1 .and. (ig0 == 9 .or. ig0 == 10) ) then
	 skhk1 = skhksfr(ig0,ib0)*(sol_spect(1)* hk1(ig0) )
	 fq = fk1o3(ig0)*238.08
	  do i=1,nv
	  tg(i) = tg(i)+ (po(i)+po(i+1)) * (pp(i+1)-pp(i))*fq
	  enddo
	
	endif

	end subroutine seijifu_ht02a_sav
!!==============================================================
