MODULE ICEDIRSFC
USE FUINPUT , only: mbsx,mccx,fi
USE FUOUTPUT ,only: foscc
implicit none
real , dimension(mbsx,mccx) :: tau_uc,tau_co


 CONTAINS 
!============================================================
   subroutine dircorrect 
     integer ic,icc,ib,isc
     real rt,taco,tauuc
     real diftau(15) 
     
     if (fi%u0 <= 0.0 ) return
     
     
     WAER: do isc  = 1,2
        
        foscc(0,isc)%dirsfc(1:15) =  foscc(0,isc)%rswdir(fi%nv+1,1:15)
        foscc(0,isc)%difsfc(1:15) =  foscc(0,isc)%rswdif(fi%nv+1,1:15)
        
        CLDCON : do icc =1 ,mccx
           
           if ( fi%fc(icc)%cldfrac > 0 ) then
              
              diftau(1:10)  = tau_co(1,icc)  -tau_uc(1,icc)
              diftau(11:15) = tau_co(2:6,icc)-tau_uc(2:6,icc)
              
              BAND: do ib = 1,15
                 
                 if ( foscc(icc,isc)%rswdir(1,ib) > 0 )then 
                    rt =  foscc(icc,isc)%rswdir(fi%nv+1      ,ib)/&
                         foscc(icc,isc)%rswdir(1,ib) ! Transmission Cloud,Gas,Rayleigh
                 else
                    rt = 0.0
                    !   print*,foscc(icc,isc)%rswdir(1,ib)
                    !   stop
                 endif
                 
                 !   taco = 128.
                 !   if ( rt > 1.0E-36 ) taco = -fi%u0 * log(rt) ! Direct Transmission to Tau
                 taco = -fi%u0 * log( max(rt,1.0e-36)) ! Direct Transmission to Tau
                 
                 tauuc = taco - diftau(ib) ! subtract of difference between tau adjusted for forward scatter peak.
                 
                 foscc(icc,isc)%dirsfc(ib) =  foscc(icc,isc)%rswdir(1,ib) * exp(- tauuc/fi%u0)
                 foscc(icc,isc)%difsfc(ib) =  foscc(icc,isc)%rswfd(fi%nv+1,ib) -  foscc(icc,isc)%dirsfc(ib)
              enddo BAND
              
           endif ! CLDFRAC >0
           
        enddo CLDCON


     enddo WAER


   end subroutine dircorrect
END MODULE ICEDIRSFC
