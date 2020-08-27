Subroutine SQ_cov(Atot,Wtot,Rtot,Ntot,Npars,DeltaP_INPOP,cov,SOS)
implicit none
integer, intent(in) :: Ntot,Npars
real(kind=8), dimension(Ntot,npars) :: Atot,AA
real(kind=8), dimension(Ntot) :: Rtot, Wtot,U
real(kind=8), dimension(npars) :: DeltaP_INPOP,dnorm,S,stdev,DeltaP,Deltazp
real(kind=8), dimension(npars,npars) :: cov,cor,V,Rxinv,Rx
real(kind=8), dimension(Npars,Ntot) :: Htot
real(kind=8), dimension(Ntot) :: Res_fit
real(kind=8) :: SOS,sum,sigma,aux,Smax,Smin,CNRx,sumb,sumc
integer :: i,j,k,IER
real(kind=8), dimension(Npars*(Npars+1)/2) :: cova
character(len=20) :: string


print*,'entree ds SQ'

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!           !                                                                C
!  BLOCK 2  !                  SQUARE-ROOT FORMULATION                       C
!           !                                                                C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CONSTRUCTION OF Rx, THE LINEAR SYSTEM MATRIX OF THE PARAMETERS ESTIMATION  C
! FORMULA                                                                    C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! OVERWRITING Atot WITH Wtot^1/2 * Atot
! OVERWRITING Rtot WITH Wtot^1/2 * Rtot

!      open(19,file='Atot_sqr.bk2.dat',status='replace')
!      write(19,111)((Atot(i,j),j=1,Npars),i=1,Ntot)
!      close(19)  

!      open(18,file='Rtot_sqr.bk2.dat')
!      write(18,114) (Rtot(i),i=1,Ntot)
!      close(18) 

! Do not comment the writing of Atot

      do i=1,Ntot
         Rtot(i)=DSQRT(Wtot(i))*Rtot(i)
         do j=1,Npars
            Atot(i,j)=DSQRT(Wtot(i))*Atot(i,j)
         enddo
      enddo   
     
!      open (17,file='Atot_sqr.dat')
!          do i=1,Ntot
!          write(17,111)(AA(i,j),j=1,Npars)
!          enddo
!      close(17)

!     open(18,file='LS/Rtot_sqr.dat')
!      read(18,114) (Rtot(i),i=1,Ntot)
!      close(18)

! UPPER TRIANGULAR MATRIX Rx(Npars,Npars)

! PERFORM Npars ORTHOGONAL TRANSFORMATIONS TO UPPER TRIANGULARIZE Atot
! THE HOUSEOLDER FACTOR U IS COMPUTED FOR EACH K-th COLUMN OF A AND THEN
! THE ORTHOGONAL TRANSFORMATION IS APPLIED.

print*,Npars,Ntot

          AA=Atot

!          do i=1,Ntot
!          write(258,111)(AA(i,j),j=1,Npars)
!          enddo

      do k=1,Npars
         do i=1,Ntot
            U(i)=0d0
         enddo
         sum=0.0D0
         do i=k,Ntot
            sum=sum+AA(i,k)**2
         enddo
         if(AA(k,k).ge.0d0)then
            sigma=SQRT(sum)
         else
            sigma=-SQRT(sum)
         endif
         aux=SQRT(1.0D0+AA(k,k)/sigma)
         sumb=1.0D0/(sigma*aux)
         U(k)=aux
         do i=k+1,Ntot
            U(i)=sumb*AA(i,k)
         enddo
         do j=k,Npars
            sumc=0.0
            do i=k,Ntot
               sumc=sumc+U(i)*AA(i,j)
            enddo
            do i=k,Ntot
               AA(i,j)=AA(i,j)-sumc*U(i)
            enddo
         enddo
      enddo

! Rx IS THE FIRST Npars ROWS AND COLUMNS OF THE TRIANGULARIZED Atot.

      do j=1,Npars
         do i=1,Npars
            Rx(i,j)=AA(i,j)
         enddo
      enddo

! FOR FURTHER COMPUTATIONS WE NEED Atot TO CONTAIN ITS ORIGINAL VALUES.

!      open (17,file='Atot_sqr.dat',status='old')
!          do i=1,Ntot
!          read(17,*)(Atot(i,j),j=1,Npars)
!          enddo
!      close(17)

  
!      open (20,file='LS/Rx.dat')
!      read(20,*) ((Rx(i,j),i=1,Npars),j=1,Npars)
!      close(20)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!           !                                                                C
!  BLOCK 3  !            INVERSION OF THE SYSTEM MATRIX                      C
!           !                                                                C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! INVERSION OF THE SYSTEM MATRIX Rx IS PERFORMED BY ITS SVD DECOMPOSITION    C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! THE INPUT MATRIX Rx IS NORMALIZED BY THE NORM OF ITS COLUMNS IN ORDER TO 
! DECREASE ITS CONDITION NUMBER.

!      call norm(Rx,Npars,dnorm)

      do j=1,Npars
         dnorm(j)=0d0
         do i=1,Npars
            dnorm(j)=dnorm(j)+Rx(i,j)**2
         enddo
         dnorm(j)=sqrt(1d0/dnorm(j))
      enddo

      do j=1,Npars
         do i=1,Npars
            Rx(i,j)=Rx(i,j)*dnorm(j)
         enddo
      enddo


! SVDCMP TRASFORMS Rx ACCORDING TO Rx=U*diag[S(i)]*Vt
! NOTE THAT U IS OVERWRITTEN ON THE INPUT MATRIX Rx.

      call SVDCMP(Rx,Npars,Npars,Npars,Npars,S,V)

! CONDITION NUMBER OF Rx

      Smax=1d-12
      Smin=1d12
      do i=1,Npars
         if(S(i).gt.Smax) Smax=S(i)
         if(S(i).lt.Smin) Smin=S(i)
      enddo

      CNRx = Smax/Smin
!      write(6,'('' largest singular value = '',e9.3)') Smax
!      write(6,'('' smallest singular value = '',e9.3)') Smin
!      write(6,'('' condition number of Rx = '',d8.3)')CNRx

! INVERSION OF Rx
! Rxinv = V*diag[1/S(i)]*U
! Since from svdcmp subroutine Rx = U and I rename V = V*diag[1/S(i)], then
! Rxinv = V*Rx^t
! SETTING THE LIMIT FOR THE SINGULAR VALUE OF Rx.
! IF THERE IS ANY S(i) EXCEEDING THE LIMIT, ITS VALUE IS SET TO ZERO
! AND RXINV IS THE PSEUDO-INVERSE OF Rx.

      Smin=Smax*1d-18
      do i=1,Npars
         if(S(i).ge.Smin) then
         S(i)=1/S(i)
         else
         S(i)=0d0
         endif
         do j=1,Npars
            V(j,i)=V(j,i)*S(i)
        enddo
     !write(259,*)i,S(i),(V(j,i),j=1,npars)
      enddo

      CALL VMULFP(V,Rx,Npars,Npars,Npars,Npars,Npars,Rxinv,Npars,IER)
      write(6,*) 'ier=',ier

! REMOVING THE NORMALIZATION FACTORS.

       do i=1,Npars
          do j=1,Npars
             Rxinv(i,j)=Rxinv(i,j)*dnorm(i)
          enddo
 !      write(257,*)(Rxinv(i,j),j=1,npars)
       enddo

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!           !                                                                C
!  BLOCK 4  !            SOLUTION OF THE LINEAR SYSTEM                       C
!           !                                                                C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ESTIMATE OF THE PARAMETERS DELTA SOLUTION VECTOR DeltaP, AND CONSTRUCTION  C
! OF THE COVARIANCE MATRIX.                                                  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! Calculation of Rxinv^t. Rxinv_output = (Rxinv_input)^t

      call VTRAN(Rxinv,Npars,Npars)
      write(6,*) 'ier=',ier

! Calculation of cov = Rxinv*Rxinv^t 
! Vector cova is the upper-packed storage mode of the covariance matrix cov.

      call VTPROF(Rxinv,Npars,Npars,Npars,cova)
      write(6,*) 'ier=',ier

! Calculation of the covariance matrix.

      do i=1,Npars
         do j=i,Npars
            cov(i,j)=cova(i+j*(j-1)/2)
            cov(j,i)=cov(i,j)
         enddo
      enddo


string='(   (e12.6,3x))'
write(string(2:4),'(I3)')Npars

!print*,"correlation matrix",Npars,string
!      open(21,file="matrix.correlation_I1")
!      write(21,string)((cor(i,j),j=1,Npars),i=1,Npars)
!     do i=1,56
!      write(21,116)(cor(i,j),j=1,56)
!     enddo
!      close(21)

!print*,"covariance matrix",Npars,string
      open(21,file="matrix.covariance_I1_full")
      write(21,string)((cov(i,j),j=1,Npars),i=1,Npars)
      close(21)

! Calculation of the correlation matrix cor.

      do i=1,Npars
         do j=1,Npars
            cor(i,j)=cov(i,j)/dsqrt(abs(cov(i,i)*cov(j,j)))
         enddo
      enddo

!string='(   (e12.6,3x))'
!write(string(2:4),'(I3)')Npars

!print*,"correlation matrix",Npars,string
!      open(21,file="matrix.correlation_I1")
!      write(21,string)((cor(i,j),j=1,Npars),i=1,Npars)
!     do i=1,100
!      write(21,116)(cor(i,j),j=1,100)
!     enddo
!      close(21)

!===============================================================
! COMPUTING THE FIRST Npars ROWS OF THE HOUSEHOLDER MATRIX Htot.

! Htot = Rxinv_output*Atot^t = (Rxinv_input)^t*Atot^t = 
! =(Rxinv_input)^t*Wtot^(-1/2)*Atot^t

call VMULFP(Rxinv,Atot,Npars,Npars,Ntot,Npars,Ntot,Htot,Npars,IER)
      write(6,*) 'ier=',ier

!do i=1,ntot
!   write(255,*)(Htot(j,i),j=1,npars)
!enddo

! COMPUTING THE DELTA SOLUTION.
! Deltazp = Htot*Wtot^(1/2)*Rtot = Htot*Rtot

CALL VMULFF(Htot,Rtot,Npars,Ntot,1,Npars,Ntot,Deltazp,Npars,IER)
      write(6,*) 'ier=',ier

!do i=1,npars
!write(257,*)Deltazp(i)
!enddo

! DeltaP = (Rxinv_output)^t*Deltazp = Rxinv_input*Deltazp

 call VMULFM(Rxinv,Deltazp,Npars,Npars,1,Npars,Npars,DeltaP,Npars,ier)
      write(6,*) 'ier=',ier

!      open (21,file='DeltaP_I1.dat')
      do i=1,Npars
!         write(21,'(e26.20)') DeltaP(i)
      enddo
!      close(21)
       
      do i=1,Npars
         stdev(i)=sqrt(cova(i+(i*(i-1)/2)))
      enddo

! Artificial inversion of the sign of DeltaP before the INPOP implementation.

!      open(421,file="correction_sigma.I1")
      do i=1,Npars
	 DeltaP_INPOP(i)=DeltaP(i)
!!	 DeltaP_INPOP(i)=DeltaP(i)
!write(421,'(i3,1x,2(d27.20,1x),f15.10)')i,DeltaP_INPOP(i),stdev(i),stdev(i)*100/DeltaP_INPOP(i)
      enddo

! OVERWRITING Wtot^1/2*Atot WITH Atot.
! OVERWRITING Wtot^1/2*Rtot WITH Rtot.

     do i=1,Ntot
         Rtot(i)=1/DSQRT(Wtot(i))*Rtot(i)
         do j=1,Npars
            Atot(i,j)=1/DSQRT(Wtot(i))*Atot(i,j)
         enddo
      enddo   
      do i=1,Ntot-Npars
         Res_fit(i)=Rtot(i)
         do j=1,Npars
            Res_fit(i)=Res_fit(i)-Atot(i,j)*DeltaP(j)
         enddo
      enddo
 
! Calculation of the postfit SOS (Sum of Squares of the weighted residuals).

      SOS=0d0
      do i=1,Ntot-Npars
         SOS=SOS+Wtot(i)*Res_fit(i)**2
      enddo
      write(6,'('' Normalized postfit SOS ='',f6.1)') SOS/(Ntot-Npars)

110   format(5x,f15.10,3x,72(d27.20,4x))
111   format(29(d27.20,4x))
112   format(f4.2)
113   format(5x,f15.10,4x,i1,3x,i3,3x,i3,4x,e27.20,2x,i5)
114   format(e27.20)
116   format((111(e12.6,3x)))  
117   format(5x,f15.10,4x,e27.20)
     
      return
end  Subroutine

