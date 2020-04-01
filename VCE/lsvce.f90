! ifort -O3 lsvce.f90 -o lsvce ; lsvce
implicit none
integer :: i,k,j,IA,IER
integer, parameter:: n1=300,m=12,n2=300,nskip=21996,nfull=40236
real(kind=8), dimension(n1,m)::A1
real(kind=8), dimension(n2,m)::A2
real(kind=8), dimension(n1+n2,m)::At
real(kind=8), dimension(n1) :: omc1
real(kind=8), dimension(n2) :: omc2
real(kind=8), dimension(m) :: AF
real(kind=8), dimension(n1+n2) :: omct,S,E,E1,E12,E21,E121,E122,tomct
real(kind=8), dimension(2*n1+2*n2) :: Wk
real(kind=8) :: dat,wpond,c1,c2,trace,omcf
real(kind=8), dimension(n1,n1) :: Q11
real(kind=8), dimension(n2,n1) :: Q10
real(kind=8), dimension(n2+n1,n1) :: Q111
real(kind=8), dimension(n2+n1,n2) :: Q01
real(kind=8), dimension(n2+n1,n2+n1) :: Q1,Q2,QY,IQY,aa,ab,mN22,mN21
real(kind=8), dimension(m,m) :: NN3,NN4
real(kind=8), dimension(n2+n1,m) :: NN2
real(kind=8), dimension(n2+n1,n2+n1) :: NN5,paortho,DNN5,matr,mq1,mq2
real(kind=8), dimension(n2+n1,n2+n1) :: mN11,mN12
real(kind=8), dimension(m,n2+n1) :: NN1,tAt,Nqy
real(kind=8), dimension(2,2) :: Nmat,INmat
real(kind=8), dimension(2) :: ll,newsig
real(kind=8), parameter :: ua=1.5e11

open(355,file="RA.out",status="old")
do i=1,nfull
read(355,*)dat,wpond,omcF,(AF(j),j=1,m)
  if(i.le.n1)then
        omc1(i)=omcF
        do j=1,m
        A1(i,j)=AF(j)
        enddo
   endif
   if(i>nskip.and.i.le.nskip+n2)then
        omc2(i-nskip)=omcF
        do j=1,m
        A2(i-nskip,j)=AF(j)
        enddo
   endif 
enddo
close(355)
omc1=omc1
omc2=omc2

print*,"======A=====",ua
do i=1,n2
!print*,(A2(i,1)),A2(i,m)
enddo

At=0d0
tAt=0d0
omct=0d0

do j=1,m
   do i=1,n1
   At(i,j)=A1(i,j)
   omct(i)=omc1(i)
   enddo
    do i=1,n2
        At(i+n1,j)=A2(i,j)
   omct(i+n1)=omc2(i)
    enddo
enddo    
   
! =======================     
! remplissages des Q ===============
!
Q1=0d0
do i=1,n1
 Q1(i,i)=1d0
enddo
!print*,'Q1====='
do i=1,n1+n2
!write(*,'(20(f10.5,1x))')(Q1(i,j),j=1,n1+n2)
        do j=1,n1+n2
!              if(Q1(i,j)==1)print*,i,j
        enddo
enddo

Q2=0d0
do i=n2+1,n1+n2
Q2(i,i)=1d0
enddo
!print*,'Q2====='
do i=1,n1+n2
!write(*,'(20(f10.5,1x))')(Q2(i,j),j=1,n1+n2)
do j=1,n1+n2
!if(Q2(i,j)==1)print*,i,j
enddo
enddo
! =======================    

QY=0d0
IQY=0d0
c1=1d0
c2=1d0
do i=1,n1+n2
   do j=1,n1+n2
        QY(i,j)=c1*Q1(i,j)+c2*Q2(i,j)
        IQY(i,j)=(1d0/c1)*Q1(i,j)+(1d0/c2)*Q2(i,j)
   enddo
enddo


tAT=TRANSPOSE(At)

NN1=MATMUL(tAt,IQY)
NN2=MATMUL(IQY,At)
NN3= MATMUL(tAt,NN2)
!print*,'NN3 ===='
!do i=1,m
!print*,(NN3(j,i),j=1,m)
!enddo

call inverse(NN3,NN4,m)
print*,'inversion 1'
Nqy=MATMUL(NN4,NN1)
print*,'NN4 ===='
!do i=1,m
!print*,(NN4(j,i),j=1,m)
!enddo


NN5=MATMUL(At,Nqy)
DNN5=0d0

do i=1,n1+n2
   DNN5(i,i)=1d0
enddo    
    
do i=1,n1+n2
  do j=1,n1+n2
        paortho(i,j)=DNN5(i,j)-NN5(i,j)
  enddo
enddo

E=MATMUL(paortho,omct)
do i=1,n1+n2
!print*,e(i),omct(i)
enddo
matr=MATMUL(IQY,paortho)
print*,'matr'


mq1=MATMUL(Q1,matr)
mq2=MATMUL(Q2,matr)
print*,'mQ'

aa=MATMUL(matr,mq1)

ab=MATMUL(matr,mq2)

print*,'mN'
mN11=MATMUL(Q1,aa)
mN12=MATMUL(Q1,ab)
mN22=MATMUL(Q2,ab)
mN21=MATMUL(Q2,aa)


print*,'Nmat-trace'
Nmat(1,1)=0.5d0*trace(mN11,n1+n2)
Nmat(1,2)=0.5d0*trace(mN12,n1+n2)
Nmat(2,2)=0.5d0*trace(mN22,n1+n2)
Nmat(2,1)=0.5d0*trace(mN21,n1+n2)
print*,'Nmat'
!print*,Nmat(1,1),Nmat(1,2),Nmat(2,1),Nmat(2,2)

!CALL LGINF(Nmat,2,2,2,0,INmat,2,S,WK,IER)
call inverse(Nmat,INmat,2)
print*,'inversion 2'
!print*,INmat(1,1),INmat(1,2),INmat(2,1),INmat(2,2)

E1=MATMUL(IQY,E)
E12=MATMUL(Q1,E1)
E21=MATMUL(Q2,E1)
E121=MATMUL(IQY,E12)
E122=MATMUL(IQY,E21)

ll=0d0
do i=1,n1+n2
ll(1)=ll(1)+e(i)*E121(i)
ll(2)=ll(2)+e(i)*E122(i)
enddo
ll=0.5*ll

newsig=MATMUL(INmat,ll)
newsig=sqrt(newsig)

print*,newsig*ua
print*,sqrt(sqrt(INmat(1,1))),sqrt(sqrt(INmat(2,2)))
End

! ===========================================

double precision function trace(N,m)
integer :: i,j,m
real(kind=8), dimension(m,m) :: N
real(kind=8) :: sum
trace=0d0
do i=1,m
trace=trace+N(i,i)
enddo
return 
end

subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

