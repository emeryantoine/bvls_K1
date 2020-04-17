! ifort -O3 lsvce.f90 -o lsvce ; lsvce
implicit none
integer :: i,j,k,ii,m
!integer, parameter:: n1=800,m=12,n2=600,nskip=21996,nfull=40236
integer, parameter:: n1=21996,n2=18240,nskip=21996,nfull=40236
integer, parameter:: niter=15
real(kind=8), parameter :: eps=1e-4
!real(kind=8), dimension(n1,m)::A1
!real(kind=8), dimension(n2,m)::A2
real(kind=8), allocatable, dimension(:,:) :: A1,A2,At
!real(kind=8), dimension(n1+n2,m)::At
!real(kind=8), dimension(n1) :: omc1
!real(kind=8), dimension(n2) :: omc2
real(kind=8), allocatable, dimension(:) :: omct,E,omc1,omc2,AF
real(kind=8) :: dat,wpond,c1,c2,trace,omcf
real(kind=8)::sum_mN11,mmN11,sum_mN12,mmN12,sum_mN21,mmN21,sum_mN22,mmN22
! real(kind=8), dimension(n2+n1,n2+n1) :: Q1,Q2,QY,IQY,aa,ab,mN22,mN21
!real(kind=8), dimension(n2+n1,n2+n1) :: NN5,paortho,DNN5,matr,mq1,mq2
!real(kind=8), dimension(n2+n1,n2+n1) :: mN11,mN12
real(kind=8), allocatable, dimension(:,:) :: Q1,Q2,QY,IQY
real(kind=8), allocatable, dimension(:,:) :: NN5,paortho
real(kind=8), allocatable, dimension(:,:) :: p1m,p2m,m11,m12

real(kind=8), allocatable, dimension(:,:) :: NN3,NN4,NN1,tAt,Nqy
!real(kind=8), dimension(m,n2+n1) :: NN1,tAt,Nqy

real(kind=8), dimension(2,2) :: Nmat,INmat
real(kind=8), dimension(2) :: ll,newsig,cc,outsol,outsolavant
real(kind=8), parameter :: ua=1.5e11

open(355,file="dimension",status="old")
read(355,*)m
close(355)

allocate(A1(n1,m))
allocate(A2(n2,m))
allocate(AF(m))
allocate(omc1(n1))
allocate(omc2(n2))


print*,"lecture",m

open(355,file="../../transfert/RA.out.mexmro.12",status="old")
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


allocate(omct(n1+n2))
allocate(At(n1+n2,m))

At=0d0
omct=0d0

print*,"construction de A et omc"

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
   
deallocate(A1)
deallocate(A2)
deallocate(omc1)
deallocate(omc2)

! =======================     
! remplissages des Q ===============
!

cc(1)=1.9
cc(2)=1.3

do ii=1,niter

print*,'ITER',ii,'====',(cc(j),j=1,2)

c1=cc(1)**2
c2=cc(2)**2

allocate(Q1(n1+n2,n1+n2))
allocate(Q2(n1+n2,n1+n2))
allocate(QY(n1+n2,n1+n2))
allocate(IQY(n1+n2,n1+n2))

Q1=0d0
do i=1,n1
 Q1(i,i)=1d0
enddo

Q2=0d0
do i=n1+1,n1+n2
Q2(i,i)=1d0
enddo

! =======================    

QY=0d0
IQY=0d0

print*,"contruction de QY et IQY"

   do j=1,n1+n2
do i=1,n1+n2
        QY(i,j)=c1*Q1(i,j)+c2*Q2(i,j)
        IQY(i,j)=(1d0/c1)*Q1(i,j)+(1d0/c2)*Q2(i,j)
   enddo
enddo
deallocate(Q1)
deallocate(Q2)
deallocate(QY)

allocate(NN3(m,m))
allocate(NN4(m,m))
allocate(NN1(m,n2+n1))
allocate(Nqy(m,n2+n1))
allocate(tAT(m,n2+n1))

tAt=0d0

print*,"contruction de NQY "

tAT=TRANSPOSE(At)
NN1=MATMUL(tAt,IQY)
NN3=MATMUL(NN1,At)
call inverse(NN3,NN4,m)
Nqy=MATMUL(NN4,NN1)

deallocate(NN3)
deallocate(NN4)
deallocate(NN1)
deallocate(IQY)
deallocate(tAT)

allocate(NN5(n1+n2,n1+n2))
!allocate(DNN5(n1+n2,n1+n2))
allocate(paortho(n1+n2,n1+n2))

NN5=MATMUL(At,Nqy)
!DNN5=0d0

!deallocate(At)
deallocate(Nqy)

print*,"contruction de paortho "

!do i=1,n1+n2
!   DNN5(i,i)=1d0
!enddo    
    
!  do j=1,n1+n2
!do i=1,n1+n2
!        paortho(i,j)=DNN5(i,j)-NN5(i,j)
!  enddo
!enddo

do j=1,n1+n2
  do i=1,n1+n2
     if(i==j)then
        paortho(i,i)=1d0-NN5(i,i)
     else
        paortho(i,j)=-NN5(i,j)
     endif
  enddo
enddo

deallocate(NN5)
!deallocate(DNN5)
!deallocate(At)
print*,'E====='
allocate(E(n1+n2))

E=MATMUL(paortho,omct)

!deallocate(omct)

ll=0d0
do i=1,n1
ll(1)=ll(1)+e(i)**2/c1**2
enddo
do i=n1+1,n2+n1
ll(2)=ll(2)+e(i)**2/c2**2
enddo
ll=ll*0.5d0
!print*,'ll =====',ll(1),ll(2),c1,c2

deallocate(e)

allocate(p1m(n1+n2,n1+n2))
allocate(p2m(n1+n2,n1+n2))
allocate(m11(n1,n1+n2))
allocate(m12(n2,n1+n2))

p1m=0d0
p2m=0d0
do j=1,n1+n2
    do i=1,n1
      p1m(i,j)=(1./c1)*paortho(i,j)
      m11(i,j)=p1m(i,j)
    enddo
   do i=n1+1,n1+n2
    p2m(i,j)=(1./c2)*paortho(i,j)
    m12(i-n1,j)=p2m(i,j)
   enddo
enddo


print*,("p1m ====")

deallocate(paortho)

!allocate(mN11(n1,n1+n2))
!allocate(mN12(n1,n1+n2))
!allocate(mN21(n2,n1+n2))
!allocate(mN22(n2,n1+n2))

print*,"mN11"

sum_mN11=0d0
do i=1,n1
           mmN11=0d0
             do k=1,n1
              mmN11=mmN11+m11(i,k)*p1m(k,i)
             enddo
             sum_mN11=sum_mN11+mmN11
enddo


!do i=1,n1
!       do j=1,n1+n2
!          mN11(i,j)=0d0
!            do k=1,n1
!             mN11(i,j)=mN11(i,j)+m11(i,k)*p1m(k,j)
!            enddo
!        enddo
!enddo  

print*,"mN12"

sum_mN12=0d0
do i=1,n1
           mmN12=0d0
             do k=n2,n1+n2
              mmN12=mmN12+m11(i,k)*p2m(k,i)
             enddo
            sum_mN12=sum_mN12+mmN12
enddo

!do i=1,n1
!       do j=1,n1+n2
!          mN12(i,j)=0d0
!            do k=n2,n1+n2
!             mN12(i,j)=mN12(i,j)+m11(i,k)*p2m(k,j)
!            enddo
!        enddo
!enddo  

print*,"mN21"

sum_mN21=0d0
do i=1,n2
           mmN21=0d0
             do k=n2,n1+n2
              mmN21=mmN21+m12(i,k)*p2m(k,i+n1)
             enddo
     sum_mN21=sum_mN21+mmN21
enddo 

!do i=1,n2
!       do j=1,n1+n2
!          mN21(i,j)=0d0
!            do k=n2,n1+n2
!             mN21(i,j)=mN21(i,j)+m12(i,k)*p2m(k,j)
!            enddo
!        enddo
!enddo 

print*,"mN22"

sum_mN22=0d0
do i=1,n2
           mmN22=0d0
             do k=1,n1
              mmN22=mmN22+m12(i,k)*p1m(k,i+n1)
             enddo
         sum_mN22=sum_mN22+mmN22
enddo

!do i=1,n2!
!       do j=1,n1+n2
!          mN22(i,j)=0d0
!            do k=1,n1
!             mN22(i,j)=mN22(i,j)+m12(i,k)*p1m(k,j)
!            enddo
!        enddo
!enddo

!print*,MATMUL(m11,p1m)

!print*,("mN22 ====")

!do i=1,n1
!print*,(mN21(i,i))
!enddo

deallocate(p1m)
deallocate(p2m)
deallocate(m11)
deallocate(m12)

!Nmat(1,1)=0.5d0*trace(mN11,n1,n1+n2,0)
!Nmat(1,2)=0.5d0*trace(mN12,n1,n1+n2,0)
!Nmat(2,1)=0.5d0*trace(mN22,n2,n1+n2,n1)
!Nmat(2,2)=0.5d0*trace(mN21,n2,n1+n2,n1)

Nmat(1,1)=0.5d0*sum_mN11
Nmat(1,2)=0.5d0*sum_mN12
Nmat(2,1)=0.5d0*sum_mN22
Nmat(2,2)=0.5d0*sum_mN21

print*,'Nmat'

!deallocate(mN11)
!deallocate(mN21)
!deallocate(mN22)
!deallocate(mN12)


!CALL LGINF(Nmat,2,2,2,0,INmat,2,S,WK,IER)
call inverse(Nmat,INmat,2)
print*,'inversion 2'
!print*,INmat(1,1),INmat(1,2),INmat(2,1),INmat(2,2)


newsig=MATMUL(INmat,ll)
print*,newsig
newsig=sqrt(newsig)

print*,'outsol',newsig*ua
print*,'cc',(cc(j)*ua,j=1,2)


if(abs(cc(1)-newsig(1))*ua<eps.and.abs(cc(2)-newsig(2))*ua<eps)then
print*,"resultat final ====",(cc(1)-newsig(1))*ua,(cc(2)-newsig(2))*ua
print*,newsig(1)*ua,newsig(2)*ua
print*,sqrt(sqrt(INmat(1,1)))*ua,sqrt(sqrt(INmat(2,2)))*ua
stop
else
print*,((newsig(j)*ua),j=1,2)

print*,'diff',(cc(1)-newsig(1))*ua,(cc(2)-newsig(2))*ua

do j=1,2
cc(j)=newsig(j)
enddo

print*,'nouveau c1,c2 ==='
print*,ii,(cc(j),j=1,2)

endif
enddo
End

! ===========================================

double precision function trace(Nmat,m,n,offset)
integer :: i,m,offset,n
real(kind=8), dimension(m,n) :: Nmat
trace=0d0
do i=1,m
trace=trace+Nmat(i,i+offset)
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


