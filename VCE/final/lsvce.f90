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
!real(kind=8)::sum_mN11,mmN11,sum_mN12,mmN12,sum_mN21,mmN21,sum_mN22,mmN22
real(kind=8)::sum_mN11,sum_mN12,sum_mN21,sum_mN22
! real(kind=8), dimension(n2+n1,n2+n1) :: Q1,Q2,QY,IQY,aa,ab,mN22,mN21
!real(kind=8), dimension(n2+n1,n2+n1) :: NN5,paortho,DNN5,matr,mq1,mq2
!real(kind=8), dimension(n2+n1,n2+n1) :: mN11,mN12
!real(kind=8), allocatable, dimension(:,:) :: Q1,Q2,QY,IQY
real(kind=8), allocatable, dimension(:) :: IQY
real(kind=8), allocatable, dimension(:,:) :: NN5, NN52
real(kind=8), allocatable, dimension(:,:) :: p1m,p2m,m11,m12

!real(kind=8), allocatable, dimension(:,:) :: NN3,NN4,NN1,tAt,Nqy
real(kind=8), allocatable, dimension(:,:) :: NN3,NN4,NN1,Nqy
!real(kind=8), dimension(m,n2+n1) :: NN1,tAt,Nqy

real(kind=8), dimension(2,2) :: Nmat,INmat
real(kind=8), dimension(2) :: ll,newsig,cc,outsol,outsolavant
real(kind=8), parameter :: ua=1.5e11
character(100) :: valuem, filetoopen

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!open(355,file="dimension",status="old")
!read(355,*)m
!close(355)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!Lecture de la taille de m a utiliser depuis les ligne de commande
!gestion d'erruer si jamais le code n'est pas utilisee correctement

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!print*,"lecture",m
if(command_argument_count() .ne. 2) then
  write(*,*)"erreur, il faut preciser 2 arguments en lignes de commande tq :"
  write(*,*) "./xxx.out <taille de m> <chemin du fichier a lire>"
  stop
endif
call get_command_argument(1, valuem)
call get_command_argument(2, filetoopen)
read(valuem, *)m
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

allocate(A1(n1,m))
allocate(A2(n2,m))
allocate(AF(m))
allocate(omc1(n1))
allocate(omc2(n2))

open(355,file=filetoopen,status="old")
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

!inutile car on va parcourir toute les cases de toute facon donc aucune valeur
!non initialisee
!At=0d0
!omct=0d0

!===construction de A et omc"===

!$OMP PARALLEL DO schedule(dynamic, 512)
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
!$OMP END PARALLEL DO
   
deallocate(A1)
deallocate(A2)
deallocate(omc1)
deallocate(omc2)

!===remplissages des Q ===

cc(1)=1.9
cc(2)=1.3

do ii=1,niter

print*,'ITER',ii,'====', cc

c1=cc(1)*cc(1)
c2=cc(2)*cc(2)


!==========================================================
!allocate(Q1(n1+n2,n1+n2))
!allocate(Q2(n1+n2,n1+n2))
!allocate(QY(n1+n2,n1+n2))
!allocate(IQY(n1+n2,n1+n2))

!Q1=0d0
!do i=1,n1
! Q1(i,i)=1d0
!enddo

!Q2=0d0
!do i=n1+1,n1+n2
!Q2(i,i)=1d0
!enddo

!QY=0d0
!IQY=0d0

!print*,"contruction de QY et IQY"

!do j=1,n1+n2
! do i=1,n1+n2
!   QY(i,j)=c1*Q1(i,j)+c2*Q2(i,j)
!   IQY(i,j)=(1d0/c1)*Q1(i,j)+(1d0/c2)*Q2(i,j)
! enddo
!enddo
!deallocate(Q1)
!deallocate(Q2)
!deallocate(QY)
!==========================================================

!supprime Q1, Q2, QY
!remplace la matrice diagonale IQY par un vecteur diagonal

!==========================================================
allocate(IQY(n1+n2))
!$OMP PARALLEL DO schedule(dynamic, 512)
do i=1, n1
  IQY(i) = 1d0/c1
end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO schedule(dynamic, 512)
do i = n1+1, n1+n2
  IQY(i) = 1d0/c2
end do
!$OMP END PARALLEL DO
!==========================================================

allocate(NN3(m,m))
allocate(NN4(m,m))
allocate(NN1(m,n2+n1))
allocate(Nqy(m,n2+n1))

!---------------------------------------------------------
!allocate(tAT(m,n2+n1))

!tAt=0d0

!print*,"contruction de NQY "

!tAT=TRANSPOSE(At)
!NN1=MATMUL(tAt,IQY)
!---------------------------------------------------------

!Pour une transposee de matrice il suffit d'inverser les indice (i,j) -> (j,i)
!Reecriture de NN1= matmul(...) pour accomoder la forme custom de IQY diagonale
!linearisee

!---------------------------------------------------------
!$OMP PARALLEL DO schedule(dynamic, 512)
do i=1, n1+n2
  do j=1, m
      NN1(j, i) = At(i, j)*IQY(i)
  enddo
enddo
!$OMP END PARALLEL DO
!---------------------------------------------------------
NN3=MATMUL(NN1,At)
call inverse(NN3,NN4,m)
Nqy=MATMUL(NN4,NN1)

deallocate(NN3)
deallocate(NN4)
deallocate(NN1)
deallocate(IQY)
!deallocate(tAT)

allocate(NN5(n1+n2,n1+n2))
!allocate(DNN5(n1+n2,n1+n2))

print*, "matmul NN5 = At * NQY"
NN5=MATMUL(At,Nqy)
print*, "sizes", size(At, 1), size(At, 2), size(Nqy, 1), size(Nqy, 2)
!call mymatmul(At, Nqy)
!time



!deallocate(At)
deallocate(Nqy)

!***************************************************
!allocate(paortho(n1+n2,n1+n2))
!do j=1,n1+n2
!  do i=1,n1+n2
!     if(i==j)then
!        paortho(i,i)=1d0-NN5(i,i)
!     else
!        paortho(i,j)=-NN5(i,j)
!     endif
!  enddo
!enddo
!****************************************************

!Au lieu de creer une nouvelle matrice, on recycle NN5 en paortho

!****************************************************
print*,"contruction de paortho "
!NN5 = -NN5

!$OMP PARALLEL DO schedule(dynamic, 512)
do i = 1, n1+n2
  do j = 1, n1+n2
    NN5(j, i) = -NN5(j, i)
  end do
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO schedule(dynamic, 512)
do i = 1, n1+n2
  NN5(i,i) = NN5(i,i) + 1d0
end do
!$OMP END PARALLEL DO
!****************************************************

!deallocate(At)
print*,'E====='
allocate(E(n1+n2))

!time
E=MATMUL(NN5,omct)

!deallocate(omct)

ll=0d0
do i=1,n1
  ll(1)=ll(1)+(e(i)*e(i))/(c1*c1)
enddo
do i=n1+1,n2+n1
  ll(2)=ll(2)+(e(i)*e(i))/(c2*c2)
enddo
ll=ll*0.5d0
!print*,'ll =====',ll(1),ll(2),c1,c2

deallocate(e)

!.......................................................
!allocate(p1m(n1+n2,n1+n2))
!allocate(p2m(n1+n2,n1+n2))
!allocate(m11(n1,n1+n2))
!allocate(m12(n2,n1+n2))
!
!p1m=0d0
!p2m=0d0
!do j=1,n1+n2
!  do i=1,n1
!    p1m(i,j)=(1./c1)*NN5(i,j)
!    m11(i,j)=p1m(i,j)
!  enddo
!  do i=n1+1,n1+n2
!    p2m(i,j)=(1./c2)*NN5(i,j)
!    m12(i-n1,j)=p2m(i,j)
!  enddo
!enddo
!......................................................

!Chnager NN5 en 1./NN5 puis remplacer les utilisations de p1m p2m m11 et m12 par
!du NN5 avec des incides differents

!......................................................
print*, "== p1m, p2m dans NN5"
!$OMP PARALLEL DO schedule(dynamic, 512)
do i = 1, n1+n2
  do j = 1, n1
    NN5(i, j) = NN5(i, j)/c1
  end do
  do j = 1+n1, n1+n2
    NN5(i, j) = NN5(i, j)/c2
  end do
end do
!$OMP END PARALLEL DO
!......................................................


!deallocate(paortho)

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Version 1
!sum_mN11=0d0
!do i=1,n1
!  mmN11=0d0
!  do k=1,n1
!    mmN11=mmN11+m11(i,k)*p1m(k,i)
!  enddo
!  sum_mN11=sum_mN11+mmN11
!enddo
!
!sum_mN12=0d0
!do i=1,n1
!  mmN12=0d0
!  do k=n2,n1+n2
!    mmN12=mmN12+m11(i,k)*p2m(k,i)
!  enddo
!  sum_mN12=sum_mN12+mmN12
!end do
!
!sum_mN21=0d0
!do i=1,n2
!  mmN21=0d0
!  do k=n2,n1+n2
!    mmN21=mmN21+m12(i,k)*p2m(k,i+n1)
!  enddo
!  sum_mN21=sum_mN21+mmN21
!enddo
!
!sum_mN22=0d0
!do i=1,n2
!  mmN22=0d0
!  do k=1,n1
!    mmN22=mmN22+m12(i,k)*p1m(k,i+n1)
!  enddo
!  sum_mN22=sum_mN22+mmN22
!enddo
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooo

!On utilise la forme changee de NN5 pour le calcul des sommes et on ne passe pas
!par des variables temporaires inutiles

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Version 2
!print*,"mN11"
!
!sum_mN11=0d0
!do i=1,n1
!  do k=1,n1
!    sum_mN11=sum_mN11+NN5(i,k)*NN5(k,i)
!  enddo
!enddo
!
!print*,"mN12"
!
!sum_mN12=0d0
!do i=1,min(n1, n2)
!  do k=n2,n1+n2
!    sum_mN12=sum_mN12+NN5(i,k)*NN5(k,i)
!  enddo
!enddo
!
!print*,"mN21"
!
!sum_mN21=0d0
!do i=1+n1,n2+n1
!  do k=n2,n1+n2
!    sum_mN21=sum_mN21+NN5(i,k)*NN5(k,i)
!  enddo
!enddo 
!
!print*,"mN22"
!
!sum_mN22=0d0
!do i=1,n2
!  do k=1,n1
!    sum_mN22=sum_mN22+NN5(i+n1,k)*NN5(k,i+n1)
!  enddo
!enddo
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooo

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Version 3
print*,"mN11"

sum_mN11=0d0
!$OMP PARALLEL DO REDUCTION(+:sum_mN11) schedule(dynamic, 512)
do i=1,n1
  do k=1,n1
    sum_mN11=sum_mN11+NN5(i,k)*NN5(k,i)
  enddo
enddo
!$OMP END PARALLEL DO

print*,"mN12"

sum_mN12=0d0
!$OMP PARALLEL DO REDUCTION(+:sum_mN12) schedule(dynamic, 512)
do i=1,min(n1, n2)
  do k=n2,n1+n2
    sum_mN12=sum_mN12+NN5(i,k)*NN5(k,i)
  enddo
enddo
!$OMP END PARALLEL DO

print*,"mN21"

sum_mN21=0d0
!$OMP PARALLEL DO REDUCTION(+:sum_mN21) schedule(dynamic, 512)
do i=1+n1,n2+n1
  do k=n2,n1+n2
    sum_mN21=sum_mN21+NN5(i,k)*NN5(k,i)
  enddo
enddo 
!$OMP END PARALLEL DO

print*,"mN22"

sum_mN22=0d0
!$OMP PARALLEL DO REDUCTION(+:sum_mN22) schedule(dynamic, 512)
do i=1,n2
  do k=1,n1
    sum_mN22=sum_mN22+NN5(i+n1,k)*NN5(k,i+n1)
  enddo
enddo
!$OMP END PARALLEL DO
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooo

!deallocate(p1m)
!deallocate(p2m)
!deallocate(m11)
!deallocate(m12)
deallocate(NN5)

Nmat(1,1)=0.5d0*sum_mN11
Nmat(1,2)=0.5d0*sum_mN12
Nmat(2,1)=0.5d0*sum_mN22
Nmat(2,2)=0.5d0*sum_mN21

print*,'Nmat'

!CALL LGINF(Nmat,2,2,2,0,INmat,2,S,WK,IER)
call inverse(Nmat,INmat,2)
print*,'inversion 2'
!print*,INmat(1,1),INmat(1,2),INmat(2,1),INmat(2,2)


newsig=MATMUL(INmat,ll)
print*,newsig
newsig=sqrt(newsig)

print*,'outsol',newsig*ua
print*,'cc',(cc(j)*ua,j=1,2)


if(abs(cc(1)-newsig(1))*ua<eps.and.abs(cc(2)-newsig(2))*ua<eps) then
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
!==========================================================
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

subroutine mymatmul(a, b, c, width_a, height_a, width_b, height_b)
!=========================================================
!calcul une multiplication matricielle avec des matrices pas foorcement carre
!Prise en charge des erreurs si les tailles de matrice ne sont pas conforme
!
!Calcul C = A*B
!
!A et B sont initialisee et rempli de valeurs qui serviront a calculer C
!C est alloue et sera remplie de 0 durant l'execution de cette subroutine
!width_a = size(A, 1)
!height_a = size(A, 2)
!tel que A est initialisee par XXX, dimension(width, height) :: A
!=========================================================


  integer :: width_a, height_a, width_b, height_b
  real(kind=8), dimension(width_a, height_a) :: a
  real(kind=8), dimension(width_b, height_b) :: b
  real(kind=8), dimension(width_a, height_b) :: c

  if(width_a /= height_b) then
    print*, "width of A incompatible with height of B"
    stop
  endif
  
  c= 0d0

  do i = 1, height_b
    do j = 1, width_a
      do k = 1, height_a
        c(j, i) = c(j, i) + a(j, k)*b(k, i)
      end do
    end do
  end do

end subroutine mymatmul
