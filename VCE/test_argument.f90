! ifort -O3 lsvce.f90 -o lsvce ; lsvce
implicit none
integer :: i,j,k,ii,m,nc,ntotal,bb,ik,xx
integer, parameter:: niter=15
real(kind=8), parameter :: eps=1e-4
real(kind=8), allocatable, dimension(:,:) :: A1,A2,At
real(kind=8), allocatable, dimension(:) :: omct,E,omc1,omc2,AF
real(kind=8) :: dat,wpond,c1,c2,omcf
real(kind=8):: norm
real(kind=8), allocatable, dimension(:) :: IQY,ll,newsig,cc
real(kind=8), allocatable, dimension(:,:) :: NN5,Nmat,INmat
real(kind=8), allocatable, dimension(:,:) :: p1m,p2m,m11,m12
real(kind=8), allocatable, dimension(:,:) :: NN3,NN4,NN1,Nqy
real(kind=8), parameter :: ua=1.5e11

character(100) :: valuem, filetoopen,valuenc
character(100), dimension(10) :: valuencn
integer, allocatable, dimension(:) :: nset

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!open(355,file="dimension",status="old")
!read(355,*)m
!close(355)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!Lecture de la taille de m a utiliser depuis les ligne de commande
!gestion d'erruer si jamais le code n'est pas utilisee correctement

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!print*,"lecture",m
!if(command_argument_count() .ne. 2) then
!  write(*,*)"erreur, il faut preciser 2 arguments en lignes de commande tq :"
!  write(*,*) "./xxx.out <taille de m> <chemin du fichier a lire>"
!  stop
!endif
call get_command_argument(1, valuem)
call get_command_argument(2, valuenc)
call get_command_argument(3, filetoopen)

read(valuem, *)m
read(valuenc, *)nc

allocate(nset(nc))

open(345,file="nset.in",status="old")
read(345,*)(nset(i),i=1,nc)
close(345)

print*,m
print*,nc
print*,filetoopen
print*,nset

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ntotal=0
do i=1,nc
ntotal=ntotal+nset(i)
enddo

print*,ntotal,m

allocate(omct(ntotal))
allocate(At(ntotal,m))

open(355,file=filetoopen,status="old")
do i=1,ntotal
read(355,*)dat,wpond,omct(i),(At(i,j),j=1,m)
enddo
close(355)

allocate(cc(nc))

do i=1,nc
cc(i)=i*2d0
enddo

do ii=1,niter

	print*,'ITER',ii,'====', cc
	
	allocate(IQY(ntotal))
	
	IQY=0d0
	
	bb=0
	do j=1,nc
		do i=1,nset(j)
		IQY(i+bb)=1d0/(cc(j)*cc(j))
		enddo
		bb=bb+nset(j)
	enddo
	
	
	allocate(NN3(m,m))
	allocate(NN4(m,m))
	allocate(NN1(m,ntotal))
	allocate(Nqy(m,ntotal))
	
  !$OMP PARALLEL DO schedule(dynamic, 512)
	do i=1, ntotal
	  do j=1, m
	      NN1(j, i) = At(i, j)*IQY(i)
	  enddo
	  !print*,(At(i,j),j=1,m)
	enddo
  !$OMP END PARALLEL DO
	!---------------------------------------------------------
	NN3=MATMUL(NN1,At)
	call inverse(NN3,NN4,m)
	Nqy=MATMUL(NN4,NN1)
	
	deallocate(NN3)
	deallocate(NN4)
	deallocate(NN1)
	
	allocate(NN5(ntotal,ntotal))
	NN5=0d0
	print*, "matmul NN5 = At * NQY"
	!NN5=MATMUL(At,Nqy)
  call matmat(At, Nqy, NN5, size(At, 1), size(At, 2), size(Nqy, 1), size(Nqy, 2))
	deallocate(Nqy)
	
	
	!****************************************************
	print*,"contruction de paortho "
	!NN5 = -NN5
	
!	!$OMP PARALLEL DO schedule(dynamic, 512)
	do i = 1, ntotal
	  do j = 1, ntotal
	    NN5(j, i) = -NN5(j, i)
	  end do
	end do
!	!$OMP END PARALLEL DO
	
!	!$OMP PARALLEL DO schedule(dynamic, 512)
	do i = 1, ntotal
	  NN5(i,i) = NN5(i,i) + 1d0
	end do
!	!$OMP END PARALLEL DO
	!****************************************************
	
	print*,'E====='
	allocate(E(ntotal))
	
	
	!E=MATMUL(NN5,omct)
	call matvect(NN5, omct, E, size(NN5, 1), size(NN5, 2), size(omct, 1))

	!deallocate(omct)
	print*,'apres E'
	
	allocate(ll(nc))
	
	ll=0d0
	bb=0
	do j=1,nc
		do i=1,nset(j)
		  ll(j)=ll(j)+(e(i+bb)*e(i+bb))/(cc(j)*cc(j))**2
!		  print*,i,i+bb,e(i+bb),cc(j)
		enddo
		bb=bb+nset(j)
		print*,'bb',bb,nset(j)
	enddo
	ll=ll*0.5d0

print*,'ll ====='!,(j,ll(j),cc(j),j=1,nc)
	
	deallocate(e)


!......................................................

print*, "== p1m, p2m dans NN5"
!$OMP PARALLEL DO schedule(dynamic, 512)
	do i = 1, ntotal
!			write(*,'(11(f10.5,1x))')IQY(i),(NN5(i,j),j=1,ntotal)
	bb=0
	   do j=1,nc
		 do k=1,nset(j)
			NN5(i,k+bb)=NN5(i,k+bb)*IQY(k+bb)
		 enddo
		 bb=bb+nset(j)
		enddo
!		write(*,'(10(f10.5,1x))')(NN5(i,j),j=1,ntotal)

	enddo
!  do j = 1, n1
!    NN5(i, j) = NN5(i, j)/c1
!  end do
!  do j = 1+n1, ntotal
!    NN5(i, j) = NN5(i, j)/c2
!  end do
!end do
!$OMP END PARALLEL DO	
	
	deallocate(IQY)
	allocate(Nmat(nc,nc))
	
	
	Nmat=0d0
	xx=0
	do ik=1,nc
		bb=0
		do j=1,nc
			do i=1,nset(ik)
				do k=1,nset(j)
				Nmat(ik,j)=Nmat(ik,j)+NN5(i+xx,k+bb)*NN5(k+bb,i+xx)
				enddo
			enddo
		bb=bb+nset(j)
		enddo
		xx=xx+nset(ik)
	enddo
	Nmat=0.5d0*Nmat
	
	deallocate(NN5)
print*,'Nmat ====='
do i=1,2
write(*,'(2(f10.5,1x))')(Nmat(i,j),j=1,2)
enddo

allocate(INmat(nc,nc))
call inverse(Nmat,INmat,nc)

deallocate(Nmat)
allocate(newsig(nc))

newsig=MATMUL(INmat,ll)
print*,'newsig',newsig
newsig=sqrt(newsig)

print*,'newsig ua',newsig*ua
print*,'cc',(cc(j)*ua,j=1,nc)

norm=0
do i=1,nc
norm=norm+(cc(i)-newsig(i))*(cc(i)-newsig(i))
enddo
print*,sqrt(norm)*ua,eps
if(sqrt(norm)*ua<eps) then
  print*,"resultat final ===="
  do i=1,nc
!	print*,i,(cc(i)-newsig(i))*ua
	print*,newsig(i)*ua,sqrt(sqrt(INmat(i,i)))*ua
  enddo
  stop
else

do j=1,nc
  cc(j)=newsig(j)
enddo

print*,'nouveau c1,c2 ==='
print*,ii,(cc(j),j=1,nc)

endif
deallocate(INmat)
deallocate(newsig)
deallocate(ll)

enddo



deallocate(nset)
deallocate(cc)
deallocate(At)
deallocate(omct)

End


! ===========================================

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

subroutine mymatmul(a, b, c, dim1_a, dim2_a, dim1_b, dim2_b)
!=========================================================
!calcul une multiplication matricielle avec des matrices pas foorcement carre
!Prise en charge des erreurs si les tailles de matrice ne sont pas conforme
!
!Calcul C = A*B
!
!A et B sont initialisee et rempli de valeurs qui serviront a calculer C
!C est alloue et sera remplie de 0 durant l'execution de cette subroutine
!dim1_a = size(A, 1)
!dim2_a = size(A, 2)
!tel que A est initialisee par XXX, dimension(width, height) :: A
!=========================================================

implicit none

integer :: dim1_a, dim2_a, dim1_b, dim2_b, OMP_GET_NUM_THREADS, i, j, k
real(kind=8), dimension(dim1_a, dim2_a) :: a
real(kind=8), dimension(dim1_b, dim2_b) :: b
real(kind=8), dimension(dim1_a, dim2_b) :: c

!print*, dim1_a, dim2_a, dim1_b, dim2_b
!print*, shape(c)
              
if(dim1_a /= dim2_b .and. dim2_b /= 1) then
  print*, "width of A incompatible with height of B"
  stop
endif

!$OMP PARALLEL
if(dim2_b > OMP_GET_NUM_THREADS()) then
  !$OMP DO schedule(dynamic, 512)
  do i = 1, dim2_b
    do j = 1, dim1_a
      do k = 1, dim2_a
        c(j, i) = c(j, i) + a(j, k)*b(k, i)
      end do
    enddo
  enddo
else
  do i = 1, dim2_b
  !$OMP DO schedule(dynamic, 512)
    do j = 1, dim1_a
      do k = 1, dim2_a
        c(j, i) = c(j, i) + a(j, k)*b(k, i)
      end do
    end do
  !$OMP END DO
  end do
endif 
!$OMP END PARALLEL

end subroutine mymatmul
