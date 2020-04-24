subroutine matmat(a, b, c, dim1_a, dim2_a, dim1_b, dim2_b)
!=========================================================
!calcul une multiplication matricielle avec des matrices pas
!foorcement carre
!Prise en charge des erreurs si les tailles de matrice ne sont pas
!conforme
!
!Calcul C = A*B
!
!A et B sont initialisee et rempli de valeurs qui serviront a
!calculer C
!C est alloue et sera initialise a 0 avant l'execution de cette
!subroutine
!dim1_a = size(A, 1)
!dim2_a = size(A, 2)
!tel que A est initialisee par XXX, dimension(width, height) :: A
!=========================================================

  use omp_lib
  include 'mpif.h'

  implicit none


  integer :: dim1_a, dim2_a, dim1_b, dim2_b, i, j, k, nbr, rank, sizze
  real(kind=8), dimension(dim1_a, dim2_a) :: a
  real(kind=8), dimension(dim1_b, dim2_b) :: b
  real(kind=8), dimension(dim1_a, dim2_b) :: c

  call MPI_INIT()
  call MPI_COMM_SIZE(MPI_COMM_WORLD, sizze)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank)

  if(rank == 0) then
   !ergvbiuboiuvefrwbiuovefwbiuoerfwbiobuebriebvfiviebfnijevbjievjni 
  end if

  
  if(dim1_a /= dim2_b) then
    print*, "width of A incompatible with height of B"
    stop
  endif
  c = 0d0
  !$OMP PARALLEL
  
  nbr = omp_get_num_threads()
  if(dim2_b > nbr) then
    !$OMP DO schedule(runtime)
    do i = 1, dim2_b
      do j = 1, dim1_a
        do k = 1, dim2_a
          c(j, i) = c(j, i) + a(j, k)*b(k, i)
        end do
      end do
    end do
    !$OMP END DO
  else
    do i = 1, dim2_b
    !$OMP DO schedule(runtime)
      do j = 1, dim1_a
        do k = 1, dim2_a
          c(j, i) = c(j, i) + a(j, k)*b(k, i)
        end do
      end do
    !$OMP END DO
    end do
  endif 
  !$OMP END PARALLEL

  call MPI_FINALIZE()

end subroutine matmat

subroutine matvect(a, b, c, dim1_a, dim2_a, dim_b)
!=========================================================
!calcul une multiplication matrice vecteur donc le resultats est un vecteur de
!meme taille que B
!Prise en charge des erreurs si les tailles de matrice ne sont pas
!conforme
!
!Calcul C = A*B
!
!A et B sont initialisee et rempli de valeurs qui serviront a
!calculer C
!C est alloue et sera initialise 0 avant l'execution de cette
!subroutine
!dim1_a = size(A, 1)
!dim2_a = size(A, 2)
!tel que A est initialisee par XXX, dimension(width, height) :: A
!=========================================================

  use omp_lib
  implicit none

  integer :: dim1_a, dim2_a, dim1_b, dim_b, i, j
  real(kind=8), dimension(dim1_a, dim2_a) :: a
  real(kind=8), dimension(dim_b) :: b, c

  if(dim1_a /= dim_b) then
    print*, "width of A incompatible with height of B"
    stop
  endif

  c = 0d0
  !$OMP PARALLEL DO schedule(runtime)
  do i = 1, dim1_a
    do j = 1, dim2_a
      c(i) = c(i) + a(i, j)*b(j)
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine matvect

subroutine matmatMPIomp(a, b, c, dim1_a, dim2_a, dim1_b, dim2_b)
  implicit none

  integer :: dim1_a, dim2_a, dim1_b, dim2_b, i, j, k, tt
  real(kind=8), dimension(dim1_a, dim2_a) :: a
  real(kind=8), dimension(dim1_b, dim2_b) :: b
  real(kind=8), dimension(dim1_a, dim2_b) :: c

  

end subroutine matmatMPIomp
