!mpifort matmul.f90 ../VCE/calculmpi.f90 -fopenmp
!mpiexec -n 4 --map-by ppr:1:socket ./a.out

subroutine matmatompi(a, b, c, dim1_a, dim2_a, dim1_b, dim2_b)
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

  !implicit none

  integer :: dim1_a, dim2_a, dim1_b, dim2_b
  integer :: i, j, k, omp_thr, rank, world
  integer :: ssdim1, ssdim2
  real(kind=8), dimension(dim1_a, dim2_a) :: a
  real(kind=8), dimension(dim1_b, dim2_b) :: b
  real(kind=8), dimension(dim1_a, dim2_b) :: c
  real(kind=8), dimension(:,:), allocatable :: tempo
  
  if(rank == 0) then
    print*, "A", A
    print*, "B", B
  endif

  if(dim1_a /= dim2_b) then
    print*, "width of A incompatible with height of B"
    stop
  endif

  call MPI_INIT(NULL, NULL)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, world)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank)
  
  ssdim1 = dim1_a

  if(rank == 0) then
    ssdim2 = dim2_b/world + modulo(dim2_b, world)
  else
    ssdim2 = dim2_b/world
  end if
  
  allocate(tempo(ssdim1, ssdim2))
  tempo = 0d0

  !$OMP PARALLEL 
  omp_thr = omp_get_num_threads()
  if(dim2_b > omp_thr) then
    !$OMP DO schedule(runtime)
    do i = 1, ssdim2
      do j = 1, ssdim1
        do k = 1, dim2_a
          tempo(j, i) = tempo(j, i) + a(j, k + ssdim2*rank)*b(k + ssdim2*rank, i)
        end do
      end do
    end do
    !$OMP END DO
  else
    do i = 1, ssdim2
    !$OMP DO schedule(runtime)
      do j = 1, ssdim1
        do k = 1, dim2_a
          tempo(j, i) = tempo(j, i) + a(j, k + ssdim2*rank)*b(k + ssdim2*rank, i)
        end do
      end do
    !$OMP END DO
    end do
  endif 
  !$OMP END PARALLEL

  !print*, rank, tempo
  !if(rank == 0) then
  !  print*, "A", A
  !  print*, "B", B
  !end if

  !reconstruct C out of all the tempo from each processuses

  call MPI_FINALIZE()
  if(rank /= 0) stop

end subroutine matmatompi


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
