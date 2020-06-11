program main

integer(kind = 4), parameter :: m = 1500, n=12
real(kind = 8), dimension(:,:), allocatable :: A
real(kind=8), dimension(:), allocatable :: B, X, R

allocate(A(m, n))
allocate(B(m))
allocate(R(m))
allocate(X(n))

R(:) = 0
X(:) = 0

do i = 1, m
  B(i) = i
  do j = 1, n
    A(i, j) = i*j
  enddo
enddo

call qr_solve(m, n, a, b, x, r)

do i = 1, m
  print*, R(i), X(i)
enddo

end program main
