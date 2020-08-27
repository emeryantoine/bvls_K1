PROGRAM MAIN

include "lapack.f90"

integer, parameter :: height = 100, width = 10
real(kind=8), dimension(:), allocatable :: p, pp
real(kind=8), dimension(:,:), allocatable :: dp
real(kind=8), dimension(:,:), allocatable :: lim
real(kind=8), dimension(:,:), allocatable :: R, R_

allocate(p(width))
allocate(pp(width))
allocate(dp(width, width))
allocate(lim(2, width))
allocate(R(width, height))
allocate(R_(width, height))

do i = 1, width
  do j = 1, height
    R(i, j)  = i+j
  end do
end do

do i = 1, width
  lim(1, i) = i
  p(i) = i+1
  lim(2, i) = i+2
end do


!Step 1 transform p_i in p_i'
do i = 1, width
  pp(i) = log(p - lim(1, i)/lim(2, i) - p(i))
end do

!Step 2 Matrice diagonale des derivées partiels
dp = 0d0
do i = 1, width
  dp(i,i) = (lim(2, i) - lim(1, i) * ( exp(pp(i)) / (exp(pp(i)) + 1)**2 ))
end do

!Step 3 contruction of R_ such a R_ = R x dp
R_ = MATMUL(R, dp)

!Step 4 apply SVP from lapack or else
!bonus, we can get the covariance and corelation matrices here
call DGESVD()


!Step 5 retransforme pp in p
do i = 1, width
  p(i) = lim(1, i) + (lim(2,i) - lim(1, i))*(exp(pp(i) + X(i))/(exp(pp(i) + X(i)) + 1))
end do

END PROGRAM
