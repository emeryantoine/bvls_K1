program matrixmul

integer, dimension(2, 2) :: A, B, C, D

integer, dimension(2, 3) :: BB
integer, dimension(3, 2) :: AA
integer, dimension(3, 3) :: CC
integer, dimension(3, 3) :: DD

A(1, 1) = 1
A(1, 2) = 2
A(2, 1) = 3
A(2, 2) = 4

B = A*2

C = matmul(A, B)

D = 0d0

do i = 1, 2
  do j = 1, 2
    do k = 1, 2
      D(j, i) = D(j, i) + A(j, k)*B(k, i)
    end do
  end do
end do

print*, A
print*, B
Print*, C
print*, D

do i = 1, 3
  do j = 1, 2
    bb(j, i) = i+j
  end do
enddo

do i = 1, 2
  do j = 1, 3
    aa(j, i) = -(i+j) +1
  end do
end do

cc = 0d0
CC = matmul(aa, bb)

print*, AA
print*, BB
Print*, CC

call mymatmul(aa, bb, DD, 3, 2, 2, 3)

print*, DD

end program matrixmul

subroutine mymatmul(a, b, c, width_a, height_a, width_b, height_b)
  
  integer :: width_a, height_a, width_b, height_b
  integer, dimension(width_a, height_a) :: a
  integer, dimension(width_b, height_b) :: b
  integer, dimension(width_a, height_b) :: c
 
  if(width_a /= height_b) then
    print*, "width of A and height of B arn't compatible"
    stop
  endif

  c = 0d0

  do i = 1, height_b
    do j = 1, width_a
      do k = 1,  height_a
        c(j, i) = c(j, i) + a(j, k)*b(k, i)
      end do
    end do
  end do

end subroutine mymatmul
