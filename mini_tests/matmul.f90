program matrixmul

integer, dimension(2, 2) :: A, B, C, D

integer, dimension(2, 3) :: AA
integer, dimension(3, 2) :: BB
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

CC = matmul(AA, BB)

call mymatmul(BB, AA, DD)

print*, AA
print*, BB
Print*, CC
print*, DD

end program matrixmul

function mymatmul(a, b, c)
  if(shape(A)(1) /= shape(B)(2)) then
    print*, "shapes of A(1) and B(2) arn't compatible"
    stop
  endif
  if(shape(A)(2) /= shape(B)(1)) then
    print*, "shape of A(2) and B(1) arn't compatible"
    stop
  endif

  c = 0d0

  do i = 1, shape(A)(2)
    do j = 1, shape(A)(1)
      do k = 1,  shape(A)(1)
        c(j, i) = c(j, i) + a(j, k)*b(k, i)
      end do
    end do
  end do

end function
