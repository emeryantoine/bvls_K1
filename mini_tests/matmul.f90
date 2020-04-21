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

do i = 1, 3
  do j = 1, 2
    aa(j, i) = i+j
  end do
enddo

do i = 1, 2
  do j = 1, 3
    bb(j, i) = -(i+j) +1
  end do
end do

cc = 0d0
CC = matmul(BB, AA)

print*, AA
print*, BB
Print*, CC

call mymatmul(bb, aa, DD, 2, 3)

print*, DD

end program matrixmul

subroutine mymatmul(b, a, c, width, height)
  
  integer :: width, height
  integer, dimension(width, height) :: a
  integer, dimension(height, width) :: b
  integer, dimension(height, height) :: c
 
  print*, shape(a)
  print*, shape(b)
  print*, shape(c)
  print*, width, height

  if(size(a, 1 ) /= width .or. size(b, 2) /= width) then
    print*, "width of A and height of B arn't compatible"
    stop
  endif
  if(size(a, 2) /= height .or. size(b, 1) /= height) then
    print*, "height of A and width of B arn't compatible"
    stop
  endif

  c = 0d0

  do i = 1, height
    do j = 1, height
      do k = 1,  width
        c(i, j) = c(i, j) + b(i, k)*a(k, j)
      end do
    end do
  end do

end subroutine mymatmul
