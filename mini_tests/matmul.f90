program matrixmul

integer, dimension(2, 2) :: A, B, C, D
integer, dimension(8,8) :: AAA, BBB, E

integer, dimension(2, 3) :: BB
integer, dimension(3, 2) :: AA
integer, dimension(3, 3) :: CC
integer, dimension(3, 3) :: DD
integer, dimension(1, 10) :: test

do i = 1, 8
  do j = 1,8 
    AAA(j, i) = i*j + i
  end do
end do

BBB = AAA*2

print*, "A", AAA
print*, "B", BBB

A(1, 1) = 1
A(1, 2) = 2
A(2, 1) = 3
A(2, 2) = 4

B = A*2

C = matmul(A, B)

call matmat(AAA, BBB, E, 8, 8, 8, 8)

D = 0d0

do i = 1, 2
  do j = 1, 2
    do k = 1, 2
      D(j, i) = D(j, i) + A(j, k)*B(k, i)
    end do
  end do
end do

print*, "A", A
print*, "B", B
Print*, "C", C
print*, "D", D

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

call matmat(aa, bb, DD, 3, 2, 2, 3)

stop

print*, DD

end program matrixmul
