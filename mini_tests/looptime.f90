program main

implicit none

integer, parameter :: sizee = 1

real, dimension(sizee) :: X, Y

integer :: start, stope, i, j, start2, stope2


call system_clock(start)
do j=0, 1
do i=1,sizee,8
X(i) = Y(i)
X(i+1) = Y(i+1)
X(i+2) = Y(i+2)
X(i+3) = Y(i+3)
X(i+4) = Y(i+4)
X(i+5) = Y(i+5)
X(i+6) = Y(i+6)
X(i+7) = Y(i+7)
end do
end do
call system_clock(stope)
print *, (stope-start)

call system_clock(start2)
do j=0, 100000
Y(1:sizee) = X(1:sizee)
end do
call system_clock(stope2)
print *, (stope2-start2)

print *, abs(-42)

end program main
