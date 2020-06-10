program main

integer, parameter :: width = 12, height = 142388
real(kind=8), dimension(width, width) :: rtilde, res1, res2, rtildet
real(kind=8), dimension(height, width) :: A
real(kind=8), dimension(width, height) :: At

open(1, file="./Rtilde", status="old", action="read")
do i = 1, width
  read(1,*) rtilde(i, :)
end do
close(1)

open(2, file="../../transfert/cas_complet/040520/RAW.out.412", status="old", &
action="read")
do i = 1, height
  read(2, *) A(i, :)
enddo
close(2)

res1(:,:) = 0
res2(:,:) = 0

do i = 1, width
  do j = 1, width
    do k = 1, height
      res1(j,i) = res1(j,i) + A(k, j)*A(k,i)
    enddo
  enddo
enddo

do i = 1, width
  do j = 1, width
    do k = 1, width
      res2(j,i) = res2(j,i) + rtilde(k, j)*rtilde(k, i)
    enddo
  enddo
enddo

!do i = 1, width
!  do j = 1, width
!    if(res1(i,j) /= res2(i,j)) print*, "diff", i, j
!  enddo
!enddo

do i = 1, width
  do j = 1, width
    print*, res1(i,j), res2(i,j)
  enddo
enddo

end program main
