PROGRAM MAIN

  integer, parameter :: width=412
  real(kind=8), dimension(100,width) :: tmp
  integer, dimension(width) :: res
  integer :: change = 0, status

status = system("rm outputA.ppm")

open(3, file="../updatedBVLS/outputA.out", status='old', action='read')
open(2,file="./outputA.ppm", status="new", action="write")
write(2,'(a)') "P2"
write(2, '(a)') "412 1420"
write(2,'(a)') "100"

do i = 1, 1420
  do k = 1, 100
    read(3,*) tmp(k,:)
  end do

  res(:) = 0d0

  do j = 1, width
    do k = 1, 100
      if(abs(tmp(k,j)) > 0) res(j) = res(j) + 1
    end do
  end do

  do k = 1, width
    write(2,*) res(k)
  end do

end do
close(2)
close(3)


END PROGRAM MAIN
