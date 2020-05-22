PROGRAM MAIN

  integer, parameter :: width=415
  real(kind=8), dimension(100,width) :: tmp
  integer, dimension(width-3) :: res
  integer :: change = 0, status

status = system("rm origin.ppm")

open(3, file="./RAW.out.412.sort", status='old', action='read')
open(2,file="./origin.ppm", status="new", action="write")
write(2,'(a)') "P2"
write(2, '(a)') "412 1420"
write(2,'(a)') "100"

do i = 1, 1420
  do k = 1, 100
    read(3,*) tmp(k,:)
    counter = counter + 1
  end do

  res(:) = 0d0

  do j = 4, width
    do k = 1, 100
      if(abs(tmp(k,j)) > 5e-3) res(j) = res(j) + 1
    end do
  end do

  do k = 1, width-3
    write(2,*) res(k)
  end do

end do
close(2)
close(3)


END PROGRAM MAIN
