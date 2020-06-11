PROGRAM MAIN

  integer, parameter :: height=142388, width=415
  real(kind=8), dimension(width) :: tmp, cnt

  cnt(:) = 0
open(3, file="./RAW.out.412.sort", status='old', action='read')

do i = 1, height
  read(3, *) tmp(:)
  cnt1 = 0
  do j = 4, width
    if(abs(tmp(j)) > 5e-3) cnt(j) = cnt(j) + 1
  end do
end do

do i = 1, width
  print*, cnt(i)
end do

close(3)

END PROGRAM MAIN
