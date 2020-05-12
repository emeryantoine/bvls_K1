PROGRAM MAIN

  integer, parameter :: height=141900, width=415
  real(kind=8), dimension(100,width) :: tmp
  integer, dimension(width-3) :: zeros, neg, smalls
  real(kind=8), dimension(width) :: tmp_RA
  integer :: change = 0

open(3, file="../../transfert/cas_complet/040520/RAW.out.412", status='old', action='read')
open(2, file="./zero.ppm", status="new", action="write")
write(2,'(a)')"P2"
write(2,'(a)')"412 1419"
write(2,'(a)')"100"

do i = 1, height, 100
  do j = 1, 100
    read(3, *) tmp(j, :)
  end do

  zeros(:) = 0
  do x = 4, width
    do y = 1, 100
      if (abs(tmp(y,x)) < 5e-3) then 
        zeros(x) = zeros(x) + 1
        change = change + 1
      endif
      if (abs(tmp(y,x)) .eq. 0d0) then
        change = change + 1
      endif
    end do
  enddo

  do x = 4,  width
    write(2, *) zeros(x)
  enddo
end do


print*, change, 141900*412

close(3)
close(2)

END PROGRAM MAIN
