PROGRAM MAIN

  integer, parameter :: height=25013, width=415
  real(kind=8), dimension(width) :: tmp

open(3, file="../../transfert/cas_complet/040520/RAW.out.412", status='old', action='read')
open(2, file="./zero.ppm", status="new", action="write")
write(2,'(a)')"P2"
write(2,'(a)')"412 142388"
write(2,'(a)')"100"

do i = 1, height
  read(3, *) tmp(:)
  
  do j = 4, width
    if(abs(tmp(j)) > 5e-3) then
      write(2,*) 100
    else
      write(2,*) 0
    end if
  end do
end do


print*, change, 141900*412

close(3)
close(2)

END PROGRAM MAIN
