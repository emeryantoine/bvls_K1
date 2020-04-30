PROGRAM MAIN

  integer, parameter :: height=141900, width=77
  real(kind=8), dimension(100,width) :: tmp
  integer, dimension(width-3) :: vals
  real(kind=8), dimension(width) :: tmp_RA
  integer :: zero=0, small=0, supone=0, total=0

  open(1, file="../../transfert/RA.out", status='old', action='read')
  do i = 1, height
    read(1, *) tmp_RA(:)
    do j = 4, width
      total = total + 1
      if (tmp_RA(j) == 0) then
        zero = zero + 1
      elseif (tmp_RA(j) > 1) then
        supone = supone + 1
      elseif (tmp_RA(j) < 1e-12) then
        small = small + 1
      end if
    end do
  end do
close(1)

open(1, file="../../transfert/RA.out", status='old', action='read')
open(2, file="./zero.ppm", status="new", action="write")
write(2,'(a)')"P2"
write(2,'(a)')"74 14190"
write(2,'(a)')"100"

do i = 1, height, 100
  do j = 1, 100
    read(1, *) tmp(j, :)
  end do
  
  vals = 0
  do x = 4, width
    do y = 1, 100
      if (tmp(y, x) == 0) vals(y) = vals(y) + 1
    end do
  enddo

  do x = 1, 100
    write(2, *) vals(x)
  enddo

end do

close(1)
close(2)

print*, "zero   small   supone   total-small&zero   total"
print*, zero, small, supone, total - small - zero, total
print*, real(zero)/real(total)*100, real(small)/real(total)*100, &
  real(supone)/real(total)*100, real(total - small - zero)/real(total)*100, "100%"

END PROGRAM MAIN
