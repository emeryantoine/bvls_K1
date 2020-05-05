PROGRAM MAIN

  integer, parameter :: height=141900, width=415
  real(kind=8), dimension(100,width) :: tmp
  integer, dimension(width-3) :: zeros, neg, smalls
  real(kind=8), dimension(width) :: tmp_RA
  integer :: zero=0, small=0, supone=0, total=0, change = 0

  open(1, file="../../transfert/cas_complet/RA.out", status='old', action='read')
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

open(3, file="../../transfert/cas_complet/RA.out", status='old', action='read')
open(20, file="./zero.ppm", status="new", action="write")
open(21, file="./neg.ppm", status="new", action="write")
open(22, file="./small.ppm", status="new", action="write")
write(20,'(a)')"P2"
write(20,'(a)')"412 1419"
write(20,'(a)')"100"
write(21,'(a)')"P2"
write(21,'(a)')"412 1419"
write(21,'(a)')"100"
write(22,'(a)')"P2"
write(22,'(a)')"412 1419"
write(22,'(a)')"100"

do i = 1, height, 100
  do j = 1, 100
    read(3, *) tmp(j, :)
  end do

  zeros(:) = 0
  neg(:) = 0
  smalls(:) = 0
  do x = 4, width
    do y = 1, 100
      if (tmp(y, x) == 0) zeros(x) = zeros(x) + 1
      if (tmp(y,x) < 0) neg(x) = neg(x) + 1
      if (abs(tmp(y,x)) < 10e-3) then 
        smalls(x) = smalls(x) + 1
        change = change + 1
      endif
    end do
  enddo

  do x = 4,  width
    write(20, *) zeros(x)
    write(21, *) neg(x)
    write(22, *) smalls(x)
  enddo

end do

print*, change, 141900*74

close(3)
close(20)
close(21)
close(22)

print*, "zero   small   supone   total-small&zero   total"
print*, zero, small, supone, total - small - zero, total
print*, real(zero)/real(total)*100, real(small)/real(total)*100, &
  real(supone)/real(total)*100, real(total - small - zero)/real(total)*100, "100%"

END PROGRAM MAIN
