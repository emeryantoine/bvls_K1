PROGRAM MAIN

  integer, parameter :: height=141970, width=77
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
print*, "zero   small   supone   total-small&zero   total"
print*, zero, small, supone, total - small - zero, total
print*, real(zero)/real(total)*100, real(small)/real(total)*100, &
  real(supone)/real(total)*100, real(total - small - zero)/real(total)*100, "100%"

  open(2, file="../../transfert/fort.299", status='old', action='read')

END PROGRAM MAIN
