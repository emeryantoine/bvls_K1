PROGRAM MAIN

  integer, parameter :: height=141970, width=344
  real(kind=8), dimension(height) :: W
  real(kind=8), dimension(height,width) :: J
  real(kind=8), dimension(width+3) :: reaad
  real(kind=8), dimension(width) :: tmp1
  real(kind=8), dimension(width, width) :: P


  open(1, file='../../transfert/RA.out', status='old', action='read')
  do i = 1, height
    read(1, *) reaad(:)
    W(i) = reaad(4)
    J(i,:) = reaad(4:347)
  end do

  tmp1 = matmul(transpose(J), W)
  P = matmul(transpose(tmp1), J)



END PROGRAM MAIN
