program main

real(kind=16) :: val = 1, res, tmp

res = exp(val)

do

  val = val + 1
  tmp = res
  res = exp(val)
  print*, val, res

  if(tmp .eq. res) EXIT

end do

end program
