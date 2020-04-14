program main

integer, parameter :: limit = 20, height=141970, width = 77, offset = -3
real, dimension(limit) :: sizes
real(kind=8), dimension(:), allocatable :: tmp, W
real :: val, expo, maxi = -1
integer :: planet = -1, nulle = 0, supone = 0, infone = 0, z1 = 0, zz1=0, zzz1=0, zzzz1=0
integer :: zzzzz1 = 0, zzzzzz1 = 0, zzzzzzz1=0
integer :: percent = (width-1)*height
real(kind=8), dimension(16) :: try

sizes = 0
try = 0

allocate(tmp(width))
allocate(W(height))

open(1, file='../../transfert/RA.out', status='old', action='read')
  do i = 1, height
    read(1, *) tmp(:)
    W(i) = tmp(2)
    do j=1, width + offset
      val = abs(tmp(j - offset))
      incr: do k=1, limit
      expo = k
      if (expo == limit) then
        if(val == 0) then
          sizes(k) = sizes(k) + 1
        endif
      end if
      if (val < 10**(-expo)) then
        sizes(k) = sizes(k) + 1
      else
        exit incr
      endif
        end do incr
    end do

    if(i >= 69346) then
      do j = 1, width + offset
        if(j >= 63) then
          do k=1, 16
            expo = k
            val = abs(tmp(j))
            if(val < 10**(-expo)) then
              try(k) = try(k) + 1
            end if
          end do
        end if
      end do
    end if
  end do

do i=1, 15
  try(i) = try(i) - try(i+1)
end do

do i =1, 16
  print*, try(i), "between 10^", -i, "and 10^", -(i+1)
end do
print*, "max value", maxi

close(1)
!print *, "there's", nonzerol, nonzeror, "non zero element in the bottom right section of matrix A"
do i = 1, limit
  print *, sizes(i), "elements are < 10^",-i, " --> ", (sizes(i)/(height*width))*100, "%"
end do

!tmp = realloc(tmp, 7)

open(2, file='../../transfert/fort.299', status='old', action='read')
do i=1, height
  read(2, *) tmp(1:7)
  if (tmp(2) /= planet) then
    planet = tmp(2)
    print *, "planet : ", planet, "starts at line : ", i
  end if
end do

close (2)

open(3, file='../../transfert/RA.out', status='old', action='read')
do i=1, height
 read(3, *) tmp(:)
 do j=2, width
  if(tmp(j) == 0) nulle = nulle + 1
  val = abs(tmp(j))
  if(val >= 1) supone = supone+1
  if(val < 1 .and. val >= 0.1) z1 = z1 + 1
  if(val < 0.1 .and. val >= 0.01) zz1 = zz1 + 1
  if(val < 0.01 .and. val >= 0.001) zzz1 = zzz1 + 1
  if(val < 0.001 .and. val >= 0.0001) zzzz1 = zzzz1 + 1
  if(val < 0.0001 .and. val >= 0.00001) zzzzz1 = zzzzz1 + 1
  if(val < 0.00001 .and. val >= 0.000001) zzzzzz1 = zzzzzz1 + 1
  if(val < 0.000001 .and. val >= 0.0000001) zzzzzzz1 = zzzzzzz1 + 1
  if(val > maxi) maxi = val 
 end do
end do
close(3)
print*, nulle,"zero elements"
print *, real(nulle)/real(percent) * 100, "%"
print*, supone, "elements higher than one and"
print*, real(supone)/real(percent)*100, "%"
print*, z1, "element betwen 1 ans 0.1"
print*, real(z1)/real(percent) * 100, "%"
print*, zz1, "element betwen 0.1 ans 0.01"
print*, real(zz1)/real(percent) * 100, "%"
print*, zzz1, "element betwen 0.01 ans 0.001"
print*, real(zzz1)/real(percent) * 100, "%"
print*, zzzz1, "element betwen 0.001 ans 0.0001"
print*, real(zzzz1)/real(percent) * 100, "%"
print*, zzzzz1, "element betwen 0.0001 ans 0.00001"
print*, real(zzzzz1)/real(percent) * 100, "%"
print*, zzzzzz1, "element betwen 0.00001 ans 0.000001"
print*, real(zzzzzz1)/real(percent) * 100, "%"
print*, zzzzzzz1, "element betwen 0.000001 ans 0.0000001"
print*, real(zzzzzzz1)/real(percent) * 100, "%"
print*, "max value in abs : ", maxi
val = percent - supone - z1-zz1-zzz1-zzzz1-zzzzz1-zzzzzz1-zzzzzzz1 - nulle
print *, val,"positive elements smaller than 0.0000001"
print *, real(val)/real(percent) * 100, "%"

end program main
