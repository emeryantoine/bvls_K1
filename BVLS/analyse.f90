program main

integer, parameter :: limit = 20, height=141970, width = 77, offset = -3
real, dimension(limit) :: sizes
real(kind=8), dimension(:), allocatable :: tmp
real :: val, expo
integer :: planet = -1
real(kind=8), dimension(16) :: try

sizes = 0
try = 0

allocate(tmp(width))

open(1, file='../../transfert/RA.out', status='old', action='read')
  do i = 1, height
    read(1, *) tmp(:)
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

print *, try(:)

close(1)
print *, "there's", nonzerol, nonzeror, "non zero element in the bottom right section of matrix A"
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
end program main
