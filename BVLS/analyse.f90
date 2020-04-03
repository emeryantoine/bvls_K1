program main

integer, parameter :: limit = 20, height=141970, width = 77, offset = -3
real, dimension(limit) :: sizes
real(kind=8), dimension(:), allocatable :: tmp
real :: val, expo

sizes = 0

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
  end do

do i = 1, limit
  print *, sizes(i), "elements are < 10^", -i, " --> ", (sizes(i)/(height*width))*100, "%"
  end do
end program main
