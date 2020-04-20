program main

integer :: start, stop

call system_clock(start)

call omphello()

call system_clock(stop)

print *, "time in tick : ", (stop - start)

end program main
