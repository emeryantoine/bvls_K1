program test
include 'mpif.h'

integer :: rank, world

call MPI_INIT()
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, world)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank)
print*, "hello world", rank, world
CALL MPI_FINALIZE()

end program test
