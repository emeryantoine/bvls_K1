program test

integer :: rank, world

CALL MPI_COM_SIZE(MPI_COM_WORLD, world)
call MPI_COM_RANK(MPI_COM_WORLD, rank)
call MPI_INIT()
print*, "hello world", rank, world
CALL MPI_FINALIZE()

end program test
