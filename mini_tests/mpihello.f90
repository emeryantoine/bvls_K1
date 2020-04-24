program main

  use mpi
  integer :: id, nbr, err

  call MPI_INIT(err)
  call MPI_COMM_RANK(MPI_COMM_WORLD, id,  err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbr, err)
  !$OMP PARALLEL
    print*, "Hello world !", id, nbr
    call sleep(5)
  !$OMP END PARALLEL
  call MPI_FINALIZE(err)

end program main
