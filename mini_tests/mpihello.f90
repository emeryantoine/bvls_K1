program mpi
  
  include "mpif.h"

  integer :: id, nbr, err

  call MPI_INIT(err)
  call MPI_COMM_SIZE(MPI_COM_WORLD, nbr, err)
  call MPI_COM_RANK(MPI_COM_WORLD, id, err)
  print*, "Hello world !", id, nbr
  call MPI_FINALIZE(err)

end program mpi
