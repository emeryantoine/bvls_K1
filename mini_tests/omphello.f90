subroutine omphello()
  implicit none

  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS, a, i

!$OMP PARALLEL
      print *, "hello world !", OMP_GET_THREAD_NUM(),"/", OMP_GET_NUM_THREADS()

      !$OMP DO      
      do i=1, 1000000
      a = a +1
      end do
      !$OMP END DO
      
      print *, a
      
      
!$OMP END PARALLEL
      print *, OMP_GET_MAX_THREADS()

end subroutine
