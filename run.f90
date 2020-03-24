PROGRAM main

INTERFACE
   SUBROUTINE BVLS (A, B, BND, X, RNORM, NSETP, W, INDEX, IERR)
    REAL(KIND(1E0)) A(:,:), B(:), BND(:,:), X(:), RNORM, W(:)
    INTEGER NSETP, INDEX(:), IERR
   END SUBROUTINE
END INTERFACE

logical :: debug = .FALSE.
integer, parameter :: width = 74
integer, parameter :: height = 141970

!todo : A, B
!raw -> 74 x 141970 
real :: tmp1, tmp2

real :: RNORM, tmp3
integer :: NSETP,IERR
real :: BND(2,width)
real, dimension(width, height) :: A
real, dimension(height) :: B, X, tmp4
real, dimension(width) :: W
integer, dimension(width) :: INDEX
integer :: err

open(1, file='BND.out.20', status='old')

  do i = 1,width
    read (1, *) BND(1, i), BND(2,i)
  end do

close(1)

if (debug .eqv. .TRUE.) then
  do i = 1, width
    print *, BND(1, i), BND(2,i)
    exit
  end do
end if

open(2, file='RAW.out.20', status='old', iostat=err)
if (err > 0) then
  print *, "an error has occur trying to oen RAW", err
end if

if (debug .eqv. .TRUE.) then
  print *, err
end if

do i=1,height
  if (debug .eqv. .TRUE.) then
    print *, "step : ", i
    read(1, *, iostat=err) tmp1, tmp2, B(i), A(:,i)
  end if
end do
close(2)

!lire les fichiers RAW et BND pour pouvoir lancer la suite
call bvls(A, B, BND, X, RNORM, NSETP, W, INDEX, IEER) 

END PROGRAM main
