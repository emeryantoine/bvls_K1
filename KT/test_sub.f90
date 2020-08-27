program main

integer, parameter :: height = 20, width = 5
real(kind=8) :: chi2
real(kind=8), dimension(height, width) :: A, logA
real(kind=8), dimension(height) :: W, R
real(kind=8), dimension(width) :: sol, xdx, logsol
real(kind=8), dimension(width, width) :: logcov
real(kind=8) :: tmp1
real(kind=8), dimension(100) :: tmp2

open(1, file='RA.out', status='old')
do i=1, height
  read(1,*) tmp1, W(i), R(i), (A(i, j),j=1, width), tmp2(:)
end do

call SQ(A, W, R, height, width, sol, xdx, chi2)
print*, chi2

logA(:,:) = log10(A(:,:))

call SQ_cov(logA, W, R, height, width, logsol, logcov, chi2)
print*, chi2

end program
