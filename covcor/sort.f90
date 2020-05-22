program sort

integer, parameter :: height = 142388, width = 415
real(kind=8), dimension(width) :: lu
integer, dimension(width-3) :: zero
integer, dimension(height) :: id, val
integer :: tmpval, tmpid, prout
intrinsic :: FINDLOC
integer :: x(1)
real(kind=8), dimension(height, width) :: mat
real(kind=8), dimension(2,height) :: bndsort
real(kind=8), dimension(2,1) :: lubnd

open(1, file="../../transfert/cas_complet/040520/RAW.out.412", status="old", action = "read")
do i = 1, height
  read(1,*) lu
  id(i) = i
  val(i) = 0
  do j = 4, width
    if(abs(lu(j)) < 5e-3) val(i) = val(i) + 1
  enddo
end do
close(1)

do i = 1, height
  do j = 1, height - 1
    if(val(j) > val(j+1)) then
      tmpval = val(j)
      tmpid = id(j)
      val(j) = val(j+1)
      id(j) = id(j+1)
      val(j+1) = tmpval
      id(j+1) = tmpid
    endif
  end do
end do

do i = 1, height
  print*, val(i)
enddo

open(1, file="../../transfert/cas_complet/040520/RAW.out.412", status="old", action="read")
do i = 1, height
  read(1,*) lu(:)
  x = FINDLOC(id, value=i)
  mat(x(1), :) = lu(:)
end do
close(1)

open(1, file="../../transfert/cas_complet/040520/BND.out", status="old",&
action="write")
do i =1, height
  read(1,*) lubnd(:,1)
  x = findloc(id, value=i)
  bndsort(:, x(1)) = lubnd(1,:)
end do
close(1)

open(1, file="../../transfert/cas_complet/040520/BND.out.sort")
do i = 1, height
  write(1,*) bndsort(:,i)
end do
close(1)

status = system("rm -f total.ppm")
open(1, file="./total.ppm", status="new", action="write")
write(1,'(a)') "P2"
write(1,'(a)') "412 1420"
write(1,'(a)') "100"
do i = 1, 142000, 100
  zero(:) = 0
  do j = 4, width
    do k = 0, 99
      if(abs(mat(i+k,j)) < 5e-3) zero(j-3) = zero(j-3) + 1
    end do
  end do
  write(1, *) zero(:)
end do
close(1)

do i = 1, height
  tmpval = 0
  do j = 4, width
    if(abs(mat(i,j)) < 5e-3) tmpval = tmpval + 1
  end do
  print*, tmpval
end do

prout = system("rm -f ./RAW.out.412.sort")
open(1, file="../../transfert/cas_complet/040520/RAW.out.412.sort", status ="new", action="write")
do i = 1, height
  write(1, *) mat(i,:)
end do
close(1)



end program sort
