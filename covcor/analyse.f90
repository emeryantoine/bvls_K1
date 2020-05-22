PROGRAM MAIN

  integer, parameter :: width=415
  real(kind=8), dimension(6,width) :: tmp
  integer, dimension(width-3) :: zeros, neg, smalls
  real(kind=8), dimension(width) :: tmp_RA
  integer :: change = 0, status, counter=0

status = system("rm *.ppm")

open(3, file="./RAW.out.412", status='old', action='read')
open(2, file="./1.ppm", status="new", action="write")
write(2,'(a)')"P2"
write(2,'(a)')"412 6253"
write(2,'(a)')"100"

do i = 1, 25012/4
  do k = 1, 4
    read(3, *) tmp(k,:)
    counter = counter + 1
  enddo
  do k = 2, 4
    tmp(1,:) = tmp(1, :) + tmp(k, :)
  end do

  do j = 4, width
    if(abs(tmp(1,j)) > 4*(5e-3)) then
      write(2, *) 0
    else
      write(2, *) 100
    endif
  end do
end do
close(2)

do i = 1, 1
  read(3,*) tmp(1,:)
  counter = counter + 1
end do

print*, counter

open(2, file="./2.ppm", status="new", action="write")
write(2, '(a)')"P2"
write(2,'(a)')"412 1570"
write(2,'(a)')"100"

do i = 25014, 26582
  read(3,*)tmp(1,:)
  
    counter = counter + 1
  do j = 4, width
    if(abs(tmp(1,j)) > 5e-3) then
      write(2,*) 0
    else
      write(2,*) 100
    endif
  end do
end do
close(2)

print*, counter
open(2, file="./3.ppm", status="new", action="write")
write(2,'(a)')"P2"
write(2,'(a)')"412 7127"
write(2,'(a)')"100"

do i = 1, 42762/6
  do k = 1, 6
    read(3,*)tmp(k,:)
    counter = counter + 1
  end do
  do k = 2, 6
    tmp(1,:) = tmp(1,:) + tmp(k,:)
  end do
  do j = 4, width
    if(abs(tmp(1,j)) > 6*(5e-3)) then
      write(2,*) 0
    else
      write(2,*) 100
    endif
  end do
end do
close(2)
do i = 1, 3
  read(3,*) tmp(1,:)
  counter = counter + 1
end do
print*, counter
open(2, file="./4.ppm", status="new", action="write")
write(2, '(a)')"P2"
write(2, '(a)')"412 2539"
write(2, '(a)')"100"

do i = 1, 2539
  do k = 1, 5
    read(3,*) tmp(k, :)
    counter = counter + 1
  end do
  do k = 2, 5
    tmp(1,:) = tmp(1,:) + tmp(k, :)
  end do

  do j = 4, width
    if(abs(tmp(1,j)) > 5*(5e-3)) then
      write(2,*) 0
    else
      write(2,*) 100
    endif
  end do
end do
close(2)

read(3,*) tmp(1,:)
print*, counter
open(2, file="./5.ppm", status="new", action="write")
write(2,'(a)') "P2"
write(2,'(a)') "412 4082"
write(2,'(a)') "100"

do i = 1, 4082
  do k = 1, 4
    read(3,*) tmp(k,:)
    counter = counter + 1
  end do
  do k = 2, 4
    tmp(1,:) = tmp(1,:) + tmp(k,:)
  end do

  do j = 4, width
    if(abs(tmp(1,j)) > 4*(5e-3)) then
      write(2,*) 0
    else
      write(2,*) 100
    endif
  end do
end do
close(2)

print*, counter
open(2, file="./6.ppm", status="new", action="write")
write(2,'(a)') "P2"
write(2,'(a)') "412 6550"
write(2,'(a)') "100"

do i = 1, 6550
  do k = 1, 4
    read(3,*) tmp(k,:)
    counter = counter + 1
  end do
  do k = 2, 4
    tmp(1,:) = tmp(1,:) + tmp(k,:)
  end do

  do j = 4, width
    if(abs(tmp(1,j)) > 4*(5e-3)) then
      write(2,*) 0
    else
      write(2,*) 100
    endif
  end do
end do
close(2)

print*, counter
open(2, file="./7.ppm", status="new", action="write")
write(2, '(a)') "P2"
write(2,'(a)') "412 2826"
write(2,'(a)') "100"

do i = 1, 2826
  do k = 1, 4
    read(3,*) tmp(k,:)
    counter = counter + 1
  end do
  do k = 2, 4
    tmp(1,:) = tmp(1,:) + tmp(k,:)
  end do

  do j = 4, width
    if(abs(tmp(1,j)) > 4*(3e-5)) then
      write(2,*) 0
    else
      write(2,*) 100
    endif
  end do
end do
close(2)

print*, counter
open(2,file="./8.ppm", status="new", action="write")
write(2,'(a)') "P2"
write(2, '(a)') "412 1524"
write(2,'(a)') "100"

do i = 1, 1524
  do k = 1, 4
    read(3,*) tmp(k,:)
    counter = counter + 1
  end do
  do k = 2, 4
    tmp(1,:) = tmp(1,:) + tmp(k,:)
  end do

  do j = 4, width
    if(abs(tmp(1,j)) > 4*(5e-3)) then
      write(2,*) 0
    else
      write(2,*) 100
    endif
  end do
end do
close(2)
close(3)

print*, counter

END PROGRAM MAIN
