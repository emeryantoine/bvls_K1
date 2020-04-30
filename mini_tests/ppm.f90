program main


  open(1, file='./img.ppm', status="new", action="write")
  write(1, '(a)')"P6"
  write(1, *)"5  5"
  write(1,*)"255"
  do i=1, 25
    write(1, *) i*10
  end do
  close(1)

end program main
