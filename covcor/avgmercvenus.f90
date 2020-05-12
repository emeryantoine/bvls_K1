program main

  integer, parameter :: height = 26582, width = 415
  real(kind=8), dimension(width)  :: tmp
  real(kind=8) :: small=0, smallavg=0, bigavg=0, big=0

  open(1, file="../../transfert/cas_complet/RA.out", status="old", action="read")
  do i = 1, height
    read(1,*) tmp(:)
    do j = 62, width
     ! if(abs(tmp(j)) > 1) then
     !   big = big + 1
     !   bigavg = bigavg + tmp(j)
     ! else
     !   small = small + 1
     !   smallavg = smallavg + tmp(j)
     ! endif
      small = small + tmp(j)   
    end do
  end do
  !smallavg = smallavg/small
  !bigavg = bigavg/big
        
 ! print*, "big : ", bigavg
 ! print*, "small : ", smallavg
  print*, small/(height*(width-62))

end program main
