Subroutine select_eop(dat)
	use var_obs

	implicit none
	integer :: i
	double precision, intent(in) :: dat
	double precision :: XINT,YINT,TINT,datmjd

	datmjd=dat-2400000.5d0
!print*,'date eop',Nbrcoeff_eop
!print*,datmjd,mjd(1),mjd(Nbrcoeff_eop)

        if(datmjd>mjd(1).and.datmjd<mjd(Nbrcoeff_eop))then

	call INTERP_EOP(mjd,Xp,Yp,ut1_utc,Nbrcoeff_eop,datmjd,XINT,YINT,TINT)	

do i=1,Nbrcoeff_eop
	if(abs(datmjd-mjd(i))<1d0.and.(datmjd-mjd(i))>0d0)then
	  pole(4)=(dpsi(i)*(1d0-(datmjd-mjd(i)))+dpsi(i+1)*(1d0-(datmjd-mjd(i+1))))/2d3
	  pole(5)=(deps(i)*(1d0-(datmjd-mjd(i)))+deps(i+1)*(1d0-(datmjd-mjd(i+1))))/2d3
	endif
enddo

	pole(1)=XINT
	pole(2)=YINT
	pole(3)=TINT

!pole=0d0
!plate=0d0

!pole(1)=0.0841220994d0
!pole(2)=0.4894010471d0
!pole(3)=32d0-31.971890d0
!pole(4)=-0.050576d0
!pole(5)=-0.0048915318d0

        else

           pole=0d0
        
        endif

!write(*,'(a4,5(f20.10,1x))')'Pole',(Pole(i),i=1,5)
	return
End Subroutine select_eop
