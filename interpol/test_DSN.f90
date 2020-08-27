! ifort test_DSN.f90 eop_interp.f -o testeop
implicit none

        double precision :: lod,datmjd,XINT,YINT,TINT,dat
        integer :: i,Imjd,k
        real(kind=8), dimension(6) :: tab_err
        integer, dimension(3) :: peu
   integer, parameter :: Nbrcoeff_eop=21336 !20404
   double precision, dimension(Nbrcoeff_eop) :: mjd,xp,yp,ut1_utc,dpsi,deps
   double precision, dimension(6) :: pole


open(350,file='/projets/PLANETO/data/INPOP/datasupp/eopc04.62-now',status='old')
        
  do i=1,Nbrcoeff_eop

  read(350,'(3(I4),I7,2(F11.6),2(F12.7),2(F11.6),2(F11.6),2(F11.7),2F12.6)') &
        (peu(k),k=1,3),Imjd,xp(i),yp(i),ut1_utc(i),lod,dpsi(i),deps(i),(tab_err(k),k=1,6)
        mjd(i)=Imjd
  enddo

close(350)


dat=2457463.17d0

        datmjd=dat-2400000.5d0

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

     else
        pole=0d0
endif
do i=1,3
print*,pole(i),pole(i+3)
enddo
End 
