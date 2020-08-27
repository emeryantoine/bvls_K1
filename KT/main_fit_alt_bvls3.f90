!ifort main_fit_alt_bvls3.f90 bvlsDP.f90 subf90_SQ_cov.f90  subf90_SQ.f90
!svdfit.f90 -L./ -limsl -o bvls

implicit none

    interface

    subroutine bvls ( atilde, btilde, bndtilde, soltilde, chi2, nsetp,w,  &
                itilde,loopA)
      real ( kind ( 1d0 ) ) atilde(:,:)
      real ( kind ( 1d0 ) ) btilde(:)
      real ( kind ( 1d0 ) ) bndtilde(:,:)
      real ( kind ( 1d0 ) ) soltilde(:)
      real ( kind ( 1d0 ) ) chi2
      integer nsetp
      real ( kind ( 1d0 ) ) w(:)
      integer itilde(:)
      integer loopA
    end subroutine

  end interface

integer, parameter :: ixt2=58,iwk_div=0,iwk_biais=0
integer, parameter :: npla=9
integer, parameter :: ntotal3=141976
real(kind=8), dimension(2,ixt2) :: BND
real(kind=8), dimension(ixt2) :: solnr_sol

real(kind=8), dimension(:,:), allocatable :: A2tot,A2b,tA2tot,A2org
real(kind=8), dimension(:),  allocatable :: R2tot, W2tot,R2b,W2b
character(len=20) :: string

real(kind=8), dimension(ntotal3,ixt2) :: Atot,logA2
real(kind=8), dimension(ixt2,ixt2) :: cor,Rt,tRt,tRtinv,Rinv,sigmamat
real(kind=8), dimension(ixt2,ixt2) :: P1,P2,logcov,cov
real(kind=8), dimension(ixt2) :: ww,xs,xs0,xdx,sol,logsolnr_sol,vconv
real(kind=8), dimension(ixt2) :: cditot,logdeltasolnr_sol,vqi
real(kind=8), dimension(54) :: cdi0
real(kind=8), dimension(357) :: cdi
real(kind=8), dimension(iwk_div) :: prc_mas,prc_div
real(kind=8), dimension(ntotal3) :: Wtot,Rtot,zz,datjjf2,unite
real(kind=8) :: chi2,sos,scale
integer :: i,j,k,info,rank,ik,ntot,nsetp,loopA,lim,b
integer, dimension(ixt2) :: istate

integer :: num2,tabsetf,Tnumf
integer :: start, finish

print*,'entree dans fit',ixt2,ntotal3

!lim=ixt2-iwk_div-iwk_biais

ntot=ntotal3

allocate(A2tot(ntot,ixt2))
allocate(A2org(ntot,ixt2))
allocate(tA2tot(ixt2,ntot))
allocate(R2tot(ntot))
allocate(W2tot(ntot))
allocate(A2b(ntot,ixt2))
allocate(R2b(ntot))
allocate(W2b(ntot))

A2tot(:,:)=0d0
R2tot(:)=0d0
W2tot(:)=0d0


! sert a la comparaison avec R
! on trouve les memes resultats

b=0
!open(556,file='RA.out.mexmro.12.extrait',status='old')
open(556,file='RA.out',status='old')
do i=1,ntot
read(556,*)datjjf2(i),W2tot(i),R2tot(i),(A2org(i,j),j=1,ixt2)
enddo
close(556)

open(557,file='BND.KT01',status='old')
do i=1,ixt2
read(557,*)BND(1,i),BND(2,i)
enddo
close(557)


A2tot(:,:)=A2org(:,:)
W2b=W2tot

print*,'-------------------------------------------'
print *, 'FIT SQ simple ss contrainte ====='
print*,'-------------------------------------------'

CALL system_clock(start)
call SQ(A2tot,W2b,R2tot,ntot,ixt2,solnr_sol,xdx,chi2)
call system_clock(finish)
print*, "SQ : ", finish - start

print*,'chi2',chi2

open(525,file='svd.out',status='replace')
do i=1,ixt2
write(525,*)solnr_sol(i),xdx(i)!,xdx(j)*chi2/(ntot-ixt2)
enddo
close(525)


print*,'-------------------------------------------'
print*,'ALT KT -------------------------------------------'
print*,'-------------------------------------------'

  cditot=0
  
open(458,file='cdi',status='old')
read(458,*)(cdi0(j+6),j=1,6)
read(458,*)(cdi0(j),j=1,6)
read(458,*)(cdi0(8*6+j),j=1,6)
do i=3,8
read(458,*)(cdi0((i-1)*6+j),j=1,6)
enddo


open(457,file='fort.457',status='old')
do j=1,ixt2-6*npla
read(457,*)cdi(j)
enddo
close(437)


  cditot=0
do i=1,6*npla
        cditot(i)=cdi0(i)
enddo
do i=6*npla+1,ixt2
        cditot(i)=cdi(i-6*npla)
enddo
  
print*,'avant conversion ============'

!cditot=0d0

CALL system_clock(start)
call logAtot(cditot,bnd,ixt2,vconv,vqi,npla)

do i=1,ntot
  do j=1,ixt2
        logA2(i,j)=vconv(j)*A2tot(i,j)
  enddo
!write(*,'(21(E10.2,x))')(logA2(i,j),j=1,ixt2)
enddo
chi2 = 0d0
W2b(:)=W2tot(:)
!call SQ_cov(logA2,W2b,R2tot,ntot,ixt2,logdeltasolnr_sol,logcov,chi2)
CALL system_clock(finish)
print*,"SQ_cov : ", finish-start 
print*,'chi2',chi2
!
print*,'SOL variable alternative KT01 -------------------------------------------'
do i=1,ixt2
!print*,logdeltasolnr_sol(i)+cditot(i)*vconv(i)
logsolnr_sol(i)=logdeltasolnr_sol(i)+vqi(i)
enddo

!call logAinv(logsolnr_sol,bnd,ixt2,solnr_sol,cditot,npla)

open(525,file='KT01.out',status='replace')
print*,'resultats apres conversion ==='
do i=1,ixt2
   do j=1,ixt2
    cov(i,j)=logcov(i,j)*vconv(i)*vconv(j)
   enddo
xdx(i)=sqrt(cov(i,i))
!write(525,'(5(e20.10,1x))')solnr_sol(i),xdx(i),solnr_sol(i)+cditot(i),bnd(1,i),bnd(2,i)!,xdx(j)*chi2/(ntot-ixt2)
write(525,'(6(e20.10,1x))')solnr_sol(i),solnr_sol(i)-cditot(i),cditot(i),xdx(i),bnd(1,i),bnd(2,i)!,xdx(j)*chi2/(ntot-ixt2)
!print*,solnr_sol(i),solnr_sol(i)+cditot(i)
!if(solnr_sol(i)+cditot(i).lt.bnd(1,i).or.solnr_sol(i)+cditot(i).gt.bnd(2,i))print*,'ATTENTION:',i,solnr_sol(i)+cditot(i),bnd(1,i),bnd(2,i)
if(abs(bnd(1,i)-bnd(2,i))>0d0)then
if(solnr_sol(i).lt.bnd(1,i).or.solnr_sol(i).gt.bnd(2,i))print*,'ATTENTION:',i,solnr_sol(i)+cditot(i),bnd(1,i),bnd(2,i)
endif
enddo
close(525)


print*,' BVLS -------------------------------------------'
print*,'avant wtot',ixt2,b,b*100.0/(ixt2*ntot)

      do i=1,Ntot
         R2tot(i)=DSQRT(W2tot(i))*R2tot(i)
         do j=1,ixt2
            A2tot(i,j)=DSQRT(W2tot(i))*A2tot(i,j)
         enddo
      enddo

open(557,file='BND.out',status='old')
do i=1,ixt2
read(557,*)BND(1,i),BND(2,i)
!BND(2,i)=1d0
!BND(1,i)=-1d0
enddo
!close(557)

solnr_sol(:)=0d0
CALL system_clock(start)
!call bvls(A2tot,R2tot,BND,solnr_sol,chi2,nsetp,ww,istate,loopA)
CALL system_clock(finish)
print*, "bvls : ", finish-start
print*,'apres wtot',loopA,chi2,nsetp,chi2/(ntot-ixt2)
print*,'chi2',chi2/(ntot-ixt2)
scale=chi2/(ntot-ixt2)
if(loopA>0)print*,'---------- WARNING !!!!--------',loopA

open(525,file='bvls.out',status='replace')
print*,'SOL BVLS -------------------------------------------'
do i=1,ixt2
write(525,*)solnr_sol(i),solnr_sol(i)+cditot(i)
enddo

close(525)


111   format(29(d27.20,4x))
!112   format(f4.2)
113   format(5x,f15.10,4x,i1,3x,i3,3x,i3,4x,e27.20,2x,i5)
114   format(e27.20)

End

subroutine logAtot(cdi,bnd,ixt2,vconv,vqi)
implicit none
integer :: ixt2,i,j
!real(kind=8), dimension(ntot,ixt2) :: A2tot,logA2
real(kind=8), dimension(ixt2) :: cdi,cdi0
real(kind=8), dimension(2,ixt2) :: bnd
real(kind=8), dimension(ixt2) :: vconv,vqi
real(kind=8) :: rap,qi
vconv=0d0

print*,'ds logA2'

cdi0=cdi

do i=1,ixt2
  if(bnd(1,i)==bnd(2,i))then
   vqi(i)=cdi(i)
   vconv(i)=1d0
  else
   rap=(cdi(i)-(bnd(1,i)))/(bnd(2,i)-cdi(i))
   vqi(i)=log10(rap)
   vconv(i)=(bnd(2,i)-bnd(1,i))*exp(qi)/((exp(qi)+1)**2)
  endif
enddo
cdi=cdi0

end subroutine


subroutine logAinv(qi,bnd,ixt2,pi,cdi)
implicit none
integer :: ixt2,i,j
real(kind=8), dimension(ixt2) :: cdi,qi,pi !logsolnr_sol
real(kind=8), dimension(2,ixt2) :: bnd
real(kind=8) :: rap

do i=1,ixt2
  if(bnd(1,i)==bnd(2,i))then
      pi(i)=qi(i)
  else 
rap=exp(qi(i))/((exp(qi(i))+1))
pi(i) = (bnd(1,i)) + (bnd(2,i)-bnd(1,i))*rap
  endif
enddo
end subroutine


