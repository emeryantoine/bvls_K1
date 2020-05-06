! Ajout de l'a priori sigma dans le fit
implicit none

    interface

    subroutine bvls ( atilde, btilde, bndtilde, soltilde, chi2, nsetp, w, itilde,  loopA)
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

integer, parameter :: ixt2=412,iwk_div=347,iwk_biais=11
integer, parameter :: ntotal3=141976
!integer :: ixt2,ntotal3,iwk_div,npla,iwk_biais

!real(kind=8), dimension(ntotal3,ixt2), intent(in) :: PD2
real(kind=8), dimension(2,ixt2) :: BND
!real(kind=8), dimension(ntotal3), intent(in) :: sig2,omc3
real(kind=8), dimension(ixt2) :: solnr_sol
!real(kind=8), dimension(ntotal3),intent(out) :: LA_res

real(kind=8), dimension(:,:), allocatable :: A2tot,A2b
real(kind=8), dimension(:),  allocatable :: R2tot, W2tot,R2b,W2b
character(len=20) :: string

real(kind=8), dimension(ntotal3,ixt2) :: Atot
real(kind=8), dimension(ixt2,ixt2) :: cor
real(kind=8), dimension(ixt2) :: ww,xs,xs0,xdx,sol
real(kind=8), dimension(iwk_div) :: prc_mas,prc_div
real(kind=8), dimension(ntotal3) :: Wtot,Rtot,zz,datjjf2,unite
real(kind=8) :: chi2,sos
integer :: i,j,k,info,rank,ik,ntot,nsetp,loopA,lim
integer, dimension(ixt2) :: istate
real :: reaad, min1=10e-3, min2=10e-30
integer :: num2,tabsetf,Tnumf
print*,'entree dans fit',ixt2,ntotal3

lim=ixt2-iwk_div-iwk_biais

ntot=ntotal3+ixt2

allocate(A2tot(ntot,ixt2))
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

open(556,file='../../transfert/cas_complet/040520/RAW.out.412',status='old')
do i=1,ntot
read(556,*)datjjf2(i),W2tot(i),R2tot(i),(A2tot(i,j),j=1,ixt2)
enddo
close(556)

if(.true.) then
  do i = 1, 74
    do j = 1,ixt2
      if (A2tot(i,j) < min1) A2tot(i,j) = 0d0
    end do
  end do
endif
if(.false.) then
  do i = 75,ntot
    do j = 1,ixt2
      if(A2tot(i,j) < min2) A2tot(i,j) = 0d0
    end do
  end do
end if

open(557,file='../../transfert/cas_complet/040520/BND.out',status='old')
do i=1,ixt2
read(557,*)BND(1,i),BND(2,i)
enddo
close(557)

print*,'avant wtot',ixt2
solnr_sol(:)=0d0
call bvls(A2tot,R2tot,BND,solnr_sol,chi2,nsetp,ww,istate,loopA)
!solnr_sol=-solnr_sol
print*,'apres wtot',loopA,chi2,nsetp,chi2/(ntot-ixt2)
print*,'chi2',chi2/(ntot-ixt2)
if(loopA>0)print*,'---------- WARNING !!!!--------',loopA


!open(556,file='RA.out.postfit',status='replace')
!do i=1,ntotal3
!write(556,'(6(E30.20,x))')(A2tot(i,j)/W2tot(i),j=1,6)
!enddo
!close(556)

print*,'SOL BVLS-------------------------------------------'
OPEN(42, file="reference.out", status="old")
print*, "result, reference, res/ref"
do j=1,ixt2
  read(42,*) reaad
  print*,solnr_sol(j), reaad, solnr_sol(j)/reaad
enddo

111   format(29(d27.20,4x))
!112   format(f4.2)
113   format(5x,f15.10,4x,i1,3x,i3,3x,i3,4x,e27.20,2x,i5)
114   format(e27.20)

End
