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

integer, parameter :: iwk_div=347, iwk_biais=11, ixt2=412 !iwk_div=347,iwk_biais=11, ixt2=412
integer, parameter :: ntotal3 = 141976 !ntotal3=141976
!integer :: ixt2,ntotal3,iwk_div,npla,iwk_biais

!real(kind=8), dimension(ntotal3,ixt2), intent(in) :: PD2
real(kind=8), dimension(2,ixt2) :: BND
  !real(kind=8), dimension(ntotal3), intent(in) :: sig2,omc3
real(kind=8), dimension(ixt2) :: solnr_sol
!real(kind=8), dimension(ntotal3),intent(out) :: LA_res
real(kind=8) :: read_sol

real(kind=8), dimension(:,:), allocatable :: A2tot,A2b
real(kind=8), dimension(:),  allocatable :: R2tot, W2tot,R2b,W2b
character(len=20) :: string

real(kind=8), dimension(ntotal3,ixt2) :: Atot
real(kind=8), dimension(ixt2,ixt2) :: cor
real(kind=8), dimension(ixt2) :: ww,xs,xs0,xdx,sol
real(kind=8), dimension(iwk_div) :: prc_mas,prc_div
real(kind=8), dimension(ntotal3) :: Wtot,Rtot,zz,datjjf2,unite
real(kind=8) :: chi2,sos
integer :: i,j,k,info,rank,ik,ntot,nsetp,loopA,lim, m, n
integer, dimension(ixt2) :: istate
real :: reaad, min1=5e-3, small = 1d0, val = 0
integer :: num2,tabsetf,Tnumf, zero = 0, summ = 0
integer, dimension(10) :: inf
integer, dimension(ixt2) :: nonnul
!print*,'entree dans fit',ixt2,ntotal3
integer, dimension(:), allocatable :: test
real(kind=8), dimension(:,:), allocatable :: R, Rt, At, RtR, AtA

lim=ixt2-iwk_div-iwk_biais

ntot=ntotal3+ixt2

allocate(test(ntot))
allocate(A2tot(ntot,ixt2))
allocate(At(ixt2,ntot))
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

!open(556,file='../covcor/RAW.out.412.sort',status='old')
!open(556, file="./RAW.out.412", status="old")
open(556, file="../../transfert/cas_complet/040520/RAW.out.412", status="old")
!open(556, file="../../transfert/RAW.out.412.sort", status="old")
do i=1,ntot
read(556,*)datjjf2(i),W2tot(i),R2tot(i),(A2tot(i,j),j=1,ixt2)
enddo
close(556)

!if(.true.) then
!  do i = 1, 74
!    do j = 1,ixt2
!      if (A2tot(i,j) < min1) A2tot(i,j) = 0d0
!    end do
!  end do
!endif
!if(.false.) then
!  do i = 75,ntot
!    do j = 1,ixt2
!      if(A2tot(i,j) < min2) A2tot(i,j) = 0d0
!      if (A2tot(i,j) /= 0) then
!        if(abs(A2tot(i,j)) < small) small = abs(A2tot(i,j))
!      end if
!    end do
!  end do
!end if


if(.true.) then
  do i = 1, ntot
    do j = 1, ixt2
      if(abs(A2tot(i,j)) < min1) then
      !if(.true.) then
        A2tot(i,j) = 0d0
        zero = zero + 1
      end if
      !if(A2tot(i,j) .eq. 0) zero = zero + 1
      !if(j .eq. 58) then
      !  A2tot(i,j) = 0d0
      !endif
    end do
  end do
endif

!print*, "zero :",zero,"/",ntot*ixt2, ":", 100*(real(zero)/real(ntot*ixt2)) ,"%"

nonnul(:) = 0
do j = 1, ixt2
  do i = 1,ntot
    if(A2tot(i,j) /= 0) nonnul(j) = nonnul(j) + 1
  end do
end do

!print*, "colonnes nulles ?:"
!do i = 1, ixt2
!  print*, i, nonnul(i)
!enddo

!summ = SUM(nonnul)
!print*, "comparaison",summ, zero, summ+zero, ntot*ixt2

!print*, "smallest number", small
open(557,file='../../transfert/cas_complet/040520/BND.out',status='old')
!open(557,file='../../transfert/cas_complet/040520/BND.out',status='old')
do i=1,ixt2
  read(557,*)BND(1,i),BND(2,i)
enddo
close(557)

if(.false.) then
  BND(:,:) = BND(:,:) * 100
endif

!print*,'avant wtot',ixt2
solnr_sol(:)=0d0
call bvls(A2tot,R2tot,BND,solnr_sol,chi2,nsetp,ww,istate,loopA)
!print *, A2tot(1,:)
!solnr_sol=-solnr_sol
!print*,'apres wtot',loopA,chi2,nsetp,chi2/(ntot-ixt2)
!print*,'chi2',chi2/(ntot-ixt2)
if(loopA>0)print*,'---------- WARNING !!!!--------',loopA

allocate(R(nsetp, nsetp))
allocate(Rt(nsetp,nsetp))
allocate(RtR(nsetp, nsetp))
allocate(AtA(nsetp,nsetp))
print*, "alloc done"

print*, "NSETP = ", nsetp

do i = 1, nsetp
  do j = 1, nsetp
    R(i,j) = A2tot(i,j)
  enddo
enddo

!print*, "R tilde fait"
!m = size(Rt, 1)
!n = size(Rt, 2)
!print*, "Rt", m, n
!m = size(R, 1)
!n = size(R, 2)
!print*, "R", m, n
!
!m = size(A2tot, 1)
!n = size(A2tot, 2)
!print*, "A", m, n
!m = size(At, 1)
!n = size(At, 2)
!print*, "At", m, n


!Rt = transpose(R)
!At = transpose(A2tot)
!print*, "transposee faites"

!RtR = matmul(Rt, R)
!AtA = matmul(At, A2tot)
!print*, "matmul fait"

!do i = 1, nsetp
!  do j = 1, nsetp
!    print*, RtR(i,j), AtA(i,j)
!  enddo
!enddo

if(.false.) then
  test(:) = 0
  do i = 1,ntot
    do j = 1, ixt2
      if(A2tot(i,j) /= 0) test(i) = test(i) + 1
    end do
  end do
endif


!open(556,file='RA.out.postfit',status='replace')
!do i=1,ntotal3
!write(556,'(6(E30.20,x))')(A2tot(i,j)/W2tot(i),j=1,6)
!enddo
!close(556)


!open(123, file="./fast.txt", status="old", ation="read")
if(.true.) then
  !print*,'SOL BVLS-------------------------------------------'
  !OPEN(42, file="ref.out", status="old")
  !print*, "result, reference, res/ref"
  open(420, file='./output.txt', status='old')
  do j=1,ixt2
    !read(42,*) reaad
    !print*, solnr_sol(j), reaad, real(real(solnr_sol(j))/real(reaad))
    READ(420,*) read_sol
    print*, read_sol, solnr_sol(j), read_sol-solnr_sol(j)
    !write(420, '(6(E30.20,x))') solnr_sol(j)
  enddo
  close(42)
endif


111   format(29(d27.20,4x))
!112   format(f4.2)
113   format(5x,f15.10,4x,i1,3x,i3,3x,i3,4x,e27.20,2x,i5)
114   format(e27.20)

End
!open(314, file="./outputA.out", status="new", action="write")
