! Ajout de l'a priori sigma dans le fit
Subroutine fit_planetes_bvls(npla,ixt2,ntotal3,PD2,BND,omc3,sig2,datjjf2,solnr_sol,LA_res,prc_mas,iwk_div,prc_div,iwk_biais &
)
use var_fit

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

integer :: ixt2,ntotal3,iwk_div,npla,iwk_biais

real(kind=8), dimension(ntotal3,ixt2), intent(in) :: PD2
real(kind=8), dimension(2,ixt2),intent(in) :: BND
real(kind=8), dimension(ntotal3), intent(in) :: sig2,omc3
real(kind=8), dimension(ixt2),intent(out) :: solnr_sol
real(kind=8), dimension(ntotal3),intent(out) :: LA_res

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
do i=1,ntotal3
   R2tot(i)=omc3(i)
enddo


DO J=1,ixt2
      DO I=1,ntotal3
      A2tot(i,j)=PD2(I,J)
      enddo
   A2tot(j+ntotal3,j)=1d0
enddo


call sig_total(npla,W2tot,sig2,ntotal3,ntot,iwk_div,prc_div)

print*,'-------------------------------------------'
! -------------------------------------------

      do i=1,Ntot
         R2tot(i)=DSQRT(W2tot(i))*R2tot(i)
         do j=1,ixt2
            A2tot(i,j)=DSQRT(W2tot(i))*A2tot(i,j)
         enddo
      enddo


! sert a la comparaison avec R
! on trouve les memes resultats

string='(   (e30.20,3x))'
write(string(2:4),'(I3)')ixt2+3

open(556,file='RAW.out',status='replace')
do i=1,ntotal3
write(556,string)datjjf2(i),W2tot(i),R2tot(i),(A2tot(i,j),j=1,ixt2)
enddo
close(556)

open(557,file='BND.out',status='replace')
do i=1,ixt2
write(557,*)BND(1,i),BND(2,i)
enddo
close(557)

print*,'avant wtot',ixt2
solnr_sol(:)=0d0
call bvls(A2tot,R2tot,BND,solnr_sol,chi2,nsetp,ww,istate,loopA)
print*,'apres wtot',loopA,chi2,nsetp,chi2/(ntot-ixt2)
print*,'chi2',chi2/(ntot-ixt2)
if(loopA>0)print*,'---------- WARNING !!!!--------',loopA


!print*,'apres wtot',iwk_div,lim,ixt2
! ================================
! ajustement par LS simple 

! ================================

open(456,file='resultats/fort.440.sig',status='replace')
print*,' solution LS simple avec sigma ================ > resultats/fort.440.sig'
do j=1,iwk_div
write(456,*)prc_mas(j)+solnr_sol(j+lim),solnr_sol(j+lim),xdx(j+lim)
enddo
close(456)
111   format(29(d27.20,4x))
!112   format(f4.2)
113   format(5x,f15.10,4x,i1,3x,i3,3x,i3,4x,e27.20,2x,i5)
114   format(e27.20)

End subroutine fit_planetes_bvls
