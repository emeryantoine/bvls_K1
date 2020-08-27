 
!##############################################################################
 
!==============================================================================!
! Calcul de la SVD de la matrice a(1:m,1:n) : A = U.W.Vt                       !
! (issu de Numerical Recipes in Fortran ; traduit en F90)                      !
!                                                                              !
! Entree : a              matrice dont on veut la SVD                          !
!          m, n, mp, np   dimensions physiques des matrices                    !
!                                                                              !
! Sortie : a, v, w        en temps que U, V et W                               !
!                                                                              !
!==============================================================================!
! Librairie EPROC :                                                            !
!    Ephemerides, Predictions et Reductions pour les Corps Celestes.           !
! Conditions d'utilisation : lire LICENSE                                      !
! (C) Copr. 1986-92 Numerical Recipes Software.                                !
! Copyright (C) 1998, Moloko-X.                                                !
!==============================================================================!
 
!##############################################################################
 
subroutine svdcmp (a, m, n, mp, np, w, v)
 
   implicit none
   integer, intent(in) :: m, mp, n, np
   real (kind=8), dimension(mp,np), intent(inout) :: a
   real (kind=8), dimension(np,np), intent(out) :: v
   real (kind=8), dimension(np), intent(out) :: w
 
   integer, parameter :: nmax=500
   integer :: i, its, j, jj, k, l, nm,b
   real (kind=8) :: anorm, c, f, g, h, s
   real (kind=8) :: scale, x, y, z, pythag
   real (kind=8), dimension(nmax) :: rv1
 
   g = 0.d0
   scale = 0.d0
   anorm = 0.d0
   do i = 1, n
      l = i+1
      rv1(i) = scale*g
      g = 0.d0
      s = 0.d0
      scale = 0.d0
      if (i <= m) then
         do k = i, m
            scale = scale + abs(a(k,i))
         end do
         if (scale /= 0.d0) then
            do k = i, m
               a(k,i) = a(k,i)/scale
               s = s + a(k,i)*a(k,i)
            end do
            f = a(i,i)
            g = -sign(sqrt(s),f)
            h = f*g-s
            a(i,i) = f-g
            do j = l, n
               s = 0.d0
               do k = i, m
                  s = s + a(k,i)*a(k,j)
               end do
               f = s/h
               do k = i, m
                  a(k,j) = a(k,j)+f*a(k,i)
               end do
            end do
            do k = i, m
              a(k,i) = scale*a(k,i)
            end do
         end if
      end if
      w(i) = scale *g
      g = 0.d0
      s = 0.d0
      scale = 0.d0
      if ((i <= m) .and. (i /= n)) then
         do k = l, n
            scale = scale + abs(a(i,k))
         end do
         if (scale /= 0.d0) then
            do k = l, n
               a(i,k) = a(i,k)/scale
               s = s + a(i,k)*a(i,k)
            end do
            f = a(i,l)
            g = -sign(sqrt(s),f)
            h = f*g-s
            a(i,l) = f-g
            do k = l, n
               rv1(k) = a(i,k)/h
            end do
            do j = l, m
               s = 0.d0
               do k = l, n
                  s = s + a(j,k)*a(i,k)
               end do
               do k = l, n
                  a(j,k) = a(j,k)+s*rv1(k)
               end do
            end do
            do k = l, n
               a(i,k) = scale*a(i,k)
            end do
         end if
      end if
      anorm = max(anorm,(abs(w(i))+abs(rv1(i))))
   end do

   do i = n, 1, -1
      if (i < n) then
         if (g /= 0.d0) then
            do j = l, n
               v(j,i) = (a(i,j)/a(i,l))/g
            end do
            do j = l, n
               s = 0.d0
               do k = l, n
                  s = s + a(i,k)*v(k,j)
               end do
               do k = l, n
                  v(k,j) = v(k,j)+s*v(k,i)
               end do
            end do
         end if
         do j = l, n
            v(i,j) = 0.d0
            v(j,i) = 0.d0
         end do
      end if
      v(i,i) = 1.d0
      g = rv1(i)
      l = i
   end do

   do i = min(m,n), 1, -1
      l = i+1
      g = w(i)
      do j = l, n
         a(i,j) = 0.d0
      end do
      if (g /= 0.d0) then
         g = 1.d0/g
         do j = l, n
            s = 0.d0
            do k = l, m
               s = s + a(k,i)*a(k,j)
            end do
            f = (s/a(i,i))*g
            do k = i, m
               a(k,j) = a(k,j)+f*a(k,i)
            end do
         end do
         do j = i, m
            a(j,i) = a(j,i)*g
         end do
      else
         do j= i, m
            a(j,i) = 0.d0
         end do
      endif
      a(i,i) = a(i,i)+1.d0
   end do

   do k = n, 1, -1
      do its = 1, 30
         do l = k, 1, -1
            nm = l-1
            if ((abs(rv1(l))+anorm) == anorm)  goto 2
            if ((abs(w(nm))+anorm) == anorm)  goto 1
         end do
1        c = 0.d0
         s = 1.d0
         do i = l, k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            if ((abs(f)+anorm) == anorm) goto 2
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            h = 1.d0/h
            c =  (g*h)
            s = -(f*h)
            do j = 1, m
              y = a(j,nm)
              z = a(j,i)
              a(j,nm) = (y*c)+(z*s)
              a(j,i) = -(y*s)+(z*c)
            end do
         end do
2        z = w(k)
         if (l == k) then
            if (z < 0.d0) then
               w(k) = -z
               do j = 1, n
                  v(j,k) = -v(j,k)
               end do
            end if
            goto 3
         end if
         if (its == 30) stop '*** no convergence in svdcmp'
         x = w(l)
         nm = k-1
         y = w(nm)
         g = rv1(nm)
         h = rv1(k)
         f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.d0*h*y)
         g = pythag(f,1.d0)
         f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
         c = 1.d0
         s = 1.d0
         do j = l, nm
            i = j+1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g
            z = pythag(f,h)
            rv1(j) = z
            c = f/z
            s = h/z
            f =  (x*c)+(g*s)
            g = -(x*s)+(g*c)
            h = y*s
            y = y*c
            do jj = 1, n
               x = v(jj,j)
               z = v(jj,i)
               v(jj,j) =  (x*c)+(z*s)
               v(jj,i) = -(x*s)+(z*c)
            end do
            z = pythag(f,h)
            w(j) = z
            if (z /= 0.d0) then
              z = 1.d0/z
              c = f*z
              s = h*z
            end if
            f =  (c*g)+(s*y)
            x = -(s*g)+(c*y)
            do jj = 1, m
              y = a(jj,j)
              z = a(jj,i)
              a(jj,j) =  (y*c)+(z*s)
              a(jj,i) = -(y*s)+(z*c)
            end do
         end do
         rv1(l) = 0.d0
         rv1(k) = f
         w(k) = x
      end do
3     continue
   end do
 
   return
end subroutine svdcmp

!=== FIN ======================================================================!
 
!##############################################################################
 
!==============================================================================!
! Resolution de A.X = B, A donnee par les matrices u, v, w                     !
! (issu de Numerical Recipes in Fortran ; traduit en F90)                      !
!                                                                              !
! Entree : u, v, w        matrices fournies par svdcmp                         !
!          m, n, mp, np   dimensions physiques des matrices                    !
!          b              second membre de l'equation                          !
!                                                                              !
! Sortie : x              vecteur solution                                     !
!                                                                              !
!==============================================================================!
! Librairie EPROC :                                                            !
!    Ephemerides, Predictions et Reductions pour les Corps Celestes.           !
! Conditions d'utilisation : lire LICENSE                                      !
! (C) Copr. 1986-92 Numerical Recipes Software.                                !
! Copyright (C) 1998, Moloko-X.                                                !
!==============================================================================!
 
!##############################################################################
 
subroutine svbksb (u, w, v, m, n, mp, np, b, x)
 
   implicit none
   integer, intent(in) :: m, mp, n, np
   real (kind=8), dimension(mp,np), intent(in) :: u
   real (kind=8), dimension(np,np), intent(in) :: v
   real (kind=8), dimension(np), intent(in) :: w
   real (kind=8), dimension(mp), intent(in) :: b
   real (kind=8), dimension(np), intent(out) :: x
 
   integer :: i, j, jj
   integer, parameter :: nmax=500
   real (kind=8) :: s
   real (kind=8), dimension(nmax) :: tmp

   do j = 1, n
      s = 0.d0
      if (w(j) /= 0.d0) then
         do i = 1, m
           s = s + u(i,j)*b(i)
         end do
         s = s/w(j)
      end if
      tmp(j) = s
   end do
   do j = 1, n
     s = 0.d0
     do jj = 1, n
       s = s + v(j,jj)*tmp(jj)
     end do
     x(j) = s
   end do

   return
end subroutine svbksb

!=== FIN ======================================================================!
subroutine svdfit1501(y,ind,TOL,sig,ndata,a,ma,u,v,w,mp,np,xdx,cor,yy,fd,chisq)
  implicit none
      INTEGER :: ma,mp,ndata,np,i,k,j,bw,ind,ib
      real(kind=8) :: chisq,xds,ss,sum,thresh,tmp,wmax,TOL
      real(kind=8) , dimension(ma) :: a,xdx
      real(kind=8) , dimension(ndata) :: sig,y,yy
      real(kind=8) , dimension(ma,ma) :: cor,impc,cov2
      real(kind=8) , dimension(np) :: w,wti
      real(kind=8) , dimension(np,np) :: v
      real(kind=8),dimension(mp,np) :: u, fd
      real(kind=8), dimension(mp) :: b(mp)
      integer, dimension(ma) :: ktst
      integer, parameter :: NMAX=300000,MMAX=5000
!      real(kind=8),parameter :: TOL=5.d-17  ! d-5 au depart



! SVD ---------------  
   print*,ndata,ma,mp,np
        chisq=0d0
	u=0d0
        v=0d0
        w=0d0
      do i=1,ndata
        tmp=1./sig(i)
        do j=1,ma
          u(i,j)=fd(i,j)*tmp
        enddo
        b(i)=y(i)*tmp
     enddo

      call svdcmp(u,ndata,ma,mp,np,w,v)
      wmax=0.
      do j=1,ma
        if(w(j).gt.wmax)wmax=w(j)
      enddo
      thresh=TOL*wmax
print*,'limite de nettoyage',thresh,TOL,wmax
      ib=0
      do j=1,ma
        if(w(j).lt.thresh.and.ind.ne.0)then
           w(j)=0.  
           ib=ib+1
         endif
      enddo
        print*,'nbr de valeurs propres annulees',ib 

      call svbksb(u,w,v,ndata,ma,mp,np,b,a)

      chisq=0.
      do i=1,ndata
        sum=0.
        do j=1,ma
          sum=sum+a(j)*fd(i,j)
        enddo
        chisq=chisq+((y(i)-sum)/sig(i))**2
      enddo


if(ind.ne.0)print*,'TOL/CHI',tol,chisq

! ================================================
! calcul des residus ---------------



yy=0d0
           do j=1,ndata
              yy(j)=y(j)
              do i=1,ma
                 YY(j)=yy(j)-fd(j,i)*a(i)
              enddo
           enddo

xds=0d0
           do i=1,ndata
              xds=xds+yy(i)*yy(i)/sig(i)
           enddo

!print*,'xds',xds,ndata-ma

           Ktst=0
           do i=1,ma
              if(w(i).ne.0d0)Ktst(i)=1
           enddo
! ================================================
! calcul de la matrice de covariance ---------------
! 
! cov=tmp*v*w-1*w-1*vt avec vt=v-1
!
      bw=0
       do i=1,ma
        wti(i)=0.
        if(w(i).ne.0.) then
           wti(i)=1d0/(w(i)*w(i))
           bw=bw+1
        endif
       enddo

!print*,ma,bw

       do i=1,ma
        do j=1,i
          sum=0.
          do k=1,ma
            sum=sum+v(i,k)*v(j,k)*wti(k)
          enddo
!          cov2(i,j)=sum*tmp**2
!          cov2(j,i)=sum*tmp**2
          cov2(i,j)=sum
          cov2(j,i)=sum
        enddo
!write(*,'(6(f20.10,x))')(cov2(j,i),j=1,ma)
      enddo

        
! ================================================

xdx=0d0

      do i = 1, ma     
         xdx(i)=sqrt(cov2(i,i))
      do j = 1, ma
         cor(i,j) = abs(cov2(i,j)/sqrt(abs(cov2(i,i)*cov2(j,j))))
      end do
   end do


do i=1,ma
!write(*,*)a(i),xdx(i)
enddo
      return
END subroutine svdfit1501


function pythag (a,b)

   implicit none 
   real (kind=8) :: a, b, pythag
   real (kind=8) :: absa, absb

   absa = abs(a)
   absb = abs(b)
   if (absa > absb) then
      pythag = absa*sqrt(1.d0+(absb/absa)**2)
   else
      if (absb == 0.d0) then
         pythag = 0.d0
      else
         pythag = absb*sqrt(1.d0+(absa/absb)**2)
      end if
   end if
 
end function

