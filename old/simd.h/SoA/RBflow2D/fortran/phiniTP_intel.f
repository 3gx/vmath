************************************************************************
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine fftqua
      use param
      integer :: n2mh,n2mp,j
cm      pi=2.*asin(1.)
      n2mh=n2m/2+1
      n2mp=n2mh+1
c
c     wave number definition
c
      do 18 j=1,n2mh
        ap(j)=(j-1)*2.d0*pi
   18 continue
      do 19 j=n2mp,n2m
        ap(j)=-(n2m-j+1)*2.d0*pi
   19 continue
c      call fftfax(nx2fft,ifx1,trigx1)
c
c   modified wave number

c    
      do 28 j=1,n2m
#ifdef SINGLE
        ak2(j)=2.d0*(1.d0-cos(ap(j)/n2m))*(float(n2m)/rext)**2
#else
        ak2(j)=2.d0*(1.d0-dcos(ap(j)/n2m))*(float(n2m)/rext)**2
#endif
   28 continue

#if 1
! egaburov: modification to make sure code works with vectorized version
      do j=1,n2m
       ak2(j)= ak2(jmhv(j))
      end do
#endif
c
      return
      end
c=====================================================
      subroutine phiniTP
      use param
      implicit none
      real, dimension(n2m) :: xr
      complex, dimension(m2mh) :: xa
      integer FFTW_EXHAUSTIVE
      parameter(FFTW_EXHAUSTIVE=64)

#ifdef SINGLE
      call sfftw_plan_dft_r2c_1d(fwd_plan,n2m,xr,xa,FFTW_EXHAUSTIVE)
      call sfftw_plan_dft_c2r_1d(bck_plan,n2m,xa,xr,FFTW_EXHAUSTIVE)
#else
      call dfftw_plan_dft_r2c_1d(fwd_plan,n2m,xr,xa,FFTW_EXHAUSTIVE)
      call dfftw_plan_dft_c2r_1d(bck_plan,n2m,xa,xr,FFTW_EXHAUSTIVE)
#endif
      
      call tridiag_matrices

      return
      end
      
      
c=======================================================================
      subroutine tridiag_matrices
      use param
      implicit none
      integer  :: jc,jp,kc,km,kp
      real :: ugmmm,a22icc,a22icp,ac2,a33icc,a33icp
c
      call fftqua
c
c   tridiagonal matrix coefficients at each k and i
c   x1 and x3 cartesian coordinates
c
      do jc=1,n2m
        jp=jpv(jc)
        a22icc=jmc(jc)*dx2q/g2rc(jc)
        a22icp=jpc(jc)*dx2q/g2rc(jp)
        ac2=-(a22icc+a22icp)
        ugmmm=1.0d0/(g2rm(jc))
        amphj(jc)=a22icc*ugmmm
        apphj(jc)=a22icp*ugmmm
        acphjj(jc)=-(amphj(jc)+apphj(jc))
      enddo
      do kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        a33icc=kmc(kc)*dx3q/g3rc(kc)
        a33icp=kpc(kc)*dx3q/g3rc(kp)
        ugmmm=1.0d0/g3rm(kc)
        amphk(kc)=a33icc*ugmmm
        apphk(kc)=a33icp*ugmmm
        acphk(kc)=-(amphk(kc)+apphk(kc))
      enddo
c
      end subroutine tridiag_matrices
c==================================================
