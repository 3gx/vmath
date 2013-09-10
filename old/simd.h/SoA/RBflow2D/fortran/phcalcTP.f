************************************************************************
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 and x1to use the real fourier transform
c
      subroutine phcalcTP
      use param
      use local_arrays, only: dph
      implicit none
      real, dimension(m2m) :: xr
      complex, dimension(m2mh) :: xa
      real, dimension(m3m-2)   :: appph
      real, dimension(m3m-1)   :: amphT, apphT
      real, dimension(m3m)     :: acphT,drhs
      real, dimension(m3m)     :: apph, amph
      integer, dimension(m3m)  :: phpiv

      integer :: jmh,n2mh
      integer :: j,k,info
      real :: acphT_b


c           
c
c   fft applied to the x2 direction to the
c   complex coeff. from cos fft
c   from physical to wave number space
      n2mh=n2m/2+1
!$OMP  PARALLEL DO 
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,dph,n2mh,fwd_plan)
!$OMP$ PRIVATE(k,j,xr,xa)
      do k=1,n3m
        do j=1,n2m
         xr(j)=dph(j,k)
        end do
#ifdef SINGLE
        call sfftw_execute_dft_r2c(fwd_plan,xr,xa)
#else
        call dfftw_execute_dft_r2c(fwd_plan,xr,xa)
#endif

        do j=1,n2mh
         dph(j,k)=real(xa(j)/n2m)
        end do
 
        do j=1,n2mh
         dph(j+n2mh,k)=aimag(xa(j)/n2m)
        end do 
      end do
!$OMP  END PARALLEL DO
        
cm==================================================================
c     inversion of the matrix in the x2 and x3 directions (real part)
c
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2,n3m,jmhv,ak2,acphk,apphk,amphk,dph)
!$OMP$ PRIVATE(acphT_b,jmh,apph,amph,drhs,acphT,amphT,apphT,info,phpiv)
!$OMP$ PRIVATE(appph)
      do j=1,n2+1
       jmh=jmhv(j)

       do k = 1,n3m
        acphT_b=acphk(k)-ak2(jmh)
        apph(k)=apphk(k)/acphT_b
        amph(k)=amphk(k)/acphT_b
        drhs(k)=dph(j,k)/acphT_b
        acphT(k)=1.0d0
       enddo

       amphT=amph(2:n3m)
       apphT=apph(1:(n3m-1))

#ifdef SINGLE
       call sgttrf(n3m, amphT, acphT, apphT, appph, phpiv, info)
       call sgttrs('N', n3m, 1, amphT, acphT, apphT, appph, phpiv, drhs,
     %               n3m, info)
#else
       call dgttrf(n3m, amphT, acphT, apphT, appph, phpiv, info)

       call dgttrs('N', n3m, 1, amphT, acphT, apphT, appph, phpiv, drhs,
     %               n3m, info)
#endif

        do k=1,n3m
          dph(j,k) = drhs(k)
        enddo

       enddo
!$OMP  END PARALLEL DO
c
c================================================================
c   inverse fft applied to the phi x1 direction
c   from wave number space to physical space
c    
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,n2mh,dph,bck_plan)
!$OMP$ PRIVATE(k,j,xa,xr)
      do k=1,n3m
       do j=1,n2mh
        xa(j)=cmplx(dph(j,k),dph(j+n2mh,k))
       end do


#ifdef SINGLE
      call sfftw_execute_dft_c2r(bck_plan,xa,xr)
#else
      call dfftw_execute_dft_c2r(bck_plan,xa,xr)
#endif


       do j=1,n2m
        dph(j,k)=xr(j)
       end do
      end do
!$OMP  END PARALLEL DO

      return
      end
c============================================================
c
      subroutine tridag(a,b,c,d,u,n,nn1)
      implicit none
      integer n,nn1
      real gam(n),a(n),b(n),c(n),d(n),u(n)
      integer j
      real bet

        if(b(1).eq.0.) then
         write(6,*) 'The solution does not exist'
        endif
        bet=b(1)
        u(1)=d(1)/bet
        do j=2,nn1
          gam(j)=c(j-1)/bet
          bet=b(j)-a(j)*gam(j)
          if(bet.eq.0.) then
            write(6,*) 'The solution does not exist'
          end if
          u(j)=(d(j)-a(j)*u(j-1))/bet
        end do
        do  j=NN1-1,1,-1
          u(j)=u(j)-gam(j+1)*u(j+1)
        end do
      return
      end
