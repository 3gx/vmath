
c
c
************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x3
c
      subroutine solrok
      use param
      use tridiag
      use local_arrays, only : dens,rhs
      implicit none
      real, dimension(m3) :: amkl,apkl,ackl, fkl
      integer :: jc,kc,info,ipkv(m3m)
      real :: betadx,ackl_b
      real amkT(m3m-1),ackT(m3m),apkT(m3m-1),appk(m3-3)

cm      dimension amkl(m3),apkl(m3),ackl(m3),fkl(m1,m3) 
c  ********* compute the dq3* sweeping in the x3 direction
c
c
      betadx=0.5d0*al*dt/pec
     
      do kc=1,n3m
        ackl_b=1.-ac3ssk(kc)*betadx
        amkl(kc)=-am3ssk(kc)*betadx/ackl_b
        ackl(kc)=1
        apkl(kc)=-ap3ssk(kc)*betadx/ackl_b
      enddo

      amkT=amkl(2:n3m)
      apkT=apkl(1:(n3m-1))
      ackT=ackl(1:n2m)

#ifdef SINGLE
      call sgttrf(n3m,amkT,ackT,apkT,appk,ipkv,info)
#else
      call dgttrf(n3m,amkT,ackT,apkT,appk,ipkv,info)
#endif

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,betadx,ac3ssk,rhs,dens,amkT,ackT,apkT)
!$OMP$ SHARED(appk,ipkv)
!$OMP$ PRIVATE(jc,kc,ackl_b,amkl,ackl,apkl,fkl,info)
      do jc=1,n2m
          do kc=1,n3m
            ackl_b=1.-ac3ssk(kc)*betadx
            fkl(kc)=rhs(jc,kc)/ackl_b
          end do
          
#ifdef SINGLE
          call sgttrs('N',n3m,1,amkT,ackT,apkT,appk,ipkv,fkl,n3m,info)
#else
          call dgttrs('N',n3m,1,amkT,ackT,apkT,appk,ipkv,fkl,n3m,info)
#endif
          
          do kc=1,n3m
            dens(jc,kc)=dens(jc,kc) + fkl(kc)
          end do
      end do
!$OMP  END PARALLEL DO
c
c
      return
      end
