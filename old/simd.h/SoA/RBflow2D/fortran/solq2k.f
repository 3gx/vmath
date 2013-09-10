
c
************************************************************************
c   this subroutine performs the inversion of the q2 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q2
c   are treated implicitly
c   in the first part the rhs is calculated
c        direction x3
c
      subroutine solq2k
      use param
      use tridiag
      use local_arrays, only : q2,rhs
      implicit none
      real, dimension(m3) :: amkl,apkl,ackl, fkl
      real amkT(m3-1),ackT(m3),apkT(m3-1),appk(m3-2)
      integer :: jc,kc,info,ipkv(m3)
      real :: betadx,ackl_b
c  ************ compute dq2 sweeping along the x3 direction      periodic
c
      betadx=beta*al
c
      do kc=1,n3m
        ackl_b=1.-ac3sk(kc)*betadx
        amkl(kc)=-am3sk(kc)*betadx/ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3sk(kc)*betadx/ackl_b
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
!$OMP$ SHARED(n2m,n3m,betadx,ac3sk,rhs,q2,amkT,ackT,apkT)
!$OMP$ SHARED(appk,ipkv)
!$OMP$ PRIVATE(jc,kc,ackl_b,amkl,ackl,apkl,fkl,info)
      do jc=1,n2m
          do kc=1,n3m
            ackl_b=1.-ac3sk(kc)*betadx
            fkl(kc)=rhs(jc,kc)/ackl_b
          end do
          
#ifdef SINGLE
          call sgttrs('N',n3m,1,amkT,ackT,apkT,appk,ipkv,fkl,n3m,info)
#else
          call dgttrs('N',n3m,1,amkT,ackT,apkT,appk,ipkv,fkl,n3m,info)
#endif
          
          do kc=1,n3m
            q2(jc,kc)=q2(jc,kc) + fkl(kc)
          end do
      end do
!$OMP  END PARALLEL DO
c
c
      return
      end
