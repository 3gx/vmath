
c
c
c
************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x3
c
      subroutine solq3k
      use param
      use tridiag
      use local_arrays, only : q3,rhs
      implicit none
      real, dimension(m3) :: amkl,apkl,ackl, fkl
      real amkT(m3-1),ackT(m3),apkT(m3-1),appk(m3-2)
      integer :: jc,kc,info,ipkv(m3)
      real :: betadx,ackl_b
cm    dimension amkl(m3),apkl(m3),ackl(m3),fkl(m1,m3)
c  ********* compute the dq3* sweeping in the x3 direction
c
      betadx=beta*al
c
      do kc=1,n3
      if (kc.eq.1) then
        amkl(1)=0.d0
        apkl(1)=0.d0
        ackl(1)=1.d0
        
      elseif (kc.ge.2 .and. kc.le.n3m) then 
        ackl_b=1.-ac3ck(kc)*betadx
        amkl(kc)=-am3ck(kc)*betadx/ackl_b
        ackl(kc)=1.
        apkl(kc)=-ap3ck(kc)*betadx/ackl_b
      
      else 
        amkl(n3)=0.d0
        apkl(n3)=0.d0
        ackl(n3)=1.d0
      endif
      enddo

      amkT=amkl(2:n3)
      apkT=apkl(1:(n3-1))
      ackT=ackl

#ifdef SINGLE
      call sgttrf(n3,amkT,ackT,apkT,appk,ipkv,info)
#else
      call dgttrf(n3,amkT,ackT,apkT,appk,ipkv,info)
#endif


!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3,n3m,betadx,q3,rhs,amkT,apkT,ackT)
!$OMP$ SHARED(ipkv,appk,ac3ck)
!$OMP$ PRIVATE(jc,kc,fkl,ackl_b,info)
      do jc=1,n2m
            fkl(1)= 0.d0
          do kc=2,n3m
            ackl_b=1.-ac3ck(kc)*betadx
            fkl(kc) = rhs(jc,kc)/ackl_b
          enddo
            fkl(n3)= 0.d0
          
#ifdef SINGLE
          call sgttrs('N',n3,1,amkT,ackT,apkT,appk,ipkv,fkl,n3,info)
#else
          call dgttrs('N',n3,1,amkT,ackT,apkT,appk,ipkv,fkl,n3,info)
#endif
          do kc=1,n3
            q3(jc,kc)=q3(jc,kc) + fkl(kc)
          end do
      end do
!$OMP  END PARALLEL DO
      return
      end
