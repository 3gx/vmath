
***********************************************************************
c************************************************************************
c   this subroutine performs the inversion of the q2 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q2
c   are treated implicitly
c   in the first part the rhs is calculated
c         direction x2
c
      subroutine solq2j
      use param
      use local_arrays, only : rhs
      implicit none
      real, dimension(m2m):: amjl,apjl,acjl,fjl
      integer :: jc,kc
      real betadx

      betadx=beta*al

c     amjl(1)=0.
c     apjl(1)=0.
c     acjl(1)=1.
c     amjl(n2)=0.
c     apjl(n2)=0.
c     acjl(n2)=1.
c   thread spawning might be too costly here
!$OMP  PARALLEL IF(n2m > 1000)
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,ap2j,ac2j,am2j,betadx,apjl,acjl,amjl)
!$OMP$ PRIVATE(jc)
!$OMP  DO
      do jc=1,n2m
       apjl(jc)=-ap2j(jc)*betadx
       acjl(jc)=1.-ac2j(jc)*betadx
       amjl(jc)=-am2j(jc)*betadx
      end do
!$OMP  END DO
!$OMP  END PARALLEL

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,rhs,amjl,acjl,apjl)
!$OMP$ PRIVATE(kc,jc,fjl)
      do kc=1,n3m 

          do jc=1,n2m
            fjl(jc)=rhs(jc,kc)
          end do

          call tripvmy_line(amjl,acjl,apjl,fjl,1,n2m,m2)
 
          do jc=1,n2m
            rhs(jc,kc) = fjl(jc)  
          end do
      end do 
!$OMP  END PARALLEL DO
      return
      end
