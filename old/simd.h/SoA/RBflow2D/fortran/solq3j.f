
c
c
************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x2
c
      subroutine solq3j
      use param
      use tridiag
      use local_arrays, only : rhs
      implicit none
      real, dimension(m2m) :: amjl,apjl,acjl,fjl
      integer :: jc,kc
      real betadx

      betadx=beta*al
!$OMP  PARALLEL IF(n2m > 1000)
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(betadx,n2m,ap3j,am3j,ac3j,apjl,amjl,acjl)
!$OMP$ PRIVATE(jc)
!$OMP  DO
      do jc=1,n2m
        apjl(jc)=-ap3j(jc)*betadx
        amjl(jc)=-am3j(jc)*betadx
        acjl(jc)=1.-ac3j(jc)*betadx
      end do
!$OMP  END DO
!$OMP  END PARALLEL

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,rhs,acjl,apjl,amjl)
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
c
      return
      end
