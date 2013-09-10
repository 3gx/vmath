
************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x2
c
      subroutine solroj
      use param
      use tridiag
      use local_arrays, only : rhs
      implicit none
      real, dimension(m2m):: amjl,apjl,acjl,fjl
      integer :: jc,kc
      real betadx

      betadx=0.5d0*al*dt/pec



!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,rhs,apscj,amscj,acscj,betadx)
!$OMP$ PRIVATE(kc,jc,fjl,apjl,amjl,acjl)
      do kc=1,n3m
          do jc=1,n2m
            apjl(jc)=-apscj(jc)*betadx
            amjl(jc)=-amscj(jc)*betadx
            acjl(jc)=1.d0-acscj(jc)*betadx       
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
