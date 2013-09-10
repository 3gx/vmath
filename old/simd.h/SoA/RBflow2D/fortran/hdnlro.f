
c***********************************************************************
      subroutine hdnlro
      use param
      use local_arrays, only: q2,q3,hro,dens
      implicit none
      integer :: jc,kc
      integer :: kpp,km,kp,jp,jmm,jpp
      real    :: h32,h33 
c
c     h term for the q3 momentum equation at i+1/2,j+1/2,k
c
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,kmv,kpv,jmv,jpv,q2,q3,udx2m,udx3m,dens,hro)
!$OMP$ PRIVATE(jc,kc,km,kpp,kp,jmm,jpp,jp,h32,h33)
      do kc=1,n3m
      km=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do jc=1,n2m
      jp=jc+1
      jmm=jmv(jc)
      jpp=jpv(jc)
c
c
c    rho q2 term
c
c
c                d  rho q_r 
c             -----------
c                d   r      
c
      h32=(q2(jp,kc)*(dens(jpp,kc)+dens(jc,kc))-
     %     q2(jc,kc)*(dens(jc,kc)+dens(jmm,kc))
     %    )*udx2m(jc)*0.5d0
c
c    rho q3 term
c
c
c                 d  rho q_x 
c                -----------
c                 d   x      
c
      h33=(q3(jc,kp)*(dens(jc,kpp)+dens(jc,kc))-
     %     q3(jc,kc)*(dens(jc,kc)+dens(jc,km))
     %    )*udx3m(kc)*0.5d0
c
      hro(jc,kc)=-(h32+h33)
c     
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
