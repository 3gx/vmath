
c***********************************************************************
      subroutine hdnl3
      use param
      use local_arrays, only: q2,q3,qcap,dens
      implicit none
      integer :: jc,kc
      integer :: km,kp,jp,jmm,jpp
      real    :: h32,h33 
      real    :: densit
c
c     h term for the q3 momentum equation at i+1/2,j+1/2,k
c
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,kmv,jmv,jpv,q2,q3,udx2c,udx3c,dens,qcap)
!$OMP$ PRIVATE(jc,kc,km,kp,jmm,jpp,jp,h32,h33,densit)
      do kc=2,n3m
      km=kmv(kc)
      kp=kc+1
      do jc=1,n2m
      jmm=jmv(jc)
      jpp=jpv(jc)
      jp=jc+1
c
c
c    q3 q2 term
c
c
c                d  q_x q_r 
c             -----------
c                d   r      
c
      h32=(((q2(jp,kc)+q2(jp,km))
     %     *(q3(jpp,kc)+q3(jc,kc)))
     %    -((q2(jc,kc)+q2(jc,km))
     %     *(q3(jc,kc)+q3(jmm,kc))))*udx2c(jc)*0.25d0
c
c    q3 q3 term
c
c
c                 d  q_x q_x 
c                -----------
c                 d   x      
c
      h33=((q3(jc,kp)+q3(jc,kc))*(q3(jc,kp)+q3(jc,kc))
     %    -(q3(jc,kc)+q3(jc,km))*(q3(jc,kc)+q3(jc,km))
     %    )*udx3c(kc)*0.25d0
c
c  add the buoyancy term
c
            densit=(dens(jc,kc)+dens(jc,km))*0.5d0
c
c
      qcap(jc,kc)=-(h32+h33)+densit
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
c
