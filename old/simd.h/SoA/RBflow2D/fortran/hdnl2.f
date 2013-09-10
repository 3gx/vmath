
c***********************************************************************
      subroutine hdnl2
      use param
      use local_arrays, only: q2,q3,dph
      implicit none
      integer :: kc,kp,jp,jm,jc
      integer :: kmm,kpp
      real    :: h22,h23 

c
c     h term for the q2 momentum equation at i+1/2,j,k+1/2
c
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,kmv,kpv,jmv,jpv,q2,q3,udx2c,udx3m,dph)
!$OMP$ PRIVATE(jc,kc,kmm,kpp,kp,jm,jp,h22,h23)
      do kc=1,n3m
      kmm=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      
c     q2 q2 term
c
c
c                 d  q_r q_r 
c                ------------
c                 d   r      
c
      h22=( (q2(jp,kc)+q2(jc,kc))
     %     *(q2(jp,kc)+q2(jc,kc))
     %     -(q2(jc,kc)+q2(jm,kc))
     %     *(q2(jc,kc)+q2(jm,kc))
     %    )*udx2c(jc)*0.25d0
c
c     q2 q3 term
c
c
c                 d  q_x q_r 
c                -----------
c                 d   x      
c
      h23=((q3(jc,kp)+q3(jm,kp))*(q2(jc,kpp)+q2(jc,kc))
     %    -(q3(jc,kc)+q3(jm,kc))*(q2(jc,kc)+q2(jc,kmm))
     %    )*udx3m(kc)*0.25d0
c
      dph(jc,kc)=-(h22+h23)
c
      enddo
      enddo
!$OMP  END PARALLEL DO
      
      return
      end
c
