
c***********************************************************************
c
c     this subroutine calculates the solenoidal vel field.
c       q(n+1)=qhat-grad(dph)*dt ,  pr=dph
c    third order runge-kutta is used.
c
      subroutine updvp
      use param
      use local_arrays, only: q2,q3,dq,dph
      implicit none
      integer :: jc,jm,kc,km
      real    :: usurm, usukm


      do kc=1,n3m
        km=kmv(kc)
        usukm = al*dt*udx3c(kc)

!RO  Compute the q2 and q3 velocity components
!RO  by updating with the gradient of dph

        do jc=1,n2m
          jm=jmv(jc)
          usurm = al*dt*udx2c(jc)
          q2(jc,kc)=q2(jc,kc)-
     %      (dph(jc,kc)-dph(jm,kc))*usurm
          q3(jc,kc)=q3(jc,kc)-
     %      (dph(jc,kc)-dph(jc,km))*usukm

        enddo 

      enddo

      q2(n2,:)=q2(1,:)
      q3(:,1) =0.0d0
      q3(:,n3) =0.0d0
      return
      end

