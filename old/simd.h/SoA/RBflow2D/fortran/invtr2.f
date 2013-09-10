
************************************************************************
c                       SUBROUTINE INVTR2
c   This subroutine performs the computation of Q~~ for the q2 momentum 
c   equation (radial direction) by a factored implicit scheme.
c   For details see the introduction of INVTR1
c   
c
      subroutine invtr2
      use param
      use local_arrays, only: q2,pr,dq,rhs,dph,ru2
      implicit none
      integer :: jc,kc,km,kp,jp,jm
      real    :: udx2
      real    :: dcq2,dpx22
      real    :: d22q2,d33q2
      real    :: alre
#ifdef DEBUG
      real    :: cksum2
#endif      
      
cm
      alre=al/ren
c
c  compute the rhs of the factored equation
c  everything at i+1/2,j,k+1/2
c
c    points inside the flowfield
c
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,kmv,kpv,jmv,udx2c,ap2je,ac2je,am2je)
!$OMP$ SHARED(q2,ap3sk,ac3sk,am3sk,pr,dph,ru2,alre)
!$OMP$ SHARED(ga,ro,dt,al,rhs,jpv)
!$OMP$ PRIVATE(jc,kc,km,kp,jm,udx2,d22q2,d33q2,dcq2,dpx22,jp)
        do kc=1,n3m
          km=kmv(kc)
          kp=kpv(kc)
          do jc=1,n2m
           jm=jmv(jc)
           jp=jpv(jc)
           udx2 = al*udx2c(jc)

c
c   22 second derivative of q2
c
            d22q2= q2(jp,kc)*ap2je(jc)
     %            +q2(jc,kc)*ac2je(jc)
     %            +q2(jm,kc)*am2je(jc)
c
c   33 second derivative of q2
c
            d33q2=q2(jc,kp)*ap3sk(kc)
     %           +q2(jc,kc)*ac3sk(kc)
     %           +q2(jc,km)*am3sk(kc)

c
c    viscid terms
c
            dcq2=d22q2+d33q2
c
c
c   component of grad(pr) along 2 direction
c
            dpx22=(pr(jc,kc)-pr(jm,kc))*udx2

c
            rhs(jc,kc)=(ga*dph(jc,kc)+ro*ru2(jc,kc)
     %                    +alre*dcq2-dpx22)*dt

cm===========================================================
c
            ru2(jc,kc)=dph(jc,kc)
       enddo
      enddo
!$OMP  END PARALLEL DO

#ifdef DEBUG        
        cksum2=0.0d0
        do jc=1,n2m
        do kc=1,n3m
        cksum2=cksum2+rhs(jc,kc)
        enddo
        enddo
        write(*,*) 'rhscksum= ',cksum2
#endif

#ifdef DEBUG
      write(*,*) 'starting solq2j'
#endif
      call solq2j
#ifdef DEBUG
      write(*,*) 'starting solq2k'
#endif
      call solq2k
      return
      end
c
