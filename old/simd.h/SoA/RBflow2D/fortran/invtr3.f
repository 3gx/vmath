
c************************************************************************
c                       SUBROUTINE INVTR3
c   This subroutine performs the computation of Q~~ for the q3 momentum 
c   equation (axial direction) by a factored implicit scheme.
c   Viscous terms are treated implicitly, nonlinear terms explicitly.
c   For details see the introduction of INVTR1 
c
      subroutine invtr3
      use param
      use local_arrays, only: q3,dq,qcap, pr,ru3,rhs
      implicit none
      integer :: jc,kc,km,kp,jpp,jmm
      real    :: udx3
      real    :: dq32,dq33,dcq3,dpx33 
      real    :: aap,aac,aam
      real    :: alre
#ifdef DEBUG
      real :: cksum2
      integer :: j,k
#endif
cm      
      alre=al/ren
c
c  compute the rhs of the factored equation
c  everything at i+1/2,j+1/2,k
c
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,kmv,jmv,jpv,ap3j,am3j,ac3j,q3)
!$OMP$ SHARED(udx3c,ap3ck,ac3ck,am3ck,pr,rhs,ru3,alre)
!$OMP$ SHARED(ga,ro,al,qcap,dt)
!$OMP$ PRIVATE(jc,kc,km,kp,udx3,jmm,jpp,aap,aam,aac)
!$OMP$ PRIVATE(dq32,dq33,dcq3,dpx33)
      do kc=2,n3m
        km=kmv(kc)
        kp=kc+1
        udx3 = al*udx3c(kc)
        do jc=1,n2m
          jmm=jmv(jc)
          jpp=jpv(jc)
          aap=ap3j(jc)
          aam=am3j(jc)
          aac=ac3j(jc)
c
c   22 second derivatives of q3
c
            dq32=aam*q3(jmm,kc)+aac*q3(jc,kc)+aap*q3(jpp,kc)
c
c   33 second derivatives of q3
c
            dq33=q3(jc,kp)*ap3ck(kc)
     %          +q3(jc,kc)*ac3ck(kc)
     %          +q3(jc,km)*am3ck(kc)
c
c   viscous terms
c
            dcq3=dq32+dq33
c
c  component of grad(pr) along x3 direction
c
            dpx33=(pr(jc,kc)-pr(jc,km))*udx3
cm=======================================================     
            rhs(jc,kc)=(ga*qcap(jc,kc)+ro*ru3(jc,kc)
     %                    +alre*dcq3-dpx33)*dt 
cm=======================================================
c
c  updating of the non-linear terms
c
            ru3(jc,kc)=qcap(jc,kc)
       enddo
      enddo
!$OMP  END PARALLEL DO

c       if(kc.eq.1) then
c        do jc=1,n2m
c         rhs(jc,kc)=0.0d0
c        enddo
c       endif

#ifdef DEBUG        
        cksum2=0.0d0
        do j=1,n2m
        do k=1,n3m
        cksum2=cksum2+rhs(j,k)
        enddo
        enddo
        write(*,*) 'rhs cksum= ',cksum2
        write(*,*) 'starting solq3j'
#endif
      call solq3j
#ifdef DEBUG        
        cksum2=0.0d0
        do j=1,n2m
        do k=1,n3m
        cksum2=cksum2+rhs(j,k)
        enddo
        enddo
        write(*,*) 'rhs cksum= ',cksum2
        write(*,*) 'starting solq3k'
#endif
      call solq3k
      return
      end
c
c
