
************************************************************************
c
c     this subroutine calculates divg(q).
c     q are the "fluxes"
c      q2=v(r)     q3=v(zeta)
c
      subroutine divg
      use param
      use local_arrays, only: q2,q3,dph
      implicit none
      integer :: jc,jp,kc,kp
      real    :: usdtal,dqcap   

      usdtal = 1.d0/(dt*al)
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,q2,q3,udx2m,udx3m,usdtal,jpv,dph)
!$OMP$ PRIVATE(jc,kc,kp,jp,dqcap)
      do kc=1,n3m
        kp=kc+1
        do jc=1,n2m
          jp=jpv(jc)
              dqcap=(q2(jp,kc)-q2(jc,kc))*udx2m(jc)
     %              +(q3(jc,kp)-q3(jc,kc))*udx3m(kc)
              dph(jc,kc)=dqcap*usdtal
      enddo
      enddo
!$OMP  END PARALLEL DO
      return
      end
c
c
c
