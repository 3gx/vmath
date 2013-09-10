
************************************************************************
c   this subroutine performs the calculation of the pressure.
c   this depends on the fractional step
c
      subroutine prcalc
      use param
      use local_arrays, only: pr,dph
      implicit none
      integer :: kp,km,jm,jp,jc,kc
      real    :: be
c
c    the pressure is evaluated at the center of the box.
c
c     p^{n+1} = p^{n} + phi^{n+1} - b * Nabla^2 phi^{n+1}
c
      be=al*beta
      do kc=1,n3m
        kp=kpv(kc)
        km=kmv(kc)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
              pr(jc,kc)=pr(jc,kc)+dph(jc,kc)-be*(
     %        (dph(jc,kp)*apphk(kc)
     %        +dph(jc,kc)*acphk(kc)
     %        +dph(jc,km)*amphk(kc)) +
     %        (dph(jp,kc)*apphj(jc)
     %        +dph(jc,kc)*acphjj(jc)
     %        +dph(jm,kc)*amphj(jc))                    )
      enddo
      enddo
      
      return
      end
c
c
