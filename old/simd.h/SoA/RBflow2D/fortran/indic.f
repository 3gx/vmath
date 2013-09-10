
************************************************************************
c
c     in this subroutine the indices ip,im,jp,jm,kp,km are calculated.
c     the indices jpc,jup,jmc,jum,kpc,kup,kmc,kum are indices 
c     for the walls
c
      subroutine indic
      use param
      implicit none
      integer :: jc,kc
c
c   streamwise direction
c
      do 4 kc=1,n3m
        kmv(kc)=kc-1
        kpv(kc)=kc+1
        if(kc.eq.1) kmv(kc)=kc
        if(kc.eq.n3m) kpv(kc)=kc
    4 continue
c
c     direction normal to non-slip walls
c
      do kc=1,n3m
        kpc(kc)=kpv(kc)-kc
        kmc(kc)=kc-kmv(kc)
        kup(kc)=1-kpc(kc)
        kum(kc)=1-kmc(kc)
      enddo

c
c   indices for radial direction
c   different for periodic BC
c
      do jc=1,n2m
        jmv(jc)=jc-1
        jpv(jc)=jc+1
        if(jc.eq.1) jmv(jc)=n2m
        if(jc.eq.n2m) jpv(jc)=1
      enddo
      do jc = 1,n2+1
       jmhv(jc) = mod(jc,n2m/2+1)
       if(jmhv(jc).eq.0) jmhv(jc) = n2m/2 + 1
      enddo
c
c   indices for radial direction
c   different for periodic BC
c
      do jc=1,n2m
        jpc(jc)=jpv(jc)-jc
        jmc(jc)=jc-jmv(jc)
        if (jc.eq.1) jmc(jc) = 1
        if (jc.eq.n2m) jpc(jc) = 1
      enddo
c

      return
      end
c
