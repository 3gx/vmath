
c***********************************************************************
c***********************************************************************
c                                                                      *
c   ************** WSON ************************                       *
c                                                                      *
c***********************************************************************
      subroutine wson(time)
      use param
      use local_arrays, only: q3,dens
      implicit none
      real :: time
      integer :: jc,kc,kp
      real :: anusin,my_anusin,vol,q3cen,fac2
    
      my_anusin=0.d0 
      anusin=0.d0    
      vol = 1.d0/(r0*alx3*dx2*dx3)

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,g2rm,g3rm,q3,dens)
!$OMP$ PRIVATE(kc,jc,kp,fac2,q3cen)
!$OMP$ REDUCTION(+:my_anusin)
      do kc=1,n3m
        kp = kc + 1
        do jc=1,n2m
          fac2 = g2rm(jc)*g3rm(kc)
            q3cen = (q3(jc,kc)+q3(jc,kp))*0.5d0
            my_anusin=my_anusin+dens(jc,kc)*q3cen*fac2
      enddo
      enddo
!$OMP  END PARALLEL DO
     
      anusin = my_anusin
           
#ifdef SINGLE
      anusin=1.d0 + sqrt(pra*ray)*anusin*vol
#else
      anusin=1.d0 + dsqrt(pra*ray)*anusin*vol
#endif
      
      open(95,file='nusse.out',status='unknown',access='sequential',
     %  position='append')
      write(95,*) time, anusin
      close(95)
      
      return
      end     
