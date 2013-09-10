
c***********************************************************************

c***********************************************************************
      subroutine densmc
      use param
      use local_arrays, only: dens
      implicit none
      integer :: j,jc,kc
      real :: my_densm,my_denmax,my_denmin
      real ::  dden,vol,fac2,surf
      real ::  my_dden,my_vol, my_anussupp
      real :: del1,del,deln,del1n,udel1q,udelq,udel1qn,udelqn
      real :: grtlow,grtupp,fcder,fcdern
      
      dden = 0.0d0
      vol = 0.0d0
c
      densm=0.d0    
      denmax=-10.d0*10.d0
      denmin=10.d0*10.d0
cm============================
      my_dden=0.d0
      my_vol=0.d0
      my_densm=0.d0    
      my_denmax=-100.d0
      my_denmin=100.d0
cm============================      

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,dens,g2rm,g3rm)
!$OMP$ PRIVATE(jc,kc,fac2)
!$OMP$ REDUCTION(min:my_denmin)
!$OMP$ REDUCTION(max:my_denmax)
!$OMP$ REDUCTION(+:my_dden,my_vol) 
      do kc=1,n3m
        do jc=1,n2m
          fac2 = g2rm(jc)*g3rm(kc)
            my_denmax = max(my_denmax,dens(jc,kc))
            my_denmin = min(my_denmin,dens(jc,kc))
            my_dden   = my_dden+(dens(jc,kc)*fac2) 
            my_vol    = my_vol+(fac2)
cm============================================================          
      enddo
      enddo
!$OMP END PARALLEL DO 
cm
     
      denmax = my_denmax
      denmin = my_denmin
      
      dden = my_dden
      vol  = my_vol
     
      densm=(dden/vol)
  
c
c     COMPUTATION OF THE NUSSELT NUMBER AT THE 
c     LOWER AND UPPER WALLS
c

      anusslow = 0.d0
      anussupp = 0.d0
cm      
      my_anussupp = 0.d0
cm      
      del1 = zm(1)-zz(1)
      del  = zm(2)-zz(1)
      udel1q = 1.d0/del1**2
      udelq = 1.d0/del**2
      fcder = 1.d0/(1.d0/del1 - 1.d0/del)
      del1n = -zm(n3m)+zz(n3)
      deln  = -zm(n3m-1)+zz(n3)
      udel1qn = 1.d0/del1n**2
      udelqn = 1.d0/deln**2
      fcdern = 1.d0/(1.d0/del1n - 1.d0/deln)
      surf = r0*alx3
      
c     serial might be faster due to thread spawning time
!$OMP  PARALLEL IF(n2m > 1000)
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(g2rm,n2m,dx2,dens,udel1q,udelq,fcder,udelqn,udel1qn)
!$OMP$ SHARED(fcdern,denbn,denbs,n3m)
!$OMP$ PRIVATE(j,fac2,grtlow,grtupp)
!$OMP$ REDUCTION(+:anusslow,anussupp)
!$OMP  DO

      do j=1,n2m
        fac2 = g2rm(j)/dx2
        
         grtlow = (dens(j,1)*udel1q-dens(j,2)*udelq
     %              -denbs(j)*(udel1q-udelq) )
     %              *fcder
           grtupp = (dens(j,n3m)*udel1qn-dens(j,n3m-1)*udelqn
     %              -denbn(j)*(udel1qn-udelqn) )
     %              *fcdern
           anusslow = anusslow + grtlow*fac2
           anussupp = anussupp + grtupp*fac2
           
cm======================================================================   
      end do
cm=====================================================================
!$OMP  END DO
!$OMP  END PARALLEL
      anusslow = anusslow / surf
      anussupp = anussupp / surf
      
      return         
      end                                                               
