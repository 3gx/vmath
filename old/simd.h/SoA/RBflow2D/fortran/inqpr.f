
c***********************************************************************
c                                                                      *
c                       INITIAL CONDITION                              *
c                                                                      *
c***********************************************************************
      subroutine inqpr
      use param
      use local_arrays, only: q2,q3,dens
      implicit none
      integer :: j,k
      real :: xxx,yyy,eps

      eps=0.01d0
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2,n3,q2,q3,zm,rc,zz,rm,eps)
!$OMP$ PRIVATE(k,j,yyy,xxx)                     
      do k=1,n3
        do j=1,n2
         yyy=zm(k) 
         xxx=rc(j)            
         q2(j,k)=(2.0d0*yyy-6.0d0*yyy**2+4.0d0*yyy**3)*sin(3*xxx)*eps
         yyy=zz(k)          
         xxx=rm(j)
         q3(j,k)=-yyy**2*(1.0d0-yyy)**2*cos(3.1*xxx)*eps
c         q2(j,k)=0.0d0
c         q3(j,k)=0.0d0
       enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,dens,denbs,zm,denbn)
!$OMP$ PRIVATE(k,j)                     
       do k=1,n3m 
        do j=1,n2m    
             dens(j,k)= denbs(j) - (denbs(j) - denbn(j))
     %                   *zm(k)
          end do 
        end do
!$OMP  END PARALLEL DO
      return                                                            
      end                                                               
c
