
c***********************************************************************
                                          
c***********************************************************************
      subroutine cfl(cflm)
      use param
      use local_arrays, only: q2,q3
      implicit none
      real    :: cflm, my_cflm
      integer :: j,k,jp,kp
      real :: qcf
      
      my_cflm=-100.d0
c                                                                       
      cflm=0.d0

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(q2,q3,n2m,n3m,udx2m,udx3m)
!$OMP$ PRIVATE(k,j,kp,jp,qcf)      
!$OMP$ REDUCTION(max:my_cflm)
      do k=1,n3m
        kp=k+1
        do j=1,n2m
          jp=j+1
            qcf= +abs((q2(j,k)+q2(jp,k))*0.5d0*udx2m(j))
     %           +abs((q3(j,k)+q3(j,kp))*0.5d0*udx3m(k))

cm===============================================================                
cm                cflm = max(cflm,qcf)
cm===============================================================
            my_cflm = max(my_cflm,qcf)
cm===============================================================
      enddo
      enddo
!$OMP END PARALLEL DO
            
      cflm = my_cflm

      return  
      end                                                               
