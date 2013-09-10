
c                                                                          
c***********************************************************************
      subroutine vmaxv
      use param
      use local_arrays, only: q2,q3
      implicit none
      real    :: my_vmax2,my_vmax3
      integer :: j,k
      
        my_vmax2=-100.d0
        my_vmax3=-100.d0
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n2m,n3m,q2,q3)
!$OMP$ PRIVATE(k,j)
!$OMP$ REDUCTION(max:my_vmax2,my_vmax3)
        do k=1,n3m
          do j=1,n2m
              my_vmax2 = max(my_vmax2,abs(q2(j,k)))
              my_vmax3 = max(my_vmax3,abs(q3(j,k)))
       enddo
       enddo
!$OMP  END PARALLEL DO 

      vmax(2) = my_vmax2
      vmax(3) = my_vmax3
      
      return   
      end     
