
c***********************************************************************
c
c                                                                       
c***********************************************************************
      subroutine divgck(qmax,qtot)
      use param
      use local_arrays, only: q2,q3
      implicit none
      real    :: qmax,qtot
      integer :: jc,kc,kp,jp
      real    :: dqcap,dvol,my_qmax
        
cm==================
      my_qmax =-100.d0
cm==================                                                                       
      dvol= 1.d0/(dx2*dx3)
      qtot=0.d0                                               
      qmax=0.d0                                                     
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(q2,q3,udx2m,udx3m,n2m,n3m)
!$OMP$ PRIVATE(jc,kc,kp,jp,dqcap)
!$OMP$ REDUCTION(max:my_qmax)
      do kc=1,n3m
        kp=kc+1
        do jc=1,n2m
          jp=jc+1
              dqcap=(q2(jp,kc)-q2(jc,kc))*udx2m(jc)
     %              +(q3(jc,kp)-q3(jc,kc))*udx3m(kc)
              my_qmax = max(abs(dqcap),my_qmax)          
c         write(*,*) jc,kc,dqcap
      enddo
      enddo
!$OMP  END PARALLEL DO 
     
      qmax = my_qmax

      qqmax=qmax
      qtot=qtot*dvol
      
      return     
      end         
