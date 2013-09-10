
c***********************************************************************
      subroutine initia
      use param
      use local_arrays
      implicit none
      integer :: j,k
      
cm=================================
      q2=0.d0
      q3=0.d0
      dens=1.d0
      pr=0.d0
      dph=0.d0
      rhs=0.d0
      dq=0.d0
      ru1=0.d0
      ru2=0.d0
      ru3=0.d0
      ruro=0.d0
      qcap=0.d0
      hro=0.d0
cm==================================         

      do j=1,n2
        amscj(j) = 0.d0
        acscj(j) = 0.d0
        apscj(j) = 0.d0
        rc(j) = 0.d0
        rm(j) = 0.d0
        g2rc(j) = 0.d0
        g2rm(j) = 0.d0
      end do
      do k=1,n3
        am3ssk(k) = 0.d0
        ac3ssk(k) = 0.d0
        ap3ssk(k) = 0.d0
        zz(k) = 0.d0
        zm(k) = 0.d0
        g3rc(k) = 0.d0
        g3rm(k) = 0.d0
      end do
      return 
      end   
