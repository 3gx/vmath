
c**************************************************************
c
      subroutine cordin
      use param                                                 
      implicit none
      integer  :: j,k
      real, dimension(1:m3) :: etaz
      real :: tstr3, z2dp
      real :: x2,x3

c
c     RADIAL COORDINATE DEFINITION
c
      rint = 0.d0
      rext = r0
c
c     UNIFORM GRID
c
      if (istr.eq.0) then
        do  j=1,n2
         x2=real(j-1)/real(n2m)
         rc(j)= r0*x2
        end do
      end if
c
c     CLUSTERING AT THE EXTERNAL RADIAL WALL
c
c
c     STAGGERED COORDINATES AND
c     METRIC QUANTITIES
c
      do j=1,n2m
        rm(j)=(rc(j)+rc(j+1))*0.5d0
        g2rm(j)=(rc(j+1)-rc(j))*dx2
      end do
      do j=2,n2m
        g2rc(j)=(rc(j+1)-rc(j-1))*dx2*0.5d0
      end do
      g2rc(1)=(rc(2)-rc(1))*dx2
      g2rc(n2)=(rc(n2)-rc(n2m))*dx2
c
c     AXIAL COORDINATE DEFINITION
c
c
c     UNIFORM GRID
c
c       write(6,*) istr3,alx3

      tstr3=tanh(str3)

      if (istr3.eq.0) then
        do k=1,n3
          x3=real(k-1)/real(n3m)
          etaz(k)=alx3*x3
          zz(k)=etaz(k)
c          write(6,*)k,zz(k)
        enddo
      endif


c       
c      CLUSTERING AT THE EXTERNAL RADIAL WALL 
c                       and  
c             CLUSTERING AT THE AXIS 
c      

        if (istr3.eq.4) then
         zz(1)=0.0d0
         do k=2,n3
          z2dp=float(2*k-n3-1)/float(n3m)
          zz(k)=(1+tanh(str3*z2dp)/tstr3)*0.5*alx3
          if(zz(j).lt.0.or.zz(j).gt.alx3)then
           write(*,*)'Forza la griglia: ','zc(',k,')=',zz(k)
           stop
          endif
         end do

        end if


      
cm-----------------------------------------
c
c     STAGGERED COORDINATES AND
c     METRIC QUANTITIES
c
      do k=1,n3m
        zm(k)=(zz(k)+zz(k+1))*0.5d0
        g3rm(k)=(zz(k+1)-zz(k))*dx3
      enddo
      do k=2,n3m
        g3rc(k)=(zz(k+1)-zz(k-1))*dx3*0.5d0
      enddo
      g3rc(1)=(zz(2)-zz(1))*dx3
      g3rc(n3)= (zz(n3)-zz(n3m))*dx3
c
c     WRITE GRID INFORMATION
c
cm====================================================
      open(unit=98,file='radcor.out',status='unknown')
      do j=1,n2
        write(98,345) j,rc(j),rm(j),g2rc(j),g2rm(j)
      end do
      close(98)
      open(unit=78,file='axicor.out',status='unknown')
      do k=1,n3
        write(78,345) k,zz(k),zm(k),g3rc(k),g3rm(k)
      end do
      close(78)
 345  format(i4,4(2x,e23.15))
cm===================================================
c
c     QUANTITIES FOR DERIVATIVES
c
      open(unit=78,file='fact3.out',status='unknown')
      do k=1,n3m
        udx3m(k) = dx3/g3rm(k)
        udx3c(k) = dx3/g3rc(k)
        write(78,*) k,udx3m(k),udx3c(k)
      end do
      udx3c(n3) = dx3/g3rc(n3)
        write(78,*) n3,udx3m(n3m),udx3c(n3)
        close(78)
      open(unit=78,file='fact2m.out',status='unknown')
      do j=1,n2m
        usg2rc(j) = 1.d0/g2rc(j)
        udx2m(j) = dx2/g2rm(j)
        udx2c(j) = dx2/g2rc(j)
        write(78,*) j,udx2m(j),udx2c(j),usg2rc(j)
      end do
      usg2rc(n2) = 1.d0/g2rc(n2)
      udx2c(n2) = dx2/g2rc(n2)
        write(78,*) n2,udx2c(j)
        close(78)
      return                                                            
      end                                                               
