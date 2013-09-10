
      subroutine radialgrid
      use param, only: pi,strr,n2, rc
      implicit none
      real, dimension(1:2*n2) :: rr, xr
      integer n2n,j,jj
      real a,b,eta,rgl,dru

      n2n=2*(n2-1)
      a =-1.d0
      b =1.d0
      dru= (b-a)/real(n2n)
      do j =1,n2n+1
        eta = a + real(j-1) * dru
        rgl= dsin(0.5d0*pi*eta)
        
        xr(j) = strr* rgl + (1.d0-strr)*eta
        
        write(22,*)j, xr(j) 
      enddo
      !write(*,*)xr
      a=0.
      b=1.
      do j=1,n2n+1
        rr(j) = 0.5*(b-a)*xr(j) + 0.5*(a+b)
        write(23,*)j,rr(j)
      enddo
      
      do j= 1,n2
        jj = n2+j-1
        rc(j) = rr(jj)-0.5
        write(24,*)j,rc(j)
      enddo

      return
      end                                                             
