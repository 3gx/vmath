
      subroutine ntrvb_line(am,ac,ap,f,nn,n)
      implicit none
      integer :: nn,n
      real, dimension(nn) :: am, ac, ap, f
      integer :: i,nm,ii
c
cm     dimension am(nn),ac(nn),ap(nn),f(mm,nn)
c
c  ******** reduction of trid. matrix to an upper rigth matrix
c
      do i=2,n
        ac(i)=ac(i)*ac(i-1)-ap(i-1)*am(i)
        ap(i)=ap(i)*ac(i-1)
        f(i)=f(i)*ac(i-1)-f(i-1)*am(i)
      enddo
c  ******** calculation of the unknown by backward elimination
c
      if(ac(n).eq.0.) then
        write(6,*) 'Error: Singular tridiagonal matrix 1 at n=',n
        stop
      endif

      f(n)=f(n)/ac(n)
      nm=n-1
      
      do ii=1,nm
        i=n-ii
      
        if(ac(i).eq.0.) then
          write(6,*) 'Error: Singular tridiagonal matrix 2 at i=',i
          stop
        endif
      
        f(i)=(f(i)-ap(i)*f(i+1))/ac(i)
      
      enddo
   
      return
      end     
