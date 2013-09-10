
#if 0
c======================================================================
      subroutine tripvmy_line(ami,aci,api,rrr,n1i,n1f,m1)
      implicit none
      integer :: n1i,n1f,m1
      real, dimension(m1) :: ami,aci,api,rrr
      real, dimension(m1) :: q,s,fei
      real    :: fn,p
      integer :: ia,ii,i,l
c                                                                       
c     vectorized for right hand side and coefficients                   
c                                                                                                                
cm      dimension ami(m1),aci(m1),api(m1),rrr(m1)
cm      dimension q(m1),s(m1),fei(m1)
      
      ia = n1i + 1
      ii = n1i + n1f 
c
c   COEFFICIENTS FOR TRIDIAGONAL INVERSION
c

c
c  THE INVERSION STARTS
c
      q(n1i) = -api(n1i)/aci(n1i) 
      s(n1i) = -ami(n1i)/aci(n1i)
      fn = rrr(n1f)
      rrr(n1i) = rrr(n1i)/aci(n1i)                   
c                                                                       
c     forward elimination sweep                                         
c                                                                       
      do 10 i=ia,n1f                                                
        p =1./( aci(i) + ami(i)*q(i-1))
        q(i) = - api(i)*p                    
        s(i) = - ami(i)*s(i-1)*p
        rrr(i) = ( rrr(i) - ami(i)*rrr(i-1))*p        
   10 continue             
c                                                                       
c     backward pass                                                     
c                                                                           
      s(n1f) = 1.   
      fei(n1f) = 0.
               
      do 11 l=ia,n1f       
        i = ii - l         
        s(i) = s(i) + q(i)*s(i+1)        
        fei(i) = rrr(i) + q(i)*fei(i+1)                     
   11 continue                     
              
      rrr(n1f)=(fn-api(i)*fei(n1i) - 
     %      ami(i)*fei(n1f-1))/(api(i)*s(n1i) +
     %      ami(i)*s(n1f-1)+aci(i))                  
c                                                                       
c     backward elimination pass                                         
c                                                                       
      do 12 l=ia,n1f         
        i = ii -l               
        rrr(i) = rrr(n1f)*s(i) + fei(i)                                
   12 continue
                                   
      return                                 
      end 
c======================================================================

#else

      subroutine tripvmy_line(ami,aci,api,rrr,n1i,n1f,m1)                   
      dimension ami(m1),aci(m1),api(m1),rrr(m1)
      dimension amiT(n1f-1), apiT(n1f-1), aciT(n1f)
      real appi(n1f-2), u(n1f), gamma
      integer ipiv(m1), info


c     Solves for a vector rrr(1:n1f) the “cyclic” set of linear equations

c     Choose a pivot as -b(1)
c      gamma=aci(1)/2
      gamma=-aci(1)
c      gamma=aci(1)

c     Set up the diagonal of the modified tridiagonal system.
c     alpha = api(n1f), beta = ami(1)

      aciT(1)=aci(1)-gamma
      aciT(2:(n1f-1)) = aci(2:(n1f-1))
      aciT(n1f)=aci(n1f)-api(n1f)*ami(1)/gamma

      amiT=ami(2:n1f)
      apiT=api(1:(n1f-1))

c      Factorize matrix

#ifdef SINGLE
       call sgttrf(n1f, amiT, aciT, apiT, appi, ipiv, info1)
#else
       call dgttrf(n1f, amiT, aciT, apiT, appi, ipiv, info1)
#endif

c      Solve A · x = r.

#ifdef SINGLE
       call sgttrs('N', n1f, 1, amiT, aciT, apiT, appi, ipiv, rrr,
     %    n1f, info1)
#else
       call dgttrs('N', n1f, 1, amiT, aciT, apiT, appi, ipiv, rrr,
     %    n1f, info1)
#endif


c      Set up the vector u.
       u(1)=gamma
       do i=2,n1f-1
        u(i) = 0.0
       enddo
       u(n1f)=api(n1f)

c      Solve A · z = u.

#ifdef SINGLE
       call sgttrs('n', n1f, 1, amiT, aciT, apiT, appi, ipiv, u,
     %    n1f, info2)
#else
       call dgttrs('n', n1f, 1, amiT, aciT, apiT, appi, ipiv, u,
     %    n1f, info2)
#endif

c      Form v · x/(1 + v · z).

       fact=(rrr(1)+ami(1)*rrr(n1f)/gamma)/
     %      (1.+u(1)+u(n1f)*ami(1)/gamma)

c      Now get the solution vector x.
       
       do i=1,n1f
        rrr(i)=rrr(i)-fact*u(i)
       enddo

      return
      end
#endif
