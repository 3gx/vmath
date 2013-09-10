
c************************************************************************
c                                                                       *
c ****************************** subrout coetar  ********************** *
c                                                                       *
c    this subroutine calculates the coefficients for the                *
c    integration in the radial direction with non-uniform coor. trasf.  *
c                                                                       *
c************************************************************************
      subroutine coetar
      use param
      implicit none
      integer :: jm,jc,jp,km,kc,kp
      real:: a22,a22p,a22m
      real:: a33,a33m,a33p
      real:: ugmm
c
c  ***********  coefficients for q2   
c

c     am2j(1)=0.d0
c     ap2j(1)=0.d0
c     ac2j(1)=1.d0
c     am2j(n2)=0.d0
c     ap2j(n2)=0.d0
c     ac2j(n2)=1.d0
c     do jc=2,n2m
c     jm=jc-1
c     jp=jc+1
c     a22=dx2q/g2rc(jc)
c     a22p=1.d0/g2rm(jc)
c     a22m=1.d0/g2rm(jm)
c     ap2j(jc)=a22*a22p
c     am2j(jc)=a22*a22m
c     ac2j(jc)=-(a22*a22p+a22*a22m)
c     enddo

cm    PERIODIC
#ifdef DEBUG
      write(*,*) 'coetar: a*2j'
#endif
      do jc=2,n2m
      if(jc.eq.1) then
      jm=n2
      jp=jc+1
      elseif(jc.eq.n2) then
      jm=jc-1
      jp=1
      else
      jm=jc-1
      jp=jc+1
      endif
      a22=dx2q/g2rc(jc)
      a22p=1.d0/g2rm(jc)
      a22m=1.d0/g2rm(jm)
      ap2j(jc)=a22*a22p
      am2j(jc)=a22*a22m
      ac2j(jc)=-(a22*a22p+a22*a22m)
      enddo

c
c  ***********  coefficients for q2  (explicit part)
c
c     do jc=2,n2m
c        jm=jc-1
c        jp=jc+1
c        ugmm = dx2q/g2rc(jc)
c        a22p=1.d0/g2rm(jc)
c        a22m=1.d0/g2rm(jm)
c        a22= 1.d0/g2rm(jc) + 1.d0/g2rm(jm)
c        ap2je(jc)=a22p*ugmm
c        am2je(jc)=a22m*ugmm
c        ac2je(jc)=-a22*ugmm
c     end do

cm    PERIODIC
#ifdef DEBUG
      write(*,*) 'coetar: a*2je'
#endif
      do jc=1,n2m
      if(jc.eq.1) then
      jm=n2m
      else
      jm=jc-1
      endif
         ugmm = dx2q/g2rc(jc)
         a22p=1.d0/g2rm(jc)
         a22m=1.d0/g2rm(jm)
         a22= 1.d0/g2rm(jc) + 1.d0/g2rm(jm)
         ap2je(jc)=a22p*ugmm
         am2je(jc)=a22m*ugmm
         ac2je(jc)=-a22*ugmm
      end do
c
c
c  ***********  coefficients for q3   inner points
c  *********** are equal to those for dens
c
c     do jc=2,n2m
c     jp=jc+1
c     jm=jc-1
c     a22=dx2q/g2rm(jc)
c     a22p= a22/g2rc(jp)
c     a22m= a22/g2rc(jc)
c     ap3j(jc)=a22p
c     am3j(jc)=a22m
c     ac3j(jc)=-(a22p+a22m)
c     apscj(jc)=a22p
c     amscj(jc)=a22m
c     acscj(jc)=-(a22p+a22m)
c     enddo

cm    PERIODIC
#ifdef DEBUG
      write(*,*) 'coetar: a*scj'
#endif
      do jc=1,n2m
      if(jc.eq.n2m) then
      jp=1
      else
      jp=jc+1
      endif
      a22=dx2q/g2rm(jc)
      a22p= a22/g2rc(jp)
      a22m= a22/g2rc(jc)
      ap3j(jc)=a22p
      am3j(jc)=a22m
      ac3j(jc)=-(a22p+a22m)
      apscj(jc)=a22p
      amscj(jc)=a22m
      acscj(jc)=-(a22p+a22m)
      enddo
c     
c    external bound. conditions  q3
c     
c      jc=n2m
c      jp=jc+1
c      jm=jc-1
c       ugmm2=dx2q/g2rm(jc)/rm(jc)
c       am3j(jc)=rc(jc)/g2rc(jc)*ugmm2
c       ap3j(jc)=0.d0
c      if(inslwr.eq.1) then
cc                 no-slip wall
c       ac3j(jc)= - (ugmm2*rc(jc)/g2rc(jc) +
c     %              ugmm2*rc(jp)/g2rc(jp)*2.d0 )
c      else
cc                 dq3/dr=0 has been assumed at the external boundary
c       ac3j(jc)= - ugmm2*rc(jc)/g2rc(jc)
c      endif
c     
c    external bound. conditions  dens
c     
c     a22=dx2q/g2rm(jc)/rm(jc)
c     a22p=2.d0*a22*rc(jp)/g2rc(jp)
c     a22m= +a22*rc(jc)/g2rc(jc)
c     apscj(jc)=0.d0
c     amscj(jc)=a22m
c     acscj(jc)=-a22m
c     
c    internal boundary conditions   q3
c     
c     jc=1
c     jp=jc+1
c     a22=dx2q/g2rm(jc)/rm(jc)
c     a22p= +a22*rc(jp)/g2rc(jp)
c     ap3j(jc)=a22p
c     am3j(jc)=0.d0
c     ac3j(jc)=-a22p
c     apscj(jc)=a22p
c     amscj(jc)=0.d0
c     acscj(jc)=-a22p
c
c
c  ***********  coefficients for q3   for x3 differentation
c  c means centered that is at k location
c
#ifdef DEBUG
      write(*,*) 'coetar: a*3ck'
#endif
      am3ck(1)=0.d0
      ap3ck(1)=0.d0
      ac3ck(1)=1.d0
      am3ck(n3)=0.d0
      ap3ck(n3)=0.d0
      ac3ck(n3)=1.d0
      do kc=2,n3m
      km=kc-1
      kp=kc+1
      a33=dx3q/g3rc(kc)
      a33p=1.d0/g3rm(kc)
      a33m=1.d0/g3rm(km)
      ap3ck(kc)=a33*a33p
      am3ck(kc)=a33*a33m
      ac3ck(kc)=-(ap3ck(kc)+am3ck(kc))
      enddo
c
c
c  **coefficients for q1, q2 ap3sk,am3sk,ac3sk
c   ap3ssk,am3ssk,ac3ssk, psc   for x3 differentation
c  s means staggered that is at k+1/2 location
c
#ifdef DEBUG
      write(*,*) 'coetar: a*2sk'
#endif
      do kc=2,n3m-1
      kp=kc+1
      km=kc-1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      ap3sk(kc)=a33p
      am3sk(kc)=a33m
      ac3sk(kc)=-(ap3sk(kc)+am3sk(kc))
      ap3ssk(kc)=ap3sk(kc)
      am3ssk(kc)=am3sk(kc)
      ac3ssk(kc)=ac3sk(kc)
      enddo
c    
c    lower wall  bound. conditions  indicated by inslws
c    differemtiation of sec. derivative at 1+1/2
c    
      kc=1
      kp=kc+1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      ap3sk(kc)=a33p
      am3sk(kc)=0.d0
      ac3sk(kc)=-(a33p+inslws*a33m*2.d0)
      ap3ssk(kc)=ap3sk(kc)
      am3ssk(kc)=am3sk(kc)
      ac3ssk(kc)=-(a33p+a33m*2.d0)
c     ac3ssk(kc)=-(a33p)
c    
c    upper wall  bound. conditions  indicated by inslws
c    differentiation of sec. derivative at n3-1/2
c    
      kc=n3m
      kp=kc+1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      am3sk(kc)=a33m
      ap3sk(kc)=0.d0
      ac3sk(kc)=-(a33m+inslwn*a33p*2.d0)
      ap3ssk(kc)=ap3sk(kc)
      am3ssk(kc)=am3sk(kc)
      ac3ssk(kc)=-(a33m+a33p*2.d0)
c     ac3ssk(kc)=-(a33m)
C
      return
      end
c

