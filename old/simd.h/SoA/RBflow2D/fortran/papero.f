c     this code is made for simulating three-dimensional flows in polar
c     ( cilindrical ) coordinates.
c     boundary condition are slip-walls ( =0 ) or non-slip wall ( =1)       
c     in the radial and the vertical direction.                  
c                                                                       
c     this code allows , by changing suitable indices , to introduce
c     the buoiancy term ( ibuo=1 )
c
c      navier-stokes equations are solved by a fractional step method   
c     ( Kim and Moin ) with the pressure in the first step.              
c            
c     The time advancement of the solution is obtained by a
c     Runge-Kutta 3rd order low storage scheme (Wray) or a 2nd   
c     order Adams-Bashfort scheme.                                      
c
c     The Poisson  equation for the pressure is solved directly                
c     introducing FFT in the azimutal and vertical direction.
c     Because of the presence of the walls in the vertical direction
c     the cosFFT in that direction is required
c                                                                       
c                 Roberto Verzicco and Paolo Orlandi                       
c                 dipartimento di meccanica ed aeronautica              
c                 universita' la sapienza di roma                 
c                                                                       
c                                                                       
c                                                                       
c     All variables are calculated in a staggered grid:
c   
c        Instead of velocities, "flux" variables are introduced to avoid
c        the singularity for r=0
c        q1=vtheta, q2= r*vr, q3= vz       
c
c        dq1,dq2,dq3 : velocity correction                              
c
c        qcap :divergence of the  non free divergent velocity field    
c
c        non linear terms:
c
c        ru1, ru2, ru3, ruro : old step      
c        h1,h2,h3, hro : new step           
c        dens : density
c        pr : pressure                                                  
c        dph : pressure correction                                      
c       pressure solver is in the package press.f                      
c       non linear terms are calculated in hdnl routines                
c       the invertions of momentum equations is performed in invtr      
c       routines                                             
c 
*======================================================================
*     MPI parallelization by Mohammad Emran, TU Ilmenau, Germnay.
*======================================================================                                                                     
      program papero
      use param
      use local_arrays
      implicit none
      character(len=4) :: dummy
      integer :: tfini,tfin,n,ns
      double precision :: ts,te
      INTEGER :: nt, i
c      double precision :: omp_get_wtime
      double precision :: rtc, sumt, sumt1
      real, dimension(1:m2, 1:m3, 4), target :: data1

      q2   => data1(1:m2, 1:m3, 1)
      q3   => data1(1:m2, 1:m3, 2)
      pr   => data1(1:m2, 1:m3, 3)
      dens => data1(1:m2, 1:m3, 4)

      nt =0

      !write(*,*) 'My thread ID:',TID
      nt = nt+1
      
c      if (myid .eq. 1)write(6,*) 'MPI tasks=', numtasks
c      if (myid .eq. 1)write(6,*) 'OMP thraeds per task=', nt
c      if (myid .eq. 1)write(6,*) 'No. of processors=', nt*numtasks
      
      !if (myid .eq. 1)write(*,*) 'Hello ierr', ierr
      !write(*,*) 'Hello myid numtasks', myid, numtasks
c*******************************************************
c******* Read input file bou.in by all processes********
c*******************************************************
c
      open(unit=15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,*) n2,n3,nsst,nwrit,nread
        read(15,301) dummy
        read(15,*) n2p,n3p
        read(15,301) dummy
        read(15,*) ntst,tprint,tpin,tmax,ireset
        read(15,301) dummy
        read(15,*) alx3,istr3,str3,rmed31,etdp3,strb3
        read(15,301) dummy
        read(15,*)strr,rext,istr,rmed,etdp,strb
        read(15,301) dummy
        read(15,*) ray,pra,dt,resid,cflmax
        read(15,301) dummy
        read(15,*) inslws,inslwn,inslwr,ibuo
        read(15,301) dummy
        read(15,*) istat,timeav,timerms,iresta
        read(15,301) dummy         
        read(15,*) nlev,dism,nn   
        read(15,301) dummy       
        read(15,*) idtv,dtmax,cfllim  
        read(15,301) dummy           
        read(15,*) nini,nfin,nstri
        read(15,301) dummy
        read(15,*) nson
301     format(a4)                
      close(15)
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
cm============================================
cm************ End of input file**************
cm============================================
      
      if( n2>m2 ) then
          write(6,*) 'Error: n2 must be <= m2'
          write(6,*) 'For optimum memory, set: n2 = m2'
          write(6,*) 'Check m2 in the beginning of param.f'
          write(6,*) 'Check n2 in the beginning of bou.in'
        STOP
c       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
      
      if( n3>m3 ) then
          write(6,*) 'Error: n3 must be <= m3'
          write(6,*) 'For optimum memory, set: n3 = m3'
          write(6,*) 'Check m3 in the beginning of param.f'
          write(6,*) 'Check n3 in the beginning of bou.in'
        STOP
c      call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
cm============================================
      r0 = rext
      rint = 0.d0
cm============================================     
cm      if(nson.ne.0) then
cm      open(unit=15,file='sonda.in',status='old')
cm        do i = 1,nson
cm          read(15,*) (coson(i,j),j=1,3)
cm        end do
cm      close(15)
cm      end if
cm===========================================     
c
c     DEFINITIONS FOR THE NATURAL CONVECTION
c
#ifdef SINGLE
      ren = sqrt(ray/pra)
      pec = sqrt(pra*ray)
#else
      ren = dsqrt(ray/pra)
      pec = dsqrt(pra*ray)
#endif
c                                                                       
c
cm==========================================
cm      call openfi
cm==========================================    
      call openfi
cm==========================================
c
      pi=2.d0*dasin(1.d0)                          
      tfini=dt*ntst                                                     
      tfin=tfini                                                        
      n2m=n2-1                                                          
      n3m=n3-1
cm======================                                                         
cm      n3mh=n3m/2+1 
cm======================                                                          
c
cm====================================================
cm====================================================                                                                             
      write(6,112)2.d0*r0/alx3
      write(32,112)2.d0*r0/alx3
  112 format(//,20x,'R A Y B E N ',//,10x,
     % '2D Cell with aspect-ratio:  D/H = ',f4.2)
      write(6,202) ray,pra
      write(32,202) ray,pra
  202 format(/,5x,'Parameters: ',' Ra=',e10.3,' Pr= ',e10.3)
      if(idtv.eq.1) then
         write(6,204) cflmax
         write(32,204) cflmax
  204 format(/,5x,'Variable dt and fixed cfl= ',
     % e11.4,/ )            
      else 
         write(6,205) dtmax,cfllim
         write(32,205) dtmax,cfllim
  205 format(/,5x,'Fixed dt= ',e11.4,' and maximum cfl=',
     %  e11.4,/ )            
cm====================================================    
      endif
cm====================================================      
c                                                                       
c     assign coefficients for time marching schemes                     
c
      if(nsst.gt.1) then   
        gam(1)=8.d0/15.d0
        gam(2)=5.d0/12.d0
        gam(3)=3.d0/4.d0
        rom(1)=0.d0
        rom(2)=-17.d0/60.d0
        rom(3)=-5.d0/12.d0
cm======================================================
        write(6,100) (gam(n),n=1,nsst),(rom(n),n=1,nsst)
  100   format(/,5x,'The time scheme is a III order Runge-Kutta'
     %  ,4x,'gam= ',3f8.3,4x,'ro= ',3f8.3)
cm======================================================
      else                                                              
        gam(1)=1.5d0
        gam(2)=0.d0
        gam(3)=0.d0
        rom(1)=-0.5d0
        rom(2)=0.d0
        rom(3)=0.d0
cm======================================================
        write(6,110) gam(1),rom(1)
  110   format(/,5x,'The time scheme is the Adams-Bashfort',4x,
     %   'gam= ',f8.3,4x,'ro= ',f8.3)
     
cm======================================================                                 
      endif                                                             
      do 10 ns=1,nsst
        alm(ns)=(gam(ns)+rom(ns))
   10 continue

cm======================================================
      
c     call mem_alloc

cm======================================================

c
c     the solution of the problem starts
c
cm======================================================

      ts=rtc()     

      call gcurv
      
      te=rtc()
      
        open(27,file="Total_time.out")
        write(27,*)"Total simulation time in sec.: ", te-ts
        close(27)
cm==================================================
c     call mem_dealloc
      
c     call MPI_FINALIZE(ierr)
cm==================================================


      sumt = 0
      do i=1,10
        print *, i, mytimer(i)
        sumt = sumt + mytimer(i)
      end do
      print  *, 'total: ', sumt

      do i=20,27
        print *, i, mytimer(i)
      end do

      print *, 'sum: ', sumt+sumt1

      print *, '40: ', mytimer(40) 
      print *, '41: ', mytimer(41) 
      print *, '50: ', mytimer(50) 
      print *, '51: ', mytimer(51) 
      
      print *, '61: ', mytimer(61) 
      print *, '62: ', mytimer(62) 

      print *, '81: ', mytimer(81) 
      print *, '82: ', mytimer(82) 
      

      stop                                                              
      end                                                               
c                                                                       
c 
c
c       

c                                                                                       
