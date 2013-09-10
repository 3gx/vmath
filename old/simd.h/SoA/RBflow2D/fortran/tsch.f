c
************************************************************************
c
c           SUBROUTINE  TSCHEM
c
c   This subroutine manages the whole integration scheme.
c   The following equations are solved:          
c   
c    ~~     n
c   Q  -  Q                n         n       n-1   alp       2  ~~   n 
c  --------- = -alp*grad (P ) + gam*H + rho*H   + ----- nabla ( Q + Q )
c    d t                                          2 Re
c
c          i                           i               i
c   where H  are the nonlinear terms, P  the pressure Q  the velocities
c       ~~
c   and Q  the provisional non solenoidal velocity field.
c   The superscripts (~~, n, n-1) indicate the time step level. 
c                        n
c   The nonlinear terms H  are computed in the routines HDNL*, while
c   in the routines INVTR* are computed the remaining terms, updated
c   the non linear terms and inverted the equation to find the provisional
c   field at the new time step.
c       ~~
c   The Q  velocity field is projected onto a solenoidal field by a 
c   scalar Phi computed through the equation
c
c                         2            1          ~~
c                    nabla (Phi ) =  ------ div ( Q  )
c                                    alp dt
c
c   The right hand side of this equation is computed in the routine
c   DIVG, while the equation is solved in PHCALC.
c
c   In the routine UPDVP the solenoidal velocity field at the new time
c   step is then computed through
c
c                n+1  ~~
c               Q   = Q  - alt*dt grad (Phi)
c
c   Finally in the routine PRCALC is updated the pressure field
c
c                n+1   n        alp dt      2
c               P   = P + Phi - ------ nabla (Phi)
c                                2 Re
c
c   When the scalar field is computed (density, concentration,
c   temperature) the routines HDNLRO and INVTRRO are used. The same
c   strategy at the velocity field is used, except that the scalar
c   field does not need any correction.
c
c   All variables are located on a staggered grid with the velocities
c   on the faces of the computational cell and all the scalars at the
c   centre. This is important when terms belonging to different equations
c   are avaluated.
c
c   Further details of the scheme can be found in the paper
c   "A finite-difference scheme for three-dimensional incompressible
c    flows in cylindrical coordinates" by R. Verzicco and P. Orlandi
c    J. of Comp. Phys. 1996.
c
c
      subroutine tschem(time)
      use param
      use local_arrays
      implicit none
      real    :: time
      integer :: ns
      integer, save :: val = 0
c#define DEBUG
#ifdef DEBUG
      double precision :: cksum2,cksum3
      integer :: j,k
#endif
      integer, parameter :: inc = 1
      double precision :: t0, t1, rtc
    
       val = val + 1
      
c 
cm      REAL timef !etime_,t(2)
c
c   TIME INTEGRATION : implicit viscous, 3rd order RK (Adams Bashfort)  
c
      do ns=1,nsst                                                 
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)
        
#ifdef DEBUG
        print *, val, '0--------------------------0'
        cksum2=0.0d0
        cksum3=0.0d0
        do j=1+inc,n2m-inc
        do k=1+inc,n3m-inc
        cksum2=cksum2+q2(j,k)**2
        cksum3=cksum3+q3(j,k)**2
        enddo
        enddo
        write(*,*) 'q2cksum= ',cksum2
        write(*,*) 'q3cksum= ',cksum3
        write(*,*) 'starting hdnl2'
#endif

        t0 = rtc()

        call hdnl2                 !!  "    "      "
        t1 = rtc()
        mytimer(1) = mytimer(1) + t1 - t0
        t0 = t1

#ifdef DEBUG        
        cksum2=0.0d0
        do j=1+inc,n2m-inc
        do k=1+inc,n3m-inc
        cksum2=cksum2+dph(j,k)**2
        enddo
        enddo
        write(*,*) 'dph= ',cksum2
        write(*,*) 'starting hdnl3'
#endif
        
        call hdnl3                          !!  "    "      "
        t1 = rtc()
        mytimer(2) = mytimer(2) + t1 - t0
        t0 = t1
#ifdef DEBUG        
        cksum2=0.0d0
        do j=1+inc,n2m-inc
        do k=1+inc,n3m-inc
        cksum2=cksum2+qcap(j,k)**2
        enddo
        enddo
        write(*,*) 'qcap= ',cksum2
        write(*,*) 'starting hdnlro'
#endif

        
        call hdnlro                         !!  "    "      "
        t1 = rtc()
        mytimer(3) = mytimer(3) + t1 - t0
        t0 = t1
#ifdef DEBUG        
        cksum2=0.0d0
        do j=1+inc,n2m-inc
        do k=1+inc,n3m-inc
        cksum2=cksum2+hro(j,k)**2
        enddo
        enddo
        write(*,*) 'q2= ',cksum2
        write(*,*) 'starting invtr2'
#endif
        
         
        call invtr2                     !! DQ2HAT=Q2HAT-Q2(N)
        t1 = rtc()
        mytimer(4) = mytimer(4) + t1 - t0
        t0 = t1
#ifdef DEBUG
        cksum2=0.0d0
        do j=1+inc,n2m-inc
        do k=1+inc,n3m-inc
        cksum2=cksum2+q2(j,k)**2
        enddo
        enddo
        write(*,*) 'q2= ',cksum2
        write(*,*) 'starting invtr3'
#endif

        call invtr3     !! DQ3HAT=Q3HAT-Q3(N)
        t1 = rtc()
        mytimer(5) = mytimer(5) + t1 - t0
        t0 = t1
#ifdef DEBUG        
        cksum2 = 0.0
        do j=1+inc,n2m-inc
        do k=1+inc,n3m-inc
        cksum2=cksum2+q3(j,k)**2
        enddo
        enddo
        write(*,*) 'q3= ',cksum2
        write(*,*) 'starting divg'
#endif
        


        call divg !here only q3 ghost(up) cell are needed
        t1 = rtc()
        mytimer(6) = mytimer(6) +t1 - t0
        t0 = t1
#ifdef DEBUG        
        cksum2=0.0d0
        do j=1+inc,n2m-inc
        do k=1+inc,n3m-inc
        cksum2=cksum2+dph(j,k)**2
        enddo
        enddo
        write(*,*) 'dph= ',cksum2
        write(*,*) 'starting phcalcTP'
#endif

       
        call phcalcTP
        t1 = rtc()
        mytimer(7) = mytimer(7) + t1 - t0
        t0 = t1
#ifdef DEBUG        
        cksum2=0.0d0
        do j=1+inc,n2m-inc
        do k=1+inc,n3m-inc
        cksum2=cksum2+dph(j,k)**2
        enddo
        enddo
        write(*,*) 'dph= ',cksum2
        write(*,*) 'starting updvp'
#endif
       
        call updvp                 !! SOLENOIDAL VEL FIELD
        t1 = rtc()
        mytimer(8) = mytimer(8) + t1 - t0
        t0 = t1
#ifdef DEBUG        
        cksum2=0.0d0
        cksum3=0.0d0
        do j=1+inc,n2m-inc
        do k=1+inc,n3m-inc
        cksum2=cksum2+q2(j,k)**2
        cksum3=cksum3+q3(j,k)**2
        enddo
        enddo
        write(*,*) 'Q2*cksum= ',cksum2
        write(*,*) 'Q3*cksum= ',cksum3
        write(*,*) 'starting invtrro'
#endif
       
        call invtrro                              !! DENSITY  ==>>> 
        t1 = rtc()
        mytimer(9) = mytimer(9) + t1 - t0
        t0 = t1
#ifdef DEBUG        
        cksum2=0.0d0
        do j=1+inc,n2m-inc
        do k=1+inc,n3m-inc
        cksum2=cksum2+dens(j,k)**2
        enddo
        enddo
        write(*,*) 'dens= ',cksum2
        write(*,*) 'starting prcalc'
#endif


        call prcalc                         !! PRESSURE FIELD
        t1 = rtc()
        mytimer(10) = mytimer(10) + t1 - t0
        t0 = t1

#ifdef DEBUG        
        cksum2=0.0d0
        do j=1+inc,n2m-inc
        do k=1+inc,n3m-inc
        cksum2=cksum2+pr(j,k)**2
        enddo
        enddo
        write(*,*) 'pr= ',cksum2
#endif
      enddo

      return                                                            
      end                                                               
c
