c#define _USEFORTRAN_
************************************************************************
      subroutine gcurv
      use param
      use local_arrays
c, only: q2,q3,dens                  
      implicit none
      integer :: ntstf, ntii,l
      real    :: time,cflm,dmax,dtot
      double precision :: ti(2), tin(3)
c      double precision    :: omp_get_wtime
      double precision :: rtc, t0, t1
      integer , parameter :: ndt = 10
c
c     Code for the computation of three-dimensional incompressible flows    
c     in cylindrical polar coordinates.                                             
c                                                                       
c     This code solves flow fields bounded in the x3 (axial) and 
c     x2 (radial) directions. All boundaries can be no-slip or free-slip
c     by setting the appropriate indices in the input file.
c     The geometry  (a cylindrical can) is given in the subroutine cordi.
c     The discretization is uniform in the axial (3) and azimuthal (1)
c     directions, while can be non-uniform in the radial direction.
c
c     The equations for the following variables
c                                                                       
c      q1=v(theta)    q2=v(r)*r     q3=v(zeta)                          
c                                                                       
c     are discretized by finite-difference schemes.                    
c     The introduction of the variable q2 is necessary to avoid the
c     problem of the singularity of the equations at the axis of symmetry
c     (r=0).
c
c     All spatial derivatives are discretized by central second-order
c     accurate finite-difference schemes including the non linear terms.                                   
c
c     The nonlinear terms are treated explicitly, while the viscous terms
c     are computed implicitly. This would lead to the inversion of a
c     large banded matrix, that however is avoided by introducing
c     a factored scheme bringing to the solution of three tridiagonal
c     matrices for each velocity component (subroutine INVTR*).
c                              
c     In time a fractional-step procedure is used in the version of 
c     Nagi Mansour introducing the pressure in the first step.                         
c
c     The non-linear terms and the cross derivatives of the viscous terms
c     are discretized by explicit  Adams-Bashfort or 3rd order Runge-Kutta
c     method (A. Wray, personal communication).                      
c
c     The scalar quantity Phi, which projects the provisional velocity field
c     onto  a divergence free field, is solved by a direct method. 
c     For the axial and azimuthal directions modified wave numbers coupled 
c     with trigonometric expansions (FFTs) are used. The equation is then
c     solved by simply inverting a tridiagonal matrix for the radial direction.
c     No explicit boundary conditions are necessary for this Poisson equation.      
c                                                                       
c     Other details of the scheme are given in the introduction of the  
c     subroutine TSCHEM
c
c     timings
c                                                                       
c
c
      tin(1) = rtc()
      
c
#ifdef DEBUG
      write(*,*) 'starting initia'
#endif
      call initia
C**************************
cm       imovie = 0
cm       tframe = 0.05
C**************************
c                                                                       
c     grid information                                                 
c                                                                       
#ifdef DEBUG
      write(*,*) 'starting meshes'
#endif
      call meshes
#ifdef DEBUG
      write(*,*) 'starting indic'
#endif
      call indic                                                        
#ifdef DEBUG
      write(*,*) 'starting cordin'
#endif
      call cordin
cm===================================                                                      
cm      if(imovie.eq.1) call pricmov
cm===================================
cm===================================
      write(6,754)n2,n3                                              
      write(32,754)n2,n3                                             
  754 format(/,5x,'numero di punti griglia: ',' n2= ',i5,
     % ' n3= ',i5/)                       
      write(6,755) 1.d0/dx2,1.d0/dx3,dt,ntst                  
      write(32,755) 1.d0/dx2,1.d0/dx3,dt,ntst               
  755 format(/,2x,' dx2=',e10.3,' dx3=',e10.3,' dt=',e10.3,
     % ' ntst=',i7,/)

cm===================================
cm===================================     
      
      time=0.d0
      ntii=0                                                            
c                                                                       
c   read or create initial fields                                       
c                                                                       
cm      call phini ! Fishpack version
#ifdef DEBUG
      write(*,*) 'starting phiniTP'
#endif
      call phiniTP                                                    
c
        do 22 l=2,ndv                                                   
           vmax(l)=0.d0
   22   continue                                                        
c
#ifdef DEBUG
      write(*,*) 'starting densbo'
#endif
      call densbo
c
c      create the initial conditions
c
      if(nread.eq.0) then
cm================================
        write(6,*)' nread=0 ---> cond. iniz. create nel code'     
        write(32,*)' nread=0 ---> cond. iniz. create nel code'
cm================================            
        ntii=0                                                          
        ntime=0                                                         
        time=0.d0
        cflm=0.d0
        
#ifdef DEBUG                                                               
        write(*,*) 'starting inqpr'
#endif
        call inqpr
        
cm======================================================================  
#ifdef DEBUG                                                               
        write(*,*) 'starting divgck' 
#endif
        call divgck(dmax,dtot)
        write(6,*)' initial divg dmax,dtot  ',dmax,dtot                 
        write(32,*)' initial divg dmax,dtot  ',dmax,dtot                
cm================================ 
cm
      else
cm
cm================================
       write(6,*)' nread=1 ---> cond. iniz. lette da file'      
       write(32,*)' nread=1 ---> cond. iniz. lette da file'
cm================================ 
cm      
#ifdef DEBUG                                                               
        write(*,*) 'starting inirea' 
#endif
       call inirea(ntii,time)
cm
cm=====================================================
       
#ifdef DEBUG                                                               
        write(*,*) 'starting divgck' 
#endif
       call divgck(dmax,dtot)
       write(6,*)' initial divg dmax,dtot  ',dmax,dtot                 
       write(32,*)' initial divg dmax,dtot  ',dmax,dtot       
cm==============================               
      endif                                                             
c
      ntstf=ntii+ntst                                                   
      ntii=ntii+1
cm================================
      write(6,711)tprint,ntii,ntstf,tpin
      write(32,711)tprint,ntii,ntstf,tpin
  711 format(3x,'check in cond : tprint =',f8.4,'  ntii =',i8,
     1       '  ntstf =',i8,2x,'tpin =',f8.4//)
cm================================ 
c       call densmc
            if(idtv.eq.1) then
            
cm================================
      write(6,*)ntime,time,vmax(2),vmax(3),
     % dt,dmax,densm,denmax,denmin
      write(32,*)ntime,time,vmax(2),vmax(3),
     % dt,dmax,densm,denmax,denmin
cm================================  
            else
cm================================
      write(6,*)ntime,time,vmax(2),vmax(3),
     % cflm,dmax, densm,denmax,denmin
      write(32,*)ntime,time,vmax(2),vmax(3),
     % cflm,dmax, densm,denmax,denmin
cm================================     
     
        cflm=cflm*dt
            endif
c
#ifdef DEBUG                                                               
        write(*,*) 'starting coetar' 
#endif
         call coetar
c      
c
c  ********* starts the time dependent calculation ***

      tin(2) = rtc()
      
       write(6,*) 'Initialization Time = ', tin(2) -tin(1), ' sec.'
c                                                                       
          ti(1) = rtc()

c**************************************
c*       opening the library
c**************************************
#ifndef _USEFORTRAN_
      mytimer(:) = 0;
      call lib_open(
     &  dx2, dx3,
     &  rc, rm, g2rc, g2rm,
     &  zz, zm, g3rc, g3rm,
c
     &  udx2rm, udx2c, udx2m,
     &  usg2rc, udx3c, udx3m,
c
     &  jmv, jmc, jpv, jpc, jmhv,
     &  kmc, kpc, kmv, kpv, kup, kum,
c
     &  ap2j,  ac2j,  am2j,
     &  ap2je, ac2je, am2je,
     &  ap3j,  ac3j,  am3j,
     &  apscj, acscj, amscj,
c
     &  ap3ck,  ac3ck,  am3ck,
     &  ap3sk,  ac3sk,  am3sk,
     &  ap3ssk, ac3ssk, am3ssk,
c
     &  ak2,
     &  amphj, apphj, acphjj,
     &  amphk, acphk, apphk,
c
     &  n2,  n3,
     &  ren, pec,
     &  alm,  gam, rom,
     &  rhs1, rhs2, rhs3
     &  )
#endif
c**************************************
c*       opening the library
c**************************************


      do 350 ntime=ntii,ntstf                                           
c
c     the calculation stops if the velocities are diverging for numerical
c     stability conditions (courant number restrictions)                
c
c          ti(1) = rtc()
c
#ifdef DEBUG                                                               
        write(*,*) 'starting cfl' 
#endif
            t0 = rtc()
            call cfl(cflm)
            t1 = rtc()
            mytimer(20) = mytimer(20) + t1 - t0
                  
            
            if(idtv.eq.1) then
              if(ntime.ne.1) then
                 dt=cflmax/cflm
                 if(dt.gt.dtmax) dt=dtmax
              endif
              if(dt.lt. 0.000001d0) go to 166
            else
              cflm=cflm*dt
              if(cflm.gt.cfllim) go to 165
            endif
            beta=dt/ren*0.5d0

            t0 = rtc();
#ifdef _USEFORTRAN_
            call tschem(time)
#else
            call lib_tschem(dt, beta, nsst, 
     &          denbs, denbn,
     &          q2, q3,
     &          pr, dens, hro, rhs,
     &          ru1, ru2, ru3, ruro,
     &          dq, qcap,
     &          dph)
#endif
            mytimer(21) = mytimer(21) + rtc() - t0


            if (mod(ntime,ndt) .eq. 0) then
              t0 = rtc();
              write(6,*) ' ---------------------------------------- '
              write(6,*) ' T = ',time,' NTIME = ',ntime,' DT = ',dt
              call outh(time, dt);
              ti(2) = rtc()
              write(6,*) 'Iteration Time= ', (ti(2)-ti(1))/ndt, ' sec.'
              ti(1) = ti(2)
              mytimer(27) = mytimer(27) + rtc() - t0
            end if

            time=time+dt
            
            ! STOP using goto statements in this simple conditionals.
            ! very bad and bug prone style
            if ((mod(time,tpin).lt.dt) .or. (ntime.eq.1)) then
              t0 = rtc();
              call vmaxv
              t1 = rtc()
              mytimer(22) = mytimer(22) + t1 - t0
              t0 = t1

              if(vmax(2).gt. 1000.d0) go to 266

              call cfl(cflm)
              t1 = rtc()
              mytimer(23) = mytimer(23) + t1 - t0
              t0 = t1

              call divgck(dmax,dtot)
              t1 = rtc()
              mytimer(24) = mytimer(24) + t1 - t0
              t0 = t1

              call densmc
              t1 = rtc()
              mytimer(25) = mytimer(25) + t1 - t0
              t0 = t1

              if(.not.(idtv.eq.1)) cflm=cflm*dt
              if(dmax.gt.resid) go to 169
            end if
            
            if(nson.ne.0) then
              t0 = rtc()
              if(mod(time,tpin).lt.dt) call wson(time)
              mytimer(26) = mytimer(26) + rtc() - t0
            end if

            !  write the flow field                                              
            if(mod(time,tprint/10.d0).lt.dt) then                 
                t0 = rtc();
                if(nwrit.eq.1) then 
                  call outpf(time)
                endif 
                mytimer(21) = mytimer(21) + rtc() - t0
            endif
          
            if(time.gt.tmax) go to 167
            if( (ti(2) -tin(1)) .gt. 85500.) then
              call continua(time)
              write(6,*) 'Restart file updated at exit'
              write(6,*) 'exit time=',ti(2) -tin(1),'sec'
               stop 'Exit after time limit'
            endif
c================================================           
c
  350 continue
  
        tin(3) = rtc()
        write(6,*) 'Total Iteration Time = ',tin(3) -tin(2),' sec.'
      call continua(time)
      go to 167                                                         
  165 continue 
cm======================================
      write(6,164)                                                      
      write(32,164)                                                     
  164 format(10x,'cfl too large  ')
cm======================================                                  
      go to 167                                                         
  166 continue
cm======================================
      write(6,168) dt 
      write(32,168) dt                                                  
  168 format(10x,'dt too small, DT= ',e14.7)
cm======================================   
      go to 167                                                         
  266 continue
cm======================================
      write(6,268)                                                      
      write(32,268)                                                     
  268 format(10x,'velocities diverged')
cm======================================                                 
      go to 167                                                         
  169 continue
  
cm======================================
      write(6,178) dmax                                 
      write(32,178) dmax                                 
  178 format(10x,'too large local residue for mass conservation : '     
     1       ,e12.5,' at ',1x,i3 )
     
cm======================================

  167 continue                                                          

#ifndef _USEFORTRAN_
      call lib_close()
#endif

c
      return                                                            
      end                                                               
c
