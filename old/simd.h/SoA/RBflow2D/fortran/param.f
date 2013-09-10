#include "../dims.h"

c===========================================================
c Declaration of global variables
c***********************************************************      
      module param
        implicit none
        integer :: m2,m3
        integer :: m2m,m3m,m2mh
c       Grid size        
        parameter (m2=__M2__,m3=__M3__)
        parameter (m2m=m2-1,m3m=m3-1,m2mh=m2m/2+1)
c==========================================================			
c       read from input file bou.in
c==========================================================
        integer   :: n2, n3
        integer   :: n2p, n3p
        integer   :: nsst, nwrit, nread, ntst, ireset
        real      :: tprint,tpin,tmax
        real      :: alx3,str3,rmed31,etdp3,strb3
        integer   :: istr,istr3
        real      :: strr,rext,rmed,etdp,strb
        real      :: ray,pra,dt,resid,cflmax
        integer   :: inslws,inslwn,inslwr,ibuo
        integer   :: istat,iresta
        real      :: timeav,timerms
        integer   :: nlev, nn
        real      :: dism
        real      :: dtmax,cfllim 
        integer   :: idtv
        integer   :: nini,nfin,nstri
        integer   :: nson
        double precision :: mytimer(100)
c=================================================
c       end of input file
c=================================================
c******* Grid parameters**************************
        real :: dx2,dx3
        real :: dx2q,dx3q
c        
        real, dimension(1:m2) :: rc,rm,g2rc,g2rm
        real, dimension(1:m3) :: zz,zm,g3rc,g3rm
c====================================================
c******* QUANTITIES FOR DERIVATIVES******************
        real, dimension(1:m2) :: udx2rm,udx2c
        real, dimension(1:m2) :: udx2m
        real, dimension(1:m2) :: usg2rc
        real, dimension(1:m3) :: udx3c,udx3m
c==========================================================
c******* Grid indices**************************************
        integer, dimension(1:m2) :: jmv,jmc,jpv,jpc
        integer, dimension(1:m2+1) :: jmhv
        integer, dimension(1:m3) :: kmc,kpc,kmv,kpv,kup,kum
c===========================================================
c******* Metric coefficients *******************************
        real, dimension(1:m2) :: ap2j,ac2j,am2j
        real, dimension(1:m2) :: ap2je,ac2je,am2je
        real, dimension(1:m2) :: ap3j,ac3j,am3j
        real, dimension(1:m2) :: apscj,acscj,amscj

        real, dimension(1:m3) :: ap3ck,ac3ck,am3ck
        real, dimension(1:m3) :: ap3sk,ac3sk,am3sk
        real, dimension(1:m3) :: ap3ssk,ac3ssk,am3ssk   
c============================================================
c******* Variables for FFT and Poisson solver****************
c        real w(m2*m3) ! This is for fishpack
        real, dimension(13) :: ifx1
        real, dimension(3*m2/2+1) :: trigx1
        real, dimension(1:m2) :: ak2,ap
        real, dimension(1:m2) :: amphj,apphj,acphjj
        real, dimension(1:m3) :: amphk,acphk,apphk
        integer*8 :: fwd_plan,bck_plan
        
c===========================================================
c******* Other variables ***********************************
        integer  :: n2m, n3m
        integer  :: iaxsy
        real :: rint, r0
        real :: ren, pec
        real :: pi
        real :: al,ga,ro
        real :: beta
        real :: qqmax,qqtot
        real :: re
        real :: anusslow,anussupp
        real :: denmax,denmin,densm
        integer :: ntime
        integer, parameter:: ndv=3
        real, dimension(2:ndv) :: vmax
        real, dimension(1:3) :: gam,rom,alm
        real, dimension(1:m2) :: denbs,denbn
              
      end module param
      
c************* End of param module******************************
c===============================================================
c******* 3D arrays, dynamically allocated by each process*******
      module local_arrays
      use param
        implicit none
#if 0      
        real, dimension(1:m2,1:m3) :: q2,q3, pr,dens
#else
        real, pointer :: q2(:,:), q3(:,:), pr(:,:), dens(:,:) !allocated in  papero.f
#endif

        real, dimension(1:m2,1:m3) :: hro,rhs
        real, dimension(1:m2,1:m3) :: rhs1, rhs2, rhs3
        real, dimension(1:m2,1:m3) :: ru1,ru2,ru3,ruro
        real, dimension(1:m2,1:m3) :: dq,qcap
  
        real, dimension(1:m2+16,1:m3) :: dph
        
      end module local_arrays
c===============================================
      module tridiag
      implicit none
        
        interface 
        
        subroutine ntrvb(am,ac,ap,f,mm,nn,n,m)
          implicit none
          integer :: mm,nn,n,m
          real, pointer, dimension(:):: am,ac,ap
          real, pointer, dimension(:,:):: f
        
        end subroutine ntrvb
        
        subroutine ntrvb_test(am,ac,ap,f,mm,nn,n,m)
          implicit none
          integer :: mm,nn,n,m
          real, pointer, dimension(:):: am,ac,ap
          real, pointer, dimension(:,:):: f
        
        end subroutine ntrvb_test
        
        end interface
      end module tridiag
c====================================================
      module lapack_blas
        use param, only: m3
        implicit none
        double precision, dimension(m3) :: wr
        double precision, dimension(m3,m3) :: zmm,zmmt
 
      end module lapack_blas
