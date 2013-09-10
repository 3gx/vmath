
c***********************************************************************
c******** This subroutine dumps snapshots ******************************
c***********************************************************************
      subroutine outpf(time)
      use param
      use local_arrays, only: q2,q3,dens
      implicit none
      real :: time
      integer :: itime, j,k, n2pp,n3pp
      character(32) :: namfi3
      character(5)  :: ipfi
c    
      itime=10.*time+0.5
      write(ipfi,99) itime
   99 format(i5.5)

      namfi3='boum'//ipfi//'.dat'
cm======================================================================
      n2pp=(n2-1)/n2p+1                                                 
      n3pp=(n3-1)/n3p+1                                                 
      open(99,file=namfi3,form='unformatted',status='unknown')
        rewind(99)
        write(99) n2pp,n3pp 
        write(99) re,time 
        write(99)   ((q2(j,k),j=1,n2,n2p),k=1,n3,n3p) 
        write(99)   ((q3(j,k),j=1,n2,n2p),k=1,n3,n3p) 
        write(99)   ((dens(j,k),j=1,n2,n2p),k=1,n3,n3p) 
      close(99)
cm======================================================================
      return                                                            
      end                                                               
