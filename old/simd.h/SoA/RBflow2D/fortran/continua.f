
c***********************************************************************
      subroutine continua(time)
      use param
      use local_arrays, only: q2,q3,dens
      implicit none
      real :: time
      integer :: j,k
      character(12) :: filcnw
      
      write(6,201)time
  201 format(10x,'At t=',e10.3,' restart file is updated')
cm==============================================================
      filcnw = 'continua.dat'      
      open(13,file=filcnw,form='unformatted',status='unknown')
       rewind(13)                                                      
        write(13)n2,n3
        write(13) re,time 
        write(13) ((dens(j,k),j=1,n2),k=1,n3)
        write(13) ((q2(j,k),j=1,n2),k=1,n3)
        write(13) ((q3(j,k),j=1,n2),k=1,n3)
      close(13)
cm==============================================================
      return                                                            
      end                                                               
