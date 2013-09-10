
c***********************************************************************
      subroutine inirea(ntii,time)
      use param
      use local_arrays, only: q2, q3, dens
      implicit none
      integer :: ntii,j,k
      real    :: time
      character(12) :: filcnw
      
      write(6,201)time
  201 format(10x,'At t=',e10.3,' reading from the restart file')

cm==============================================================      
      filcnw = 'continua.dat'
      open(13,file=filcnw,form='unformatted',status='unknown')
        rewind(13)                                                      
        read(13)n2,n3
        read(13) re,time 
        read(13) ((dens(j,k),j=1,n2),k=1,n3)
        read(13) ((q2(j,k),j=1,n2),k=1,n3)
        read(13) ((q3(j,k),j=1,n2),k=1,n3)
      close(13)
cm=============================================================

      if (ireset.eq.1) then                                             
cm      ihist=0                                                           
      ntii=0                                                            
      time=0.d0                                        
      endif                                                             
      return                                                            
      end                                                               
c                                                                       
