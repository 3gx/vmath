
c***********************************************************************
c***********************************************************************
      subroutine outh(time,cflm)
      use param
      implicit none
      real ::  time,cflm
      character*14 filth,filte,filba,filen,filso,filve,filnu,filbaper
      common/names/filth,filte,filba,filen,filso,filve,filnu,filbaper
c                                                                       

      write(6,159)ntime,time,
     % vmax(2),
     % vmax(3),
     % cflm,qqmax,densm,denmax,denmin
       open(32,file=filth,status='unknown',access='sequential',
     % position='append')
      write(32,159)ntime,time,
     % vmax(2),
     % vmax(3),
     % cflm,qqmax,densm,denmax,denmin
       close(32)
  159 format(1x,i5,2x,e10.4,4x,/,4x,
     % 3(e9.3,1x),/,4x,e9.3,
     % 2x,e9.3,1x,3(e12.6,1x))
c
       open(97,file=filte,status='unknown',access='sequential',
     % position='append')
       write(97,546) time, densm, anusslow, anussupp
 546   format(4(1x,e14.6))
       close(97)
c
       open(96,file=filve,status='unknown')
       write(96,*) time, vmax(2), vmax(3)
       close(96)
      return                                                            
      end                                                               
