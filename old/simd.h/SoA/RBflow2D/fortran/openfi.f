
      subroutine openfi
      implicit none
      character*14 filth,filte,filba,filen,filso,filve,filnu,filbaper
      common/names/filth,filte,filba,filen,filso,filve,filnu,filbaper
      filbaper = 'filbaper.out'
      open(46,file='nfbou',status='old')
      read(46,'(a)')filth
      read(46,'(a)')filte
      read(46,'(a)')filba
      read(46,'(a)')filen
      read(46,'(a)')filso
      open(32,file=filth,status='unknown',access='sequential',
     % position='append')
      open(34,file=filba,status='unknown')
      open(35,file=filbaper,status='unknown')
      open(39,file=filen,status='unknown')
      !open(97,file=filte,status='unknown')
      open(97,file=filte,status='unknown',access='sequential',
     % position='append')
      filnu='nusse.out'
      filve='maxvel.out'
      !open(95,file=filnu,status='unknown')
       open(95,file=filnu,status='unknown',access='sequential',
     % position='append')
      open(96,file=filve,status='unknown')
      !rewind(32)
      rewind(34)
      rewind(35)
      rewind(39)
      !rewind(95)
      !rewind(97)
      close(46)
      return
      end   
