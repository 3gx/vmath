
c***********************************************************************
c                                                                      *
c                       CONDIZIONI AL CONTORNO                         *
c                                                                      *
c***********************************************************************
      subroutine densbo
      use param
      implicit none
      integer :: j
      do j=1,n2
              denbn(j)=0.d0
              denbs(j)=1.d0
      enddo
      return
      end
C
