
************************************************************************
c                       SUBROUTINE INVTRRO
c   This subroutine performs the computation of he scalar field.
c   For details see the introduction of INVTR1
c
      subroutine invtrro
      use param
      use local_arrays, only: dens,dq,hro,ruro,rhs
      implicit none
      integer :: jc,kc,km,kp,jp,jpp,jmm
      real    :: dq32,dq33,dcq3 
      real    :: aap,aac,aam
      real    :: alpec
      real    :: del1,del2, fcder
c
      alpec=al/pec
c
c  compute the rhs of the factored equation
c  everything at i+1/2,j+1/2,k+1/2
c
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n3m,n2m,kmv,kpv,jmv,jpv,apscj,amscj,acscj)
!$OMP$ SHARED(dens,ap3ssk,am3ssk,ac3ssk,alpec,hro,ruro)
!$OMP$ SHARED(dt,ga,ro,rhs,zm,zz,denbn,denbs,n3)
!$OMP$ PRIVATE(jc,kc,km,kp,jmm,jpp,jp,aap,aam,aac,dq32)
!$OMP$ PRIVATE(dq33,dcq3,del1,del2,fcder)
      do kc=1,n3m
        if( (kc .ge. 2) .and. (kc .le. n3m-1) ) then
        km=kmv(kc)
        kp=kpv(kc)
        do jc=1,n2m
          jmm=jmv(jc)
          jpp=jpv(jc)
          jp=jc+1
          aap=apscj(jc)
          aam=amscj(jc)
          aac=acscj(jc)
c
c   22 second derivatives of dens
c
            dq32=dens(jpp,kc)*aap+dens(jc,kc)*aac
     %          +dens(jmm,kc)*aam
c
c   33 second derivatives of dens
c
               dq33=dens(jc,kp)*ap3ssk(kc)
     %              +dens(jc,kc)*ac3ssk(kc)
     %              +dens(jc,km)*am3ssk(kc)
c
            dcq3=dq32+dq33
c
c    right hand side of the density equation
c
cm===========================================================
            rhs(jc,kc)=(ga*hro(jc,kc)+ro*ruro(jc,kc)
     %              +alpec*dcq3)*dt
cm===========================================================
c
c    updating of the non-linear terms
c
            ruro(jc,kc)=hro(jc,kc)
        enddo
        endif
c
c
c       LOWER HOT WALL
c
cm
      
      del1 = zm(1)-zz(1)
      del2 = zm(2)-zm(1)
      fcder = 2.d0/(del1*del2*(del1+del2))
        if (kc .eq. 1) then
        kp = kc + 1
            do jc=1,n2m
              jmm=jmv(jc)
              jpp=jpv(jc)
              jp = jc + 1
          aap=apscj(jc)
          aam=amscj(jc)
          aac=acscj(jc)
c
c   22 second derivatives of dens
c
            dq32=dens(jpp,kc)*aap+dens(jc,kc)*aac
     %          +dens(jmm,kc)*aam
c
c   33 second derivatives of dens
c
            dq33=(dens(jc,kp)*del1-dens(jc,kc)*(del1+del2)
     %            +denbs(jc)*del2)
     %              *fcder
c
            dcq3=dq32+dq33
c
c    right hand side of the density equation
c
cm=======================================================
            rhs(jc,kc)=(ga*hro(jc,kc)+ro*ruro(jc,kc)
     %              +alpec*dcq3)*dt
cm=======================================================
c
c    updating of the non-linear terms
c
            ruro(jc,kc)=hro(jc,kc)
c       
      enddo
      endif
c
c       UPPER COLD WALL
c     
cm      csd=4./3.
      del1 = zz(n3)-zm(n3m)
      del2 = zm(n3m)-zm(n3m-1)
      fcder = 2.d0/(del1*del2*(del1+del2))
!m    kc = n3m
      if (kc .eq. n3m) then
        km = kc - 1
            do jc=1,n2m
              jmm=jmv(jc)
              jpp=jpv(jc)
              jp = jc + 1
          aap=apscj(jc)
          aam=amscj(jc)
          aac=acscj(jc)
c
c   22 second derivatives of dens
c
            dq32=dens(jpp,kc)*aap+dens(jc,kc)*aac
     %          +dens(jmm,kc)*aam
c
c   33 second derivatives of dens
c
            dq33=(dens(jc,km)*del1-dens(jc,kc)*(del1+del2)
     %            +denbn(jc)*del2)
     %              *fcder
c
            dcq3=dq32+dq33
c
c    right hand side of the density equation
c
c
c
cm========================================================
            rhs(jc,kc)=(ga*hro(jc,kc)+ro*ruro(jc,kc)
     %              +alpec*dcq3)*dt
cm========================================================
c
c    updating of the non-linear terms
c
            ruro(jc,kc)=hro(jc,kc)
        enddo
        endif
        enddo
!$OMP  END PARALLEL DO

      call solroj

      call solrok
c

      return
      end
