c---------------------------------------------------------------------
      subroutine drivels
c
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
c      
      if(istep.eq.0)then
         if(ifld_tls.gt.0)ifredist(ifld_tls-1) = .true.
         if(ifld_cls.gt.0)ifclsredist(ifld_cls-1) = .true.
         dt_tls = dt
         dt_cls = dt

!     Turn off internal solvers
         idpss(ifld_tls-1) = -1
         idpss(ifld_cls-1) = -1
         return
      endif

      
      if(ifld_tls .gt. 0)then
         if(idpss(ifld_tls-1).eq.-1)then
            if(nid.eq.0)write(*,*)"Solving externally ifield =",ifld_tls
            call nekls_advance(ifld_tls,nsteps_tls) !internal solver must be off
         endif
      endif

      if(ifld_cls .gt. 0)then
         if(idpss(ifld_cls-1).eq.-1)then
            if(nid.eq.0)write(*,*)"Solving externally ifield =",ifld_cls
            call nekls_advance(ifld_cls,nsteps_cls) !internal solver must be off
         endif
      endif
      return
      end
c---------------------------------------------------------------------
      subroutine nekls_advance(ifld,nstepsls)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
c
      real simdt,simtime
      integer simstep
      integer idpss2(ldimt)

c     store physical time
      simstep = istep
      simdt = dt
      simtime = time

      istep = 0
      time = 0.
      if(ifld.eq.ifld_tls)dt = dt_tls
      if(ifld.eq.ifld_cls)dt = dt_cls

      call copy(idpss2,idpss,ldimt)
      do i=1,ldimt
         idpss(i)=-1  !Turn off everything
      enddo
      idpss(ifld-1)=0           !Helmholtz
      
      do i=1,nstepsls
         istep = istep+1
         call settime
         do igeom=1,ngeom
            call heat(igeom)
         enddo
      enddo
      call copy(idpss,idpss2,ldimt)
      
      
      if(ifld.eq.ifld_tls)dt_tls = dt
      if(ifld.eq.ifld_cls)dt_cls = dt
c     restore physical time
      dt = simdt
      istep = simstep
      time = simtime
      return
      end
