C----------------------------------------------------------------------     
      include "experimental/lshmholtz.f"
C----------------------------------------------------------------------     
      subroutine ls_init(nsteps_in,eps_in,dt_in,
     $                   ifld_cls_in,ifld_clsr_in,
     $                   ifld_tls_in,ifld_tlsr_in,
     $                    ifdebug)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      real eps_in,dt_in
      integer nsteps_in
      integer ifld_cls_in, ifld_tls_in
      integer ifld_clsr_in, ifld_tlsr_in
      integer ntot, ifdebug
      
      ! multiple of element length
      eps_cls = eps_in
      
      nsteps_cls = nsteps_in

      ifld_cls = ifld_cls_in
      ifld_clsr = ifld_clsr_in
      ifld_tls = ifld_tls_in
      ifld_tlsr = ifld_tlsr_in

      dt_cls = dt_in

      ifls_debug = ifdebug

      ntot = lx1*ly1*lz1*nelv
      !no natural BCs for LS fields
      !need to set tmasks to one
      !also turn off internal solvers
      if(ifld_cls_in.ne.0)then 
        call rone(tmask(1,1,1,1,ifld_cls_in-1),ntot)
        idpss(ifld_cls_in-1) = -1
      endif
      if(ifld_clsr_in.ne.0)then 
        call rone(tmask(1,1,1,1,ifld_clsr_in-1),ntot)
        idpss(ifld_clsr_in-1) = -1
      endif
      if(ifld_tls_in.ne.0)then 
        call rone(tmask(1,1,1,1,ifld_tls_in-1),ntot)
        idpss(ifld_tls_in-1) = -1
      endif
      if(ifld_tlsr_in.ne.0)then 
        call rone(tmask(1,1,1,1,ifld_tlsr_in-1),ntot)
        idpss(ifld_tlsr_in-1) = -1
      endif

      if(nio.eq.0)write(*,*)"Initialized Level-Set"
      
      if(nio.eq.0)write(*,*)"Debug mode",ifls_debug

      return
      end
C----------------------------------------------------------------------     
      subroutine ls_drive
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer istep_save, ifld_save
      real dt_save, time_save
      integer i,ntot

      !replace internal istep
      istep_save = ISTEP
      dt_save = dt
      time_save = time
      ifld_save = ifield

      ISTEP = 0
      dt = dt_cls
      time = 0.0
      ifield = ifld_clsr

      if(ifls_debug.eq.1 .and. nio.eq.0)then
        write(*,*) "Field", ifield
        write(*,*) "istep", istep
        write(*,*) "time step", dt_cls
        write(*,*) "Max iteration count", nsteps_cls
      endif

      if(ifls_debug.eq.1)
     $ call lsmonitor(t(1,1,1,1,ifld_tls-1),'TLS  ')

      if(ifls_debug.eq.1)
     $ call lsmonitor(t(1,1,1,1,ifld_clsr-1),'CLSr ')

      !Convert normal vector to rst space
      !Note that normals do not change over re-dist steps
      call cls_normals(clsnx,clsny,clsnz,ifld_tls)
      call vector_to_rst(clsnx,clsny,clsnz,
     $                   clsnr,clsns,clsnt)

      if(ifls_debug.eq.1)then
        ntot = lx1*ly1*lz1*nelv
        call copy(vx,clsnr,ntot)
        call copy(vy,clsns,ntot)
        if(if3d)call copy(vz,clsnt,ntot)
      endif
      if(ifls_debug.eq.1)call lsmonitor(clsnx,'xnorm')
      if(ifls_debug.eq.1)call lsmonitor(clsnr,'rnorm')

      do i=1,nsteps_cls
        istep = istep + 1
        call ls_advance
      enddo

      !replace back
      ISTEP = istep_save
      dt = dt_save
      time = time_save
      ifield = ifld_save

      return
      end
C----------------------------------------------------------------------     
      subroutine ls_advance
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer igeom

      call nekgsync()

      call settime_cls

      call setprop_cls

      if(ifls_debug.eq.1 .and. nio.eq.0)then 
        write(*,*)"ngeom: ",ngeom
      endif

      if (.not.iftmsh(ifield)) imesh = 1
      if (     iftmsh(ifield)) imesh = 2

      do igeom = 1,ngeom
        call unorm
        !diffusion array must be filled out here - pending
        call settolt
        call cdcls(igeom)
      enddo


      return
      end
C----------------------------------------------------------------------     
      subroutine settime_cls
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'LVLSET'

      integer irst, ilag

      irst = param(46)

      do ILAG=10,2,-1
        DTLAG(ILAG) = DTLAG(ILAG-1)
      enddo

      CALL setdt_cls
      DTLAG(1) = DT
      IF (ISTEP.EQ.1 .and. irst.le.0) DTLAG(2) = DT

      TIMEF    = TIME
      TIME     = TIME+DT

      CALL SETORDBD
      if (irst.gt.0) nbd = nbdinp
      CALL RZERO (BD,10)
      CALL SETBD (BD,DTLAG,NBD)
      if (PARAM(27).lt.0) then
        NAB = NBDINP
      else
        NAB = 3
      endif
      IF (ISTEP.lt.NAB.and.irst.le.0) NAB = ISTEP
      CALL RZERO   (AB,10)
      CALL SETABBD (AB,DTLAG,NAB,NBD)

      if(ifls_debug.eq.1 .and. nio.eq.0)then
        write(*,*)"BDF/EXT order",nbd,nab
      endif

      return
      end
C----------------------------------------------------------------------     
      subroutine setdt_cls
      implicit none
      include 'SIZE'
      include 'LVLSET'

      real umax

      call compute_cfl(umax,clsnx,clsny,clsnz,1)

      ! worry about adjust dt based on CFL later

      return
      end
C----------------------------------------------------------------------     
      subroutine cdcls(igeom)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
      include 'ORTHOT'

      common /lsscratch/ ta(lx1,ly1,lz1,lelt),
     $                   tb(lx1,ly1,lz1,lelt) 
      real ta,tb

      COMMON /SCRVH/ H1(LX1,LY1,LZ1,LELT),
     $               H2(lx1,ly1,lz1,lelt) 

      real h1,h2

      integer igeom,n,iter
      logical ifconv
      integer ifld1,isd,intype

      n = lx1*ly1*lz1*nelv

      ifld1 = ifield-1
      napproxt(1,ifld1) = laxtt


      if(ifls_debug.eq.1 .and. nio.eq.0)
     $  write(*,*)"in cdcls",igeom,ifls_debug

      if(igeom.eq.1)then
        call makeq_cls
        call lagscal
      else
        write(name4t,'(A4)')"CLSR"

        if(ifls_debug.eq.1)call lsmonitor(bq(1,1,1,1,ifield-1),'bq   ')

        isd = 1
        do iter=1,nmxnl
          intype = 0
          if(iftran) intype = -1
          call sethlm(h1,h2,intype)
          ! call bcneusc (ta,-1)
          ! call add2 (h2,ta,n)
          !following is divergence term
          call add2 (h2,adq(1,1,1,1,ifield-1),n)
          ! call bcdirsc (t(1,1,1,1,ifield-1))
          call axhelm_cls2(ta,t(1,1,1,1,ifield-1),h1,h2,imesh,isd) 
          call sub3(tb,bq(1,1,1,1,ifield-1),ta,n)
          ! call bcneusc (ta,1)
          ! call add2(tb,ta,n)
          
          if(ifls_debug.eq.1)call lsmonitor(tb,'tbrhs')

          call hsolve_cls(name4t,ta,tb,h1,h2,
     $                tmask(1,1,1,1,ifield-1),
     $                tmult(1,1,1,1,ifield-1),
     $                imesh,tolht(ifield),nmxt(ifield-1),1,
     $                approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)

          if(ifls_debug.eq.1) call lsmonitor(ta,'phidt')
          call add2(t(1,1,1,1,ifield-1),ta,n)
          call cvgnlps (ifconv)
          if (ifconv) exit
        enddo
      endif

      return
      end
C----------------------------------------------------------------------     
      subroutine axhelm_cls(au,u,helm1,helm2,imsh,isd)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      real au(lx1,ly1,lz1,1)
      real u(lx1,ly1,lz1,1)
      real helm1(lx1,ly1,lz1,1)
      real helm2(lx1,ly1,lz1,1)

      integer imsh,isd

      COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
     $  ,             DUDS  (LX1,LY1,LZ1)
     $  ,             DUDT  (LX1,LY1,LZ1)
     $  ,             TMP1  (LX1,LY1,LZ1)
     $  ,             TMP2  (LX1,LY1,LZ1)
     $  ,             TMP3  (LX1,LY1,LZ1)
      real dudr,duds,dudt,tmp1,tmp2,tmp3

      real tm1(lx1,ly1,lz1)
      real tm2(lx1,ly1,lz1)
      real tm3(lx1,ly1,lz1)
      equivalence (dudr,tm1),(duds,tm2),(dudt,tm3)

      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV

      integer ntot,nxy,nyz,nxz,nxyz
      integer e,iz

      common /ltmp/ ur(lx1,ly1,lz1),
     $              us(lx1,ly1,lz1),
     $              ut(lx1,ly1,lz1)
      real ur,us,ut

      NXY=lx1*ly1
      NYZ=ly1*lz1
      NXZ=lx1*lz1
      NXYZ=lx1*ly1*lz1
      NTOT=NXYZ*NELV

      call rzero(au,ntot)

      do e = 1,nelv
        if(.not.if3d)then
          call col3(ur,clsnr(1,1,1,e),u(1,1,1,e),NXYZ)
          call col3(us,clsns(1,1,1,e),u(1,1,1,e),NXYZ)
          call mxm  (dxm1,lx1,ur,lx1,dudr,nyz)
          call mxm  (us,lx1,dytm1,ly1,duds,ly1)
          call col3 (tmp1,dudr,g1m1(1,1,1,e),nxyz)
          call col3 (tmp2,duds,g2m1(1,1,1,e),nxyz)
          if (ifdfrm(e)) then
            call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
            call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
          endif
          call col2 (tmp1,helm1(1,1,1,e),nxyz)
          call col2 (tmp2,helm1(1,1,1,e),nxyz)
          call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
          call mxm  (tmp2,lx1,dym1,ly1,tm2,ly1)
          call col2(tm1,clsnr(1,1,1,e),NXYZ)
          call col2(tm2,clsns(1,1,1,e),NXYZ)
          call add2 (au(1,1,1,e),tm1,nxyz)
          call add2 (au(1,1,1,e),tm2,nxyz)
        else
          call col3(ur,clsnr(1,1,1,e),u(1,1,1,e),NXYZ)
          call col3(us,clsns(1,1,1,e),u(1,1,1,e),NXYZ)
          call col3(ut,clsnt(1,1,1,e),u(1,1,1,e),NXYZ)
          call mxm(dxm1,lx1,ur,lx1,dudr,nyz)
          do iz=1,lz1
            call mxm(us,lx1,dytm1,ly1,duds(1,1,iz),ly1)
          enddo
          call mxm     (ut,nxy,dztm1,lz1,dudt,lz1)
          call col3    (tmp1,dudr,g1m1(1,1,1,e),nxyz)
          call col3    (tmp2,duds,g2m1(1,1,1,e),nxyz)
          call col3    (tmp3,dudt,g3m1(1,1,1,e),nxyz)
          if (ifdfrm(e)) then
            call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
            call addcol3 (tmp1,dudt,g5m1(1,1,1,e),nxyz)
            call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
            call addcol3 (tmp2,dudt,g6m1(1,1,1,e),nxyz)
            call addcol3 (tmp3,dudr,g5m1(1,1,1,e),nxyz)
            call addcol3 (tmp3,duds,g6m1(1,1,1,e),nxyz)
          endif
          call col2 (tmp1,helm1(1,1,1,e),nxyz)
          call col2 (tmp2,helm1(1,1,1,e),nxyz)
          call col2 (tmp3,helm1(1,1,1,e),nxyz)
          call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
          do iz=1,lz1
            call mxm(tmp2(1,1,iz),lx1,dym1,ly1,tm2(1,1,iz),ly1)
          enddo
          call mxm  (tmp3,nxy,dzm1,lz1,tm3,lz1)
          call col2(tm1,clsnr(1,1,1,e),NXYZ)
          call col2(tm2,clsns(1,1,1,e),NXYZ)
          call col2(tm3,clsnt(1,1,1,e),NXYZ)
          call add2 (au(1,1,1,e),tm1,nxyz)
          call add2 (au(1,1,1,e),tm2,nxyz)
          call add2 (au(1,1,1,e),tm3,nxyz)
        endif
      enddo


      call addcol4 (au,helm2,bm1,u,ntot)

      if(ifsvv(ifield-1))call axhelm_svv(au,u,imsh,isd)
      !lets worry about axisymmetry later

      if(ifls_debug.eq.1 .and. nio.eq.0)
     $ write(*,*)"SVV status",ifsvv(ifield-1)
      if(ifls_debug.eq.1) call lsmonitor(au,'Diff ')

      return
      end
C----------------------------------------------------------------------     
      subroutine setprop_cls
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer n,i
      real deltael

      n = lx1*ly1*lz1*lelv

      call cfill(VTRANS(1,1,1,1,ifield),1.0,n)

      do i=1,n
        VDIFF(i,1,1,1,ifield) = deltael(i,1,1,1) * eps_cls
      enddo

      if(ifls_debug.eq.1)then
        call lsmonitor(vtrans(1,1,1,1,ifield),'rho  ')
        call lsmonitor(vdiff(1,1,1,1,ifield),'diff ')
      endif

      return
      end
C----------------------------------------------------------------------     
      subroutine vector_to_rst(ux,uy,uz,ur,us,ut)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real ux(1),uy(1),uz(1)
      real ur(1),us(1),ut(1)

      integer i,ntot
      real xr,yr,zr,xs,ys,zs,xt,yt,zt

      integer icalld
      save icalld
      data icalld /0/

      common /lvlset_x/ xrm1(lx1,ly1,lz1,lelt),
     $                  yrm1(lx1,ly1,lz1,lelt),
     $                  zrm1(lx1,ly1,lz1,lelt),
     $                  xsm1(lx1,ly1,lz1,lelt),
     $                  ysm1(lx1,ly1,lz1,lelt),
     $                  zsm1(lx1,ly1,lz1,lelt),
     $                  xtm1(lx1,ly1,lz1,lelt),
     $                  ytm1(lx1,ly1,lz1,lelt),
     $                  ztm1(lx1,ly1,lz1,lelt)
      real xrm1,yrm1,zrm1 
      real xsm1,ysm1,zsm1
      real xtm1,ytm1,ztm1

      if(icalld.eq.0)then
        call XYZRST(XRM1,YRM1,ZRM1,
     $               XSM1,YSM1,ZSM1,
     $                XTM1,YTM1,ZTM1,ifaxis)
        icalld = 1
      endif
      ntot = lx1*ly1*lz1*nelv

      if(if3d)then
        do i=1,ntot
          xr = xrm1(i,1,1,1)
          yr = yrm1(i,1,1,1)
          zr = zrm1(i,1,1,1)
          xs = xsm1(i,1,1,1)
          ys = ysm1(i,1,1,1)
          zs = zsm1(i,1,1,1)
          xt = xtm1(i,1,1,1)
          yt = ytm1(i,1,1,1)
          zt = ztm1(i,1,1,1)
          ur(i) = xr*ux(i) + yr*uy(i) + zr*uz(i)
          us(i) = xs*ux(i) + ys*uy(i) + zs*uz(i)
          ut(i) = xt*ux(i) + yt*uy(i) + zt*uz(i)
        enddo
      else
        do i=1,ntot
          xr = xrm1(i,1,1,1)
          yr = yrm1(i,1,1,1)
          xs = xsm1(i,1,1,1)
          ys = ysm1(i,1,1,1)
          ur(i) = xr*ux(i) + yr*uy(i)
          us(i) = xs*ux(i) + ys*uy(i)
        enddo
      endif

      return
      end
C----------------------------------------------------------------------     
      subroutine makeq_cls
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer ntot

      common /lsscratch/ ta(lx1,ly1,lz1,lelt),
     $                   tb(lx1,ly1,lz1,lelt) 
      real ta,tb

      integer i

      ntot = lx1*ly1*lz1*nelv

      !this is to add divergence. Currently acts as dummy
      call makeq_aux

      !(1-psi)*psi
      call copy(tb,t(1,1,1,1,ifield-1),ntot)
      call cmult(tb,-1.0,ntot)
      call cadd(tb,1.0,ntot)
      call col2(tb,t(1,1,1,1,ifield-1),ntot)

      if(ifls_debug.eq.1) call lsmonitor(tb,'convp')

      !convop
      call rzero(ta,ntot)
      call convect_cons(ta,tb,.false.,clsnx,clsny,clsnz,.false.)
      call invcol2(ta,bm1,ntot)

      if(ifls_debug.eq.1) call lsmonitor(ta,'advec')

      !convab
      do i=1,ntot
        bq(i,1,1,1,ifield-1) = -bm1(i,1,1,1)*ta(i,1,1,1)
      enddo

      if(ifls_debug.eq.1) call lsmonitor(bq(1,1,1,1,ifield-1),'bqarr')

      call makeabq

      if(ifls_debug.eq.1) call lsmonitor(bq(1,1,1,1,ifield-1),'bqabq')

      call makebdq

      if(ifls_debug.eq.1)call lsmonitor(bq(1,1,1,1,ifield-1),'bqmke')

      return
      end
c---------------------------------------------------------------
      real function deltael(ix,iy,iz,iel)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer ix,iy,iz,iel

      real dx(lx1,ly1,lz1,lelt)
      save dx 

      integer icalld
      data icalld /0/
      save icalld

      integer nxyz,n,ie
      real dd,dinv
      real dxmin_e
      real dxmax_e

      nxyz = nx1*ny1*nz1
      n    = nxyz*nelv

      if (icalld.eq.0 .or. ifmvbd) then
         dinv = 1./ldim
         do ie = 1,nelv
            dd = dxmax_e(ie)
            call cfill(dx(1,1,1,ie),dd,nxyz) 
         enddo
         icalld = 1
      endif

      deltael = dx(ix,iy,iz,iel)

      return
      end 
c---------------------------------------------------------------
      real function heaviside(ix,iy,iz,iel,phi)
      include 'SIZE'
      include 'LVLSET'

      integer ix,iy,iz,iel

      real eps, deltael, phi

      eps = deltael(ix,iy,iz,iel)*eps_cls
      heaviside = 0.5*(tanh(phi/(2.0*eps))+1.0)

      return
      end
c---------------------------------------------------------------
      subroutine cls_normals(cnx,cny,cnz,ifld)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real cnx(lx1,ly1,lz1,1)
      real cny(lx1,ly1,lz1,1)
      real cnz(lx1,ly1,lz1,1)

      real cmag(lx1,ly1,lz1,lelv)

      integer ntot,ifld,i

      ntot = lx1*ly1*lz1*nelv

      !must be calc from TLS field
      call gradm1(cnx,cny,cnz,t(1,1,1,1,ifld-1))
      call opcolv(cnx,cny,cnz,bm1)
      call opdssum(cnx,cny,cnz)
      call opcolv(cnx,cny,cnz,binvm1)

      call col3(cmag,cnx,cnx,ntot)
      call addcol3(cmag,cny,cny,ntot)
      if(if3d) call addcol3(cmag,cnz,cnz,ntot)
      call vsqrt(cmag,ntot)

      do i=1,ntot
        if(cmag(i,1,1,1).gt.0.)then
          cnx(i,1,1,1) = cnx(i,1,1,1)/cmag(i,1,1,1)
          cny(i,1,1,1) = cny(i,1,1,1)/cmag(i,1,1,1)
          if(if3d)cnx(i,1,1,1) = cnx(i,1,1,1)/cmag(i,1,1,1)
        endif
      enddo

      return
      end
c---------------------------------------------------------------
      subroutine lsmonitor(u,aname)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real u(1)
      character*5 aname
      
      integer n

      real norm, amin, amax
      real gl2norm, glmax, glmin

      n = lx1*ly1*lz1*nelv

      norm = gl2norm(u,n)
      
      amin = glmin(u,n)

      amax = glmax(u,n)

      if(nio.eq.0)then
        write(6,1000)aname," norm, min, max:",
     $   norm,amin,amax
      endif

1000  format(a,10x,a,1p3E13.4)

      return
      end
c---------------------------------------------------------------
      subroutine axhelm_cls2(au,u,helm1,helm2,imsh,isd)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      real au(lx1,ly1,lz1,1)
      real u(lx1,ly1,lz1,1)
      real helm1(lx1,ly1,lz1,1)
      real helm2(lx1,ly1,lz1,1)

      integer imsh,isd

      COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
     $  ,             DUDS  (LX1,LY1,LZ1)
     $  ,             DUDT  (LX1,LY1,LZ1)
     $  ,             TMP1  (LX1,LY1,LZ1)
     $  ,             TMP2  (LX1,LY1,LZ1)
     $  ,             TMP3  (LX1,LY1,LZ1)
      real dudr,duds,dudt,tmp1,tmp2,tmp3

      real tm1(lx1,ly1,lz1)
      real tm2(lx1,ly1,lz1)
      real tm3(lx1,ly1,lz1)
      equivalence (dudr,tm1),(duds,tm2),(dudt,tm3)

      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV

      integer ntot,nxy,nyz,nxz,nxyz
      integer ie,iz

      common /ltmp/ ur(lx1,ly1,lz1),
     $              us(lx1,ly1,lz1),
     $              ut(lx1,ly1,lz1)
      real ur,us,ut

      real tmp(lx1,ly1,lz1,lelt)
      real tmpx(lx1,ly1,lz1,lelt)
      real tmpy(lx1,ly1,lz1,lelt)
      real tmpz(lx1,ly1,lz1,lelt)

      NXY=lx1*ly1
      NYZ=ly1*lz1
      NXZ=lx1*lz1
      NXYZ=lx1*ly1*lz1
      NTOT=NXYZ*NELV

      call rzero(au,ntot)

      call gradm1(tmpx,tmpy,tmpz,u)
      call opcolv(tmpx,tmpy,tmpz,bm1)
      call opdssum(tmpx,tmpy,tmpz)
      call opcolv(tmpx,tmpy,tmpz,binvm1)

      if(if3d)then
         call vdot3(tmp,tmpx,tmpy,tmpz,clsnx,clsny,clsnz,ntot)
      else
         call vdot2(tmp,tmpx,tmpy,clsnx,clsny,ntot)
      endif

      do ie = 1,nelv
        call vdot2(tmp1,rxm1(1,1,1,ie),rym1(1,1,1,ie),
     $           clsnx(1,1,1,ie),clsny(1,1,1,ie),nxyz)
        call vdot2(tmp2,sxm1(1,1,1,ie),sym1(1,1,1,ie),
     $           clsnx(1,1,1,ie),clsny(1,1,1,ie),nxyz)
        call col2(tmp1,w3m1,nxyz)
        call col2(tmp2,w3m1,nxyz)

        call col2(tmp1,tmp(1,1,1,ie),nxyz)
        call col2(tmp2,tmp(1,1,1,ie),nxyz)

        call col2(tmp1,helm1(1,1,1,ie),nxyz)
        call col2(tmp2,helm1(1,1,1,ie),nxyz)

        call mxm(dxtm1,lx1,tmp1,lx1,tm1,nyz)
        call mxm(tmp2,lx1,dym1,ly1,tm2,ly1)

        call add2(au(1,1,1,ie),tm1,nxyz)
        call add2(au(1,1,1,ie),tm2,nxyz)
      enddo


      call addcol4 (au,helm2,bm1,u,ntot)

      if(ifsvv(ifield-1))call axhelm_svv(au,u,imsh,isd)
      !lets worry about axisymmetry later

      if(ifls_debug.eq.1 .and. nio.eq.0)
     $ write(*,*)"SVV status",ifsvv(ifield-1)
      if(ifls_debug.eq.1) call lsmonitor(au,'Diff ')

      return
      end
