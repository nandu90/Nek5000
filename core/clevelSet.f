c---------------------------------------------------------------------
      real function heaviphase(ix,iy,iz,ie,phi)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      real eps
      
      eps = lsH(ix,iy,iz,ie)*eps_cls
      heaviphase = 0.5*(tanh(2.0*PI*phi/eps)+1.0)
      
      end
c---------------------------------------------------------------------
      subroutine clsconv(ifld)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
            
      real tmp(lx1,ly1,lz1,lelv)
      
      common /clstemp2/ tmp

      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv

c     psi(1-psi)
      call copy(tmp,t(1,1,1,1,ifld-1),ntot)
      call cmult(tmp,-1.,ntot)
      call cadd(tmp,1.,ntot)
      call col2(tmp,t(1,1,1,1,ifld-1),ntot)

      call convect_cons(clsadv,tmp,.false.,clsnx,clsny,clsnz,.false.)
c      call col2(clsadv,bm1,ntot)
      call invcol2(clsadv,bm1,ntot)
     
      return
      end
c---------------------------------------------------------------------
      real function q_clsconv(ix,iy,iz,ie)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      q_clsconv = -clsadv(ix,iy,iz,ie)
      
      end
c---------------------------------------------------------------------      
      subroutine clsaxhelm(ifld)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
      
      real tmpx(lx1,ly1,lz1,lelv)
      real tmpy(lx1,ly1,lz1,lelv)
      real tmpz(lx1,ly1,lz1,lelv)
      real tmp(lx1,ly1,lz1,lelv)

      real tmp1(lx1,ly1,lz1)
      real tmp2(lx1,ly1,lz1)
      real tmp3(lx1,ly1,lz1)
      real tm1(lx1,ly1,lz1)
      real tm2(lx1,ly1,lz1)
      real tm3(lx1,ly1,lz1)
      real eps(lx1,ly1,lz1)
      
      common /clstemp/ tmpx,tmpy,tmpz,tmp,
     $     tm1,tm2,tm3,tmp1,tmp2,tmp3

      integer e

c     From axhelm - needed for ifdfrm
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
      
      nxy = lx1*ly1
      nyz = ly1*lz1
      nxz = lx1*lz1
      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv

c     \nabla \psi
      call gradm1(tmpx,tmpy,tmpz,t(1,1,1,1,ifld-1))
      call opcolv(tmpx,tmpy,tmpz,bm1)
      call opdssum(tmpx,tmpy,tmpz)
      call opcolv(tmpx,tmpy,tmpz,binvm1)

c     \nabla \psi . n
      if(if3d)then
         call vdot3(tmp,tmpx,tmpy,tmpz,clsnx,clsny,clsnz,ntot)
      else
         call vdot2(tmp,tmpx,tmpy,clsnx,clsny,ntot)
      endif

c     Get the normal vector in rst space
c     This should be moved out of this routine - pending
c      call vec_xyz2rst(clsnx,clsny,clsnz,clsnr,clsns,clsnt)
c      call copy(vx,clsnr,ntot)
c      call copy(vy,clsns,ntot)
      
      call rzero(clsau,ntot)
            
c     The following is almost similar to axhelm
      do e=1,nelv
         if(ldim.eq.2)then
            call copy(tm1,clsnx(1,1,1,e),nxyz)
            call copy(tm2,clsny(1,1,1,e),nxyz)

            call col3(tmp1,tm1,g1m1(1,1,1,e),nxyz)
            call col3(tmp2,tm2,g2m1(1,1,1,e),nxyz)
            if (ifdfrm(e)) then
               call addcol3 (tmp1,tm2,g4m1(1,1,1,e),nxyz)
               call addcol3 (tmp2,tm1,g4m1(1,1,1,e),nxyz)
            endif

            call copy(eps,lsH(1,1,1,e),nxyz)
            call cmult(eps,eps_cls,nxyz)
            call col2(eps,tmp(1,1,1,e),nxyz)

            call col2(tmp1,eps,nxyz)
            call col2(tmp2,eps,nxyz)

            call mxm(dxtm1,lx1,tmp1,lx1,tm1,nyz)
            call mxm(tmp2,lx1,dym1,ly1,tm2,ly1)

            call add2 (clsau(1,1,1,e),tm1,nxyz)
            call add2 (clsau(1,1,1,e),tm2,nxyz)
         else
            call copy(tm1,clsnr(1,1,1,e),nxyz)
            call copy(tm2,clsns(1,1,1,e),nxyz)
            call copy(tm3,clsnt(1,1,1,e),nxyz)

            call col3(tmp1,tm1,g1m1(1,1,1,e),nxyz)
            call col3(tmp2,tm2,g2m1(1,1,1,e),nxyz)
            call col3(tmp3,tm3,g3m1(1,1,1,e),nxyz)
            if (ifdfrm(e)) then
               call addcol3 (tmp1,tm2,g4m1(1,1,1,e),nxyz)
               call addcol3 (tmp1,tm3,g5m1(1,1,1,e),nxyz)
               call addcol3 (tmp2,tm1,g4m1(1,1,1,e),nxyz)
               call addcol3 (tmp2,tm3,g6m1(1,1,1,e),nxyz)
               call addcol3 (tmp3,tm1,g5m1(1,1,1,e),nxyz)
               call addcol3 (tmp3,tm2,g6m1(1,1,1,e),nxyz)
            endif

            call copy(eps,lsH(1,1,1,e),nxyz)
            call cmult(eps,eps_cls,nxyz)
            call col2(eps,tmp(1,1,1,e),nxyz)

            call col2(tmp1,eps,nxyz)
            call col2(tmp2,eps,nxyz)
            call col2(tmp3,eps,nxyz)

            call mxm(dxtm1,lx1,tmp1,lx1,tm1,nyz)
            do iz=1,lz1
               call mxm(tmp2(1,1,iz),lx1,dym1,ly1,tm2(1,1,iz),ly1)
            enddo
            call mxm(tmp3,nxy,dzm1,lz1,tm3,lz1)
            
            call add2 (clsau(1,1,1,e),tm1,nxyz)
            call add2 (clsau(1,1,1,e),tm2,nxyz)
            call add2 (clsau(1,1,1,e),tm3,nxyz)
         endif
      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine vec_xyz2rst(ux,uy,uz,ur,us,ut)
c
      include 'SIZE'
      include 'TOTAL'

      real ux(1),uy(1),uz(1)
      real ur(1),us(1),ut(1)
      
      COMMON /CLS_SCRNS/ XRM1(LX1,LY1,LZ1,LELT)
     $     ,             YRM1(LX1,LY1,LZ1,LELT)
     $     ,             XSM1(LX1,LY1,LZ1,LELT)
     $     ,             YSM1(LX1,LY1,LZ1,LELT)
     $     ,             XTM1(LX1,LY1,LZ1,LELT)
     $     ,             YTM1(LX1,LY1,LZ1,LELT)
     $     ,             ZRM1(LX1,LY1,LZ1,LELT)
     $     ,             ZSM1(LX1,LY1,LZ1,LELT)
     $     ,             ZTM1(LX1,LY1,LZ1,LELT)

      integer icalld
      save icalld
      data icalld /0/

      if(icalld.eq.0)then
         CALL XYZRST(XRM1,YRM1,ZRM1,XSM1,YSM1,ZSM1,XTM1,YTM1,ZTM1,
     $        IFAXIS)
         icalld = 1
      endif
      

      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv

      do i = 1,ntot
         if(ldim .eq. 2)then
            ur(i) = ux(i)*xrm1(i,1,1,1)+uy(i)*yrm1(i,1,1,1)
            us(i) = ux(i)*xsm1(i,1,1,1)+uy(i)*ysm1(i,1,1,1)
         else
            ur(i) = ux(i)*xrm1(i,1,1,1)+uy(i)*yrm1(i,1,1,1)
     $           +uz(i)*zrm1(i,1,1,1)
            us(i) = ux(i)*xsm1(i,1,1,1)+uy(i)*ysm1(i,1,1,1)
     $           +uz(i)*zsm1(i,1,1,1)
            ut(i) = ux(i)*xtm1(i,1,1,1)+uy(i)*ytm1(i,1,1,1)
     $           +uz(i)*ztm1(i,1,1,1)
         endif
      enddo

      return
      end
c---------------------------------------------------------------------         
