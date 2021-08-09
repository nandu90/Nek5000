c-----------------------------------------------------------------------
      real function find_root_utau(u,S,visc)
      implicit real(a-h,o-z)
      parameter (n=100)

      external func_utau
      external fder_utau

      tol = 1.e-8
      utauo = 0.1
      do i=1,n

        func = func_utau(utauo,u,S,visc)
        fder = fder_utau(utauo,u,S,visc)
        utaun = utauo - func/fder
c        write(*,'(I6, 5(1X,G22.15))') i, utaun, utauo, func, fder
        if(dabs(utaun-utauo).le.tol) goto 10
        utauo = utaun

      enddo
      write(*,*) 'Could not converge for utau ', utauo, u, S
      call exitt

 10   continue
c      write(*,*) 'utau is', utaun
      find_root_utau = utaun

      return
      end
c-----------------------------------------------------------------------
      real function func_utau(utau,u,S,visc) 
      implicit real(a-h,o-z)
      real kap, S, E

      kap = 0.4
      E   = 9.0

      func_utau = utau*utau - (kap*visc*S/E)*exp(kap*u/utau) 

      return
      end
c-----------------------------------------------------------------------
      real function fder_utau(utau,u,S,visc) 
      implicit real(a-h,o-z)
      real kap, S, E

      kap = 0.4
      E   = 9.0

      fder_utau = 2.0*utau + 
     $             (kap*kap*visc*u*S/E)*exp(kap*u/utau)/utau/utau

      return
      end
c-----------------------------------------------------------------------
      subroutine getangent(st,ix,iy,iz,iside,e)

c     calculate surface normal

      include 'SIZE'
      include 'GEOM'
      include 'TOPOL'

      real st(3)
      integer e,f

      f = eface1(iside)

      if (1.le.f.and.f.le.2) then     ! "r face"
         st(1) = t1x(iy,iz,iside,e)
         st(2) = t1y(iy,iz,iside,e)
         st(3) = t1z(iy,iz,iside,e)
      elseif (3.le.f.and.f.le.4) then ! "s face"
         st(1) = t1x(ix,iz,iside,e)
         st(2) = t1y(ix,iz,iside,e)
         st(3) = t1z(ix,iz,iside,e)
      elseif (5.le.f.and.f.le.6) then ! "t face"
         st(1) = t1x(ix,iy,iside,e)
         st(2) = t1y(ix,iy,iside,e)
         st(3) = t1z(ix,iy,iside,e)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine getbitangent(sb,ix,iy,iz,iside,e)

c     calculate surface normal

      include 'SIZE'
      include 'GEOM'
      include 'TOPOL'

      real sb(3)
      integer e,f

      f = eface1(iside)

      if (1.le.f.and.f.le.2) then     ! "r face"
         sb(1) = t2x(iy,iz,iside,e)
         sb(2) = t2y(iy,iz,iside,e)
         sb(3) = t2z(iy,iz,iside,e)
      elseif (3.le.f.and.f.le.4) then ! "s face"
         sb(1) = t2x(ix,iz,iside,e)
         sb(2) = t2y(ix,iz,iside,e)
         sb(3) = t2z(ix,iz,iside,e)
      elseif (5.le.f.and.f.le.6) then ! "t face"
         sb(1) = t2x(ix,iy,iside,e)
         sb(2) = t2y(ix,iy,iside,e)
         sb(3) = t2z(ix,iy,iside,e)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine fixmask_old(c1mask,c2mask,c3mask)

c     fixes masks for SYM face corners

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'

      real   c1mask(lx1,ly1,lz1,1)
     $      ,c2mask(lx1,ly1,lz1,1)
     $      ,c3mask(lx1,ly1,lz1,1)

      common /ctmp0/ im1(lx1,ly1,lz1),im2(lx1,ly1,lz1)
      integer e,f,val,im1,im2

      character*3 cb

      n = lx1*ly1*lz1

      nface = 2*ldim

      do e=1,nelv
         call izero (im1,n)
         call izero (im2,n)
         do f=1,nface
            cb  = cbc (f,e,1)
            if (cb.eq.'SYM')  call add_iface_e(im1,f,1,lx1,ly1,lz1)
c            if (cb.eq.'SYM')  call iface_e(im2,f,2,lx1,ly1,lz1)
         enddo
c         call icol2(im2,im1,n)

         tolr = 1.e-8
         k = 1
         do j=1,ly1,ly1-1
         do i=1,lx1,lx1-1
            xx = xm1(i,j,k,e) 
            yy = ym1(i,j,k,e) 
            rr = sqrt(xx*xx+yy*yy)
c            if  ( im2(i,j,k) .eq. 2) then  ! corner of SYM & 'SYM' faces
            if  ( im1(i,j,k) .eq. 2) then  ! corner of SYM & 'SYM' faces
               c1mask(i,j,k,e) = 0.
               c2mask(i,j,k,e) = 0.
            elseif(rr.le.tolr) then
               c1mask(i,j,k,e) = 0.
               c2mask(i,j,k,e) = 0.
            endif
         enddo
         enddo
      enddo

c      do e=1,nelv
c         if ( v1mask(1,1,1,e).eq.0 .and.
c     $        v1mask(2,1,1,e).eq.0 .and.
c     $        v1mask(1,2,1,e).eq.0 ) v2mask(1,1,1,e)=0.
c
c         if ( v1mask(nx1  ,1,1,e).eq.0 .and.
c     $        v1mask(nx1-1,1,1,e).eq.0 .and.
c     $        v1mask(nx1  ,2,1,e).eq.0 ) v2mask(nx1,1,1,e)=0.
c
c         if ( v1mask(1,ny1  ,1,e).eq.0 .and.
c     $        v1mask(1,ny1-1,1,e).eq.0 .and.
c     $        v1mask(2,ny1  ,1,e).eq.0 ) v2mask(1,ny1,1,e)=0.
c
c         if ( v1mask(nx1,ny1  ,1,e).eq.0 .and.
c     $        v1mask(nx1,ny1-1,1,e).eq.0 .and.
c     $        v1mask(nx1-1,ny1,1,e).eq.0 ) v2mask(nx1,ny1,1,e)=0.
c      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine add_iface_e(a,iface,val,nx,ny,nz)

C     Assign the value VAL to face(IFACE,IE) of array A.
C     IFACE is the input in the pre-processor ordering scheme.

      include 'SIZE'
      integer a(nx,ny,nz),val
      call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)
      do 100 iz=kz1,kz2
      do 100 iy=ky1,ky2
      do 100 ix=kx1,kx2
         a(ix,iy,iz)=a(ix,iy,iz)+val
  100 continue
      return
      end

c-----------------------------------------------------------------------
      subroutine un_face_e(u,iface,val,nx,ny,nz)

C     Assign the value VAL to face(IFACE,IE) of array A.
C     IFACE is the input in the pre-processor ordering scheme.

      include 'SIZE'
      real a(nx,ny,nz),val
      call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)
      do 100 iz=kz1,kz2
      do 100 iy=ky1,ky2
      do 100 ix=kx1,kx2
         a(ix,iy,iz)=a(ix,iy,iz)+val
  100 continue
      return
      end

c-----------------------------------------------------------------------
      subroutine get_gradywd
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      common /gradywd/ ywdx(lx1,ly1,lz1,lelv),ywdy(lx1,ly1,lz1,lelv)
     $                ,ywdz(lx1,ly1,lz1,lelv),ywdc(lx1,ly1,lz1,lelv)

      ntot = lx1*ly1*lz1*lelv

      call gradm1 (ywdx,ywdy,ywdz,   ywd)
      call opcolv (ywdx,ywdy,ywdz,   bm1)
      call opdssum(ywdx,ywdy,ywdz)
      call opcolv (ywdx,ywdy,ywdz,binvm1)
      call copy   (ywdc,ywd ,ntot)

      return
      end

c-----------------------------------------------------------------------
      subroutine fixmask(cbtype,c1mask,c2mask,c3mask)

c     fixes masks for SYM face corners

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'PARALLEL'

      real   c1mask(lx1,ly1,lz1,1)
     $      ,c2mask(lx1,ly1,lz1,1)
     $      ,c3mask(lx1,ly1,lz1,1)

      real   arr(lx1,ly1,lz1,lelv),ar1(lx1,ly1,lz1,lelv)
     $      ,ar2(lx1,ly1,lz1,lelv),ar3(lx1,ly1,lz1,lelv)

      character*3 cb,cbtype

      nxyz1= lx1*ly1*lz1
      ntot1= nxyz1*nelv
      nfaces = 2*ldim
      tol  = 1.e-07

      call rzero  (arr,  ntot1)
      call oprzero(ar1,ar2,ar3)

      write(*,*) 'element faces from fixmask'
      do 1000 iel=1,nelv
      ieg = lglel(iel)
      do 100 iface=1,nfaces
         cb = cbc(iface,iel,1)
         if (cb.eq.cbtype) then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
            ia = 0
            do 10 iz=kz1,kz2
            do 10 iy=ky1,ky2
            do 10 ix=kx1,kx2
               ia =ia + 1
               arr(ix,iy,iz,iel)=arr(ix,iy,iz,iel)+1.
               ar1(ix,iy,iz,iel)=ar1(ix,iy,iz,iel)+unx(ia,1,iface,iel)
               ar2(ix,iy,iz,iel)=ar2(ix,iy,iz,iel)+uny(ia,1,iface,iel)
               if(if3d)
     $         ar3(ix,iy,iz,iel)=ar3(ix,iy,iz,iel)+unz(ia,1,iface,iel)
 10         continue
         endif
 100  continue
 1000 continue

      call dssum  (arr,nx1,ny1,nz1)
      call opdssum(ar1,  ar2,  ar3)

      do 2000 iel=1,nelv
      ieg = lglel(iel)
      do 200 iface=1,nfaces
         cb = cbc(iface,iel,1)
         if (cb.eq.cbtype) then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
            ia = 0
            do 20 iz=kz1,kz2
            do 20 iy=ky1,ky2
            do 20 ix=kx1,kx2
               ia =ia + 1
               amul = arr(ix,iy,iz,iel)
               ar1s = ar1(ix,iy,iz,iel)/amul
               ar2s = ar2(ix,iy,iz,iel)/amul
               ar3s = 0.0

               if(if3d)
     $         ar3s = ar3(ix,iy,iz,iel)/amul
               unmag = sqrt(ar1s*ar1s+ar2s*ar2s+ar3s*ar3s)

               if((1.0-abs(unmag)).ge.tol) then
                 write(*,'(A,3I5,2(2X,G14.7))') 'normal vector '
     $                             , ieg, iface, ia, unmag, amul
                 if(amul.eq.2.) then
                    c1mask(ix,iy,iz,iel) = 0.
                    c2mask(ix,iy,iz,iel) = 0.
                 endif
                 if(amul.eq.3.) then
                    c1mask(ix,iy,iz,iel) = 0.
                    c2mask(ix,iy,iz,iel) = 0.
                    c3mask(ix,iy,iz,iel) = 0.
                 endif
               endif
 20         continue
         endif
 200  continue
 2000 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine fixmask2(cbtype1,cbtype2,c1mask,c2mask,c3mask)

c     fixes masks for SYM face corners

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'PARALLEL'

      real   c1mask(lx1,ly1,lz1,1)
     $      ,c2mask(lx1,ly1,lz1,1)
     $      ,c3mask(lx1,ly1,lz1,1)

      common /scruz/ rmlt(lx1,ly1,lz1,lelv),runx(lx1,ly1,lz1,lelv)
     $              ,runy(lx1,ly1,lz1,lelv),runz(lx1,ly1,lz1,lelv)
     $              ,rt1x(lx1,ly1,lz1,lelv),rt1y(lx1,ly1,lz1,lelv)
     $              ,rt1z(lx1,ly1,lz1,lelv),rt2x(lx1,ly1,lz1,lelv)
     $              ,rt2y(lx1,ly1,lz1,lelv),rt2z(lx1,ly1,lz1,lelv)

      character*3 cb,cbtype1,cbtype2

      nxyz1= lx1*ly1*lz1
      ntot1= nxyz1*nelv
      nfaces = 2*ldim
      tol  = 1.e-07

      call rzero  (rmlt,    ntot1)
      call oprzero(runx,runy,runz)
      call oprzero(rt1x,rt1y,rt1z)
      call oprzero(rt2x,rt2y,rt2z)

      write(*,*) 'element faces from fixmask2'
      do 1000 iel=1,nelv
      ieg = lglel(iel)
      do 100 iface=1,nfaces
         cb = cbc(iface,iel,1)
         if (cb.eq.cbtype1 .or. cb.eq.cbtype2) then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
            ia = 0
            do 10 iz=kz1,kz2
            do 10 iy=ky1,ky2
            do 10 ix=kx1,kx2
              ia =ia + 1
              rmlt(ix,iy,iz,iel)=rmlt(ix,iy,iz,iel)+1.
              runx(ix,iy,iz,iel)=runx(ix,iy,iz,iel)+unx(ia,1,iface,iel)
              runy(ix,iy,iz,iel)=runy(ix,iy,iz,iel)+uny(ia,1,iface,iel)
              rt1x(ix,iy,iz,iel)=rt1x(ix,iy,iz,iel)+t1x(ia,1,iface,iel)
              rt1y(ix,iy,iz,iel)=rt1y(ix,iy,iz,iel)+t1y(ia,1,iface,iel)
              rt2x(ix,iy,iz,iel)=rt2x(ix,iy,iz,iel)+t2x(ia,1,iface,iel)
              rt2y(ix,iy,iz,iel)=rt2y(ix,iy,iz,iel)+t2y(ia,1,iface,iel)
              if(if3d) then
               runz(ix,iy,iz,iel)=runz(ix,iy,iz,iel)+unz(ia,1,iface,iel)
               rt1z(ix,iy,iz,iel)=rt1z(ix,iy,iz,iel)+t1z(ia,1,iface,iel)
               rt2z(ix,iy,iz,iel)=rt2z(ix,iy,iz,iel)+t2z(ia,1,iface,iel)
              endif
 10         continue
         endif
 100  continue
 1000 continue

      call dssum  (rmlt,nx1,ny1,nz1)
      call opdssum(runx, runy, runz)
      call opdssum(rt1x, rt1y, rt1z)
      call opdssum(rt2x, rt2y, rt2z)

      do 2000 iel=1,nelv
      ieg = lglel(iel)
      do 200 iface=1,nfaces
         cb = cbc(iface,iel,1)
         if (cb.eq.cbtype1 .or. cb.eq.cbtype2) then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
            ia = 0
            do 20 iz=kz1,kz2
            do 20 iy=ky1,ky2
            do 20 ix=kx1,kx2
               ia =ia + 1
               amul = rmlt(ix,iy,iz,iel)
               runxs= runx(ix,iy,iz,iel)/amul
               runys= runy(ix,iy,iz,iel)/amul
               runzs= 0.0
               rt1xs= rt1x(ix,iy,iz,iel)/amul
               rt1ys= rt1y(ix,iy,iz,iel)/amul
               rt1zs= 0.0
               rt2xs= rt2x(ix,iy,iz,iel)/amul
               rt2ys= rt2y(ix,iy,iz,iel)/amul
               rt2zs= 0.0

               if(if3d) then
               runzs= runz(ix,iy,iz,iel)/amul
               rt1zs= rt1z(ix,iy,iz,iel)/amul
               rt2zs= rt2z(ix,iy,iz,iel)/amul
               endif
               unmag = sqrt(runxs*runxs+runys*runys+runzs*runzs)
               t1mag = sqrt(rt1xs*rt1xs+rt1ys*rt1ys+rt1zs*rt1zs)
               t2mag = sqrt(rt2xs*rt2xs+rt2ys*rt2ys+rt2zs*rt2zs)

               if((1.0-abs(unmag)).ge.tol) then
                 write(*,'(3(1X,A),3I5,2(2X,G14.7))') 'normal vector '
     $            ,cbtype1, cbtype2, ieg, iface, ia, unmag, amul
                 if    (if3d) then
                   if(amul.eq.2.) then
                    if    ((1.0-abs(t1mag)).ge.tol) then
                      c1mask(ix,iy,iz,iel) = 0.
                      c2mask(ix,iy,iz,iel) = 0.
                    elseif((1.0-abs(t2mag)).ge.tol) then
                      c1mask(ix,iy,iz,iel) = 0.
                      c3mask(ix,iy,iz,iel) = 0.
                    endif
                   endif
                   if(amul.eq.3.) then
                      c1mask(ix,iy,iz,iel) = 0.
                      c2mask(ix,iy,iz,iel) = 0.
                      c3mask(ix,iy,iz,iel) = 0.
                   endif
                 else
                   if(amul.eq.2.) then
                      c1mask(ix,iy,iz,iel) = 0.
                      c2mask(ix,iy,iz,iel) = 0.
                   endif
                 endif
               endif

 20         continue
         endif
 200  continue
 2000 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine distf2(d,ifld,b1,b2,dmin,emin,xn,yn,zn)

c     Generate a distance function to boundary with bc "b1" or "b2".
c     This approach does not yet work with periodic boundary conditions.

c     INPUT:  ifld - field type for which distance function is to be found.
c             ifld = 1 for velocity
c             ifld = 2 for temperature, etc.

c     OUTPUT: d = distance to nearest boundary with boundary condition "b"

c     Work arrays:  dmin,emin,xn,yn,zn

      include 'SIZE'
      include 'GEOM'       ! Coordinates
      include 'INPUT'      ! cbc()
      include 'TSTEP'      ! nelfld
      include 'PARALLEL'   ! gather-scatter handle for field "ifld"

      real d(lx1,ly1,lz1,lelt)
      character*3 b1, b2

      real dmin(lx1,ly1,lz1,lelt),emin(lx1,ly1,lz1,lelt)
      real xn(lx1,ly1,lz1,lelt),yn(lx1,ly1,lz1,lelt)
      real zn(lx1,ly1,lz1,lelt)


      integer e,eg,f

      nel = nelfld(ifld)
      n = lx1*ly1*lz1*nel

      call domain_size(xmin,xmax,ymin,ymax,zmin,zmax)

      xmn = min(xmin,ymin)
      xmx = max(xmax,ymax)
      if (if3d) xmn = min(xmn ,zmin)
      if (if3d) xmx = max(xmx ,zmax)

      big = 10*(xmx-xmn)
      call cfill (d,big,n)

      call opcopy(xn,yn,zn,xm1,ym1,zm1)

      nface = 2*ldim
      do e=1,nel     ! Set d=0 on walls
      do f=1,nface
         if (cbc(f,e,1).eq.b1 .or. cbc(f,e,1).eq.b2)
     $              call facev(d,e,f,0.,lx1,ly1,lz1)
      enddo
      enddo

      nxyz = lx1*ly1*lz1

      do ipass=1,10000
         dmax    = 0
         nchange = 0
         do e=1,nel
            do k=1,lz1
            do j=1,ly1
            do i=1,lx1
              i0=max(  1,i-1)
              j0=max(  1,j-1)
              k0=max(  1,k-1)
              i1=min(lx1,i+1)
              j1=min(ly1,j+1)
              k1=min(lz1,k+1)
              do kk=k0,k1
              do jj=j0,j1
              do ii=i0,i1

               dself  = d(i,j,k,e)
               dneigh = d(ii,jj,kk,e)
               if (dneigh.lt.dself) then  ! check neighbor's nearest point
                  d2 = (xm1(i,j,k,e)-xn(ii,jj,kk,e))**2
     $               + (ym1(i,j,k,e)-yn(ii,jj,kk,e))**2
                  if (if3d) d2 = d2 + (zm1(i,j,k,e)-zn(ii,jj,kk,e))**2
                  if (d2.gt.0) d2 = sqrt(d2)
                  if (d2.lt.dself) then
                    nchange = nchange+1
                    d (i,j,k,e) = d2
                    xn(i,j,k,e) = xn(ii,jj,kk,e)
                    yn(i,j,k,e) = yn(ii,jj,kk,e)
                    zn(i,j,k,e) = zn(ii,jj,kk,e)
                    dmax = max(dmax,d(i,j,k,e))
                  endif
               endif
              enddo
              enddo
              enddo

            enddo
            enddo
            enddo

            re = lglel(e)
            call cfill(emin(1,1,1,e),re,nxyz)
            call copy (dmin(1,1,1,e),d(1,1,1,e),nxyz)

         enddo
         nchange = iglsum(nchange,1)

         call fgslib_gs_op(gsh_fld(ifld),dmin,1,3,0) ! min over all elements


         nchange2=0
         do e=1,nel
         do i=1,nxyz
          if (dmin(i,1,1,e).ne.d(i,1,1,e)) then
             nchange2 = nchange2+1
             emin(i,1,1,e) = 0  ! Flag
          endif
         enddo
         enddo
         call copy(d,dmin,n)                !   Ensure updated distance
         nchange2 = iglsum(nchange2,1)
         nchange  = nchange + nchange2
         call fgslib_gs_op(gsh_fld(ifld),emin,1,4,0) ! max over all elements

         do e=1,nel    ! Propagate nearest wall points
         do i=1,nxyz
          eg = emin(i,1,1,e)
          if (eg.ne.lglel(e)) then
             xn(i,1,1,e) = 0
             yn(i,1,1,e) = 0
             zn(i,1,1,e) = 0
          endif
         enddo
         enddo
         call fgslib_gs_op(gsh_fld(ifld),xn,1,1,0) !   Sum over all elements to
         call fgslib_gs_op(gsh_fld(ifld),yn,1,1,0) !   convey nearest point
         call fgslib_gs_op(gsh_fld(ifld),zn,1,1,0) !   to shared neighbor.

         dmax = glmax(dmax,1)
         if (nio.eq.0) write(6,1) ipass,nchange,dmax
    1    format(i9,i12,1pe12.4,' max wall distance 2')
         if (nchange.eq.0) goto 1000
      enddo
 1000 continue

c     wgt = 0.3
c     call filter_d2(d,lx1,lz1,wgt,.true.)

      return
      end
c-----------------------------------------------------------------------
