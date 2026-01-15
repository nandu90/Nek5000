c-----------------------------------------------------------------------
      subroutine hsolve_cls(name,u,r,h1,h2,vmk,vml,
     $                 imsh,tol,maxit,isd
     $                 ,approx,napprox,bi)
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'LVLSET'
c
      CHARACTER*4    NAME
      REAL           U    (LX1,LY1,LZ1,1)
      REAL           R    (LX1,LY1,LZ1,1)
      REAL           H1   (LX1,LY1,LZ1,1)
      REAL           H2   (LX1,LY1,LZ1,1)
      REAL           vmk  (LX1,LY1,LZ1,1)
      REAL           vml  (LX1,LY1,LZ1,1)
      REAL           bi   (LX1,LY1,LZ1,1)
      REAL           approx (1)
      integer        napprox(1)
      common /ctmp2/ w1   (lx1,ly1,lz1,lelt)
      common /ctmp3/ w2   (2+2*mxprev)

      logical ifstdh
      character*4  cname
      character*6  name6

      logical ifwt,ifvec

      call chcopy(cname,name,4)
      call capit (cname,4)

      ifstdh = .true. ! default is no projection
      if (ifprojfld(ifield)) then
         ifstdh = .false.
      endif

      p945 = param(94)
      if (cname.eq.'PRES') then
         ifstdh = .false.
         p945 = param(95)
      endif

      if (ifield.gt.ldimt_proj+1) ifstdh = .true.
      if (param(93).eq.0)         ifstdh = .true.
      if (p945.eq.0)              ifstdh = .true.
      if (istep.lt.p945)          ifstdh = .true.

      if(ifls_debug.eq.1 .and. nio.eq.0)write(*,*)
     $   "In hmholtz_cls",ifstdh

      if(ifls_debug.eq.1)call lsmonitor(r,'rhslv')

      if (ifstdh) then
         call hmholtz_cls(name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd)
      else

         n = lx1*ly1*lz1*nelfld(ifield)

         call col2   (r,vmk,n)
         call dssum  (r,lx1,ly1,lz1)

         call blank (name6,6)
         call chcopy(name6,name,4)
         ifwt  = .true.
         ifvec = .false.

         call project1_cls
     $       (r,n,approx,napprox,h1,h2,vmk,vml,ifwt,ifvec,name6)

         call hmhzpf_cls (name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd,bi)

         call project2_cls
     $       (u,n,approx,napprox,h1,h2,vmk,vml,ifwt,ifvec,name6)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine hmholtz_cls(name,u,rhs,h1,h2,mask,mult,
     $  imsh,tli,maxit,isd)
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'CTIMER'
      include 'LVLSET'

      CHARACTER      NAME*4
      REAL           U    (LX1,LY1,LZ1,1)
      REAL           RHS  (LX1,LY1,LZ1,1)
      REAL           H1   (LX1,LY1,LZ1,1)
      REAL           H2   (LX1,LY1,LZ1,1)
      REAL           MASK (LX1,LY1,LZ1,1)
      REAL           MULT (LX1,LY1,LZ1,1)

      logical iffdm
      character*3 nam3

      tol = abs(tli)

      iffdm = .false.
      if (ifsplit) iffdm = .true.
      if (icalld.eq.0.and.iffdm) call set_fdm_prec_h1A
      icalld = icalld+1

#ifdef TIMER
      if (name.ne.'PRES') then
        nhmhz = nhmhz + 1
        etime1 = dnekclock()
      endif
#endif


      ntot = lx1*ly1*lz1*nelfld(ifield)
      if (imsh.eq.1) ntot = lx1*ly1*lz1*nelv
      if (imsh.eq.2) ntot = lx1*ly1*lz1*nelt

      if(ifls_debug.eq.1 .and. nio.eq.0)then
        write(*,*)"In hmholtz_cls ",name,imsh
      endif

c
      call chcopy(nam3,name,3)
c
                          kfldfdm = -1
c     if (nam3.eq.'TEM' ) kfldfdm =  0
c     if (name.eq.'TEM1') kfldfdm =  0  ! hardcode for temp only, for mpaul
c     if (name.eq.'VELX') kfldfdm =  1
c     if (name.eq.'VELY') kfldfdm =  2
c     if (name.eq.'VELZ') kfldfdm =  3
      if (name.eq.'PRES') kfldfdm =  ldim+1
c     if (.not.iffdm) kfldfdm=-1
C
      call dssum   (rhs,lx1,ly1,lz1)

      if(ifls_debug.eq.1)call lsmonitor(rhs,'rhmlz')
      if(ifls_debug.eq.1)call lsmonitor(mask,'mhmlz')
      if(ifls_debug.eq.1)call lsmonitor(mult,'multz')

      call col2    (rhs,mask,ntot)
c      if (nio.eq.0.and.istep.le.10) 
c     $    write(6,*) param(22),' p22 ',istep,imsh
      if (param(22).eq.0.or.istep.le.10)
     $    call chktcg1 (tol,rhs,h1,h2,mask,mult,imsh,isd)

      if (tli.lt.0) tol=tli ! caller-specified relative tolerance

      if (imsh.eq.1) call cggo_cls
     $   (u,rhs,h1,h2,mask,mult,imsh,tol,maxit,isd,binvm1,name)
      if (imsh.eq.2) call cggo_cls
     $   (u,rhs,h1,h2,mask,mult,imsh,tol,maxit,isd,bintm1,name)

#ifdef TIMER
      if (name.ne.'PRES') thmhz=thmhz+(dnekclock()-etime1)
#endif

      return
      END
c-----------------------------------------------------------------------
      subroutine hmhzpf_cls(name,u,r,h1,h2,mask,mult,imesh,
     $   tli,maxit,isd,bi)
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'FDMH1'
      include 'CTIMER'
c
      CHARACTER*4    NAME
      REAL           U    (LX1,LY1,LZ1,1)
      REAL           R    (LX1,LY1,LZ1,1)
      REAL           H1   (LX1,LY1,LZ1,1)
      REAL           H2   (LX1,LY1,LZ1,1)
      REAL           MASK (LX1,LY1,LZ1,1)
      REAL           MULT (LX1,LY1,LZ1,1)
      REAL           bi   (LX1,LY1,LZ1,1)
      COMMON /CTMP0/ W1   (LX1,LY1,LZ1,LELT)
     $ ,             W2   (LX1,LY1,LZ1,LELT)
c
      etime1=dnekclock()
c
      IF (IMESH.EQ.1) NTOT = lx1*ly1*lz1*NELV
      IF (IMESH.EQ.2) NTOT = lx1*ly1*lz1*NELT
c
      tol = tli
      if (param(22).ne.0) tol = abs(param(22))
      CALL CHKTCG1 (TOL,R,H1,H2,MASK,MULT,IMESH,ISD)
c
c
c     Set flags for overlapping Schwarz preconditioner (pff 11/12/98)
c
                          kfldfdm = -1
c     if (name.eq.'TEMP') kfldfdm =  0
c     if (name.eq.'VELX') kfldfdm =  1
c     if (name.eq.'VELY') kfldfdm =  2
c     if (name.eq.'VELZ') kfldfdm =  3
      if (name.eq.'PRES') kfldfdm =  ldim+1

      if (ifdg) then
         call cggo_dg (u,r,h1,h2,bi,mask,name,tol,maxit)
      else
         call cggo_cls
     $      (u,r,h1,h2,mask,mult,imesh,tol,maxit,isd,bi,name)
      endif
      thmhz=thmhz+(dnekclock()-etime1)
c
c
      return
      end
c-----------------------------------------------------------------------
      subroutine cggo_cls(x,f,h1,h2,mask,mult,imsh,
     $  tin,maxit,isd,binv,name)
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'LVLSET'
 
      COMMON  /CPRINT/ IFPRINT, IFHZPC
      LOGICAL          IFPRINT, IFHZPC
 
      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv

      logical ifmcor,ifprint_hmh
 
      real x(1),f(1),h1(1),h2(1),mask(1),mult(1),binv(1)
      parameter        (lg=lx1*ly1*lz1*lelt)
      COMMON /SCRCG/ d (lg) , scalar(2)
      common /SCRMG/ r (lg) , w (lg) , p (lg) , z (lg)
c
      parameter (maxcg=900)
      common /tdarray/ diagt(maxcg),upper(maxcg)
      common /iterhm/ niterhm
      character*4 name

      if (ifsplit.and.name.eq.'PRES') then
         if (param(42).eq.0) then
           n = lx1*ly1*lz1*nelv
           call copy      (x,f,n)
           iter = maxit
           call hmh_gmres (x,h1,h2,mult,iter)
           niterhm = iter
           return
         elseif(param(42).eq.2) then 
           n = lx1*ly1*lz1*nelv
           call copy       (x,f,n)
           iter = maxit
           call hmh_flex_cg(x,h1,h2,mult,iter)
           niterhm = iter
           return
         endif
      endif

      if(ifls_debug.eq.1 .and. nio.eq.0)then
        write(*,*)"in cggo_cls"
      endif

      if(ifls_debug.eq.1)call lsmonitor(f,'fcggo')

      call rzero(diagt,maxcg)
      call rzero(upper,maxcg)
      rho = 0.00

      NXYZ   = lx1*ly1*lz1
      NEL    = NELV
      VOL    = VOLVM1
      IF (IMSH.EQ.2) NEL=NELT
      IF (IMSH.EQ.2) VOL=VOLTM1
      n      = NEL*NXYZ

      tol=abs(tin)

      if (restol(ifield).ne.0) tol=restol(ifield)
      if (name.eq.'PRES'.and.param(21).ne.0) tol=abs(param(21))

      if (tin.lt.0) tol=abs(tin)
      niter = min(maxit,maxcg)

      if (.not.ifsolv) then
         call setfast(h1,h2,imsh)
         ifsolv = .true.
      endif

      if (kfldfdm.lt.0) then
         call setprec(D,h1,h2,imsh,isd)
      elseif(param(100).ne.2) then
         call set_fdm_prec_h1b(d,h1,h2,nel)
      endif

      call copy (r,f,n)
      call rzero(x,n)
      call rzero(p,n)

      fmax = glamax(f,n)
      if (fmax.eq.0.0) return


      ifmcor = .false.
      h2max = glmax(h2  ,n)
      skmin = glmin(mask,n)
      if (skmin.gt.0.and.h2max.eq.0) ifmcor = .true.

      if (name.eq.'PRES') then
        continue
      elseif (ifmcor) then
         smean = -1./glsum(bm1,n) ! Modified 5/4/12 pff
         rmean = smean*glsc2(r,mult,n)
         call copy(x,bm1,n)
         call dssum(x,lx1,ly1,lz1)
         call add2s2(r,x,rmean,n)
         call rzero(x,n)
      endif
C
      krylov = 0
      rtz1=1.0
      niterhm = 0

      do iter=1,niter
        if (kfldfdm.lt.0) then  ! Jacobi Preconditioner
          call col3(z,r,d,n)
        else                                       ! Schwarz Preconditioner
          if (name.eq.'PRES'.and.param(100).eq.2) then
            call h1_overlap_2(z,r,mask)
            call crs_solve_h1 (w,r)  ! Currently, crs grd only for P
            call add2         (z,w,n)
          else   
            call fdm_h1(z,r,d,mask,mult,nel,ktype(1,1,kfldfdm),w)
            if (name.eq.'PRES') then 
              call crs_solve_h1 (w,r)  ! Currently, crs grd only for P
              call add2         (z,w,n)
            endif
          endif
        endif
        if (name.eq.'PRES') then
          call ortho (z)
        elseif (ifmcor) then
          rmean = smean*glsc2(z,bm1,n)
          call cadd(z,rmean,n)
        endif

        rtz2=rtz1
        scalar(1)=vlsc3 (z,r,mult,n)
        if(param(18).eq.1) then
          scalar(2)=vlsc3(r,r,mult,n)
        else 
          scalar(2)=vlsc32(r,mult,binv,n)
        endif
        call gop(scalar,w,'+  ',2)
        rtz1=scalar(1)
        rbn2=sqrt(scalar(2)/vol)
        if (iter.eq.1) rbn0 = rbn2
        if (param(22).lt.0) tol=abs(param(22))*rbn0
        if (tin.lt.0)       tol=abs(tin)*rbn0

         ifprint_hmh = .false.
         if (nio.eq.0.and.ifprint.and.param(74).ne.0) ifprint_hmh=.true.
         ! if (nio.eq.0.and.istep.eq.1)                 ifprint_hmh=.true.

         if (ifprint_hmh)
     &      write(6,3002) istep,'  Hmholtz ' // name,
     &                    iter,rbn2,h1(1),tol,h2(1),ifmcor



         IF (rbn2.LE.TOL.and.(iter.gt.1 .or. istep.le.5)) THEN
            NITER = ITER-1
            if (nio.eq.0)
     &         write(6,3000) istep,'  Hmholtz ' // name,
     &                       niter,rbn2,rbn0,tol
            goto 9999
         ENDIF
c
         beta = rtz1/rtz2
         if (iter.eq.1) beta=0.0
         call add2s1 (p,z,beta,n)
         call axhelm_cls (w,p,h1,h2,imsh,isd)
         call dssum  (w,lx1,ly1,lz1)
         call col2   (w,mask,n)
c
         rho0 = rho
         rho  = glsc3(w,p,mult,n)
         alpha=rtz1/rho
         alphm=-alpha
         call add2s2(x,p ,alpha,n)
         call add2s2(r,w ,alphm,n)
c
c        Generate tridiagonal matrix for Lanczos scheme
         if (iter.eq.1) then
            krylov = krylov+1
            diagt(iter) = rho/rtz1
         elseif (iter.le.maxcg) then
            krylov = krylov+1
            diagt(iter)    = (beta**2 * rho0 + rho ) / rtz1
            upper(iter-1)  = -beta * rho0 / sqrt(rtz2 * rtz1)
         endif
      enddo
      niter = iter-1
c
      if (nio.eq.0) write (6,3001) istep, '  Error Hmholtz ' // name,
     &                             niter,rbn2,rbn0,tol


3000  format(i11,a,1x,I7,1p4E13.4)
3001  format(i11,a,1x,I7,1p4E13.4)
3002  format(i11,a,1x,I7,1p4E13.4,l4)
9999  continue
      niterhm = niter
      ifsolv = .false.
      return
      end
c-----------------------------------------------------------------------
      subroutine project1_cls(b,n,rvar,ivar
     $                       ,h1,h2,msk,w,ifwt,ifvec,name6)

c     1. Compute the projection of x onto X

c     2. Re-orthogonalize the X basis set and corresponding B=A*X
c        vectors if A has changed.

c     Output:  b = b - projection of b onto B

c     Input:   n     = length of field (or multifields, when ifvec=true)
c              rvar  = real array of field values, including old h1,h2, etc.
c              ivar  = integer array of pointers, etc.
c              h1    = current h1, for Axhelm(.,.,h1,h2,...)
c              h2    = current h2
c              msk   = mask for Dirichlet BCs
c              w     = weight for inner products (typ. w=vmult, tmult, etc.)
c              ifwt  = use weighted inner products when ifwt=.true.
c              ifvec = are x and b vectors, or scalar fields?
c              name6 = discriminator for action of A*x

c     The idea here is to have one pair of projection routines for
c     constructing the new rhs (project1) and reconstructing the new
c     solution (x = xbar + dx) plus updating the approximation space.
c     The latter functions are done in project2.
c
c     The approximation space X and corresponding right-hand sides,
c     B := A*X are stored in rvar, as well as h1old and h2old and a
c     couple of other auxiliary arrays.

c     In this new code, we retain both X and B=A*X and we re-orthogonalize
c     at each timestep (with no extra matrix-vector products, but O(nm)
c     work.   The idea is to retain fresh vectors by injecting the most
c     recent solution and pushing the oldest off the stack, hopefully
c     keeping the number of vectors, m, small.


      include 'SIZE'   ! For nid/nio
      include 'TSTEP'  ! For istep
      include 'CTIMER'

      real b(n),rvar(n,1),h1(n),h2(n),w(n),msk(n)
      integer ivar(1)
      character*6 name6
      logical ifwt,ifvec

      etime0 = dnekclock()

      nn = n
      if (ifvec) nn = n*ldim

      call proj_get_ivar
     $   (m,mmx,ixb,ibb,ix,ib,ih1,ih2,ivar,n,ifvec,name6)

      if (m.le.0) return

      ireset=iproj_chk_ls(rvar(ih1,1),rvar(ih2,1),h1,h2,n) ! Updated matrix?

      bb4 = glsc3(b,w,b,n)
      bb4 = sqrt(bb4)


c     Re-orthogonalize basis set w.r.t. new vectors if space has changed.

      if (ireset.eq.1) then

         do j=0,m-1         ! First, set B := A*X
            jb = ib+j*nn    !Iterate through xs and bs
            jx = ix+j*nn
            call proj_matvec_cls(rvar(jb,1),rvar(jx,1),n
     $                          ,h1,h2,msk,name6)
         enddo

         if (nio.eq.0 .and. loglevel.gt.2)
     $      write(6,'(13x,A)') 'Reorthogonalize Basis'

         call proj_ortho_full
     $      (rvar(ix,1),rvar(ib,1),n,m,w,ifwt,ifvec,name6)

         ivar(2) = m ! Update number of saved vectors

      endif

c     ixb is pointer to xbar,  ibb is pointer to bbar := A*xbar

      call project1_a(rvar(ixb,1),rvar(ibb,1),b,rvar(ix,1),rvar(ib,1)
     $               ,n,m,w,ifwt,ifvec)

      baf = glsc3(b,w,b,n)
      baf = sqrt(baf)
      ratio = bb4/baf

      tproj = tproj + dnekclock() - etime0

      if (nio.eq.0) write(6,1) istep,'  Project ' // name6,
     &                         baf,bb4,ratio,m,mmx
    1 format(i11,a,6x,1p3e13.4,i4,i4)

      return
      end
c-----------------------------------------------------------------------
      subroutine project2_cls(x,n,rvar,ivar
     $                       ,h1,h2,msk,w,ifwt,ifvec,name6)

      include 'CTIMER'

      real x(n),b(n),rvar(n,1),h1(n),h2(n),w(n),msk(n)
      integer ivar(1)
      character*6 name6
      logical ifwt,ifvec

      etime0 = dnekclock()

      call proj_get_ivar(m,mmx,ixb,ibb,ix,ib,ih1,ih2,
     $                             ivar,n,ifvec,name6)

c     ix  is pointer to X,     ib  is pointer to B
c     ixb is pointer to xbar,  ibb is pointer to bbar := A*xbar

      call project2_a_cls(x,rvar(ixb,1),rvar(ix,1),rvar(ib,1)
     $              ,n,m,mmx,h1,h2,msk,w,ifwt,ifvec,name6)

      ivar(2) = m ! Update number of saved vectors

      tproj = tproj + dnekclock() - etime0

      return
      end
c-----------------------------------------------------------------------
      subroutine project2_a_cls(x,xbar,xx,bb,n,m,mmx
     $                         ,h1,h2,msk,w,ifwt,ifvec,name6)

      include 'SIZE'

      real x(n),xbar(n),xx(n,1),bb(n,1),h1(n),h2(n),w(n),msk(n)
      character*6 name6
      logical ifwt,ifvec

      nn = n
      if (ifvec) nn=ldim*n

      if (m.gt.0) call add2(x,xbar,n)      ! Restore desired solution

       !Uncomment this if using full reorthogonalization
c      if (m.eq.mmx) then ! Push old vector off the stack
c         do k=2,mmx
c            call copy (xx(1,k-1),xx(1,k),nn)
c            call copy (bb(1,k-1),bb(1,k),nn)
c         enddo
c      endif

      m = min(m+1,mmx)
      !print *, "m", m
      call copy        (xx(1,m),x,nn)   ! Update (X,B)
      call proj_matvec_cls(bb(1,m),xx(1,m),n,h1,h2,msk,name6)
      call proj_ortho  (xx,bb,n,m,w,ifwt,ifvec,name6) !Update orthogonalization
      !Uncomment the if block above if using full reorthogonalization
c      call proj_ortho_full  (xx,bb,n,m,w,ifwt,ifvec,name6) !Fully reorthogonalize

      return
      end
c-----------------------------------------------------------------------
      function iproj_chk_ls(h1old,h2old,h1,h2,n)
      include 'SIZE'
      include 'TOTAL'
c     Matrix has changed if h1/h2 differ from old values

      real h1(n),h2(n),h1old(n),h2old(n)
      integer istep_save
      save istep_save
      data istep_save /0/

      iproj_chk_ls = 0

      if (ifmvbd) then
         iproj_chk_ls = 1
         return
      endif

      if (istep.lt.istep_save) then
         iproj_chk_ls = 1
         return
      endif
      istep_save = istep

      dh1 = 0.
      dh2 = 0.
      do i=1,n
         dh1 = max(dh1,abs(h1(i)-h1old(i)))
         dh2 = max(dh2,abs(h2(i)-h2old(i)))
      enddo
      dh = max(dh1,dh2)
      dh = glmax(dh,1)  ! Max across all processors

      if (dh.gt.0) then

         call copy(h1old,h1,n)   ! Save old h1 / h2 values
         call copy(h2old,h2,n)

         iproj_chk_ls = 1      ! Force re-orthogonalization of basis

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine proj_matvec_cls(b,x,n,h1,h2,msk,name6)
      include 'SIZE'
      include 'TOTAL'
      real b(n),x(n),h1(n),h2(n),msk(n)
      character*6 name6

c     This is the default matvec for nekcem.

c     The code can later be updated to support different matvec
c     implementations, which would be discriminated by the character
c     string "name6"

      isd  = 1    ! This probably won't work for axisymmetric
      imsh = 1
      if (iftmsh(ifield)) imsh=2

      call axhelm_cls  (b,x,h1,h2,imsh,isd)       ! b = A x
      call dssum   (b,lx1,ly1,lz1)
      call col2    (b,msk,n)

      return
      end
c-----------------------------------------------------------------------
