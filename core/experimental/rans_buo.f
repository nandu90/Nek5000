c--------------------------------------------------------------      
      subroutine rans_buo_init(id, g)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
      include 'RANS_BUO'

      integer id
      real g(3)

      ifsgdh = .false.
      ifggdh = .false.
      if(id.eq.1) ifsgdh = .true.
      if(id.eq.2) ifggdh = .true.

      if(ifaxis)then
        if(nid.eq.0)write(6,*)
     &   "ERROR:Axisymmetric BC not yet supported with Buoyancy"
        call exitt
      endif

      ! if(.not.ifsgdh .or. ifggdh)then
      !   if(nid.eq.0)write(6,*)
      ! &   "ERROR: Only SGDH supported for now"
      !   call exitt
      ! endif

      if(.not.ifrans_ktau_stndrd .and.
     &   .not.ifrans_komg_stndrd)then
        if(nid.eq.0)then
          write(6,*)
     &    "ERROR: Buoyancy is only supported with:"
          write(6,*)"std. k-tau     : m_id = 4"
          write(6,*)"reg. k-omega   : m_id = 0"
        endif
        call exitt
      endif

      buo_gvec(1) = g(1)
      buo_gvec(2) = g(2)
      buo_gvec(3) = g(3)
      
      return
      end
c--------------------------------------------------------------      
      real function buo_diff(ix,iy,iz,iel)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
      include 'RANS_BUO'

      integer ix,iy,iz,iel
      real Pr_t, rans_mut


      Pr_t = coeffs(1)
      mu_t = rans_mut(ix,iy,iz,iel)

      if(.not.ifggdh)then
        buo_diff = cpfld(ifield,1) + mu_t/Pr_t
      else
        buo_diff = 1.0
      endif

      return
      end
c--------------------------------------------------------------      
      real function buo_ksrc(ix,iy,iz,iel)
      implicit none
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'
      include 'RANS_BUO'
                        
      integer ix,iy,iz,iel
      
      logical ifevalbuo
      data ifevalbuo /.true./
      common/ifbuoSrc/ifevalbuo

      if(ix*iy*iz*iel .eq. 1 .and. ifevalbuo)then
        if (nid.eq.0 .and. loglevel.gt.2) 
     $   write(*,*)'updating buoyancy source terms'
     
        if(ifrans_ktau_stndrd) call rans_ktau_stndrd_buo
        if(ifrans_komg_stndrd) call rans_komega_stndrd_buo
        ifevalbuo = .false.
      endif

      if(ifld_k .gt. ifld_omega) ifevalbuo = .true.

      buo_ksrc = ksrc_buo(ix,iy,iz,iel)
    
      return
      end
c--------------------------------------------------------------      
      real function buo_omgsrc(ix,iy,iz,iel)
      implicit none
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'
      include 'RANS_BUO'

      integer ix,iy,iz,iel
      
      logical ifevalbuo
      data ifevalbuo /.true./
      common/ifbuoSrc/ifevalbuo

      if(ix*iy*iz*iel .eq. 1 .and. ifevalbuo)then
        if (nid.eq.0 .and. loglevel.gt.2) 
     $   write(*,*)'updating buoyancy source terms'
     
        if(ifrans_ktau_stndrd) call rans_ktau_stndrd_buo
        if(ifrans_komg_stndrd) call rans_komega_stndrd_buo
        ifevalbuo = .false.
      endif

      if(ifld_omega .gt. ifld_k) ifevalbuo = .true.

      buo_omgsrc = omgsrc_buo(ix,iy,iz,iel)
      return
      end
c--------------------------------------------------------------      
      real function buo_kdiag(ix,iy,iz,iel)
      implicit none
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'
      include 'RANS_BUO'
                        
      integer ix,iy,iz,iel
      
      logical ifevalbuo
      data ifevalbuo /.true./
      common/ifbuoSrc/ifevalbuo

      if(ix*iy*iz*iel .eq. 1 .and. ifevalbuo)then
        if (nid.eq.0 .and. loglevel.gt.2) 
     $   write(*,*)'updating buoyancy source terms'
     
        if(ifrans_ktau_stndrd) call rans_ktau_stndrd_buo
        if(ifrans_komg_stndrd) call rans_komega_stndrd_buo
        ifevalbuo = .false.
      endif

      if(ifld_k .gt. ifld_omega) ifevalbuo = .true.

      buo_kdiag = kdiag_buo(ix,iy,iz,iel)
    
      return
      end
c--------------------------------------------------------------      
      real function buo_omgdiag(ix,iy,iz,iel)
      implicit none
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'
      include 'RANS_BUO'
                        
      integer ix,iy,iz,iel
      
      logical ifevalbuo
      data ifevalbuo /.true./
      common/ifbuoSrc/ifevalbuo

      if(ix*iy*iz*iel .eq. 1 .and. ifevalbuo)then
        if (nid.eq.0 .and. loglevel.gt.2) 
     $   write(*,*)'updating buoyancy source terms'
     
        if(ifrans_ktau_stndrd) call rans_ktau_stndrd_buo
        if(ifrans_komg_stndrd) call rans_komega_stndrd_buo
        ifevalbuo = .false.
      endif

      if(ifld_omega .gt. ifld_k) ifevalbuo = .true.

      buo_omgdiag = omgdiag_buo(ix,iy,iz,iel)
    
      return
      end
c--------------------------------------------------------------      
      subroutine rans_ktau_stndrd_buo
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
      include 'RANS_BUO'

      integer lxyz
      parameter(lxyz=lx1*ly1*lz1)

      real xflux(lxyz),yflux(lxyz),zflux(lxyz)
      integer e,i

      real Pr_t,alpinf_str,alp_inf
      real rho,tau,mu_t0,dotsrc
      real alp_str,alpha,gamm,G_w,G_k

      integer icalld
      save icalld
      data icalld /0/
      
      Pr_t         = coeffs(1)
      alpinf_str   = coeffs(4)
      alp_inf      = coeffs(9)

      if(.not.ifsgdh .and. .not.ifggdh)then
        ! no Buoyancy production term in k-tau eqns
        if(icalld.eq.0)then
          call rzero(ksrc_buo,lxyz*nelt)
          call rzero(kdiag_buo,lxyz*nelt)
          call rzero(omgsrc_buo,lxyz*nelt)
          call rzero(omgdiag_buo,lxyz*nelt)
          icalld=1
        endif
        return
      endif

      if(ifield.eq.3)then
        ! Evaluate vel and temp gradient once per time step
        call eval_tgrad
        if(ifggdh) call eval_vgrad
      endif

      do e=1,nelv
         if(ifsgdh)call flux_sgdh_compute(xflux,yflux,zflux,e)
         if(ifggdh)call flux_ggdh_compute(xflux,yflux,zflux,e)

         do i=1,lxyz
            rho = vtrans(i,1,1,e,1)

            k = t(i,1,1,e,ifld_k-1)
            tau = t(i,1,1,e,ifld_omega-1)

            alp_str = alpinf_str
            mu_t = rho*alp_str*k*tau
            mu_t0 = rho*alp_str* tau

            dotsrc = buo_gvec(1)*xflux(i) + buo_gvec(2)*yflux(i)
            if(if3d) dotsrc = dotsrc + buo_gvec(3)*zflux(i)
            dotsrc = dotsrc/Pr_t

            if(ifrans_diag)then
              ksrc_buo(i,1,1,e) = 0.0
              kdiag_buo(i,1,1,e) = -mu_t0 * dotsrc
            else
              ksrc_buo(i,1,1,e) = mu_t0 *dotsrc * k
              kdiag_buo(i,1,1,e) = 0.0
            endif

            alpha = alp_inf/alp_str
            gamm = alpha*alp_str

            G_w = rho * gamm*tau*dotsrc

            ! The sign of production in tau equation is opposite of what i expect
            ! With opposite sign the tau is unstable and keeps going -ve
            ! To do: Confirm this...
            if(ifrans_diag)then
              omgsrc_buo(i,1,1,e) = 0.0 
              omgdiag_buo(i,1,1,e) = -G_w
            else
              omgsrc_buo(i,1,1,e) = G_w*tau
              omgdiag_buo(i,1,1,e) = 0.0
            endif
         enddo
      enddo
      
      return
      end
c--------------------------------------------------------------      
      subroutine rans_komega_stndrd_buo
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
      include 'RANS_BUO'

      integer lxyz
      parameter(lxyz=lx1*ly1*lz1)

      real xflux(lxyz),yflux(lxyz),zflux(lxyz)
      integer e,i

      real Pr_t,alpinf_str,alp_inf
      real rho,omega,mu_t0,dotsrc
      real alp_str,alpha,gamm,G_w
      real tiny

      integer icalld
      save icalld
      data icalld /0/
      
      Pr_t         = coeffs(1)
      alpinf_str   = coeffs(4)
      alp_inf      = coeffs(9)
      tiny         = coeffs(16)

      if(.not.ifsgdh .and. .not.ifggdh)then
        ! no Buoyancy production term in k-tau eqns
        if(icalld.eq.0)then
          call rzero(ksrc_buo,lxyz*nelt)
          call rzero(kdiag_buo,lxyz*nelt)
          call rzero(omgsrc_buo,lxyz*nelt)
          call rzero(omgdiag_buo,lxyz*nelt)
          icalld=1
        endif
        return
      endif

      if(ifield.eq.3)then
        ! Evaluate vel and temp gradient once per time step
        call eval_tgrad
        if(ifggdh) call eval_vgrad
      endif

      do e=1,nelv
         if(ifsgdh)call flux_sgdh_compute(xflux,yflux,zflux,e)
         if(ifggdh)call flux_ggdh_compute(xflux,yflux,zflux,e)

         do i=1,lxyz
            rho = vtrans(i,1,1,e,1)

            k = t(i,1,1,e,ifld_k-1)
            omega = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e)

            alp_str = alpinf_str
            mu_t = rho*alp_str*k/(omega+tiny)

            dotsrc = buo_gvec(1)*xflux(i) + buo_gvec(2)*yflux(i)
            if(if3d) dotsrc = dotsrc + buo_gvec(3)*zflux(i)
            dotsrc = dotsrc/Pr_t

            if(ifrans_diag)then
              ksrc_buo(i,1,1,e) = mu_t * dotsrc
              kdiag_buo(i,1,1,e) = 0.0
            else
              ksrc_buo(i,1,1,e) = mu_t * dotsrc
              kdiag_buo(i,1,1,e) = 0.0
            endif

            alpha = alp_inf/alp_str
            gamm = alpha*alp_str

            G_w = rho*gamm*dotsrc

            if(ifrans_diag)then
              omgsrc_buo(i,1,1,e) = G_w 
              omgdiag_buo(i,1,1,e) = 0.0
            else
              omgsrc_buo(i,1,1,e) = G_w
              omgdiag_buo(i,1,1,e) = 0.0
            endif
         enddo
      enddo
      
      return
      end
c--------------------------------------------------------------      
      subroutine flux_sgdh_compute(xflux,yflux,zflux,ie)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_BUO'
      
      integer lxyz
      parameter(lxyz=lx1*ly1*lz1)

      real xflux(1),yflux(1),zflux(1)
      integer i, ie

      do i = 1, lxyz
        xflux(i) = tx_buo(i,1,1,ie)
        yflux(i) = ty_buo(i,1,1,ie)
        if(if3d)zflux(i) = tz_buo(i,1,1,ie)
      enddo

      return
      end
c--------------------------------------------------------------      
      subroutine eval_tgrad
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_BUO'

      ! if(nio.eq.0)write(*,*)"called tgrad",ifield
      call gradm1(tx_buo,ty_buo,tz_buo,t)
      call opcolv(tx_buo,ty_buo,tz_buo,bm1)
      call opdssum(tx_buo,ty_buo,tz_buo)
      call opcolv(tx_buo,ty_buo,tz_buo,binvm1)

      return
      end
c--------------------------------------------------------------      
      subroutine flux_ggdh_compute(xflux,yflux,zflux,ie)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
      include 'RANS_BUO'
      
      common /scruz/ u_x (lx1,ly1,lz1,lelv),u_y (lx1,ly1,lz1,lelv)
     $              ,u_z (lx1,ly1,lz1,lelv),v_x (lx1,ly1,lz1,lelv)
     $              ,v_y (lx1,ly1,lz1,lelv),v_z (lx1,ly1,lz1,lelv)
     $              ,w_x (lx1,ly1,lz1,lelv),w_y (lx1,ly1,lz1,lelv)
     $              ,w_z (lx1,ly1,lz1,lelv)
      real u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z

      integer lxyz
      parameter(lxyz=lx1*ly1*lz1)

      real xflux(1),yflux(1),zflux(1)
      integer i, ie

      real s11,s12,s13,s22,s23,s33
      real mu_t0,rho,rans_mut
      real gtx,gty,gtz

      do i = 1, lxyz
        mu_t = rans_mut(i,1,1,ie)
        rho = vtrans(i,1,1,ie,1)
        k = t(i,1,1,ie,ifld_k-1)
        mu_t0 = 0.0
        if(k.ne.0.)mu_t0 = mu_t/k
        
        gtx = tx_buo(i,1,1,ie)
        gty = ty_buo(i,1,1,ie)
        gtz = tz_buo(i,1,1,ie)
        if(if3d)then
          s11 = (2./3.)*rho - mu_t0*2.0*u_x(i,1,1,ie)
          s22 = (2./3.)*rho - mu_t0*2.0*v_y(i,1,1,ie)
          s33 = (2./3.)*rho - mu_t0*2.0*w_z(i,1,1,ie)
          s12 = -mu_t0*(u_y(i,1,1,ie)+v_x(i,1,1,ie))
          s23 = -mu_t0*(w_y(i,1,1,ie)+v_z(i,1,1,ie))
          s13 = -mu_t0*(u_z(i,1,1,ie)+w_x(i,1,1,ie))
          xflux(i) = s11*gtx+s12*gty+s13*gtz
          yflux(i) = s12*gtx+s22*gty+s23*gtz
          zflux(i) = s13*gtx+s23*gty+s33*gtz
        elseif(ifaxis)then  !To Do
          if(nio.eq.0)then
            write(*,*)"Axisymmetric not yet supported with Buoyancy"
          endif
          call exitt
        else
          s11 = (2./3.)*rho - mu_t0*2.0*u_x(i,1,1,ie)
          s22 = (2./3.)*rho - mu_t0*2.0*v_y(i,1,1,ie)
          s12 = -mu_t0*(u_y(i,1,1,ie)+v_x(i,1,1,ie))
          xflux(i) = s11*gtx+s12*gty
          yflux(i) = s12*gtx+s22*gty
        endif
      enddo

      return
      end
c--------------------------------------------------------------      
      subroutine store_geom
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_BUO'

      integer nxyz,ntot

      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelt

      !Memory of all geom factors is contiguous
      call copy(gms,g1m1,ntot*6)

      return
      end
c--------------------------------------------------------------      
      subroutine compute_ggdh_geom
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
      include 'RANS_BUO'

      common /storewj/ wj(lx1*ly1*lz1,lelt)
      real wj
  
      integer icalld
      save icalld
      data icalld /0/

      integer lxyz
      parameter(lxyz=lx1*ly1*lz1)

      real grs(lxyz,lelt,6)
      equivalence(grs,g1m1)

      real rxg(lxyz,lelt,3,3)
      equivalence(rxg,rxm1)

      real s_ij(lxyz,6)
      integer i,ie,j,ii,jj

      integer ir1_2d(3),ir2_2d(3)
      save ir1_2d,ir2_2d
      data ir1_2d /1,2,1/
      data ir2_2d /1,2,2/

      integer q_2d(2,2)
      save q_2d
      data q_2d /1,3,3,2/

      integer ir1,ir2,i2

      integer ntot

      if(icalld.eq.0)then
        !Store weighted jacobian
        do ie=1,nelt
          call invcol3(wj(1,ie),w3m1,jacm1(1,1,1,ie),lxyz)
        enddo
        icalld = 1
      endif

      ntot=lxyz*nelt

      call rzero(grs,ntot*6)

      do ie=1,nelt
        call compute_aniso_tensor(s_ij,ie)
        do ii=1,2
          ir1 = ir1_2d(ii)
          ir2 = ir2_2d(ii)
          i2 = ii
          if(ii.eq.3)i2=i2+1
          do j=1,2
            do i=1,2
              jj = q_2d(i,j)
              call addcol5(grs(1,ie,i2),rxg(1,ie,ir1,i),
     &          rxg(1,ie,ir2,j),wj(1,ie),s_ij(1,jj),lxyz)
            enddo
          enddo
        enddo
      enddo

      return
      end
c--------------------------------------------------------------      
      subroutine compute_aniso_tensor(s_ij,ie)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
      include 'RANS_BUO'

      integer lxyz
      parameter(lxyz=lx1*ly1*lz1)
      
      real s_ij(lxyz,6)

      common /sijoij/    sij  (lx1*ly1*lz1,6,lelv)
     $                , oij  (lx1*ly1*lz1,lelv,3)
      real sij, oij

      integer ie,i,j

      real Pr_t,mu_t0,rho,rans_mut

      real uxr,uxs,uyr,uys
      common /utmp1/ uxr(lxyz),uxs(lxyz),
     $               uyr(lxyz),uys(lxyz)  

      Pr_t = coeffs(1)
      
      call local_grad2(uxr,uxs,vx,lx1-1,ie,dxm1,dytm1)
      call local_grad2(uyr,uys,vy,lx1-1,ie,dxm1,dytm1)

      do i=1,lxyz
        mu_t = rans_mut(i,1,1,ie)
        rho = vtrans(i,1,1,ie,1)
        k = t(i,1,1,ie,ifld_k-1)
        mu_t0 = 0.0
        if(k.ne.0.)mu_t0 = mu_t/k

        if(if3d)then
          ! s_ij(i,1) = (2./3.)*rho - mu_t0*sij(i,1,ie)
          ! s_ij(i,2) = (2./3.)*rho - mu_t0*sij(i,2,ie)
          ! s_ij(i,3) = (2./3.)*rho - mu_t0*sij(i,3,ie)
          ! s_ij(i,4) = -mu_t0*sij(i,4,ie)
          ! s_ij(i,5) = -mu_t0*sij(i,5,ie)
          ! s_ij(i,6) = -mu_t0*sij(i,6,ie)
          ! call cmult(s_ij,mu_t/Pr_t,6)
          ! call cadd(s_ij,cpfld(ifield,1),6)
          continue
        else
          s_ij(i,1) = (2./3.)*rho - mu_t0*2.0*uxr(i)
          s_ij(i,2) = (2./3.)*rho - mu_t0*2.0*uys(i)
          s_ij(i,3) = -mu_t0*(uxs(i)+uyr(i))
          do j=1,3
            s_ij(i,j) = s_ij(i,j)*mu_t*0.3
          enddo
          do j=1,2
            s_ij(i,j) = s_ij(i,j)+cpfld(2,1)
          enddo
          ! call cmult(s_ij(i,1),mu_t/Pr_t,3)
          ! call cadd(s_ij(i,1),cpfld(ifield,1),3)

          ! s_ij(i,1) = cpfld(ifield,1)+mu_t/Pr_t
          ! s_ij(i,2) = cpfld(ifield,1)+mu_t/Pr_t
          ! s_ij(i,3) = 0.0 
        endif
      enddo
      return
      end
c--------------------------------------------------------------      
      subroutine addcol5(a,b,c,d,e,n)
      real a(1),b(1),c(1),d(1),e(1)
      integer n
      
      do i=1,n
        a(i) = a(i) + b(i)*c(i)*d(i)*e(i)
      enddo

      return
      end
c--------------------------------------------------------------      
      subroutine apply_buo_geom
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_BUO'

      integer nxyz,ntot
      integer ix,iy,iz,iel

      integer icalld
      save icalld
      data icalld /0/

      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelt

      if(ifggdh)then
        if(ifield.eq.2)then
          if(icalld.eq.0)then 
            call store_geom
            icalld = 1
          endif
          !apply geom factors for temperature
          call compute_ggdh_geom
          ! continue
        elseif(ifield.eq.3)then
          !replace geom factors
          call copy(g1m1,gms,ntot*6)
        endif
      endif

      return
      end
c--------------------------------------------------------------      
      subroutine eval_vgrad
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z
      common /scruz/ u_x (lx1,ly1,lz1,lelv),u_y (lx1,ly1,lz1,lelv)
     $              ,u_z (lx1,ly1,lz1,lelv),v_x (lx1,ly1,lz1,lelv)
     $              ,v_y (lx1,ly1,lz1,lelv),v_z (lx1,ly1,lz1,lelv)
     $              ,w_x (lx1,ly1,lz1,lelv),w_y (lx1,ly1,lz1,lelv)
     $              ,w_z (lx1,ly1,lz1,lelv)


      call gradm1  (u_x, u_y, u_z, vx)
      call gradm1  (v_x, v_y, v_z, vy)
      if(if3d) call gradm1  (w_x, w_y, w_z, vz)

      call opcolv  (u_x, u_y, u_z,bm1)
      call opcolv  (v_x, v_y, v_z,bm1)
      if(if3d) call opcolv  (w_x, w_y, w_z,bm1)

      call opdssum (u_x, u_y, u_z)
      call opdssum (v_x, v_y, v_z)
      if(if3d) call opdssum (w_x, w_y, w_z)

      call opcolv  (u_x, u_y, u_z,binvm1)
      call opcolv  (v_x, v_y, v_z,binvm1)
      if(if3d) call opcolv  (w_x, w_y, w_z,binvm1)
      return
      end
