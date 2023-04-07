c--------------------------------------------------------------      
      subroutine rans_buo_init(id, g)
      implicit none
      include 'SIZE'
      include 'RANS_KOMG'
      include 'RANS_BUO'

      integer id
      real g(3)

      ifsgdh = .false.
      ifggdh = .false.
      if(id.eq.1) ifsgdh = .true.
      if(id.eq.2) ifggdh = .true.

      if(.not.ifsgdh .or. ifggdh)then
        if(nid.eq.0)write(6,*)
     &   "ERROR: Only SGDH supported for now"
        call exitt
      endif

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
      include 'RANS_KOMG'
      include 'RANS_BUO'

      integer ix,iy,iz,iel
      real Pr_t, rans_mut

      Pr_t = coeffs(1)
      mu_t = rans_mut(ix,iy,iz,iel)

      if(.not.ifggdh)then
        buo_diff = mu_t/Pr_t
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
      real alp_str,alpha,gamm,G_w

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
        ! Evaluate temp gradient once per time step
        call eval_tgrad
      endif

      do e=1,nelv
         if(ifsgdh)call flux_sgdh_compute(xflux,yflux,zflux,e)

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
        ! Evaluate temp gradient once per time step
        call eval_tgrad
      endif

      do e=1,nelv
         if(ifsgdh)call flux_sgdh_compute(xflux,yflux,zflux,e)

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



