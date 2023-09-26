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

      ! if(ifaxis.and.id.ne.0)then
      !   if(nid.eq.0)write(6,*)
      ! &   "ERROR:Axisymmetric BC not yet supported with Buoyancy"
      !   call exitt
      ! endif

      if(.not.ifrans_ktau_stndrd .and.
     &   .not.ifrans_komg_stndrd .and.
     &   .not.ifrans_ktauSST_stndrd .and.
     &   .not.ifrans_komgSST_stndrd)then
        if(nid.eq.0)then
          write(6,*)
     &    "ERROR: Buoyancy is only supported with:"
          write(6,*)"reg. k-omega     : m_id = 0"
          write(6,*)"reg. k-omega SST : m_id = 2"
          write(6,*)"std. k-tau       : m_id = 4"
          write(6,*)"std. k-tau SST   : m_id = 6"
        endif
        call exitt
      endif

      buo_gvec(1) = g(1)
      buo_gvec(2) = g(2)
      buo_gvec(3) = g(3)
 
      Cs_buo = 0.3
      Prt_buo = 0.9
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
      real rans_mut


      mu_t = rans_mut(ix,iy,iz,iel)

      if(.not.ifggdh)then
        buo_diff = cpfld(ifield,1) + mu_t/Prt_buo
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
        if(ifrans_ktauSST_stndrd) call rans_ktauSST_stndrd_buo
        if(ifrans_komgSST_stndrd) call rans_komgSST_stndrd_buo
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
        if(ifrans_ktauSST_stndrd) call rans_ktauSST_stndrd_buo
        if(ifrans_komgSST_stndrd) call rans_komgSST_stndrd_buo
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
        if(ifrans_ktauSST_stndrd) call rans_ktauSST_stndrd_buo
        if(ifrans_komgSST_stndrd) call rans_komgSST_stndrd_buo
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
        if(ifrans_ktauSST_stndrd) call rans_ktauSST_stndrd_buo
        if(ifrans_komgSST_stndrd) call rans_komgSST_stndrd_buo
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

      real alpinf_str,alp_inf
      real rho,tau,mu_t0,dotsrc
      real alp_str,alpha,gamm,G_w,G_k

      integer icalld
      save icalld
      data icalld /0/
      
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
        ! if(ifggdh) call comp_sij
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
            if(.not.ifggdh)dotsrc = dotsrc/Prt_buo

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
      subroutine rans_ktauSST_stndrd_buo
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
      include 'RANS_BUO'

      integer lxyz
      parameter(lxyz=lx1*ly1*lz1)

      real xflux(lxyz),yflux(lxyz),zflux(lxyz)
      integer e,i

      real alpinf_str,alp_inf
      real rho,tau,mu_t0,dotsrc
      real alp_str,alpha,gamm,G_w,G_k

      real Pr_t,sigk1,sigom1
      real r_k,beta1,alp0_str
      real beta_str,gamma1,r_b,akk
      real alpha_0,r_w
      real alp1,beta2,sigk2,sigom2,gamma2
      
      real arg1,arg1_1,arg1_2,arg2,arg2_1,arg2_2
      real argf1,argf2,argn
      real fun1,fun2,gamma
      real tiny,tinysst,denom,st_magn
      real xk,yw,ywm1,ywm2

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)
      real St_mag2,Om_mag2,OiOjSk,DivQ

      integer icalld
      save icalld
      data icalld /0/

c Turbulent viscosity constants
        Pr_t         = coeffs( 1)
        sigk1        = coeffs( 2)
        sigom1       = coeffs( 3)

c Low Reynolds number correction constants
        alpinf_str   = coeffs( 4)
        r_k          = coeffs( 5)
        beta1        = coeffs( 6)
        alp0_str     = coeffs( 7)

c Dissipation of K constants
        beta_str     = coeffs( 8)
        gamma1       = coeffs( 9)
        r_b          = coeffs(10)
        akk          = coeffs(11)

c Production of omega constants
        alpha_0      = coeffs(12)
        r_w          = coeffs(13)

c Dissipation of omega constants
c         beta_0 = defined earlier	     

        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional SST and k and epsilon constants
        alp1         = coeffs(17)
        beta2        = coeffs(18)
        sigk2        = coeffs(19)
        sigom2       = coeffs(20)
        gamma2       = coeffs(21)

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
        ! if(ifggdh) call comp_sij
      endif

      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      do e=1,nelv
         if(ifsgdh)call flux_sgdh_compute(xflux,yflux,zflux,e)
         if(ifggdh)call flux_ggdh_compute(xflux,yflux,zflux,e)

         do i=1,lxyz
           rho = vtrans(i,1,1,e,1)
           mu = vdiff(i,1,1,e,1)
           nu = mu/rho

           k = t(i,1,1,e,ifld_k-1)
           tau = t(i,1,1,e,ifld_omega-1)

           St_magn = sqrt(St_mag2(i,e))

           yw     = ywd  (i,1,1,e)
           ywm1   = ywdm1(i,1,1,e)
           ywm2   = ywm1*ywm1
           arg2_1 =     sqrt(k) * ywm1 * tau / beta_str
           arg2_2 =    500.0*nu * ywm2 * tau
           arg2   = 2.0*arg2_1
           argF2  =     sqrt(k)*yw/(500.0*nu*beta_str)
           if(2.0*argF2 .le. 1.0) arg2   = arg2_2
           Fun2   = tanh(arg2 * arg2)

           tinySST= 1.0e-10
           arg1_1 = arg2_1
           if(    argF2 .le. 1.0) arg1_1   = arg2_2
           arg1_2 =     4.0 * rho * sigom2 * k * ywm2 / tinySST
           argF1  = tinySST * tau /(2.0 * rho * sigom2)
           if(xk .gt. argF1) arg1_2 =   2.0 * k * tau * ywm2 / xk
           arg1   = min(    arg1_1, arg1_2)
           Fun1   = tanh(arg1 * arg1 * arg1 * arg1)

           mu_t   = rho * k * tau
           mu_t0  = rho * tau
           argn   = Fun2*St_magn ! this can also be Om_magn
           if(alp1.le.(argn*tau)) then
             mu_t   = 0.0
             mu_t0  = 0.0
             if(argn.ne.0.)then
               mu_t   = rho * alp1 * k/argn
               mu_t0  = rho * alp1 /argn
             endif
             denom  = argn/ alp1
           else
             denom  = 0.
             if(tau.ne.0.) denom  = 1.0/tau
           endif

           dotsrc = buo_gvec(1)*xflux(i) + buo_gvec(2)*yflux(i)
           if(if3d) dotsrc = dotsrc + buo_gvec(3)*zflux(i)
           if(.not.ifggdh)dotsrc = dotsrc/Prt_buo

           if(ifrans_diag)then
             ksrc_buo(i,1,1,e) = 0.0
             kdiag_buo(i,1,1,e) = -mu_t0 * dotsrc
           else
             ksrc_buo(i,1,1,e) = mu_t * dotsrc
             kdiag_buo(i,1,1,e) = 0.0
           endif

           gamma = Fun1 * gamma1 + (1.0 - Fun1) * gamma2

           G_w = rho*gamma*tau*dotsrc

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

      real alpinf_str,alp_inf
      real rho,omega,mu_t0,dotsrc
      real alp_str,alpha,gamm,G_w
      real tiny

      integer icalld
      save icalld
      data icalld /0/
      
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
        ! if(ifggdh) call comp_sij
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
            if(.not.ifggdh)dotsrc = dotsrc/Prt_buo

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
      subroutine rans_komgSST_stndrd_buo
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
      include 'RANS_BUO'

      integer lxyz
      parameter(lxyz=lx1*ly1*lz1)

      real xflux(lxyz),yflux(lxyz),zflux(lxyz)
      integer e,i

      real alpinf_str,alp_inf
      real rho,omega,mu_t0,dotsrc
      real alp_str,alpha,gamm,G_w

      real Pr_t,sigk1,sigom1
      real r_k,beta1,alp0_str
      real beta_str,gamma1,r_b,akk
      real alpha_0,r_w
      real alp1,beta2,sigk2,sigom2,gamma2
      
      real arg1,arg1_1,arg1_2,arg2,arg2_1,arg2_2
      real argf1,argf2,argn
      real fun1,fun2,gamma
      real tiny,tinysst,denom,st_magn
      real xk,yw,ywm1,ywm2

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)
      real St_mag2, Om_mag2, OiOjSk, DivQ

      integer icalld
      save icalld
      data icalld /0/
      
      Pr_t         = coeffs( 1)
      sigk1        = coeffs( 2)
      sigom1       = coeffs( 3)

      alpinf_str   = coeffs( 4)
      r_k          = coeffs( 5)
      beta1        = coeffs( 6)
      alp0_str     = coeffs( 7)

      beta_str     = coeffs( 8)
      gamma1       = coeffs( 9)
      r_b          = coeffs(10)
      akk          = coeffs(11)

      alpha_0      = coeffs(12)
      r_w          = coeffs(13)

      kv_min       = coeffs(14)
      omeg_max     = coeffs(15)
      tiny         = coeffs(16)

      alp1         = coeffs(17)
      beta2        = coeffs(18)
      sigk2        = coeffs(19)
      sigom2       = coeffs(20)
      gamma2       = coeffs(21)

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
        ! if(ifggdh) call comp_sij
      endif

      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      do e=1,nelv
         if(ifsgdh)call flux_sgdh_compute(xflux,yflux,zflux,e)
         if(ifggdh)call flux_ggdh_compute(xflux,yflux,zflux,e)

         do i=1,lxyz
            rho = vtrans(i,1,1,e,1)
            mu = vdiff(i,1,1,e,1)
            nu = mu/rho

            k = t(i,1,1,e,ifld_k-1)
            omega = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e)

            St_magn = sqrt(St_mag2(i,e))

            yw     = ywd  (i,1,1,e)
            ywm1   = ywdm1(i,1,1,e)
            ywm2   = ywm1*ywm1
            arg2_1 =     sqrt(k) * ywm1 / omega / beta_str
            arg2_2 =          500.0*nu * ywm2 / omega
            arg2   = 2.0*arg2_1
            argF2  =     sqrt(k)*yw/(500.0*nu*beta_str)
            if(2.0*argF2 .le. 1.0) arg2   = arg2_2
            Fun2   = tanh(arg2 * arg2)

            tinySST= 1.0e-10
            arg1_1 = arg2_1
            if(    argF2 .le. 1.0) arg1_1   = arg2_2
            arg1_2 =     4.0 * rho * sigom2 * k * ywm2 / tinySST
            argF1  = tinySST * omega /(2.0 * rho * sigom2)
            if(xk .gt. argF1) arg1_2 =   2.0 * k * omega * ywm2 / xk
            arg1   = min(    arg1_1, arg1_2)
            Fun1   = tanh(arg1 * arg1 * arg1 * arg1)

            mu_t   = rho * k/(omega + tiny)
            argn   = Fun2*St_magn ! this can also be Om_magn
            if(omega.le.argn/alp1) then
              mu_t   = rho * alp1 * k/argn
              denom  = argn/ alp1
            else
              denom  = omega
            endif

            dotsrc = buo_gvec(1)*xflux(i) + buo_gvec(2)*yflux(i)
            if(if3d) dotsrc = dotsrc + buo_gvec(3)*zflux(i)
            if(.not.ifggdh)dotsrc = dotsrc/Prt_buo

            if(ifrans_diag)then
              ksrc_buo(i,1,1,e) = mu_t * dotsrc
              kdiag_buo(i,1,1,e) = 0.0
            else
              ksrc_buo(i,1,1,e) = mu_t * dotsrc
              kdiag_buo(i,1,1,e) = 0.0
            endif

            gamma = Fun1 * gamma1 + (1.0 - Fun1) * gamma2
            G_w = rho*gamma*dotsrc

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
          s11 = (2./3.)*rho - mu_t0*sij(i,ie,1)
          s22 = (2./3.)*rho - mu_t0*sij(i,ie,2)
          s33 = (2./3.)*rho - mu_t0*sij(i,ie,3)
          s12 = -mu_t0*sij(i,ie,4)
          s23 = -mu_t0*sij(i,ie,5)
          s13 = -mu_t0*sij(i,ie,6)
          xflux(i) = Cs_buo*(s11*gtx+s12*gty+s13*gtz)
          yflux(i) = Cs_buo*(s12*gtx+s22*gty+s23*gtz)
          zflux(i) = Cs_buo*(s13*gtx+s23*gty+s33*gtz)
        elseif(ifaxis)then  !To Do
          s11 = (2./3.)*rho - mu_t0*sij(i,ie,1)
          s22 = (2./3.)*rho - mu_t0*sij(i,ie,2)
          s12 = -mu_t0*sij(i,ie,4)
          xflux(i) = Cs_buo*(s11*gtx+s12*gty)
          yflux(i) = Cs_buo*(s12*gtx+s22*gty)
        else
          s11 = (2./3.)*rho - mu_t0*sij(i,ie,1)
          s22 = (2./3.)*rho - mu_t0*sij(i,ie,2)
          s12 = -mu_t0*sij(i,ie,3)
          xflux(i) = Cs_buo*(s11*gtx+s12*gty)
          yflux(i) = Cs_buo*(s12*gtx+s22*gty)
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
      ntot = nxyz*lelt

      !Memory of all geom factors is contiguous
      !Be careful ntot should be multiple of lelt not nelt
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

      integer ir1_3d(6),ir2_3d(6)
      save ir1_3d,ir2_3d
      data ir1_3d /1,1,2,2,3,3/
      data ir2_3d /1,2,2,3,3,1/

      integer q_3d(3,3)
      save q_3d
      data q_3d /1,4,6,4,2,5,6,5,3/

      integer ir1,ir2,i2

      integer ntot

      if(icalld.eq.0)then
        !Store weighted jacobian
        do ie=1,lelt
          call invcol3(wj(1,ie),w3m1,jacm1(1,1,1,ie),lxyz)
        enddo
        icalld = 1
      endif

      ntot=lxyz*lelt

      call rzero(grs,ntot*6)

      !TO DO: delete call from above
      call comp_sij_buo

      do ie=1,lelt
        call compute_aniso_tensor(s_ij,ie)
        if(if3d)then
          do ii=1,6
            ir1 = ir1_3d(ii)
            ir2 = ir2_3d(ii)
            i2 = ii
            do j=1,3
              do i=1,3
                jj = q_3d(i,j)
                call addcol5(grs(1,ie,i2),rxg(1,ie,ir1,i),
     &          rxg(1,ie,ir2,j),wj(1,ie),s_ij(1,jj),lxyz)
              enddo
            enddo
          enddo
        else
          do ii=1,3
            ir1 = ir1_2d(ii)
            ir2 = ir2_2d(ii)
            i2 = ii
            if(ii.eq.3)i2=i2+1 !since axhelm uses g4m1 instead of g3m1
            do j=1,2
              do i=1,2
                jj = q_2d(i,j)
                call addcol5(grs(1,ie,i2),rxg(1,ie,ir1,i),
     &          rxg(1,ie,ir2,j),wj(1,ie),s_ij(1,jj),lxyz)
              enddo
            enddo
          enddo
        endif
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

      integer ie,i,j

      real mu_t0,rho,rans_mut

      do i=1,lxyz
        mu_t = rans_mut(i,1,1,ie)
        rho = vtrans(i,1,1,ie,1)
        k = t(i,1,1,ie,ifld_k-1)
        mu_t0 = 0.0
        if(k.ne.0.)mu_t0 = mu_t/k

        if(if3d)then
          s_ij(i,1) = (2./3.)*rho - mu_t0*sij(i,ie,1)
          s_ij(i,2) = (2./3.)*rho - mu_t0*sij(i,ie,2)
          s_ij(i,3) = (2./3.)*rho - mu_t0*sij(i,ie,3)
          s_ij(i,4) = -mu_t0*sij(i,ie,4)
          s_ij(i,5) = -mu_t0*sij(i,ie,5)
          s_ij(i,6) = -mu_t0*sij(i,ie,6)
          do j=1,6
            s_ij(i,j) = s_ij(i,j)*mu_t*Cs_buo
          enddo
          do j=1,3
            s_ij(i,j) = s_ij(i,j)+cpfld(2,1)
          enddo
        else
          s_ij(i,1) = (2./3.)*rho - mu_t0*sij(i,ie,1)
          s_ij(i,2) = (2./3.)*rho - mu_t0*sij(i,ie,2)
          s_ij(i,3) = -mu_t0*sij(i,ie,3)
          do j=1,3
            s_ij(i,j) = s_ij(i,j)*mu_t*Cs_buo
          enddo
          do j=1,2
            s_ij(i,j) = s_ij(i,j)+cpfld(2,1)
          enddo
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
      ntot = nxyz*lelt

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
          !Be careful - ntot should be multiple of lelt not nelt
          call copy(g1m1,gms,ntot*6)
        endif
      endif

      return
      end
c--------------------------------------------------------------      
      subroutine comp_sij_buo
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_BUO'
      
c
      integer e
      logical iflmc
c
      integer lxyz
      parameter(lxyz=lx1*ly1*lz1)

      real ur(lxyz),us(lxyz),ut(lxyz)
     $    ,vr(lxyz),vs(lxyz),vt(lxyz)
     $    ,wr(lxyz),ws(lxyz),wt(lxyz)

      real j ! Inverse Jacobian

      integer n, nxyz, nij, i
      
      real onethird,r,trs

      iflmc = .false. !low-Mach correction - off for now

      nij = 3
      if(if3d.or.ifaxis)nij = 6

      n    = lx1-1      ! Polynomial degree
      nxyz = lx1*ly1*lz1

      if (if3d) then     ! 3D CASE
       do e=1,nelv
        call local_grad3(ur,us,ut,vx,N,e,dxm1,dxtm1)
        call local_grad3(vr,vs,vt,vy,N,e,dxm1,dxtm1)
        call local_grad3(wr,ws,wt,vz,N,e,dxm1,dxtm1)

        do i=1,nxyz

         j = jacmi(i,e)

         sij(i,e,1) = j*  ! du/dx + du/dx
     $   2*(ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e)+ut(i)*txm1(i,1,1,e))

         sij(i,e,2) = j*  ! dv/dy + dv/dy
     $   2*(vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e)+vt(i)*tym1(i,1,1,e))

         sij(i,e,3) = j*  ! dw/dz + dw/dz
     $   2*(wr(i)*rzm1(i,1,1,e)+ws(i)*szm1(i,1,1,e)+wt(i)*tzm1(i,1,1,e))

         sij(i,e,4) = j*  ! du/dy + dv/dx
     $   (ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e)+ut(i)*tym1(i,1,1,e) +
     $    vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e)+vt(i)*txm1(i,1,1,e) )

         sij(i,e,5) = j*  ! dv/dz + dw/dy
     $   (wr(i)*rym1(i,1,1,e)+ws(i)*sym1(i,1,1,e)+wt(i)*tym1(i,1,1,e) +
     $    vr(i)*rzm1(i,1,1,e)+vs(i)*szm1(i,1,1,e)+vt(i)*tzm1(i,1,1,e) )

         sij(i,e,6) = j*  ! du/dz + dw/dx
     $   (ur(i)*rzm1(i,1,1,e)+us(i)*szm1(i,1,1,e)+ut(i)*tzm1(i,1,1,e) +
     $    wr(i)*rxm1(i,1,1,e)+ws(i)*sxm1(i,1,1,e)+wt(i)*txm1(i,1,1,e) )
        enddo
       enddo

      elseif (ifaxis) then  ! AXISYMMETRIC CASE  

c
c        Notation:                       ( 2  x  Acheson, p. 353)
c                     Cylindrical
c            Nek5k    Coordinates
c
c              x          z
c              y          r
c              z          theta
c

         do e=1,nelv
            call setaxdy ( ifrzer(e) )  ! change dytm1 if on-axis
            call local_grad2(ur,us,vx,N,e,dxm1,dytm1)
            call local_grad2(vr,vs,vy,N,e,dxm1,dytm1)
            call local_grad2(wr,ws,vz,N,e,dxm1,dytm1)

            do i=1,nxyz
               j = jacmi(i,e)
               r = ym1(i,1,1,e)                              ! Cyl. Coord:

               sij(i,e,1) = j*  ! du/dx + du/dx              ! e_zz
     $           2*(ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e))

               sij(i,e,2) = j*  ! dv/dy + dv/dy              ! e_rr
     $           2*(vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e))

               if (r.gt.0) then                              ! e_@@
                  sij(i,e,3) = 2.*vy(i,1,1,e)/r  ! v / r  ! corrected AT: factor of 2, 10/30/18
               else
                  sij(i,e,3) = j*  ! L'Hopital's rule: e_@@ = 2dv/dr
     $            2*(vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e))
               endif

               sij(i,e,4) = j*  ! du/dy + dv/dx             ! e_zr
     $            ( ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e) +
     $              vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e) )

               if (r.gt.0) then                             ! e_r@
                  sij(i,e,5) = j*  ! dw/dy 
     $              ( wr(i)*rym1(i,1,1,e)+ws(i)*sym1(i,1,1,e) )
     $              - vz(i,1,1,e) / r
               else
                  sij(i,e,5) = 0
               endif

               sij(i,e,6) = j*  ! dw/dx                     ! e_@z
     $            ( wr(i)*rxm1(i,1,1,e)+ws(i)*sxm1(i,1,1,e) )
            enddo
         enddo

      else              ! 2D CASE

         do e=1,nelv
            call local_grad2(ur,us,vx,N,e,dxm1,dxtm1)
            call local_grad2(vr,vs,vy,N,e,dxm1,dxtm1)

            do i=1,nxyz
               j = jacmi(i,e)

               sij(i,e,1) = j*  ! du/dx + du/dx
     $           2*(ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e))

               sij(i,e,2) = j*  ! dv/dy + dv/dy
     $           2*(vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e))

               sij(i,e,3) = j*  ! du/dy + dv/dx
     $           (ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e) +
     $            vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e) )

            enddo
         enddo
      endif

      if(iflomach .and. iflmc) then

        onethird = 1./3.
        if(if3d .or. ifaxis) then
          do e=1,nelv
            do i=1,nxyz
             trS = sij(i,e,1) + sij(i,e,2) + sij(i,e,3) ! 2(du/dx + dv/dy + dw/dz or v/r) in axisym
             sij(i,e,1) = sij(i,e,1) - onethird*trS     ! 2S - (2/3)Q = S'-(1/3)tr(S')
            enddo
          enddo
        else
          do e=1,nelv
            do i=1,nxyz
             trS = sij(i,e,1) + sij(i,e,2)              ! 2(du/dx + dv/dy)
             sij(i,e,1) = sij(i,e,1) - onethird*trS     ! 2S - (2/3)Q = S'-(1/3)tr(S')
            enddo
           enddo
        endif

      endif

      call dssum_sij(nij)

      return
      end
c-----------------------------------------------------------------------
      subroutine dssum_sij(nij)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_BUO'

      integer nij,nxyz,ntot,i

      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv

      do i=1,nij
        call col2(sij(1,1,i),bm1,ntot)
        call dssum(sij(1,1,i),lx1,ly1,lz1)
        call col2(sij(1,1,i),binvm1,ntot)
      enddo

      return
      end
c-----------------------------------------------------------------------
