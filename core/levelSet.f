c---------------------------------------------------------------------
      subroutine tlsconv(ifld)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
c
      real tmp(lx1,ly1,lz1,lelv)

      common /tlstemp/ tls
      
      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv

c     get LS sign
      call getLSsgn(ifld,eps_redist)
      
c     Redistancing velocity
      call gradm1(tlsvx,tlsvy,tlsvz,t(1,1,1,1,ifld-1))
      call opcolv(tlsvx,tlsvy,tlsvz,bm1)
      call opdssum(tlsvx,tlsvy,tlsvz)
      call opcolv(tlsvx,tlsvy,tlsvz,binvm1)
      
      call col3(tmp,tlsvx,tlsvx,ntot)
      call addcol3(tmp,tlsvy,tlsvy,ntot)
      if(if3d)call addcol3(tmp,tlsvz,tlsvz,ntot)

      do i=1,ntot
         tmp(i,1,1,1) = max(1.0e-15,sqrt(tmp(i,1,1,1)))
      enddo

      call opicol2(tlsvx,tlsvy,tlsvz,tmp,tmp,tmp)
      call opcolv(tlsvx,tlsvy,tlsvz,signls)

c     convective term
      call convect_new(tlsadv,t(1,1,1,1,ifld-1),.false.,
     $     tlsvx,tlsvy,tlsvz,.false.)
      
      call invcol2(tlsadv,bm1,ntot)
      return
      end
c---------------------------------------------------------------------
      real function q_tlsconv(ix,iy,iz,ie)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
c
      q_tlsconv = -tlsadv(ix,iy,iz,ie)

      return
      end
c---------------------------------------------------------------------
      subroutine getElementLength
c      
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
c
      integer icalld
      save icalld
      data icalld /0/

      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv

      if(icalld .eq. 0)then
         call copy(lsH,jacm1,ntot)
         do i=1,ntot
            if(if3d)then
               lsH(i,1,1,1) = lsH(i,1,1,1)**(1.0/3.0)
            else
               lsH(i,1,1,1) = lsH(i,1,1,1)**(1.0/2.0)
            endif
         enddo
         call cmult(lsH,2.0,ntot)
         icalld = 1
      else
         return
      endif

      return
      end
c---------------------------------------------------------------------
      subroutine getLSsgn(ifld,eps)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
c
      real eps
            
      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv

      call getElementLength
      
      do i=1,ntot
         signls(i,1,1,1) = tanh(2.0*PI*t(i,1,1,1,ifld-1)
     $        /(eps*lsH(i,1,1,1)))
      enddo

      return
      end
c---------------------------------------------------------------------
      real function sgnfuncLS(ix,iy,iz,ie)
c
      include 'SIZE'
      include 'LVLSET'

      sgnfuncLS = signls(ix,iy,iz,ie)

      end
c---------------------------------------------------------------------      
      
