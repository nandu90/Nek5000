c---------------------------------------------------------------------      
      subroutine getRedistVel(ifld)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
c
      real tmp(lx1,ly1,lz1,lelt)
c
      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv

      call gradm1(vx,vy,vz,t(1,1,1,1,ifld-1))
      call opcolv(vx,vy,vz,bm1)
      call opdssum(vx,vy,vz)
      call opcolv(vx,vy,vz,binvm1)
      
      call col3(tmp,vx,vx,ntot)
      call addcol3(tmp,vy,vy,ntot)
      if(if3d)call addcol3(tmp,vz,vz,ntot)

      do i=1,ntot
         tmp(i,1,1,1) = max(1.0e-15,sqrt(tmp(i,1,1,1)))
      enddo

      call opicol2(vx,vy,vz,tmp,tmp,tmp)
      call opcolv(vx,vy,vz,signls)

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
      
