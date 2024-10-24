      subroutine prc_omesh()
        include 'SIZE'
        include 'TSTEP'
        include 'INPUT'
        include 'CTIMER'
  
        INCLUDE 'GEOM'
        INCLUDE 'SOLN'
        INCLUDE 'TOPOL'
        INCLUDE 'PRECIC'

        character cb*3 !type of bc
        integer prctn !counter for the new vertices on outer Boundary
        real tol1, tol2, tol3 !tolerance for vertices equality
        prctn = 0
        prcnve = 0
        NEL   = NELFLD(1)
        nfaces = 6
        DO IE = 1,NEL  !iterate on each face of each element
        DO IFACE=1,NFACES
           CB = CBC(IFACE,IE,1)
           if (CB.eq.'v  '.or.CB.eq.'snl') then !if the face is a boundary
           call facind (kx1,kx2,ky1,ky2,kz1, kz2,lx1,ly1,lz1,iface)
           do 100 iz=kz1,kz2 !iterate on each vertex of the face
           do 100 iy=ky1,ky2
           do 100 ix=kx1,kx2
              if(prctn.eq.0) then !case of the first vertex
              prcvrt((prctn)*3) = xm1(ix,iy,iz,ie)
              prcvrt((prctn)*3+1) = ym1(ix,iy,iz,ie)
              prcvrt((prctn)*3+2) = zm1(ix,iy,iz,ie)
              prcvid(prctn) = prctn
              mpgprc(ix,iy,iz,ie) = prctn
              prctn = prctn + 1
              prcnve = prcnve + 1
              goto 100
              endif
              do icnt = 0,prctn !parse all the already defined vertices
                 tol1 = abs(xm1(ix,iy,iz,ie) - prcvrt(icnt*3))
                 tol2 = abs(ym1(ix,iy,iz,ie) - prcvrt(icnt*3+1))
                 tol3 = abs(zm1(ix,iy,iz,ie) - prcvrt(icnt*3+2))
                 if(tol1.lt.1e-10.and.tol2.lt.1e-10.and. !if already existing, map this vertex to the one existing
     &             tol3.lt.1e-10) then
                 
                 mpgprc(ix,iy,iz,ie) = icnt
                 goto 100
                 endif
              enddo
              prcvrt(3*(prctn)) = xm1(ix,iy,iz,ie) !else, create a new boundary vertex
              prcvrt(3*(prctn)+1) = ym1(ix,iy,iz,ie)
              prcvrt(3*(prctn)+2) = zm1(ix,iy,iz,ie)
              prcvid(prctn) = prctn
              mpgprc(ix,iy,iz,ie) = prctn
              prctn = prctn + 1
              prcnve = prcnve + 1
  100  continue

           endif
        ENDDO
        ENDDO
      
      return
      end

c---------------------------------------------------------------------------------------------      
      subroutine setup_interp(nekcom)
         include 'size'
         include 'total'
         include 'precic'
        integer nekcom
        integer nxf, nyf, nzf
        real bb_t
        integer n, npt_max
        real tol
        real vmodif(0:omshdi*3-1)
        real hvpm
        integer i

        nxf = 2 * lx1
        nyf = 2 * ly1
        nzf = 2 * lz1
        bb_t = 0.01
        n = lx1*ly1*lz1*lelt
        npt_max = 128
        tol = 5e-13
        print *, np_global, 'np global'
        print *, nekcom
        call fgslib_findpts_setup(handle_u, nekcom, np_global, 3,
     &  xm1, ym1, zm1, lx1, ly1, lz1, nelt, nxf, nyf, nzf, bb_t,
     &  n, n, npt_max, tol)

      !we want to interpolate the Nek velocity onto a staggered VPM grid. The coordinates for the 
      !3 components of velocity are (with xi,yj,zk the Murphy vertex where the vorticity is stored):
      !for u: (xi, yj + h/2, zk + h/2)
      !for v: (xi+h/2, yj  , zk + h/2)
      !for u: (xi+ h/2, yj + h/2, zk )
      hvpm = 1./12.
      !modify the vertex coordinates to interpolate the x-component: adding h/2 to y and z
      do i = 0, omshdi-1
         vmodif(3*i)   = prcvr2(3*i)-1e-7
         vmodif(3*i+1) = prcvr2(3*i+1) + 0.5*hvpm
         vmodif(3*i+2) = prcvr2(3*i+2) + 0.5*hvpm
      enddo

        call fgslib_findpts(handle_u, rcodeu, 1, 
     &                    procu, 1,
     &                    elidu, 1,
     &                    rstu, 3,
     &                    distu, 1,
     &                    vmodif(0), 3,
     &                    vmodif(1), 3,
     &                    vmodif(2), 3, omshdi)
      print *, "iebushbi here"
   !    do i = 0, omshdi-1
   !       ! if(rcodeu(i).eq.2.and.(abs(vmodif(3*i)-1.125).lt.1.e-8)) then
   !       if(rcodeu(i).eq.2) then
   !          print *, vmodif(3*i),vmodif(3*i+1),vmodif(3*i+2),
   !   & distu(i)       
   !          print *, rstu(3*i), rstu(3*i+1), rstu(3*i+2)
   !       endif 
   !    enddo
      call fgslib_findpts_setup(handle_v, nekcom, np_global, 3,
     &  xm1, ym1, zm1, lx1, ly1, lz1, nelt, nxf, nyf, nzf, bb_t,
     &  n, n, npt_max, tol)

      do i = 0, omshdi-1
         vmodif(3*i)   = prcvr2(3*i) + 0.5*hvpm
         vmodif(3*i+1) = prcvr2(3*i+1) 
         vmodif(3*i+2) = prcvr2(3*i+2) + 0.5*hvpm
      enddo

        call fgslib_findpts(handle_v, rcodev, 1, 
     &                    procv, 1,
     &                    elidv, 1,
     &                    rstvv, 3,
     &                    distv, 1,
     &                    vmodif(0), 3,
     &                    vmodif(1), 3,
     &                    vmodif(2), 3, omshdi)

      
      call fgslib_findpts_setup(handle_w, nekcom, np_global, 3,
     &  xm1, ym1, zm1, lx1, ly1, lz1, nelt, nxf, nyf, nzf, bb_t,
     &  n, n, npt_max, tol)
      do i = 0, omshdi-1
         vmodif(3*i)   = prcvr2(3*i) + 0.5*hvpm
         vmodif(3*i+1) = prcvr2(3*i+1) + 0.5*hvpm
         vmodif(3*i+2) = prcvr2(3*i+2) 
      enddo

        call fgslib_findpts(handle_w, rcodew, 1, 
     &                    procw, 1,
     &                    elidw, 1,
     &                    rstw, 3,
     &                    distw, 1,
     &                    vmodif(0), 3,
     &                    vmodif(1), 3,
     &                    vmodif(2), 3, omshdi)
      return 
      end

c-----------------------------------------------------------------------------
      subroutine interpolate_u

         INCLUDE 'SIZE'
         INCLUDE 'TOTAL'
         include 'precic'
 
         
         call fgslib_findpts_eval(handle_u,prcwdt(0), 3,
     &                         rcodeu, 1, 
     &                         procu, 1, 
     &                         elidu, 1, 
     &                         rstu, 3, omshdi, 
     &                         vx)
 
         call fgslib_findpts_eval(handle_v,prcwdt(1), 3,
     &                         rcodev, 1, 
     &                         procv, 1, 
     &                         elidv, 1, 
     &                         rstvv, 3, omshdi, 
     &                         vy)

         call fgslib_findpts_eval(handle_w,prcwdt(2), 3,
     &                         rcodew, 1, 
     &                         procw, 1, 
     &                         elidw, 1, 
     &                         rstw, 3, omshdi, 
     &                         vz)
       return
       end 
c------------------------------------------------------------------------------------------------------------------------

      subroutine interpolate_u_old(omshdi, vert, wdt)

        INCLUDE 'SIZE'
        INCLUDE 'TOTAL'

        common /nekmpi/ nekcomm

        save inth_hpts
        integer nxf, nyf, nzf
        real bb_t
        integer n, npt_max
        real tol
        integer omshdi
        integer rcode(0:omshdi-1),proc(0:omshdi-1),elid (0:omshdi-1)
        real rst(0:omshdi*3-1)
        real dist(0:omshdi-1)
        real vert(0:omshdi*3-1)
        real wdt(0:omshdi*3-1)
        real vmodif(0:omshdi*3-1)
        real hvpm
        integer i

        nxf = 2 * lx1
        nyf = 2 * ly1
        nzf = 2 * lz1
        bb_t = 0.0
        n = lx1*ly1*lz1*lelt
        npt_max = 128
        tol = 5e-13
        print *, "here here"
        print *, inth_hpts

   !      call fgslib_findpts_setup(inth_hpts, nekcomm, np, 3,
   !   &  xm1, ym1, zm1, lx1, ly1, lz1, nelt, nxf, nyf, nzf, bb_t,
   !   &  n, n, npt_max, tol)



      !we want to interpolate the Nek velocity onto a staggered VPM grid. The coordinates for the 
      !3 components of velocity are (with xi,yj,zk the Murphy vertex where the vorticity is stored):
      !for u: (xi, yj + h/2, zk + h/2)
      !for v: (xi+h/2, yj  , zk + h/2)
      !for u: (xi+ h/2, yj + h/2, zk )

      hvpm = 1./12.
      !modify the vertex coordinates to interpolate the x-component: adding h/2 to y and z
      do i = 0, omshdi-1

         vmodif(3*i)   = vert(3*i)
         vmodif(3*i+1) = vert(3*i+1) + 0.5*hvpm
         vmodif(3*i+2) = vert(3*i+2) + 0.5*hvpm
      enddo

   !      call fgslib_findpts(inth_hpts, rcode, 1, 
   !   &                    proc, 1,
   !   &                    elid, 1,
   !   &                    rst, 3,
   !   &                    dist, 1,
   !   &                    vmodif(0), 3,
   !   &                    vmodif(1), 3,
   !   &                    vmodif(2), 3, omshdi)
        
   !      call fgslib_findpts_eval(inth_hpts,wdt(0), 3,
   !   &                         rcode, 1, 
   !   &                         proc, 1, 
   !   &                         elid, 1, 
   !   &                         rst, 3, omshdi, 
   !   &                         vx)
       !modify the vertex coordinates to interpolate the y-component: adding h/2 to x and z
      do i = 0, omshdi-1
         vmodif(3*i)   = vert(3*i) + 0.5*hvpm
         vmodif(3*i+1) = vert(3*i+1) 
         vmodif(3*i+2) = vert(3*i+2) + 0.5*hvpm
      enddo

   !      call fgslib_findpts(inth_hpts, rcode, 1, 
   !   &                    proc, 1,
   !   &                    elid, 1,
   !   &                    rst, 3,
   !   &                    dist, 1,
   !   &                    vmodif(0), 3,
   !   &                    vmodif(1), 3,
   !   &                    vmodif(2), 3, omshdi)

   !      call fgslib_findpts_eval(inth_hpts,wdt(1), 3,
   !   &                         rcode, 1, 
   !   &                         proc, 1, 
   !   &                         elid, 1, 
   !   &                         rst, 3, omshdi, 
   !   &                         vy)

       !modify the vertex coordinates to interpolate the z-component: adding h/2 to x and y
      do i = 0, omshdi-1
         vmodif(3*i)   = vert(3*i) + 0.5*hvpm
         vmodif(3*i+1) = vert(3*i+1) + 0.5*hvpm
         vmodif(3*i+2) = vert(3*i+2) 
      enddo

   !      call fgslib_findpts(inth_hpts, rcode, 1, 
   !   &                    proc, 1,
   !   &                    elid, 1,
   !   &                    rst, 3,
   !   &                    dist, 1,
   !   &                    vmodif(0), 3,
   !   &                    vmodif(1), 3,
   !   &                    vmodif(2), 3, omshdi)
   !      call fgslib_findpts_eval(inth_hpts,wdt(2), 3,
   !   &                         rcode, 1, 
   !   &                         proc, 1, 
   !   &                         elid, 1, 
   !   &                         rst, 3, omshdi, 
   !   &                         vz)
      return
      end 



c-----subroutine to assign x,y,z------------------------------------------------
      subroutine assign_coord(i,j,k,e,x_,y_,z_)
      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'GEOM'
      include 'DXYZ'
      include 'SOLN' 
      real x_,y_,z_ 
      integer i,j,k,e 
      x_ = xm1(i,j,k,e)
      y_ = ym1(i,j,k,e)
      z_ = zm1(i,j,k,e)
      return
      end
 
c-----temp subroutine to test condition-------------------------------------------
      subroutine cust_func(x_,y_,z_,u_,v_,w_,omx,omy,omz)

      real x_,y_,z_,u_,v_,w_,omx,omy,omz
      
      u_ = x_*y_*sin(z_)
      v_ = x_ - exp(y_)*cos(z_)
      w_ = cos(x_) +exp(x_)*y_*y_+z_ 
      omx = 2.*y_*exp(x_)-exp(y_)*sin(z_)
      omy = x_*y_*cos(z_) + sin(x_) - exp(x_)*y_*y_ 
      omz = 1. -x_*sin(z_)
      return
      end
c--------------------------------------------------------------------------------
      subroutine find_boundaries
      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'GEOM'
      include 'DXYZ'
      include 'SOLN'
      include 'PRECIC'    
      
      integer e
      common /nekcb/ cbcmesh
      character cb*3

      !variables for the new implementation
      do e = 1,nelt 
        !first, find the faces that have 'v' BC's, if any
        n_fbc = 0! number of faces that have a BC
        n_fnbc = 0! number of faces that do not have a BC

        cb = cbc(4,e,1) ! face x-
        if(cb.eq.'v  ') then 
          n_fbc = n_fbc + 1
          tmpl1(n_fbc,e) = 4
          tmpl3(4,e) = 1
        else 
          n_fnbc = n_fnbc + 1
          tmpl2(n_fnbc,e) = 4
          tmpl3(4,e) = 0
        endif 

        cb = cbc(2,e,1) ! face x+
        if(cb.eq.'v  ') then 
          n_fbc = n_fbc + 1
          tmpl1(n_fbc,e) = 2
          tmpl3(2,e) = 1
        else 
          n_fnbc = n_fnbc + 1
          tmpl2(n_fnbc,e) = 2
          tmpl3(2,e) = 0
        endif

        cb = cbc(1,e,1) ! face y-
        if(cb.eq.'v  ') then 
          n_fbc = n_fbc + 1
          tmpl1(n_fbc,e) = 1
          tmpl3(1,e) = 1
        else 
          n_fnbc = n_fnbc + 1
          tmpl2(n_fnbc,e) = 1
          tmpl3(1,e) = 0
        endif

        cb = cbc(3,e,1) ! face y+
        if(cb.eq.'v  ') then 
          n_fbc = n_fbc + 1
          tmpl1(n_fbc,e) = 3
          tmpl3(3,e) = 1
        else 
          n_fnbc = n_fnbc + 1
          tmpl2(n_fnbc,e) = 3
          tmpl3(3,e) = 0
        endif
        cb = cbc(5,e,1) ! face z-
        if(cb.eq.'v  ') then 
          n_fbc = n_fbc + 1
          tmpl1(n_fbc,e) = 5
          tmpl3(5,e) = 1
        else 
          n_fnbc = n_fnbc + 1
          tmpl2(n_fnbc,e) = 5
          tmpl3(5,e) = 0
        endif

        cb = cbc(6,e,1) ! face z+
        if(cb.eq.'v  ') then 
          n_fbc = n_fbc + 1
          tmpl1(n_fbc,e) = 6
          tmpl3(6,e) = 1
        else 
          n_fnbc = n_fnbc + 1
          tmpl2(n_fnbc,e) = 6
          tmpl3(6,e) = 0
        endif
        ! now we also put the faces w/out BC in tmpl1
        do i = 1, n_fnbc
          tmpl1(n_fbc + i,e) = tmpl2(i,e)
        enddo
        ! and in tmpl2, we put 1 if there's a BC and 0 otherwise
        do i = 1, n_fbc
          tmpl2(i,e) = 1
        enddo
        do i = 1, n_fnbc
          tmpl2(n_fbc+i,e) = 0
        enddo
        l_fbc(e) = n_fbc
        l_fnbc(e) = n_fnbc
        ! we also write the inverse maping of tmpl1 in tmpl4
        do i = 1, 2*ldim 
         tmpl4(tmpl1(i,e),e) = i 
        enddo 
   !      print *, tmpl1(1,e),tmpl1(2,e),tmpl1(3,e),tmpl1(4,e),
   !   &   tmpl1(5,e),tmpl1(6,e)
   !      print *, tmpl4(1,e),tmpl4(2,e),tmpl4(3,e),tmpl4(4,e),
   !   &   tmpl4(5,e),tmpl4(6,e)
      enddo




      return 

      end
c------------------------------------------------------------------------------------
      subroutine assign_corner
         include 'SIZE'
         include 'TSTEP'
         include 'INPUT'
         include 'CTIMER'
         include 'GEOM'
         include 'DXYZ'
         include 'SOLN'
         include 'PRECIC'

      x_corner(1) = 0.6
      x_corner(2) = 2.6
      x_corner(3) = 2.6
      x_corner(4) = 0.6
      x_corner(5) = 0.6
      x_corner(6) = 2.6
      x_corner(7) = 2.6
      x_corner(8) = 0.6
      y_corner(1) = 1.
      y_corner(2) = 1.
      y_corner(3) = 3.
      y_corner(4) = 3.
      y_corner(5) = 1.
      y_corner(6) = 1.
      y_corner(7) = 3.
      y_corner(8) = 3.
      z_corner(1) = 0.625
      z_corner(2) = 0.625
      z_corner(3) = 0.625
      z_corner(4) = 0.625
      z_corner(5) = 1.375
      z_corner(6) = 1.375
      z_corner(7) = 1.375
      z_corner(8) = 1.375
      ! x_corner(1) = 0.3
      ! x_corner(2) = 1.7
      ! x_corner(3) = 1.7
      ! x_corner(4) = 0.3
      ! x_corner(5) = 0.3
      ! x_corner(6) = 1.7
      ! x_corner(7) = 1.7
      ! x_corner(8) = 0.3
      ! y_corner(1) = 0.3
      ! y_corner(2) = 0.3
      ! y_corner(3) = 1.7
      ! y_corner(4) = 1.7
      ! y_corner(5) = 0.3
      ! y_corner(6) = 0.3
      ! y_corner(7) = 1.7
      ! y_corner(8) = 1.7
      ! z_corner(1) = 0.5
      ! z_corner(2) = 0.5
      ! z_corner(3) = 0.5
      ! z_corner(4) = 0.5
      ! z_corner(5) = 1.55
      ! z_corner(6) = 1.55
      ! z_corner(7) = 1.55
      ! z_corner(8) = 1.55

      return 
      end
      
c----------------------------------------------------------------------------------------
      subroutine compute_bc(kstep)
      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'GEOM'
      include 'DXYZ'
      include 'SOLN'
      include 'PRECIC'
      integer e,i,j,ifc,k
      integer ia, b1,b2
      real u_,v_,w_
      real omx,omy,omz
      real x_,y_,z_
      integer prcmpi

      !variables for the new implementation
      integer myface
      real dudy,dvdx,dudz,dwdx,dvdz,dwdy
      real ddx,ddy,ddz,dtol

      ! vars for the system of the corner
      real aa,bb,cc,dd,ee,ff
      real rr,ss,tt
      integer impose

      integer kstep


      dtol = 1e-12
         do e = 1,nelt! nelt 
            do ifc = 1,l_fbc(e)
               myface = tmpl1(ifc,e)
               if(myface.eq.4.or.myface.eq.2) then 
                  if(myface.eq.4) then 
                     ia = 1
                     b1 = 2
                     b2 = lx1
                  else 
                     ia = lx1 
                     b1 = 1
                     b2 = lx1-1
                  endif
                  do 210 i = 1, lx1 !loop on the y
                  do 211 k = 1, lx1 !loop on the z
                     call assign_coord(ia,i,k,e,x_,y_,z_)
                     do icorner = 1, 8
                      ddx = abs(x_-x_corner(icorner))
                      ddy = abs(y_-y_corner(icorner))
                      ddz = abs(z_-z_corner(icorner))
                      if((ddx.lt.dtol.and.ddy.lt.dtol).or.
     &                (ddx.lt.dtol.and.ddz.lt.dtol).or.
     &                (ddz.lt.dtol.and.ddy.lt.dtol)) then
                         prcmpi = mpgprc(ia,i,k,e)
                         usol(k,i,myface,e) = prcrdt(3*prcmpi)
                         vsol(k,i,myface,e) = prcrdt(3*prcmpi+1)
                         wsol(k,i,myface,e) = prcrdt(3*prcmpi+2)
                         goto 211
                      endif
                     enddo

                     dudz = 0.
                     dwdx = 0.
                     dudy = 0.
                     dvdx = 0.
                     if (kstep.gt.2) then
                     do j = b1,b2
                        dvdx = dvdx + dxm1(ia,j)*((2./1.)*vy(j,i,k,e)
     &                   -(1./1.)*vylag(j,i,k,e,1))
                        dwdx = dwdx + dxm1(ia,j)*((2./1.)*vz(j,i,k,e)
     &                   -(1./1.)*vzlag(j,i,k,e,1))
                     enddo
                  else 
                     do j = b1,b2
                        dvdx = dvdx + dxm1(ia,j)*vy(j,i,k,e)
                        dwdx = dwdx + dxm1(ia,j)*vz(j,i,k,e)
                     enddo
                  endif


                     do j = 1,lx1
                        prcmpi = mpgprc(ia,j,k,e)
                        dudy = dudy + dytm1(j,i)*prcrdt(3*prcmpi)
                        prcmpi = mpgprc(ia,i,j,e)
                        dudz = dudz + dytm1(j,k)*prcrdt(3*prcmpi)
                     enddo 
                     prcmpi = mpgprc(ia,i,k,e)
                     usol(k,i,myface,e) = prcrdt(3*prcmpi) 
                     omy = prcrdx(3*prcmpi+1)
                     omz = prcrdx(3*prcmpi+2)

                     vsol(k,i,myface,e) = (omz/jacmi((k-1)*lx1*lx1
     &              +(i-1)*lx1+ia,e)-dvdx*rxm1(ia,i,k,e) + dudy*
     &                 sym1(ia,i,k,e))/(dxm1(ia,ia)*rxm1(ia,i,k,e))
                     wsol(k,i,myface,e) = -(omy/jacmi((k-1)*lx1*lx1
     &            +(i-1)*lx1+ia,e)-dudz*tzm1(ia,i,k,e) + dwdx*
     &               rxm1(ia,i,k,e))/(dxm1(ia,ia)*rxm1(ia,i,k,e))
211            continue
210            continue

               else if(myface.eq.1.or.myface.eq.3) then 
                  if(myface.eq.1) then 
                     ia = 1
                     b1 = 2
                     b2 = lx1
                  else 
                     ia = lx1
                     b1 = 1
                     b2 = lx1-1 
                  endif
                  do 230 i = 1, lx1 !loop on the x
                  do 231 k = 1, lx1 !loop on the z
                     call assign_coord(i,ia,k,e,x_,y_,z_)
                     do icorner = 1, 8
                      ddx = abs(x_-x_corner(icorner))
                      ddy = abs(y_-y_corner(icorner))
                      ddz = abs(z_-z_corner(icorner))
                      if((ddx.lt.dtol.and.ddy.lt.dtol).or.
     &                (ddx.lt.dtol.and.ddz.lt.dtol).or.
     &                (ddz.lt.dtol.and.ddy.lt.dtol)) then
                         prcmpi = mpgprc(i,ia,k,e)
                         usol(k,i,myface,e) = prcrdt(3*prcmpi) 
                         vsol(k,i,myface,e) = prcrdt(3*prcmpi+1)
                         wsol(k,i,myface,e) = prcrdt(3*prcmpi+2)
                         goto 231
                      endif
                      
                     enddo
                     dwdy = 0.
                     dvdz = 0.
                     dudy = 0.
                     dvdx = 0.
                     if (kstep.gt.2) then
                     do j = b1,b2
                        dudy = dudy + dytm1(j,ia)*((4./3.)*vx(i,j,k,e)
     &                   -(1./3.)*vxlag(i,j,k,e,1))
                        dwdy = dwdy + dytm1(j,ia)*((4./3.)*vz(i,j,k,e)
     &                   -(1./3.)*vzlag(i,j,k,e,1))
            ! print *, vz(i,j,k,e), 2.0*vz(i,j,k,e)-vzlag(i,j,k,e,1)
                        
                     enddo
                  else
                     do j = b1,b2
                        dudy = dudy + dytm1(j,ia)*vx(i,j,k,e)
                        dwdy = dwdy + dytm1(j,ia)*vz(i,j,k,e)
            ! print *, vz(i,j,k,e), 2.0*vz(i,j,k,e)-vzlag(i,j,k,e,1)
                        
                     enddo
                  endif

                     do j = 1,lx1
                        prcmpi = mpgprc(j,ia,k,e)
                        dvdx = dvdx + dxm1(i,j)*prcrdt(3*prcmpi+1)
                        prcmpi = mpgprc(i,ia,j,e)
                        dvdz = dvdz + dytm1(j,k)*prcrdt(3*prcmpi+1)
                     enddo 
                     prcmpi = mpgprc(i,ia,k,e)
                     vsol(k,i,myface,e) = prcrdt(3*prcmpi+1)
                     omx = prcrdx(3*prcmpi)
                     omz = prcrdx(3*prcmpi+2)
                     usol(k,i,myface,e) = -(omz/jacmi((k-1)*lx1*lx1+
     & (ia-1)*lx1+i,e)-dvdx*rxm1(i,ia,k,e) + dudy*sym1(i,ia,k,e))/
     & (dytm1(ia,ia)*sym1(i,ia,k,e))

                     wsol(k,i,myface,e) = (omx/jacmi((k-1)*lx1*lx1+
     & (ia-1)*lx1+i,e)-dwdy*sym1(i,ia,k,e) + dvdz*tzm1(i,ia,k,e))/
     & (dytm1(ia,ia)*sym1(i,ia,k,e))


231            continue
230            continue
             else if(myface.eq.5.or.myface.eq.6) then 
                  if(myface.eq.5) then 
                     ia = 1
                     b1 = 2
                     b2 = lx1
                  else 
                     ia = lx1 
                     b1 = 1
                     b2 = lx1-1
                  endif
                  do 250 i = 1, lx1 !loop on the x
                  do 251 k = 1, lx1 !loop on the y
                     call assign_coord(i,k,ia,e,x_,y_,z_)
                     do icorner = 1, 8
                      ddx = abs(x_-x_corner(icorner))
                      ddy = abs(y_-y_corner(icorner))
                      ddz = abs(z_-z_corner(icorner))
                      if((ddx.lt.dtol.and.ddy.lt.dtol).or.
     &                (ddx.lt.dtol.and.ddz.lt.dtol).or.
     &                (ddz.lt.dtol.and.ddy.lt.dtol)) then
                         prcmpi = mpgprc(i,k,ia,e)
                         usol(k,i,myface,e) = prcrdt(3*prcmpi) 
                         vsol(k,i,myface,e) = prcrdt(3*prcmpi+1)
                         wsol(k,i,myface,e) = prcrdt(3*prcmpi+2)
                         goto 251
                      endif
                      
                     enddo
                     dwdy = 0.
                     dvdz = 0.
                     dudz = 0.
                     dwdx = 0.
                     if (kstep.gt.2) then
                     do j = b1,b2
                        dudz = dudz + dytm1(j,ia)*((4./3.)*vx(i,k,j,e)
     &                   -(1./3.)*vxlag(i,k,j,e,1))
                        dvdz = dvdz + dytm1(j,ia)*((4./3.)*vy(i,k,j,e)
     &                   -(1./3.)*vylag(i,k,j,e,1))
                     enddo
                  else 
                     do j = b1,b2
                        dudz = dudz + dytm1(j,ia)*vx(i,k,j,e)
                        dvdz = dvdz + dytm1(j,ia)*vy(i,k,j,e)
                     enddo
                  endif



                     do j = 1,lx1
                        prcmpi = mpgprc(i,j,ia,e)
                        dwdy = dwdy + dytm1(j,k)*prcrdt(3*prcmpi+2)
                        prcmpi = mpgprc(j,k,ia,e)
                        dwdx = dwdx + dxm1(i,j)*prcrdt(3*prcmpi+2)
                     enddo 
                     prcmpi = mpgprc(i,k,ia,e)
                     wsol(k,i,myface,e) = prcrdt(3*prcmpi+2)
                     omx = prcrdx(3*prcmpi)
                     omy = prcrdx(3*prcmpi+1)

                     usol(k,i,myface,e) = (omy/jacmi((ia-1)*lx1*lx1+
     & (k-1)*lx1+i,e)-dudz*tzm1(i,k,ia,e) + dwdx*rxm1(i,k,ia,e))/ 
     & (dytm1(ia,ia)*tzm1(i,k,ia,e))

                     vsol(k,i,myface,e) = -(omx/jacmi((ia-1)*lx1*lx1+
     & lx1*(k-1)+i,e)-dwdy*sym1(i,k,ia,e) + dvdz*tzm1(i,k,ia,e))/
     &  (dytm1(ia,ia)*tzm1(i,k,ia,e))


251            continue
250            continue

               endif

            enddo
         
         !trying to have correct BC's for the limit edge between x+ and y+:y- faces
         call assign_coord(1,lx1,1,e,x_,y_,z_)
         ddx = abs(x_-0.6)
         ddy = abs(y_-3.)
         if (ddy.lt.dtol.and.ddx.lt.dtol) then 
            i = lx1
            do k = 1, lx1 
               dudz = 0;
               dwdx = 0.
               do j = 2,lx1
                  dwdx = dwdx + dxm1(1,j)*wsol(k,j,3,e)
               enddo

               do j = 1,lx1
                  prcmpi = mpgprc(1,i,j,e)
                  dudz = dudz + dytm1(j,k)*prcrdt(3*prcmpi)
               enddo 
               prcmpi = mpgprc(1,i,k,e)
               omy = prcrdx(3*prcmpi+1)
               wsol(k,i,4,e) = -(omy/jacmi((k-1)*lx1*lx1
     &         +(i-1)*lx1+1,e)-dudz*tzm1(1,i,k,e) + dwdx*
     &        rxm1(1,i,k,e))/(dxm1(1,1)*rxm1(1,i,k,e))
              wsol(k,1,3,e) = wsol(k,i,4,e)
              usol(k,i,4,e) = prcrdt(3*prcmpi) 
              vsol(k,i,4,e) = prcrdt(3*prcmpi+1)
              usol(k,1,3,e) = prcrdt(3*prcmpi) 
              vsol(k,1,3,e) = prcrdt(3*prcmpi+1)

              !just try to impose a mean value... 
              usol(k,i,4,e) = 0.5*(usol(k,i-1,4,e)+usol(k,2,3,e))
              vsol(k,i,4,e) = 0.5*(vsol(k,i-1,4,e)+vsol(k,2,3,e))
              wsol(k,i,4,e) = 0.5*(wsol(k,i-1,4,e)+wsol(k,2,3,e))

              usol(k,1,3,e) = usol(k,i,4,e)
              vsol(k,1,3,e) = vsol(k,i,4,e)
              wsol(k,1,3,e) = wsol(k,i,4,e)

              
            enddo
         endif



         call assign_coord(1,1,1,e,x_,y_,z_)
         ddx = abs(x_-0.6)
         ddy = abs(y_-1.)
         if (ddy.lt.dtol.and.ddx.lt.dtol) then 
            i = 1
            do k = 1, lx1 

              !just try to impose a mean value... 
              usol(k,i,4,e) = 0.5*(usol(k,i+1,4,e)+usol(k,2,1,e))
              vsol(k,i,4,e) = 0.5*(vsol(k,i+1,4,e)+vsol(k,2,1,e))
              wsol(k,i,4,e) = 0.5*(wsol(k,i+1,4,e)+wsol(k,2,1,e))

              usol(k,1,1,e) = usol(k,i,4,e)
              vsol(k,1,1,e) = vsol(k,i,4,e)
              wsol(k,1,1,e) = wsol(k,i,4,e)

            enddo
         endif

         !doing it also for the 2 lateral edges...
         call assign_coord(1,1,1,e,x_,y_,z_)
         ddx = abs(x_-0.6)
         ddz = abs(z_-0.625)
         if (ddz.lt.dtol.and.ddx.lt.dtol) then 
            k = 1
            do i = 1, lx1 

              !just try to impose a mean value... 
              usol(k,i,4,e) = 0.5*(usol(k+1,i,4,e)+usol(i,2,5,e))
              vsol(k,i,4,e) = 0.5*(vsol(k+1,i,4,e)+vsol(i,2,5,e))
              wsol(k,i,4,e) = 0.5*(wsol(k+1,i,4,e)+wsol(i,2,5,e))

              usol(i,1,5,e) = usol(k,i,4,e)
              vsol(i,1,5,e) = vsol(k,i,4,e)
              wsol(i,1,5,e) = wsol(k,i,4,e)

            enddo
         endif
         call assign_coord(1,1,lx1,e,x_,y_,z_)
         ddx = abs(x_-0.6)
         ddz = abs(z_-1.375)
         if (ddz.lt.dtol.and.ddx.lt.dtol) then 
            k = lx1
            do i = 1, lx1 

              !just try to impose a mean value... 
              usol(k,i,4,e) = 0.5*(usol(k-1,i,4,e)+usol(i,2,6,e))
              vsol(k,i,4,e) = 0.5*(vsol(k-1,i,4,e)+vsol(i,2,6,e))
              wsol(k,i,4,e) = 0.5*(wsol(k-1,i,4,e)+wsol(i,2,6,e))

              usol(i,1,6,e) = usol(k,i,4,e)
              vsol(i,1,6,e) = vsol(k,i,4,e)
              wsol(i,1,6,e) = wsol(k,i,4,e)

            enddo
         endif

         ! and doing it also on the other side of the cube
         call assign_coord(lx1,lx1,1,e,x_,y_,z_)
         ddx = abs(x_-2.6)
         ddy = abs(y_-3.)
         if (ddy.lt.dtol.and.ddx.lt.dtol) then 
            i = lx1
            do k = 1, lx1 
              !just try to impose a mean value... 
              usol(k,i,2,e) = 0.5*(usol(k,i-1,2,e)+usol(k,lx1-1,3,e))
              vsol(k,i,2,e) = 0.5*(vsol(k,i-1,2,e)+vsol(k,lx1-1,3,e))
              wsol(k,i,2,e) = 0.5*(wsol(k,i-1,2,e)+wsol(k,lx1-1,3,e))

              usol(k,lx1,3,e) = usol(k,i,2,e)
              vsol(k,lx1,3,e) = vsol(k,i,2,e)
              wsol(k,lx1,3,e) = wsol(k,i,2,e)

              
            enddo
         endif



         call assign_coord(lx1,1,1,e,x_,y_,z_)
         ddx = abs(x_-2.6)
         ddy = abs(y_-1.)
         if (ddy.lt.dtol.and.ddx.lt.dtol) then 
            i = 1
            do k = 1, lx1 

              !just try to impose a mean value... 
              usol(k,i,2,e) = 0.5*(usol(k,i+1,2,e)+usol(k,lx1-1,1,e))
              vsol(k,i,2,e) = 0.5*(vsol(k,i+1,2,e)+vsol(k,lx1-1,1,e))
              wsol(k,i,2,e) = 0.5*(wsol(k,i+1,2,e)+wsol(k,lx1-1,1,e))

              usol(k,lx1,1,e) = usol(k,i,2,e)
              vsol(k,lx1,1,e) = vsol(k,i,2,e)
              wsol(k,lx1,1,e) = wsol(k,i,2,e)

            enddo
         endif

         !doing it also for the 2 lateral edges...
         call assign_coord(lx1,1,1,e,x_,y_,z_)
         ddx = abs(x_-2.6)
         ddz = abs(z_-0.625)
         if (ddz.lt.dtol.and.ddx.lt.dtol) then 
            k = 1
            do i = 1, lx1 

              !just try to impose a mean value... 
              usol(k,i,2,e) = 0.5*(usol(k+1,i,2,e)+usol(i,lx1-1,5,e))
              vsol(k,i,2,e) = 0.5*(vsol(k+1,i,2,e)+vsol(i,lx1-1,5,e))
              wsol(k,i,2,e) = 0.5*(wsol(k+1,i,2,e)+wsol(i,lx1-1,5,e))

              usol(i,lx1,5,e) = usol(k,i,2,e)
              vsol(i,lx1,5,e) = vsol(k,i,2,e)
              wsol(i,lx1,5,e) = wsol(k,i,2,e)

            enddo
         endif
         call assign_coord(lx1,1,lx1,e,x_,y_,z_)
         ddx = abs(x_-2.6)
         ddz = abs(z_-1.375)
         if (ddz.lt.dtol.and.ddx.lt.dtol) then 
            k = lx1
            do i = 1, lx1 

              !just try to impose a mean value... 
              usol(k,i,2,e) = 0.5*(usol(k-1,i,2,e)+usol(i,lx1-1,6,e))
              vsol(k,i,2,e) = 0.5*(vsol(k-1,i,2,e)+vsol(i,lx1-1,6,e))
              wsol(k,i,2,e) = 0.5*(wsol(k-1,i,2,e)+wsol(i,lx1-1,6,e))

              usol(i,lx1,6,e) = usol(k,i,2,e)
              vsol(i,lx1,6,e) = vsol(k,i,2,e)
              wsol(i,lx1,6,e) = wsol(k,i,2,e)

            enddo
         endif
         

         call assign_coord(1,lx1,lx1,e,x_,y_,z_)
         ddy = abs(y_-3.)
         ddz = abs(z_-0.625)
         if (ddz.lt.dtol.and.ddy.lt.dtol) then 
            k = lx1
            do i = 1, lx1 

              !just try to impose a mean value... 
              usol(k,i,6,e) = 0.5*(usol(k-1,i,6,e)+usol(lx1-1,i,3,e))
              vsol(k,i,6,e) = 0.5*(vsol(k-1,i,6,e)+vsol(lx1-1,i,3,e))
              wsol(k,i,6,e) = 0.5*(wsol(k-1,i,6,e)+wsol(lx1-1,i,3,e))

              usol(k,i,3,e) = usol(k,i,6,e)
              vsol(k,i,3,e) = vsol(k,i,6,e)
              wsol(k,i,3,e) = wsol(k,i,6,e)

            enddo
         endif

         call assign_coord(1,1,lx1,e,x_,y_,z_)
         ddy = abs(y_-1.)
         ddz = abs(z_-1.375)
         if (ddz.lt.dtol.and.ddy.lt.dtol) then 
            k = 1
            do i = 1, lx1 

              !just try to impose a mean value... 
              usol(k,i,6,e) = 0.5*(usol(k+1,i,6,e)+usol(lx1-1,i,1,e))
              vsol(k,i,6,e) = 0.5*(vsol(k+1,i,6,e)+vsol(lx1-1,i,1,e))
              wsol(k,i,6,e) = 0.5*(wsol(k+1,i,6,e)+wsol(lx1-1,i,1,e))

              usol(lx1,i,1,e) = usol(k,i,6,e)
              vsol(lx1,i,1,e) = vsol(k,i,6,e)
              wsol(lx1,i,1,e) = wsol(k,i,6,e)

            enddo
         endif

         call assign_coord(1,lx1,1,e,x_,y_,z_)
         ddy = abs(y_-3.)
         ddz = abs(z_-0.625)
         if (ddz.lt.dtol.and.ddy.lt.dtol) then 
            k = lx1
            do i = 1, lx1 

              !just try to impose a mean value... 
              usol(k,i,5,e) = 0.5*(usol(k-1,i,5,e)+usol(2,i,3,e))
              vsol(k,i,5,e) = 0.5*(vsol(k-1,i,5,e)+vsol(2,i,3,e))
              wsol(k,i,5,e) = 0.5*(wsol(k-1,i,5,e)+wsol(2,i,3,e))

              usol(1,i,3,e) = usol(k,i,5,e)
              vsol(1,i,3,e) = vsol(k,i,5,e)
              wsol(1,i,3,e) = wsol(k,i,5,e)

            enddo
         endif

         call assign_coord(1,1,1,e,x_,y_,z_)
         ddy = abs(y_-1.)
         ddz = abs(z_-0.625)
         if (ddz.lt.dtol.and.ddy.lt.dtol) then 
            k = 1
            do i = 1, lx1 

              !just try to impose a mean value... 
              usol(k,i,5,e) = 0.5*(usol(k+1,i,5,e)+usol(2,i,1,e))
              vsol(k,i,5,e) = 0.5*(vsol(k+1,i,5,e)+vsol(2,i,1,e))
              wsol(k,i,5,e) = 0.5*(wsol(k+1,i,5,e)+wsol(2,i,1,e))

              usol(1,i,1,e) = usol(k,i,5,e)
              vsol(1,i,1,e) = vsol(k,i,5,e)
              wsol(1,i,1,e) = wsol(k,i,5,e)

            enddo
         endif

   !       ! do the corner! yes yes! its crazy
   !       call assign_coord(1,lx1,lx1,e,x_,y_,z_)
   !       ddy = abs(y_-3.)
   !       ddz = abs(z_-0.875)
   !       ddx = abs(x_-1.0)
   !       if (ddz.lt.dtol.and.ddy.lt.dtol.and.ddx.lt.dtol) then 
   !          aa = - dytm1(lx1,lx1)*tzm1(1,lx1,lx1,e)*1.0001
   !          bb = dytm1(lx1,lx1)*sym1(1,lx1,lx1,e)
   !          cc = dytm1(lx1,lx1)*tzm1(1,lx1,lx1,e)
   !          dd = - dxm1(1,1)*rxm1(1,lx1,lx1,e)
   !          ee = - dytm1(lx1,lx1)*sym1(1,lx1,lx1,e)
   !          ff = dxm1(1,1)*rxm1(1,lx1,lx1,e)

   !          dwdy = 0.
   !          dvdz = 0.
   !          dudz = 0.
   !          dwdx = 0.
   !          dvdx = 0.
   !          dudy = 0.

   !          do j = 1, lx1-1 
   !    dwdy  = dwdy + sym1(1,lx1,lx1,e)*dytm1(j,lx1)*wsol(lx1,j,4,e)
   !    dudy  = dudy + sym1(1,lx1,lx1,e)*dytm1(j,lx1)*usol(lx1,j,4,e)

   !    dvdz  = dvdz + tzm1(1,lx1,lx1,e)*dytm1(j,lx1)*vsol(j,lx1,4,e)
   !    dudz  = dudz + tzm1(1,lx1,lx1,e)*dytm1(j,lx1)*usol(j,lx1,4,e)

   !          enddo 
   !          do j = 2, lx1
   !    dwdx  = dwdx + rxm1(1,lx1,lx1,e)*dxm1(1,j)*wsol(lx1,j,6,e)
   !    dvdx  = dvdx + rxm1(1,lx1,lx1,e)*dxm1(1,j)*vsol(lx1,j,6,e)

   !          enddo 
   !          prcmpi = mpgprc(1,lx1,lx1,e)
   !          omx = prcrdx(3*prcmpi)
   !          omy = prcrdx(3*prcmpi+1)
   !          omz = prcrdx(3*prcmpi+2)

   !    rr= omx/jacmi((lx1-1)*lx1*lx1+(lx1-1)*lx1+1,e)-dwdy+dvdz
   !    ss= omy/jacmi((lx1-1)*lx1*lx1+(lx1-1)*lx1+1,e)-dudz+dwdx
   !    tt= omz/jacmi((lx1-1)*lx1*lx1+(lx1-1)*lx1+1,e)-dvdx+dudy 

   !    vsol(lx1,lx1,4,e) = (cc*bb*tt+dd*ee*rr-ss*ee*bb)
   !   & /(dd*ee*aa+bb*cc*ff)
   !    usol(lx1,lx1,4,e) = (tt-ff*vsol(lx1,lx1,4,e))/ee 
   !    wsol(lx1,lx1,4,e) = (rr-aa*vsol(lx1,lx1,4,e))/bb 

   !    ! just try pure dirichlet...
   !    usol(lx1,lx1,4,e) = prcrdt(3*prcmpi)
   !    vsol(lx1,lx1,4,e) = prcrdt(3*prcmpi+1)
   !    wsol(lx1,lx1,4,e) = prcrdt(3*prcmpi+2)

   !    usol(lx1,1,6,e) = usol(lx1,lx1,4,e)
   !    usol(lx1,1,3,e) = usol(lx1,lx1,4,e)

   !    vsol(lx1,1,6,e) = vsol(lx1,lx1,4,e)
   !    vsol(lx1,1,3,e) = vsol(lx1,lx1,4,e)

   !    wsol(lx1,1,6,e) = wsol(lx1,lx1,4,e)
   !    wsol(lx1,1,3,e) = wsol(lx1,lx1,4,e)



   !    ! we can also check that the vorticity is correclty imposed: 
   !    dwdy = 0.
   !    dvdz = 0.
   !    dudz = 0.
   !    dwdx = 0.
   !    dvdx = 0.
   !    dudy = 0.

   !          do j = 1, lx1 
   !    dwdy  = dwdy + sym1(1,lx1,lx1,e)*dytm1(j,lx1)*wsol(lx1,j,4,e)
   !    dudy  = dudy + sym1(1,lx1,lx1,e)*dytm1(j,lx1)*usol(lx1,j,4,e)

   !    dvdz  = dvdz + tzm1(1,lx1,lx1,e)*dytm1(j,lx1)*vsol(j,lx1,4,e)
   !    dudz  = dudz + tzm1(1,lx1,lx1,e)*dytm1(j,lx1)*usol(j,lx1,4,e)

   !          enddo 
   !          do j = 1, lx1
   !    dwdx  = dwdx + rxm1(1,lx1,lx1,e)*dxm1(1,j)*wsol(lx1,j,6,e)
   !    dvdx  = dvdx + rxm1(1,lx1,lx1,e)*dxm1(1,j)*vsol(lx1,j,6,e)
   !          enddo

   !    print *, 'om_x:', omx, (dwdy-dvdz)*jacmi((lx1-1)
   !   & *lx1*lx1+(lx1-1)*lx1+1,e) 
   !    print *, 'om_y:', omy, (dudz-dwdx)*jacmi((lx1-1)
   !   & *lx1*lx1+(lx1-1)*lx1+1,e)
   !    print *, 'om_z:', omz, (dvdx-dudy)*jacmi((lx1-1)
   !   & *lx1*lx1+(lx1-1)*lx1+1,e)
   !    print *, usol(lx1,lx1,4,e)
   !    print *, vsol(lx1,lx1,4,e)
   !    print *, wsol(lx1,lx1,4,e)

   !    !! check the system:
   !    print *, aa,bb,cc,dd,ee,ff 
   !    print *, rr, ss, tt



            

            


   !          usol(lx1,lx1,4,e) = (1./3.)*(usol(lx1-1,lx1,4,e)+ 
   !   &       usol(lx1,lx1-1,4,e) + usol(lx1-1,1,6,e))

   !          vsol(lx1,lx1,4,e) = (1./3.)*(vsol(lx1-1,lx1,4,e)+ 
   !   &       vsol(lx1,lx1-1,4,e) + vsol(lx1-1,1,6,e))

   !          wsol(lx1,lx1,4,e) = (1./3.)*(wsol(lx1-1,lx1,4,e)+ 
   !   &       wsol(lx1,lx1-1,4,e) + wsol(lx1-1,1,6,e))

   !          usol(lx1,1,6,e) = usol(lx1,lx1,4,e)
   !          usol(lx1,1,3,e) = usol(lx1,lx1,4,e)

   !          vsol(lx1,1,6,e) = vsol(lx1,lx1,4,e)
   !          vsol(lx1,1,3,e) = vsol(lx1,lx1,4,e)

   !          wsol(lx1,1,6,e) = wsol(lx1,lx1,4,e)
   !          wsol(lx1,1,3,e) = wsol(lx1,lx1,4,e)
           

         ! endif

         
      enddo

      return 
      end


       