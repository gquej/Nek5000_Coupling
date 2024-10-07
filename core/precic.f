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
c------------------------------------------------------------------------------------------------------------------------

      subroutine interpolate_u(omshdi, prcvr2, prcwdt)

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
        real prcvr2(0:omshdi*3-1)
        real prcwdt(0:omshdi*3-1)
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

        call fgslib_findpts_setup(inth_hpts, nekcomm, np, 3,
     &  xm1, ym1, zm1, lx1, ly1, lz1, nelt, nxf, nyf, nzf, bb_t,
     &  n, n, npt_max, tol)



      !we want to interpolate the Nek velocity onto a staggered VPM grid. The coordinates for the 
      !3 components of velocity are (with xi,yj,zk the Murphy vertex where the vorticity is stored):
      !for u: (xi, yj + h/2, zk + h/2)
      !for v: (xi+h/2, yj  , zk + h/2)
      !for u: (xi+ h/2, yj + h/2, zk )

      hvpm = 1./24.
      !modify the vertex coordinates to interpolate the x-component: adding h/2 to y and z
      do i = 0, omshdi-1

         vmodif(3*i)   = prcvr2(3*i)
         vmodif(3*i+1) = prcvr2(3*i+1) + 0.5*hvpm
         vmodif(3*i+2) = prcvr2(3*i+2) + 0.5*hvpm
      enddo

        call fgslib_findpts(inth_hpts, rcode, 1, 
     &                    proc, 1,
     &                    elid, 1,
     &                    rst, 3,
     &                    dist, 1,
     &                    vmodif(0), 3,
     &                    vmodif(1), 3,
     &                    vmodif(2), 3, omshdi)
        
        call fgslib_findpts_eval(inth_hpts,prcwdt(0), 3,
     &                         rcode, 1, 
     &                         proc, 1, 
     &                         elid, 1, 
     &                         rst, 3, omshdi, 
     &                         vx)
       !modify the vertex coordinates to interpolate the y-component: adding h/2 to x and z
      do i = 0, omshdi-1
         vmodif(3*i)   = prcvr2(3*i) + 0.5*hvpm
         vmodif(3*i+1) = prcvr2(3*i+1) 
         vmodif(3*i+2) = prcvr2(3*i+2) + 0.5*hvpm
      enddo

        call fgslib_findpts(inth_hpts, rcode, 1, 
     &                    proc, 1,
     &                    elid, 1,
     &                    rst, 3,
     &                    dist, 1,
     &                    vmodif(0), 3,
     &                    vmodif(1), 3,
     &                    vmodif(2), 3, omshdi)

        call fgslib_findpts_eval(inth_hpts,prcwdt(1), 3,
     &                         rcode, 1, 
     &                         proc, 1, 
     &                         elid, 1, 
     &                         rst, 3, omshdi, 
     &                         vy)

       !modify the vertex coordinates to interpolate the z-component: adding h/2 to x and y
      do i = 0, omshdi-1
         vmodif(3*i)   = prcvr2(3*i) + 0.5*hvpm
         vmodif(3*i+1) = prcvr2(3*i+1) + 0.5*hvpm
         vmodif(3*i+2) = prcvr2(3*i+2) 
      enddo

        call fgslib_findpts(inth_hpts, rcode, 1, 
     &                    proc, 1,
     &                    elid, 1,
     &                    rst, 3,
     &                    dist, 1,
     &                    vmodif(0), 3,
     &                    vmodif(1), 3,
     &                    vmodif(2), 3, omshdi)
        call fgslib_findpts_eval(inth_hpts,prcwdt(2), 3,
     &                         rcode, 1, 
     &                         proc, 1, 
     &                         elid, 1, 
     &                         rst, 3, omshdi, 
     &                         vz)
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
        print *, tmpl1(1,e),tmpl1(2,e),tmpl1(3,e),tmpl1(4,e),
     &   tmpl1(5,e),tmpl1(6,e)
        print *, tmpl4(1,e),tmpl4(2,e),tmpl4(3,e),tmpl4(4,e),
     &   tmpl4(5,e),tmpl4(6,e)
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

      x_corner(1) = -0.5
      x_corner(2) = 0.5
      x_corner(3) = 0.5
      x_corner(4) = -0.5
      x_corner(5) = -0.5
      x_corner(6) = 0.5
      x_corner(7) = 0.5
      x_corner(8) = -0.5
      y_corner(1) = -0.5
      y_corner(2) = -0.5
      y_corner(3) = 0.5
      y_corner(4) = 0.5
      y_corner(5) = -0.5
      y_corner(6) = -0.5
      y_corner(7) = 0.5
      y_corner(8) = 0.5
      z_corner(1) = -0.5
      z_corner(2) = -0.5
      z_corner(3) = -0.5
      z_corner(4) = -0.5
      z_corner(5) = 0.5
      z_corner(6) = 0.5
      z_corner(7) = 0.5
      z_corner(8) = 0.5

      return 
      end
      
c----------------------------------------------------------------------------------------
      subroutine compute_bc
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

                     do j = b1,b2
                        dvdx = dvdx + dxm1(ia,j)*vy(j,i,k,e)
                        dwdx = dwdx + dxm1(ia,j)*vz(j,i,k,e)
                     enddo

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

                     do j = b1,b2
                        dudy = dudy + dytm1(j,ia)*vx(i,j,k,e)
                        dwdy = dwdy + dytm1(j,ia)*vz(i,j,k,e)
                     enddo
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

                     do j = b1,b2
                        dudz = dudz + dytm1(j,ia)*vx(i,k,j,e)
                        dvdz = dvdz + dytm1(j,ia)*vy(i,k,j,e)
                     enddo

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

         enddo

      return 
      end


       