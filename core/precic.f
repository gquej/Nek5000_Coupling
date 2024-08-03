      subroutine prc_omesh
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
           if (CB.eq.'v  ') then !if the face is a boundary
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