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
        integer rcode(0:omshdi-1), proc(0:omshdi-1), elid (0:omshdi-1)
        real rst(0:omshdi*3-1)
        real dist(0:omshdi-1)
        real prcvr2(0:omshdi*3-1)
        real prcwdt(0:omshdi*3-1)


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

        call fgslib_findpts(inth_hpts, rcode, 1, 
     &                    proc, 1,
     &                    elid, 1,
     &                    rst, 3,
     &                    dist, 1,
     &                    prcvr2(0), 3,
     &                    prcvr2(1), 3,
     &                    prcvr2(2), 3, omshdi)
        
        call fgslib_findpts_eval(inth_hpts,prcwdt(0), 3,
     &                         rcode, 1, 
     &                         proc, 1, 
     &                         elid, 1, 
     &                         rst, 3, omshdi, 
     &                         vx)
        call fgslib_findpts_eval(inth_hpts,prcwdt(1), 3,
     &                         rcode, 1, 
     &                         proc, 1, 
     &                         elid, 1, 
     &                         rst, 3, omshdi, 
     &                         vy)
        call fgslib_findpts_eval(inth_hpts,prcwdt(2), 3,
     &                         rcode, 1, 
     &                         proc, 1, 
     &                         elid, 1, 
     &                         rst, 3, omshdi, 
     &                         vz)
      return
      end 


       