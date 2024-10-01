c-----------------------------------------------------------------------
      subroutine nek_init(comm)
c
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
c
      include 'OPCTR'
      include 'CTIMER'

C     used scratch arrays
C     NOTE: no initial declaration needed. Linker will take 
c           care about the size of the CBs automatically
c
c      COMMON /CTMP1/ DUMMY1(LCTMP1)
c      COMMON /CTMP0/ DUMMY0(LCTMP0)
c
c      COMMON /SCRNS/ DUMMY2(LX1,LY1,LZ1,LELT,7)
c      COMMON /SCRUZ/ DUMMY3(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCREV/ DUMMY4(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRVH/ DUMMY5(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRMG/ DUMMY6(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCRCH/ DUMMY7(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRSF/ DUMMY8(LX1,LY1,LZ1,LELT,3)
c      COMMON /SCRCG/ DUMM10(LX1,LY1,LZ1,LELT,1)

      integer comm
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
  
      common /rdump/ ntdump

      real kwave2
      logical ifemati

      real rtest
      integer itest
      integer*8 itest8
      character ctest
      logical ltest 

      common /c_is1/ glo_num(lx1 * ly1 * lz1, lelt)
      common /ivrtx/ vertex((2 ** ldim) * lelt)
      integer*8 glo_num, ngv
      integer*8 vertex

      ! set word size for REAL
      wdsize = sizeof(rtest)
      ! set word size for INTEGER
      isize = sizeof(itest)
      ! set word size for INTEGER*8
      isize8 = sizeof(itest8) 
      ! set word size for LOGICAL
      lsize = sizeof(ltest) 
      ! set word size for CHARACTER
      csize = sizeof(ctest)

      call setupcomm(comm,newcomm,newcommg,'','')
      intracomm   = newcomm   ! within a session
      nekcomm     = newcomm
      iglobalcomm = newcommg  ! across all sessions
      call iniproc()

      if (nid.eq.nio) call printHeader

      etimes = dnekclock()
      istep  = 0

      call opcount(1)

      call initdim         ! Initialize / set default values.
      call initdat
      call files

      call readat          ! Read .rea +map file

      if (nio.eq.0) then
         write(6,12) 'nelgt/nelgv/lelt:',nelgt,nelgv,lelt
         write(6,12) 'lx1/lx2/lx3/lxd: ',lx1,lx2,lx3,lxd
 12      format(1X,A,4I12)
         write(6,*)
      endif

      call setvar          ! Initialize most variables

      instep=1             ! Check for zero steps
      if (nsteps.eq.0 .and. fintim.eq.0.) instep=0

      igeom = 2
      call setup_topo      ! Setup domain topology  

      call genwz           ! Compute GLL points, weights, etc.

      if(nio.eq.0) write(6,*) 'call usrdat'
      call usrdat
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat' 

      call gengeom(igeom)  ! Generate geometry, after usrdat 

      if (ifmvbd) call setup_mesh_dssum ! Set mesh dssum (needs geom)

      if(nio.eq.0) write(6,*) 'call usrdat2'
      call usrdat2
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat2' 

      call count_bdry   ! count the number of faces with assigned BCs
      call fix_geom

      call vrdsmsh          ! verify mesh topology
      call mesh_metrics     ! print some metrics

      call setlog(.true.)   ! Initalize logical flags

      if (ifneknekc) call neknek_setup

      call bcmask  ! Set BC masks for Dirichlet boundaries.

      if (fintim.ne.0.0 .or. nsteps.ne.0) 
     $   call geneig(igeom) ! eigvals for tolerances

      call dg_setup ! Setup DG, if dg flag is set.

      if (ifflow.and.iftran) then ! Init pressure solver 
         if (fintim.ne.0 .or. nsteps.ne.0) call prinit
      endif

      if(ifcvode) call cv_setsize

      if(nio.eq.0) write(6,*) 'call usrdat3'
      call usrdat3
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat3'

      call setics
      call setprop

      if (instep.ne.0) then
         if (ifneknekc) call neknek_exchange
         if (ifneknekc) call chk_outflow

         if (nio.eq.0) write(6,*) 'call userchk'
         call userchk
         if(nio.eq.0) write(6,'(A,/)') ' done :: userchk' 
      endif

      call setprop      ! call again because input has changed in userchk

      if (ifcvode .and. nsteps.gt.0) call cv_init

      call comment
      call sstest (isss) 

      call dofcnt

      jp = 0  ! Set perturbation field count to 0 for baseline flow
      p0thn = p0th

      call in_situ_init()

      call time00       !     Initalize timers to ZERO
      call opcount(2)

      ntdump=0
      if (timeio.ne.0.0) ntdump = int( time/timeio )

      tinit = dnekclock_sync() - etimes
      if (nio.eq.0) then
        write (6,*) ' '
        if (time.ne.0.0) write (6,'(a,e14.7)') ' Initial time:',time
        write (6,'(a,g13.5,a)') 
     &     ' Initialization successfully completed ', tinit, ' sec'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_solve

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'GEOM'
      include 'DXYZ'
      integer e,i,j,ifc
      common /nekcb/ cbcmesh
      character cb*3

      real nx,ny
      real mat_sn(2*lx1,2*lx1,lelt)
      real uart, vart, omart
      real temp1, temp2

      real u_,v_,p_,om_,du_
      real x_,y_

      !variables for the new implementation
      integer tmpl1(2*ldim,nelv), tmpl2(2*ldim,nelv)
      integer tmpl3(2*ldim,nelv), tmpl4(2*ldim,nelv)
      integer n_fbc, n_fnbc
      integer l_fbc(nelv), l_fnbc(nelv)

      call nekgsync()

      if (instep.eq.0) then
        if(nid.eq.0) write(6,'(/,A,/,A,/)') 
     &     ' nsteps=0 -> skip time loop',
     &     ' running solver in post processing mode'
      else
        if(nio.eq.0) write(6,'(/,A,/)') 'Starting time loop ...'
      endif

      isyc  = 0
      if(ifsync) isyc=1
      itime = 0
#ifdef TIMER
      itime = 1
#endif

      ! start measurements
      dtmp = dnekgflops()

      istep  = 0
      msteps = 1

      irstat = int(param(120))

      !assemble the matrix to impose vorticity BC
      nx = -1.0
      ny = 0.

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
          tmpl3(0,e) = 0
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
        print *, tmpl1(1,e),tmpl1(2,e),tmpl1(3,e),tmpl1(4,e)
        print *, tmpl4(1,e),tmpl4(2,e),tmpl4(3,e),tmpl4(4,e)
      enddo

      ! now, we build the matrix, side by side
      ! do e = 1, nelv 
      !    ! first, init the matrix to 0
      !    do i = 1, 2*lx1*2*ldim 
      !       do j = 1, 2*lx1*2*ldim 
      !          mat_sn(i,j,e) = 0.
      !       enddo 
      !    enddo 

      !    do ifc = 1, l_fbc(e)
      !       if(tmpl1(ifc,e).eq.4) then
      !          !first point: might impact variables of side 1 if they are unknowns
      !          if(tmpl2(tmpl4(1,e),e).eq.1) then 
      !             ehexto8hex
      !          else 
      !             ehexto8hex
      !          endif
      !          !then do loop on point 2, lx1-1
      !          do i = 2, lx1-1
      !             e 
      !          enddo
      !          !last point: might impact variables of side 3 if they are unknowns
      !          last point
      !       else if () then 

      !       else if () then 

      !       else if () then 

      !       endif


      !    enddo

      ! enddo




      do e = 1,nelt 
         do i = 1,lx1
            do j = 1, lx1 
               mat_sn(i,j,e) = -sym1(1,i,1,e)*dytm1(j,i)
               mat_sn(i,j+lx1,e) = sxm1(1,i,1,e)*dytm1(j,i)
            enddo
            mat_sn(i,i,e)    = -rym1(1,i,1,e)*dxm1(1,1) - 
     &       sym1(1,i,1,e)*dytm1(i,i)
            mat_sn(i,i+lx1,e) = rxm1(1,i,1,e)*dxm1(1,1) +
     &        sxm1(1,i,1,e)*dytm1(i,i)
            mat_sn(i+lx1,i,e) = nx
            mat_sn(i+lx1,i+lx1,e) = ny
         enddo

         call LUdec(mat_sn(1,1,e),2*lx1, persnl(1,e),1e-12)
         ! call invmat(mat_sn(1,1,e),2*lx1, matsnl(1,1,e), persnl(1,e))
      enddo
      !creating an uniform rhs with u = 1, v = 1
      uart = 1.
      vart = 1.8792
      omart = 0.

      do e = 1,nelt
         do i = 1,lx1
            temp1 = 0.
            temp2 = 0.
            do j = 2,lx1
               x_ = xm1(j,i,1,e)
               y_ = ym1(j,i,1,e)
               call GVLFD(x_,y_,u_,v_, p_, om_, du_, du_, du_, du_)
               temp1 = temp1 + dxm1(1,j)* v_
               temp2 = temp2 + dxm1(1,j)* u_
            enddo
            x_ = xm1(1,i,1,e)
            y_ = ym1(1,i,1,e)
            call GVLFD(x_,y_,u_,v_, p_, om_, du_, du_, du_, du_)
            rhssnl(i,e) = om_/jacmi((i-1)*lx1+1,e) - 
     &       rxm1(1,i,1,e)*temp1 + rym1(1,i,1,e)*temp2
            rhssnl(i+lx1,e) = -u_
         enddo
      enddo 

      do e = 1,nelt 
         ! call solaxb(matsnl(1,1,e),2*lx1, rhssnl(1,e),xsol, persnl(1,e))
         call LUsolv(mat_sn(1,1,e),2*lx1, rhssnl(1,e), persnl(1,e),
     &    xsol(1,e))
      enddo

      do e = 1, nelv
        if(abs(xm1(1,1,1,e)+0.5).lt.1e-12) then 
        print *, 'Element ', e, xm1(1,1,1,e), ym1(1,1,1,e)
        do i = 1, LX1
            x_ = xm1(1,i,1,e)
            y_ = ym1(1,i,1,e)
            call GVLFD(x_,y_,u_,v_, p_, om_, du_, du_, du_, du_)
            print *, abs(xsol(i,e)-u_), abs(xsol(i+lx1,e)-v_)
            ! print *, xsol(i,e),u_, xsol(i+lx1,e), v_
        enddo
      endif
      enddo

      print *, 'finding the boundaries'

      do e = 1,nelv 
         do ifc = 1,4
            cb = cbc(ifc,e,1)
            print *, 'Element ', e, 'Face ', ifc, 'BC ', cb
         enddo
      enddo







      do kstep=1,nsteps,msteps
         call nek__multi_advance(kstep,msteps)
         if(kstep.ge.nsteps) lastep = 1
         call check_ioinfo  
         call set_outfld
         etime1 = dnekclock()
         call userchk
         tuchk = tuchk + dnekclock()-etime1
         call prepost (ifoutfld,'his')
         call in_situ_check()
         if (mod(kstep,irstat).eq.0 .and. lastep.eq.0) call runstat 
         if (lastep .eq. 1) goto 1001
      enddo
 1001 lastep=1

      call comment

c     check for post-processing mode
      if (instep.eq.0) then
         nsteps=0
         istep=0
         if(nio.eq.0) write(6,*) 'call userchk'
         call userchk
         if(nio.eq.0) write(6,*) 'done :: userchk'
         call prepost (.true.,'his')
      else
         if (nio.eq.0) write(6,'(/,A,/)') 
     $      'end of time-step loop' 
      endif


      RETURN
      END

c-----------------------------------------------------------------------
      subroutine nek_advance

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      common /cgeom/ igeom

      ntot = lx1*ly1*lz1*nelv

      call nekgsync

      call setup_convect(2) ! Save conv vel
      print *, "before settime", time
      if (iftran) call settime
      print *, "after settime", time
      if (ifmhd ) call cfl_check
      call setsolv
      call comment

      if (ifsplit) then   ! PN/PN formulation
         do igeom=1,ngeom

         if (ifneknekc .and. igeom.gt.2) then
            if (ifneknekm.and.igeom.eq.3) call neknek_setup
            call neknek_exchange
         endif

         ! call here before we overwrite wx 
         if (ifheat .and. ifcvode) call heat_cvode (igeom)   

         if (ifgeom) then
            call gengeom (igeom)
            call geneig  (igeom)
         endif

         if (ifheat) call heat (igeom)

         if (igeom.eq.2) then  
            call setprop
            call rzero(qtl,ntot)
            if (iflomach) call qthermal
         endif

         if (ifflow)          call fluid    (igeom)
         if (ifmvbd)          call meshv    (igeom)
         if (igeom.eq.ngeom.and.filterType.eq.1)
     $                        call q_filter(param(103))

         enddo

      else                ! PN-2/PN-2 formulation
         call setprop
         do igeom=1,ngeom

            if (ifneknekc .and. igeom.gt.2) then
              if (ifneknekm.and.igeom.eq.3) call neknek_setup
              call neknek_exchange
            endif

            ! call here before we overwrite wx 
            if (ifheat .and. ifcvode) call heat_cvode (igeom)   

            if (ifgeom) then
               if (.not.ifrich) call gengeom (igeom)
               call geneig  (igeom)
            endif

            if (ifmhd) then
               if (ifheat)      call heat     (igeom)
                                call induct   (igeom)
            elseif (ifpert) then
               if (ifbase.and.ifheat)  call heat          (igeom)
               if (ifbase.and.ifflow)  call fluid         (igeom)
               if (ifflow)             call fluidp        (igeom)
               if (ifheat)             call heatp         (igeom)
            else  ! std. nek case
               if (ifheat)             call heat          (igeom)
               if (ifflow)             call fluid         (igeom)
               if (ifmvbd)             call meshv         (igeom)
            endif
            if (igeom.eq.ngeom.and.filterType.eq.1)
     $         call q_filter(param(103))
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine nek_end

      include 'SIZE'
      include 'TOTAL'
      include 'DPROCMAP'

      if(instep.ne.0) call runstat

c      if (ifstrs) then
c         call fgslib_crs_free(xxth_strs) 
c      else
c         call fgslib_crs_free(xxth(1))
c      endif

#ifdef DPROCMAP
#ifdef MPI
      call MPI_Win_free(dProcmapH, ierr)
#endif
#endif 
      call in_situ_end()
      call exitt0()

      return
      end
c-----------------------------------------------------------------------
      subroutine nek__multi_advance(kstep,msteps)

      include 'SIZE'
      include 'TOTAL'

      do i=1,msteps
         istep = istep+i
         call nek_advance

         if (ifneknekc) then 
            call neknek_exchange
            call bcopy
            call chk_outflow
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
