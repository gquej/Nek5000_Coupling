c-----------------------------------------------------------------------
      subroutine nek_init(comm)
c
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
c
      include 'OPCTR'
      include 'CTIMER'
      include 'PRECIC'

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
      integer nel
      real bbbx1, bbbx2, bbby1, bbby2, bbbz1, bbbz2
      real bbbox_tol

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

      parnm = 'Nek'
      config = '../../../../Coupling_dir/precice-config.xml'
      meshnm = 'Nek-Mesh'
      omeshn = 'Murphy-Mesh'
      rdDtNm = 'Murphy_u'
      wrDtNm = 'Nek_u'
      prcdim = 3
      !we'd like to find region spanned by this rank, let's do it in a dirty way:
      bbbx1 = 1.e30
      bbbx2 = -1.e30
      bbby1 = 1.e30
      bbby2 = -1.e30
      bbbz1 = 1.e30
      bbbz2 = -1.e30
      bbbox_tol = 0.02
      nel = nelfld(1)
      do ie = 1,nel 
         do k = 1, lx1 
         do j = 1, lx1 
         do i = 1, lx1 
         if(xm1(i,j,k,ie).lt.bbbx1) then 
            bbbx1 = xm1(i,j,k,ie)
         endif 
         if(xm1(i,j,k,ie).gt.bbbx2) then 
            bbbx2 = xm1(i,j,k,ie)
         endif

         if(ym1(i,j,k,ie).lt.bbby1) then 
            bbby1 = ym1(i,j,k,ie)
         endif 
         if(ym1(i,j,k,ie).gt.bbby2) then 
            bbby2 = ym1(i,j,k,ie)
         endif

         if(zm1(i,j,k,ie).lt.bbbz1) then 
            bbbz1 = zm1(i,j,k,ie)
         endif 
         if(zm1(i,j,k,ie).gt.bbbz2) then 
            bbbz2 = zm1(i,j,k,ie)
         endif
      enddo
      enddo
      enddo


      enddo
      print *, bbbx1, bbbx2
      print *, bbby1, bbby2
      print *, bbbz1, bbbz2

      prcbox(0) = bbbx1 -bbbox_tol
      prcbox(1) = bbbx2 +bbbox_tol
      prcbox(2) = bbby1 -bbbox_tol
      prcbox(3) = bbby2 +bbbox_tol
      prcbox(4) = bbbz1 -bbbox_tol
      prcbox(5) = bbbz2 +bbbox_tol

      ! prcbox(0) = 0.6
      ! prcbox(1) = 2.6
      ! prcbox(2) = 1.
      ! prcbox(3) = 3.
      ! prcbox(4) = 0.5
      ! prcbox(5) = 1.4
      call precicef_create(parnm,config,nid_global,np_global)
      call prc_omesh
      print *, "number of vert", prcnve 
      call precicef_set_vertices(meshnm,prcnve, prcvrt, prcvid)
      call precicef_set_mesh_access_region(omeshn, prcbox)
      call precicef_initialize()
      
      call precicef_get_mesh_vertex_size(omeshn, omshdi)
      call precicef_get_mesh_vertex_ids_and_coordinates(omeshn,
     & omshdi,prcvi2, prcvr2)
      print *, nid,'other mesh dim', omshdi
      

      call assign_corner
      call find_boundaries
      print *, 'ok here'
      call setup_interp(nekcomm)
      print *, 'ok here2'

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_solve

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      INCLUDE 'PRECIC'
      INCLUDE 'PARALLEL'
      integer i,j
      real*8 min_dt
      real*8 precice_dt
      real*8 solver_dt

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

      do kstep=1,nsteps,msteps
         ! solver_dt = 0.00005

         ! min_dt = 1e-10
         ! call precicef_get_max_time_step_size(precice_dt)
         ! if((precice_dt-solver_dt).lt.min_dt) then
         !    prcdt = precice_dt
         ! else 
         !    prcdt = solver_dt
         ! ENDIF 
         ! print *, 'solver dt ', solver_dt
         ! prcdt = solver_dt



         rdDtNm = 'Murphy_u'
         rdDtNa = 'Murphy_w'
         ! rdDtNb = 'Data1_gy'
         ! rdDtNc = 'Data1_gz'
         
         
         call settime  !compute dt
            
         call precicef_read_data(meshnm, rdDtNm, prcnve, 
     &    prcvid ,dt, prcrdt)
         call precicef_read_data(meshnm, rdDtNa, prcnve, 
     &    prcvid ,dt, prcrdx)
         call compute_bc(kstep)
         
   !       call precicef_read_data(meshnm, rdDtNb, prcnve, 
   !   &    prcvid ,prcdt, prcrdy)
   !       call precicef_read_data(meshnm, rdDtNc, prcnve, 
   !   &    prcvid ,prcdt, prcrdz)
!          do 120 i = 0,1000
!             print *, prcrdx(i), prcrdy(i), prcrdz(i)
!  120        continue
   !       do 23 i = 0, prcnve
   !          print *, prcvrt(i*3), prcvrt(i*3+1), prcvrt(i*3+2)
   !          print *, prcvrt(i*3)-prcrdt(i*3), prcvrt(i*3+1)-
   !   &       prcrdt(i*3+1), prcvrt(i*3+2)-prcrdt(i*3+2)
   !          print *, "   "
   ! 23        continue
         
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
         print *, "ici 1"
         call interpolate_u
         print *, "ici2"   
         print *, omshdi, 'omshdi'
         call precicef_write_data(omeshn, wrDtNm, omshdi, prcvi2 ,
     &     prcwdt)
         print *, 'ici3'
         call precicef_advance(dt)
         print *, 'ici4'
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

      !if (iftran) call settime
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
      call precicef_finalize()

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
