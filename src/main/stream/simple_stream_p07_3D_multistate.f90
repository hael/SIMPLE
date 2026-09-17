!@descr: stream pipeline stage 7 — multistate 3D reconstruction/refinement of pooled particles
!==============================================================================
! MODULE: simple_stream_p07_3D_multistate
!
! PURPOSE:
!   Drives the continuous multistate 3D reconstruction/refinement loop for the
!   streaming pipeline. Watches for completed 3D-export sets (written by
!   stream_p06_pool2D) and imports them into a growing pool project.
!
!   For now this stage only imports newly available sets into the pool; no
!   3D reconstruction/refinement is performed yet.
!
! ENTRY POINT:
!   stream_p07_3D_multistate%execute(cline) — called by the stream master
!
! INTERNAL SUBROUTINES:
!   import_sets_into_pool — read new exported sets into the pool
!   sigterm_handler       — SIGTERM handler: sets l_terminate for graceful exit
!
! DEPENDENCIES:
!   simple_stream_api, unix
!==============================================================================
module simple_stream_p07_3D_multistate
use unix,                    only: SIGTERM
use simple_commanders_cavgs, only: commander_model_cavgs_rejection
use simple_gui_utils,        only: mrc2jpeg_tiled
use simple_qsys_env,         only: qsys_env
use simple_stream_api
implicit none

public :: stream_p07_3D_multistate
private
#include "simple_local_flags.inc"

integer, parameter       :: NSTATES3D  = 3                 ! number of classes for abinitio3D
integer, parameter       :: NSTAGES3D  = 5                 ! number of stages for abinitio3D

type, extends(commander_base) :: stream_p07_3D_multistate
  contains
    procedure :: execute => exec_stream_p07_3D_multistate
end type stream_p07_3D_multistate

contains

    ! Manages multistate 3D reconstruction/refinement
    subroutine exec_stream_p07_3D_multistate( self, cline )
        class(stream_p07_3D_multistate), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)          :: params
        type(rec_list)            :: setslist
        type(stream_watcher)      :: project_buff
        type(sp_project)          :: spproj_glob
        type(qsys_env)            :: qenv
        type(string), allocatable :: projects(:)
        integer                   :: i, nprojects, nimported, nptcls_glob, abinitio_stage, refine_stage
        integer                   :: envlen, refine_it, nptcls_at_last_refine
        character(len=STDLEN)     :: preproc_part_env
        logical                   :: l_terminate, l_pause_ingestion
        volatile :: l_terminate ! set asynchronously by sigterm_handler
        l_terminate           = .false.
        abinitio_stage        = 0
        refine_stage          = 0
        nptcls_glob           = 0
        refine_it             = 0
        nptcls_at_last_refine = 0
        call signal(SIGTERM, sigterm_handler)   ! graceful shutdown on SIGTERM
        call cline%set('oritype', 'mic')
        call cline%set('mkdir',   'yes')
        ! generate own project file if projfile isnt set
        if( .not.cline%defined('projfile') )then
            call cline%set('projname', 'stream_3Dmultistate')
            call cline%set('projfile', 'stream_3Dmultistate.simple')
            call spproj_glob%update_projinfo(cline)
            call spproj_glob%update_compenv(cline)
            call spproj_glob%write
        endif
        call create_stream_project(spproj_glob, cline, string('3Dmultistate'))
        ! master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call get_environment_variable(SIMPLE_STREAM_PREPROC_PARTITION, preproc_part_env, envlen)
        if(envlen > 0) then
            call qenv%new(params, 1,stream=.true.,qsys_partition=string(trim(preproc_part_env)))
        else
            call qenv%new(params, 1,stream=.true.)
        end if
        ! wait if dir_target doesn't exist yet
        call wait_for_folder2(params%dir_target)
        call wait_for_folder2(params%dir_target//'/spprojs_completed')
        ! master project file
        call spproj_glob%read( params%projfile )
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_3D_multistate must start from an empty project (eg from root project folder)')
        ! project watcher
        project_buff = stream_watcher(LONGTIME, params%dir_target//'/'//DIR_STREAM_COMPLETED, spproj=.true., nretries=10)
        ! Infinite loop
        nprojects = 0 ! # of projects per iteration
        nimported = 0 ! # of sets per iteration
        l_pause_ingestion = .false.
        do
            if( file_exists(TERM_STREAM) .or. l_terminate ) then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            ! detection of new projects
            call project_buff%watch(nprojects, projects)
            if( nprojects > 0 )then
                ! zero-padded names sort lexicographically in export order
                call lex_sort(projects)
                ! memoize detected projects
                call project_buff%add2history(projects)
                do i = 1,nprojects
                    call setslist%push2chunk_list(projects(i), setslist%size() + 1, .true.)
                enddo
            endif
            ! Import new particles, paused while abinitio3D or refine3D is running
            if( .not. l_pause_ingestion ) call import_sets_into_pool( nimported )
            ! abinitio stage
            if( abinitio_stage < 2 .and. spproj_glob%os_ptcl2D%get_noris() /= 0 ) then
                if( abinitio_stage == 0 ) then
                    l_pause_ingestion = .true.
                    ! start abinitio 3D
                    call start_abinitio3D(spproj_glob, string('abinitio3D'), 5000)
                    abinitio_stage = 1
                else if( abinitio_stage == 1 ) then
                  ! Test for abinitio 3D completion
                  if( file_exists(string('abinitio3D')//'/'//TASK_FINISHED) ) then
                      ! stage complete
                      call finish_abinitio3D(spproj_glob, string('abinitio3D'))
                      abinitio_stage    = 2
                      l_pause_ingestion = .false.
                  end if
                end if
            end if
            ! refine stage
            if( abinitio_stage == 2 .and. spproj_glob%os_ptcl2D%get_noris() /= 0 ) then
                if( refine_stage == 0 ) then
                    ! only enter if the particle count has grown since the last refine stage
                    if( spproj_glob%os_ptcl2D%get_noris() > nptcls_at_last_refine ) then
                        refine_it             = refine_it + 1
                        nptcls_at_last_refine = spproj_glob%os_ptcl2D%get_noris()
                        write(logfhandle,'(A,I0)')'>>> ENTERING REFINE STAGE ', refine_it
                        l_pause_ingestion = .true.
                        ! start refine 3D
                     !   call start_refine3D(spproj_glob, string('refine3D/it_')//int2str(refine_it), 5000)
                        call spproj_glob%write(string('refine_') // int2str(refine_it) // METADATA_EXT)
                        refine_stage = 1
                    end if
                else if( refine_stage == 1 ) then
                    ! Test for refine stage completion
                !    if( file_exists(string('refine3D/it_')//int2str(refine_it)//'/'//TASK_FINISHED) ) then
                !        ! stage complete
                !        call finish_refine3D(spproj_glob, string('refine3D')//int2str(refine_it))
                !        refine_stage = 0
                !    end if
                    l_pause_ingestion = .false.
                    refine_stage = 0 ! for testing purposes
                end if
            end if
            ! Wait
            call sleep(WAITTIME)
        enddo
        ! Cleanup and final project
        call spproj_glob%kill
        call qsys_cleanup(params)
        ! end gracefully
        call simple_end('**** SIMPLE_STREAM_3D_MULTISTATE NORMAL STOP ****')
        contains

            ! imports new sets of exported particles into the pool
            subroutine import_sets_into_pool( nimported )
                integer,           intent(out) :: nimported
                type(sp_project),  allocatable :: spprojs(:)
                type(rec_iterator)             :: it
                type(chunk_rec)                :: crec
                logical, allocatable :: l_processed(:), l_imported(:), l_include_now(:)
                integer :: nsets2import, iset, nptcls2import, nmics2import, pool_nmics, nptcls
                integer :: i, fromp, imic, ind, iptcl, jptcl, jmic, nptcls_sel_tot, first_iset
                type(cmdline)                         :: cline_quality
                type(commander_model_cavgs_rejection) :: xmodel_cavgs_rejection
                type(string)                          :: quality_dir, quality_projfile, cwd_before, selected_jpeg, rejected_jpeg
                nimported = 0
                if( setslist%size() == 0 ) return
                l_processed = setslist%get_processed_flags()
                l_imported  = setslist%get_included_flags()
                nsets2import = count(l_processed(:).and.(.not.l_imported(:)))
                if( nsets2import == 0 ) return
                allocate(l_include_now(setslist%size()), source=.false.)
                ! read sets in
                allocate(spprojs(setslist%size()))
                nptcls2import = 0
                nmics2import  = 0
                it            = setslist%begin()
                do iset = 1,setslist%size()
                    call it%get(crec)
                    if( crec%included .or. (.not.crec%processed .or. crec%busy) )then
                        ! move iterator
                        call it%next()
                        cycle
                    endif
                    quality_dir      = string('quality_selection/') // get_fbody(basename(crec%projfile), 'simple')
                    quality_projfile = quality_dir // '/' // basename(crec%projfile)
                    call simple_mkdir('quality_selection')
                    call simple_mkdir(quality_dir)
                    call simple_copy_file(crec%projfile, quality_projfile)
                    ! class-average quality selection, applied to the copy before import
                    call cline_quality%set('projfile',       basename(crec%projfile))
                    call cline_quality%set('mskdiam',        params%mskdiam)
                    call cline_quality%set('quality_mode',   'apply')
                    call cline_quality%set('rejection_type', 'pool')
                    call cline_quality%set('mkdir',          'no')
                    call simple_getcwd(cwd_before)
                    call simple_chdir(quality_dir)
                    call xmodel_cavgs_rejection%execute(cline_quality)
                    call mrc2jpeg_tiled(string('quality_selected_cavgs')//MRC_EXT, string('quality_selected_cavgs')//JPG_EXT)
                    call mrc2jpeg_tiled(string('quality_rejected_cavgs')//MRC_EXT, string('quality_rejected_cavgs')//JPG_EXT)
                    selected_jpeg = cwd_before//'/'//quality_dir//'/quality_selected_cavgs'//JPG_EXT
                    if( file_exists(selected_jpeg) ) then
                        write(logfhandle, '(A)')   '>>> QUALITY SELECTED CLASS AVERAGES'
                        write(logfhandle, '(A,A)') '>>> JPEG ', selected_jpeg%to_char()
                    endif
                    rejected_jpeg = cwd_before//'/'//quality_dir//'/quality_rejected_cavgs'//JPG_EXT
                    if( file_exists(rejected_jpeg) ) then
                        write(logfhandle, '(A)')   '>>> QUALITY REJECTED CLASS AVERAGES'
                        write(logfhandle, '(A,A)') '>>> JPEG ', rejected_jpeg%to_char()
                    endif
                    call simple_chdir(cwd_before)
                    ! read the entire project file for the current set
                    call spprojs(iset)%read(quality_projfile)
                    nmics2import  = nmics2import  + spprojs(iset)%os_mic%get_noris()
                    nptcls2import = nptcls2import + spprojs(iset)%os_ptcl3D%get_noris()
                    l_include_now(iset) = .true.
                    call it%next()
                enddo
                ! reallocations
                pool_nmics = spproj_glob%os_mic%get_noris()
                nptcls     = spproj_glob%os_ptcl3D%get_noris()
                if( pool_nmics == 0 )then
                    call spproj_glob%os_mic%new(nmics2import, is_ptcl=.false.)
                    call spproj_glob%os_stk%new(nmics2import, is_ptcl=.false.)
                    call spproj_glob%os_ptcl2D%new(nptcls2import, is_ptcl=.true.)
                    call spproj_glob%os_ptcl3D%new(nptcls2import, is_ptcl=.true.)
                    ! cls2D is only transferred once, to seed the pool project on the initial import
                    first_iset = findloc(l_include_now, .true., dim=1)
                    if( first_iset > 0 ) spproj_glob%os_cls2D = spprojs(first_iset)%os_cls2D
                    fromp = 1
                else
                    call spproj_glob%os_mic%reallocate(pool_nmics+nmics2import)
                    call spproj_glob%os_stk%reallocate(pool_nmics+nmics2import)
                    call spproj_glob%os_ptcl2D%reallocate(nptcls+nptcls2import)
                    call spproj_glob%os_ptcl3D%reallocate(nptcls+nptcls2import)
                    fromp = spproj_glob%os_stk%get_top(pool_nmics)+1
                endif
                ! parameters transfer
                imic           = pool_nmics
                nptcls_sel_tot = 0
                it             = setslist%begin()
                do iset = 1,setslist%size()
                    call it%get(crec)
                    if( crec%included .or. (.not.crec%processed .or. crec%busy) .or. .not.l_include_now(iset) )then
                        ! move iterator
                        call it%next()
                        cycle
                    endif
                    ind = 1
                    do jmic = 1,spprojs(iset)%os_mic%get_noris()
                        imic = imic + 1
                        ! micrograph
                        call spproj_glob%os_mic%transfer_ori(imic, spprojs(iset)%os_mic, jmic)
                        ! stack
                        call spproj_glob%os_stk%transfer_ori(imic, spprojs(iset)%os_stk, jmic)
                        nptcls = spprojs(iset)%os_stk%get_int(jmic,'nptcls')
                        call spproj_glob%os_stk%set(imic, 'fromp', fromp)
                        call spproj_glob%os_stk%set(imic, 'top',   fromp+nptcls-1)
                        ! particles
                        do i = 1,nptcls
                            iptcl = fromp + i - 1
                            jptcl = ind   + i - 1
                            call spproj_glob%os_ptcl2D%transfer_ori(iptcl, spprojs(iset)%os_ptcl2D, jptcl)
                            call spproj_glob%os_ptcl2D%set_stkind(iptcl, imic)
                            call spproj_glob%os_ptcl3D%transfer_ori(iptcl, spprojs(iset)%os_ptcl3D, jptcl)
                            call spproj_glob%os_ptcl3D%set_stkind(iptcl, imic)
                        enddo
                        ind   = ind   + nptcls
                        fromp = fromp + nptcls
                    enddo
                    nptcls_sel_tot = nptcls_sel_tot + spprojs(iset)%os_ptcl3D%get_noris(consider_state=.true.)
                    ! display
                    write(logfhandle,'(A,I6,A,I6)')'>>> TRANSFERRED ',spprojs(iset)%os_ptcl3D%get_noris(consider_state=.true.),' PARTICLES FROM SET ',crec%id
                    call flush(logfhandle)
                    ! global list update
                    crec%included = .true.
                    call setslist%replace_iterator(it, crec)
                    ! move iterator
                    call it%next()
                enddo
                nimported   = count(l_include_now)
                nptcls_glob = nptcls_glob + nptcls_sel_tot
                ! cleanup
                do iset = 1,setslist%size()
                    call spprojs(iset)%kill
                enddo
                if( allocated(l_imported)     ) deallocate(l_imported)
                if( allocated(l_processed)    ) deallocate(l_processed)
                if( allocated(l_include_now)  ) deallocate(l_include_now)
                deallocate(spprojs)
            end subroutine import_sets_into_pool

            ! Run ab-initio 3D classification 
            subroutine start_abinitio3D( spproj_stage, outdir, mskdiam_in )
                type(sp_project),   intent(inout) :: spproj_stage
                type(string),          intent(in) :: outdir
                integer,               intent(in) :: mskdiam_in
                type(cmdline)                     :: cline_abinitio3D
                type(string)                     :: cwd, cwd_abinitio3D, server_address
                call simple_getcwd(cwd)
                call simple_mkdir('abinitio3D')
                call simple_mkdir(outdir)
                call simple_chdir(outdir)
                call simple_getcwd(cwd_abinitio3D)
                CWD_GLOB       = cwd_abinitio3D%to_char()
                server_address = qenv%get_persistent_worker_server_address()
                call spproj_stage%write(string('abinitio3D.simple'))
                call cline_abinitio3D%kill()
                call cline_abinitio3D%set('prg',              'abinitio3D')
                call cline_abinitio3D%set('mkdir',                    'no')
                call cline_abinitio3D%set('pgrp',                     'c1')
                call cline_abinitio3D%set('nstates',             NSTATES3D)
                call cline_abinitio3D%set('lpstart',                    50) ! 20
                call cline_abinitio3D%set('lpstop',                     10) ! 6
                call cline_abinitio3D%set('force_lp_range',          'yes')
                call cline_abinitio3D%set('mskdiam',            mskdiam_in)
                call cline_abinitio3D%set('nparts',                      4)
                call cline_abinitio3D%set('nthr',                       16)
                call cline_abinitio3D%set('nstages',             NSTAGES3D)
                call cline_abinitio3D%set('projfile',  'abinitio3D.simple')
                call cline_abinitio3D%printline()
                call qenv%exec_simple_prg_in_queue_async( cline_abinitio3D, string('./distr_abinitio3D'), string('simple_log_abinitio3D'), exec_bin=string('simple_exec') )
                call simple_chdir(cwd)
                CWD_GLOB = cwd%to_char()
            end subroutine start_abinitio3D

            subroutine finish_abinitio3D( spproj_stage, outdir )
                type(sp_project), intent(inout) :: spproj_stage
                type(string),        intent(in) :: outdir
                type(string)                    :: cwd
                if( .not. file_exists(outdir) ) THROW_HARD('Output directory does not exist: ')
                call simple_getcwd(cwd)
                call simple_chdir(outdir)
                call spproj_stage%kill()
                call spproj_stage%read(string('abinitio3D.simple')) ! read the project with abinitio3D output
                call simple_chdir(cwd)
            end subroutine finish_abinitio3D

            ! Run refine 3D classification 
            subroutine start_refine3D( spproj_stage, outdir, mskdiam_in )
                type(sp_project),   intent(inout) :: spproj_stage
                type(string),          intent(in) :: outdir
                integer,               intent(in) :: mskdiam_in
                type(cmdline)                     :: cline_refine3D
                type(string)                      :: cwd, cwd_refine3D, server_address
                type(string)                      :: vol_fname
                real                               :: vol_smpd
                integer                           :: vol_box, state, nstates_stage
                call simple_getcwd(cwd)
                call simple_mkdir('refine3D')
                call simple_mkdir(outdir)
                call simple_chdir(outdir)
                call simple_getcwd(cwd_refine3D)
                CWD_GLOB       = cwd_refine3D%to_char()
                server_address = qenv%get_persistent_worker_server_address()
                nstates_stage  = spproj_stage%os_ptcl3D%get_n('state')
                call spproj_stage%write(string('refine3D.simple'))
                call cline_refine3D%kill()
                call cline_refine3D%set('prg',                'refine3D')
                call cline_refine3D%set('mkdir',                    'no')
                call cline_refine3D%set('balance',                  'no')
                call cline_refine3D%set('frac_best',                 1.0)
                call cline_refine3D%set('fillin',                   'no')
                call cline_refine3D%set('update_frac',               1.0)
                call cline_refine3D%set('trail_rec',                'np')
                call cline_refine3D%set('maxits',                      1)
                call cline_refine3D%set('refine',               'greedy')
                call cline_refine3D%set('greedy_sampling',         'yes')
                call cline_refine3D%set('update_missing',          'yes')
                call cline_refine3D%set('pgrp',                     'c1')
                call cline_refine3D%set('nstates',         nstates_stage)
                call cline_refine3D%set('lpstart',                    20)
                call cline_refine3D%set('lpstop',                      6)
                call cline_refine3D%set('force_lp_range',          'yes')
                call cline_refine3D%set('mskdiam',            mskdiam_in)
                call cline_refine3D%set('nparts',                      4)
                call cline_refine3D%set('nthr',                       16)
                call cline_refine3D%set('nstages',                     1)
                call cline_refine3D%set('sigma_est',            'global')
                call cline_refine3D%set('projfile',    'refine3D.simple')
                ! state volumes from the abinitio3D output segment
                do state = 1, nstates_stage
                    if( .not. spproj_stage%isthere_in_osout('vol', state) )then
                        THROW_HARD('missing state volume in project out segment; start_refine3D')
                    endif
                    call spproj_stage%get_vol('vol', state, vol_fname, vol_smpd, vol_box)
                    call cline_refine3D%set('vol'//int2str(state), vol_fname)
                enddo
                call cline_refine3D%printline()
                call qenv%exec_simple_prg_in_queue_async( cline_refine3D, string('./distr_refine3D'), string('simple_log_refine3D'), exec_bin=string('simple_exec') )
                call simple_chdir(cwd)
                CWD_GLOB = cwd%to_char()
                call vol_fname%kill
            end subroutine start_refine3D

            subroutine finish_refine3D( spproj_stage, outdir )
                type(sp_project), intent(inout) :: spproj_stage
                type(string),        intent(in) :: outdir
                type(string)                    :: cwd
                if( .not. file_exists(outdir) ) THROW_HARD('Output directory does not exist: ')
                call simple_getcwd(cwd)
                call simple_chdir(outdir)
                call spproj_stage%kill()
                call spproj_stage%read(string('refine3D.simple')) ! read the project with refine3D output
                call simple_chdir(cwd)
            end subroutine finish_refine3D

            ! Called asynchronously on SIGTERM. Exits immediately after logging.
            subroutine sigterm_handler()
                write(logfhandle, '(A)') 'SIGTERM RECEIVED'
                l_terminate = .true.
            end subroutine sigterm_handler

    end subroutine exec_stream_p07_3D_multistate

end module simple_stream_p07_3D_multistate
