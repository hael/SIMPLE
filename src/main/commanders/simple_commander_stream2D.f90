! concrete commander: cluster2D_stream for streaming 2D alignment and clustering of single-particle images
module simple_commander_stream2D
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_commander_base,     only: commander_base
use simple_parameters,         only: parameters
use simple_sp_project,         only: sp_project
use simple_guistats,           only: guistats
use simple_stream_utils
use simple_qsys_funs
use simple_qsys_env
use simple_commander_cluster2D_stream
use simple_moviewatcher
use simple_progress
use simple_nice
implicit none

private
#include "simple_local_flags.inc"

public :: commander_stream_sieve_cavgs

type, extends(commander_base) :: commander_stream_sieve_cavgs
  contains
    procedure :: execute => exec_sieve_cavgs
end type commander_stream_sieve_cavgs

! module constants
character(len=STDLEN), parameter :: DIR_STREAM_COMPLETED = trim(PATH_HERE)//'spprojs_completed/' ! location for projects processed
integer,               parameter :: LONGTIME  = 60                                               ! time lag after which a movie/project is processed
integer,               parameter :: WAITTIME  = 10                                               ! movie folder watched every WAITTIME seconds
integer(kind=dp),      parameter :: FLUSH_TIMELIMIT = 900 ! time (secs) after which leftover particles join the pool IF the 2D analysis is paused

contains

    subroutine exec_sieve_cavgs( self, cline )
        class(commander_stream_sieve_cavgs), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        character(len=STDLEN),     parameter   :: micspproj_fname = './streamdata.simple'
        integer,                   parameter   :: PAUSE_NITERS    = 5   ! # of iterations after which 2D analysis is paused
        integer,                   parameter   :: PAUSE_TIMELIMIT = 600 ! time (secs) after which 2D analysis is paused
        type(projrecord),          allocatable :: projrecords(:)
        type(parameters)                       :: params
        type(qsys_env)                         :: qenv
        type(simple_nice_communicator)         :: nice_communicator
        type(projs_list)                       :: chunkslist, setslist
        type(guistats)                         :: gui_stats
        type(oris)                             :: moldiamori, chunksizeori
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj_glob
        type(json_core)                        :: json
        type(json_value),          pointer     :: snapshot_json
        character(len=LONGSTRLEN), allocatable :: projects(:)
        character(len=LONGSTRLEN)              :: cwd_job
        character(len=STDLEN)                  :: chunk_part_env
        real                                   :: moldiam
        integer(kind=dp)                       :: time_last_import
        integer                                :: nchunks_glob, nchunks_imported, nprojects, iter
        integer                                :: n_imported, n_imported_prev, n_added, nptcls_glob, n_failed_jobs
        integer                                :: i, envlen
        logical                                :: l_pause, l_params_updated
        nullify(snapshot_json)
        call cline%set('oritype',      'mic')
        call cline%set('mkdir',        'yes')
        call cline%set('autoscale',    'yes')
        call cline%set('reject_mics',  'no')
        call cline%set('kweight_chunk','default')
        call cline%set('prune',        'no')
        call cline%set('wiener',       'full')
        call cline%set('reject_cls',   'no')
        call cline%set('remove_chunks','no')
        call cline%set('refine',       'snhc_smpl')
        call cline%set('ml_reg',       'no')
        call cline%set('objfun',       'euclid')
        if( .not. cline%defined('dir_target')   ) THROW_HARD('DIR_TARGET must be defined!')
        if( .not. cline%defined('walltime')     ) call cline%set('walltime',     29*60) ! 29 minutes
        if( .not. cline%defined('dynreslim')    ) call cline%set('dynreslim',    'no')
        if( .not. cline%defined('nchunksperset')) call cline%set('nchunksperset', 2)
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ! sanity check for restart
        if( cline%defined('dir_exec') )then
            if( .not.file_exists(cline%get_carg('dir_exec')) )then
                THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
            endif
        endif      
        ! master parameters
        call cline%set('numlen', 5)
        call cline%set('stream','yes')
        call params%new(cline)
        params%nthr2D       = params%nthr
        params%ml_reg_chunk = trim(params%ml_reg)
        call simple_getcwd(cwd_job)
        call cline%set('mkdir', 'no')
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        ! restart
        if( cline%defined('dir_exec') )then
            call cline%delete('dir_exec')
            call del_file(micspproj_fname)
            call cleanup_root_folder
        endif
        ! mskdiam
        if( .not. cline%defined('mskdiam') ) then
            ! nice communicator status
            nice_communicator%stat_root%stage = "waiting for mask diameter"
            call nice_communicator%cycle()
            write(logfhandle,'(A,F8.2)')'>>> WAITING UP TO 5 MINUTES FOR '//trim(STREAM_MOLDIAM)
            do i=1, 30
                if(file_exists(trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM))) exit
                call sleep(10)
            end do
            if( .not. file_exists(trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM))) THROW_HARD('either mskdiam must be given or '// trim(STREAM_MOLDIAM) // ' exists in target_dir')
            ! read mskdiam from file
            call moldiamori%new(1, .false.)
            call moldiamori%read( trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM) )
            if( .not. moldiamori%isthere(1, "moldiam") ) THROW_HARD( 'moldiam missing from ' // trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM) )
            moldiam = moldiamori%get(1, "moldiam")
           
            ! write acopy for stream 3d
            call moldiamori%write(1, trim(STREAM_MOLDIAM))
            call moldiamori%kill
            params%mskdiam = moldiam * 1.2
            call cline%set('mskdiam', params%mskdiam)
            write(logfhandle,'(A,F8.2)')'>>> MASK DIAMETER SET TO', params%mskdiam
        endif
        ! Computing environment
        call get_environment_variable(SIMPLE_STREAM_CHUNK_PARTITION, chunk_part_env, envlen)
        if(envlen > 0) then
            call qenv%new(1, exec_bin='simple_exec', qsys_partition=trim(chunk_part_env))
        else
            call qenv%new(1, exec_bin='simple_exec')
        end if
        ! Resolution based class rejection
        call set_lpthres_type("off")
        ! Number of particles per class
        if(params%nptcls_per_cls == 0) write(logfhandle,'(A)')'>>> # PARTICLES PER CLASS WILL BE AUTO DETERMINED AFTER 100 IMPORTED MICROGRAPHS'
        ! initialise progress monitor
        call progressfile_init()
        ! master project file
        call spproj_glob%read( params%projfile )
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_cluster2D must start from an empty project (eg from root project folder)')
        ! movie watcher init
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true.)
        ! Infinite loop
        nptcls_glob      = 0       ! global number of particles
        nchunks_glob     = 0       ! global number of completed chunks
        n_imported       = 0       ! global number of imported processed micrographs
        n_imported_prev  = 0
        nprojects        = 0
        iter             = 0       ! global number of infinite loop iterations
        n_failed_jobs    = 0
        l_pause          = .false.
        time_last_import = huge(time_last_import)      ! used for flushing unprocessed particles
        ! guistats init
        call gui_stats%init(.true.)
        call gui_stats%set('particles', 'particles_imported',          0,            primary=.true.)
        call gui_stats%set('2D',        'iteration',                   0,            primary=.true.)
        ! Joe: nparts is not an input, also see project_buff below
        call gui_stats%set('compute',   'compute_in_use',      int2str(0) // '/' // int2str(params%nparts), primary=.true.)
        ! nice
        nice_communicator%stat_root%stage = "importing particles"
        call nice_communicator%update_cls2D(particles_imported=0)
        call nice_communicator%cycle()
        do
            if( file_exists(trim(TERM_STREAM)) .or. nice_communicator%exit)then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            if( nice_communicator%stop .or. test_repick() ) then
                if(test_repick()) call write_repick_refs("../repick_refs.mrc")
                ! termination
                write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                call spproj_glob%kill
                call qsys_cleanup
                call nice_communicator%terminate(stop=.true.)
                call simple_end('**** SIMPLE_STREAM2D USER STOP ****')
                call EXIT(0)
            endif
            iter = iter + 1
            ! detection of new projects
            call project_buff%watch(nprojects, projects, max_nmovies=10*params%nparts)
            ! update global records
            if( nprojects > 0 )then
                call update_records_with_project(projects, n_imported )
                call project_buff%add2history(projects)
            endif
            ! project update
            if( nprojects > 0 )then
                n_imported = size(projrecords)
                write(logfhandle,'(A,I6,I8)') '>>> # MICROGRAPHS / PARTICLES IMPORTED : ',n_imported, nptcls_glob
                ! guistats
                call gui_stats%set('particles', 'particles_imported', int2commastr(nptcls_glob), primary=.true.)
                call gui_stats%set_now('particles', 'last_particles_imported')
                call nice_communicator%update_cls2D(particles_imported=nptcls_glob, last_particles_imported=.true.)
                ! update progress monitor
                call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                time_last_import = time8()
                if( n_imported < 1000 )then
                    call update_user_params(cline)
                else if( n_imported > n_imported_prev + 100 )then
                    call update_user_params(cline)
                    n_imported_prev = n_imported
                endif
            endif
            ! Chunk book-keeping section
            nice_communicator%view_cls2D%mskdiam  = params%mskdiam
            nice_communicator%view_cls2D%boxsizea = get_boxa()
            call update_user_params2D(cline, l_params_updated, nice_communicator%update_arguments)
            if( l_params_updated ) l_pause = .false.    ! resuming
            call update_chunks
            call memoize_chunks(chunkslist, nchunks_imported)
            call update_user_params2D(cline, l_params_updated, nice_communicator%update_arguments)
            if( l_params_updated ) l_pause = .false.    ! resuming
            if( nchunks_imported > 0 )then
                nchunks_glob         = nchunks_glob + nchunks_imported
                l_pause              = .false.   ! resuming
                ! build sets
                call generate_sets(chunkslist, setslist)
            endif
            ! Sets analysis section
            if( setslist%n > 0 )then
                ! processing with cluster_cavgs / match_cavgs
                if( setslist%processed(1) .and. (setslist%n > 1) )then
                    do i = 2,setslist%n
                        call is_set_processed(i)
                    enddo
                    call submit_match_cavgs
                else
                    call is_set_processed(1)
                    call submit_cluster_cavgs
                endif
            endif
            if( l_pause )then
                ! Whether to flush particles
                if( (time8()-time_last_import > FLUSH_TIMELIMIT) .and. all_chunks_available() )then
                    ! Remaining unclassified particles: TBD
                endif
            endif
            ! 2D analyses
            call analyze2D_new_chunks(projrecords)
            call sleep(WAITTIME)
            ! guistats
            call gui_stats%write_json
            ! nice
            if(l_pause) then
                nice_communicator%stat_root%stage = "paused chunk 2D analysis"
                nice_communicator%stat_root%user_input = .true.
            else if(get_nchunks() > 0) then
                nice_communicator%stat_root%stage = "classifying chunks"
                nice_communicator%stat_root%user_input = .false.
            end if
            call nice_communicator%cycle()
        end do
        ! termination
        nice_communicator%stat_root%stage = "terminating"
        call nice_communicator%cycle()
        call terminate_stream2D( projrecords )
        call update_user_params(cline) ! Joe: bit late for this?
        ! final stats
        call gui_stats%hide('compute', 'compute_in_use')
        call gui_stats%deactivate_section('compute')
        call gui_stats%write_json
        call gui_stats%kill
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
        call nice_communicator%terminate()
        call simple_end('**** SIMPLE_STREAM2D NORMAL STOP ****')
        contains

            ! updates global records
            subroutine update_records_with_project( projectnames, n_imported )
                character(len=LONGSTRLEN), allocatable, intent(in)  :: projectnames(:)
                integer,                                intent(out) :: n_imported
                type(sp_project),     allocatable :: spprojs(:)
                type(projrecord),     allocatable :: old_records(:)
                character(len=:),     allocatable :: fname, abs_fname
                real    :: avgmicptcls, nptcls_per_cls
                integer :: iproj, n_spprojs, n_old, irec, n_completed, nptcls, nmics, imic, n_ptcls, first
                n_imported = 0
                n_ptcls    = 0
                if( .not.allocated(projectnames) ) return
                n_spprojs  = size(projectnames)
                if( n_spprojs == 0 )return
                n_old = 0 ! on first import
                if( allocated(projrecords) ) n_old = size(projrecords)
                allocate(spprojs(n_spprojs))
                ! because pick_extract purges state=0 and nptcls=0 mics,
                ! all mics can be assumed associated with particles
                nmics = 0
                first = 0
                do iproj = 1,n_spprojs
                    call spprojs(iproj)%read_segment('mic', trim(projectnames(iproj)))
                    nmics = nmics + spprojs(iproj)%os_mic%get_noris()
                    if( (first == 0) .and. (nmics > 0) ) first = iproj
                enddo
                if( nmics == 0 ) return
                ! import micrographs
                n_completed = n_old + nmics
                n_imported  = nmics
                ! reallocate records
                if( n_old == 0 )then
                    allocate(projrecords(nmics))
                else
                    call move_alloc(projrecords, old_records)
                    allocate(projrecords(n_completed))
                    projrecords(1:n_old) = old_records(:)
                    deallocate(old_records)
                endif
                ! update global records and some global variables
                irec = n_old
                do iproj = 1,n_spprojs
                    do imic = 1,spprojs(iproj)%os_mic%get_noris()
                        irec      = irec + 1
                        nptcls    = spprojs(iproj)%os_mic%get_int(imic,'nptcls')
                        n_ptcls   = n_ptcls + nptcls ! global update
                        fname     = trim(projectnames(iproj))
                        abs_fname = simple_abspath(fname, errmsg='stream_cluster2D :: update_projects_list 1')
                        projrecords(irec)%projname   = trim(abs_fname)
                        projrecords(irec)%micind     = imic
                        projrecords(irec)%nptcls     = nptcls
                        projrecords(irec)%nptcls_sel = nptcls
                        projrecords(irec)%included   = .false.
                    enddo
                enddo
                nptcls_glob = nptcls_glob + n_ptcls ! global update
                ! Updates global parameters once and init 2D
                if(params%nptcls_per_cls == 0) then
                    if(size(projrecords) .gt. 100) then
                        avgmicptcls = nptcls_glob / size(projrecords)
                        avgmicptcls = ceiling(avgmicptcls / 10) * 10.0
                        ! these parameters may need tweaking
                        nptcls_per_cls = 1000 * (20 + (0.15 * avgmicptcls))
                        nptcls_per_cls = nptcls_per_cls / real(params%ncls_start)
                        nptcls_per_cls = ceiling(nptcls_per_cls / 100) * 100.0
                        write(logfhandle,'(A,I6)')   '>>> AVERAGE # PARTICLES PER MICROGRAPH : ', int(avgmicptcls)
                        write(logfhandle,'(A,I6,A)') '>>> USING ', int(nptcls_per_cls), ' PARTICLES PER CLASS'
                        params%nptcls_per_cls = int(nptcls_per_cls)
                        call cline%set('nptcls_per_cls', nptcls_per_cls)
                        params%smpd = spprojs(first)%os_mic%get(1,'smpd')
                        call spprojs(first)%read_segment('stk', trim(projectnames(first)))
                        params%box  = nint(spprojs(first)%os_stk%get(1,'box'))
                        call init_cluster2D_stream(cline, spproj_glob, micspproj_fname)
                        call cline%delete('ncls')
                        ! write out for stream3d to pick up
                        call chunksizeori%new(1, .false.)
                        call chunksizeori%set(1, 'nptcls_per_cls', params%nptcls_per_cls)
                        call chunksizeori%write(1, trim(STREAM_CHUNKSIZE))
                        call chunksizeori%kill
                    end if
                else if( n_old == 0 )then
                    params%smpd = spprojs(first)%os_mic%get(1,'smpd')
                    call spprojs(first)%read_segment('stk', trim(projectnames(first)))
                    params%box  = nint(spprojs(first)%os_stk%get(1,'box'))
                    call init_cluster2D_stream(cline, spproj_glob, micspproj_fname)
                    call cline%delete('ncls')
                endif
                ! cleanup
                do iproj = 1,n_spprojs
                    call spprojs(iproj)%kill
                enddo
                deallocate(spprojs)
            end subroutine update_records_with_project

            subroutine generate_sets( chunks, sets )
                class(projs_list), intent(in)    :: chunks
                class(projs_list), intent(inout) :: sets
                type(sp_project)              :: spproj
                character(len=:), allocatable :: tmpl
                integer :: navail_chunks, n, iset, ic_start, ic_end
                navail_chunks = chunks%n - sets%n * params%nchunksperset
                n = floor(real(navail_chunks) / real(params%nchunksperset))
                if( n < 1 )return
                do iset = 1,n
                    ! merge chunks project into designated folder
                    ic_start = sets%n*params%nchunksperset + 1
                    ic_end   = ic_start + params%nchunksperset - 1
                    tmpl     = 'set_'//int2str(sets%n+1)
                    call simple_mkdir(tmpl)
                    call merge_chunks(chunks%projfiles(ic_start:ic_end), tmpl, spproj, projname_out=tmpl)
                    ! update global list, also increments sets%n
                    call sets%append(tmpl//'/'//tmpl//trim(METADATA_EXT), sets%n+1, .false.)
                enddo
                call spproj%kill
            end subroutine generate_sets

            subroutine submit_cluster_cavgs
                type(cmdline)              :: cline_cluster_cavgs
                character(len=XLONGSTRLEN) :: cwd
                if( setslist%n < 1 )       return  ! no sets generated yet
                if( setslist%busy(1) )      return  ! ongoing
                if( setslist%processed(1) ) return  ! already done
                call cline_cluster_cavgs%set('prg',          'cluster_cavgs')
                call cline_cluster_cavgs%set('projfile',     basename(setslist%projfiles(1)))
                call cline_cluster_cavgs%set('mkdir',        'no')
                call cline_cluster_cavgs%set('nthr',         params%nthr)
                call cline_cluster_cavgs%set('mskdiam',      params%mskdiam)
                call cline_cluster_cavgs%set('verbose_exit', 'yes')
                call chdir(stemname(setslist%projfiles(1)))
                call simple_getcwd(cwd)
                cwd_glob = trim(cwd)
                call qenv%exec_simple_prg_in_queue_async(cline_cluster_cavgs,&
                    &'cluster_cavgs_script', 'cluster_cavgs.log')
                call chdir('..')
                call simple_getcwd(cwd_glob)
                setslist%busy(1)      = .true.
                setslist%processed(1) = .false.
                call cline_cluster_cavgs%kill
            end subroutine submit_cluster_cavgs

            subroutine submit_match_cavgs
                type(cmdline)                 :: cline_match_cavgs
                character(len=:), allocatable :: path
                character(len=XLONGSTRLEN)    :: cwd
                integer :: iset
                if( setslist%n < 2 ) return    ! not enough sets generated
                ! any unprocessed and not being processed?
                if( .not.any((.not.setslist%processed(2:)) .and. (.not.setslist%busy(2:))) ) return
                call cline_match_cavgs%set('prg',          'match_cavgs')
                call cline_match_cavgs%set('projfile',     simple_abspath(setslist%projfiles(1) ,"submit_match_cavgs"))
                call cline_match_cavgs%set('mkdir',        'no')
                call cline_match_cavgs%set('nthr',         params%nthr)
                call cline_match_cavgs%set('mskdiam',      params%mskdiam)
                call cline_match_cavgs%set('verbose_exit', 'yes')
                call cline_match_cavgs%delete('nparts')
                call cline_match_cavgs%delete('nparts_chunk')
                do iset = 2,setslist%n
                    if( setslist%processed(iset) ) cycle ! already done
                    if( setslist%busy(iset) )      cycle ! ongoing
                    call cline_match_cavgs%set('projfile_target', basename(setslist%projfiles(iset)))
                    path = stemname(setslist%projfiles(iset))
                    call chdir(path)
                    call simple_getcwd(cwd)
                    cwd_glob = trim(cwd)
                    call qenv%exec_simple_prg_in_queue_async(cline_match_cavgs,&
                        &'match_cavgs_script', 'match_cavgs.log')
                    call chdir('..')
                    call simple_getcwd(cwd_glob)
                    setslist%busy(iset)      = .true.   ! ongoing
                    setslist%processed(iset) = .false.  ! not complete
                enddo
                call cline_match_cavgs%kill
            end subroutine submit_match_cavgs

            ! Check for status of individual sets
            subroutine is_set_processed( i )
                integer, intent(in) :: i
                character(len=:), allocatable :: fname
                if( setslist%n < 1        ) return  ! no sets generated yet
                if( setslist%processed(i) ) return  ! already done
                fname = trim(stemname(setslist%projfiles(i)))//'/'//trim(TASK_FINISHED)
                if( file_exists(fname) )then
                    setslist%busy(i)      = .false.
                    setslist%processed(i) = .true.  ! now ready for pool import
                else
                    setslist%processed(i) = .false.
                endif
            end subroutine is_set_processed

            subroutine import_sets_into_pool( nimported )
                integer,          intent(out) :: nimported
                type(sp_project), allocatable :: spprojs(:)
                class(sp_project),    pointer :: pool
                integer :: nsets2import, iset, nptcls2import, nmics2import, nmics, nptcls
                integer :: i, fromp, fromp_prev, imic, ind, iptcl, jptcl, jmic, nptcls_sel
                nimported = 0
                if( setslist%n == 0 )          return
                if( .not.is_pool_available() ) return
                nsets2import = count(setslist%processed(:).and.(.not.setslist%imported(:)))
                if( nsets2import == 0 )        return
                ! read sets in
                allocate(spprojs(setslist%n))
                nptcls2import = 0
                nmics2import  = 0
                do iset = 1,setslist%n
                    if( setslist%imported(iset) )       cycle   ! already imported
                    if( .not.setslist%processed(iset) ) cycle   ! not processed yet
                    if( setslist%busy(iset) )           cycle   ! being processed
                    call spprojs(iset)%read_non_data_segments(setslist%projfiles(iset))
                    call spprojs(iset)%read_segment('mic',    setslist%projfiles(iset))
                    call spprojs(iset)%read_segment('stk',    setslist%projfiles(iset))
                    call spprojs(iset)%read_segment('ptcl2D', setslist%projfiles(iset))
                    nmics2import  = nmics2import  + spprojs(iset)%os_mic%get_noris()
                    nptcls2import = nptcls2import + spprojs(iset)%os_ptcl2D%get_noris()
                    ! global list update
                    setslist%imported(iset) = .true.
                enddo
                ! reallocations
                call get_pool_ptr(pool)
                nmics  = pool%os_mic%get_noris()
                nptcls = pool%os_ptcl2D%get_noris()
                if( nmics == 0 )then
                    call pool%os_mic%new(nmics2import, is_ptcl=.false.)
                    call pool%os_stk%new(nmics2import, is_ptcl=.false.)
                    call pool%os_ptcl2D%new(nptcls2import, is_ptcl=.true.)
                    fromp = 1
                else
                    call pool%os_mic%reallocate(nmics+nmics2import)
                    call pool%os_stk%reallocate(nmics+nmics2import)
                    call pool%os_ptcl2D%reallocate(nptcls+nptcls2import)
                    fromp = pool%os_stk%get_top(nmics)+1
                endif
                ! parameters transfer
                imic = nmics
                do iset = 1,nsets2import
                    fromp_prev = fromp
                    ind = 1
                    do jmic = 1,spprojs(iset)%os_mic%get_noris()
                        imic  = imic+1
                        ! micrograph
                        call pool%os_mic%transfer_ori(imic, spprojs(iset)%os_mic, jmic)
                        ! stack
                        call pool%os_stk%transfer_ori(imic, spprojs(iset)%os_stk, jmic)
                        nptcls = spprojs(iset)%os_stk%get_int(jmic,'nptcls')
                        call pool%os_stk%set(imic, 'fromp', fromp)
                        call pool%os_stk%set(imic, 'top',   fromp+nptcls-1)
                        ! particles
                        !$omp parallel do private(i,iptcl,jptcl) default(shared) proc_bind(close)
                        do i = 1,nptcls
                            iptcl = fromp + i - 1
                            jptcl = ind   + i - 1
                            call pool%os_ptcl2D%transfer_ori(iptcl, spprojs(iset)%os_ptcl2D, jptcl)
                            call pool%os_ptcl2D%set_stkind(iptcl, imic)
                            call pool%os_ptcl2D%set(iptcl, 'updatecnt', 0)  ! new particle
                            call pool%os_ptcl2D%set(iptcl, 'frac',      0.) ! new particle
                        enddo
                        !$omp end parallel do
                        ind   = ind   + nptcls
                        fromp = fromp + nptcls
                    enddo
                    nptcls_sel = spprojs(iset)%os_ptcl2D%get_noris(consider_state=.true.)
                    ! TODO
                    ! storing sigmas as per stack individual documents
                    ! if( l_update_sigmas ) call converged_chunks(ichunk)%split_sigmas_into(SIGMAS_DIR)
                    ! display
                    write(logfhandle,'(A,I6,A,I6)')'>>> TRANSFERRED ',nptcls_sel,' PARTICLES FROM SET ',setslist%ids(iset)
                    call flush(logfhandle)
                enddo
                nimported = nsets2import
                ! cleanup
                do iset = 1,setslist%n
                    call spprojs(iset)%kill
                enddo
                nullify(pool)
            end subroutine import_sets_into_pool

    end subroutine exec_sieve_cavgs

end module simple_commander_stream2D
