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
use simple_stream_communicator
use simple_progress
use simple_imgproc
use simple_nice
implicit none

private
#include "simple_local_flags.inc"

public :: commander_stream_sieve_cavgs, commander_stream_abinitio2D

type, extends(commander_base) :: commander_stream_sieve_cavgs
  contains
    procedure :: execute => exec_sieve_cavgs
end type commander_stream_sieve_cavgs

type, extends(commander_base) :: commander_stream_abinitio2D
  contains
    procedure :: execute => exec_stream_abinitio2D
end type commander_stream_abinitio2D

! module constants
character(len=STDLEN), parameter :: DIR_STREAM_COMPLETED = trim(PATH_HERE)//'spprojs_completed/' ! location for processed projects
character(len=STDLEN), parameter :: micspproj_fname = './streamdata.simple'
integer,               parameter :: LONGTIME        = 60    ! time lag after which a movie/project is processed
integer,               parameter :: WAITTIME        = 10    ! movie folder watched every WAITTIME seconds
integer(kind=dp),      parameter :: FLUSH_TIMELIMIT = 900   ! time (secs) after which leftover particles join the pool IF the 2D analysis is paused
integer,               parameter :: PAUSE_NITERS    = 5     ! # of iterations after which 2D analysis is paused
integer,               parameter :: PAUSE_TIMELIMIT = 600   ! time (secs) after which 2D analysis is paused

contains

    ! Manages individual chunks/sets classification, matching & rejection
    ! TODO: handling of un-classified particles
    subroutine exec_sieve_cavgs( self, cline )
        class(commander_stream_sieve_cavgs), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(projrecord),          allocatable :: projrecords(:)
        type(parameters)                       :: params
        type(qsys_env)                         :: qenv
        type(simple_nice_communicator)         :: nice_communicator
        type(stream_http_communicator)         :: http_communicator
        type(projs_list)                       :: chunkslist, setslist
        type(guistats)                         :: gui_stats
        type(oris)                             :: moldiamori, chunksizeori
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj_glob
        type(json_value),          pointer     :: accepted_cls2D, rejected_cls2D
        character(len=LONGSTRLEN), allocatable :: projects(:)
        character(len=STDLEN)                  :: chunk_part_env
        character(len=:),          allocatable :: selection_jpeg
        integer,                   allocatable :: accepted_cls_ids(:), rejected_cls_ids(:)
        real             :: moldiam
        integer(kind=dp) :: time_last_import
        integer          :: nchunks_glob, nchunks_imported, nprojects, iter, i, envlen
        integer          :: n_imported, n_imported_prev, n_added, nptcls_glob, n_failed_jobs
        integer          :: n_accepted, n_rejected, jpg_ntiles, jpg_nxtiles, jpg_nytiles, xtile, ytile
        logical          :: l_params_updated, l_wait_for_user, selection_jpeg_created, found
        call cline%set('oritype',      'mic')
        call cline%set('mkdir',        'yes')
        call cline%set('autoscale',    'yes')
        call cline%set('reject_mics',  'no')
        call cline%set('reject_cls',   'no') ! refers to previous implementation
        call cline%set('kweight_chunk','default')
        call cline%set('prune',        'no')
        call cline%set('wiener',       'full')
        call cline%set('refine',       'snhc_smpl')
        call cline%set('ml_reg',       'no')
        call cline%set('objfun',       'euclid')
        call cline%set('sigma_est',    'global')
        if( .not.cline%defined('walltime')     ) call cline%set('walltime',     29*60) ! 29 minutes
        if( .not.cline%defined('dynreslim')    ) call cline%set('dynreslim',    'no')
        if( .not.cline%defined('nchunksperset')) call cline%set('nchunksperset', 2)
        if( .not.cline%defined('remove_chunks')) call cline%set('remove_chunks','yes')
        if( .not.cline%defined('center')       ) call cline%set('center',       'yes')
        if( .not.cline%defined('algorithm')    ) call cline%set('algorithm',    'abinitio2D')
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ! restart
        call cleanup4restart
        ! generate own project file if projfile isnt set
        if(cline%get_carg('projfile') .eq. '') then 
            call cline%set('projname', 'sieve_cavgs')
            call cline%set('projfile', 'sieve_cavgs.simple')
            call spproj_glob%update_projinfo(cline)
            call spproj_glob%update_compenv(cline)
            call spproj_glob%write
        endif
        ! master parameters
        call cline%set('numlen', 5)
        call cline%set('stream','yes')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        selection_jpeg_created = .false.
        l_wait_for_user        = .false.
        if(params%interactive == 'yes') l_wait_for_user = .true.
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, "")
        call nice_communicator%cycle()
        ! http communicator init
        call http_communicator%create(params%niceprocid, params%niceserver, "sieve_cavgs")
        call communicator_init()
        call http_communicator%send_jobstats()
        ! wait if dir_target doesn't exist yet
        if(.not. dir_exists(trim(params%dir_target))) then
            write(logfhandle, *) ">>> WAITING FOR ", trim(params%dir_target), " TO BE GENERATED"
            do i=1, 360
                if(dir_exists(trim(params%dir_target))) then
                    write(logfhandle, *) ">>> ", trim(params%dir_target), " FOUND"
                    exit
                endif
                call sleep(10)
            end do
        endif
        ! mskdiam
        if( .not. cline%defined('mskdiam') )then
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
        ! project watcher
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true., nretries=10)
        call simple_mkdir(trim(PATH_HERE)//trim(DIR_STREAM_COMPLETED))
        ! Infinite loop
        nptcls_glob      = 0       ! global number of particles
        nchunks_glob     = 0       ! global number of completed chunks
        n_imported       = 0       ! global number of imported processed micrographs
        n_imported_prev  = 0
        nprojects        = 0
        iter             = 0       ! global number of infinite loop iterations
        n_failed_jobs    = 0
        n_accepted       = 0       ! global number of accepted particles
        n_rejected       = 0       ! global number of rejected particles
        time_last_import = huge(time_last_import)   ! used for flushing unprocessed particles
        ! guistats init
        !call gui_stats%init(.true.)
        !call gui_stats%set('particles', 'particles_imported',          0,            primary=.true.)
        !call gui_stats%set('2D',        'iteration',                   0,            primary=.true.)
        ! Joe: nparts is not an input, also see project_buff below
        !all gui_stats%set('compute',   'compute_in_use',      int2str(0) // '/' // int2str(params%nparts), primary=.true.)
        ! nice
        nice_communicator%stat_root%stage = "importing particles"
        call nice_communicator%update_cls2D(particles_imported=0)
        call nice_communicator%cycle()

        do
            ! termination
            if( file_exists(trim(TERM_STREAM)) .or. nice_communicator%exit ) exit
            if( nice_communicator%stop .or. test_repick() ) then
                if(test_repick()) call write_repick_refs("../repick_refs.mrc")
                ! termination
                write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                call spproj_glob%kill
                call qsys_cleanup
                call nice_communicator%terminate(stop=.true.)
                call simple_end('**** SIMPLE_STREAM_SIEVE_CAVGS USER STOP ****')
                call EXIT(0)
            endif
            iter = iter + 1
            ! http stats
            call http_communicator%json%update(http_communicator%job_json, "stage",               "finding and processing new particles", found)    
            call http_communicator%json%update(http_communicator%job_json, "particles_accepted",  n_accepted,                             found)
            call http_communicator%json%update(http_communicator%job_json, "particles_rejected",  n_rejected,                             found)
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
            !    call gui_stats%set('particles', 'particles_imported', int2commastr(nptcls_glob), primary=.true.)
             !   call gui_stats%set_now('particles', 'last_particles_imported')
             !   call nice_communicator%update_cls2D(particles_imported=nptcls_glob, last_particles_imported=.true.)
                ! update progress monitor
              !  call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                ! http stats
                call http_communicator%json%update(http_communicator%job_json, "particles_imported",       nptcls_glob,      found)
                call http_communicator%json%update(http_communicator%job_json, "last_particles_imported",  stream_datestr(), found)
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
            call update_chunks
            call memoize_chunks(chunkslist, nchunks_imported)
            call update_user_params2D(cline, l_params_updated, nice_communicator%update_arguments)
            if( nchunks_imported > 0 )then
                nchunks_glob = nchunks_glob + nchunks_imported
                ! build sets
                call generate_sets(chunkslist, setslist)
            endif
            ! Sets analysis section
            if( setslist%n > 0 )then
                if( setslist%processed(1) .and. (setslist%n > 1) .and. .not. l_wait_for_user) then
                    ! all sets but the first employ match_cavgs
                    do i = 2,setslist%n
                        call is_set_processed(i)
                    enddo
                    call submit_match_cavgs
                else
                    if(l_wait_for_user .and. setslist%processed(1)) then
                        call http_communicator%json%update(http_communicator%job_json, "user_input", .true., found)
                        if(.not. selection_jpeg_created) then
                            call generate_selection_jpeg()
                            xtile = 0
                            ytile = 0
                            call http_communicator%json%remove(accepted_cls2D, destroy=.true.)
                            call http_communicator%json%create_array(accepted_cls2D, "accepted_cls2D")
                            call http_communicator%json%add(http_communicator%job_json, accepted_cls2D)
                            call http_communicator%json%remove(rejected_cls2D, destroy=.true.)
                            call http_communicator%json%create_array(rejected_cls2D, "rejected_cls2D")
                            call http_communicator%json%add(http_communicator%job_json, rejected_cls2D)
                            if(allocated(accepted_cls_ids) .and. allocated(rejected_cls_ids)) then
                                do i=0, jpg_ntiles - 1
                                    if(any( accepted_cls_ids == i + 1)) then
                                        call communicator_add_cls2D_accepted(selection_jpeg, i + 1, xtile * (100 / (jpg_nxtiles - 1)), ytile * (100 / (jpg_nytiles - 1)), 100 * jpg_nytiles, 100 * jpg_nxtiles)
                                    else if(any( rejected_cls_ids == i + 1)) then
                                        call communicator_add_cls2D_rejected(selection_jpeg, i + 1, xtile * (100 / (jpg_nxtiles - 1)), ytile * (100 / (jpg_nytiles - 1)), 100 * jpg_nytiles, 100 * jpg_nxtiles)                              
                                    endif
                                    xtile = xtile + 1
                                    if(xtile .eq. jpg_nxtiles) then
                                        xtile = 0
                                        ytile = ytile + 1
                                    endif
                                end do
                            endif
                        endif
                        call http_communicator%json%get(http_communicator%update_arguments, 'accepted_cls2D', accepted_cls2D, found)
                        if(found) then
                            call http_communicator%json%get(http_communicator%update_arguments, 'rejected_cls2D', rejected_cls2D, found)
                            if(found) then
                                call http_communicator%json%update(http_communicator%job_json, "user_input", .false., found)
                                l_wait_for_user = .false.
                            endif
                        endif
                    endif
                    ! first set uses cluster_cavgs
                    call is_set_processed(1)
                    call submit_cluster_cavgs
                endif
                ! make completed sets available to abinitio2D_stream
                call flag_complete_sets
            endif
            ! 2D analyses
            call analyze2D_new_chunks(projrecords)
            call sleep(WAITTIME)
            ! nice
            if(get_nchunks() > 0) then
                nice_communicator%stat_root%stage = "classifying chunks"
                nice_communicator%stat_root%user_input = .false.
            end if
            ! http stats send
            call http_communicator%send_jobstats()
        end do
        ! termination
        write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
        nice_communicator%stat_root%stage = "terminating"
        call nice_communicator%cycle()
        call terminate_chunks
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
        call simple_end('**** SIMPLE_STREAM_SIEVE_CAVGS NORMAL STOP ****')
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
                        call init_chunk_clustering(cline, spproj_glob)
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
                    call init_chunk_clustering(cline, spproj_glob)
                    call cline%delete('ncls')
                endif
                ! cleanup
                do iproj = 1,n_spprojs
                    call spprojs(iproj)%kill
                enddo
                deallocate(spprojs)
            end subroutine update_records_with_project

            subroutine generate_sets( chunks, sets )
                use simple_euclid_sigma2, only: average_sigma2_groups
                class(projs_list), intent(in)    :: chunks
                class(projs_list), intent(inout) :: sets
                type(sp_project)                       :: spproj
                character(len=LONGSTRLEN), allocatable :: starfiles(:)
                character(len=:),          allocatable :: tmpl
                integer :: navail_chunks, n, iset, i, ic, ic_start, ic_end
                navail_chunks = chunks%n - sets%n * params%nchunksperset
                n = floor(real(navail_chunks) / real(params%nchunksperset))
                if( n < 1 )return
                do iset = 1,n
                    ! merge chunks project into designated folder
                    ic_start = sets%n*params%nchunksperset + 1
                    ic_end   = ic_start + params%nchunksperset - 1
                    tmpl     = trim(DIR_SET)//int2str(sets%n+1)
                    call simple_mkdir(tmpl)
                    call merge_chunks(chunks%projfiles(ic_start:ic_end), tmpl, spproj, projname_out=tmpl)
                    ! average and stash sigma2
                    allocate(starfiles(params%nchunksperset))
                    do i = 1,params%nchunksperset
                        ic = ic_start + i - 1
                        starfiles(i) = trim(SIGMAS_DIR)//'/chunk_'//int2str(chunks%ids(ic))//trim(STAR_EXT)
                    enddo
                    call average_sigma2_groups(tmpl//'/'//tmpl//trim(STAR_EXT), starfiles)
                    deallocate(starfiles)
                    ! update global list and increment sets%n
                    call sets%append(tmpl//'/'//tmpl//trim(METADATA_EXT), sets%n+1, .false.)
                    ! remove imported chunk
                    if( trim(params%remove_chunks).eq.'yes' )then
                        do ic = ic_start,ic_end
                            call simple_rmdir(stemname(chunks%projfiles(ic)))
                        enddo
                    endif
                enddo
                call spproj%kill
            end subroutine generate_sets

            subroutine submit_cluster_cavgs
                type(cmdline)              :: cline_cluster_cavgs
                character(len=XLONGSTRLEN) :: cwd
                if( setslist%n < 1 )        return  ! no sets generated yet
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

            ! make completed project files visible to the watcher of the next application
            subroutine flag_complete_sets
                type(sp_project)              :: spproj_imported
                character(len=:), allocatable :: destination
                integer :: iset, n_state_nonzero
                do iset = 1,setslist%n
                    if( setslist%imported(iset) ) cycle
                    if( setslist%processed(iset) )then
                        destination = trim(DIR_STREAM_COMPLETED)//trim(DIR_SET)//int2str(iset)//trim(METADATA_EXT)
                        call simple_copy_file(setslist%projfiles(iset), destination)
                        setslist%imported(iset) = .true.
                        ! update particle counts
                        call spproj_imported%read_segment("ptcl2D", destination)
                        n_state_nonzero = spproj_imported%os_ptcl2D%count_state_gt_zero()
                        n_accepted = n_accepted + n_state_nonzero
                        n_rejected = n_rejected + spproj_imported%os_ptcl2D%get_noris() - n_state_nonzero 
                        call spproj_imported%kill()
                    endif
                enddo
            end subroutine flag_complete_sets
            
            subroutine communicator_init()
                call http_communicator%json%add(http_communicator%job_json, "stage",                   "initialising")
                call http_communicator%json%add(http_communicator%job_json, "particles_imported ",     0)
                call http_communicator%json%add(http_communicator%job_json, "particles_accepted",      0)
                call http_communicator%json%add(http_communicator%job_json, "particles_rejected",      0)
                call http_communicator%json%add(http_communicator%job_json, "user_input",              .false.)
                call http_communicator%json%add(http_communicator%job_json, "last_particles_imported", "")
                call http_communicator%json%create_array(accepted_cls2D, "accepted_cls2D")
                call http_communicator%json%add(http_communicator%job_json, accepted_cls2D)
                call http_communicator%json%create_array(rejected_cls2D, "rejected_cls2D")
                call http_communicator%json%add(http_communicator%job_json, rejected_cls2D)
            end subroutine communicator_init

            ! Remove previous files from folder to restart
            subroutine cleanup4restart
                character(len=STDLEN), allocatable :: folders(:)
                integer :: i
                if( .not.cline%defined('dir_exec') )then
                    ! nothing to do
                else
                    if( .not.file_exists(cline%get_carg('dir_exec')) )then
                        THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
                    endif
                    call cline%delete('dir_exec')
                    call del_file(TERM_STREAM)
                    call del_file(USER_PARAMS2D)
                    call simple_rmdir(SIGMAS_DIR)
                    call simple_rmdir(DIR_STREAM_COMPLETED)
                    folders = simple_list_dirs('.')
                    if( allocated(folders) )then
                        do i = 1,size(folders)
                            if( str_has_substr(folders(i),trim(DIR_CHUNK)).or.&
                                &str_has_substr(folders(i),trim(DIR_SET)) )then
                                call simple_rmdir(folders(i))
                            endif
                        enddo
                    endif
                endif
            end subroutine cleanup4restart

            subroutine generate_selection_jpeg()
                type(sp_project)              :: set1_proj
                integer                       :: ncls, icls
                real                          :: smpd
                call set1_proj%read_segment('cls2D', setslist%projfiles(i))
                call set1_proj%read_segment('out',   setslist%projfiles(i))
                call set1_proj%get_cavgs_stk(selection_jpeg, ncls, smpd)
                call mrc2jpeg_tiled(selection_jpeg, swap_suffix(selection_jpeg, JPG_EXT, params%ext), ntiles=jpg_ntiles, n_xtiles=jpg_nxtiles, n_ytiles=jpg_nytiles)
                allocate(accepted_cls_ids(0))
                allocate(rejected_cls_ids(0))
                do icls=1, set1_proj%os_cls2D%get_noris()
                    if(set1_proj%os_cls2D%isthere(icls, 'accept')) then
                        if(set1_proj%os_cls2D%get(icls, 'accept') .gt. 0.0) then
                            accepted_cls_ids = [accepted_cls_ids, icls]
                        else
                            rejected_cls_ids = [rejected_cls_ids, icls]
                        endif
                    endif
                enddo
                call set1_proj%kill()
                selection_jpeg_created = .true.
                selection_jpeg = swap_suffix(selection_jpeg, JPG_EXT, params%ext)
            end subroutine generate_selection_jpeg

            subroutine communicator_add_cls2D_accepted(path, idx, spritex, spritey, spriteh, spritew)
                character(*),     intent(in)  :: path
                integer,          intent(in)  :: spritex, spritey, spriteh, spritew, idx
                type(json_value), pointer     :: template
                call http_communicator%json%create_object(template, "")
                call http_communicator%json%add(template, "path",    path)
                call http_communicator%json%add(template, "spritex", spritex)
                call http_communicator%json%add(template, "spritey", spritey)
                call http_communicator%json%add(template, "spriteh", spriteh)
                call http_communicator%json%add(template, "spritew", spritew)
                call http_communicator%json%add(template, "idx",     idx)
                call http_communicator%json%add(accepted_cls2D,      template)
            end subroutine communicator_add_cls2D_accepted

            subroutine communicator_add_cls2D_rejected(path, idx, spritex, spritey, spriteh, spritew)
                character(*),     intent(in)  :: path
                integer,          intent(in)  :: spritex, spritey, spriteh, spritew, idx
                type(json_value), pointer     :: template
                call http_communicator%json%create_object(template, "")
                call http_communicator%json%add(template, "path",    path)
                call http_communicator%json%add(template, "spritex", spritex)
                call http_communicator%json%add(template, "spritey", spritey)
                call http_communicator%json%add(template, "spriteh", spriteh)
                call http_communicator%json%add(template, "spritew", spritew)
                call http_communicator%json%add(template, "idx",     idx)
                call http_communicator%json%add(rejected_cls2D,      template)
            end subroutine communicator_add_cls2D_rejected

    end subroutine exec_sieve_cavgs

    ! Manages Global 2D Clustering
    ! TODO: handling of un-classified particles
    subroutine exec_stream_abinitio2D( self, cline )
        class(commander_stream_abinitio2D), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        character(len=STDLEN),     parameter   :: micsspproj_fname = './streamdata.simple'
        type(parameters)                       :: params
        type(simple_nice_communicator)         :: nice_communicator
        type(projs_list)                       :: setslist
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj_glob
        character(len=LONGSTRLEN), allocatable :: projects(:)
        integer(kind=dp) :: time_last_import, time_last_iter
        integer :: i, iter, nprojects, nimported, nptcls_glob, nsets_imported, pool_iter, iter_last_import
        logical :: l_pause, l_params_updated
        call cline%set('oritype',      'mic')
        call cline%set('mkdir',        'yes')
        call cline%set('autoscale',    'yes')
        call cline%set('reject_mics',  'no')
        call cline%set('reject_cls',   'no') ! refers to previous implementation
        call cline%set('kweight_pool', 'default')
        call cline%set('wiener',       'full')
        call cline%set('refine',       'snhc_smpl')
        call cline%set('ml_reg',       'no')
        call cline%set('objfun',       'euclid')
        call cline%set('sigma_est',    'global')
        call cline%set('cls_init',     'rand')
        ! if( .not. cline%defined('walltime')  ) call cline%set('walltime',  29*60) ! 29 minutes
        if( .not. cline%defined('dynreslim') ) call cline%set('dynreslim', 'yes')
        if( .not.cline%defined('center')     ) call cline%set('center',    'yes')
        if( .not.cline%defined('algorithm')  ) call cline%set('algorithm', 'abinitio2D')
        if( .not.cline%defined('ncls')       ) call cline%set('ncls',       200)
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ! restart
        call cleanup4restart
        ! master parameters
        call cline%set('numlen', 5)
        call cline%set('stream','yes')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        ! initialise progress monitor
        call progressfile_init()
        ! master project file
        call spproj_glob%read( params%projfile )
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_cluster2D must start from an empty project (eg from root project folder)')
        ! project watcher
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true., nretries=10)
        ! Infinite loop
        iter             = 0
        nprojects        = 0        ! # of projects per iteration
        nimported        = 0        ! # of sets per iteration
        nptcls_glob      = 0        ! global # of particles present in the pool
        l_pause          = .false.  ! pause clustering
        nsets_imported   = 0        ! Global # of sets imported
        pool_iter        = 0
        time_last_import = huge(time_last_import)
        iter_last_import = -1       ! Pool iteration # when set(s) last imported
        do
            ! termination
            if( file_exists(trim(TERM_STREAM)) .or. nice_communicator%exit ) exit
            iter = iter + 1
            ! detection of new projects
            call project_buff%watch(nprojects, projects)
            if( nprojects > 0 )then
                ! memoize detected projects
                call project_buff%add2history(projects)
                do i = 1,nprojects
                    call setslist%append(projects(i), setslist%n+1, .true.)
                enddo
            endif
            ! check on progress, updates particles & alignement parameters
            ! TODO: class remapping
            if( l_pause )then
                call generate_pool_stats
                ! TODO Flush raw particles
            else
                ! progress
                call update_pool_status
                ! updates particles if iteration completes
                call update_pool
            endif
            ! Import new particles & clustering init
            call import_sets_into_pool( nimported )
            if( nimported > 0 )then
                time_last_import = time8()
                iter_last_import = get_pool_iter()
                call unpause_pool
            endif
            nsets_imported = count(setslist%imported)
            call update_user_params2D(cline, l_params_updated, nice_communicator%update_arguments)
            if( l_params_updated ) call unpause_pool
            ! pause?
            if( (pool_iter >= iter_last_import+PAUSE_NITERS) .or.&
                & (time8()-time_last_import>PAUSE_TIMELIMIT) )then
                if( .not.l_pause )then
                    l_pause = is_pool_available()
                    if( l_pause ) write(logfhandle,'(A)')'>>> PAUSING 2D ANALYSIS'
                endif
            endif
            ! Performs clustering iteration
            if( l_pause )then
                ! skip iteration
            else
                ! optionally updates alignement parameters
                call update_pool_aln_params
                ! initiates new iteration
                pool_iter = get_pool_iter()
                call iterate_pool
                if( get_pool_iter() > pool_iter )then
                    nice_communicator%stat_root%user_input = .true.
                    time_last_iter = time8()
                endif
            endif
            ! Wait
            call sleep(WAITTIME)
        enddo
        ! Cleanup and final project
        call terminate_stream2D
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
        call nice_communicator%terminate()
        call simple_end('**** SIMPLE_STREAM_ABINITIO2D NORMAL STOP ****')
        contains

            ! resumes with new sets or analysis parameters have been updated
            subroutine unpause_pool
                if( l_pause ) write(logfhandle,'(A)')'>>> RESUMING 2D ANALYSIS'
                l_pause = .false.
            end subroutine unpause_pool

            ! imports new sets of pre-classified particles into the pool
            ! and initialize the clustering module
            subroutine import_sets_into_pool( nimported )
                use simple_euclid_sigma2, only: average_sigma2_groups, sigma2_star_from_iter
                integer,                   intent(out) :: nimported
                type(sp_project),          allocatable :: spprojs(:)
                character(len=LONGSTRLEN), allocatable :: sigmas(:)
                class(sp_project),             pointer :: pool
                integer :: nsets2import, iset, nptcls2import, nmics2import, nmics, nptcls
                integer :: i, fromp, fromp_prev, imic, ind, iptcl, jptcl, jmic, nptcls_sel
                nimported = 0
                if( setslist%n == 0 ) return
                if( nptcls_glob == 0 )then
                    ! first import
                else
                    ! at other times only import when the pool is free
                    if( .not.is_pool_available() ) return
                endif
                nsets2import = count(setslist%processed(:).and.(.not.setslist%imported(:)))
                if( nsets2import == 0 ) return
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
                do iset = 1,setslist%n
                    if( setslist%imported(iset) )       cycle   ! already imported
                    if( .not.setslist%processed(iset) ) cycle   ! not processed yet
                    if( setslist%busy(iset) )           cycle   ! being processed
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
                            ! new particle
                            call pool%os_ptcl2D%set(iptcl, 'updatecnt', 0)
                            call pool%os_ptcl2D%set(iptcl, 'frac',      0.)
                            call pool%os_ptcl2D%delete_2Dclustering(iptcl, keepshifts=.true.)
                        enddo
                        !$omp end parallel do
                        ind   = ind   + nptcls
                        fromp = fromp + nptcls
                    enddo
                    nptcls_sel = spprojs(iset)%os_ptcl2D%get_noris(consider_state=.true.)
                    ! display
                    write(logfhandle,'(A,I6,A,I6)')'>>> TRANSFERRED ',nptcls_sel,' PARTICLES FROM SET ',setslist%ids(iset)
                    call flush(logfhandle)
                    ! global list update
                    setslist%imported(iset) = .true.
                enddo
                nimported = nsets2import
                ! average all previously imported sigmas
                allocate(sigmas(count(setslist%imported(:))))
                i = 0
                do iset = 1,setslist%n
                    if( .not.setslist%imported(iset) ) cycle
                    i         = i+1
                    ind       = setslist%ids(iset)
                    sigmas(i) = trim(params%dir_target)//'/set_'//int2str(ind)//'/set_'//int2str(ind)//trim(STAR_EXT)
                enddo
                call average_sigma2_groups(sigma2_star_from_iter(get_pool_iter()+1), sigmas)
                deallocate(sigmas)
                ! initialize fundamental parameters & clustering
                if( nptcls_glob == 0 )then
                    params%smpd = pool%get_smpd()
                    params%box  = pool%get_box()
                    call init_pool_clustering(cline, spproj_glob, micsspproj_fname, reference_generation=.false.)
                endif
                ! global count
                nptcls_glob = nptcls_glob + nptcls_sel
                ! cleanup
                do iset = 1,setslist%n
                    call spprojs(iset)%kill
                enddo
                deallocate(spprojs)
                nullify(pool)
            end subroutine import_sets_into_pool

            ! Remove previous files from folder to restart
            subroutine cleanup4restart
                if( .not.cline%defined('dir_exec') )then
                    ! nothing to do
                else
                    if( .not.file_exists(cline%get_carg('dir_exec')) )then
                        THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
                    endif
                    call cline%delete('dir_exec')
                    call del_file(micspproj_fname)
                    call cleanup_root_folder
                endif
            end subroutine cleanup4restart


    end subroutine exec_stream_abinitio2D

end module simple_commander_stream2D
