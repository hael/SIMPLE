! concrete commander: cluster2D_stream for streaming 2D alignment and clustering of single-particle images
module simple_commanders_stream2D
include 'simple_lib.f08'
use simple_binimage,           only: binimage
use simple_cmdline,            only: cmdline
use simple_commander_base,     only: commander_base
use simple_ctf,                only: ctf
use simple_ctf_estimate_iter,  only: ctf_estimate_iter
use simple_guistats,           only: guistats
use simple_parameters,         only: parameters
use simple_particle_extractor, only: ptcl_extractor
use simple_sp_project,         only: sp_project
use simple_commanders_abinitio2D
use simple_commanders_cluster2D_stream
use simple_commanders_imgproc
use simple_commanders_preprocess
use simple_commanders_project
use simple_gui_utils
use simple_micproc
use simple_moviewatcher
use simple_nice
use simple_progress
use simple_qsys_env
use simple_qsys_funs
use simple_segmentation
use simple_stack_io
use simple_strategy2D_utils
use simple_stream_communicator
use simple_stream_utils
implicit none

public :: commander_stream_sieve_cavgs, commander_stream_abinitio2D, stream_test_sieve_cavgs
private
#include "simple_local_flags.inc"

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
character(len=STDLEN), parameter :: micspproj_fname      = './streamdata.simple'
character(len=STDLEN), parameter :: REJECTED_CLS_STACK   = './rejected_cls.mrc'
integer,               parameter :: LONGTIME             = 60    ! time lag after which a movie/project is processed
integer,               parameter :: WAITTIME             = 10    ! movie folder watched every WAITTIME seconds
integer(kind=dp),      parameter :: FLUSH_TIMELIMIT      = 900   ! time (secs) after which leftover particles join the pool IF the 2D analysis is paused
integer,               parameter :: PAUSE_NITERS         = 5     ! # of iterations after which 2D analysis is paused
integer,               parameter :: PAUSE_TIMELIMIT      = 600   ! time (secs) after which 2D analysis is paused

contains

    ! Manages individual chunks/sets classification, matching & rejection
    ! TODO: handling of un-classified particles
    subroutine exec_sieve_cavgs( self, cline )
        class(commander_stream_sieve_cavgs), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(projrecord),          allocatable :: projrecords(:)
        type(parameters)                       :: params
        type(qsys_env)                         :: qenv
        type(stream_http_communicator)         :: http_communicator
        type(projs_list)                       :: chunkslist, setslist
        type(oris)                             :: moldiamori, chunksizeori
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj_glob
        type(json_value),          pointer     :: accepted_cls2D, rejected_cls2D, latest_accepted_cls2D, latest_rejected_cls2D
        character(len=LONGSTRLEN), allocatable :: projects(:)
        character(len=STDLEN)                  :: chunk_part_env
        character(len=:),          allocatable :: selection_jpeg
        integer,                   allocatable :: accepted_cls_ids(:), rejected_cls_ids(:), jpg_cls_map(:)
        real,                      allocatable :: cls_res(:), cls_pop(:)
        real             :: mskdiam
        integer(kind=dp) :: time_last_import
        integer          :: nchunks_glob, nchunks_imported, nprojects, iter, i, envlen
        integer          :: n_imported, n_imported_prev, nptcls_glob, n_failed_jobs
        integer          :: n_accepted, n_rejected, jpg_ntiles, jpg_nxtiles, jpg_nytiles, xtile, ytile
        integer          :: latest_processed_set, latest_displayed_set
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
        if( .not.cline%defined('nmics')        ) call cline%set('nmics',        100)
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
        ! http communicator init
        call http_communicator%create(params%niceprocid, params%niceserver, "sieve_cavgs")
        call communicator_init()
        call http_communicator%send_jobstats()
        ! wait if dir_target doesn't exist yet
        call wait_for_folder(http_communicator, params%dir_target, '**** SIMPLE_STREAM_SIEVE_CAVGS USER STOP ****')
        call wait_for_folder(http_communicator, trim(params%dir_target)//'/spprojs', '**** SIMPLE_STREAM_SIEVE_CAVGS USER STOP ****')
        call wait_for_folder(http_communicator, trim(params%dir_target)//'/spprojs_completed', '**** SIMPLE_STREAM_SIEVE_CAVGS USER STOP ****')
        ! mskdiam
        if( .not. cline%defined('mskdiam') )then
            ! obtain mskdiam from moldiam ori file written by pick_extract (generated by generate_pickrefs)
            write(logfhandle,'(A,F8.2)')'>>> WAITING UP TO 5 MINUTES FOR '//trim(STREAM_MOLDIAM)
            do i=1, 30
                if(file_exists(trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM))) exit
                call sleep(10)
                call http_communicator%send_jobstats()
                if( http_communicator%exit )then
                    ! termination
                    write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                    call http_communicator%term()
                    call simple_end('**** SIMPLE_STREAM_PICK_EXTRACT USER STOP ****')
                    call EXIT(0)
                endif
            end do
            if( .not. file_exists(trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM))) THROW_HARD('either mskdiam must be given or '// trim(STREAM_MOLDIAM) // ' exists in target_dir')
            ! read mskdiam from file
            call moldiamori%new(1, .false.)
            call moldiamori%read( trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM) )
            if( .not. moldiamori%isthere(1, "mskdiam") ) THROW_HARD( 'mskdiam missing from ' // trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM) )
            mskdiam = moldiamori%get(1, "mskdiam")
            ! write a copy for stream 2D downstream
            if(file_exists(trim(STREAM_MOLDIAM))) call del_file(trim(STREAM_MOLDIAM))
            call moldiamori%write(1, trim(STREAM_MOLDIAM))
            call moldiamori%kill
            params%mskdiam = mskdiam
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
        if(params%nptcls_per_cls == 0) write(logfhandle, '(A,I6,A)')'>>> # PARTICLES PER CLASS WILL BE AUTO DETERMINED AFTER', params%nmics ,' IMPORTED MICROGRAPHS'
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
        latest_displayed_set = 0   !
        time_last_import = huge(time_last_import)   ! used for flushing unprocessed particles
        do
            ! termination
            if( file_exists(trim(TERM_STREAM)) .or. http_communicator%exit ) exit
            if( http_communicator%stop .or. test_repick() ) then
                if(test_repick()) call write_repick_refs("../repick_refs.mrc")
                ! termination
                write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                call spproj_glob%kill
                call qsys_cleanup
            !    call nice_communicator%terminate(stop=.true.)
                call simple_end('**** SIMPLE_STREAM_SIEVE_CAVGS USER STOP ****')
                call EXIT(0)
            endif
            iter = iter + 1
            ! http stats
            call http_communicator%json%update(http_communicator%job_json, "stage",               "finding and sieving particles", found)    
            call http_communicator%json%update(http_communicator%job_json, "particles_accepted",  n_accepted,                      found)
            call http_communicator%json%update(http_communicator%job_json, "particles_rejected",  n_rejected,                      found)
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
                write(logfhandle,'(A,I6,I8)') '>>> # MICROGRAPHS / PARTICLES IMPORTED : ', n_imported, nptcls_glob
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
            !call update_user_params2D(cline, l_params_updated, nice_communicator%update_arguments)
            call update_user_params2D(cline, l_params_updated)
            call update_chunks
            call memoize_chunks(chunkslist, nchunks_imported)
            !call update_user_params2D(cline, l_params_updated, nice_communicator%update_arguments)
            call update_user_params2D(cline, l_params_updated)
            if( nchunks_imported > 0 )then
                nchunks_glob = nchunks_glob + nchunks_imported
                ! build sets
                call generate_sets(chunkslist, setslist)
            endif
            ! Sets analysis section
            if( setslist%n > 0 )then
                if( setslist%processed(1) .and. (setslist%n > 1) .and. (.not.l_wait_for_user)) then
                    ! all sets but the first employ match_cavgs
                    latest_processed_set = 0
                    do i = 2,setslist%n
                        call is_set_processed(i)
                        if( setslist%processed(i) ) latest_processed_set = i
                    enddo
                    call submit_match_cavgs
                    ! http stats
                    call http_communicator%json%remove(accepted_cls2D, destroy=.true.)
                    call http_communicator%json%remove(rejected_cls2D, destroy=.true.)
                    if(latest_processed_set > 0 .and. latest_displayed_set .ne. latest_processed_set) then
                        latest_displayed_set = latest_processed_set
                        call generate_set_jpeg()
                        xtile = 0
                        ytile = 0
                        call http_communicator%json%remove(latest_accepted_cls2D, destroy=.true.)
                        call http_communicator%json%create_array(latest_accepted_cls2D, "latest_accepted_cls2D")
                        call http_communicator%json%add(http_communicator%job_json, latest_accepted_cls2D)
                        call http_communicator%json%remove(latest_rejected_cls2D, destroy=.true.)
                        call http_communicator%json%create_array(latest_rejected_cls2D, "latest_rejected_cls2D")
                        call http_communicator%json%add(http_communicator%job_json, latest_rejected_cls2D)
                        if(allocated(accepted_cls_ids) .and. allocated(rejected_cls_ids)) then
                            do i=0, size(jpg_cls_map) - 1
                                if(any( accepted_cls_ids == jpg_cls_map(i + 1))) then
                                    call communicator_add_cls2D_accepted(selection_jpeg,&
                                        &i + 1,&
                                        &xtile * (100.0 / (jpg_nxtiles - 1)),&
                                        &ytile * (100.0 / (jpg_nytiles - 1)),&
                                        &100 * jpg_nytiles,&
                                        &100 * jpg_nxtiles,&
                                        &latest=.true.,&
                                        &res=cls_res(jpg_cls_map(i+1)),&
                                        &pop=nint(cls_pop(jpg_cls_map(i+1))))
                                else if(any( rejected_cls_ids == jpg_cls_map(i + 1))) then
                                    call communicator_add_cls2D_rejected(selection_jpeg,&
                                        &i + 1,&
                                        &xtile * (100.0 / (jpg_nxtiles - 1)),&
                                        &ytile * (100.0 / (jpg_nytiles - 1)),&
                                        &100 * jpg_nytiles,&
                                        &100 * jpg_nxtiles,&
                                        &latest=.true.,&
                                        &res=cls_res(jpg_cls_map(i+1)),&
                                        &pop=nint(cls_pop(jpg_cls_map(i+1))))
                                endif
                                xtile = xtile + 1
                                if(xtile .eq. jpg_nxtiles) then
                                    xtile = 0
                                    ytile = ytile + 1
                                endif
                            end do
                        endif
                    endif
                else
                    ! first set uses cluster_cavgs
                    call is_set_processed(1)
                    call submit_cluster_cavgs
                    ! interactive selection
                    if( l_wait_for_user .and. setslist%processed(1) ) then
                        call http_communicator%json%update(http_communicator%job_json, "user_input", .true., found)
                        call http_communicator%json%update(http_communicator%job_json, "stage", "waiting for user selection", found)
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
                                do i=0, size(jpg_cls_map) - 1
                                    if(any( accepted_cls_ids == jpg_cls_map(i + 1))) then
                                        call communicator_add_cls2D_accepted(selection_jpeg,&
                                            &jpg_cls_map(i + 1),&
                                            &xtile * (100.0 / (jpg_nxtiles - 1)),&
                                            &ytile * (100.0 / (jpg_nytiles - 1)),&
                                            &100 * jpg_nytiles,&
                                            &100 * jpg_nxtiles,&
                                            &res=cls_res(jpg_cls_map(i+1)),&
                                            &pop=nint(cls_pop(jpg_cls_map(i+1))))
                                    else if(any( rejected_cls_ids == jpg_cls_map(i + 1))) then
                                        call communicator_add_cls2D_rejected(selection_jpeg,&
                                            &jpg_cls_map(i + 1),&
                                            &xtile * (100.0 / (jpg_nxtiles - 1)),&
                                            &ytile * (100.0 / (jpg_nytiles - 1)),&
                                            &100 * jpg_nytiles,&
                                            &100 * jpg_nxtiles,&
                                            &res=cls_res(jpg_cls_map(i+1)),&
                                            &pop=nint(cls_pop(jpg_cls_map(i+1))))
                                    endif
                                    xtile = xtile + 1
                                    if(xtile .eq. jpg_nxtiles) then
                                        xtile = 0
                                        ytile = ytile + 1
                                    endif
                                end do
                            endif
                        endif
                    endif
                endif
            endif
            if( l_wait_for_user )then
                ! nothing for abinitio2D_stream until the first set has been selected
                if(allocated(setslist%processed)) then
                    if(size(setslist%processed) .gt. 0 .and. setslist%processed(1)) then
                        call http_communicator%json%get(http_communicator%update_arguments, 'accepted_cls2D', accepted_cls_ids, found) ! accepted_cls2D now contains user selection
                        if(found) then
                            call http_communicator%json%get(http_communicator%update_arguments, 'rejected_cls2D', rejected_cls_ids, found) ! rejected_cls2D now contains user selection
                            if(found) then
                                call http_communicator%json%update(http_communicator%job_json, "user_input", .false., found)
                                ! apply interactive selection
                                call report_interactive_selection( rejected_cls_ids )
                                write(logfhandle,'(A)') '>>> RECEIVED USER SELECTIONS'
                                ! http stats
                                call http_communicator%json%remove(accepted_cls2D, destroy=.true.)
                                call http_communicator%json%remove(rejected_cls2D, destroy=.true.)
                                xtile = 0
                                ytile = 0
                                call http_communicator%json%remove(latest_accepted_cls2D, destroy=.true.)
                                call http_communicator%json%create_array(latest_accepted_cls2D, "latest_accepted_cls2D")
                                call http_communicator%json%add(http_communicator%job_json, latest_accepted_cls2D)
                                call http_communicator%json%remove(latest_rejected_cls2D, destroy=.true.)
                                call http_communicator%json%create_array(latest_rejected_cls2D, "latest_rejected_cls2D")
                                call http_communicator%json%add(http_communicator%job_json, latest_rejected_cls2D)
                                if(allocated(accepted_cls_ids) .and. allocated(rejected_cls_ids)) then
                                    do i=0, size(jpg_cls_map) - 1
                                        if(any( accepted_cls_ids == jpg_cls_map(i + 1))) then
                                            call communicator_add_cls2D_accepted(selection_jpeg,&
                                                &jpg_cls_map(i + 1),&
                                                &xtile * (100.0 / (jpg_nxtiles - 1)),&
                                                &ytile * (100.0 / (jpg_nytiles - 1)),&
                                                &100 * jpg_nytiles,&
                                                &100 * jpg_nxtiles,&
                                                &latest=.true.,&
                                                &res=cls_res(i+1),&
                                                &pop=nint(cls_pop(i+1)))
                                        else if(any( rejected_cls_ids == jpg_cls_map(i + 1))) then
                                            call communicator_add_cls2D_rejected(selection_jpeg,&
                                                &jpg_cls_map(i + 1),&
                                                &xtile * (100.0 / (jpg_nxtiles - 1)),&
                                                &ytile * (100.0 / (jpg_nytiles - 1)),&
                                                &100 * jpg_nytiles,&
                                                &100 * jpg_nxtiles,&
                                                &latest=.true.,&
                                                &res=cls_res(i+1),&
                                                &pop=nint(cls_pop(i+1)))                      
                                        endif
                                        xtile = xtile + 1  
                                        if(xtile .eq. jpg_nxtiles) then
                                            xtile = 0
                                            ytile = ytile + 1
                                        endif
                                    end do
                                endif
                                l_wait_for_user = .false.
                            endif
                        endif
                    endif
                endif
            else
                ! make completed sets available to abinitio2D_stream
                call flag_complete_sets
            endif
            ! 2D analyses
            call analyze2D_new_chunks(projrecords)
            call sleep(WAITTIME)
            ! http stats send
            call http_communicator%send_jobstats()
        end do
        ! termination
        write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
        call terminate_chunks
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
        call http_communicator%term()
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
                n_spprojs = size(projectnames)
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
                    if(size(projrecords) .gt. params%nmics) then
                        avgmicptcls    = nptcls_glob / size(projrecords)
                        avgmicptcls    = ceiling(avgmicptcls / 10) * 10.0
                        nptcls_per_cls = ceiling(params%nmics * avgmicptcls / params%ncls)
                        write(logfhandle,'(A,I6)')   '>>> AVERAGE # PARTICLES PER MICROGRAPH : ', int(avgmicptcls)
                        write(logfhandle,'(A,I6,A)') '>>> USING ', int(nptcls_per_cls), ' PARTICLES PER CLASS'
                        params%nptcls_per_cls = int(nptcls_per_cls)
                        call cline%set('nptcls_per_cls', nptcls_per_cls)
                        params%smpd = spprojs(first)%os_mic%get(1,'smpd')
                        call spprojs(first)%read_segment('stk', trim(projectnames(first)))
                        params%box  = nint(spprojs(first)%os_stk%get(1,'box'))
                        if(params%mskdiam .eq. 0.0) then
                            params%mskdiam = 0.9 * ceiling(params%box * spprojs(first)%os_stk%get(1,'smpd'))
                            call cline%set('mskdiam', params%mskdiam)
                            write(logfhandle,'(A,F8.2)')'>>> MASK DIAMETER SET TO', params%mskdiam
                        endif
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
                    if(params%mskdiam .eq. 0.0) then
                        params%mskdiam = 0.9 * ceiling(params%box * spprojs(first)%os_stk%get(1,'smpd'))
                        call cline%set('mskdiam', params%mskdiam)
                        write(logfhandle,'(A,F8.2)')'>>> MASK DIAMETER SET TO', params%mskdiam
                    endif
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

            ! apply user-inputted selection on the first set
            subroutine report_interactive_selection( cls2reject )
                integer, allocatable, intent(in) :: cls2reject(:)
                integer, allocatable             :: states(:)
                type(sp_project)                 :: spproj
                integer                          :: ncls2reject, i
                if( .not.allocated(cls2reject) ) return
                ncls2reject = size(cls2reject)
                ! read all fields
                call spproj%read(setslist%projfiles(1))
                ! selection
                allocate(states(spproj%os_cls2D%get_noris()),source=1)
                ! state=0 for all rejected
                do i=1, ncls2reject
                    states(cls2reject(i)) = 0
                enddo
                ! maintain state=0 for junk 
                do i=1, spproj%os_cls2D%get_noris()
                    if(.not. spproj%os_cls2D%isthere(i, 'cluster'))  states(i) = 0
                enddo
                call spproj%map_cavgs_selection(states)
                call spproj%write(setslist%projfiles(1))
                call spproj%kill
                deallocate(states)
            end subroutine report_interactive_selection

            ! make completed project files visible to the watcher of the next application
            subroutine flag_complete_sets
                use simple_image, only:image
                type(sp_project)              :: spproj
                type(image)                   :: img
                character(len=:), allocatable :: destination, stk
                real    :: smpd
                integer :: ldim(3), icls, iset, n_state_nonzero, nimgs, ncls
                do iset = 1,setslist%n
                    if( setslist%imported(iset) ) cycle
                    if( setslist%processed(iset) )then
                        destination = trim(DIR_STREAM_COMPLETED)//trim(DIR_SET)//int2str(iset)//trim(METADATA_EXT)
                        call simple_rename(setslist%projfiles(iset), destination)
                        setslist%projfiles(iset) = destination ! relocation
                        setslist%imported(iset)  = .true.
                        write(logfhandle,'(A,I3)')'>>> COMPLETED SET ',setslist%ids(iset)
                        ! update particle counts
                        call spproj%read_segment('ptcl2D', destination)
                        n_state_nonzero = spproj%os_ptcl2D%count_state_gt_zero()
                        n_accepted      = n_accepted + n_state_nonzero
                        n_rejected      = n_rejected + spproj%os_ptcl2D%get_noris() - n_state_nonzero
                        ! updates stack of rejected classes
                        call spproj%read_segment('cls2D', destination)
                        call spproj%read_segment('out', destination)
                        call spproj%get_cavgs_stk(stk, ncls, smpd)
                        nimgs = 0
                        if( file_exists(REJECTED_CLS_STACK) )then
                            call find_ldim_nptcls(REJECTED_CLS_STACK, ldim, nimgs)
                        else
                            call find_ldim_nptcls(stk, ldim, ncls)
                        endif
                        if( .not.img%exists() )then
                            ldim(3) = 1
                            call img%new(ldim,smpd)
                        endif
                        do icls = 1,ncls
                            if( spproj%os_cls2D%get_state(icls) == 0 )then
                                nimgs = nimgs+1
                                call img%read(stk,icls)
                                call img%write(REJECTED_CLS_STACK, nimgs)
                            endif
                        enddo
                        call spproj%kill
                    endif
                enddo
                call img%kill
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
                call http_communicator%json%create_array(latest_accepted_cls2D, "latest_accepted_cls2D")
                call http_communicator%json%add(http_communicator%job_json, latest_accepted_cls2D)
                call http_communicator%json%create_array(latest_rejected_cls2D, "latest_rejected_cls2D")
                call http_communicator%json%add(http_communicator%job_json, latest_rejected_cls2D)
            end subroutine communicator_init

            ! Remove previous files from folder to restart
            subroutine cleanup4restart
                character(len=STDLEN), allocatable :: folders(:)
                character(len=XLONGSTRLEN)         :: cwd_restart
                logical :: l_restart
                integer :: i
                call simple_getcwd(cwd_restart)
                l_restart = .false.
                if(cline%defined('outdir') .and. dir_exists(trim(cline%get_carg('outdir')))) then
                    l_restart = .true.
                    call chdir(trim(cline%get_carg('outdir')))
                endif
                if(cline%defined('dir_exec')) then
                    if( .not.file_exists(cline%get_carg('dir_exec')) )then
                        THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
                    endif
                    l_restart = .true.
                endif
                if( l_restart ) then
                    write(logfhandle, *) ">>> RESTARTING EXISTING JOB", trim(cwd_restart)
                    if(cline%defined('dir_exec')) call cline%delete('dir_exec')
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
                call chdir(trim(cwd_restart))
            end subroutine cleanup4restart

            subroutine generate_selection_jpeg()
                type(sp_project)              :: set1_proj
                integer                       :: ncls, icls
                real                          :: smpd
                call set1_proj%read_segment('cls2D', setslist%projfiles(1))
                call set1_proj%read_segment('out',   setslist%projfiles(1))
                call set1_proj%get_cavgs_stk(selection_jpeg, ncls, smpd)
                call mrc2jpeg_tiled(selection_jpeg, swap_suffix(selection_jpeg, JPG_EXT, params%ext), ntiles=jpg_ntiles, n_xtiles=jpg_nxtiles, n_ytiles=jpg_nytiles)
                if(allocated(cls_res))     deallocate(cls_res)
                if(allocated(cls_pop))     deallocate(cls_pop)
                if(allocated(jpg_cls_map)) deallocate(jpg_cls_map)
                allocate(jpg_cls_map(0))
                allocate(accepted_cls_ids(0))
                allocate(rejected_cls_ids(0))
                do icls=1, set1_proj%os_cls2D%get_noris()
                    if(set1_proj%os_cls2D%isthere(icls, 'accept') .and. set1_proj%os_cls2D%isthere(icls, 'res') .and. set1_proj%os_cls2D%isthere(icls, 'pop')) then
                        if(set1_proj%os_cls2D%get(icls, 'accept') .gt. 0.0) then
                            accepted_cls_ids = [accepted_cls_ids, icls]
                        else
                            rejected_cls_ids = [rejected_cls_ids, icls]
                        endif
                    endif
                    jpg_cls_map = [jpg_cls_map, icls]
                enddo
                cls_res = set1_proj%os_cls2D%get_all("res")
                cls_pop = set1_proj%os_cls2D%get_all("pop")
                call set1_proj%kill()
                selection_jpeg_created = .true.
                selection_jpeg = swap_suffix(selection_jpeg, JPG_EXT, params%ext)
            end subroutine generate_selection_jpeg

            subroutine generate_set_jpeg()
                type(sp_project)              :: set_proj
                integer                       :: ncls, icls
                real                          :: smpd
                call set_proj%read_segment('cls2D', setslist%projfiles(latest_processed_set))
                call set_proj%read_segment('out',   setslist%projfiles(latest_processed_set))
                call set_proj%get_cavgs_stk(selection_jpeg, ncls, smpd)
                call mrc2jpeg_tiled(selection_jpeg, swap_suffix(selection_jpeg, JPG_EXT, params%ext), ntiles=jpg_ntiles, n_xtiles=jpg_nxtiles, n_ytiles=jpg_nytiles)
                if(allocated(accepted_cls_ids)) deallocate(accepted_cls_ids)
                if(allocated(rejected_cls_ids)) deallocate(rejected_cls_ids)
                if(allocated(cls_res))          deallocate(cls_res)
                if(allocated(cls_pop))          deallocate(cls_pop)
                if(allocated(jpg_cls_map))      deallocate(jpg_cls_map)
                allocate(jpg_cls_map(0))
                allocate(accepted_cls_ids(0))
                allocate(rejected_cls_ids(0))
                do icls=1, set_proj%os_cls2D%get_noris()
                    if(set_proj%os_cls2D%isthere(icls, 'accept') .and. set_proj%os_cls2D%isthere(icls, 'res') .and. set_proj%os_cls2D%isthere(icls, 'pop')) then
                        if(set_proj%os_cls2D%get(icls, 'accept') .gt. 0.0) then
                            accepted_cls_ids = [accepted_cls_ids, icls]
                        else
                            rejected_cls_ids = [rejected_cls_ids, icls]
                        endif
                    endif
                    jpg_cls_map = [jpg_cls_map, icls]
                enddo
                cls_res = set_proj%os_cls2D%get_all("res")
                cls_pop = set_proj%os_cls2D%get_all("pop")
                call set_proj%kill()
                selection_jpeg = swap_suffix(selection_jpeg, JPG_EXT, params%ext)
            end subroutine generate_set_jpeg

            subroutine communicator_add_cls2D_accepted(path, idx, spritex, spritey, spriteh, spritew, latest, res, pop)
                character(*),      intent(in) :: path
                real,              intent(in) :: spritex, spritey
                integer,           intent(in) :: spriteh, spritew, idx
                logical, optional, intent(in) :: latest
                integer, optional, intent(in) :: pop
                real,    optional, intent(in) :: res
                type(json_value),  pointer    :: template
                logical                       :: l_latest = .false.
                if(present(latest)) l_latest = latest
                call http_communicator%json%create_object(template, "")
                call http_communicator%json%add(template, "path",    path)
                call http_communicator%json%add(template, "spritex", dble(spritex))
                call http_communicator%json%add(template, "spritey", dble(spritey))
                call http_communicator%json%add(template, "spriteh", spriteh)
                call http_communicator%json%add(template, "spritew", spritew)
                call http_communicator%json%add(template, "idx",     idx)
                if(present(res)) call http_communicator%json%add(template, "res",     dble(res))
                if(present(pop)) call http_communicator%json%add(template, "pop",     pop)
                if(l_latest) then
                    call http_communicator%json%add(latest_accepted_cls2D, template)
                else
                    call http_communicator%json%add(accepted_cls2D, template)
                endif
            end subroutine communicator_add_cls2D_accepted

            subroutine communicator_add_cls2D_rejected(path, idx, spritex, spritey, spriteh, spritew, latest, res, pop)
                character(*),      intent(in) :: path
                real,              intent(in) :: spritex, spritey
                integer,           intent(in) :: spriteh, spritew, idx
                logical, optional, intent(in) :: latest
                integer, optional, intent(in) :: pop
                real,    optional, intent(in) :: res
                type(json_value),  pointer    :: template
                logical                       :: l_latest = .false.
                if(present(latest)) l_latest = latest
                call http_communicator%json%create_object(template, "")
                call http_communicator%json%add(template, "path",    path)
                call http_communicator%json%add(template, "spritex", dble(spritex))
                call http_communicator%json%add(template, "spritey", dble(spritey))
                call http_communicator%json%add(template, "spriteh", spriteh)
                call http_communicator%json%add(template, "spritew", spritew)
                call http_communicator%json%add(template, "idx",     idx)
                if(present(res)) call http_communicator%json%add(template, "res",     dble(res))
                if(present(pop)) call http_communicator%json%add(template, "pop",     pop)
                if(l_latest) then
                    call http_communicator%json%add(latest_rejected_cls2D, template)
                else
                    call http_communicator%json%add(rejected_cls2D, template)
                endif
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
        type(stream_http_communicator)         :: http_communicator
        type(json_value),          pointer     :: latest_cls2D
        type(projs_list)                       :: setslist
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj_glob
        type(oris)                             :: moldiamori             
        character(kind=CK,len=:),  allocatable :: snapshot_filename
        character(len=LONGSTRLEN), allocatable :: projects(:)
        integer(kind=dp) :: time_last_import, time_last_iter
        integer :: i, iter, nprojects, nimported, nptcls_glob, nsets_imported, pool_iter, iter_last_import
        integer :: xtile, ytile, mskdiam_update, extra_pause_iters
        logical :: l_pause, l_params_updated, found
        real    :: mskdiam
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
        ! restart
        call cleanup4restart
        ! generate own project file if projfile isnt set
        if(cline%get_carg('projfile') .eq. '') then 
            call cline%set('projname', 'stream_abinitio2D')
            call cline%set('projfile', 'stream_abinitio2D.simple')
            call spproj_glob%update_projinfo(cline)
            call spproj_glob%update_compenv(cline)
            call spproj_glob%write
        endif
        ! master parameters
        call cline%set('numlen', 5)
        call cline%set('stream','yes')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! nice communicator init
        call nice_communicator%init(params%niceprocid,'')
        call nice_communicator%cycle()
        ! http communicator init
        call http_communicator%create(params%niceprocid, params%niceserver, "classification_2D")
        call communicator_init()
        call http_communicator%send_jobstats()
        ! wait if dir_target doesn't exist yet
        call wait_for_folder(http_communicator, params%dir_target, '**** SIMPLE_STREAM_ABINITIO2D NORMAL STOP ****')
        call wait_for_folder(http_communicator, trim(params%dir_target)//'/spprojs_completed', '**** SIMPLE_STREAM_ABINITIO2D NORMAL STOP ****')
        ! wait for and retrieve mskdiam from sieving
        if( .not. cline%defined('mskdiam') )then
            write(logfhandle,'(A,F8.2)')'>>> WAITING UP TO 5 MINUTES FOR '//trim(STREAM_MOLDIAM)
            do i=1, 30
                if(file_exists(trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM))) exit
                call sleep(10)
                call http_communicator%send_jobstats()
                if( http_communicator%exit )then
                    ! termination
                    write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                    call http_communicator%term()
                    call simple_end('**** SIMPLE_STREAM_PICK_EXTRACT USER STOP ****')
                    call EXIT(0)
                endif
            end do
            if( .not. file_exists(trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM))) THROW_HARD('either mskdiam must be given or '// trim(STREAM_MOLDIAM) // ' exists in target_dir')
            ! read mskdiam from file
            call moldiamori%new(1, .false.)
            call moldiamori%read( trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM) )
            if( .not. moldiamori%isthere(1, "mskdiam") ) THROW_HARD( 'mskdiam missing from ' // trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM) )
            mskdiam = moldiamori%get(1, "mskdiam")
            params%mskdiam = mskdiam
            call cline%set('mskdiam', params%mskdiam)
            write(logfhandle,'(A,F8.2)')'>>> MASK DIAMETER SET TO', params%mskdiam
        endif
        ! initialise progress monitor
        call progressfile_init()
        ! master project file
        call spproj_glob%read( params%projfile )
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_cluster2D must start from an empty project (eg from root project folder)')
        ! project watcher
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true., nretries=10)
        ! Infinite loop
        iter              = 0
        nprojects         = 0        ! # of projects per iteration
        nimported         = 0        ! # of sets per iteration
        nptcls_glob       = 0        ! global # of particles present in the pool
        l_pause           = .false.  ! pause clustering
        nsets_imported    = 0        ! Global # of sets imported
        pool_iter         = 0
        time_last_import  = huge(time_last_import)
        iter_last_import  = -1       ! Pool iteration # when set(s) last imported
        extra_pause_iters = 0        ! extra iters before pause
        do
            ! termination
            if( file_exists(trim(TERM_STREAM)) .or. http_communicator%exit ) exit
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
                ! http stats
                call http_communicator%json%update(http_communicator%job_json, "last_particles_imported",  stream_datestr(), found)
                time_last_import = time8()
                iter_last_import = get_pool_iter()
                call unpause_pool
            endif
            nsets_imported = count(setslist%imported)
            call update_user_params2D(cline, l_params_updated, nice_communicator%update_arguments)
            if( l_params_updated ) call unpause_pool
            ! pause?
            if( (pool_iter >= iter_last_import+PAUSE_NITERS+extra_pause_iters) .or.&
                & (time8()-time_last_import>PAUSE_TIMELIMIT) )then
                if( .not.l_pause )then
                    l_pause = is_pool_available()
                    extra_pause_iters = 0
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
            ! http stats
            call http_communicator%json%update(http_communicator%job_json, "stage",               "finding and classifying particles", found)
            call http_communicator%json%update(http_communicator%job_json, "particles_imported",  nptcls_glob,                         found)
            call http_communicator%json%update(http_communicator%job_json, "particles_accepted",  get_pool_assigned(),                 found)
            call http_communicator%json%update(http_communicator%job_json, "particles_rejected",  get_pool_rejected(),                 found)
            call http_communicator%json%update(http_communicator%job_json, "iteration",           last_complete_iter,                  found) ! -1 as get_pool_iter returns currently running iteration
            if(get_pool_iter() > 1) then
                call http_communicator%json%update(http_communicator%job_json, "user_input", .true., found)
                xtile = 0
                ytile = 0
                call http_communicator%json%remove(latest_cls2D, destroy=.true.)
                call http_communicator%json%create_array(latest_cls2D, "latest_cls2D")
                call http_communicator%json%add(http_communicator%job_json, latest_cls2D)
                if(allocated(pool_jpeg_map)) then
                    do i=0, size(pool_jpeg_map) - 1
                        call communicator_add_cls2D(&
                            &trim(get_pool_cavgs_jpeg()),&
                            &trim(cwd_glob) // '/' // trim(get_pool_cavgs_mrc()),&
                            &pool_jpeg_map(i + 1),&
                            &xtile * (100.0 / (get_pool_cavgs_jpeg_ntilesx() - 1)),&
                            &ytile * (100.0 / (get_pool_cavgs_jpeg_ntilesy() - 1)),&
                            &100 * get_pool_cavgs_jpeg_ntilesy(),&
                            &100 * get_pool_cavgs_jpeg_ntilesx(),&
                            &pop=pool_jpeg_pop(i + 1),&
                            &res=pool_jpeg_res(i + 1))
                        xtile = xtile + 1
                        if(xtile .eq. get_pool_cavgs_jpeg_ntilesx()) then
                            xtile = 0
                            ytile = ytile + 1
                        endif
                    end do
                endif
            endif
            if(associated(http_communicator%update_arguments) .and. last_complete_iter .gt. 0) then
                ! project snapshot if requested
                call http_communicator%json%get(http_communicator%update_arguments, "snapshot_iteration", snapshot_iteration, found) 
                if(found) then
                    call http_communicator%json%get(http_communicator%update_arguments, "snapshot_selection", snapshot_selection, found)
                    if(found) then
                        call http_communicator%json%get(http_communicator%update_arguments, "snapshot_filename", snapshot_filename, found)
                        if(found) then
                            call write_project_stream2D(&
                                &snapshot_projfile=trim(cwd_glob) // '/' // DIR_SNAPSHOT // '/' // swap_suffix(snapshot_filename, "", ".simple") // '/' //snapshot_filename,&
                                &snapshot_starfile_base=trim(cwd_glob) // '/' // DIR_SNAPSHOT // '/' // swap_suffix(snapshot_filename, "", ".simple") // '/' // swap_suffix(snapshot_filename, "", ".simple"),&
                                &optics_dir=trim(params%optics_dir))
                            call http_communicator%json%add(http_communicator%job_json, "snapshot_filename",  snapshot_filename)
                            call http_communicator%json%add(http_communicator%job_json, "snapshot_nptcls",    last_snapshot_nptcls)
                            call http_communicator%json%add(http_communicator%job_json, "snapshot_time",      stream_datestr())
                        endif
                    endif
                endif
                ! update mskdiam if requested
                call http_communicator%json%get(http_communicator%update_arguments, "mskdiam", mskdiam_update, found) 
                if(found) then
                    call update_mskdiam(mskdiam_update)
                    if( pool_iter > iter_last_import) extra_pause_iters = PAUSE_NITERS
                    time_last_import = time8()
                    call unpause_pool()
                endif
                call http_communicator%json%destroy(http_communicator%update_arguments)
                nullify(http_communicator%update_arguments)
            endif
            call http_communicator%send_jobstats()
            ! ensure snapshot info is only returned once
            call http_communicator%json%remove_if_present(http_communicator%job_json, "snapshot_filename")
            call http_communicator%json%remove_if_present(http_communicator%job_json, "snapshot_nptcls")
            call http_communicator%json%remove_if_present(http_communicator%job_json, "snapshot_time")
            ! Wait
            call sleep(WAITTIME)
        enddo
        ! Cleanup and final project
        call terminate_stream2D(optics_dir=trim(params%optics_dir))
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully 
        call http_communicator%term()
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
                    if(params%mskdiam <= 0.0) then
                        params%mskdiam = 0.9 * ceiling(params%box * params%smpd)
                        call cline%set('mskdiam', params%mskdiam)
                        write(logfhandle,'(A,F8.2)')'>>> MASK DIAMETER SET TO', params%mskdiam
                    endif
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
                character(len=XLONGSTRLEN) :: cwd_restart
                logical :: l_restart
                call simple_getcwd(cwd_restart)
                l_restart = .false.
                if(cline%defined('outdir') .and. dir_exists(trim(cline%get_carg('outdir')))) then
                    l_restart = .true.
                    call chdir(trim(cline%get_carg('outdir')))
                endif
                if(cline%defined('dir_exec')) then
                    if( .not.file_exists(cline%get_carg('dir_exec')) )then
                        THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
                    endif
                    l_restart = .true.
                endif
                if( l_restart ) then
                    write(logfhandle, *) ">>> RESTARTING EXISTING JOB", trim(cwd_restart)
                    if(cline%defined('dir_exec')) call cline%delete('dir_exec')
                    call del_file(micspproj_fname)
                    call cleanup_root_folder
                endif
                call chdir(trim(cwd_restart))
            end subroutine cleanup4restart

            subroutine communicator_init()
                call http_communicator%json%add(http_communicator%job_json, "stage",                   "initialising")
                call http_communicator%json%add(http_communicator%job_json, "particles_imported ",     0)
                call http_communicator%json%add(http_communicator%job_json, "particles_accepted",      0)
                call http_communicator%json%add(http_communicator%job_json, "particles_rejected",      0)
                call http_communicator%json%add(http_communicator%job_json, "iteration",               0)
                call http_communicator%json%add(http_communicator%job_json, "user_input",              .false.)
                call http_communicator%json%add(http_communicator%job_json, "last_particles_imported", "")
                call http_communicator%json%create_array(latest_cls2D, "latest_cls2D")
                call http_communicator%json%add(http_communicator%job_json, latest_cls2D)
            end subroutine communicator_init

            subroutine communicator_add_cls2D(path, mrcpath, mrc_idx, spritex, spritey, spriteh, spritew, pop, res)
                character(*),      intent(in)  :: path, mrcpath
                real,              intent(in)  :: spritex, spritey
                integer,           intent(in)  :: spriteh, spritew, mrc_idx
                integer, optional, intent(in)  :: pop
                real,    optional, intent(in)  :: res
                type(json_value),  pointer     :: template
                call http_communicator%json%create_object(template, "")
                call http_communicator%json%add(template, "path",    path)
                call http_communicator%json%add(template, "mrcpath", mrcpath)
                call http_communicator%json%add(template, "mrcidx",  mrc_idx)
                call http_communicator%json%add(template, "spritex", dble(spritex))
                call http_communicator%json%add(template, "spritey", dble(spritey))
                call http_communicator%json%add(template, "spriteh", spriteh)
                call http_communicator%json%add(template, "spritew", spritew)
                call http_communicator%json%add(template, "mskdiam",    nint(params%mskdiam))
                call http_communicator%json%add(template, "mskscale",   dble(params%box * params%smpd))
                if(present(pop)) call http_communicator%json%add(template, "pop", pop)
                if(present(res)) call http_communicator%json%add(template, "res", dble(res))
                call http_communicator%json%add(latest_cls2D,        template)
            end subroutine communicator_add_cls2D

    end subroutine exec_stream_abinitio2D

    subroutine stream_test_sieve_cavgs( cline )
        use simple_commanders_preprocess, only: commander_make_pickrefs, commander_extract
        use simple_commanders_project,    only: commander_new_project
        use simple_picker_utils,          only: exec_refpick
        class(cmdline),          intent(inout) :: cline
        real,                      parameter   :: LPSTART = 30., LPSTOP = 8.
        integer,                   parameter   :: MAX_NMICS      = 100
        integer,                   parameter   :: MAXPOP_DEFAULT = 100000
        character(len=*),          parameter   :: DIR_CTF    = 'ctf/'
        character(len=*),          parameter   :: DIR_PICKER = 'picker/'
        character(len=*),          parameter   :: DIR_PTCLS  = 'particles/'
        character(len=*),          parameter   :: PROJNAME = 'project', PROJFILE = 'project.simple'
        character(len=LONGSTRLEN), allocatable :: micnames(:)
        type(image),               allocatable :: pickrefs(:)
        type(parameters)                       :: params
        type(ctfparams)                        :: ctfvars
        type(ctf_estimate_iter)                :: ctfiter
        type(ori)                              :: omic
        type(oris)                             :: os_mic
        type(cmdline)                          :: cline_new_proj, cline_extract
        type(cmdline)                          :: cline_make_pickrefs, cline_sieve_cavgs
        type(commander_new_project)            :: xnew_project
        type(commander_make_pickrefs)          :: xmakepickrefs
        type(commander_stream_sieve_cavgs)     :: xsieve_cavgs
        type(sp_project)                       :: spproj
        character(len=LONGSTRLEN) :: boxfile, thumb_den
        integer                   :: nmics,nboxes,imic,nptcls,nselmics,jmic,to,iproj,np
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir',           'yes')
        ! CTF
        if( .not. cline%defined('hp')             ) call cline%set('hp',                30.)
        if( .not. cline%defined('lp')             ) call cline%set('lp',                 5.)
        if( .not. cline%defined('dfmin')          ) call cline%set('dfmin',   DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')          ) call cline%set('dfmax',   DFMAX_DEFAULT)
        if( .not. cline%defined('ctfpatch')       ) call cline%set('ctfpatch',        'yes')
        ! # of particles
        if( .not. cline%defined('maxpop')         ) call cline%set('maxpop', MAXPOP_DEFAULT)
        ! Compute
        call cline%set('nparts', 1)
        ! command-line parsing
        call params%new(cline)
        ! project creation
        call cline_new_proj%set('dir',      PATH_HERE)
        call cline_new_proj%set('projname', PROJNAME)
        call xnew_project%execute_safe(cline_new_proj)
        call spproj%read(PROJFILE)
        ! Micrographs list
        call read_filetable(params%filetab, micnames)
        nmics = size(micnames)
        write(*,*) 'Initial # of micrographs: ', nmics
        ! Picking references
        call cline_make_pickrefs%set('prg',     'make_pickrefs')
        call cline_make_pickrefs%set('mkdir',   'no')
        call cline_make_pickrefs%set('projfile', PROJFILE)
        call cline_make_pickrefs%set('nthr',     params%nthr)
        call cline_make_pickrefs%set('pickrefs', params%pickrefs)
        call cline_make_pickrefs%set('smpd',     params%smpd)
        call xmakepickrefs%execute_safe(cline_make_pickrefs)
        params%pickrefs = trim(PICKREFS_FBODY)//params%ext
        pickrefs        = read_cavgs_into_imgarr(params%pickrefs)
        write(*,*) 'Prepared picking references'
        ! Prep for CTF estimation, picking & extraction
        call simple_mkdir(DIR_CTF)
        call simple_mkdir(DIR_PICKER)
        call simple_mkdir(DIR_PTCLS)
        call simple_mkdir('spprojs/')
        call simple_mkdir('spprojs_completed/')
        call cline_extract%set('prg',     'extract')
        call cline_extract%set('mkdir',   'no')
        call cline_extract%set('nthr',     params%nthr)
        call cline_extract%set('stream',  'yes')            ! micrograph pruning
        call cline_extract%set('smpd',     params%smpd)     ! stream=yes
        call cline_extract%set('dir',      PATH_HERE)
        if( cline%defined('box') ) call cline_extract%set('box', params%box)
        ctfvars%smpd    = params%smpd
        ctfvars%kv      = params%kv
        ctfvars%cs      = params%cs
        ctfvars%fraca   = params%fraca
        ctfvars%ctfflag = CTFFLAG_YES
        call os_mic%new(nmics, is_ptcl=.false.)
        call os_mic%set_all2single('state',0)
        nboxes = 0
        iproj  = 0
        nptcls = 0
        do imic = 1, nmics, STREAM_NMOVS_SET
            to = min(nmics,imic+STREAM_NMOVS_SET-1)
            if( to == imic ) exit
            do jmic = imic,to
                ! CTF estimation
                call omic%new(is_ptcl=.false.)
                call omic%set_state(1)
                call ctfiter%iterate(ctfvars, micnames(jmic), omic, trim(DIR_CTF), l_gen_thumb=.true.)
                call omic%set('intg',    micnames(jmic)) ! for future project import
                call omic%set('imgkind', 'mic')
                call omic%set('smpd',    params%smpd)
                call omic%set('cs',      params%cs)
                call omic%set('kv',      params%kv)
                call omic%set('fraca',   params%fraca)
                call omic%set('ctf',     'yes')
                call os_mic%set_ori(jmic, omic)
                call omic%kill
                ! CTFres/icefrac selection goes here
                if( os_mic%get_state(jmic)==0 ) cycle
                ! Picking
                call exec_refpick( micnames(jmic), boxfile, thumb_den, params%smpd, np, pickrefs, DIR_PICKER )
                if( np > 0 )then
                    call os_mic%set(jmic, 'boxfile',   simple_abspath(boxfile))
                    call os_mic%set(jmic, 'thumb_den', simple_abspath(thumb_den))
                    nboxes   = nboxes   + np
                    nselmics = nselmics + 1
                else
                    call os_mic%set_state(jmic, 0)
                    cycle
                endif
            enddo
            ! Extract
            call extract2project( imic, to, np )
            nptcls = nptcls + np
            write(logfhandle,*)'>>> imic, #mics, #picks #extracted', imic, nselmics, nboxes, nptcls
            ! limits
            if( nselmics > MAX_NMICS )     exit
            if( nptcls   > params%maxpop ) exit
        end do
        write(*,*)'Picked micrographs:  ', nmics
        write(*,*)'Picked particles:    ', nboxes
        write(*,*)'Extracted particles: ', nptcls
        call dealloc_imgarr(pickrefs)
        call os_mic%kill
        ! Align & Sieve
        call cline_sieve_cavgs%set('prg',           'sieve_cavgs')
        call cline_sieve_cavgs%set('ncls',           params%ncls)
        call cline_sieve_cavgs%set('nchunks',        1)
        call cline_sieve_cavgs%set('nchunksperset',  1)
        call cline_sieve_cavgs%set('nparts',         params%nparts)
        call cline_sieve_cavgs%set('nthr',           params%nthr)
        call cline_sieve_cavgs%set('dir_target',     params%cwd)
        call cline_sieve_cavgs%set('projfile',       PROJFILE)
        call cline_sieve_cavgs%set('remove_chunks',  'no')
        if( cline%defined('nptcls_per_cls') )then
            call cline_sieve_cavgs%set('nptcls_per_cls', params%nptcls_per_cls)
        else
            call cline_sieve_cavgs%set('nptcls_per_cls', floor(real(nptcls)/real(params%ncls)))
        endif
        ! optional mskdiam
        if( cline%defined('mskdiam') ) call cline_sieve_cavgs%set('mskdiam', params%mskdiam)
        call xsieve_cavgs%execute_safe(cline_sieve_cavgs)
        ! cleanup
        call spproj%kill
        contains

            ! extract & packages into lots of particles
            subroutine extract2project( fromm, tomm, nptcls_extracted )
                integer,        intent(in)    :: fromm, tomm
                integer,        intent(inout) :: nptcls_extracted
                type(sp_project)              :: proj
                type(commander_extract)       :: xextract
                character(len=:), allocatable :: fname
                integer                       :: nm, ns
                iproj = iproj+1
                fname = 'proj_'//int2str_pad(iproj,5)//trim(METADATA_EXT)
                call simple_chdir(DIR_PTCLS)
                proj%os_mic  = os_mic%extract_subset(fromm, tomm)
                proj%compenv = spproj%compenv
                call proj%update_projinfo(fname)
                call proj%write(fname)
                call cline_extract%set('fromp',    1)
                call cline_extract%set('top',      tomm-fromm+1)
                call cline_extract%set('projfile', fname)
                call xextract%execute_safe(cline_extract)
                call proj%read_data_info(fname, nm, ns, nptcls_extracted)
                call simple_rename(fname, '../spprojs_completed/'//fname)
                call del_files(JOB_FINISHED_FBODY, 1)
                call simple_chdir('..')
                call proj%kill
            end subroutine extract2project

    end subroutine stream_test_sieve_cavgs

end module simple_commanders_stream2D
