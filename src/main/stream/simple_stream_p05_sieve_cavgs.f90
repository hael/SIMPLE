module simple_stream_p05_sieve_cavgs
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_euclid_sigma2,  only: average_sigma2_groups
use simple_guistats,       only: guistats
use simple_linked_list,    only: list_iterator
use simple_parameters,     only: parameters
use simple_projfile_utils, only: merge_chunk_projfiles
use simple_sp_project,     only: sp_project
use simple_gui_utils
use simple_progress
use simple_qsys_env
use simple_qsys_funs
use simple_rec_list
use simple_stack_io
use simple_strategy2D_utils
use simple_stream_chunk2D_utils
use simple_stream_cluster2D_utils
use simple_stream_communicator
use simple_stream_pool2D_utils
use simple_stream_utils
use simple_stream_watcher
use json_kinds
implicit none

public :: stream_p05_sieve_cavgs
private

#include "simple_local_flags.inc"
type, extends(commander_base) :: stream_p05_sieve_cavgs
  contains
    procedure :: execute => exec_stream_p05_sieve_cavgs
end type stream_p05_sieve_cavgs

contains

    ! Manages individual chunks/sets classification, matching & rejection
    ! TODO: handling of un-classified particles
    subroutine exec_stream_p05_sieve_cavgs( self, cline )
        class(stream_p05_sieve_cavgs), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(rec_list)                 :: project_list, chunk_list, set_list
        type(string),      allocatable :: projects(:), completed_projfiles(:)
        integer,           allocatable :: accepted_cls_ids(:), rejected_cls_ids(:), jpg_cls_map(:)
        real,              allocatable :: cls_res(:), cls_pop(:)
        logical,           allocatable :: l_processed(:)
        type(parameters)               :: params
        type(qsys_env)                 :: qenv
        type(stream_http_communicator) :: http_communicator
        type(oris)                     :: moldiamori, chunksizeori, nmicsori
        type(stream_watcher)           :: project_buff
        type(chunk_rec)                :: crec
        type(rec_iterator)             :: it
        type(sp_project)               :: spproj_glob
        type(json_value), pointer      :: accepted_cls2D, rejected_cls2D, latest_accepted_cls2D, latest_rejected_cls2D
        type(json_core)                :: json
        character(len=STDLEN)          :: chunk_part_env
        type(string)     :: selection_jpeg, mapfileprefix
        real             :: mskdiam
        integer(kind=dp) :: time_last_import
        integer          :: nchunks_glob, nchunks_imported, nprojects, iter, i, envlen
        integer          :: n_imported, n_imported_prev, nptcls_glob, n_failed_jobs
        integer          :: n_accepted, n_rejected, jpg_ntiles, jpg_nxtiles, jpg_nytiles, xtile, ytile
        integer          :: latest_processed_set, latest_displayed_set, optics_map_id
        logical          :: l_params_updated, l_wait_for_user, selection_jpeg_created, found
        call cline%set('oritype',      'mic')
        call cline%set('mkdir',        'yes')
        call cline%set('autoscale',    'yes')
        call cline%set('reject_mics',  'no')
        call cline%set('reject_cls',   'no') ! refers to previous implementation
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
        if( .not.cline%defined('projfile') ) then
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
        call http_communicator%create(params%niceprocid, params%niceserver%to_char(), "sieve_cavgs")
        call communicator_init()
        call http_communicator%send_jobstats()
        ! wait if dir_target doesn't exist yet
        call wait_for_folder(http_communicator, params%dir_target, '**** SIMPLE_STREAM_SIEVE_CAVGS USER STOP ****')
        call wait_for_folder(http_communicator, params%dir_target//'/spprojs', '**** SIMPLE_STREAM_SIEVE_CAVGS USER STOP ****')
        call wait_for_folder(http_communicator, params%dir_target//'/spprojs_completed', '**** SIMPLE_STREAM_SIEVE_CAVGS USER STOP ****')
        ! mskdiam
        if( .not. cline%defined('mskdiam') )then
            ! obtain mskdiam from moldiam ori file written by pick_extract (generated by generate_pickrefs)
            write(logfhandle,'(A,F8.2)')'>>> WAITING UP TO 24 HOURS FOR '//STREAM_MOLDIAM
            do i=1, 8640
                if(file_exists(params%dir_target//'/'//STREAM_MOLDIAM)) exit
                call sleep(10)
                call http_communicator%send_jobstats()
                if( http_communicator%exit_status() )then
                    ! termination
                    write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                    call http_communicator%term()
                    call simple_end('**** SIMPLE_STREAM_PICK_EXTRACT USER STOP ****')
                    call EXIT(0)
                endif
            end do
            if( .not. file_exists(params%dir_target//'/'//STREAM_MOLDIAM)) THROW_HARD('either mskdiam must be given or '// STREAM_MOLDIAM // ' exists in target_dir')
            ! read mskdiam from file
            call moldiamori%new(1, .false.)
            call moldiamori%read(params%dir_target//'/'//STREAM_MOLDIAM)
            if( .not. moldiamori%isthere(1, "mskdiam") ) THROW_HARD( 'mskdiam missing from ' // params%dir_target%to_char()//'/'//STREAM_MOLDIAM)
            mskdiam = moldiamori%get(1, "mskdiam")
            ! write a copy for stream 2D downstream
            if(file_exists(STREAM_MOLDIAM)) call del_file(STREAM_MOLDIAM)
            call moldiamori%write(1, string(STREAM_MOLDIAM))
            call moldiamori%kill
            params%mskdiam = mskdiam
            call cline%set('mskdiam', params%mskdiam)
            write(logfhandle,'(A,F8.2)')'>>> MASK DIAMETER SET TO', params%mskdiam
            ! read nmics from file if present
            if( file_exists( params%dir_target//'/'//STREAM_NMICS) ) then
                call nmicsori%new(1, .false.)
                call nmicsori%read(params%dir_target//'/'//STREAM_NMICS )
                if( nmicsori%isthere(1, "nmics" ) ) then
                    params%nmics = nmicsori%get_int(1, "nmics")
                    write(logfhandle,'(A,I6)')'>>> NMICS SET TO', params%nmics
                endif
            endif
        endif
        ! Computing environment
        call get_environment_variable(SIMPLE_STREAM_CHUNK_PARTITION, chunk_part_env, envlen)
        if(envlen > 0) then
            call qenv%new(1, exec_bin=string('simple_exec'), qsys_partition=string(trim(chunk_part_env)))
        else
            call qenv%new(1, exec_bin=string('simple_exec'))
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
        project_buff = stream_watcher(LONGTIME, params%dir_target//'/'//DIR_STREAM_COMPLETED, spproj=.true., nretries=10)
        call simple_mkdir(PATH_HERE//DIR_STREAM_COMPLETED)
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
            if( file_exists(TERM_STREAM) .or. http_communicator%exit_status() ) exit
            if( http_communicator%stop_status() .or. test_repick() ) then
                if(test_repick()) call write_repick_refs(string("../repick_refs.mrc"))
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
            call http_communicator%update_json("stage",               "finding and sieving particles", found)    
            call http_communicator%update_json("particles_accepted",  n_accepted,                      found)
            call http_communicator%update_json("particles_rejected",  n_rejected,                      found)
            ! detection of new projects
            call project_buff%watch(nprojects, projects, max_nmovies=10*params%nparts)
            ! update global records
            if( nprojects > 0 )then
                call update_records_with_project(projects, n_imported )
                call project_buff%add2history(projects)
            endif
            ! project update
            if( nprojects > 0 )then
                n_imported = project_list%size()
                write(logfhandle,'(A,I6,I8)') '>>> # MICROGRAPHS / PARTICLES IMPORTED : ', n_imported, nptcls_glob
                ! http stats
                call http_communicator%update_json("particles_imported",       nptcls_glob,      found)
                call http_communicator%update_json("last_particles_imported",  stream_datestr(), found)
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
            call memoize_chunks(chunk_list, nchunks_imported)
            !call update_user_params2D(cline, l_params_updated, nice_communicator%update_arguments)
            call update_user_params2D(cline, l_params_updated)
            if( nchunks_imported > 0 )then
                nchunks_glob = nchunks_glob + nchunks_imported
                ! build sets
                call generate_sets(chunk_list, set_list)
            endif
            ! Sets analysis section
            if( set_list%size() > 0 )then
                call set_list%at(1, crec)
                if( crec%processed .and. crec%included .and. (set_list%size() > 1) .and. (.not.l_wait_for_user)) then
                    l_processed = set_list%get_processed_flags()
                    ! all sets but the first employ match_cavgs
                    latest_processed_set = 0
                    it = set_list%begin()
                    do i = 1,set_list%size()
                        if( i == 1 )then
                            ! move iterator
                            call it%next()
                            cycle
                        endif
                        call is_set_processed(it)
                        if( l_processed(i) ) latest_processed_set = i
                        ! move iterator
                        call it%next()
                    enddo
                    call submit_match_cavgs
                    ! http stats
                    call json%remove(accepted_cls2D, destroy=.true.)
                    call json%remove(rejected_cls2D, destroy=.true.)
                    if(latest_processed_set > 0 .and. latest_displayed_set .ne. latest_processed_set) then
                        latest_displayed_set = latest_processed_set
                        call generate_set_jpeg()
                        xtile = 0
                        ytile = 0
                        call json%remove(latest_accepted_cls2D, destroy=.true.)
                        call json%create_array(latest_accepted_cls2D, "latest_accepted_cls2D")
                        call http_communicator%add_to_json(latest_accepted_cls2D)
                        call json%remove(latest_rejected_cls2D, destroy=.true.)
                        call json%create_array(latest_rejected_cls2D, "latest_rejected_cls2D")
                        call http_communicator%add_to_json(latest_rejected_cls2D)
                        if(allocated(accepted_cls_ids) .and. allocated(rejected_cls_ids)) then
                            do i=0, size(jpg_cls_map) - 1
                                if(any( accepted_cls_ids == jpg_cls_map(i + 1))) then
                                    call add_cls2D_accepted_to_json(selection_jpeg%to_char(),&
                                        &i + 1,&
                                        &xtile * (100.0 / (jpg_nxtiles - 1)),&
                                        &ytile * (100.0 / (jpg_nytiles - 1)),&
                                        &100 * jpg_nytiles,&
                                        &100 * jpg_nxtiles,&
                                        &latest=.true.,&
                                        &res=cls_res(jpg_cls_map(i+1)),&
                                        &pop=nint(cls_pop(jpg_cls_map(i+1))))
                                else if(any( rejected_cls_ids == jpg_cls_map(i + 1))) then
                                    call add_cls2D_rejected_to_json(selection_jpeg%to_char(),&
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
                    it = set_list%begin()
                    call is_set_processed(it)
                    call submit_cluster_cavgs
                    ! interactive selection
                    call set_list%at(1, crec)
                    if( l_wait_for_user .and. crec%processed ) then
                        call http_communicator%update_json("user_input", .true., found)
                        call http_communicator%update_json("stage", "waiting for user selection", found)
                        if(.not. selection_jpeg_created) then
                            call generate_selection_jpeg()
                            xtile = 0
                            ytile = 0
                            call json%remove(accepted_cls2D, destroy=.true.)
                            call json%create_array(accepted_cls2D, "accepted_cls2D")
                            call http_communicator%add_to_json( accepted_cls2D)
                            call json%remove(rejected_cls2D, destroy=.true.)
                            call json%create_array(rejected_cls2D, "rejected_cls2D")
                            call http_communicator%add_to_json( rejected_cls2D)
                            if(allocated(accepted_cls_ids) .and. allocated(rejected_cls_ids)) then
                                do i=0, size(jpg_cls_map) - 1
                                    if(any( accepted_cls_ids == jpg_cls_map(i + 1))) then
                                        call add_cls2D_accepted_to_json(selection_jpeg%to_char(),&
                                            &jpg_cls_map(i + 1),&
                                            &xtile * (100.0 / (jpg_nxtiles - 1)),&
                                            &ytile * (100.0 / (jpg_nytiles - 1)),&
                                            &100 * jpg_nytiles,&
                                            &100 * jpg_nxtiles,&
                                            &res=cls_res(jpg_cls_map(i+1)),&
                                            &pop=nint(cls_pop(jpg_cls_map(i+1))))
                                    else if(any( rejected_cls_ids == jpg_cls_map(i + 1))) then
                                        call add_cls2D_rejected_to_json(selection_jpeg%to_char(),&
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
                if( set_list%size() > 0 ) then
                    call set_list%at(1, crec)
                    if( crec%processed ) then
                        call http_communicator%get_json_arg('accepted_cls2D', accepted_cls_ids, found) ! accepted_cls2D now contains user selection
                        if(found) then
                            call http_communicator%get_json_arg('rejected_cls2D', rejected_cls_ids, found) ! rejected_cls2D now contains user selection
                            if(found) then
                                call http_communicator%update_json("user_input", .false., found)
                                ! apply interactive selection
                                call report_interactive_selection( rejected_cls_ids )
                                write(logfhandle,'(A)') '>>> RECEIVED USER SELECTIONS'
                                ! http stats
                                call json%remove(accepted_cls2D, destroy=.true.)
                                call json%remove(rejected_cls2D, destroy=.true.)
                                xtile = 0
                                ytile = 0
                                call json%remove(latest_accepted_cls2D, destroy=.true.)
                                call json%create_array(latest_accepted_cls2D, "latest_accepted_cls2D")
                                call http_communicator%add_to_json( latest_accepted_cls2D)
                                call json%remove(latest_rejected_cls2D, destroy=.true.)
                                call json%create_array(latest_rejected_cls2D, "latest_rejected_cls2D")
                                call http_communicator%add_to_json( latest_rejected_cls2D)
                                if(allocated(accepted_cls_ids) .and. allocated(rejected_cls_ids)) then
                                    do i=0, size(jpg_cls_map) - 1
                                        if(any( accepted_cls_ids == jpg_cls_map(i + 1))) then
                                            call add_cls2D_accepted_to_json(selection_jpeg%to_char(),&
                                                &jpg_cls_map(i + 1),&
                                                &xtile * (100.0 / (jpg_nxtiles - 1)),&
                                                &ytile * (100.0 / (jpg_nytiles - 1)),&
                                                &100 * jpg_nytiles,&
                                                &100 * jpg_nxtiles,&
                                                &latest=.true.,&
                                                &res=cls_res(i+1),&
                                                &pop=nint(cls_pop(i+1)))
                                        else if(any( rejected_cls_ids == jpg_cls_map(i + 1))) then
                                            call add_cls2D_rejected_to_json(selection_jpeg%to_char(),&
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
            endif
            ! make completed sets available to abinitio2D_stream
            if ( params%interactive == 'no' ) then
                ! for non-iteractive
                call flag_complete_sets
            else if ( .not. l_wait_for_user) then
                ! interactive -> l_wait_for_user must be false (ie user selection performed) before completing any sets
                call flag_complete_sets
            endif
            ! 2D analyses
            call analyze2D_new_chunks(project_list)
            call sleep(WAITTIME)
            ! http stats send
            call http_communicator%send_jobstats()
        end do
        ! termination
        write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
        call terminate_chunks
        ! merge sets and write project
        it = set_list%begin()
        do i=1, set_list%size()
            call it%get(crec)
            if(.not. crec%busy .and. crec%processed .and. crec%included) then
                if(.not. allocated(completed_projfiles))then
                    allocate(completed_projfiles(1))
                    completed_projfiles(1) = crec%projfile
                else
                    completed_projfiles = [completed_projfiles, crec%projfile]
                endif
            endif
            ! move iterator
            call it%next()
        enddo
        if(allocated(completed_projfiles)) then
            call merge_chunk_projfiles(completed_projfiles, string('./'), spproj_glob, projname_out=string("tmp"), write_proj=.false.)
            call spproj_glob%update_projinfo(cline%get_carg('projfile')) ! update projinfo with projfile name as modified by merge_chunks
            deallocate(completed_projfiles)
        endif
        ! add optics
        if(cline%defined('optics_dir')) then
            optics_map_id = get_latest_optics_map_id(params%optics_dir)
            if(optics_map_id .gt. 0) then
                mapfileprefix = params%optics_dir // '/' // OPTICS_MAP_PREFIX // int2str(optics_map_id)
                call spproj_glob%import_optics_map(mapfileprefix)
            endif
        endif
        ! write project and star files (just in case you want ot import these particles/micrographs elsewhere)
        call spproj_glob%write
        call spproj_glob%write_mics_star(string("micrographs.star"))
        call spproj_glob%write_ptcl2D_star(string("particles.star"))
        ! cleanup
        if( allocated(l_processed) ) deallocate(l_processed)
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
        call http_communicator%term()
        call simple_end('**** SIMPLE_STREAM_SIEVE_CAVGS NORMAL STOP ****')
        
        contains

            ! updates global records
            subroutine update_records_with_project( projectnames, n_imported )
                type(string), allocatable, intent(in)  :: projectnames(:)
                integer,                   intent(out) :: n_imported
                type(sp_project), allocatable :: spprojs(:)
                type(project_rec) :: prec
                type(string)      :: fname, abs_fname
                real              :: avgmicptcls, nptcls_per_cls
                integer :: iproj, n_spprojs, n_old, irec, n_completed, nptcls, nmics, imic, n_ptcls, first
                n_imported = 0
                n_ptcls    = 0
                if( .not.allocated(projectnames) ) return
                n_spprojs = size(projectnames)
                if( n_spprojs == 0 )return
                n_old = 0 ! on first import
                allocate(spprojs(n_spprojs))
                ! because pick_extract purges state=0 and nptcls=0 mics,
                ! all mics can be assumed associated with particles
                nmics = 0
                first = 0
                do iproj = 1,n_spprojs
                    call spprojs(iproj)%read_segment('mic', projectnames(iproj))
                    nmics = nmics + spprojs(iproj)%os_mic%get_noris()
                    if( (first == 0) .and. (nmics > 0) ) first = iproj
                enddo
                if( nmics == 0 ) return
                ! import micrographs
                n_completed = n_old + nmics
                n_imported  = nmics
                ! update global records and some global variables
                irec = n_old
                do iproj = 1,n_spprojs
                    do imic = 1,spprojs(iproj)%os_mic%get_noris()
                        irec      = irec + 1
                        nptcls    = spprojs(iproj)%os_mic%get_int(imic,'nptcls')
                        n_ptcls   = n_ptcls + nptcls ! global update
                        fname     = projectnames(iproj)
                        abs_fname = simple_abspath(fname)
                        prec%projname   = abs_fname
                        prec%micind     = imic
                        prec%nptcls     = nptcls
                        prec%nptcls_sel = nptcls
                        prec%included   = .false.
                        call project_list%push_back(prec)
                    enddo
                enddo
                nptcls_glob = nptcls_glob + n_ptcls ! global update
                ! Updates global parameters once and init 2D
                if(params%nptcls_per_cls == 0) then
                    if( project_list%size() .gt. params%nmics) then
                        avgmicptcls    = nptcls_glob / project_list%size()
                        avgmicptcls    = ceiling(avgmicptcls / 10) * 10.0
                        nptcls_per_cls = ceiling(params%nmics * avgmicptcls / params%ncls)
                        write(logfhandle,'(A,I6)')   '>>> AVERAGE # PARTICLES PER MICROGRAPH : ', int(avgmicptcls)
                        write(logfhandle,'(A,I6,A)') '>>> USING ', int(nptcls_per_cls), ' PARTICLES PER CLASS'
                        params%nptcls_per_cls = int(nptcls_per_cls)
                        call cline%set('nptcls_per_cls', nptcls_per_cls)
                        params%smpd = spprojs(first)%os_mic%get(1,'smpd')
                        call spprojs(first)%read_segment('stk', projectnames(first))
                        params%box  = nint(spprojs(first)%os_stk%get(1,'box'))
                        if(params%mskdiam == 0.0) then
                            params%mskdiam = 0.9 * ceiling(params%box * spprojs(first)%os_stk%get(1,'smpd'))
                            call cline%set('mskdiam', params%mskdiam)
                            write(logfhandle,'(A,F8.2)')'>>> MASK DIAMETER SET TO', params%mskdiam
                        endif
                        call init_chunk_clustering(cline, spproj_glob)
                        call cline%delete('ncls')
                        ! write out for stream3d to pick up
                        call chunksizeori%new(1, .false.)
                        call chunksizeori%set(1, 'nptcls_per_cls', params%nptcls_per_cls)
                        call chunksizeori%write(1, string(STREAM_CHUNKSIZE))
                        call chunksizeori%kill
                    end if
                else if( n_old == 0 )then
                    params%smpd = spprojs(first)%os_mic%get(1,'smpd')
                    call spprojs(first)%read_segment('stk', projectnames(first))
                    params%box  = nint(spprojs(first)%os_stk%get(1,'box'))
                    if(params%mskdiam == 0.0) then
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
                class(rec_list), intent(in)    :: chunks
                class(rec_list), intent(inout) :: sets
                type(string),      allocatable :: starfiles(:), projfiles(:)
                integer,           allocatable :: ids(:)
                type(sp_project) :: spproj
                type(string)     :: tmpl
                integer          :: navail_chunks, n, iset, i, ic, ic_start, ic_end
                navail_chunks = chunks%size() - sets%size() * params%nchunksperset
                n = floor(real(navail_chunks) / real(params%nchunksperset))
                if( n < 1 )return
                do iset = 1,n
                    ! merge chunks project into designated folder
                    ic_start = sets%size() * params%nchunksperset + 1
                    ic_end   = ic_start + params%nchunksperset - 1
                    tmpl     = DIR_SET//int2str(sets%size() + 1)
                    call simple_mkdir(tmpl)
                    projfiles = chunks%get_projfiles([ic_start,ic_end])
                    call merge_chunk_projfiles(projfiles, tmpl, spproj, projname_out=tmpl)
                    ! average and stash sigma2
                    allocate(starfiles(params%nchunksperset))
                    ids = chunks%get_ids()
                    do i = 1,params%nchunksperset
                        ic = ic_start + i - 1
                        starfiles(i) = SIGMAS_DIR//'/chunk_'//int2str(ids(ic))//STAR_EXT
                    enddo
                    call average_sigma2_groups(tmpl//'/'//tmpl//STAR_EXT, starfiles)
                    ! update global list
                    call sets%push2chunk_list(tmpl//'/'//tmpl//METADATA_EXT, sets%size() + 1, .false.)
                    ! remove imported chunk
                    if( trim(params%remove_chunks).eq.'yes' )then
                        do ic = ic_start,ic_end
                            call simple_rmdir(stemname(projfiles(ic)))
                        enddo
                    endif
                    call projfiles%kill
                    call starfiles%kill
                    if( allocated(ids) ) deallocate(ids)
                enddo
                call spproj%kill
            end subroutine generate_sets

            subroutine submit_cluster_cavgs
                type(cmdline)   :: cline_cluster_cavgs
                type(string)    :: cwd
                type(chunk_rec) :: crec
                if( set_list%size() < 1 ) return  ! no sets generated yet
                call set_list%at(1, crec)
                if( crec%busy )           return  ! ongoing
                if( crec%processed )      return  ! already done
                call cline_cluster_cavgs%set('prg',          'cluster_cavgs')
                call cline_cluster_cavgs%set('projfile',     basename(crec%projfile))
                call cline_cluster_cavgs%set('mkdir',        'no')
                call cline_cluster_cavgs%set('nthr',         params%nthr)
                call cline_cluster_cavgs%set('mskdiam',      params%mskdiam)
                call cline_cluster_cavgs%set('verbose_exit', 'yes')
                call simple_chdir(stemname(crec%projfile))
                call simple_getcwd(cwd)
                CWD_GLOB = cwd%to_char()
                call qenv%exec_simple_prg_in_queue_async(cline_cluster_cavgs,&
                    &string('cluster_cavgs_script'), string('cluster_cavgs.log'))
                call simple_chdir('..')
                call simple_getcwd(cwd)
                CWD_GLOB = cwd%to_char()
                crec%busy      = .true.
                crec%processed = .false.
                call set_list%replace_at(1, crec)
                call cline_cluster_cavgs%kill
            end subroutine submit_cluster_cavgs

            subroutine submit_match_cavgs
                type(cmdline)        :: cline_match_cavgs
                type(string)         :: path, cwd
                type(chunk_rec)      :: crec
                integer              :: iset
                type(rec_iterator)   :: it
                logical, allocatable :: l_processed(:), l_busy(:)
                if( set_list%size() < 2 ) return ! not enough sets generated
                l_processed = set_list%get_processed_flags()
                l_busy      = set_list%get_busy_flags()
                call set_list%at(1, crec)
                ! any unprocessed and not being processed?
                if ( all(l_processed(2:) .or. l_busy(2:)) ) return
                call cline_match_cavgs%set('prg',          'match_cavgs')
                call cline_match_cavgs%set('projfile',     simple_abspath(crec%projfile))
                call cline_match_cavgs%set('mkdir',        'no')
                call cline_match_cavgs%set('nthr',         params%nthr)
                call cline_match_cavgs%set('mskdiam',      params%mskdiam)
                call cline_match_cavgs%set('verbose_exit', 'yes')
                call cline_match_cavgs%delete('nparts')
                it = set_list%begin()
                do iset = 1,set_list%size()
                    call it%get(crec)
                    if( iset == 1 .or. (crec%processed .or. crec%busy) )then
                        call it%next()
                        cycle
                    endif
                    call cline_match_cavgs%set('projfile_target', basename(crec%projfile))
                    path = stemname(crec%projfile)
                    call simple_chdir(path)
                    call simple_getcwd(cwd)
                    CWD_GLOB = cwd%to_char()
                    call qenv%exec_simple_prg_in_queue_async(cline_match_cavgs,&
                        &string('match_cavgs_script'), string('match_cavgs.log'))
                    call simple_chdir('..')
                    call simple_getcwd(cwd)
                    CWD_GLOB = cwd%to_char()
                    crec%busy      = .true.   ! ongoing
                    crec%processed = .false.  ! not complete
                    call set_list%replace_iterator(it, crec)
                    ! move iterator
                    call it%next()
                enddo
                if( allocated(l_processed) ) deallocate(l_processed)
                if( allocated(l_busy)      ) deallocate(l_busy)
                call cline_match_cavgs%kill
            end subroutine submit_match_cavgs

            ! Check for status of individual sets
            subroutine is_set_processed( it )
                class(rec_iterator), intent(inout) :: it
                type(string)    :: fname
                type(chunk_rec) :: crec
                if( set_list%size() < 1 ) return  ! no sets generated yet
                call it%get(crec)
                if( crec%processed )      return  ! already done
                fname = stemname(crec%projfile)//'/'//TASK_FINISHED
                if( file_exists(fname) )then
                    crec%busy      = .false.
                    crec%processed = .true.       ! now ready for pool import
                else
                    crec%processed = .false.
                endif
                call set_list%replace_iterator(it, crec)
            end subroutine is_set_processed

            ! apply user-inputted selection on the first set
            subroutine report_interactive_selection( cls2reject )
                integer, allocatable, intent(in) :: cls2reject(:)
                integer, allocatable             :: states(:)
                type(sp_project)                 :: spproj
                integer, allocatable :: clusters_accepted(:)
                type(chunk_rec) :: crec
                integer :: i, n, n_clusters, cluster
                logical :: l_class_selection
                if( .not.allocated(cls2reject) )then
                    write(logfhandle, *) ">>> cls2reject not allocated"
                    return  ! gui must return an non empty vector
                endif
                ! read all fields
                call set_list%at(1, crec)
                call spproj%read(crec%projfile)
                n = spproj%os_cls2D%get_noris()
                ! undo the default particles cluster_cavgs selection
                call spproj%os_ptcl2D%set_all2single('state', 1)
                call spproj%os_ptcl3D%set_all2single('state', 1)
                ! selection
                allocate(states(n), source=1)
                l_class_selection = .true.
                if( size(cls2reject) == 1 )then
                    ! no rejection applied, all classes are to be selected
                    l_class_selection = cls2reject(1) /= 0
                endif
                if( l_class_selection )then
                    ! get all clusters and calculate their populations
                    n_clusters = maxval(spproj%os_cls2D%get_all_asint('cluster'))
                    allocate(clusters_accepted(n_clusters), source=0)
                    do i=1, spproj%os_cls2D%get_noris()
                        ! set clusters_accepted(cluster) to 1 if any child class is selected
                        if( any(cls2reject == i) ) cycle
                        cluster = spproj%os_cls2D%get_int(i, 'cluster')
                        if( cluster > 0 ) clusters_accepted(cluster) = 1
                    enddo
                    do i=1, size(clusters_accepted)
                        if( clusters_accepted(i) > 0) then
                            write(logfhandle, *) ">>> CLUSTER ", i, " HAS BEEN SELECTED"
                        else
                            write(logfhandle, *) ">>> CLUSTER ", i, " HAS BEEN DESELECTED"
                        endif
                    enddo
                    do i=1, spproj%os_cls2D%get_noris()
                        cluster = spproj%os_cls2D%get_int(i, 'cluster')
                        if( cluster == 0 ) then
                            states(i) = 0
                        else
                            states(i) = clusters_accepted(cluster)
                        endif
                    enddo
                endif
                ! maintain state=0 for junk
                do i = 1, n
                    if(.not. spproj%os_cls2D%isthere(i, 'cluster')) states(i) = 0
                enddo
                ! report selection to particles
                call spproj%map_cavgs_selection(states)
                call spproj%write(crec%projfile)
                ! cleanup
                call spproj%kill
                if( allocated(clusters_accepted) ) deallocate(clusters_accepted)
                deallocate(states)
            end subroutine report_interactive_selection

            ! make completed project files visible to the watcher of the next application
            subroutine flag_complete_sets
                use simple_image, only:image
                type(sp_project)   :: spproj
                type(image)        :: img
                type(string)       :: destination, stk
                type(rec_iterator) :: it
                type(chunk_rec)    :: crec
                real    :: smpd
                integer :: ldim(3), icls, iset, n_state_nonzero, nimgs, ncls
                it = set_list%begin()
                do iset = 1,set_list%size()
                    call it%get(crec)
                    if( crec%included )then
                        ! move iterator
                        call it%next() 
                        cycle
                    endif
                    if( crec%processed )then
                        destination = DIR_STREAM_COMPLETED//DIR_SET//int2str(iset)//METADATA_EXT
                        call simple_rename(crec%projfile, destination)
                        crec%projfile = destination ! relocation
                        crec%included = .true.
                        write(logfhandle,'(A,I3)')'>>> COMPLETED SET ', crec%id
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
                            call find_ldim_nptcls(string(REJECTED_CLS_STACK), ldim, nimgs)
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
                                call img%write(string(REJECTED_CLS_STACK), nimgs)
                            endif
                        enddo
                        call spproj%kill
                    endif
                    ! move iterator
                    call it%next() 
                enddo
                call img%kill
            end subroutine flag_complete_sets

            subroutine communicator_init()
                call http_communicator%add_to_json( "stage",                   "initialising")
                call http_communicator%add_to_json( "particles_imported ",     0)
                call http_communicator%add_to_json( "particles_accepted",      0)
                call http_communicator%add_to_json( "particles_rejected",      0)
                call http_communicator%add_to_json( "user_input",              .false.)
                call http_communicator%add_to_json( "last_particles_imported", "")
                call json%create_array(accepted_cls2D, "accepted_cls2D")
                call http_communicator%add_to_json( accepted_cls2D)
                call json%create_array(rejected_cls2D, "rejected_cls2D")
                call http_communicator%add_to_json( rejected_cls2D)
                call json%create_array(latest_accepted_cls2D, "latest_accepted_cls2D")
                call http_communicator%add_to_json( latest_accepted_cls2D)
                call json%create_array(latest_rejected_cls2D, "latest_rejected_cls2D")
                call http_communicator%add_to_json( latest_rejected_cls2D)
            end subroutine communicator_init

            ! Remove previous files from folder to restart
            subroutine cleanup4restart
                type(string), allocatable :: folders(:)
                type(string) :: cwd_restart, str_dir
                logical      :: l_restart
                integer      :: i
                call simple_getcwd(cwd_restart)
                l_restart = .false.
                if(cline%defined('outdir') .and. dir_exists(cline%get_carg('outdir'))) then
                    l_restart = .true.
                    call simple_chdir(cline%get_carg('outdir'))
                endif
                if(cline%defined('dir_exec')) then
                    if( .not.file_exists(cline%get_carg('dir_exec')) )then
                        str_dir = cline%get_carg('dir_exec')
                        THROW_HARD('Previous directory does not exists: '//str_dir%to_char())
                        call str_dir%kill
                    endif
                    l_restart = .true.
                endif
                if( l_restart ) then
                    write(logfhandle,'(A)') ">>> RESTARTING EXISTING JOB", cwd_restart%to_char()
                    if(cline%defined('dir_exec')) call cline%delete('dir_exec')
                    call del_file(TERM_STREAM)
                    call del_file(USER_PARAMS2D)
                    call simple_rmdir(SIGMAS_DIR)
                    call simple_rmdir(DIR_STREAM_COMPLETED)
                    folders = simple_list_dirs('.')
                    if( allocated(folders) )then
                        do i = 1,size(folders)
                            if( folders(i)%has_substr(DIR_CHUNK).or.folders(i)%has_substr(DIR_SET) )then
                                call simple_rmdir(folders(i))
                            endif
                        enddo
                    endif
                endif
                call simple_chdir(cwd_restart)
            end subroutine cleanup4restart

            subroutine generate_selection_jpeg()
                type(sp_project) :: set1_proj
                integer          :: ncls, icls
                real             :: smpd, stkbox
                type(chunk_rec)  :: crec
                call set_list%at(1, crec)
                call set1_proj%read_segment('cls2D', crec%projfile)
                call set1_proj%read_segment('out',   crec%projfile)
                call set1_proj%read_segment('stk',   crec%projfile)
                call set1_proj%get_cavgs_stk(selection_jpeg, ncls, smpd, box=stkbox)
                call mrc2jpeg_tiled(selection_jpeg, swap_suffix(selection_jpeg, JPG_EXT, params%ext%to_char()), ntiles=jpg_ntiles, n_xtiles=jpg_nxtiles, n_ytiles=jpg_nytiles, mskdiam_px=ceiling((params%mskdiam * stkbox) / (smpd * set1_proj%get_box())))
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
                selection_jpeg = swap_suffix(selection_jpeg, JPG_EXT, params%ext%to_char())
            end subroutine generate_selection_jpeg

            subroutine generate_set_jpeg()
                type(sp_project) :: set_proj
                integer          :: ncls, icls
                real             :: smpd, stkbox
                type(chunk_rec)  :: crec
                call set_list%at(latest_processed_set, crec)
                call set_proj%read_segment('cls2D', crec%projfile)
                call set_proj%read_segment('out',   crec%projfile)
                call set_proj%read_segment('stk',   crec%projfile)
                call set_proj%get_cavgs_stk(selection_jpeg, ncls, smpd, box=stkbox)
                call mrc2jpeg_tiled(selection_jpeg, swap_suffix(selection_jpeg, JPG_EXT, params%ext%to_char()), ntiles=jpg_ntiles, n_xtiles=jpg_nxtiles, n_ytiles=jpg_nytiles, mskdiam_px=ceiling((params%mskdiam * stkbox) / (smpd * set_proj%get_box())))
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
                selection_jpeg = swap_suffix(selection_jpeg, JPG_EXT, params%ext%to_char())
            end subroutine generate_set_jpeg

            subroutine add_cls2D_accepted_to_json(path, idx, spritex, spritey, spriteh, spritew, latest, res, pop)
                character(*),      intent(in) :: path
                real,              intent(in) :: spritex, spritey
                integer,           intent(in) :: spriteh, spritew, idx
                logical, optional, intent(in) :: latest
                integer, optional, intent(in) :: pop
                real,    optional, intent(in) :: res
                type(json_value),  pointer    :: template
                logical                       :: l_latest = .false.
                if(present(latest)) l_latest = latest
                call json%create_object(template, "")
                call json%add(template, "path",    path)
                call json%add(template, "spritex", dble(spritex))
                call json%add(template, "spritey", dble(spritey))
                call json%add(template, "spriteh", spriteh)
                call json%add(template, "spritew", spritew)
                call json%add(template, "idx",     idx)
                if(present(res)) call json%add(template, "res",     dble(res))
                if(present(pop)) call json%add(template, "pop",     pop)
                if(l_latest) then
                    call json%add(latest_accepted_cls2D, template)
                else
                    call json%add(accepted_cls2D, template)
                endif
            end subroutine add_cls2D_accepted_to_json

            subroutine add_cls2D_rejected_to_json(path, idx, spritex, spritey, spriteh, spritew, latest, res, pop)
                character(*),      intent(in) :: path
                real,              intent(in) :: spritex, spritey
                integer,           intent(in) :: spriteh, spritew, idx
                logical, optional, intent(in) :: latest
                integer, optional, intent(in) :: pop
                real,    optional, intent(in) :: res
                type(json_value),  pointer    :: template
                logical                       :: l_latest = .false.
                if(present(latest)) l_latest = latest
                call json%create_object(template, "")
                call json%add(template, "path",    path)
                call json%add(template, "spritex", dble(spritex))
                call json%add(template, "spritey", dble(spritey))
                call json%add(template, "spriteh", spriteh)
                call json%add(template, "spritew", spritew)
                call json%add(template, "idx",     idx)
                if(present(res)) call json%add(template, "res",     dble(res))
                if(present(pop)) call json%add(template, "pop",     pop)
                if(l_latest) then
                    call json%add(latest_rejected_cls2D, template)
                else
                    call json%add(rejected_cls2D, template)
                endif
            end subroutine add_cls2D_rejected_to_json

    end subroutine exec_stream_p05_sieve_cavgs

end module simple_stream_p05_sieve_cavgs