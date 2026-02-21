!@descr: task 5 in the stream pipeline: chunk-based 2D clustering and automatic selection of high-quality class averages (sieving)
module simple_stream_p05_sieve_cavgs
use simple_stream_api
use simple_stream_pool2D_utils, only: set_lpthres_type
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
        type(string),                    allocatable :: projects(:), completed_projfiles(:)
        integer,                         allocatable :: accepted_cls_ids(:), rejected_cls_ids(:), jpg_cls_map(:)
        real,                            allocatable :: cls_res(:), cls_pop(:)
        logical,                         allocatable :: l_processed(:)
        integer,                           parameter :: NCLS_MAX = 100, NCLS_MIN=10, MIN_PTCLS_PER_CLASS=200
        real,                              parameter :: CHUNK_MULTIPLIER = 1.5 ! accounts for more particles found using template pick vs initial pick
        type(json_value),                    pointer :: accepted_cls2D, rejected_cls2D, latest_accepted_cls2D, latest_rejected_cls2D
        type(rec_list)                               :: project_list, chunk_list, set_list
        type(parameters)                             :: params
        type(qsys_env)                               :: qenv
        type(stream_http_communicator)               :: http_communicator
        type(oris)                                   :: moldiamori, nmicsori
        type(stream_watcher)                         :: project_buff
        type(chunk_rec)                              :: crec
        type(rec_iterator)                           :: it
        type(sp_project)                             :: spproj_glob
        type(json_core)                              :: json
        type(string)                                 :: selection_jpeg, mapfileprefix
        character(len=STDLEN)                        :: chunk_part_env
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
            call qenv%new(params, 1, exec_bin=string('simple_exec'), qsys_partition=string(trim(chunk_part_env)))
        else
            call qenv%new(params, 1, exec_bin=string('simple_exec'))
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
                ! update projects and start chunks
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
            call update_user_params2D(cline, l_params_updated)
            call update_chunks
            call memoize_chunks(chunk_list, nchunks_imported)
            call update_user_params2D(cline, l_params_updated)
            if( nchunks_imported > 0 ) nchunks_glob = nchunks_glob + nchunks_imported
            ! create sets and cluster/match cavgs when required
            call manage_sets()
            ! visualization section
            if( set_list%size() > 0 )then
                call set_list%at(1, crec)
                if( crec%processed .and. (set_list%size() > 1) .and. (.not.l_wait_for_user)) then
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
                        if( l_processed(i) ) latest_processed_set = i
                        ! move iterator
                        call it%next()
                    enddo
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
                    ! interactive selection
                    call set_list%at(1, crec)
                    if( l_wait_for_user .and. crec%waiting ) then
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
                    if( crec%waiting ) then
                        call http_communicator%get_json_arg('accepted_cls2D', accepted_cls_ids, found) ! accepted_cls2D now contains user selection
                        if(found) then
                            call http_communicator%get_json_arg('rejected_cls2D', rejected_cls_ids, found) ! rejected_cls2D now contains user selection
                            if(found) then
                                call http_communicator%update_json("user_input", .false., found)
                                ! apply interactive selection
                                write(logfhandle,'(A)') '>>> RECEIVED USER SELECTIONS', rejected_cls_ids
                                call user_select_reference_set( rejected_cls_ids )
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
            ! Create new chunks and 
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
                real              :: avgmicptcls, nptcls_per_cls, chunk_size
                integer           :: iproj, n_spprojs, n_recs_prev, irec, n_completed, nptcls, ncls
                integer           :: nmics, imic, n_ptcls, first
                n_imported = 0
                n_ptcls    = 0
                if( .not.allocated(projectnames) ) return
                n_spprojs = size(projectnames)
                if( n_spprojs == 0 )return
                n_recs_prev = project_list%size()
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
                n_completed = n_recs_prev + nmics
                n_imported  = nmics
                ! update global records and some global variables
                irec = n_recs_prev
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
                ! Updates global parameters once and init chunk 2D clustering
                if(params%nptcls_per_cls == 0) then
                    if( project_list%size() .gt. params%nmics) then
                        ! nptcls_per_cls is calculated after params%nmics processed micrographs
                        avgmicptcls    = nptcls_glob / project_list%size()
                        avgmicptcls    = ceiling(avgmicptcls / 10) * 10.0
                        chunk_size     = ceiling((CHUNK_MULTIPLIER * avgmicptcls * params%nmics) / 1000.0 ) * 1000
                        ncls           = min(NCLS_MAX, max(NCLS_MIN, ceiling(chunk_size/MIN_PTCLS_PER_CLASS)))
                        nptcls_per_cls = ceiling(chunk_size /(ncls * 10.0)) * 10.0
                        write(logfhandle,'(A,I6)')   '>>> AVERAGE # PARTICLES PER MICROGRAPH : ', int(avgmicptcls)
                        write(logfhandle,'(A,I6,A,I6,A)') '>>> USING ', ncls, ' CLASSES WITH ',int(nptcls_per_cls), ' PARTICLES PER CLASS'
                        params%nptcls_per_cls = int(nptcls_per_cls)
                        params%ncls           = ncls
                        call cline%set('nptcls_per_cls', nptcls_per_cls)
                        call cline%set('ncls',           ncls)
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
                    end if
                else if( n_recs_prev == 0 )then
                    ! nptcls_per_cls is provided on command line: initialize upon first set import
                    params%smpd = spprojs(first)%os_mic%get(1,'smpd')
                    call spprojs(first)%read_segment('stk', projectnames(first))
                    params%box  = nint(spprojs(first)%os_stk%get(1,'box'))
                    if( params%mskdiam < 0.5 ) then
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
                    write(logfhandle,'(A,A)') ">>> RESTARTING EXISTING JOB IN ", cwd_restart%to_char()
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

            subroutine user_select_reference_set( my_cls2reject )
                integer, allocatable, intent(in) :: my_cls2reject(:)
                type(sp_project)                 :: my_spproj
                type(chunk_rec)                  :: my_ref_srec
                integer, allocatable             :: my_states(:)
                integer                          :: my_i, my_ncls
                if( .not.allocated(my_cls2reject) )then
                    write(logfhandle, *) ">>> cls2reject not allocated"
                    return  ! gui must return an non empty vector
                endif
                ! get reference set record
                call set_list%at(1, my_ref_srec)
                ! read project -> only need 2d so can improve
                call my_spproj%read(my_ref_srec%projfile)
                my_ncls = my_spproj%os_cls2D%get_noris()
                if( my_ncls == 0 )then
                    write(logfhandle, *) ">>> no cls2D information in project file"
                    return 
                endif
                ! allocate states
                allocate( my_states(my_ncls) )
                ! set all state flags to 1
                my_states = 1
                call my_spproj%os_cls2D%set_all2single('state', 1)
                ! set all rejected cls state flags to 0
                do my_i=1,size(my_cls2reject)
                    my_states( my_cls2reject(my_i) ) = 0
                    call my_spproj%os_cls2D%set(my_cls2reject(my_i), 'state', 0)
                enddo
                ! write project file -> state particle mapping and junk rejection handled by cluster_cavgs_selection
                call my_spproj%write(my_ref_srec%projfile)
                ! cleanup
                deallocate(my_states)
                call my_spproj%kill()
                ! submit cluster_cavgs_selection
                call submit_cluster_cavgs_selection()
            end subroutine user_select_reference_set

            subroutine submit_aggregate_chunks( my_set_it, my_chunk_it, my_cavgs_max )
                class(rec_iterator), intent(inout) :: my_set_it, my_chunk_it
                integer,             intent(in)    :: my_cavgs_max
                type(cmdline)                      :: my_cline_aggregate_chunks
                type(string)                       :: my_cwd
                type(chunk_rec)                    :: my_crec, my_srec
                call my_chunk_it%get(my_crec)
                call my_set_it%get(my_srec)
                call my_cline_aggregate_chunks%set('prg',                             'aggregate_chunks')
                call my_cline_aggregate_chunks%set('projfile',                basename(my_srec%projfile))
                call my_cline_aggregate_chunks%set('mkdir',                                         'no')
                call my_cline_aggregate_chunks%set('nthr',                                   params%nthr)
                call my_cline_aggregate_chunks%set('mskdiam',                             params%mskdiam)
                call my_cline_aggregate_chunks%set('ncls',                                  my_cavgs_max)
                call my_cline_aggregate_chunks%set('projfile_target',    string('../')//my_crec%projfile)
                call my_cline_aggregate_chunks%set('verbose_exit',                                 'yes')
                call my_cline_aggregate_chunks%set('verbose_exit_fname',                      TARGET_MET)
                call simple_chdir(stemname(my_srec%projfile))
                call del_file(TASK_FINISHED)
                call simple_getcwd(my_cwd)
                CWD_GLOB = my_cwd%to_char()
                call qenv%exec_simple_prg_in_queue_async(my_cline_aggregate_chunks,&
                    &string('aggregate_chunks_script'), string('aggregate_chunks.log'), string('simple_private_exec'))
                call simple_chdir('..')
                call simple_getcwd(my_cwd)
                CWD_GLOB = my_cwd%to_char()
                my_srec%busy      = .true.
                my_srec%processed = .false.
                call set_list%replace_iterator(my_set_it, my_srec)
                call my_cline_aggregate_chunks%kill
                call my_cwd%kill
            end subroutine submit_aggregate_chunks

            subroutine submit_cluster_cavgs( my_set_it )
                class(rec_iterator), intent(inout) :: my_set_it
                type(cmdline)                      :: my_cline_cluster_cavgs
                type(string)                       :: my_cwd
                type(chunk_rec)                    :: my_srec
                call my_set_it%get(my_srec)
                call my_cline_cluster_cavgs%set('prg',                           'cluster_cavgs')
                call my_cline_cluster_cavgs%set('projfile',           basename(my_srec%projfile))
                call my_cline_cluster_cavgs%set('mkdir',                                    'no')
                call my_cline_cluster_cavgs%set('nthr',                              params%nthr)
                call my_cline_cluster_cavgs%set('mskdiam',                        params%mskdiam)
                call my_cline_cluster_cavgs%set('verbose_exit',                            'yes')
                call my_cline_cluster_cavgs%set('verbose_exit_fname',        CLUSTERING_COMPLETE)
                call simple_chdir(stemname(my_srec%projfile))
                call del_file(TASK_FINISHED)
                call simple_getcwd(my_cwd)
                CWD_GLOB = my_cwd%to_char()
                call qenv%exec_simple_prg_in_queue_async(my_cline_cluster_cavgs,&
                    &string('cluster_cavgs_script'), string('cluster_cavgs.log'))
                call simple_chdir('..')
                call simple_getcwd(my_cwd)
                CWD_GLOB = my_cwd%to_char()
                my_srec%busy      = .true.
                my_srec%processed = .false.
                call set_list%replace_iterator(my_set_it, my_srec)
                call my_cline_cluster_cavgs%kill
                call my_cwd%kill
                write(logfhandle, '(A,I4,A)') '>>> CLUSTERING CAVGS IN REFERENCE SET (', my_srec%id, ')'
            end subroutine submit_cluster_cavgs

            subroutine submit_match_cavgs( my_set_it )
                class(rec_iterator), intent(inout) :: my_set_it
                type(rec_iterator)                 :: my_refset_it
                type(cmdline)                      :: my_cline_match_cavgs
                type(string)                       :: my_cwd
                type(chunk_rec)                    :: my_srec, my_srec_refs
                my_refset_it = set_list%begin()
                call my_set_it%get(my_srec)
                call my_refset_it%get(my_srec_refs)
                ! exit if refset not processed
                if( .not.my_srec_refs%processed ) return
                call my_cline_match_cavgs%set('prg',                                       'match_cavgs')
                call my_cline_match_cavgs%set('projfile',                     basename(my_srec%projfile))
                call my_cline_match_cavgs%set('projfile_ref',       string('../')//my_srec_refs%projfile)
                call my_cline_match_cavgs%set('mkdir',                                              'no')
                call my_cline_match_cavgs%set('nthr',                                        params%nthr)
                call my_cline_match_cavgs%set('mskdiam',                                  params%mskdiam)
                call my_cline_match_cavgs%set('verbose_exit',                                      'yes')
                call my_cline_match_cavgs%set('verbose_exit_fname',                  CLUSTERING_COMPLETE)
                call simple_chdir(stemname(my_srec%projfile))
                call del_file(TASK_FINISHED)
                call simple_getcwd(my_cwd)
                CWD_GLOB = my_cwd%to_char()
                call qenv%exec_simple_prg_in_queue_async(my_cline_match_cavgs,&
                    &string('match_cavgs_script'), string('match_cavgs.log'))
                call simple_chdir('..')
                call simple_getcwd(my_cwd)
                CWD_GLOB = my_cwd%to_char()
                my_srec%busy      = .true.
                my_srec%processed = .false.
                call set_list%replace_iterator(my_set_it, my_srec)
                call my_cline_match_cavgs%kill
                call my_cwd%kill
                write(logfhandle, '(A,I4,A,I4,A)') '>>> MATCHING CAVGS IN SET', my_srec%id, ' TO REFERENCE SET (', my_srec_refs%id, ')'
            end subroutine submit_match_cavgs

            subroutine submit_cluster_cavgs_selection()
                type(rec_iterator)                 :: my_refset_it
                type(cmdline)                      :: my_cline_cluster_cavgs_selection
                type(string)                       :: my_cwd
                type(chunk_rec)                    :: my_srec_refs
                my_refset_it = set_list%begin()
                call my_refset_it%get(my_srec_refs)
                call my_cline_cluster_cavgs_selection%set('prg',                     'cluster_cavgs_selection')
                call my_cline_cluster_cavgs_selection%set('projfile',          basename(my_srec_refs%projfile))
                call my_cline_cluster_cavgs_selection%set('mkdir',                                        'no')
                call my_cline_cluster_cavgs_selection%set('nthr',                                  params%nthr)
                call my_cline_cluster_cavgs_selection%set('mskdiam',                            params%mskdiam)
                call my_cline_cluster_cavgs_selection%set('verbose_exit',                                'yes')
                call my_cline_cluster_cavgs_selection%set('verbose_exit_fname',                  USER_SELECTED)
                call simple_chdir(stemname(my_srec_refs%projfile))
                call del_file(TASK_FINISHED)
                call simple_getcwd(my_cwd)
                CWD_GLOB = my_cwd%to_char()
                call qenv%exec_simple_prg_in_queue_async(my_cline_cluster_cavgs_selection,&
                    &string('cluster_cavgs_selection_script'), string('cluster_cavgs_selection.log'))
                call simple_chdir('..')
                call simple_getcwd(my_cwd)
                CWD_GLOB = my_cwd%to_char()
                my_srec_refs%busy      = .true.
                my_srec_refs%processed = .false.
                call set_list%replace_iterator(my_refset_it, my_srec_refs)
                call my_cline_cluster_cavgs_selection%kill
                call my_cwd%kill
                write(logfhandle, '(A,I4,A)') '>>> SELECTING CAVGS IN REFERENCE SET (', my_srec_refs%id, ')'
            end subroutine submit_cluster_cavgs_selection

            subroutine create_new_set( my_chunk_it )
                class(rec_iterator), intent(inout) :: my_chunk_it
                type(chunk_rec)                    :: my_crec
                type(string)                       :: my_tmpl, my_set_proj
                my_tmpl     = DIR_SET//int2str(set_list%size() + 1)
                my_set_proj = my_tmpl//'/'//my_tmpl//METADATA_EXT
                call my_chunk_it%get(my_crec)
                call simple_mkdir(my_tmpl)
                call simple_copy_file(my_crec%projfile, my_set_proj)
                call simple_touch(my_tmpl//'/'//TASK_FINISHED)
                call set_list%push2chunk_list(my_set_proj, set_list%size() + 1, .false.)
                my_crec%processed = .true.
                call chunk_list%replace_iterator(my_chunk_it, my_crec)
                write(logfhandle, '(A,I4,A,I4)') '>>> CREATED SET', set_list%size(), ' AND APPENDED CHUNK', my_crec%id
            end subroutine create_new_set

            subroutine append_to_set( my_set_it, my_chunk_it )
                class(rec_iterator), intent(inout) :: my_set_it, my_chunk_it
                type(chunk_rec)                    :: my_crec, my_srec
                integer                            :: my_cavgs_max
                my_cavgs_max = SIEVING_MATCH_CAVGS_MAX                                   ! match set size
                if( my_set_it == set_list%begin() ) my_cavgs_max = SIEVING_REF_CAVGS_MAX ! reference set size
                call my_chunk_it%get(my_crec)
                call my_set_it%get(my_srec)
                my_crec%processed = .true.
                call chunk_list%replace_iterator(my_chunk_it, my_crec)
                call submit_aggregate_chunks(my_set_it, my_chunk_it, my_cavgs_max)
                write(logfhandle, '(A,I4,A,I4)') '>>> APPENDING CHUNK', my_crec%id, ' TO SET', my_srec%id
            end subroutine append_to_set

            subroutine update_set_status( my_set_it, my_l_is_ref )
                class(rec_iterator), intent(inout) :: my_set_it
                logical,             intent(in)    :: my_l_is_ref
                type(string)                       :: my_fname1, my_fname2, my_fname3, my_fname4
                type(chunk_rec)                    :: my_srec
                if( set_list%size() < 1 ) return  ! no sets generated yet
                call my_set_it%get(my_srec)
                if( my_srec%processed )      return  ! already done
                my_fname1 = stemname(my_srec%projfile)//'/'//TASK_FINISHED
                my_fname2 = stemname(my_srec%projfile)//'/'//CLUSTERING_COMPLETE
                my_fname3 = stemname(my_srec%projfile)//'/'//TARGET_MET
                my_fname4 = stemname(my_srec%projfile)//'/'//USER_SELECTED
                ! busy is set if a process is running
                if( file_exists(my_fname1) )then
                    my_srec%busy      = .false.
                else
                    my_srec%busy      = .true.
                endif
                ! included is set once the number of non-junk cavgs in set >= cavgs_max
                if( file_exists(my_fname3) )then
                    ! now ready for cavgs clustering/matching
                    my_srec%included = .true. 
                else
                    my_srec%included = .false.
                endif
                if( my_l_is_ref ) then
                    ! processed is set only when user cavg selection has been completed
                    if( file_exists(my_fname4) )then
                        ! now ready for pool import
                        my_srec%processed = .true.
                        my_srec%waiting   = .false.
                        call simple_copy_file(my_srec%projfile, string(DIR_STREAM_COMPLETED) // basename(my_srec%projfile))
                    else if( file_exists(my_fname2) ) then
                        ! waiting for user selection
                        my_srec%processed = .false.
                        my_srec%waiting   = .true.
                    else
                        my_srec%processed = .false.
                    endif
                else
                    ! processed is set only when cavgs clustering/matching has been completed
                    if( file_exists(my_fname2) )then
                        ! now ready for pool import
                        my_srec%processed = .true.
                        call simple_copy_file(my_srec%projfile, string(DIR_STREAM_COMPLETED) // basename(my_srec%projfile)) 
                    else
                        my_srec%processed = .false.
                    endif
                endif
                call set_list%replace_iterator(my_set_it, my_srec)
            end subroutine update_set_status

            subroutine manage_sets() ! note: this causes submit_cluster_cavgs/submit_match_cavgs to be run once the following chunk is complete -> can be improved
                type(rec_iterator) :: my_chunk_it, my_set_it
                type(chunk_rec)    :: my_srec, my_crec
                integer            :: my_ichunk
                if( chunk_list%size() > 0 )then
                    ! we have completed chunks
                    if( set_list%size() > 0 )then
                        ! At least 1 set created. Get last set. No new sets created until last is complete so no need to test all
                        my_set_it = set_list%end()
                        ! update final set complete status
                        call update_set_status( my_set_it, my_set_it == set_list%begin() )
                        ! get last set record
                        call my_set_it%get(my_srec)
                        ! iterate through chunks
                        my_chunk_it = chunk_list%begin()
                        do my_ichunk=1, chunk_list%size()
                            ! get chunk record
                            call my_chunk_it%get(my_crec)
                            ! get latest copy of set record
                            call my_set_it%get(my_srec)
                            ! test if chunk has been added to a set
                            if( .not.my_crec%processed )then
                                ! chunk has not been added to a set
                                if( my_srec%processed )then
                                    ! last set complete -> create new set
                                    call create_new_set(my_chunk_it)
                                    ! update l_set_it to new set
                                    my_set_it = set_list%end()
                                else if( .not.my_srec%busy .and. .not.my_srec%waiting )then
                                    if( my_srec%included ) then
                                        ! set reached target size
                                        if( my_set_it == set_list%begin() ) then
                                            ! run cluster cavgs on reference set
                                            call submit_cluster_cavgs(my_set_it)
                                        else
                                            ! run match cavgs
                                            call submit_match_cavgs(my_set_it)
                                        endif
                                    else
                                        ! last chunk added but set not reached target size -> add next chunk
                                        call append_to_set(my_set_it, my_chunk_it)
                                    endif
                                endif
                            endif
                            ! increment iterator
                            call my_chunk_it%next()
                        enddo                 
                    else
                        ! generate 1st(reference) set
                        my_chunk_it = chunk_list%begin()
                        call create_new_set(my_chunk_it)
                    endif
                endif
            end subroutine manage_sets

    end subroutine exec_stream_p05_sieve_cavgs

end module simple_stream_p05_sieve_cavgs