!@descr: task 6 in the stream pipeline: global 2D refinement of pooled particles from the sieving
module simple_stream_p06_pool2D_new
use unix,                        only: SIGTERM
use simple_stream_api
use simple_stream_pool2D_utils
use simple_stream_mq_defs,       only: mq_stream_master_in, mq_stream_master_out
use simple_gui_metadata_api,     only: gui_metadata_cavg2D,                                    &
                                       gui_metadata_stream_pool2D, &
                                       sprite_sheet_pos, &
                                       GUI_METADATA_STREAM_POOL2D_TYPE,              &
                                       GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE
implicit none

public :: stream_p06_pool2D
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_p06_pool2D
  contains
    procedure :: execute => exec_stream_p06_pool2D
end type stream_p06_pool2D

contains

    ! Manages Global 2D Clustering
    ! //TODO- handling of un-classified particles
    subroutine exec_stream_p06_pool2D( self, cline )
        class(stream_p06_pool2D), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        logical,                      parameter :: SIGMAS_IN_PROJECT = .true. ! temporary flag to turn on/off getting sigma file location friom project whilst in development
        character(len=:),           allocatable :: meta_buffer
        type(parameters)                      :: params
        type(json_value), pointer             :: latest_cls2D
        type(json_core)                       :: json
        type(rec_list)                        :: setslist
        type(stream_watcher)                  :: project_buff
        type(sp_project)                      :: spproj_glob
        type(oris)                            :: moldiamori
        type(gui_metadata_cavg2D)             :: meta_cavg2D        
        type(gui_metadata_stream_pool2D)      :: meta_pool2D
        character(kind=CK,len=:), allocatable :: snapshot_filename
        type(string),             allocatable :: projects(:)
        logical,                  allocatable :: l_imported(:)
        type(string)                          :: str_avgs_jpeg, str_avgs_mrc
        integer(kind=dp) :: time_last_import, time_last_iter
        integer :: i, iter, nprojects, nimported, nptcls_glob=0, nsets_imported, pool_iter, iter_last_import
        integer :: xtile, ytile, mskdiam_update, extra_pause_iters, last_sent_iter
        logical :: l_pause, l_params_updated, found, l_terminate=.false.
        real    :: mskdiam
        call signal(SIGTERM, sigterm_handler)   ! graceful shutdown on SIGTERM
        call cline%set('oritype',      'mic')
        call cline%set('mkdir',        'yes')
        call cline%set('autoscale',    'yes')
        call cline%set('reject_mics',  'no')
        call cline%set('wiener',       'full')
        call cline%set('refine',       'snhc_smpl')
        call cline%set('ml_reg',       'no')
        call cline%set('objfun',       'euclid')
        call cline%set('sigma_est',    'global')
        call cline%set('cls_init',     'rand')
        call cline%set('numlen',       5)
        ! if( .not. cline%defined('walltime')  ) call cline%set('walltime',  29*60) ! 29 minutes
        if( .not. cline%defined('dynreslim') ) call cline%set('dynreslim', 'yes')
        if( .not.cline%defined('center')     ) call cline%set('center',    'yes')
        if( .not.cline%defined('ncls')       ) call cline%set('ncls',       200)
        ! restart
        call cleanup4restart
        ! generate own project file if projfile isnt set
        if( .not.cline%defined('projfile') )then
            call cline%set('projname', 'stream_abinitio2D')
            call cline%set('projfile', 'stream_abinitio2D.simple')
            call spproj_glob%update_projinfo(cline)
            call spproj_glob%update_compenv(cline)
            call spproj_glob%write
        endif
        ! initialise metadata
        call meta_pool2D%new(GUI_METADATA_STREAM_POOL2D_TYPE)
        call meta_cavg2D%new(GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE)
        call send_meta(string('initialising'))
        call create_stream_project(spproj_glob, cline, string('pool2D'))
        ! master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! wait if dir_target doesn't exist yet
        call send_meta(string('waiting on particle sieving'))
        call wait_for_folder2(params%dir_target)
        call wait_for_folder2(params%dir_target//'/spprojs_completed')
        ! wait for and retrieve mskdiam from sieving
        if( .not. cline%defined('mskdiam') )then
            write(logfhandle,'(A,F8.2)')'>>> WAITING UP TO 24 HOURS FOR '//STREAM_MOLDIAM
            do i=1, 8640
                if(file_exists(params%dir_target//'/'//STREAM_MOLDIAM)) exit
                call sleep(10)
            end do
            if( .not. file_exists(params%dir_target//'/'//STREAM_MOLDIAM)) THROW_HARD('either mskdiam must be given or '// STREAM_MOLDIAM // ' exists in target_dir')
            ! read mskdiam from file
            call moldiamori%new(1, .false.)
            call moldiamori%read( params%dir_target//'/'//STREAM_MOLDIAM )
            if( .not. moldiamori%isthere(1, "mskdiam") ) THROW_HARD('mskdiam missing from '//params%dir_target%to_char()//'/'//STREAM_MOLDIAM)
            mskdiam = moldiamori%get(1, "mskdiam")
            params%mskdiam = mskdiam
            call cline%set('mskdiam', params%mskdiam)
            write(logfhandle,'(A,F8.2)')'>>> MASK DIAMETER SET TO', params%mskdiam
        endif
        ! master project file
        call spproj_glob%read( params%projfile )
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_cluster2D must start from an empty project (eg from root project folder)')
        ! project watcher
        project_buff = stream_watcher(LONGTIME, params%dir_target//'/'//DIR_STREAM_COMPLETED, spproj=.true., nretries=10)
        ! Infinite loop
        iter              = 0
        last_sent_iter    = 0
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
            if( file_exists(TERM_STREAM) .or. l_terminate ) then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            iter = iter + 1
            ! detection of new projects
            call project_buff%watch(nprojects, projects)
            if( nprojects > 0 )then
                ! memoize detected projects
                call project_buff%add2history(projects)
                do i = 1,nprojects
                    call setslist%push2chunk_list(projects(i), setslist%size() + 1, .true.)
                enddo
            endif
            ! check on progress, updates particles & alignment parameters
            ! TODO: class remapping
            if( l_pause )then
                call generate_pool_stats(params)
                ! TODO Flush raw particles
            else
                ! progress
                call update_pool_status
                ! updates particles if iteration completes
                call update_pool(params)
            endif
            ! Import new particles & clustering init
            call import_sets_into_pool( nimported )
            if( nimported > 0 )then
                time_last_import = time8()
                iter_last_import = get_pool_iter()
                call unpause_pool
            endif
            l_imported = setslist%get_included_flags()
            nsets_imported = count(l_imported)
           ! call update_user_params2D(params, cline, l_params_updated, nice_comm%update_arguments)
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
                call send_meta(string('paused whilst awaiting new particles'))
            else
                ! optionally updates alignment parameters
                call update_pool_aln_params
                ! initiates new iteration
                pool_iter = get_pool_iter()
                call iterate_pool(params)
                if( get_pool_iter() > pool_iter )then
               !     nice_comm%stat_root%user_input = .true.
                    time_last_iter = time8()
                endif
                call send_meta(string('finding and classifying particles'))
            endif
            ! http stats
            if(get_pool_iter() > 1) then
                call meta_pool2D%set_user_input(.true.)
                if( last_sent_iter /= last_complete_iter ) then
                    call send_cavgs2D()
                    last_sent_iter = last_complete_iter
                endif
            endif
            !//TODO-snapshot generation and mask update
            ! if( http_communicator%arg_associated() .and. last_complete_iter .gt. 0) then
            !     ! project snapshot if requested
            !     call http_communicator%get_json_arg("snapshot_iteration", snapshot_iteration, found) 
            !     if(found) then
            !         call http_communicator%get_json_arg("snapshot_selection", snapshot_selection, found)
            !         if(found) then
            !             call http_communicator%get_json_arg("snapshot_filename", snapshot_filename, found)
            !             if(found) then
            !                 call write_project_stream2D(params,&
            !                     &snapshot_projfile=string(CWD_GLOB) // '/' // DIR_SNAPSHOT // '/' // swap_suffix(snapshot_filename, "", ".simple") // '/' //snapshot_filename,&
            !                     &snapshot_starfile_base=string(CWD_GLOB) // '/' // DIR_SNAPSHOT // '/' // swap_suffix(snapshot_filename, "", ".simple") // '/' // swap_suffix(snapshot_filename, "", ".simple"),&
            !                     &optics_dir=params%optics_dir)
            !                 call http_communicator%add_to_json("snapshot_filename", snapshot_filename)
            !                 call http_communicator%add_to_json("snapshot_nptcls",   snapshot_last_nptcls)
            !                 call http_communicator%add_to_json("snapshot_time",     stream_datestr())
            !             endif
            !         endif
            !     endif
            !     ! update mskdiam if requested
            !     call http_communicator%get_json_arg("mskdiam", mskdiam_update, found) 
            !     if(found) then
            !         call update_mskdiam(params, mskdiam_update)
            !         if( pool_iter > iter_last_import) extra_pause_iters = PAUSE_NITERS
            !         time_last_import = time8()
            !         call unpause_pool()
            !     endif
            !     call http_communicator%destroy_arg()
            ! endif
         !   call http_communicator%send_jobstats()
            ! Wait
            call sleep(WAITTIME)
        enddo
        ! Cleanup and final project
        if( allocated(l_imported) ) deallocate(l_imported)
        call terminate_stream2D(params, optics_dir=params%optics_dir)
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup(params)
        ! end gracefully 
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
                integer,           intent(out) :: nimported
                type(sp_project),  allocatable :: spprojs(:)
                type(string),      allocatable :: sigmas(:)
                class(sp_project), pointer     :: pool
                type(rec_iterator)             :: it
                type(chunk_rec)                :: crec
                logical, allocatable :: l_processed(:), l_imported(:)
                integer :: nsets2import, iset, nptcls2import, nmics2import, nmics, nptcls
                integer :: i, fromp, fromp_prev, imic, ind, iptcl, jptcl, jmic, nptcls_sel
                nimported = 0
                if( setslist%size()== 0 ) return
                if( nptcls_glob == 0 )then
                    ! first import
                else
                    ! at other times only import when the pool is free
                    if( .not.is_pool_available() ) return
                endif
                l_processed = setslist%get_processed_flags()
                l_imported  = setslist%get_included_flags()
                nsets2import = count(l_processed(:).and.(.not.l_imported(:)))
                if( nsets2import == 0 ) return
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
                    call spprojs(iset)%read_non_data_segments(crec%projfile)
                    call spprojs(iset)%read_segment('mic',    crec%projfile)
                    call spprojs(iset)%read_segment('stk',    crec%projfile)
                    call spprojs(iset)%read_segment('ptcl2D', crec%projfile)
                    nmics2import  = nmics2import  + spprojs(iset)%os_mic%get_noris()
                    nptcls2import = nptcls2import + spprojs(iset)%os_ptcl2D%get_noris()
                    call it%next()
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
                it   = setslist%begin()
                do iset = 1,setslist%size()
                    call it%get(crec)
                    if( crec%included .or. (.not.crec%processed .or. crec%busy) )then
                        ! move iterator
                        call it%next()
                        cycle
                    endif
                    fromp_prev = fromp
                    ind = 1
                    do jmic = 1,spprojs(iset)%os_mic%get_noris()
                        imic = imic + 1
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
                    write(logfhandle,'(A,I6,A,I6)')'>>> TRANSFERRED ',nptcls_sel,' PARTICLES FROM SET ',crec%id
                    call flush(logfhandle)
                    ! global list update
                    crec%included = .true.
                    call setslist%replace_iterator(it, crec)
                    ! move iterator
                    call it%next()
                enddo
                nimported = nsets2import
                ! average all previously imported sigmas
                l_imported  = setslist%get_included_flags()
                allocate(sigmas(count(l_imported)))
                i  = 0
                it = setslist%begin()
                do iset = 1,setslist%size()
                    call it%get(crec)
                    if( .not.l_imported(iset) )then
                        ! move iterator
                        call it%next()
                        cycle
                    endif
                    i         = i+1
                    ind       = crec%id
                    if(SIGMAS_IN_PROJECT) then
                        call spprojs(iset)%read_segment('out', crec%projfile)
                        call spprojs(iset)%get_sigma2(sigmas(i))
                    else
                        sigmas(i) = params%dir_target//'/set_'//int2str(ind)//'/set_'//int2str(ind)//STAR_EXT
                    endif
                    ! move iterator
                    call it%next()
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
                    call init_pool_clustering(params, cline, spproj_glob, string(MICSPPROJ_FNAME), reference_generation=.false.)
                endif
                ! global count
                nptcls_glob = nptcls_glob + nptcls_sel
                ! cleanup
                do iset = 1,setslist%size()
                    call spprojs(iset)%kill
                enddo
                if( allocated(l_imported)  ) deallocate(l_imported)
                if( allocated(l_processed) ) deallocate(l_processed)
                deallocate(spprojs)
                nullify(pool)
            end subroutine import_sets_into_pool

            ! Remove previous files from folder to restart
            subroutine cleanup4restart
                type(string) :: cwd_restart, str_dir
                logical :: l_restart
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
                    call del_file(micspproj_fname)
                    call cleanup_root_folder
                endif
                call simple_chdir(cwd_restart)
            end subroutine cleanup4restart

            ! Broadcast progress to the GUI.
            subroutine send_meta( my_stage )
                type(string), intent(in) :: my_stage
                call meta_pool2D%set(                             &
                stage              = my_stage,                    &
                iteration          = last_complete_iter,          &
                particles_imported = nptcls_glob,                 &
                particles_accepted = get_pool_assigned(),         &
                particles_rejected = get_pool_rejected())
                ! call json%add(template, "mskdiam",  nint(params%mskdiam))
               ! call json%add(template, "mskscale", dble(params%box * params%smpd))
                if( meta_pool2D%assigned() .and. mq_stream_master_in%is_active() ) then
                    call meta_pool2D%serialise(meta_buffer)
                    call mq_stream_master_in%send(meta_buffer)
                endif
            end subroutine send_meta

            subroutine send_cavg2D_meta( my_path, my_stk, my_i, my_i_max, my_xtile, my_ytile )
                type(string), intent(in) :: my_path, my_stk
                integer,      intent(in) :: my_i, my_i_max, my_xtile, my_ytile
                integer                  :: my_idx, my_xtiles, my_ytiles
                my_idx    = pool_jpeg_map(my_i)
                my_xtiles = get_pool_cavgs_jpeg_ntilesx()
                my_ytiles = get_pool_cavgs_jpeg_ntilesy()
                call meta_cavg2D%set(                                    &
                    path    = my_path,                                   &
                    mrcpath = my_stk,                                    &
                    i       = my_i,                                      &
                    i_max   = my_i_max,                                  &
                    res     = pool_jpeg_res(my_i),                       &
                    pop     = pool_jpeg_pop(my_i),                       &
                    idx     = my_idx,                                    &
                    sprite  = sprite_sheet_pos(                                                           &
                                    x = merge(0.0, my_xtile * (100.0 / (my_xtiles - 1)), my_xtiles == 1), &
                                    y = merge(0.0, my_ytile * (100.0 / (my_ytiles - 1)), my_ytiles == 1), &
                                    h = 100 * my_ytiles,                                                  &
                                    w = 100 * my_xtiles))
                if( meta_cavg2D%assigned() .and. mq_stream_master_in%is_active() ) then
                    call meta_cavg2D%serialise(meta_buffer)
                    call mq_stream_master_in%send(meta_buffer)
                endif
            end subroutine send_cavg2D_meta

            subroutine send_cavgs2D()
                type(string) :: my_path, my_stk, my_cwd
                integer :: my_i, my_xtile, my_ytile, my_n
                integer :: my_xtiles
                call simple_getcwd(my_cwd)
                my_xtiles = get_pool_cavgs_jpeg_ntilesx()
                if(allocated(pool_jpeg_map)) then
                    my_path  = my_cwd // '/' // get_pool_cavgs_jpeg()
                    my_stk   = my_cwd // '/' // get_pool_cavgs_mrc()
                    my_xtile = 0
                    my_ytile = 0
                    my_n     = size(pool_jpeg_map)
                    do my_i = 1, my_n
                        call send_cavg2D_meta(my_path, my_stk, my_i, my_n, my_xtile, my_ytile)
                        my_xtile = my_xtile + 1
                        if( my_xtile == my_xtiles ) then
                            my_xtile = 0
                            my_ytile = my_ytile + 1
                        endif
                    end do
                end if
            end subroutine send_cavgs2D   

            ! Called asynchronously on SIGTERM. Exits immediately after logging.
            subroutine sigterm_handler()
                write(logfhandle, '(A)') 'SIGTERM RECEIVED'
                l_terminate = .true.
            end subroutine sigterm_handler
            
    end subroutine exec_stream_p06_pool2D

end module simple_stream_p06_pool2D_new
