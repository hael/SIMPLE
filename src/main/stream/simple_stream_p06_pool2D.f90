module simple_stream_p06_pool2D
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_euclid_sigma2,  only: average_sigma2_groups, sigma2_star_from_iter
use simple_parameters,     only: parameters
use simple_sp_project,     only: sp_project
use simple_gui_utils
use simple_nice
use simple_progress
use simple_qsys_funs
use simple_rec_list
use simple_stream2D_state
use simple_stream_cluster2D_utils
use simple_stream_communicator
use simple_stream_pool2D_utils
use simple_stream_utils
use simple_stream_watcher
use json_kinds
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
    ! TODO: handling of un-classified particles
    subroutine exec_stream_p06_pool2D( self, cline )
        class(stream_p06_pool2D), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(parameters)                      :: params
        type(simple_nice_communicator)        :: nice_communicator
        type(stream_http_communicator)        :: http_communicator
        type(json_value), pointer             :: latest_cls2D
        type(json_core)                       :: json
        type(rec_list)                        :: setslist
        type(stream_watcher)                  :: project_buff
        type(sp_project)                      :: spproj_glob
        type(oris)                            :: moldiamori             
        character(kind=CK,len=:), allocatable :: snapshot_filename
        type(string),             allocatable :: projects(:)
        logical,                  allocatable :: l_imported(:)
        type(string)                          :: str_avgs_jpeg, str_avgs_mrc
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
        if( .not.cline%defined('projfile') )then
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
        call http_communicator%create(params%niceprocid, params%niceserver%to_char(), "classification_2D")
        call communicator_init()
        call http_communicator%send_jobstats()
        ! wait if dir_target doesn't exist yet
        call wait_for_folder(http_communicator, params%dir_target, '**** SIMPLE_STREAM_ABINITIO2D NORMAL STOP ****')
        call wait_for_folder(http_communicator, params%dir_target//'/spprojs_completed', '**** SIMPLE_STREAM_ABINITIO2D NORMAL STOP ****')
        ! wait for and retrieve mskdiam from sieving
        if( .not. cline%defined('mskdiam') )then
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
            call moldiamori%read( params%dir_target//'/'//STREAM_MOLDIAM )
            if( .not. moldiamori%isthere(1, "mskdiam") ) THROW_HARD('mskdiam missing from '//params%dir_target%to_char()//'/'//STREAM_MOLDIAM)
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
        project_buff = stream_watcher(LONGTIME, params%dir_target//'/'//DIR_STREAM_COMPLETED, spproj=.true., nretries=10)
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
            if( file_exists(TERM_STREAM) .or. http_communicator%exit_status() ) exit
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
                call http_communicator%update_json("last_particles_imported",  stream_datestr(), found)
                time_last_import = time8()
                iter_last_import = get_pool_iter()
                call unpause_pool
            endif
            l_imported = setslist%get_included_flags()
            nsets_imported = count(l_imported)
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
            call http_communicator%update_json("stage",               "finding and classifying particles", found)
            call http_communicator%update_json("particles_imported",  nptcls_glob,                         found)
            call http_communicator%update_json("particles_accepted",  get_pool_assigned(),                 found)
            call http_communicator%update_json("particles_rejected",  get_pool_rejected(),                 found)
            call http_communicator%update_json("iteration",           last_complete_iter,                  found) ! -1 as get_pool_iter returns currently running iteration
            if(get_pool_iter() > 1) then
                call http_communicator%update_json("user_input", .true., found)
                xtile = 0
                ytile = 0
                call json%remove(latest_cls2D, destroy=.true.)
                call json%create_array(latest_cls2D, "latest_cls2D")
                call http_communicator%add_to_json(latest_cls2D)
                if(allocated(pool_jpeg_map)) then
                    do i=0, size(pool_jpeg_map) - 1
                        str_avgs_jpeg = get_pool_cavgs_jpeg()
                        str_avgs_mrc  = string(CWD_GLOB) // '/' // get_pool_cavgs_mrc()
                        call add_cls2D_to_json(&
                            &str_avgs_jpeg%to_char(),&
                            &str_avgs_mrc%to_char(),&
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
                        call str_avgs_jpeg%kill
                        call str_avgs_mrc%kill
                    end do
                endif
            endif
            if( http_communicator%arg_associated() .and. last_complete_iter .gt. 0) then
                ! project snapshot if requested
                call http_communicator%get_json_arg("snapshot_iteration", snapshot_iteration, found) 
                if(found) then
                    call http_communicator%get_json_arg("snapshot_selection", snapshot_selection, found)
                    if(found) then
                        call http_communicator%get_json_arg("snapshot_filename", snapshot_filename, found)
                        if(found) then
                            call write_project_stream2D(&
                                &snapshot_projfile=string(CWD_GLOB) // '/' // DIR_SNAPSHOT // '/' // swap_suffix(snapshot_filename, "", ".simple") // '/' //snapshot_filename,&
                                &snapshot_starfile_base=string(CWD_GLOB) // '/' // DIR_SNAPSHOT // '/' // swap_suffix(snapshot_filename, "", ".simple") // '/' // swap_suffix(snapshot_filename, "", ".simple"),&
                                &optics_dir=params%optics_dir)
                            call http_communicator%add_to_json("snapshot_filename", snapshot_filename)
                            call http_communicator%add_to_json("snapshot_nptcls",   snapshot_last_nptcls)
                            call http_communicator%add_to_json("snapshot_time",     stream_datestr())
                        endif
                    endif
                endif
                ! update mskdiam if requested
                call http_communicator%get_json_arg("mskdiam", mskdiam_update, found) 
                if(found) then
                    call update_mskdiam(mskdiam_update)
                    if( pool_iter > iter_last_import) extra_pause_iters = PAUSE_NITERS
                    time_last_import = time8()
                    call unpause_pool()
                endif
                call http_communicator%destroy_arg()
            endif
            call http_communicator%send_jobstats()
            ! ensure snapshot info is only returned once
            call http_communicator%remove_from_json_if_present("snapshot_filename")
            call http_communicator%remove_from_json_if_present("snapshot_nptcls")
            call http_communicator%remove_from_json_if_present("snapshot_time")
            ! Wait
            call sleep(WAITTIME)
        enddo
        ! Cleanup and final project
        if( allocated(l_imported) ) deallocate(l_imported)
        call terminate_stream2D(optics_dir=params%optics_dir)
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
                    sigmas(i) = params%dir_target//'/set_'//int2str(ind)//'/set_'//int2str(ind)//STAR_EXT
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
                    call init_pool_clustering(cline, spproj_glob, string(MICSPPROJ_FNAME), reference_generation=.false.)
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

            subroutine communicator_init()
                call http_communicator%add_to_json("stage",                   "initialising")
                call http_communicator%add_to_json("particles_imported ",     0)
                call http_communicator%add_to_json("particles_accepted",      0)
                call http_communicator%add_to_json("particles_rejected",      0)
                call http_communicator%add_to_json("iteration",               0)
                call http_communicator%add_to_json("user_input",              .false.)
                call http_communicator%add_to_json("last_particles_imported", "")
                call json%create_array(latest_cls2D, "latest_cls2D")
                call http_communicator%add_to_json(latest_cls2D)
            end subroutine communicator_init

            subroutine add_cls2D_to_json(path, mrcpath, mrc_idx, spritex, spritey, spriteh, spritew, pop, res)
                character(*),      intent(in) :: path, mrcpath
                real,              intent(in) :: spritex, spritey
                integer,           intent(in) :: spriteh, spritew, mrc_idx
                integer, optional, intent(in) :: pop
                real,    optional, intent(in) :: res
                type(json_value),  pointer    :: template
                call json%create_object(template, "")
                call json%add(template, "path",     path)
                call json%add(template, "mrcpath",  mrcpath)
                call json%add(template, "mrcidx",   mrc_idx)
                call json%add(template, "spritex",  dble(spritex))
                call json%add(template, "spritey",  dble(spritey))
                call json%add(template, "spriteh",  spriteh)
                call json%add(template, "spritew",  spritew)
                call json%add(template, "mskdiam",  nint(params%mskdiam))
                call json%add(template, "mskscale", dble(params%box * params%smpd))
                if(present(pop)) call json%add(template, "pop", pop)
                if(present(res)) call json%add(template, "res", dble(res))
                call json%add(latest_cls2D, template)
            end subroutine add_cls2D_to_json

    end subroutine exec_stream_p06_pool2D

end module simple_stream_p06_pool2D
