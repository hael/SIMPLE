!@descr: task 3 in the stream pipeline: the first 2D analyis from segmentation picked particles used for initial screening and generation of picking references
module simple_stream_p03_opening2D_new
use unix,                         only: SIGTERM, c_time
use simple_stream_api
use simple_stream_mq_defs,        only: mq_stream_master_in, mq_stream_master_out
use simple_commanders_pick,       only: commander_extract
use simple_commanders_cavgs,      only: commander_shape_rank_cavgs
use simple_commanders_abinitio2D, only: commander_abinitio2D
use simple_mini_stream_utils,     only: segdiampick_mics
use simple_gui_metadata_api

implicit none

public :: stream_p03_opening2D
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_p03_opening2D
  contains
    procedure :: execute => exec_stream_p03_opening2D
end type stream_p03_opening2D

contains

    subroutine exec_stream_p03_opening2D( self, cline )
        implicit none
        class(stream_p03_opening2D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)                 :: params
        type(stream_http_communicator)   :: http_communicator, http_gen_pickrefs_communicator
        type(stream_watcher)             :: project_buff
        type(sp_project)                 :: spproj, spproj_part
        type(cmdline)                    :: cline_extract, cline_abinitio2D, cline_shape_rank
        type(commander_extract)          :: xextract
        type(commander_abinitio2D)       :: xabinitio2D
        type(commander_shape_rank_cavgs) :: xshape_rank
        type(oris)                       :: nmics_ori
        type(json_value),    pointer     :: latest_picked_micrographs, latest_cls2D, selected_refs
        type(json_core)                  :: json
        type(gui_metadata_stream_update)          :: meta_update
        type(gui_metadata_stream_initial_picking) :: meta_initial_picking
        type(gui_metadata_stream_opening2D)       :: meta_opening2D
        type(gui_metadata_micrograph)             :: meta_micrograph
        type(string),        allocatable :: projects(:)
        character(len=:),    allocatable :: buffer
        character(len=:),    allocatable :: meta_buffer
        type(string)                     :: final_selection_source, cavgsstk, mapfileprefix, str_dir, str_thumb, str_box, projfile
        integer,             allocatable :: cavg_inds(:)
        integer,             parameter   :: NCLS_MIN = 10, NCLS_MAX = 100, NPARTS2D = 4, NTHUMB_MAX = 10
        real,                parameter   :: LPSTOP = 8.
        integer,             allocatable :: final_selection(:)
        integer                          :: nprojects, i, nptcls, ncls, nthr2D, box_in_pix=0, box_for_pick, box_for_extract
        integer                          :: ithumb, xtiles, ytiles, xtile, ytile, ncls_stk, iori, i_max, optics_map_id
        integer                          :: n_non_zero, nmics
        logical                          :: found, increase_nmics = .false.
        logical                          :: restart_requested, l_terminate=.false.
        real                             :: mskdiam_estimate=0, smpd_stk
        call signal(SIGTERM, sigterm_handler)
        if( .not. cline%defined('dir_target')     ) THROW_HARD('DIR_TARGET must be defined!')
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir',            'yes')
        if( .not. cline%defined('nptcls_per_cls') ) call cline%set('nptcls_per_cls',     200)
        if( .not. cline%defined('pick_roi')       ) call cline%set('pick_roi',         'yes')
        if( .not. cline%defined('outdir')         ) call cline%set('outdir',              '')
        if( .not. cline%defined('nmics')          ) call cline%set('nmics',              100)
        ! sanity check for restart
        if(cline%defined('outdir') .and. dir_exists(cline%get_carg('outdir'))) then
            write(logfhandle,'(A)') ">>> RESTARTING EXISTING JOB"
            call del_file(TERM_STREAM)
        endif
        ! generate own project file if projfile doesn't exist
        projfile = cline%get_carg('projfile')
        if( .not.file_exists(projfile) )then
            call cline%set('projfile', projfile)
            call spproj%update_projinfo(cline)
            call spproj%update_compenv(cline)
            call spproj%write()
        endif
        ! master parameters
        call params%new(cline)
        ! initialise metadata
        call meta_update%new(                        GUI_METADATA_STREAM_UPDATE_TYPE)
        call meta_initial_picking%new(      GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
        call meta_opening2D%new(                  GUI_METADATA_STREAM_OPENING2D_TYPE)
        call meta_micrograph%new(GUI_METADATA_STREAM_INITIAL_PICKING_MICROGRAPH_TYPE)
        
        ! ---------- Restartable block: loop instead of GOTO ----------
        restart_requested = .false.
        do
            restart_requested = .false.
            if(file_exists(STREAM_NMICS)) call del_file(STREAM_NMICS)
            call nmics_ori%new(1,          .false.)
            call nmics_ori%set(1, 'nmics', params%nmics)
            call nmics_ori%write(1, string(STREAM_NMICS))
            call nmics_ori%kill
            ! master project file
            call spproj%read( params%projfile )
            if( spproj%os_mic%get_noris() /= 0 ) call spproj%os_mic%new(1, .false.) !!!!!!!?????
            ! wait if dir_target doesn't exist yet
            call wait_for_folder(http_communicator, params%dir_target, '**** SIMPLE_GEN_PICKREFS USER STOP ****')
            call wait_for_folder(http_communicator, params%dir_target//'/spprojs', '**** SIMPLE_GEN_PICKREFS USER STOP ****')
            call wait_for_folder(http_communicator, params%dir_target//'/spprojs_completed', '**** SIMPLE_GEN_PICKREFS USER STOP ****')
            ! movie watcher init
            project_buff = stream_watcher(LONGTIME, params%dir_target//'/'//DIR_STREAM_COMPLETED, spproj=.true., nretries=10)
            ! import at least params%nmics micrographs - default 100
            write(logfhandle, '(A,I6,A)') '>>> IMPORTING AT LEAST', params%nmics, ' MICROGRAPHS'
            ! send generate pickrefs display info to gui
            call send_meta(string('waiting for micrographs'))
            call send_meta2D(string('waiting for particles'), box_in_pix)
            call micimporter( params%nmics )
            call send_meta(string('picking particles'))                              
            ! segmentation-based picking
            nmics = params%nmics ! because of Susan's bug report
            call segdiampick_mics(spproj, params%pcontrast, nmics, params%moldiam_max, box_in_pix, mskdiam_estimate)
            call send_meta(string('extracting particles'))
            if(spproj%os_mic%isthere('thumb_den') .and. spproj%os_mic%isthere('xdim') .and. spproj%os_mic%isthere('ydim') &
            .and. spproj%os_mic%isthere('smpd') .and. spproj%os_mic%isthere('boxfile') ) then
                i_max = min(spproj%os_mic%get_noris(), NTHUMB_MAX)
                do iori=1, i_max
                    ithumb = spproj%os_mic%get_noris() - i_max + iori
                    call send_micrograph_meta(iori, i_max, ithumb)
                end do
            endif     
            ! extract
            call cline_extract%set('prg',        'extract')
            call cline_extract%set('mkdir',           'no')
            call cline_extract%set('nparts',   params%nthr)
            call cline_extract%set('nthr',               1)
            call cline_extract%set('projfile',    projfile)
            call xextract%execute(cline_extract)
            call spproj%read(projfile)
            ! send generate pickrefs display info to gui
            call send_meta(string('complete'))
            call send_meta2D(string('classifying particles'), box_in_pix) 
            ! 2D analysis
            nptcls = spproj%os_ptcl2D%get_noris()
            ncls   = min(NCLS_MAX,max(NCLS_MIN,nptcls/params%nptcls_per_cls))
            nthr2D = max(1,floor(real(params%nthr)/4.))
            call cline_abinitio2D%set('prg',               'abinitio2D')
            call cline_abinitio2D%set('mkdir',                     'no')
            call cline_abinitio2D%set('ncls',                      ncls)
            call cline_abinitio2D%set('sigma_est',             'global')
            call cline_abinitio2D%set('center',                   'yes')
            call cline_abinitio2D%set('autoscale',                'yes')
            call cline_abinitio2D%set('lpstop',                  LPSTOP)
            call cline_abinitio2D%set('mskdiam',       mskdiam_estimate)
            call cline_abinitio2D%set('nthr',                    nthr2D)
            call cline_abinitio2D%set('nparts',                NPARTS2D)
            call cline_abinitio2D%set('projfile',              projfile)
            call xabinitio2D%execute(cline_abinitio2D)
            call send_meta2D(string('shape ranking particles'), box_in_pix)
            ! shape rank cavgs
            call cline_shape_rank%set('nthr',               params%nthr)
            call cline_shape_rank%set('projfile',              projfile)
            call xshape_rank%execute(cline_shape_rank)
            call spproj%read(projfile)
            call spproj%shape_ranked_cavgs2jpg(cavg_inds, string("shape_ranked_")//int2str(params%nmics)//JPG_EXT,&
            &xtiles, ytiles, mskdiam_px=ceiling(mskdiam_estimate * spproj%get_smpd()))
            call spproj%get_cavgs_stk(cavgsstk, ncls_stk, smpd_stk)
            ! send generate pickrefs display info to gui
            xtile = 0
            ytile = 0
            n_non_zero = spproj%os_ptcl2D%count_state_gt_zero()
            call http_gen_pickrefs_communicator%add_to_json("particles_accepted",  n_non_zero)
            call http_gen_pickrefs_communicator%add_to_json("particles_rejected",  spproj%os_ptcl2D%get_noris() - n_non_zero)
            call json%remove(latest_cls2D, destroy=.true.)
            call json%create_array(latest_cls2D, "latest_cls2D")
            call http_gen_pickrefs_communicator%add_to_json(latest_cls2D)
            call http_gen_pickrefs_communicator%add_to_json("user_input", .true.)
            call meta_opening2D%set_user_input(.true.)
            call send_meta2D(string('waiting for user selection'), box_in_pix)
            if(allocated(cavg_inds)) then
                do i=0, size(cavg_inds) - 1
                    call add_cls2D_to_json(trim(CWD_GLOB) // '/' // "shape_ranked_" // int2str(params%nmics) // JPG_EXT,&
                        cavgsstk%to_char(),&
                        cavg_inds(i + 1),&
                        xtile * (100.0 / (xtiles - 1)),&
                        ytile * (100.0 / (ytiles - 1)),&
                        100 * ytiles,&
                        100 * xtiles,&
                        scale=box_in_pix * spproj%get_smpd(),&
                        pop=spproj%os_cls2D%get_int(cavg_inds(i + 1), 'pop'),&
                        res=spproj%os_cls2D%get(cavg_inds(i + 1), 'res')&
                    )
                    xtile = xtile + 1
                    if(xtile .eq. xtiles) then
                        xtile = 0
                        ytile = ytile + 1
                    endif
                end do
            endif
            ! wait for user interaction
            write(logfhandle, *) ">>> WAITING FOR USER TO SELECT REFERENCES"
            do 
                call send_jobstats()
                call http_gen_pickrefs_communicator%get_json_arg('increase_nmics', increase_nmics, found)
                if (found) then
                    if(increase_nmics) then
                        increase_nmics = .false.
                        call cline%set('nmics', params%nmics + NMICS_DELTA)
                        params%nmics = params%nmics + NMICS_DELTA
                        call cleanup_previous
                        call meta_opening2D%set_user_input(.false.)
                        call send_meta2D(string('waiting for particles'), box_for_extract)
                        restart_requested = .true.
                        exit  ! exit inner wait-loop to restart outer loop
                    endif
                endif
                call http_gen_pickrefs_communicator%get_json_arg('final_selection', final_selection, found)
                if(found) then
                    call http_gen_pickrefs_communicator%get_json_arg('final_selection_source', buffer, found)
                    if(found) then
                        final_selection_source = buffer
                        call process_selected_refs(params, final_selection_source, spproj%get_smpd(), final_selection, mskdiam_estimate, box_for_pick, box_for_extract, xtiles, ytiles)
                        call http_gen_pickrefs_communicator%add_to_json("mask_diam", nint(mskdiam_estimate))
                        call http_gen_pickrefs_communicator%add_to_json("mskscale",  dble(box_for_extract * spproj%get_smpd()))
                        call http_gen_pickrefs_communicator%add_to_json("box_size",  box_for_extract)
                        call meta_opening2D%set_user_input(.false.)
                        call send_meta2D(string('applying user selection'), box_for_extract)
                        xtile = 0
                        ytile = 0
                        n_non_zero = 0
                        do i=1, size(final_selection)
                            call add_selected_refs_to_json(trim(CWD_GLOB) // '/' // STREAM_SELECTED_REFS // JPG_EXT,&
                                &xtile * (100.0 / (xtiles - 1)),&
                                &ytile * (100.0 / (ytiles - 1)),&
                                &100 * ytiles,&
                                &100 * xtiles,&
                                &pop=spproj%os_cls2D%get_int(final_selection(i), 'pop'),&
                                &res=spproj%os_cls2D%get(final_selection(i), 'res')&
                            )
                            n_non_zero = n_non_zero + spproj%os_cls2D%get_int(final_selection(i), 'pop')
                            xtile = xtile + 1
                            if(xtile .eq. xtiles) then
                                xtile = 0
                                ytile = ytile + 1
                            endif
                        end do
                        call http_gen_pickrefs_communicator%add_to_json("particles_accepted",  n_non_zero)
                        call http_gen_pickrefs_communicator%add_to_json("particles_rejected",  spproj%os_ptcl2D%get_noris() - n_non_zero)
                        if( allocated(buffer) ) deallocate(buffer)               
                    endif
                    exit
                endif
                call sleep(WAITTIME) ! may want to increase as 3s default
            enddo
            if(restart_requested) cycle    ! repeat the restartable block
            exit                           ! normal completion of restart loop
        end do ! end restart loop
        call send_jobstats()
        ! termination
        call http_communicator%update_json("stage", "terminating", found)
        call http_gen_pickrefs_communicator%update_json("stage", "terminating", found)
        call send_jobstats()
        if( allocated(projects)  ) deallocate(projects)
        if( allocated(cavg_inds) ) deallocate(cavg_inds)
        ! add optics
        if(cline%defined('optics_dir')) then
            optics_map_id = get_latest_optics_map_id(params%optics_dir)
            if(optics_map_id .gt. 0) then
                mapfileprefix = params%optics_dir//'/'//OPTICS_MAP_PREFIX//int2str(optics_map_id)
                call spproj%import_optics_map(mapfileprefix)
            endif
        endif
        ! write project and star files (just in case you want ot import these particles/micrographs elsewhere)
        call spproj%write
        call spproj%write_mics_star(string("micrographs.star"))
        call spproj%write_ptcl2D_star(string("particles.star"))
        ! cleanup
        call spproj%kill
        ! end gracefully
        call http_communicator%term()
        call http_gen_pickrefs_communicator%term()
        call simple_end('**** SIMPLE_GEN_PICKREFS NORMAL STOP ****')

        contains

            subroutine send_jobstats()
                call http_communicator%send_jobstats()
                call http_gen_pickrefs_communicator%send_jobstats()
            end subroutine send_jobstats

            subroutine background_heartbeats()
                call http_communicator%background_heartbeat()
                call http_gen_pickrefs_communicator%background_heartbeat()
            end subroutine background_heartbeats

            subroutine join_background_heartbeats()
                call http_communicator%join_background_heartbeat()
                call http_gen_pickrefs_communicator%join_background_heartbeat()
                if(http_communicator%exit_status() .or. http_gen_pickrefs_communicator%exit_status()) then
                    write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                    call spproj%kill
                    call http_communicator%term()
                    call http_gen_pickrefs_communicator%term()
                    call simple_end('**** SIMPLE_GEN_PICKREFS USER STOP ****')
                    call EXIT(0)
                endif
            end subroutine join_background_heartbeats

            subroutine cleanup_previous()
                type(string), allocatable :: list(:)
                integer :: i
                write(logfhandle,'(A)')'>>> CLEANING UP'
                call spproj%kill
                call cline%set('projfile', '')
                call simple_list_files('*', list)
                do i=1, size(list)
                    call del_file(list(i))
                enddo
              !  call cline%set('projname', PROJNAME_GEN_PICKREFS)
               ! call cline%set('projfile', PROJFILE_GEN_PICKREFS)
                call spproj%update_projinfo(cline)
                call spproj%update_compenv(cline)
                call spproj%write
            end subroutine cleanup_previous

            subroutine micimporter( nmics )
                integer, intent(in) :: nmics
                integer :: n_imported, n_new_oris, iproj, iori, imic
                n_imported= 0
                do
                    if( file_exists(TERM_STREAM) .or. http_communicator%exit_status()) then
                        ! termination
                        write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                        call spproj%kill
                        call qsys_cleanup(params)
                        call http_communicator%term()
                        call http_gen_pickrefs_communicator%term()
                        call simple_end('**** SIMPLE_GEN_PICKREFS USER STOP ****')
                        call EXIT(0)
                    endif
                    ! detection of new projects
                    call project_buff%watch( nprojects, projects, max_nmovies=50 )
                    ! append projects to processing stack
                    if( nprojects > 0 )then
                        n_imported= spproj%os_mic%get_noris()
                        if(n_imported> 0) then
                            n_new_oris  =  n_imported+ nprojects * STREAM_NMOVS_SET
                            call spproj%os_mic%reallocate(n_new_oris)
                        else
                            n_new_oris = nprojects * STREAM_NMOVS_SET
                            call spproj%os_mic%new(n_new_oris, .false.)
                        end if
                        do iproj = 1, nprojects 
                            call project_buff%add2history(projects(iproj)) ! I got this one, so please don't give it to me again
                            call spproj_part%read(projects(iproj))
                            do iori = 1, STREAM_NMOVS_SET
                                n_imported = n_imported + 1
                                ! set state=0 mics to state=-1
                                if(spproj_part%os_mic%get(iori, 'state') < 1) call spproj_part%os_mic%set(iori, 'state', -1)
                                call spproj%os_mic%transfer_ori(n_imported, spproj_part%os_mic, iori)
                            end do
                            call spproj_part%kill()
                        enddo
                        write(logfhandle,'(A,I4,A,A)')'>>> ', nprojects * STREAM_NMOVS_SET, ' NEW MICROGRAPHS DETECTED; ', cast_time_char(simple_gettime())
                    else
                        call sleep(WAITTIME) ! may want to increase as 3s default
                    endif
                    ! micrograph rejection
                    do imic = 1,spproj%os_mic%get_noris()
                        if( spproj%os_mic%get(imic, 'state') < 0 ) cycle
                        call spproj%os_mic%set(imic, 'state', 1) ! set all to state 1 pre rejection
                        if( spproj%os_mic%isthere(imic, 'ctfres') ) then
                            if( spproj%os_mic%get(imic,'ctfres') > (params%ctfresthreshold-0.001) ) call spproj%os_mic%set(imic, 'state', 0)
                        end if
                        if( spproj%os_mic%isthere(imic, 'icefrac') ) then
                            if( spproj%os_mic%get(imic,'icefrac') > (params%icefracthreshold-0.001) ) call spproj%os_mic%set(imic, 'state', 0)
                        end if
                        if( spproj%os_mic%isthere(imic, 'astig') ) then
                            if( spproj%os_mic%get(imic,'astig') > (params%astigthreshold-0.001) ) call spproj%os_mic%set(imic, 'state', 0)
                        end if
                    enddo
                    ! update thresholds if sent from gui
            !        call update_user_params(params, cline, http_communicator)
                    if( spproj%os_mic%count_state_gt_zero() >= nmics ) then
                        ! set state=-1 mics back to 1
                        do imic = 1,spproj%os_mic%get_noris()
                            if( spproj%os_mic%get(imic, 'state') < 0 ) call spproj%os_mic%set(imic, 'state', 0)
                        enddo
                        return
                    endif
                end do
            end subroutine micimporter

            subroutine communicator_init_initial_picking()
                call http_communicator%add_to_json("stage",       "initialising")
                call http_communicator%add_to_json("micrographs_imported",     0)
                call http_communicator%add_to_json("micrographs_accepted",     0)
                call http_communicator%add_to_json("micrographs_rejected",     0)
                call http_communicator%add_to_json("particles_extracted",      0)
                call http_communicator%add_to_json("particles_per_mic",        0)
                call http_communicator%add_to_json("user_input",         .false.)
                call http_communicator%add_to_json("last_micrograph_imported", "")
                call json%create_array(latest_picked_micrographs, "latest_picked_micrographs")
                call http_communicator%add_to_json(latest_picked_micrographs)
            end subroutine communicator_init_initial_picking

            subroutine communicator_gen_pickrefs_init()
                call http_gen_pickrefs_communicator%add_to_json("stage",              "initialising")
                call http_gen_pickrefs_communicator%add_to_json("particles_imported", 0)
                call http_gen_pickrefs_communicator%add_to_json("particles_accepted", 0)
                call http_gen_pickrefs_communicator%add_to_json("particles_rejected", 0)
                call http_gen_pickrefs_communicator%add_to_json("mask_diam",          0)
                call http_gen_pickrefs_communicator%add_to_json("box_size",           0)
                call http_gen_pickrefs_communicator%add_to_json("mskscale",           dble(0.0))
                call http_gen_pickrefs_communicator%add_to_json("user_input",         .false.)
                call json%create_array(latest_cls2D, "latest_cls2D")
                call http_gen_pickrefs_communicator%add_to_json(latest_cls2D)
                call json%create_array(selected_refs, "selected_references")
                call http_gen_pickrefs_communicator%add_to_json(selected_refs)
            end subroutine communicator_gen_pickrefs_init

            subroutine add_micrograph_to_json( path, xdim, ydim, boxfile_path )
                character(*),     intent(in)  :: path, boxfile_path
                integer,          intent(in)  :: xdim, ydim
                type(nrtxtfile)               :: boxfile
                type(json_value), pointer     :: micrograph, boxes, box
                real,             allocatable :: boxdata(:,:)
                integer                       :: i, x, y
                call json%create_object(micrograph, "")
                call json%add(micrograph, "path",    path)
                call json%add(micrograph, "xdim"   , xdim)
                call json%add(micrograph, "ydim",    ydim)
                call json%create_array(boxes, "boxes")
                call boxfile%new(string(boxfile_path), 1)
                allocate(boxdata(boxfile%get_nrecs_per_line(), boxfile%get_ndatalines()))
                if(boxfile%get_nrecs_per_line() >= 4) then
                    do i=1, boxfile%get_ndatalines()
                        call boxfile%readNextDataLine(boxdata(:,i))
                        call json%create_object(box, "")
                        x = nint(boxdata(1,i) + boxdata(3,i)/2)
                        y = nint(boxdata(2,i) + boxdata(4,i)/2)
                        call json%add(box, "x",    x)
                        call json%add(box, "y",    y)
                        call json%add(boxes, box)
                    enddo
                endif
                call boxfile%kill()
                deallocate(boxdata)
                call json%add(micrograph, boxes)
                call json%add(latest_picked_micrographs, micrograph)
            end subroutine add_micrograph_to_json
            
            subroutine add_cls2D_to_json( path, mrcpath, mrc_idx, spritex, spritey, spriteh, spritew, res, pop, scale )
                character(*),      intent(in) :: path, mrcpath
                real,              intent(in) :: spritex, spritey
                integer,           intent(in) :: spriteh, spritew, mrc_idx
                integer, optional, intent(in) :: pop
                real,    optional, intent(in) :: res, scale
                type(json_value),  pointer    :: template
                call json%create_object(template, "")
                call json%add(template, "path",     path)
                call json%add(template, "mrcpath",  mrcpath)
                call json%add(template, "mrcidx",   mrc_idx)
                call json%add(template, "spritex",  dble(spritex))
                call json%add(template, "spritey",  dble(spritey))
                call json%add(template, "spriteh",  spriteh)
                call json%add(template, "spritew",  spritew)
                if(present(scale)) call json%add(template, "mskscale", dble(scale))
                if(present(res))   call json%add(template, "res",      dble(res))
                if(present(pop))   call json%add(template, "pop",      pop)
                call json%add(latest_cls2D, template)
             end subroutine add_cls2D_to_json

            subroutine send_meta( my_stage )
                type(string), intent(in) :: my_stage
                call meta_initial_picking%set(stage                = my_stage,                           & 
                                              micrographs_imported = spproj%os_mic%get_noris(),          &
                                              micrographs_accepted = spproj%os_mic%count_state_gt_zero(),&
                                              particles_extracted  = spproj%os_ptcl2D%get_noris())
                if( meta_initial_picking%assigned() ) then
                    if( mq_stream_master_in%is_active() ) then
                        write(logfhandle,*) "SENDING"
                        call meta_initial_picking%serialise(meta_buffer)
                        call mq_stream_master_in%send(meta_buffer)
                    endif
                endif
            end subroutine send_meta

            subroutine send_meta2D( my_stage, box_size )
                type(string), intent(in) :: my_stage
                integer,      intent(in) :: box_size
                real                     :: smpd = 0.0
                if( spproj%os_ptcl2D%get_noris() > 0 ) smpd = spproj%get_smpd()
                call meta_opening2D%set(stage              = my_stage,                              & 
                                        particles_imported = spproj%os_ptcl2D%get_noris(),          &
                                        particles_accepted = spproj%os_ptcl2D%count_state_gt_zero(),&
                                        mask_diam          = nint(mskdiam_estimate),                &
                                        mask_scale         = box_size * smpd,                       &
                                        box_size           = box_size)
                if( meta_opening2D%assigned() ) then
                    if( mq_stream_master_in%is_active() ) then
                        write(logfhandle,*) "SENDING"
                        call meta_opening2D%serialise(meta_buffer)
                        call mq_stream_master_in%send(meta_buffer)
                    endif
                endif
            end subroutine send_meta2D

            subroutine send_micrograph_meta( my_iori, my_i_max, my_ithumb )
                integer,         intent(in) :: my_iori, my_i_max, my_ithumb
                type(nrtxtfile)             :: my_boxfile
                type(string)                :: my_boxpath
                real,           allocatable :: boxdata(:)
                integer                     :: my_i, my_x, my_y, my_nrecs, my_xdim, my_ydim
                call meta_micrograph%set(path   = spproj%os_mic%get_str(my_ithumb, "thumb"),&
                                         dfx    = spproj%os_mic%get(my_ithumb,       "dfx"),&
                                         dfy    = spproj%os_mic%get(my_ithumb,       "dfy"),&
                                         ctfres = spproj%os_mic%get(my_ithumb,    "ctfres"),&
                                         i_max  = my_i_max                                 ,&
                                         i      = my_iori)
                call meta_micrograph%clear_coordinates()
                my_boxpath = spproj%os_mic%get_str(my_ithumb, "boxfile")
                if( my_boxpath%strlen() > 0 ) then
                    if( file_exists(my_boxpath) ) then
                        call my_boxfile%new(my_boxpath, 1)
                        my_nrecs = my_boxfile%get_nrecs_per_line()
                        my_xdim  = spproj%os_mic%get(my_ithumb, "xdim")
                        my_ydim  = spproj%os_mic%get(my_ithumb, "ydim")
                        if(my_nrecs >= 4) then
                            allocate(boxdata(my_nrecs))
                            do my_i=1, my_boxfile%get_ndatalines()
                                call my_boxfile%readNextDataLine(boxdata(:))
                                my_x = nint(boxdata(1) + boxdata(3)/2)
                                my_y = nint(boxdata(2) + boxdata(4)/2)
                                call meta_micrograph%set_coordinate(my_i, my_x, my_y, my_xdim, my_ydim)
                            enddo
                            deallocate(boxdata)
                        endif
                        call my_boxfile%kill()
                        if( meta_micrograph%assigned() .and. mq_stream_master_in%is_active() ) then
                            call meta_micrograph%serialise(meta_buffer)
                            call mq_stream_master_in%send(meta_buffer)
                        endif
                    endif
                endif
            end subroutine send_micrograph_meta

             subroutine add_selected_refs_to_json( path, spritex, spritey, spriteh, spritew, res, pop )
                character(*),      intent(in) :: path
                real,              intent(in) :: spritex, spritey
                integer,           intent(in) :: spriteh, spritew
                integer, optional, intent(in) :: pop
                real,    optional, intent(in) :: res
                type(json_value),  pointer    :: template
                call json%create_object(template, "")
                call json%add(template, "path",    path)
                call json%add(template, "spritex", dble(spritex))
                call json%add(template, "spritey", dble(spritey))
                call json%add(template, "spriteh", spriteh)
                call json%add(template, "spritew", spritew)
                if(present(res)) call json%add(template, "res",     dble(res))
                if(present(pop)) call json%add(template, "pop",     pop)
                call json%add(selected_refs, template)
            end subroutine add_selected_refs_to_json
            
            subroutine sigterm_handler()
                write(logfhandle, '(A)') 'SIGTERM RECEIVED'
                l_terminate = .true.
            end subroutine sigterm_handler

    end subroutine exec_stream_p03_opening2D

end module simple_stream_p03_opening2D_new
