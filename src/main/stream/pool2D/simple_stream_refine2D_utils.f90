!@descr: utilities for running 2D clustering in the stream
module simple_stream_refine2D_utils
use simple_core_module_api
use simple_stream2D_state
use simple_class_frcs,           only: class_frcs
use simple_cmdline,              only: cmdline
use simple_commanders_refine2D, only: commander_rank_cavgs
use simple_image,                only: image
use simple_parameters,           only: parameters
use simple_qsys_funs,            only: qsys_cleanup
use simple_rec_list,             only: rec_list
use simple_sp_project,           only: sp_project
use simple_stack_io,             only: stack_io
use simple_starproject,          only: starproject
use simple_starproject_stream,   only: starproject_stream
use simple_optics_maps,          only: import_latest_optics_map
use simple_syslib,               only: get_current_rss_bytes, get_peak_rss_bytes
implicit none

public :: cleanup_root_folder
public :: setup_downscaling
public :: terminate_chunks
public :: terminate_stream2D
public :: test_repick
public :: tidy_2Dstream_iter
public :: write_project_stream2D
public :: write_pool_snapshot
public :: write_repick_refs
public :: build_pool_publication
public :: publish_pool_state
public :: delete_pool_publication
private
#include "simple_local_flags.inc"

type(starproject_stream) :: starproj_stream
logical, parameter       :: DEBUG_HERE = .false.
integer, allocatable     :: repick_selection(:)  ! selection for selecting classes for re-picking
integer(timer_int_kind)  :: t
integer                  :: repick_iteration = 0 ! iteration to select classes from for re-picking

contains

    ! Deselects the particles of the classes @p selection does not name (all classes are kept
    ! when it names none in range).
    subroutine apply_snapshot_selection( snapshot_projfile, selection )
        type(sp_project), intent(inout) :: snapshot_projfile
        integer,          intent(in)    :: selection(:)
        logical, allocatable :: cls_mask(:)
        integer              :: nptcls_rejected, ncls_rejected, iptcl
        integer              :: icls, jcls, i
        if( snapshot_projfile%os_cls2D%get_noris() == 0 ) return
        allocate(cls_mask(ncls_glob), source=.false.)
        do i = 1, size(selection)
            icls = selection(i)
            if( icls < 1 ) cycle
            if( icls > ncls_glob ) cycle
            cls_mask(icls) = .true.
        enddo
        if( count(cls_mask) == 0 ) return
        ncls_rejected   = 0
        do icls = 1,ncls_glob
            if(cls_mask(icls) ) cycle
            nptcls_rejected = 0
            !$omp parallel do private(iptcl,jcls) reduction(+:nptcls_rejected) proc_bind(close)
            do iptcl = 1,snapshot_projfile%os_ptcl2D%get_noris()
                if( snapshot_projfile%os_ptcl2D%get_state(iptcl) == 0 )cycle
                jcls = snapshot_projfile%os_ptcl2D%get_class(iptcl)
                if( jcls == icls )then
                    call snapshot_projfile%os_ptcl2D%reject(iptcl)
                    call snapshot_projfile%os_ptcl2D%delete_2Dclustering(iptcl)
                    nptcls_rejected = nptcls_rejected + 1
                endif
            enddo
            !$omp end parallel do
            if( nptcls_rejected > 0 )then
                ncls_rejected = ncls_rejected + 1
                call snapshot_projfile%os_cls2D%set_state(icls,0)
                call snapshot_projfile%os_cls2D%set(icls,'pop',0.)
                call snapshot_projfile%os_cls2D%set(icls,'corr',-1.)
                call snapshot_projfile%os_cls2D%set(icls,'prev_pop_even',0.)
                call snapshot_projfile%os_cls2D%set(icls,'prev_pop_odd', 0.)
                write(logfhandle,'(A,I6,A,I4)')'>>> USER REJECTED FROM SNAPSHOT: ',nptcls_rejected,' PARTICLE(S) IN CLASS ',icls
            endif
        enddo
    end subroutine apply_snapshot_selection

    ! Remove previous files from folder to restart
    subroutine cleanup_root_folder( all )
        logical, optional, intent(in)  :: all
        type(string), allocatable :: files(:), folders(:)
        integer :: i
 !       call simple_rmdir(DIR_SNAPSHOT) ! need snapshots keeping 
        call del_file(USER_PARAMS2D)
        call del_file(POOL_PROJFILE)
        call del_file(POOL_DISTR_EXEC_FNAME)
        call del_file(TERM_STREAM)
        call simple_list_files_regexp(string('.'), '\.mrc$|\.mrcs$|\.txt$|\.star$|\.eps$|\.jpeg$|\.jpg$|\.dat$|\.bin$', files)
        if( allocated(files) )then
            do i = 1,size(files)
                call del_file(files(i))
            enddo
        endif
        folders = simple_list_dirs('.')
        if( allocated(folders) )then
            do i = 1,size(folders)
                if( folders(i)%has_substr(DIR_CHUNK) ) call simple_rmdir(folders(i))
            enddo
        endif
        if( present(all) )then
            if( all )then
                call del_file(POOL_LOGFILE)
                call del_file(REFINE2D_FINISHED)
                call del_file('simple_script_single')
            endif
        endif
    end subroutine cleanup_root_folder

    subroutine debug_print( str )
        character(len=*), intent(in) :: str
        if( DEBUG_HERE )then
            write(logfhandle,*) trim(str)
            call flush(logfhandle)
        endif
    end subroutine debug_print

    integer function get_box()
        get_box = int(pool_proj%get_box())
    end function get_box

    integer function get_boxa()
        if(pool_proj%os_stk%get_noris() .gt. 0) then
            get_boxa = ceiling(real(pool_proj%get_box()) * pool_proj%get_smpd())
        else
            get_boxa = 0
        end if
    end function get_boxa

    ! For ranking class-averages
    subroutine rank_cavgs
        type(commander_rank_cavgs) :: xrank_cavgs
        type(cmdline)              :: cline_rank_cavgs
        type(string)               :: refs_ranked, stk
        refs_ranked = add2fbody(refs_glob, MRC_EXT ,'_ranked')
        call cline_rank_cavgs%set('projfile', orig_projfile)
        stk = string(POOL_DIR)//refs_glob
        call cline_rank_cavgs%set('stk',            stk)
        call cline_rank_cavgs%set('outstk', refs_ranked)
        call xrank_cavgs%execute(cline_rank_cavgs)
        call cline_rank_cavgs%kill
    end subroutine rank_cavgs

    ! Final rescaling of references
    subroutine rescale_cavgs( params, src, dest )
        class(parameters), intent(in) :: params
        class(string), intent(in) :: src, dest
        integer, allocatable :: cls_pop(:)
        type(image)    :: img, img_pad
        type(stack_io) :: stkio_r, stkio_w
        type(string)   :: dest_here
        integer        :: ldim(3),icls, ncls_here
        call debug_print('in rescale_cavgs '//src%to_char()//' -> '//dest%to_char())
        if(src == dest)then
            dest_here = 'tmp_cavgs.mrc'
        else
            dest_here = dest
        endif
        call img%new([pool_dims%box,pool_dims%box,1],pool_dims%smpd)
        call img_pad%new([params%box,params%box,1],params%smpd)
        cls_pop = nint(pool_proj%os_cls2D%get_all('pop'))
        call find_ldim_nptcls(src,ldim,ncls_here)
        call stkio_r%open(src, pool_dims%smpd, 'read', bufsz=ncls_here)
        call stkio_r%read_whole
        call stkio_w%open(dest_here, params%smpd, 'write', box=params%box, bufsz=ncls_here)
        do icls = 1,ncls_here
            if( cls_pop(icls) > 0 )then
                call img%zero_and_unflag_ft
                call stkio_r%get_image(icls, img)
                call img%fft
                call img%pad(img_pad, backgr=0., antialiasing=.false.)
                call img_pad%ifft
            else
                img_pad = 0.
            endif
            call stkio_w%write(icls, img_pad)
        enddo
        call stkio_r%close
        call stkio_w%close
        if ( src == dest ) call simple_rename('tmp_cavgs.mrc',dest)
        call img%kill
        call img_pad%kill
        call debug_print('end rescale_cavgs')
    end subroutine rescale_cavgs

    subroutine set_dimensions( params )
        class(parameters), intent(inout) :: params
        call setup_downscaling( params )
        pool_dims%smpd  = params%smpd_crop
        pool_dims%box   = params%box_crop
        pool_dims%boxpd = 2 * round2even(KBALPHA * real(params%box_crop/2)) ! logics from parameters
        pool_dims%msk   = params%msk_crop
        chunk_dims = pool_dims ! chunk & pool have the same dimensions to start with
        ! Scaling-related command lines update
        call cline_refine2D_chunk%set('smpd_crop', chunk_dims%smpd)
        call cline_refine2D_chunk%set('box_crop',   chunk_dims%box)
        call cline_refine2D_chunk%set('msk_crop',   chunk_dims%msk)
        call cline_refine2D_chunk%set('box',        params%box)
        call cline_refine2D_chunk%set('smpd',       params%smpd)
        call cline_refine2D_pool%set('smpd_crop',   pool_dims%smpd)
        call cline_refine2D_pool%set('box_crop',    pool_dims%box)
        call cline_refine2D_pool%set('msk_crop',    pool_dims%msk)
        call cline_refine2D_pool%set('box',         params%box)
        call cline_refine2D_pool%set('smpd',        params%smpd)
    end subroutine set_dimensions

    ! private routine for resolution-related updates to command-lines
    subroutine set_resolution_limits( params )
        class(parameters), intent(inout) :: params
        lpstart = max(lpstart, 2.0*params%smpd_crop)
        if( l_no_chunks )then
            params%lpstop = lpstop
        else
            if( master_cline%defined('lpstop') )then
                params%lpstop = max(2.0*params%smpd_crop,params%lpstop)
            else
                params%lpstop = 2.0*params%smpd_crop
            endif
            call cline_refine2D_chunk%delete('lp')
            call cline_refine2D_chunk%set('lpstart', lpstart)
            call cline_refine2D_chunk%set('lpstop',   lpstart)
        endif
        call cline_refine2D_pool%set('lpstart',   lpstart)
        call cline_refine2D_pool%set('lpstop',    params%lpstop)
        if( .not.master_cline%defined('cenlp') )then
            call cline_refine2D_chunk%set('cenlp', lpcen)
            call cline_refine2D_pool%set( 'cenlp', lpcen)
        else
            call cline_refine2D_chunk%set('cenlp', params%cenlp)
            call cline_refine2D_pool%set( 'cenlp', params%cenlp)
        endif
        ! Will use resolution update scheme from solve2D
        if( .not.l_no_chunks )then
            if( master_cline%defined('lpstop') )then
                ! already set above
            else
                call cline_refine2D_chunk%delete('lpstop')
            endif
        endif
        write(logfhandle,'(A,F5.1)') '>>> POOL STARTING LOW-PASS LIMIT (IN A): ', lpstart
        write(logfhandle,'(A,F5.1)') '>>> POOL   HARD RESOLUTION LIMIT (IN A): ', params%lpstop
        write(logfhandle,'(A,F5.1)') '>>> CENTERING     LOW-PASS LIMIT (IN A): ', lpcen
    end subroutine set_resolution_limits

    ! Determines dimensions for downscaling
    subroutine setup_downscaling( params )
        class(parameters), intent(inout) :: params
        real    :: SMPD_TARGET = MAX_SMPD  ! target sampling distance
        real    :: smpd, scale_factor
        integer :: box
        if( params%box == 0 ) THROW_HARD('FATAL ERROR')
        scale_factor          = 1.0
        params%smpd_crop = params%smpd
        params%box_crop  = params%box
        if( l_scaling .and. params%box >= CHUNK_MINBOXSZ )then
            call autoscale(params%box, params%smpd, SMPD_TARGET, box, smpd, scale_factor, minbox=CHUNK_MINBOXSZ)
            l_scaling = box < params%box
            if( l_scaling )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',box
                params%smpd_crop = smpd
                params%box_crop  = box
            endif
        endif
        params%msk_crop = round2even(params%mskdiam / params%smpd_crop / 2.)
    end subroutine setup_downscaling

    ! ends chunks processing
    subroutine terminate_chunks( params )
        class(parameters), intent(in) :: params
        integer :: ichunk
        do ichunk = 1,params%nchunks
            call chunks(ichunk)%terminate_chunk
        enddo
    end subroutine terminate_chunks

    ! ends processing, generates project & cleanup
    subroutine terminate_stream2D( params, project_list, optics_dir)
        class(parameters),         intent(inout) :: params
        class(rec_list), optional, intent(inout) :: project_list
        class(string),   optional, intent(in)    :: optics_dir
        integer      :: ipart, lastmap
        call terminate_chunks( params )
        if( pool_iter <= 0 )then
            ! no 2D yet
            call write_raw_project
        else
            if( .not.l_pool_available )then
                pool_iter = pool_iter-1 ! iteration pool_iter not complete so fall back on previous iteration
                if( pool_iter <= 0 )then
                    ! no 2D yet
                    call write_raw_project
                else
                    refs_glob = CAVGS_ITER_FBODY//int2str_pad(pool_iter,3)//MRC_EXT
                    ! tricking the asynchronous master process to come to a hard stop
                    call simple_touch(POOL_DIR//TERM_STREAM)
                    do ipart = 1,params%nparts_pool
                        call simple_touch(POOL_DIR//JOB_FINISHED_FBODY//int2str_pad(ipart,numlen))
                    enddo
                    call simple_touch(POOL_DIR//'CAVGASSEMBLE_FINISHED')
                endif
            endif
            if( pool_iter >= 1 )then
                call write_project_stream2D(params, write_star=.true., clspath=.true., optics_dir=optics_dir)
                call rank_cavgs
            endif
        endif
        ! cleanup
        call del_file(POOL_DIR//POOL_PROJFILE)
        call del_file(projfile4gui)
        if( .not. DEBUG_HERE )then
            call qsys_cleanup(params)
        endif

        contains

            ! no pool clustering performed, all available info is written down
            subroutine write_raw_project
                if( present(project_list) )then
                    if( project_list%size() > 0 )then
                        ! the groups of the newest optics map, when there is one
                        if( present(optics_dir) ) lastmap = import_latest_optics_map(pool_proj, optics_dir)
                        call pool_proj%projrecords2proj(project_list)
                        call starproj_stream%copy_micrographs_optics(pool_proj, verbose=DEBUG_HERE)
                        call starproj_stream%stream_export_micrographs(params, pool_proj, params%outdir, optics_set=.true.)
                        call starproj_stream%stream_export_particles_2D(params, pool_proj, params%outdir, optics_set=.true.)
                        call pool_proj%write(orig_projfile)
                    endif
                endif
            end subroutine write_raw_project

    end subroutine terminate_stream2D

    logical function test_repick()
        test_repick = allocated(repick_selection)
    end function test_repick

    ! Removes some unnecessary files
    !> Removes the files of pool iteration @p iter (none when it is below 1): its class averages
    !! (and their even and odd halves and JPEG), its class STAR file and its FRCs copy.
    subroutine tidy_2Dstream_iter( iter )
        integer, intent(in) :: iter
        type(string) :: prefix
        if( iter < 1 ) return
        prefix = POOL_DIR//CAVGS_ITER_FBODY//int2str_pad(iter,3)
        call del_file(prefix//JPG_EXT)
        call del_file(prefix//'_even'//MRC_EXT)
        call del_file(prefix//'_odd'//MRC_EXT)
        call del_file(prefix//MRC_EXT)
        call del_file(POOL_DIR//CLS2D_STARFBODY//'_iter'//int2str_pad(iter,3)//STAR_EXT)
        call del_file(string(POOL_DIR)//swap_suffix(FRCS_FILE,"_iter"//int2str_pad(iter, 3)//".bin",".bin"))
    end subroutine tidy_2Dstream_iter


    !> Writes the pool's project (the stage's final project): the class averages and FRCs at the
    !! native sampling, the newest optics map's groups (with @p optics_dir), ptcl3D prepared for
    !! 3D, the class STAR file and, with @p write_star, the micrograph and particle STAR files.
    subroutine write_project_stream2D( params, write_star, clspath, optics_dir )
        class(parameters), intent(inout) :: params
        logical,       optional, intent(in) :: write_star
        logical,       optional, intent(in) :: clspath
        class(string), optional, intent(in) :: optics_dir
        type(class_frcs)   :: frcs
        type(oris)         :: os_backup
        type(string)       :: projfile, projfname, cavgsfname, frcsfname, pool_refs
        integer            :: lastmap
        logical            :: l_write_star, l_clspath
        l_write_star = .false.
        l_clspath    = .false.
        if(present(write_star)) l_write_star = write_star
        if(present(clspath))    l_clspath    = clspath
        ! file naming
        projfname  = get_fbody(orig_projfile, METADATA_EXT, separator=.false.)
        cavgsfname = get_fbody(refs_glob, MRC_EXT, separator=.false.)
        frcsfname  = get_fbody(FRCS_FILE, BIN_EXT, separator=.false.)
        call pool_proj%projinfo%set(1,'projname', projfname)
        projfile   = projfname//METADATA_EXT
        call pool_proj%projinfo%set(1,'projfile', projfile)
        cavgsfname = cavgsfname//MRC_EXT
        frcsfname  = frcsfname//BIN_EXT
        pool_refs  = string(POOL_DIR)//refs_glob
        lastmap    = 0
        write(logfhandle,'(A,A,A,A)')'>>> WRITING PROJECT ', projfile%to_char(), ' AT: ',cast_time_char(simple_gettime())
        if( present(optics_dir) ) lastmap = import_latest_optics_map(pool_proj, optics_dir)
        if( l_scaling )then
            os_backup = pool_proj%os_cls2D
            call rescale_refs( params, cavgsfname )
            call pool_proj%os_out%kill
            call pool_proj%add_cavgs2os_out(cavgsfname, params%smpd, 'cavg', clspath=l_clspath)
            pool_proj%os_cls2D = os_backup
            call os_backup%kill
            ! rescale frcs
            call frcs%read(string(POOL_DIR)//FRCS_FILE)
            call frcs%pad(params%smpd, params%box)
            call frcs%write(frcsfname)
            call frcs%kill
            call pool_proj%add_frcs2os_out(frcsfname, 'frc2D')
        else
            call pool_proj%os_out%kill
            call pool_proj%add_cavgs2os_out(cavgsfname, params%smpd, 'cavg', clspath=l_clspath)
            call pool_proj%add_frcs2os_out(frcsfname, 'frc2D')
        endif
        ! the 3D field as the STAR files have it: 2D clustering removed, shifts kept
        pool_proj%os_ptcl3D = pool_proj%os_ptcl2D
        call pool_proj%os_ptcl3D%delete_2Dclustering
        call pool_proj%os_ptcl3D%clean_entry('updatecnt', 'sampled')
        call pool_proj%write(projfile)
        ! write starfiles
        call starproj%export_cls2D(pool_proj)
        if(l_write_star) then
            if( DEBUG_HERE ) t = tic()
            call starproj_stream%stream_export_micrographs(params, pool_proj, params%outdir, optics_set=.true.)
            if( DEBUG_HERE ) print *,'ms_export  : ', toc(t); call flush(6); t = tic()
            call starproj_stream%stream_export_particles_2D(params, pool_proj, params%outdir, optics_set=.true.)
            if( DEBUG_HERE ) print *,'ptcl_export  : ', toc(t); call flush(6)
        end if
        call pool_proj%os_ptcl3D%kill
        call pool_proj%os_cls2D%delete_entry('stk')

        contains

        ! rescale classes to original scale
        subroutine rescale_refs( params, cavgs_fname )
            class(parameters), intent(in) :: params
            class(string), intent(in) :: cavgs_fname
            type(string) :: source, destination
            call rescale_cavgs(params, pool_refs, cavgs_fname)
            source  = add2fbody(pool_refs, MRC_EXT, '_even')
            destination = add2fbody(cavgs_fname, MRC_EXT, '_even')
            call rescale_cavgs(params, source, destination)
            source  = add2fbody(pool_refs, MRC_EXT,'_odd')
            destination = add2fbody(cavgs_fname, MRC_EXT,'_odd')
            call rescale_cavgs(params, source, destination)
        end subroutine rescale_refs

    end subroutine write_project_stream2D

    !> Writes a snapshot of pool iteration @p iteration as @p projfile: the classes @p selection
    !! names (the others' particles deselected), with the class averages and FRCs of that
    !! iteration at the pool's sampling, the newest optics map's groups (when @p optics_dir is not
    !! empty) and micrograph and particle STAR files (@p starfile_base, optics group ids offset by
    !! @p optics_offset). The iteration is the current one or one of the history
    !! (POOL_NHISTORY). @p nptcls is the number of selected particles written; 0 when the
    !! iteration is no longer kept or its files are missing, and then nothing is written.
    !! @p cavgs_* report the selected class averages' sprite sheet for the GUI.
    subroutine write_pool_snapshot( iteration, selection, projfile, starfile_base, optics_dir, optics_offset,&
            &nptcls, cavgs_jpeg, cavgs_mrc, cavgs_ntilesx, cavgs_ntilesy, cavgs_idx, cavgs_pop, cavgs_res )
        integer,              intent(in)  :: iteration
        integer,              intent(in)  :: selection(:)
        class(string),        intent(in)  :: projfile, starfile_base, optics_dir
        integer,              intent(in)  :: optics_offset
        integer,              intent(out) :: nptcls
        type(string),         intent(out) :: cavgs_jpeg, cavgs_mrc
        integer,              intent(out) :: cavgs_ntilesx, cavgs_ntilesy
        integer, allocatable, intent(out) :: cavgs_idx(:), cavgs_pop(:)
        real,    allocatable, intent(out) :: cavgs_res(:)
        type(sp_project) :: snapshot_proj
        type(string)     :: dir, stk, frcs, cavgsfname, frcsfname
        real             :: smpd
        integer          :: ncls, islot, lastmap
        logical          :: l_found
        nptcls        = 0
        cavgs_jpeg    = ''
        cavgs_mrc     = ''
        cavgs_ntilesx = 0
        cavgs_ntilesy = 0
        allocate(cavgs_idx(0), cavgs_pop(0), cavgs_res(0))
        write(logfhandle,'(A,I4,A,A,A,A)') '>>> WRITING SNAPSHOT FROM ITERATION ', iteration, ': ', projfile%to_char(),&
            &' AT: ', cast_time_char(simple_gettime())
        call log_rss('snapshot/before copy')
        ! the iteration: the current one, or one the history keeps
        l_found = .false.
        if( iteration == pool_iter )then
            call snapshot_proj%copy(pool_proj)
            call snapshot_proj%add_frcs2os_out(string(POOL_DIR)//FRCS_FILE, 'frc2D')
            l_found = .true.
        else if( iteration >= 1 )then
            islot   = pool_history_slot(iteration)
            l_found = pool_history_iter(islot) == iteration
            if( l_found ) call snapshot_proj%copy(pool_proj_history(islot))
        endif
        ! and its files
        if( l_found )then
            call snapshot_proj%get_cavgs_stk(stk, ncls, smpd, fail=.false.)
            call snapshot_proj%get_frcs(frcs, 'frc2D', fail=.false.)
            l_found = ncls > 0
            if( l_found ) l_found = file_exists(stk) .and. file_exists(frcs)
            if( l_found ) l_found = file_exists(add2fbody(stk, MRC_EXT, '_even'))
            if( l_found ) l_found = file_exists(add2fbody(stk, MRC_EXT, '_odd'))
        endif
        if( .not. l_found )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> WARNING: ITERATION ', iteration, ' IS NOT KEPT (THE POOL KEEPS THE LAST ',&
                &POOL_NHISTORY, '); SNAPSHOT NOT WRITTEN'
            call snapshot_proj%kill
            return
        endif
        call log_rss('snapshot/after copy')
        dir = stemname(projfile)
        if( .not. file_exists(stemname(dir)) ) call simple_mkdir(stemname(dir))
        if( .not. file_exists(dir) )           call simple_mkdir(dir)
        if( optics_dir%strlen() > 0 ) lastmap = import_latest_optics_map(snapshot_proj, optics_dir)
        call apply_snapshot_selection(snapshot_proj, selection)
        call log_rss('snapshot/after selection')
        cavgsfname = dir//'/cavgs'//STK_EXT
        frcsfname  = dir//'/'//FRCS_FILE
        call simple_copy_file(stk, cavgsfname)
        call simple_copy_file(add2fbody(stk, MRC_EXT, '_even'), dir//'/cavgs_even'//STK_EXT)
        call simple_copy_file(add2fbody(stk, MRC_EXT, '_odd'),  dir//'/cavgs_odd'//STK_EXT)
        call simple_copy_file(frcs, frcsfname)
        call snapshot_proj%os_out%kill
        ! copied as they are, at the pool's (possibly cropped) sampling
        call snapshot_proj%add_cavgs2os_out(cavgsfname, smpd, 'cavg')
        call snapshot_proj%add_frcs2os_out(frcsfname, 'frc2D')
        call snapshot_proj%set_cavgs_thumb(projfile)
        call snapshot_cavgs_meta(snapshot_proj, cavgsfname, cavgs_jpeg, cavgs_mrc, cavgs_ntilesx, cavgs_ntilesy,&
            &cavgs_idx, cavgs_pop, cavgs_res)
        snapshot_proj%os_ptcl3D = snapshot_proj%os_ptcl2D
        call snapshot_proj%write(projfile)
        call snapshot_proj%write_mics_star(starfile_base//"_micrographs.star", optics_offset=optics_offset)
        call snapshot_proj%write_ptcl2D_star(starfile_base//"_particles.star", optics_offset=optics_offset)
        nptcls = snapshot_proj%os_ptcl2D%count_state_gt_zero()
        call log_rss('snapshot/after write')
        call snapshot_proj%kill
    end subroutine write_pool_snapshot

    ! The sprite sheet set_cavgs_thumb wrote for the selected class averages of @p proj
    ! (thumb, thumbnx, thumbny, thumbidx in its cls2D), so the GUI gets the snapshot's class
    ! averages without recomputing their geometry; @p cavgsfname is their stack.
    subroutine snapshot_cavgs_meta( proj, cavgsfname, jpeg_path, mrc_path, ntilesx, ntilesy, cls_idx, cls_pop, cls_res )
        class(sp_project),    intent(in)  :: proj
        class(string),        intent(in)  :: cavgsfname
        type(string),         intent(out) :: jpeg_path, mrc_path
        integer,              intent(out) :: ntilesx, ntilesy
        integer, allocatable, intent(out) :: cls_idx(:), cls_pop(:)
        real,    allocatable, intent(out) :: cls_res(:)
        integer :: iori, noris, n_selected, thumbidx
        jpeg_path = ''
        mrc_path  = ''
        ntilesx   = 0
        ntilesy   = 0
        noris     = proj%os_cls2D%get_noris()
        n_selected = 0
        do iori = 1,noris
            if( proj%os_cls2D%isthere(iori,'thumbidx') ) n_selected = n_selected + 1
        enddo
        allocate(cls_idx(n_selected), cls_pop(n_selected), source=0)
        allocate(cls_res(n_selected), source=0.0)
        if( noris == 0 ) return
        jpeg_path = proj%os_cls2D%get_str(1,'thumb')
        mrc_path  = cavgsfname
        ntilesx   = proj%os_cls2D%get_int(1,'thumbnx')
        ntilesy   = proj%os_cls2D%get_int(1,'thumbny')
        do iori = 1,noris
            if( .not. proj%os_cls2D%isthere(iori,'thumbidx') ) cycle
            thumbidx = proj%os_cls2D%get_int(iori,'thumbidx')
            if( thumbidx < 1 .or. thumbidx > n_selected ) cycle
            cls_idx(thumbidx) = iori
            cls_pop(thumbidx) = nint(proj%os_cls2D%get(iori,'pop'))
            cls_res(thumbidx) = proj%os_cls2D%get(iori,'res')
        enddo
    end subroutine snapshot_cavgs_meta

    ! The process's resident memory, logged at @p phase.
    subroutine log_rss( phase )
        use, intrinsic :: iso_c_binding, only: c_int64_t, c_double
        character(len=*), intent(in) :: phase
        integer(c_int64_t) :: current_rss, peak_rss
        real(c_double)     :: current_mib, peak_mib
        current_rss = get_current_rss_bytes()
        peak_rss    = get_peak_rss_bytes()
        if( current_rss < 0_c_int64_t .or. peak_rss < 0_c_int64_t ) return
        current_mib = real(current_rss, c_double) / 1048576.0_c_double
        peak_mib    = real(peak_rss,    c_double) / 1048576.0_c_double
        write(logfhandle,'(A,A,A,F10.1,A,F10.1,A)') '>>> RSS ', trim(phase), ': current=', current_mib, ' MiB peak=', peak_mib, ' MiB'
        call flush(logfhandle)
    end subroutine log_rss

    !> The pool's classified state as a project for 3D: the stacks whose particles have been
    !! through an iteration (updatecnt > 0), with their micrographs and particles (2D and 3D, the
    !! 3D a copy of the 2D), and the classes. Stacks are renumbered from 1 and their particle
    !! ranges and stack indices follow, so the project is self-consistent; a particle keeps its
    !! image index in its stack. Particles never classified are left out. @p nstks is the number
    !! of stacks published (0: nothing to publish, @p pub is empty).
    subroutine build_pool_publication( src, pub, nstks )
        class(sp_project), intent(inout) :: src
        class(sp_project), intent(inout) :: pub
        integer,           intent(out)   :: nstks
        logical, allocatable :: l_stk(:)
        integer :: nstks_src, istk, jstk, iptcl, jptcl, fromp, top, nptcls, stkind_src, ind_in_stk
        call pub%kill
        nstks     = 0
        nstks_src = src%os_stk%get_noris()
        if( nstks_src == 0 ) return
        if( src%os_mic%get_noris() /= nstks_src ) THROW_HARD('build_pool_publication: # micrographs /= # stacks')
        allocate(l_stk(nstks_src), source=.false.)
        nptcls = 0
        do istk = 1,nstks_src
            fromp = src%os_stk%get_fromp(istk)
            top   = src%os_stk%get_top(istk)
            do iptcl = fromp,top
                if( src%os_ptcl2D%get_updatecnt(iptcl) > 0 )then
                    l_stk(istk) = .true.
                    exit
                endif
            enddo
            if( l_stk(istk) ) nptcls = nptcls + top - fromp + 1
        enddo
        nstks = count(l_stk)
        if( nstks == 0 ) return
        pub%projinfo = src%projinfo
        pub%compenv  = src%compenv
        pub%jobproc  = src%jobproc
        pub%os_optics = src%os_optics
        pub%os_cls2D  = src%os_cls2D
        call pub%os_mic%new(nstks,    is_ptcl=.false.)
        call pub%os_stk%new(nstks,    is_ptcl=.false.)
        call pub%os_ptcl2D%new(nptcls, is_ptcl=.true.)
        jstk  = 0
        jptcl = 0
        do istk = 1,nstks_src
            if( .not. l_stk(istk) ) cycle
            jstk  = jstk + 1
            fromp = src%os_stk%get_fromp(istk)
            top   = src%os_stk%get_top(istk)
            call pub%os_mic%transfer_ori(jstk, src%os_mic, istk)
            call pub%os_stk%transfer_ori(jstk, src%os_stk, istk)
            call pub%os_stk%set(jstk, 'fromp', jptcl + 1)
            call pub%os_stk%set(jstk, 'top',   jptcl + top - fromp + 1)
            do iptcl = fromp,top
                jptcl = jptcl + 1
                ! the image index in the stack, as the source maps it, before the ranges change
                call src%map_ptcl_ind2stk_ind('ptcl2D', iptcl, stkind_src, ind_in_stk)
                call pub%os_ptcl2D%transfer_ori(jptcl, src%os_ptcl2D, iptcl)
                call pub%os_ptcl2D%set_stkind(jptcl, jstk)
                call pub%os_ptcl2D%set(jptcl, 'indstk', ind_in_stk)
            enddo
        enddo
        pub%os_ptcl3D = pub%os_ptcl2D
        deallocate(l_stk)
    end subroutine build_pool_publication

    !> Publishes the pool's classified state for 3D as @p projfile (build_pool_publication), with
    !! the pool's class averages (and their even and odd halves) and FRCs at the native sampling
    !! beside it, registered with the pool's mask diameter. The class averages and FRCs are written
    !! first and the project last, under a temporary name renamed into place, so a reader that
    !! finds the project finds it complete. @p nstks is the number of stacks published (0: nothing
    !! was written).
    subroutine publish_pool_state( params, projfile, nstks )
        class(parameters), intent(in)  :: params
        class(string),     intent(in)  :: projfile
        integer,           intent(out) :: nstks
        type(sp_project) :: pub
        type(class_frcs) :: frcs
        type(string)     :: pool_refs, cavgsfname, frcsfname
        call build_pool_publication(pool_proj, pub, nstks)
        if( nstks == 0 ) return
        call pool_publication_names(projfile, cavgsfname, frcsfname)
        pool_refs = string(POOL_DIR)//refs_glob
        if( l_scaling )then
            call rescale_cavgs(params, pool_refs, cavgsfname)
            call rescale_cavgs(params, add2fbody(pool_refs, MRC_EXT, '_even'), add2fbody(cavgsfname, MRC_EXT, '_even'))
            call rescale_cavgs(params, add2fbody(pool_refs, MRC_EXT, '_odd'),  add2fbody(cavgsfname, MRC_EXT, '_odd'))
            call frcs%read(string(POOL_DIR)//FRCS_FILE)
            call frcs%pad(params%smpd, params%box)
            call frcs%write(frcsfname)
            call frcs%kill
        else
            call simple_copy_file(pool_refs, cavgsfname)
            call simple_copy_file(add2fbody(pool_refs, MRC_EXT, '_even'), add2fbody(cavgsfname, MRC_EXT, '_even'))
            call simple_copy_file(add2fbody(pool_refs, MRC_EXT, '_odd'),  add2fbody(cavgsfname, MRC_EXT, '_odd'))
            call simple_copy_file(string(POOL_DIR)//FRCS_FILE, frcsfname)
        endif
        call pub%os_out%kill
        call pub%add_cavgs2os_out(cavgsfname, params%smpd, 'cavg', mskdiam=params%mskdiam)
        call pub%add_frcs2os_out(frcsfname, 'frc2D')
        call pub%write(projfile, tempfile=.true.)
        write(logfhandle,'(A,A,A,I8,A,I8,A)') '>>> PUBLISHED THE POOL FOR 3D: ', projfile%to_char(), ', ',&
            &nstks, ' STACK(S), ', pub%os_ptcl2D%count_state_gt_zero(), ' PARTICLE(S)'
        call pub%kill
    end subroutine publish_pool_state

    !> Removes the publication @p projfile and the class averages and FRCs written beside it.
    subroutine delete_pool_publication( projfile )
        class(string), intent(in) :: projfile
        type(string) :: cavgsfname, frcsfname
        call pool_publication_names(projfile, cavgsfname, frcsfname)
        call del_file(projfile)
        call del_file(cavgsfname)
        call del_file(add2fbody(cavgsfname, MRC_EXT, '_even'))
        call del_file(add2fbody(cavgsfname, MRC_EXT, '_odd'))
        call del_file(frcsfname)
    end subroutine delete_pool_publication

    ! the class averages and FRCs files of the publication @p projfile, beside it
    subroutine pool_publication_names( projfile, cavgsfname, frcsfname )
        class(string), intent(in)  :: projfile
        type(string),  intent(out) :: cavgsfname, frcsfname
        type(string) :: fbody
        fbody      = get_fbody(projfile, METADATA_EXT, separator=.false.)
        cavgsfname = fbody//'_cavgs'//STK_EXT
        frcsfname  = fbody//'_'//FRCS_FILE
    end subroutine pool_publication_names

    ! Handles user inputted class rejection
    subroutine write_repick_refs(refsout)
        class(string), intent(in) :: refsout
        type(image)  :: img        
        type(string) :: refsin
        integer      :: icls, i
        if( .not. l_stream2D_active ) return
        if( pool_proj%os_cls2D%get_noris() == 0 ) return
        if(repick_iteration .lt. 1) return
        refsin = CAVGS_ITER_FBODY//int2str_pad(repick_iteration,3)//MRC_EXT
        if(.not. file_exists(string(POOL_DIR)//refsin)) return
        if(file_exists(refsout) ) call del_file(refsout)
        call img%new([pool_dims%box,pool_dims%box,1], pool_dims%smpd)
        img = 0.
        if(allocated(repick_selection)) then
            do i = 1, size(repick_selection)
                icls = repick_selection(i)
                if( icls <= 0 ) cycle
                if( icls > ncls_glob ) cycle
                call img%read(string(POOL_DIR)//refsin,icls)
                call img%write(refsout,i)
            end do
        end if    
        call update_stack_nimgs(refsout, size(repick_selection))
        call img%kill
    end subroutine write_repick_refs

end module simple_stream_refine2D_utils
