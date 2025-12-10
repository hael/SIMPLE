module simple_stream_cluster2D_utils
include 'simple_lib.f08'
use simple_class_frcs,         only: class_frcs
use simple_cmdline,            only: cmdline
use simple_euclid_sigma2,      only: sigma2_star_from_iter
use simple_image,              only: image
use simple_parameters,         only: params_glob
use simple_sp_project,         only: sp_project
use simple_stack_io,           only: stack_io
use simple_starproject,        only: starproject
use simple_starproject_stream, only: starproject_stream
use simple_commanders_cluster2D
use simple_rec_list
use simple_stream2D_state
implicit none

public :: apply_snapshot_selection
public :: cleanup_root_folder
public :: consolidate_sigmas
public :: debug_print
public :: get_box
public :: get_boxa
public :: rank_cavgs
public :: rescale_cavgs
public :: set_dimensions
public :: set_resolution_limits
public :: setup_downscaling
public :: terminate_chunks
public :: terminate_stream2D
public :: test_repick
public :: tidy_2Dstream_iter
public :: transfer_cavg
public :: update_user_params2D
public :: write_project_stream2D
public :: write_repick_refs
private
#include "simple_local_flags.inc"

type(starproject_stream) :: starproj_stream
logical, parameter       :: DEBUG_HERE = .false.
integer, allocatable     :: repick_selection(:)  ! selection for selecting classes for re-picking
integer(timer_int_kind)  :: t
integer                  :: snapshot_jobid   = 0 ! nice job id for snapshot
integer                  :: repick_iteration = 0 ! iteration to select classes from for re-picking
character(16)            :: snapshot_last = ""   ! timestamp string for last snapshot

contains

    subroutine apply_snapshot_selection(snapshot_projfile)
        type(sp_project),     intent(inout) :: snapshot_projfile
        logical, allocatable :: cls_mask(:)
        integer              :: nptcls_rejected, ncls_rejected, iptcl
        integer              :: icls, jcls, i
        if( snapshot_projfile%os_cls2D%get_noris() == 0 ) return
        if(allocated(snapshot_selection)) then
            allocate(cls_mask(ncls_glob), source=.false.)
            do i = 1, size(snapshot_selection)
                icls = snapshot_selection(i)
                if( icls == 0 ) cycle
                if( icls > ncls_glob ) cycle
                cls_mask(icls) = .true.
            enddo
            deallocate(snapshot_selection)
        else
            allocate(cls_mask(ncls_glob), source=.true.)
        endif
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
 !       call qsys_cleanup(nparts=params_glob%nparts_pool) ! cyril - fails on stream 2d classification restart
        call simple_rmdir(SIGMAS_DIR)
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
                call del_file(CLUSTER2D_FINISHED)
                call del_file('simple_script_single')
            endif
        endif
    end subroutine cleanup_root_folder

    ! Private utility to aggregate sigma2
    subroutine consolidate_sigmas( project, nstks )
        use simple_euclid_sigma2, only: consolidate_sigma2_groups, average_sigma2_groups
        type(sp_project),           intent(in) :: project
        integer,                    intent(in) :: nstks
        type(string) :: stack_fname, ext, fbody
        type(string), allocatable :: sigma_fnames(:)
        integer :: i, istk
        if( l_update_sigmas )then
            if( trim(params_glob%sigma_est).eq.'group' )then
                allocate(sigma_fnames(nstks))
                do istk = 1,nstks
                    call project%os_stk%getter(istk,'stk',stack_fname)
                    stack_fname = basename(stack_fname)
                    ext         = fname2ext(stack_fname)
                    fbody       = get_fbody(stack_fname, ext)
                    sigma_fnames(istk) = SIGMAS_DIR//'/'//fbody%to_char()//STAR_EXT
                enddo
                call consolidate_sigma2_groups(sigma2_star_from_iter(pool_iter), sigma_fnames)
                deallocate(sigma_fnames)
            else
                ! sigma_est=global & first iteration
                if( pool_iter==1 )then
                    allocate(sigma_fnames(glob_chunk_id))
                    do i = 1,glob_chunk_id
                        sigma_fnames(i) = SIGMAS_DIR//'/chunk_'//int2str(i)//STAR_EXT
                    enddo
                    call average_sigma2_groups(sigma2_star_from_iter(pool_iter), sigma_fnames)
                    deallocate(sigma_fnames)
                endif
            endif
            do i = 1,params_glob%nparts_pool
                call del_file(SIGMA2_FBODY//int2str_pad(i,numlen)//'.dat')
            enddo
        endif
    end subroutine consolidate_sigmas

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
        refs_ranked = add2fbody(refs_glob, params_glob%ext ,'_ranked')
        call cline_rank_cavgs%set('projfile', orig_projfile)
        if( l_wfilt )then
            stk = string(POOL_DIR)//add2fbody(refs_glob,params_glob%ext,WFILT_SUFFIX)
        else
            stk = string(POOL_DIR)//refs_glob
        endif
        call cline_rank_cavgs%set('stk',            stk)
        call cline_rank_cavgs%set('outstk', refs_ranked)
        call xrank_cavgs%execute_safe(cline_rank_cavgs)
        call cline_rank_cavgs%kill
    end subroutine rank_cavgs

    ! Final rescaling of references
    subroutine rescale_cavgs( src, dest )
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
        call img_pad%new([params_glob%box,params_glob%box,1],params_glob%smpd)
        cls_pop = nint(pool_proj%os_cls2D%get_all('pop'))
        call find_ldim_nptcls(src,ldim,ncls_here)
        call stkio_r%open(src, pool_dims%smpd, 'read', bufsz=ncls_here)
        call stkio_r%read_whole
        call stkio_w%open(dest_here, params_glob%smpd, 'write', box=params_glob%box, bufsz=ncls_here)
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

    subroutine set_dimensions
        call setup_downscaling
        pool_dims%smpd  = params_glob%smpd_crop
        pool_dims%box   = params_glob%box_crop
        pool_dims%boxpd = 2 * round2even(params_glob%alpha * real(params_glob%box_crop/2)) ! logics from parameters
        pool_dims%msk   = params_glob%msk_crop
        chunk_dims = pool_dims ! chunk & pool have the same dimensions to start with
        ! Scaling-related command lines update
        call cline_cluster2D_chunk%set('smpd_crop', chunk_dims%smpd)
        call cline_cluster2D_chunk%set('box_crop',  chunk_dims%box)
        call cline_cluster2D_chunk%set('msk_crop',  chunk_dims%msk)
        call cline_cluster2D_chunk%set('box',       params_glob%box)
        call cline_cluster2D_chunk%set('smpd',      params_glob%smpd)
        call cline_cluster2D_pool%set('smpd_crop',  pool_dims%smpd)
        call cline_cluster2D_pool%set('box_crop',   pool_dims%box)
        call cline_cluster2D_pool%set('msk_crop',   pool_dims%msk)
        call cline_cluster2D_pool%set('box',        params_glob%box)
        call cline_cluster2D_pool%set('smpd',       params_glob%smpd)
    end subroutine set_dimensions

    ! private routine for resolution-related updates to command-lines
    subroutine set_resolution_limits
        lpstart = max(lpstart, 2.0*params_glob%smpd_crop)
        if( l_no_chunks )then
            params_glob%lpstop = lpstop
        else
            if( master_cline%defined('lpstop') )then
                params_glob%lpstop = max(2.0*params_glob%smpd_crop,params_glob%lpstop)
            else
                params_glob%lpstop = 2.0*params_glob%smpd_crop
            endif
            call cline_cluster2D_chunk%delete('lp')
            call cline_cluster2D_chunk%set('lpstart', lpstart)
            call cline_cluster2D_chunk%set('lpstop',  lpstart)
        endif
        call cline_cluster2D_pool%set('lpstart',  lpstart)
        call cline_cluster2D_pool%set('lpstop',   params_glob%lpstop)
        if( .not.master_cline%defined('cenlp') )then
            call cline_cluster2D_chunk%set('cenlp', lpcen)
            call cline_cluster2D_pool%set( 'cenlp', lpcen)
        else
            call cline_cluster2D_chunk%set('cenlp', params_glob%cenlp)
            call cline_cluster2D_pool%set( 'cenlp', params_glob%cenlp)
        endif
        ! Will use resolution update scheme from abinitio2D
        if( l_abinitio2D .and. (.not.l_no_chunks))then
            if( master_cline%defined('lpstop') )then
                ! already set above
            else
                call cline_cluster2D_chunk%delete('lpstop')
            endif
        endif
        write(logfhandle,'(A,F5.1)') '>>> POOL STARTING LOW-PASS LIMIT (IN A): ', lpstart
        write(logfhandle,'(A,F5.1)') '>>> POOL   HARD RESOLUTION LIMIT (IN A): ', params_glob%lpstop
        write(logfhandle,'(A,F5.1)') '>>> CENTERING     LOW-PASS LIMIT (IN A): ', lpcen
    end subroutine set_resolution_limits

    ! Determines dimensions for downscaling
    subroutine setup_downscaling
        real    :: SMPD_TARGET = MAX_SMPD  ! target sampling distance
        real    :: smpd, scale_factor
        integer :: box
        if( params_glob%box == 0 ) THROW_HARD('FATAL ERROR')
        scale_factor          = 1.0
        params_glob%smpd_crop = params_glob%smpd
        params_glob%box_crop  = params_glob%box
        if( l_scaling .and. params_glob%box >= CHUNK_MINBOXSZ )then
            call autoscale(params_glob%box, params_glob%smpd, SMPD_TARGET, box, smpd, scale_factor, minbox=CHUNK_MINBOXSZ)
            l_scaling = box < params_glob%box
            if( l_scaling )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params_glob%box,'/',box
                params_glob%smpd_crop = smpd
                params_glob%box_crop  = box
            endif
        endif
        params_glob%msk_crop = round2even(params_glob%mskdiam / params_glob%smpd_crop / 2.)
    end subroutine setup_downscaling

    ! ends chunks processing
    subroutine terminate_chunks
        integer :: ichunk
        do ichunk = 1,params_glob%nchunks
            call chunks(ichunk)%terminate_chunk
        enddo
    end subroutine terminate_chunks

    ! ends processing, generates project & cleanup
    subroutine terminate_stream2D( project_list, optics_dir)
        class(rec_list), optional, intent(inout) :: project_list
        class(string),   optional, intent(in)    :: optics_dir
        type(string) :: mapfileprefix
        integer      :: ipart, lastmap
        call terminate_chunks
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
                    refs_glob = CAVGS_ITER_FBODY//int2str_pad(pool_iter,3)//params_glob%ext%to_char()
                    ! tricking the asynchronous master process to come to a hard stop
                    call simple_touch(POOL_DIR//TERM_STREAM)
                    do ipart = 1,params_glob%nparts_pool
                        call simple_touch(POOL_DIR//JOB_FINISHED_FBODY//int2str_pad(ipart,numlen))
                    enddo
                    call simple_touch(POOL_DIR//'CAVGASSEMBLE_FINISHED')
                endif
            endif
            if( pool_iter >= 1 )then
                call write_project_stream2D(write_star=.true., clspath=.true., optics_dir=optics_dir)
                call rank_cavgs
            endif
        endif
        ! cleanup
        call simple_rmdir(SIGMAS_DIR)
        call del_file(POOL_DIR//POOL_PROJFILE)
        call del_file(projfile4gui)
        if( .not. DEBUG_HERE )then
            call qsys_cleanup
        endif

        contains

            ! no pool clustering performed, all available info is written down
            subroutine write_raw_project
                if( present(project_list) )then
                    if( project_list%size() > 0 )then
                        lastmap = 0
                        if(present(optics_dir)) then
                            call get_latest_optics_map()
                            if(lastmap .gt. 0) then
                                mapfileprefix = optics_dir//'/'//OPTICS_MAP_PREFIX//int2str(lastmap)
                                call pool_proj%import_optics_map(mapfileprefix)
                            endif
                        endif
                        call pool_proj%projrecords2proj(project_list)
                        call starproj_stream%copy_micrographs_optics(pool_proj, verbose=DEBUG_HERE)
                        call starproj_stream%stream_export_micrographs(pool_proj, params_glob%outdir, optics_set=.true.)
                        call starproj_stream%stream_export_particles_2D(pool_proj, params_glob%outdir, optics_set=.true.)
                        call pool_proj%write(orig_projfile)
                    endif
                endif
            end subroutine write_raw_project

            subroutine get_latest_optics_map()
                type(string), allocatable :: map_list(:)
                type(string) :: map_i_str, map_str
                integer      :: imap, prefix_len, testmap
                lastmap = 0
                if(optics_dir .ne. "") then
                    if(dir_exists(optics_dir)) call simple_list_files(optics_dir%to_char()//'/'//OPTICS_MAP_PREFIX//'*'// TXT_EXT, map_list)
                endif
                if(allocated(map_list)) then
                    prefix_len = len(optics_dir%to_char() // '/' // OPTICS_MAP_PREFIX) + 1
                    do imap=1, size(map_list)
                        map_str   = map_list(imap)%to_char([prefix_len,map_list(imap)%strlen_trim()])
                        map_i_str = swap_suffix(map_str, "", TXT_EXT)
                        testmap   = map_i_str%to_int()
                        if(testmap > lastmap) lastmap = testmap
                        call map_str%kill
                        call map_i_str%kill
                    enddo
                    deallocate(map_list)
                endif
            end subroutine get_latest_optics_map

    end subroutine terminate_stream2D

    logical function test_repick()
        test_repick = allocated(repick_selection)
    end function test_repick

    ! Removes some unnecessary files
    subroutine tidy_2Dstream_iter
        type(string) :: prefix
        if( pool_iter > 5 )then
            prefix = POOL_DIR//CAVGS_ITER_FBODY//int2str_pad(pool_iter-5,3)
            call del_file(prefix//JPG_EXT)
            call del_file(prefix//'_even'//params_glob%ext%to_char())
            call del_file(prefix//'_odd'//params_glob%ext%to_char())
            call del_file(prefix//params_glob%ext%to_char())
            if( l_wfilt )then
                call del_file(prefix//WFILT_SUFFIX//'_even'//params_glob%ext%to_char())
                call del_file(prefix//WFILT_SUFFIX//'_odd'//params_glob%ext%to_char())
                call del_file(prefix//WFILT_SUFFIX//params_glob%ext%to_char())
            endif
            call del_file(POOL_DIR//CLS2D_STARFBODY//'_iter'//int2str_pad(pool_iter-5,3)//STAR_EXT)
            call del_file(prefix // '.jpg')
            if( l_update_sigmas ) call del_file(string(POOL_DIR)//sigma2_star_from_iter(pool_iter-5))
        endif
    end subroutine tidy_2Dstream_iter

    !>  Transfers references & partial arrays from a chunk to the pool
    subroutine transfer_cavg( refs_in, dir, indin, refs_out, indout )
        class(string), intent(in) :: refs_in, dir, refs_out
        integer,       intent(in) :: indin, indout
        type(image)  :: img, img2
        type(string) :: stkout, stkin, refs_out_here, refs_in_here
        integer      :: ipart
        call debug_print('in transfer_cavg '//int2str(indin)//' '//int2str(indout))
        refs_in_here = refs_in
        call img%new([chunk_dims%box,chunk_dims%box,1], chunk_dims%smpd)
        call img2%new([pool_dims%box,pool_dims%box,1], pool_dims%smpd)
        ! making sure we are writing to the correct folder
        refs_out_here = string(POOL_DIR)//refs_out
        ! merged class
        call read_pad_write(refs_in_here, indin, refs_out_here, indout)
        if( l_wfilt )then
            stkout = add2fbody(refs_out_here,params_glob%ext,WFILT_SUFFIX)
            if( pool_dims%box > chunk_dims%box )then
                call img2%write(stkout,indout)
            else
                call img%write(stkout,indout)
            endif
        endif
        ! e/o
        stkin  = add2fbody(refs_in_here, params_glob%ext,'_even')
        stkout = add2fbody(refs_out_here,params_glob%ext,'_even')
        call read_pad_write(stkin, indin, stkout, indout)
        stkin  = add2fbody(refs_in_here,params_glob%ext,'_odd')
        stkout = add2fbody(refs_out_here,params_glob%ext,'_odd')
        call read_pad_write(stkin, indin, stkout, indout)
        ! temporary matrices, logics from chunk%read
        call img%new([chunk_dims%boxpd,chunk_dims%boxpd,1], chunk_dims%smpd)
        call img%zero_and_flag_ft
        if( pool_dims%box > chunk_dims%box ) call img2%new([pool_dims%boxpd,pool_dims%boxpd,1], pool_dims%smpd)
        call write_inside_ftstack('/cavgs_even_part',     'cavgs_even_part',     'cavgs_even_wfilt_part')
        call write_inside_ftstack('/cavgs_odd_part',      'cavgs_odd_part',      'cavgs_odd_wfilt_part')
        call write_inside_ftstack('/ctfsqsums_even_part', 'ctfsqsums_even_part', 'ctfsqsums_even_wfilt_part')
        call write_inside_ftstack('/ctfsqsums_odd_part',  'ctfsqsums_odd_part',  'ctfsqsums_odd_wfilt_part')
        ! cleanup
        call img%kill
        call img2%kill
        contains

            subroutine read_pad_write(strin, iin, strout, iout)
                class(string), intent(in) :: strin, strout
                integer,       intent(in) :: iin, iout
                call img%read(strin, iin)
                if( pool_dims%box > chunk_dims%box )then
                    call img2%zero_and_flag_ft
                    call img%fft
                    call img%pad(img2, antialiasing=.false.)
                    call img2%ifft
                    call img2%write(strout,iout)
                    call img%zero_and_unflag_ft
                else
                    call img%write(strout,iout)
                endif
            end subroutine read_pad_write

            subroutine write_inside_ftstack(tmplin, tmplout, tmplout_wfilt)
                character(len=*), intent(in) :: tmplin, tmplout, tmplout_wfilt
                stkin = dir//trim(tmplin)//params_glob%ext%to_char()
                call img%read(stkin, indin)
                if( pool_dims%box > chunk_dims%box )then
                    call img%pad(img2, antialiasing=.false.)
                    do ipart = 1,params_glob%nparts_pool
                        stkout = POOL_DIR//trim(tmplout)//int2str_pad(ipart,numlen)//params_glob%ext%to_char()
                        call img2%write(stkout,indout)
                        if( l_wfilt )then
                            stkout = POOL_DIR//trim(tmplout_wfilt)//int2str_pad(ipart,numlen)//params_glob%ext%to_char()
                            call img2%write(stkout,indout)
                        endif
                    enddo
                else
                    do ipart = 1,params_glob%nparts_pool
                        stkout = POOL_DIR//trim(tmplout)//int2str_pad(ipart,numlen)//params_glob%ext%to_char()
                        call img%write(stkout,indout)
                        if( l_wfilt )then
                            stkout = POOL_DIR//trim(tmplout_wfilt)//int2str_pad(ipart,numlen)//params_glob%ext%to_char()
                            call img%write(stkout,indout)
                        endif
                    enddo
                endif
            end subroutine write_inside_ftstack

    end subroutine transfer_cavg

    !> Updates current parameters with user input
    subroutine update_user_params2D( cline_here, updated, update_arguments)
        type(cmdline), intent(inout) :: cline_here
        logical,       intent(out)   :: updated
        type(json_value), pointer, optional, intent(inout) :: update_arguments
        character(kind=CK,len=:), allocatable :: snapshot
        integer,                  allocatable :: prune_selection(:) ! selection for pruning classes on the fly, JSON-carried selection object / string
        type(oris)            :: os
        type(json_core)       :: json
        type(string)          :: val
        real                  :: lpthres, ndev
        integer               :: lpthres_int, mskdiam_int
        logical               :: found
        updated = .false.
        call os%new(1, is_ptcl=.false.)
        if( file_exists(USER_PARAMS2D) )then
            call os%read(string(USER_PARAMS2D))
            ! class resolution threshold for rejection
            if( os%isthere(1,'lpthres') )then
                lpthres = os%get(1,'lpthres')
                if( abs(lpthres-params_glob%lpthres) > 0.001 )then
                    if( lpthres < 3.0*chunk_dims%smpd )then
                        write(logfhandle,'(A,F8.2)')'>>> REJECTION lpthres TOO LOW: ',lpthres
                    else
                        params_glob%lpthres = lpthres
                        write(logfhandle,'(A,F8.2)')'>>> REJECTION lpthres UPDATED TO: ',params_glob%lpthres
                        updated = .true.
                    endif
                endif
            endif
            ! class standard deviation of resolution threshold for rejection
            if( os%isthere(1,'ndev2D') )then
                ndev = os%get(1,'ndev2D')
                if( abs(ndev-params_glob%ndev) > 0.001 )then
                    if( ndev < 0.1 )then
                        write(logfhandle,'(A,F8.2)')'>>> REJECTION NDEV TOO LOW: ',ndev
                    else
                        params_glob%ndev = ndev
                        write(logfhandle,'(A,F8.2)')'>>> REJECTION NDEV   UPDATED TO: ',params_glob%ndev
                        updated = .true.
                    endif
                endif
            endif
            ! class rejection
            if( os%isthere(1,'reject_cls') )then
                val = os%get_str(1, 'reject_cls')
                if( val .ne. trim(params_glob%reject_cls) )then
                    select case(val%to_char())
                        case('yes')
                            write(logfhandle,'(A)')'>>> ACTIVATING CLASS REJECTION'
                            params_glob%reject_cls = val%to_char()
                            updated = .true.
                        case('no')
                            write(logfhandle,'(A)')'>>> DE-ACTIVATING CLASS REJECTION'
                            params_glob%reject_cls = val%to_char()
                            updated = .true.
                        case('old')
                            write(logfhandle,'(A)')'>>> DE-ACTIVATING IMAGE MOMENTS-BASED CLASS REJECTION'
                            params_glob%reject_cls = val%to_char()
                            updated = .true.
                        case DEFAULT
                            THROW_WARN('Unknown flag for class rejection: '//val%to_char())
                    end select
                endif
            endif
            ! remove once processed
            call del_file(USER_PARAMS2D)
        endif
        ! nice
        if(present(update_arguments)) then
            if(associated(update_arguments)) then
                ! snapshot
                if(snapshot_jobid .eq. 0) then
                    call json%get(update_arguments, 'snapshot', snapshot, found)
                    if(found) then
                        if(allocated(snapshot_selection)) deallocate(snapshot_selection)
                        allocate(snapshot_selection(0))
                        params_glob%snapshot = snapshot
                        params_glob%updated  = 'yes'
                        write(logfhandle,'(A,A)')'>>> SNAPSHOT REQUESTED: ', snapshot
                    end if
                    ! snapshot selection
                    call json%get(update_arguments, 'snapshot_selection', snapshot_selection, found)
                    if(found) then
                        write(logfhandle,'(A,A)')'>>> SNAPSHOT SELECTION RECEIVED'
                    end if
                    ! snapshot iteration
                    call json%get(update_arguments, 'snapshot_iteration', snapshot_iteration, found)
                    if(found) then
                        write(logfhandle,'(A,A)')'>>> SNAPSHOT ITERATION RECEIVED'
                    end if
                    ! snapshot job id
                    call json%get(update_arguments, 'snapshot_jobid', snapshot_jobid, found)
                    if(found) then
                        write(logfhandle,'(A,A)')'>>> SNAPSHOT JOB ID RECEIVED'
                    end if
                end if
                ! prune selection
                call json%get(update_arguments, 'prune_selection', prune_selection, found)
                if(found) then
                    write(logfhandle,'(A,A)')'>>> PRUNE SELECTION RECEIVED'
                end if
                ! repick selection
                call json%get(update_arguments, 'repick_selection', repick_selection, found)
                if(found) then
                    write(logfhandle,'(A,A)')'>>> REPICK SELECTION RECEIVED'
                end if
                ! repick iteration
                call json%get(update_arguments, 'repick_iteration', repick_iteration, found)
                if(found) then
                    write(logfhandle,'(A,A)')'>>> REPICK ITERATION RECEIVED'
                end if
                ! lpthres
                call json%get(update_arguments, 'lpthres', lpthres_int, found)
                if(found) then
                    write(logfhandle,'(A,A)')'>>> LPTHRES UPDATE RECEIVED'
                    if( real(lpthres_int) < 1.0)then
                        lpthres_type = "auto"
                        call mskdiam2streamresthreshold(params_glob%mskdiam, params_glob%lpthres)
                        write(logfhandle,'(A,F8.2)')'>>> REJECTION lpthres SET TO AUTO: ', params_glob%lpthres
                        updated = .true.
                    else if( real(lpthres_int) .gt. LOWRES_REJECT_THRESHOLD)then
                        lpthres_type = "off"
                        params_glob%lpthres = real(lpthres_int)
                        write(logfhandle,'(A,F8.2)')'>>> REJECTION lpthres SET TO OFF: ', params_glob%lpthres
                        updated = .true.
                    else
                        lpthres_type = "manual"
                        params_glob%lpthres = real(lpthres_int)
                        write(logfhandle,'(A,F8.2)')'>>> REJECTION lpthres UPDATED TO: ',params_glob%lpthres
                        updated = .true.
                    endif 
                end if
                ! mskdiam
                call json%get(update_arguments, 'mskdiam', mskdiam_int, found)
                if(found) then
                    write(logfhandle,'(A,A)')'>>> MASK DIAMETER UPDATE RECEIVED'
                    if( real(mskdiam_int) .gt. 49.0)then
                        params_glob%mskdiam = real(mskdiam_int)
                        call cline_cluster2D_chunk%set('mskdiam',   params_glob%mskdiam)
                        call cline_cluster2D_pool%set('mskdiam',   params_glob%mskdiam)
                        write(logfhandle,'(A,F8.2)')'>>> MASK DIAMETER UPDATED TO: ', params_glob%mskdiam
                        updated = .true.
                    endif 
                end if
                call json%destroy(update_arguments)
            end if
        end if
        call os%kill
    end subroutine update_user_params2D

    !> produces consolidated project
    subroutine write_project_stream2D( write_star, clspath, snapshot_projfile, snapshot_starfile_base, optics_dir)
        logical,           optional, intent(in)    :: write_star
        logical,           optional, intent(in)    :: clspath
        class(string),     optional, intent(in)    :: snapshot_projfile, snapshot_starfile_base, optics_dir
        integer                        :: snapshot_complete_jobid = 0 ! nice job id for completed snapshot
        type(class_frcs)               :: frcs
        type(oris)                     :: os_backup
        type(sp_project)               :: snapshot_proj
      !  type(simple_nice_communicator) :: snapshot_comm
        type(json_core)                :: json
        type(string)                   :: projfile, projfname, cavgsfname, frcsfname, src, dest, mapfileprefix
        type(string)                   :: pool_refs, l_stkname, l_frcsname
        logical                        :: l_write_star, l_clspath, l_snapshot, snapshot_proj_found = .false.
        logical,     parameter         :: DEBUG_HERE = .false.
        real                           :: l_smpd
        integer                        :: l_ncls, lastmap
        l_write_star = .false.
        l_clspath    = .false.
        l_snapshot   = .false.
        if(present(write_star)) l_write_star = write_star
        if(present(clspath))    l_clspath    = clspath
        ! file naming
        projfname  = get_fbody(orig_projfile, METADATA_EXT, separator=.false.)
        cavgsfname = get_fbody(refs_glob, params_glob%ext, separator=.false.)
        frcsfname  = get_fbody(FRCS_FILE, BIN_EXT, separator=.false.)
        call pool_proj%projinfo%set(1,'projname', projfname)
        projfile   = projfname//METADATA_EXT
        call pool_proj%projinfo%set(1,'projfile', projfile)
        cavgsfname = cavgsfname//params_glob%ext
        frcsfname  = frcsfname//BIN_EXT
        pool_refs  = string(POOL_DIR)//refs_glob
        if(present(snapshot_projfile)) l_snapshot = snapshot_iteration > 0
        lastmap = 0
        if( l_snapshot ) then
            write(logfhandle, '(A,I4,A,A,A,I0,A)') ">>> WRITING SNAPSHOT FROM ITERATION ", snapshot_iteration, snapshot_projfile%to_char(),&
            &' AT: ',cast_time_char(simple_gettime()), snapshot_jobid, params_glob%niceserver%to_char()
            call json%destroy(snapshot_json)
            nullify(snapshot_json)
            if(snapshot_iteration .eq. pool_iter) then
                call snapshot_proj%copy(pool_proj)
                call snapshot_proj%add_frcs2os_out(string(POOL_DIR)//FRCS_FILE, 'frc2D')
                snapshot_proj_found = .true.
            else
                call snapshot_proj%copy(pool_proj_history(snapshot_iteration))
                snapshot_proj_found = .true.
            end if
            if(snapshot_proj_found) then
                if(.not. file_exists(stemname(stemname(snapshot_projfile)))) call simple_mkdir(stemname(stemname(snapshot_projfile)))
                if(.not. file_exists(stemname(snapshot_projfile)))           call simple_mkdir(stemname(snapshot_projfile))
                if(present(optics_dir)) then
                    call get_latest_optics_map()
                    if(lastmap .gt. 0) then
                        mapfileprefix = optics_dir // '/' // OPTICS_MAP_PREFIX // int2str(lastmap)
                        call snapshot_proj%import_optics_map(mapfileprefix)
                    endif
                endif
                call apply_snapshot_selection(snapshot_proj)
                call snapshot_proj%get_cavgs_stk(l_stkname, l_ncls, l_smpd)
                call snapshot_proj%get_frcs(l_frcsname, 'frc2D')
                call simple_copy_file(l_stkname, stemname(snapshot_projfile) // "/cavgs" // STK_EXT)
                call simple_copy_file(add2fbody(l_stkname, params_glob%ext,'_even'), stemname(snapshot_projfile) // "/cavgs_even" // STK_EXT)
                call simple_copy_file(add2fbody(l_stkname, params_glob%ext,'_odd'),  stemname(snapshot_projfile) // "/cavgs_odd"  // STK_EXT)
                call simple_copy_file(l_frcsname, stemname(snapshot_projfile) // '/' // FRCS_FILE)
                call snapshot_proj%os_out%kill
                call snapshot_proj%add_cavgs2os_out(stemname(snapshot_projfile) // "/cavgs" // STK_EXT, l_smpd, 'cavg')
                call snapshot_proj%add_frcs2os_out(stemname(snapshot_projfile) // '/' // FRCS_FILE, 'frc2D')
                call snapshot_proj%set_cavgs_thumb(snapshot_projfile)
                snapshot_proj%os_ptcl3D = snapshot_proj%os_ptcl2D
                call snapshot_proj%write(snapshot_projfile)
                if(present(snapshot_starfile_base)) then
                    if( DEBUG_HERE ) t = tic()
                    call snapshot_proj%write_mics_star(snapshot_starfile_base // "_micrographs.star")
                    if( DEBUG_HERE ) print *,'ms_export  : ', toc(t); call flush(6); t = tic()
                    call snapshot_proj%write_ptcl2D_star(snapshot_starfile_base // "_particles.star")
                    if( DEBUG_HERE ) print *,'ptcl_export  : ', toc(t); call flush(6)
                endif
                call set_snapshot_time()
                snapshot_last_nptcls = snapshot_proj%os_ptcl2D%count_state_gt_zero()
                call snapshot_proj%kill
            else
                write(logfhandle, '(A)') ">>> FAILED TO WRITE SNAPSHOT"
            end if
            snapshot_complete_jobid = snapshot_jobid
            snapshot_jobid = 0
            snapshot_iteration = 0
        else
            write(logfhandle,'(A,A,A,A)')'>>> WRITING PROJECT ', projfile%to_char(), ' AT: ',cast_time_char(simple_gettime())
            if(present(optics_dir)) then
                call get_latest_optics_map()
                if(lastmap .gt. 0) then
                    mapfileprefix = optics_dir // '/' // OPTICS_MAP_PREFIX // int2str(lastmap)
                    call pool_proj%import_optics_map(mapfileprefix)
                endif
            endif
            if( l_scaling )then
                os_backup = pool_proj%os_cls2D
                ! rescale classes
                if( l_wfilt )then
                    src  = add2fbody(pool_refs, params_glob%ext,WFILT_SUFFIX)
                    dest = add2fbody(cavgsfname,params_glob%ext,WFILT_SUFFIX)
                    call rescale_cavgs(src, dest)
                    src  = add2fbody(pool_refs, params_glob%ext,WFILT_SUFFIX//'_even')
                    dest = add2fbody(cavgsfname,params_glob%ext,WFILT_SUFFIX//'_even')
                    call rescale_cavgs(src, dest)
                    src  = add2fbody(pool_refs, params_glob%ext,WFILT_SUFFIX//'_odd')
                    dest = add2fbody(cavgsfname,params_glob%ext,WFILT_SUFFIX//'_odd')
                    call rescale_cavgs(src, dest)
                endif
                call rescale_cavgs(pool_refs, cavgsfname)
                src  = add2fbody(pool_refs, params_glob%ext,'_even')
                dest = add2fbody(cavgsfname,params_glob%ext,'_even')
                call rescale_cavgs(src, dest)
                src  = add2fbody(pool_refs, params_glob%ext,'_odd')
                dest = add2fbody(cavgsfname,params_glob%ext,'_odd')
                call rescale_cavgs(src, dest)
                call pool_proj%os_out%kill
                call pool_proj%add_cavgs2os_out(cavgsfname, params_glob%smpd, 'cavg', clspath=l_clspath)
                if( l_wfilt )then
                    src = add2fbody(cavgsfname,params_glob%ext,WFILT_SUFFIX)
                    call pool_proj%add_cavgs2os_out(src, params_glob%smpd, 'cavg'//WFILT_SUFFIX, clspath=l_clspath)
                endif
                pool_proj%os_cls2D = os_backup
                call os_backup%kill
                ! rescale frcs
                call frcs%read(string(POOL_DIR)//FRCS_FILE)
                call frcs%pad(params_glob%smpd, params_glob%box)
                call frcs%write(frcsfname)
                call frcs%kill
                call pool_proj%add_frcs2os_out(frcsfname, 'frc2D')
            else
                call pool_proj%os_out%kill
                call pool_proj%add_cavgs2os_out(cavgsfname, params_glob%smpd, 'cavg', clspath=l_clspath)
                if( l_wfilt )then
                    src = add2fbody(cavgsfname,params_glob%ext,WFILT_SUFFIX)
                    call pool_proj%add_cavgs2os_out(src, params_glob%smpd, 'cavg'//WFILT_SUFFIX)
                endif
                call pool_proj%add_frcs2os_out(frcsfname, 'frc2D')
            endif
            pool_proj%os_ptcl3D = pool_proj%os_ptcl2D ! test addition
            call pool_proj%write(projfile)
        endif
        ! 3d field
        pool_proj%os_ptcl3D = pool_proj%os_ptcl2D
        call pool_proj%os_ptcl3D%delete_2Dclustering
        call pool_proj%os_ptcl3D%clean_entry('updatecnt', 'sampled')
        ! write starfiles
        call starproj%export_cls2D(pool_proj)
        if(l_write_star) then
            call starproj_stream%copy_micrographs_optics(pool_proj, verbose=DEBUG_HERE)
            if( DEBUG_HERE ) t = tic()
            call starproj_stream%stream_export_micrographs(pool_proj, params_glob%outdir, optics_set=.true.)
            if( DEBUG_HERE ) print *,'ms_export  : ', toc(t); call flush(6); t = tic()
            call starproj_stream%stream_export_particles_2D(pool_proj, params_glob%outdir, optics_set=.true.)
            if( DEBUG_HERE ) print *,'ptcl_export  : ', toc(t); call flush(6)
        end if
        call pool_proj%os_ptcl3D%kill
        call pool_proj%os_cls2D%delete_entry('stk')

        contains
            
        subroutine set_snapshot_time()
            character(8)  :: date
            character(10) :: time
            character(5)  :: zone
            integer,dimension(8) :: values
            ! using keyword arguments
            call date_and_time(date,time,zone,values)
            call date_and_time(DATE=date,ZONE=zone)
            call date_and_time(TIME=time)
            call date_and_time(VALUES=values)
            write(snapshot_last, '(I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') values(1), '/', values(2), '/', values(3), '_', values(5), ':', values(6)
        end subroutine set_snapshot_time

        subroutine get_latest_optics_map()
            type(string), allocatable :: map_list(:)
            type(string) :: map_i_str, map_str
            integer      :: imap, prefix_len, testmap
            lastmap = 0
            if(optics_dir .ne. "") then
                if(dir_exists(optics_dir)) call simple_list_files(optics_dir%to_char()//'/'//OPTICS_MAP_PREFIX //'*'// TXT_EXT, map_list)
            endif
            if(allocated(map_list)) then
                prefix_len = len(optics_dir%to_char() // '/' // OPTICS_MAP_PREFIX) + 1
                do imap=1, size(map_list)
                    map_str   = map_list(imap)%to_char([prefix_len,map_list(imap)%strlen_trim()])
                    map_i_str = swap_suffix(map_str, "", TXT_EXT)
                    testmap   = map_i_str%to_int()
                    if(testmap > lastmap) lastmap = testmap
                    call map_str%kill
                    call map_i_str%kill
                enddo
                deallocate(map_list)
            endif
        end subroutine get_latest_optics_map
        
    end subroutine write_project_stream2D

    ! Handles user inputted class rejection
    subroutine write_repick_refs(refsout)
        class(string), intent(in) :: refsout
        type(image)  :: img        
        type(string) :: refsin
        integer      :: icls, i
        if( .not. l_stream2D_active ) return
        if( pool_proj%os_cls2D%get_noris() == 0 ) return
        if(repick_iteration .lt. 1) return
        refsin = CAVGS_ITER_FBODY//int2str_pad(repick_iteration,3)//params_glob%ext%to_char()
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

end module simple_stream_cluster2D_utils
