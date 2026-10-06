!@descr: stateless helpers of the stream's 2D pool and chunks: folder clean-up, downscaling, iteration files, set appending, class draws, publication building and naming, snapshot sprite sheets
module simple_stream_refine2D_utils
use simple_core_module_api
use simple_parameters,           only: parameters
use simple_sp_project,           only: sp_project
use simple_optics_maps,          only: import_latest_optics_map
use simple_syslib,               only: get_current_rss_bytes, get_peak_rss_bytes
implicit none

public :: cleanup_root_folder
public :: setup_downscaling
public :: tidy_2Dstream_iter
public :: build_pool_publication
public :: delete_pool_publication
public :: pool_publication_names
public :: snapshot_cavgs_meta
public :: log_rss
public :: draw_new_classes
public :: append_project_sets
private
#include "simple_local_flags.inc"

contains

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

    ! Determines dimensions for downscaling; @p l_scaling, on entry whether to downscale, on return
    ! whether the box was downscaled
    subroutine setup_downscaling( params, l_scaling )
        class(parameters), intent(inout) :: params
        logical,           intent(inout) :: l_scaling
        real    :: SMPD_TARGET = MAX_SMPD  ! target sampling distance
        real    :: smpd, scale_factor, msk_max
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
        ! the mask within the box (D40)
        msk_max = (real(params%box_crop) - COSMSKHALFWIDTH) / 2.
        if( params%msk_crop > msk_max )then
            write(logfhandle,'(A,F8.2,A,F7.1,A)') '>>> MASK DIAMETER ', params%mskdiam, ' A EXCEEDS THE BOX; MASK RADIUS CLAMPED TO ',&
                &real(floor(msk_max)), ' PIXELS'
            params%msk_crop = floor(msk_max)
        endif
    end subroutine setup_downscaling

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
    !! image index in its stack. Stacks never classified are left out; within a published stack,
    !! a particle never updated (a fractional update did not sample it) keeps its row but is
    !! published deselected, since its class is the random one it was given on import. The 3D
    !! particles are prepared as the final project's (2D clustering removed, shifts kept, no update
    !! counts), and with @p optics_dir the newest optics map's groups and optics table are applied
    !! (the pool keeps no optics table of its own). @p nstks is the number of stacks published
    !! (0: nothing to publish, @p pub is empty).
    subroutine build_pool_publication( src, pub, nstks, optics_dir )
        class(sp_project),       intent(inout) :: src
        class(sp_project),       intent(inout) :: pub
        integer,                 intent(out)   :: nstks
        class(string), optional, intent(in)    :: optics_dir
        logical, allocatable :: l_stk(:)
        integer :: nstks_src, istk, jstk, iptcl, jptcl, fromp, top, nptcls, stkind_src, ind_in_stk, lastmap
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
                if( src%os_ptcl2D%get_updatecnt(iptcl) == 0 ) call pub%os_ptcl2D%set_state(jptcl, 0)
            enddo
        enddo
        pub%os_ptcl3D = pub%os_ptcl2D
        call pub%os_ptcl3D%delete_2Dclustering
        call pub%os_ptcl3D%clean_entry('updatecnt', 'sampled')
        if( present(optics_dir) )then
            if( optics_dir%strlen() > 0 ) lastmap = import_latest_optics_map(pub, optics_dir)
        endif
        deallocate(l_stk)
    end subroutine build_pool_publication

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

    !> Appends the sieve sets @p sets to the pool project @p proj, in order: their micrographs, their
    !! stacks renumbered after the project's (the particle ranges follow), and their particles as new
    !! rows that point at their stack in the project, keep their shifts and have no 2D parameters,
    !! update count or fractional-update weight. Sets without micrographs (the sieve's empty final
    !! set) add nothing.
    subroutine append_project_sets( proj, sets )
        class(sp_project), intent(inout) :: proj
        class(sp_project), intent(inout) :: sets(:)
        integer :: nmics_new, nptcls_new, pool_nmics, pool_nptcls, fromp, imic, iset, ind, jmic, nptcls, i, iptcl, jptcl
        nmics_new  = 0
        nptcls_new = 0
        do iset = 1,size(sets)
            nmics_new  = nmics_new  + sets(iset)%os_mic%get_noris()
            nptcls_new = nptcls_new + sets(iset)%os_ptcl2D%get_noris()
        enddo
        if( nmics_new == 0 ) return
        ! room in the project
        pool_nmics  = proj%os_mic%get_noris()
        pool_nptcls = proj%os_ptcl2D%get_noris()
        if( pool_nmics == 0 )then
            call proj%os_mic%new(nmics_new,     is_ptcl=.false.)
            call proj%os_stk%new(nmics_new,     is_ptcl=.false.)
            call proj%os_ptcl2D%new(nptcls_new, is_ptcl=.true.)
            fromp = 1
        else
            call proj%os_mic%reallocate(pool_nmics + nmics_new)
            call proj%os_stk%reallocate(pool_nmics + nmics_new)
            call proj%os_ptcl2D%reallocate(pool_nptcls + nptcls_new)
            fromp = proj%os_stk%get_top(pool_nmics) + 1
        endif
        ! the sets, in order
        imic = pool_nmics
        do iset = 1,size(sets)
            ind = 1
            do jmic = 1,sets(iset)%os_mic%get_noris()
                imic = imic + 1
                call proj%os_mic%transfer_ori(imic, sets(iset)%os_mic, jmic)
                call proj%os_stk%transfer_ori(imic, sets(iset)%os_stk, jmic)
                nptcls = sets(iset)%os_stk%get_int(jmic, 'nptcls')
                call proj%os_stk%set(imic, 'fromp', fromp)
                call proj%os_stk%set(imic, 'top',   fromp + nptcls - 1)
                !$omp parallel do private(i,iptcl,jptcl) default(shared) proc_bind(close)
                do i = 1,nptcls
                    iptcl = fromp + i - 1
                    jptcl = ind   + i - 1
                    call proj%os_ptcl2D%transfer_ori(iptcl, sets(iset)%os_ptcl2D, jptcl)
                    call proj%os_ptcl2D%set_stkind(iptcl, imic)
                    call proj%os_ptcl2D%set(iptcl, 'updatecnt', 0)
                    call proj%os_ptcl2D%set(iptcl, 'frac',      0.)
                    call proj%os_ptcl2D%delete_2Dclustering(iptcl, keepshifts=.true.)
                enddo
                !$omp end parallel do
                ind   = ind   + nptcls
                fromp = fromp + nptcls
            enddo
        enddo
    end subroutine append_project_sets

    ! Classes for the never-updated particles among the first @p nptcls of @p spproj, drawn among
    ! its populated classes (of the first @p ncls) from a generator seeded with the iteration
    ! @p iter, in one thread, so a run is reproducible; the process's generator state is restored
    ! after. With no populated class they keep the class they have.
    subroutine draw_new_classes( spproj, nptcls, iter, ncls )
        use iso_fortran_env, only: int64
        class(sp_project), intent(inout) :: spproj
        integer,           intent(in)    :: nptcls, iter, ncls
        integer, allocatable :: clspops(:), populated(:), saved_seed(:), seed(:)
        integer :: ncls_drawn, iptcl, nseed, i
        clspops    = spproj%os_cls2D%get_all_asint('pop')
        ncls_drawn = min(ncls, size(clspops))
        if( ncls_drawn < 1 ) return
        populated  = pack([(i, i=1,ncls_drawn)], clspops(:ncls_drawn) > 0)
        if( size(populated) == 0 ) return
        call random_seed(size=nseed)
        allocate(saved_seed(nseed), seed(nseed))
        call random_seed(get=saved_seed)
        do i = 1,nseed
            seed(i) = int(modulo(int(iter,int64) * 7919_int64 + 104729_int64 * int(i-1,int64), int(huge(0)-1,int64)) + 1_int64)
        enddo
        call random_seed(put=seed)
        do iptcl = 1,min(nptcls, spproj%os_ptcl2D%get_noris())
            if( spproj%os_ptcl2D%get_state(iptcl) == 0 ) cycle
            if( spproj%os_ptcl2D%get_updatecnt(iptcl) /= 0 ) cycle
            call spproj%os_ptcl2D%set_class(iptcl, populated(irnd_uni(size(populated))))
        enddo
        call random_seed(put=saved_seed)
    end subroutine draw_new_classes

end module simple_stream_refine2D_utils
