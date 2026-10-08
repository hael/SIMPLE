!@descr: stateless helpers of the stream's 2D pool: folder clean-up, iteration files, set appending, class draws, publication building and naming, snapshot sprite sheets
module simple_stream_refine2D_utils
use simple_core_module_api
use simple_sp_project,           only: sp_project
use simple_optics_maps,          only: import_latest_optics_map
use simple_class_frcs,           only: class_frcs
use simple_image,                only: image
implicit none

public :: cleanup_root_folder
public :: tidy_2Dstream_iter
public :: build_pool_publication
public :: pool_publication_nselected
public :: build_sieve_publication
public :: combine_sieve_classes
public :: delete_pool_publication
public :: pool_publication_names
public :: snapshot_cavgs_meta
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
    !> The particles a publication of @p src selects (build_pool_publication): those an iteration
    !! has updated and the pool keeps (state > 0).
    integer function pool_publication_nselected( src ) result( n )
        class(sp_project), intent(in) :: src
        integer :: iptcl
        n = 0
        do iptcl = 1,src%os_ptcl2D%get_noris()
            if( src%os_ptcl2D%get_state(iptcl) <= 0 ) cycle
            if( src%os_ptcl2D%get_updatecnt(iptcl) > 0 ) n = n + 1
        enddo
    end function pool_publication_nselected

    !> The first publication for 3D built from the sieve's own 2D (sieve_ini3D), from the sieve's
    !! sets @p sets, each read whole (micrographs, stacks, particles, classes): their stacks
    !! renumbered from 1 and their particles with their image index in the stack (indstk), as
    !! build_pool_publication sets them; the sieve's 2D parameters and selection, each set's class
    !! labels offset by the classes of the sets before it (a particle without a class of its set is
    !! deselected); and the sets' class tables concatenated in set order, the sieve's states kept and
    !! each population counted from the selected particles. The 3D particles are prepared as
    !! build_pool_publication prepares them. @p nstks is the number of stacks published (0:
    !! nothing), @p ncls the classes. The class averages go with it through combine_sieve_classes.
    subroutine build_sieve_publication( sets, pub, nstks, ncls )
        class(sp_project), intent(inout) :: sets(:)
        class(sp_project), intent(inout) :: pub
        integer,           intent(out)   :: nstks, ncls
        integer, allocatable :: pops(:)
        integer :: iset, ifirst, istk, jstk, iptcl, jptcl, nptcls, fromp, top, stkind, ind_in_stk
        integer :: offset, ncls_set, icls
        call pub%kill
        nstks  = 0
        ncls   = 0
        nptcls = 0
        ifirst = 0
        do iset = 1,size(sets)
            ncls = ncls + sets(iset)%os_cls2D%get_noris()
            if( sets(iset)%os_stk%get_noris() == 0 ) cycle
            if( sets(iset)%os_mic%get_noris() /= sets(iset)%os_stk%get_noris() ) THROW_HARD('build_sieve_publication: # micrographs /= # stacks')
            if( ifirst == 0 ) ifirst = iset
            nstks  = nstks  + sets(iset)%os_stk%get_noris()
            nptcls = nptcls + sets(iset)%os_ptcl2D%get_noris()
        enddo
        if( nstks == 0 )then
            ncls = 0
            return
        endif
        pub%projinfo  = sets(ifirst)%projinfo
        pub%compenv   = sets(ifirst)%compenv
        pub%jobproc   = sets(ifirst)%jobproc
        pub%os_optics = sets(ifirst)%os_optics
        call pub%os_mic%new(nstks,    is_ptcl=.false.)
        call pub%os_stk%new(nstks,    is_ptcl=.false.)
        call pub%os_ptcl2D%new(nptcls, is_ptcl=.true.)
        if( ncls > 0 ) call pub%os_cls2D%new(ncls, is_ptcl=.false.)
        jstk   = 0
        jptcl  = 0
        offset = 0
        do iset = 1,size(sets)
            ncls_set = sets(iset)%os_cls2D%get_noris()
            do icls = 1,ncls_set
                call pub%os_cls2D%transfer_ori(offset + icls, sets(iset)%os_cls2D, icls)
            enddo
            do istk = 1,sets(iset)%os_stk%get_noris()
                jstk  = jstk + 1
                fromp = sets(iset)%os_stk%get_fromp(istk)
                top   = sets(iset)%os_stk%get_top(istk)
                call pub%os_mic%transfer_ori(jstk, sets(iset)%os_mic, istk)
                call pub%os_stk%transfer_ori(jstk, sets(iset)%os_stk, istk)
                call pub%os_stk%set(jstk, 'fromp', jptcl + 1)
                call pub%os_stk%set(jstk, 'top',   jptcl + top - fromp + 1)
                do iptcl = fromp,top
                    jptcl = jptcl + 1
                    ! the image index in the stack, as the set maps it, before the ranges change
                    call sets(iset)%map_ptcl_ind2stk_ind('ptcl2D', iptcl, stkind, ind_in_stk)
                    call pub%os_ptcl2D%transfer_ori(jptcl, sets(iset)%os_ptcl2D, iptcl)
                    call pub%os_ptcl2D%set_stkind(jptcl, jstk)
                    call pub%os_ptcl2D%set(jptcl, 'indstk', ind_in_stk)
                    icls = pub%os_ptcl2D%get_class(jptcl)
                    if( icls >= 1 .and. icls <= ncls_set )then
                        call pub%os_ptcl2D%set_class(jptcl, offset + icls)
                    else
                        call pub%os_ptcl2D%set_class(jptcl, 0)
                        call pub%os_ptcl2D%set_state(jptcl, 0)
                    endif
                enddo
            enddo
            offset = offset + ncls_set
        enddo
        ! the populations of the selected particles
        if( ncls > 0 )then
            allocate(pops(ncls), source=0)
            do jptcl = 1,nptcls
                if( pub%os_ptcl2D%get_state(jptcl) <= 0 ) cycle
                icls = pub%os_ptcl2D%get_class(jptcl)
                if( icls >= 1 .and. icls <= ncls ) pops(icls) = pops(icls) + 1
            enddo
            do icls = 1,ncls
                call pub%os_cls2D%set(icls, 'class', icls)
                call pub%os_cls2D%set(icls, 'pop',   pops(icls))
            enddo
        endif
        pub%os_ptcl3D = pub%os_ptcl2D
        call pub%os_ptcl3D%delete_2Dclustering
        call pub%os_ptcl3D%clean_entry('updatecnt', 'sampled')
    end subroutine build_sieve_publication

    !> The sieve sets' class averages, their even and odd halves and their FRCs, concatenated in set
    !! order into @p cavgsfname (its halves beside it) and @p frcsfname, the files of the publication
    !! build_sieve_publication builds. @p stks(i) is set i's class-average stack ('' for a set
    !! without classes) with its halves beside it (_even, _odd), @p frcs(i) its FRCs and @p ncls(i)
    !! its classes; @p smpd is the class averages' pixel size. .false., with nothing written, when a
    !! file is missing, a stack does not hold its set's classes, or the sets' boxes or FRC sizes
    !! differ.
    function combine_sieve_classes( stks, frcs, ncls, smpd, cavgsfname, frcsfname ) result( l_ok )
        class(string), intent(in) :: stks(:), frcs(:)
        integer,       intent(in) :: ncls(:)
        real,          intent(in) :: smpd
        class(string), intent(in) :: cavgsfname, frcsfname
        logical :: l_ok
        type(class_frcs)  :: frcs_set, frcs_all
        type(image)       :: img
        type(string)      :: src(3), dst(3), ext, dst_ext
        real, allocatable :: frc(:)
        integer :: iset, ldim(3), n, ntot, box, frc_box, filtsz, ihalf, icls, k
        l_ok    = .false.
        ntot    = 0
        box     = 0
        frc_box = 0
        filtsz  = 0
        do iset = 1,size(stks)
            if( stks(iset)%strlen() == 0 )then
                if( ncls(iset) > 0 ) return
                cycle
            endif
            ext = string('.')//fname2ext(stks(iset))
            if( .not. file_exists(stks(iset)) ) return
            if( .not. file_exists(add2fbody(stks(iset), ext%to_char(), '_even')) ) return
            if( .not. file_exists(add2fbody(stks(iset), ext%to_char(), '_odd'))  ) return
            if( .not. file_exists(frcs(iset)) ) return
            call find_ldim_nptcls(stks(iset), ldim, n)
            if( n /= ncls(iset) ) return
            if( box == 0 ) box = ldim(1)
            if( ldim(1) /= box ) return
            call frcs_set%read(frcs(iset))
            if( frcs_set%get_ncls() /= n ) return
            if( frc_box == 0 )then
                frc_box = frcs_set%get_box()
                filtsz  = frcs_set%get_filtsz()
            endif
            if( frcs_set%get_box() /= frc_box .or. frcs_set%get_filtsz() /= filtsz ) return
            call frcs_set%kill
            ntot = ntot + n
        enddo
        if( ntot == 0 ) return
        ! the stacks, main and halves, in set order
        dst_ext = string('.')//fname2ext(cavgsfname)
        dst(1)  = cavgsfname
        dst(2)  = add2fbody(cavgsfname, dst_ext%to_char(), '_even')
        dst(3)  = add2fbody(cavgsfname, dst_ext%to_char(), '_odd')
        do ihalf = 1,3
            if( file_exists(dst(ihalf)) ) call del_file(dst(ihalf))
        enddo
        call img%new([box, box, 1], smpd)
        k = 0
        do iset = 1,size(stks)
            if( stks(iset)%strlen() == 0 ) cycle
            ext    = string('.')//fname2ext(stks(iset))
            src(1) = stks(iset)
            src(2) = add2fbody(stks(iset), ext%to_char(), '_even')
            src(3) = add2fbody(stks(iset), ext%to_char(), '_odd')
            do icls = 1,ncls(iset)
                do ihalf = 1,3
                    call img%read(src(ihalf), icls)
                    call img%write(dst(ihalf), k + icls)
                enddo
            enddo
            k = k + ncls(iset)
        enddo
        call img%kill
        ! the FRCs, class by class, at the sets' FRC box (its pixel size follows from the averages')
        call frcs_all%new(ntot, frc_box, smpd * real(box) / real(frc_box))
        allocate(frc(filtsz))
        k = 0
        do iset = 1,size(stks)
            if( stks(iset)%strlen() == 0 ) cycle
            call frcs_set%read(frcs(iset))
            do icls = 1,ncls(iset)
                call frcs_set%frc_getter(icls, frc)
                call frcs_all%set_frc(k + icls, frc)
            enddo
            call frcs_set%kill
            k = k + ncls(iset)
        enddo
        call frcs_all%write(frcsfname)
        call frcs_all%kill
        l_ok = .true.
    end function combine_sieve_classes

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
