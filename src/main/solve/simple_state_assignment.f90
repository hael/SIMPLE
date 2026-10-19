!@descr: multistate 3D by class-average assignment, shared by stream p07 and solve3D_pyramid: the rules, the reference check, the matching's mapping onto the particles, the local refinement's seeds and result, and the jobs' command lines
!==============================================================================
! MODULE: simple_state_assignment
!
! PURPOSE:
!   The method of stream p07 (simple_stream_stage_solve3D) without its stage, so the batch program
!   solve3D_pyramid (simple_commanders_solve3D_pyramid) runs the same steps on the same rules
!   (doc/policies/stream/stream_3D_ingestion_policy.md, doc/policies/3D/solve3D_pyramid_policy.md):
!     - the route to the reference volumes: solve2D of a first set into NCLS2D_FIRST classes, the
!       pool model's selection of them (select_route_classes), cls_expansion into expansion_budget
!       sub-classes, solve3D_cavgs on the sub-class averages;
!     - the check of solve3D_cavgs' volumes (check_state_refs, refs_usable, docked_pair_cc), with
!       the volumes, their record and the class averages per state written into one folder;
!     - the assignment: assign_cavgs_states' matching mapped onto the particles (map_assignment:
!       reject_by_fit, flag_ambiguous, map2ptcls), the matched class averages per state
!       (write_matched_cavgs); for solve3D_pyramid, a solve3D_cavgs of each state's class averages
!       alone (select_state_classes), its volume docked onto the state's reference
!       (docked_cc_common_grid) and its poses mapped onto the state's particles (map_state_poses);
!     - the local refinement: its seeds (set_refine_seeds), its project with the states renumbered
!       (write_refine_project) and its result mapped back (map_refine_result);
!     - the command lines of the jobs (set_*_cline): the keys that make the method, set alike by
!       both callers, which add the resources (threads, parts) themselves.
!   Free procedures without state: the callers own the projects, the folders and the jobs.
!==============================================================================
module simple_state_assignment
use simple_defs,                   only: logfhandle, COSMSKHALFWIDTH
use simple_defs_fname,             only: MRC_EXT, JPG_EXT
use simple_error,                  only: simple_exception
use simple_string,                 only: string
use simple_string_utils,           only: int2str, int2str_pad
use simple_fileio,                 only: add2fbody, del_file, file_exists, fname2ext, get_fpath, simple_abspath,&
                                        &simple_copy_file
use simple_syslib,                 only: simple_mkdir
use simple_math,                   only: round2even
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_imghead,                only: find_ldim_nptcls, find_img_smpd
use simple_refine3D_fnames,        only: refine3D_oris_heatmap_fname, refine3D_reprojs_fname
use simple_oris,                   only: oris
use simple_sym,                    only: sym
use simple_cmdline,                only: cmdline
use simple_sp_project,             only: sp_project
use simple_eo_partition,           only: eo_transaction
use simple_image,                  only: image
use simple_imgarr_utils,           only: dealloc_imgarr, read_cavgs_into_imgarr
use simple_cavg_quality_model,     only: cavg_quality_model, CAVG_QUALITY_MODEL_POOL_DEFAULT
use simple_cavg_quality_types,     only: cavg_quality_result
use simple_cavg_quality_selection, only: score_project_cavgs, write_cavg_selection_stacks
use simple_dock_vols,              only: dock_vols
use simple_gui_utils,              only: mrc2jpeg_tiled
implicit none

! the method's constants
public :: MIN_PTCLS_PER_STATE, NCLS2D_FIRST, NCLS_EXPANDED, LPSTART_CAVGS3D, LPSTOP_CAVGS3D
public :: EO_ANGDIST_MAX, EO_OVERLAP_MIN, STATE_FRAC_MIN, NSTATES_MIN, REFS_PAIR_CC_MAX, DOCK_HP_REFS
public :: LP_ASSIGN, NSPACE_ASSIGN, FIT_NMAD, FIT_REJECT_MAX, MARGIN_AMBIGUOUS
! the jobs' folders and project names, and the files written
public :: SOLVE2D_DIR, CLSEXP_DIR, CAVGS3D_DIR, ASSIGN_LABEL, REFINE_DIR, REFS_DIR, REFS_RECORD
public :: QUALITY_DIR, SELECTED_CAVGS, REJECTED_CAVGS
! the values of a check and of an assignment
public :: refs_check, assign_stats
! rules
public :: refs_usable, reject_by_fit, flag_ambiguous, compact_states, expansion_budget, count_active_classes
public :: fit_mskdiam, median_of, usable_job_2D, cut_whole_stacks
! steps
public :: score_cavgs_with_model, select_route_classes, model_rejected_rows
public :: check_state_refs, docked_pair_cc, docked_cc_common_grid, read_refs_record
public :: map_assignment, write_matched_cavgs, select_state_classes, map_state_poses
public :: set_refine_seeds, write_refine_project, map_refine_result
public :: write_state_cavgs, cavgs_fbody, log_jpeg
! the jobs' command lines
public :: set_solve2D_route_cline, set_clsexp_cline, set_cavgs3D_cline, set_state_cavgs3D_cline, set_assign_cline
public :: set_refine_cont_cline
private
#include "simple_local_flags.inc"

! a job needs this many selected particles per state
integer, parameter :: MIN_PTCLS_PER_STATE = 5
! the route: solve2D into NCLS2D_FIRST classes, then cls_expansion of the classes the pool model
! keeps into NCLS_EXPANDED sub-classes in all (each assignment's expansion into as many); never
! fewer than the classes split, which cls_expansion refuses (expansion_budget)
integer, parameter :: NCLS2D_FIRST        = 150
integer, parameter :: NCLS_EXPANDED       = 500
! solve3D_cavgs' low-pass ramp (A), fixed (force_lp_range) rather than derived from the class FRCs
real,    parameter :: LPSTART_CAVGS3D     = 20.
real,    parameter :: LPSTOP_CAVGS3D      = 12.
! the check of the reference volumes (provisional, until a validation run sets them)
real,    parameter :: EO_ANGDIST_MAX      = 50.  ! the class averages' even/odd angular distance at most (degrees)
real,    parameter :: EO_OVERLAP_MIN      = 70.  ! their even/odd state agreement at least (%)
real,    parameter :: STATE_FRAC_MIN      = 0.05 ! a state counts when it holds this fraction of the class averages
integer, parameter :: NSTATES_MIN         = 2    ! states that count, at least
real,    parameter :: REFS_PAIR_CC_MAX    = 0.97 ! two states' correlation after docking, at most (distinct states)
real,    parameter :: DOCK_HP_REFS        = 100. ! the docking band: from this high-pass to the matching's limit (A)
! the assignment (provisional)
real,    parameter :: LP_ASSIGN           = 15.  ! the matching's low-pass limit (A), the same for every state
integer, parameter :: NSPACE_ASSIGN       = 1000 ! the matching's directions per state, about 6 degrees apart
real,    parameter :: FIT_NMAD            = 3.   ! a best correlation below the median less this many MADs is rejected
real,    parameter :: FIT_REJECT_MAX      = 0.30 ! of the (sub-)classes matched, rejected by fit at most
real,    parameter :: MARGIN_AMBIGUOUS    = 0.02 ! a margin below this flags a (sub-)class ambiguous

character(len=*), parameter :: SOLVE2D_DIR      = 'solve2D'
character(len=*), parameter :: CLSEXP_DIR       = 'cls_expansion'
character(len=*), parameter :: CAVGS3D_DIR      = 'solve3D_cavgs'
character(len=*), parameter :: ASSIGN_LABEL     = 'assign_cavgs_states'
character(len=*), parameter :: REFINE_DIR       = 'refine3D'
character(len=*), parameter :: REFS_DIR         = 'state_refs'
character(len=*), parameter :: REFS_RECORD      = 'state_refs.txt'
character(len=*), parameter :: QUALITY_DIR      = 'quality_selection'
character(len=*), parameter :: SELECTED_CAVGS   = 'quality_selected_cavgs'
character(len=*), parameter :: REJECTED_CAVGS   = 'quality_rejected_cavgs'

!> the check of solve3D_cavgs' volumes (check_state_refs)
type :: refs_check
    integer :: nheld       = 0     !< the states that hold class averages, renumbered 1..nheld
    integer :: box         = 0     !< the volumes' box and pixel size (A)
    real    :: smpd        = 0.
    real    :: angdist     = -1.   !< the even/odd angular distance checked (degrees): the median over the class
                                   !! averages whose halves are in one state, else the mean over all (-1: none)
    real    :: angdist_all = -1.   !< the mean over every class (degrees)
    real    :: overlap     = -1.   !< the even/odd state agreement (%)
    real    :: cc_max      = -1.   !< the largest correlation of two docked states (-1: not compared)
    integer, allocatable :: counts(:)     !< class averages per state
    real,    allocatable :: cc_pairs(:,:) !< each pair's docked correlation, the better hand (-1: not compared)
    logical :: l_vols      = .false. !< a volume for every state that holds class averages
    logical :: l_pass      = .false. !< the check passed (refs_usable)
end type refs_check

!> the values of an assignment (map_assignment), for the log and the reports
type :: assign_stats
    integer :: ncls       = 0  !< (sub-)classes matched
    integer :: ncand      = 0  !< candidates: with a state and particles
    integer :: nrej       = 0  !< rejected by fit
    integer :: namb       = 0  !< flagged ambiguous (assigned)
    integer :: nmax       = 0  !< the cap on the rejection
    real    :: med        = 0. !< the candidates' median best correlation, its median absolute deviation and
    real    :: mad        = 0. !! the threshold below which one is rejected
    real    :: thres      = 0.
    integer :: nassigned  = 0  !< particles given a state
    integer :: ndesel     = 0  !< particles selected before and deselected
    logical :: l_model    = .false. !< the pool model's verdict was given: the next three are counted
    integer :: nboth      = 0  !< particles rejected by the fit and the model
    integer :: nmodel     = 0  !< by the model alone (not applied)
    integer :: nfit       = 0  !< by the fit alone
    integer, allocatable :: counts(:) !< particles per state
end type assign_stats

contains

    !---------------- rules ----------------

    !> .true. when reference volumes may be used: the class averages' even/odd angular distance
    !! @p angdist (degrees) at most EO_ANGDIST_MAX (not negative: measured), their even/odd state
    !! agreement @p overlap (%) at least EO_OVERLAP_MIN, the largest correlation between two docked
    !! states @p cc_max at most REFS_PAIR_CC_MAX (distinct states; negative: not compared), and at
    !! least NSTATES_MIN states each holding STATE_FRAC_MIN of the class averages (@p counts per
    !! state).
    pure logical function refs_usable( angdist, overlap, counts, cc_max )
        real,    intent(in) :: angdist, overlap, cc_max
        integer, intent(in) :: counts(:)
        integer :: ntot
        refs_usable = .false.
        if( angdist < 0. .or. angdist > EO_ANGDIST_MAX ) return
        if( overlap < EO_OVERLAP_MIN ) return
        if( cc_max > REFS_PAIR_CC_MAX ) return
        ntot = sum(counts)
        if( ntot < 1 ) return
        refs_usable = count(real(counts) >= STATE_FRAC_MIN * real(ntot)) >= NSTATES_MIN
    end function refs_usable

    !> Of the candidate classes @p l_cand, those whose best correlation @p corrs falls below the
    !! candidates' median less FIT_NMAD times their median absolute deviation are rejected
    !! (@p l_rej), the worst first, at most FIT_REJECT_MAX of the candidates (@p nmax_out). The
    !! median, the deviation and the threshold are returned for the log (0 without candidates).
    pure subroutine reject_by_fit( corrs, l_cand, l_rej, med_out, mad_out, thres_out, nmax_out )
        real,                 intent(in)    :: corrs(:)
        logical,              intent(in)    :: l_cand(:)
        logical, allocatable, intent(inout) :: l_rej(:)
        real,    optional,    intent(out)   :: med_out, mad_out, thres_out
        integer, optional,    intent(out)   :: nmax_out
        real, allocatable :: vals(:)
        real    :: med, mad, thres, worst
        integer :: ncand, nmax, i, iworst
        if( allocated(l_rej) ) deallocate(l_rej)
        allocate(l_rej(size(corrs)), source=.false.)
        if( present(med_out)   ) med_out   = 0.
        if( present(mad_out)   ) mad_out   = 0.
        if( present(thres_out) ) thres_out = 0.
        if( present(nmax_out)  ) nmax_out  = 0
        ncand = count(l_cand)
        if( ncand < 1 ) return
        vals  = pack(corrs, l_cand)
        med   = median_of(vals)
        vals  = abs(vals - med)
        mad   = median_of(vals)
        thres = med - FIT_NMAD * mad
        nmax  = int(FIT_REJECT_MAX * real(ncand))
        if( present(med_out)   ) med_out   = med
        if( present(mad_out)   ) mad_out   = mad
        if( present(thres_out) ) thres_out = thres
        if( present(nmax_out)  ) nmax_out  = nmax
        ! the worst first, while below the threshold and under the cap
        do while( count(l_rej) < nmax )
            iworst = 0
            worst  = thres
            do i = 1,size(corrs)
                if( .not. l_cand(i) .or. l_rej(i) ) cycle
                if( corrs(i) < worst )then
                    worst  = corrs(i)
                    iworst = i
                endif
            enddo
            if( iworst == 0 ) exit
            l_rej(iworst) = .true.
        enddo
    end subroutine reject_by_fit

    !> The classes @p l_kept whose @p margins (best less second-best state) fall below
    !! MARGIN_AMBIGUOUS: ambiguous, assigned as any other.
    pure function flag_ambiguous( margins, l_kept ) result( l_amb )
        real,    intent(in) :: margins(:)
        logical, intent(in) :: l_kept(:)
        logical :: l_amb(size(margins))
        l_amb = l_kept .and. margins < MARGIN_AMBIGUOUS
    end function flag_ambiguous

    !> The positive states of @p states renumbered from 1 in ascending order, in place; @p newind(s)
    !! is state s's new number (0 for a state none holds) and @p nstates its count.
    pure subroutine compact_states( states, newind, nstates )
        integer,              intent(inout) :: states(:)
        integer, allocatable, intent(inout) :: newind(:)
        integer,              intent(out)   :: nstates
        integer :: s, smax, i
        nstates = 0
        smax    = 0
        do i = 1,size(states)
            smax = max(smax, states(i))
        enddo
        if( allocated(newind) ) deallocate(newind)
        allocate(newind(max(smax, 0)), source=0)
        do s = 1,smax
            if( any(states == s) )then
                nstates   = nstates + 1
                newind(s) = nstates
            endif
        enddo
        do i = 1,size(states)
            if( states(i) > 0 ) states(i) = newind(states(i))
        enddo
    end subroutine compact_states

    !> The total number of sub-classes cls_expansion splits @p nsplit classes into (its ncls):
    !! NCLS_EXPANDED, and never fewer than the classes, which cls_expansion refuses.
    pure integer function expansion_budget( nsplit ) result( nsub )
        integer, intent(in) :: nsplit
        nsub = max(NCLS_EXPANDED, nsplit)
    end function expansion_budget

    !> The classes cls_expansion splits: the distinct positive class labels @p classes of the rows
    !! with a positive state @p states (its active particles).
    pure integer function count_active_classes( states, classes ) result( n )
        integer, intent(in) :: states(:), classes(:)
        logical, allocatable :: l_seen(:)
        integer :: i, cmax
        n    = 0
        cmax = 0
        do i = 1,min(size(states), size(classes))
            if( states(i) > 0 ) cmax = max(cmax, classes(i))
        enddo
        if( cmax < 1 ) return
        allocate(l_seen(cmax), source=.false.)
        do i = 1,min(size(states), size(classes))
            if( states(i) > 0 .and. classes(i) > 0 ) l_seen(classes(i)) = .true.
        enddo
        n = count(l_seen)
    end function count_active_classes

    !> The mask diameter (A) for class averages of @p box pixels at @p smpd: @p mskdiam, or the
    !! box default when it is not positive or too large for the box, as parameters falls back.
    pure real function fit_mskdiam( mskdiam, box, smpd )
        real,    intent(in) :: mskdiam, smpd
        integer, intent(in) :: box
        real :: msk, msk_default
        fit_mskdiam = mskdiam
        msk_default = round2even((real(box) - COSMSKHALFWIDTH) / 2.)
        if( msk_default <= 0. ) return
        msk = 0.
        if( smpd > 0. ) msk = round2even((mskdiam / smpd) / 2.)
        if( msk < 0.1 .or. msk > msk_default ) fit_mskdiam = (real(box) - COSMSKHALFWIDTH) * smpd
    end function fit_mskdiam

    !> The median of @p x (none empty).
    pure real function median_of( x )
        real, intent(in) :: x(:)
        real    :: y(size(x)), t
        integer :: n, i, j
        n = size(x)
        y = x
        do i = 2,n
            t = y(i)
            j = i - 1
            do while( j >= 1 )
                if( y(j) <= t ) exit
                y(j+1) = y(j)
                j = j - 1
            enddo
            y(j+1) = t
        enddo
        if( mod(n, 2) == 1 )then
            median_of = y((n + 1) / 2)
        else
            median_of = 0.5 * (y(n / 2) + y(n / 2 + 1))
        endif
    end function median_of

    !> .true. when the 2D project @p res (solve2D's, cls_expansion's or a classified input) can be
    !! used: @p n particles, one class average per class, and the class averages, their halves and
    !! the FRCs on disk.
    logical function usable_job_2D( res, n ) result( l_ok )
        class(sp_project), intent(inout) :: res
        integer,           intent(in)    :: n
        type(string) :: stk, frcs, ext
        real         :: smpd
        integer      :: ncls
        l_ok = .false.
        ncls = 0
        call res%get_cavgs_stk(stk, ncls, smpd, fail=.false.)
        if( res%os_ptcl2D%get_noris() /= n .or. ncls <= 0 .or. ncls /= res%os_cls2D%get_noris() ) return
        call res%get_frcs(frcs, 'frc2D', fail=.false.)
        ext  = string('.')//fname2ext(stk)
        l_ok = file_exists(stk) .and. file_exists(frcs)
        if( l_ok ) l_ok = file_exists(add2fbody(stk, ext%to_char(), '_even')) .and.&
            &file_exists(add2fbody(stk, ext%to_char(), '_odd'))
    end function usable_job_2D

    !> The first @p nmax of the rows @p selected, whole stacks in row (arrival) order: a stack
    !! (@p stkind, the rows of one stack consecutive) is taken with all its selected rows while they
    !! fit, and the stack that crosses @p nmax ends the cut, so the cut may hold fewer; when the
    !! first stack with selected rows alone exceeds @p nmax, its first @p nmax selected rows. Every
    !! selected row when they are @p nmax at most.
    pure function cut_whole_stacks( selected, stkind, nmax ) result( l_in )
        logical, intent(in) :: selected(:)
        integer, intent(in) :: stkind(:), nmax
        logical :: l_in(size(selected))
        integer :: i, j, k, ntaken, nstk
        l_in = .false.
        if( nmax < 1 ) return
        if( count(selected) <= nmax )then
            l_in = selected
            return
        endif
        ntaken = 0
        i      = 1
        do while( i <= size(selected) )
            ! rows i..j: one stack
            j = i
            do while( j < size(selected) )
                if( stkind(j + 1) /= stkind(i) ) exit
                j = j + 1
            enddo
            nstk = count(selected(i:j))
            if( ntaken + nstk > nmax )then
                if( ntaken == 0 )then
                    ! the first stack alone exceeds the count: its first nmax selected rows
                    do k = i,j
                        if( .not. selected(k) ) cycle
                        if( ntaken == nmax ) exit
                        l_in(k) = .true.
                        ntaken  = ntaken + 1
                    enddo
                endif
                return
            endif
            l_in(i:j) = selected(i:j)
            ntaken    = ntaken + nstk
            i         = j + 1
        enddo
    end function cut_whole_stacks

    !---------------- the pool model ----------------

    !> The class averages of @p proj scored with the quality model @p preset, with the mask diameter
    !! @p mskdiam (A) fitted to their box: @p cls_states is the model's verdict per class (0
    !! rejected). The selected and rejected class averages are written with JPEGs in @p dir.
    !! Nothing in @p proj changes. .false. when @p proj has no class averages.
    logical function score_cavgs_with_model( proj, preset, mskdiam, dir, cls_states ) result( l_ok )
        class(sp_project),    intent(inout) :: proj
        character(len=*),     intent(in)    :: preset
        real,                 intent(in)    :: mskdiam
        class(string),        intent(in)    :: dir
        integer, allocatable, intent(inout) :: cls_states(:)
        type(cavg_quality_model)  :: model
        type(cavg_quality_result) :: quality
        type(image), allocatable  :: cavg_imgs(:)
        type(string)              :: stk
        real                      :: smpd, box_stk
        integer                   :: ncls
        l_ok    = .false.
        box_stk = 0.
        call proj%get_cavgs_stk(stk, ncls, smpd, fail=.false., box=box_stk) ! the os_out entry's box, as a real
        if( ncls <= 0 .or. proj%os_cls2D%get_noris() == 0 ) return
        call model%init_preset(preset)
        call score_project_cavgs(proj, model, fit_mskdiam(mskdiam, nint(box_stk), smpd), cavg_imgs, quality, smpd=smpd)
        call model%kill
        if( size(cavg_imgs) /= proj%os_cls2D%get_noris() ) THROW_HARD('# class averages /= # cls2D entries in '//dir%to_char())
        cls_states = quality%states
        call quality%kill
        call simple_mkdir(dir)
        call write_cavg_selection_stacks(cavg_imgs, cls_states, dir//'/'//SELECTED_CAVGS//MRC_EXT,&
            &dir//'/'//REJECTED_CAVGS//MRC_EXT)
        call mrc2jpeg_tiled(dir//'/'//SELECTED_CAVGS//MRC_EXT, dir//'/'//SELECTED_CAVGS//JPG_EXT)
        call mrc2jpeg_tiled(dir//'/'//REJECTED_CAVGS//MRC_EXT, dir//'/'//REJECTED_CAVGS//JPG_EXT)
        call log_jpeg('QUALITY SELECTED CLASS AVERAGES', dir//'/'//SELECTED_CAVGS//JPG_EXT)
        call log_jpeg('QUALITY REJECTED CLASS AVERAGES', dir//'/'//REJECTED_CAVGS//JPG_EXT)
        write(logfhandle,'(A,A,A,I6,A,I6)') '>>> CAVG QUALITY (', preset, ') SELECTED / REJECTED : ', count(cls_states > 0), ' / ',&
            &count(cls_states <= 0)
        call dealloc_imgarr(cavg_imgs)
        l_ok = .true.
    end function score_cavgs_with_model

    !> The route's classes of @p proj (its first set, classified): with @p l_model, the pool
    !! model's selection (score_cavgs_with_model, the stacks into @p dir) mapped to the particles by
    !! class, so the particles of a rejected class are deselected; when the model cannot score the
    !! classes or keeps none, and without @p l_model, every populated class is selected and no
    !! particle deselected.
    subroutine select_route_classes( proj, l_model, mskdiam, dir )
        class(sp_project), intent(inout) :: proj
        logical,           intent(in)    :: l_model
        real,              intent(in)    :: mskdiam
        class(string),     intent(in)    :: dir
        integer, allocatable :: states2D(:), states3D(:), cls_states(:)
        logical :: l_kept
        states2D = proj%os_ptcl2D%get_all_asint('state')
        states3D = proj%os_ptcl3D%get_all_asint('state')
        l_kept   = .false.
        if( l_model )then
            l_kept = score_cavgs_with_model(proj, CAVG_QUALITY_MODEL_POOL_DEFAULT, mskdiam, dir, cls_states)
            if( l_kept ) call proj%map_cavgs_selection(cls_states)
            if( l_kept ) l_kept = any(proj%os_cls2D%get_all_asint('state') > 0)
            if( .not. l_kept ) write(logfhandle,'(A)') '>>> WARNING: THE POOL MODEL KEPT NO CLASS; EVERY POPULATED CLASS IS USED'
        else
            write(logfhandle,'(A)') '>>> NO POOL MODEL (use_model=no): EVERY POPULATED CLASS IS USED'
        endif
        if( .not. l_kept )then
            ! map_cavgs_selection only deselects: the particles' states from before it first
            call proj%os_ptcl2D%set_all('state', states2D)
            call proj%os_ptcl3D%set_all('state', states3D)
            cls_states = merge(1, 0, proj%os_cls2D%get_all_asint('pop') > 0)
            if( .not. any(cls_states > 0) ) cls_states = 1
            call proj%map_cavgs_selection(cls_states)
        endif
        write(logfhandle,'(A,I6,A,I6,A,I8,A)') '>>> THE ROUTE KEEPS ', count(proj%os_cls2D%get_all_asint('state') > 0),&
            &' OF ', proj%os_cls2D%get_noris(), ' CLASSES AND ', proj%os_ptcl3D%count_state_gt_zero(), ' PARTICLES'
    end subroutine select_route_classes

    !> The selected rows of @p proj whose class the pool model rejects, @p cls_states its verdict
    !! per class (score_cavgs_with_model): logged beside the fit's rejection (map_assignment).
    function model_rejected_rows( proj, cls_states ) result( l_rej )
        class(sp_project), intent(inout) :: proj
        integer,           intent(in)    :: cls_states(:)
        logical, allocatable :: l_rej(:)
        integer :: n, i, icls
        n = proj%os_ptcl2D%get_noris()
        allocate(l_rej(n), source=.false.)
        do i = 1,n
            if( proj%os_ptcl3D%get_state(i) <= 0 ) cycle
            icls = proj%os_ptcl2D%get_class(i)
            if( icls < 1 .or. icls > size(cls_states) ) cycle
            l_rej(i) = cls_states(icls) <= 0
        enddo
    end function model_rejected_rows

    !---------------- the reference volumes ----------------

    !> The volumes of solve3D_cavgs' result @p res checked (@p chk), under @p label in the log: the
    !! states that hold class averages (cls3D, at most @p nstates_max) renumbered from 1
    !! (compact_states); their volumes (vol_cavg entries) and the JPEGs beside them copied into
    !! @p dir with a record (REFS_RECORD: the count, box and pixel size; the angular distance
    !! checked, the state agreement, the largest correlation between two docked states and the mean
    !! angular distance over every class; the file names), and the class averages each state was
    !! made from (write_state_cavgs). The states are refined independently, so each is in a frame
    !! and hand of its own: with @p l_dock every pair is compared after docking in both hands
    !! (docked_pair_cc, the mask diameter @p mskdiam fitted to the volumes' box); without, none is
    !! (a test without the docking's time). chk%l_vols is .false., with nothing written, when no
    !! state holds class averages or one has no volume; chk%l_pass when the check passes
    !! (refs_usable).
    subroutine check_state_refs( res, nstates_max, mskdiam, dir, label, l_dock, chk )
        class(sp_project), intent(inout) :: res
        integer,           intent(in)    :: nstates_max
        real,              intent(in)    :: mskdiam
        class(string),     intent(in)    :: dir
        character(len=*),  intent(in)    :: label
        logical,           intent(in)    :: l_dock
        type(refs_check),  intent(out)   :: chk
        integer,      allocatable :: cls_states(:), newind(:)
        type(string), allocatable :: names(:)
        type(string) :: volpath, src
        real         :: smpd, cc, cc_direct, cc_mirror, msk
        integer      :: s, s2, i, box, ncopied, funit, ios
        if( res%os_cls3D%get_noris() > 0 ) cls_states = res%os_cls3D%get_all_asint('state')
        if( .not. allocated(cls_states) ) allocate(cls_states(0))
        call compact_states(cls_states, newind, chk%nheld)
        if( chk%nheld < 1 .or. chk%nheld > nstates_max )then
            write(logfhandle,'(A,A,A)') '>>> WARNING: ', label, ' LEFT NO CLASS AVERAGE WITH A STATE'
            return
        endif
        allocate(chk%counts(chk%nheld), source=0)
        do i = 1,size(cls_states)
            if( cls_states(i) > 0 ) chk%counts(cls_states(i)) = chk%counts(cls_states(i)) + 1
        enddo
        ! the even/odd agreement solve3D_cavgs records with its volumes: the same values on every
        ! vol_cavg entry (solve3D_cavgs_policy.md, section 7), so the first is read. The angular
        ! distance checked is the median over the class averages whose halves are in one state
        ! (eo_angdist_median); the mean over every class (eo_angdist), which compares poses across
        ! states' frames too, is logged, and checked only when the median is missing
        do i = 1,res%os_out%get_noris()
            if( .not. res%os_out%isthere(i, 'imgkind') ) cycle
            if( res%os_out%get_str(i, 'imgkind') /= 'vol_cavg' ) cycle
            if( res%os_out%isthere(i, 'eo_angdist') )        chk%angdist_all = res%os_out%get(i, 'eo_angdist')
            if( res%os_out%isthere(i, 'eo_angdist_median') ) chk%angdist     = res%os_out%get(i, 'eo_angdist_median')
            if( res%os_out%isthere(i, 'eo_state_overlap') )  chk%overlap     = res%os_out%get(i, 'eo_state_overlap')
            exit
        enddo
        if( chk%angdist < 0. ) chk%angdist = chk%angdist_all
        ! the volumes, renumbered
        allocate(names(chk%nheld))
        call simple_mkdir(dir)
        ncopied = 0
        do s = 1,size(newind)
            if( newind(s) < 1 ) cycle
            if( .not. res%isthere_in_osout('vol_cavg', s) ) cycle
            call res%get_vol('vol_cavg', s, volpath, smpd, box)
            if( .not. file_exists(volpath) ) cycle
            names(newind(s)) = 'vol_state'//int2str_pad(newind(s), 2)//MRC_EXT
            call simple_copy_file(volpath, dir//'/'//names(newind(s))%to_char())
            src = get_fpath(volpath)//refine3D_reprojs_fname(s)
            if( file_exists(src) ) call simple_copy_file(src, dir//'/'//refine3D_reprojs_fname(newind(s)))
            src = get_fpath(volpath)//refine3D_oris_heatmap_fname(s)
            if( file_exists(src) ) call simple_copy_file(src, dir//'/'//refine3D_oris_heatmap_fname(newind(s)))
            ncopied  = ncopied + 1
            chk%box  = box
            chk%smpd = smpd
        enddo
        if( ncopied /= chk%nheld )then
            write(logfhandle,'(A,A,A,I3,A,I3,A)') '>>> WARNING: ', label, ' HAS VOLUMES FOR ', ncopied, ' OF THE ', chk%nheld,&
                &' STATES HOLDING CLASS AVERAGES'
            return
        endif
        chk%l_vols = .true.
        ! every pair of states compared after docking: near-duplicates are no states
        allocate(chk%cc_pairs(chk%nheld,chk%nheld), source=-1.)
        if( l_dock )then
            msk = fit_mskdiam(mskdiam, chk%box, chk%smpd)
            do s = 1,chk%nheld - 1
                do s2 = s + 1,chk%nheld
                    cc = docked_pair_cc(dir//'/'//names(s)%to_char(), dir//'/'//names(s2)%to_char(), chk%smpd, msk,&
                        &cc_direct, cc_mirror)
                    write(logfhandle,'(A,I3,A,I3,A,F7.4,A,F7.4,A)') '>>>   STATES ', s, ' AND ', s2, ': CORRELATION AFTER DOCKING ',&
                        &cc_direct, ' (', cc_mirror, ' MIRRORED)'
                    chk%cc_pairs(s,s2) = cc
                    chk%cc_pairs(s2,s) = cc
                    chk%cc_max = max(chk%cc_max, cc)
                enddo
            enddo
        else
            write(logfhandle,'(A)') '>>> THE STATES ARE NOT COMPARED (NO DOCKING)'
        endif
        ! the class averages each state was made from, a stack and a JPEG per state
        call write_check_cavgs()
        open(newunit=funit, file=dir%to_char()//'/'//REFS_RECORD, status='replace', action='write', iostat=ios)
        if( ios /= 0 ) THROW_HARD('cannot write the record of '//dir%to_char())
        write(funit,'(I0,1X,I0,1X,F12.6)') chk%nheld, chk%box, chk%smpd
        write(funit,'(F12.4,1X,F12.4,1X,F12.4,1X,F12.4)') chk%angdist, chk%overlap, chk%cc_max, chk%angdist_all
        do s = 1,chk%nheld
            write(funit,'(A)') names(s)%to_char()
        enddo
        close(funit)
        chk%l_pass = refs_usable(chk%angdist, chk%overlap, chk%counts, chk%cc_max)
        write(logfhandle,'(A,I3,A,F7.2,A,F7.2,A,F7.2,A,F7.4)') '>>> '//label//': ', chk%nheld,&
            &' STATES HOLD CLASS AVERAGES; EVEN/ODD ANGULAR DISTANCE ', chk%angdist, ' DEGREES (MEDIAN, HALVES IN ONE STATE; ',&
            &chk%angdist_all, ' MEAN OVER ALL), STATE AGREEMENT ', chk%overlap, '%, LARGEST CORRELATION BETWEEN TWO DOCKED STATES ',&
            &chk%cc_max
        do s = 1,chk%nheld
            write(logfhandle,'(A,I3,A,I6,A)') '>>>   STATE ', s, ': ', chk%counts(s), ' CLASS AVERAGES'
        enddo
        if( chk%l_pass )then
            write(logfhandle,'(A)') '>>> THE VOLUMES PASS THE CHECK'
        else
            write(logfhandle,'(A)') '>>> THE VOLUMES FAIL THE CHECK'
        endif

    contains

        ! the class averages solve3D_cavgs aligned (its project's cavg entry, one per cls3D entry),
        ! grouped by their renumbered state into dir: none when they are not on disk
        subroutine write_check_cavgs()
            type(image), allocatable :: cavgs(:)
            type(string) :: stk
            real         :: smpd_stk
            integer      :: ncls_stk
            ncls_stk = 0
            call res%get_cavgs_stk(stk, ncls_stk, smpd_stk, fail=.false.)
            if( ncls_stk < 1 .or. ncls_stk /= size(cls_states) .or. .not. file_exists(stk) )then
                write(logfhandle,'(A)') '>>> NO CLASS AVERAGES OF '//label//' TO GROUP BY STATE'
                return
            endif
            cavgs = read_cavgs_into_imgarr(res)
            call write_state_cavgs(cavgs, cls_states, chk%nheld, dir, label)
            call dealloc_imgarr(cavgs)
        end subroutine write_check_cavgs

    end subroutine check_state_refs

    !> The correlation of the volumes @p vol1 and @p vol2 (pixel size @p smpd, mask diameter
    !! @p mskdiam in A) once @p vol2 is docked onto @p vol1 over rotations and shifts, from
    !! DOCK_HP_REFS to the matching's limit (LP_ASSIGN): the better of @p vol2 as it is
    !! (@p cc_direct) and mirrored (@p cc_mirror), since independently refined states may differ in
    !! hand. A correlation that is not finite counts as identical (1), with a warning.
    real function docked_pair_cc( vol1, vol2, smpd, mskdiam, cc_direct, cc_mirror ) result( cc )
        class(string), intent(in)  :: vol1, vol2
        real,          intent(in)  :: smpd, mskdiam
        real,          intent(out) :: cc_direct, cc_mirror
        type(dock_vols) :: docker
        type(image)     :: vol
        type(string)    :: mirror
        real    :: eul(3), shift(3)
        integer :: ldim(3), nimgs
        ! vol2 mirrored, beside it
        mirror = add2fbody(vol2, MRC_EXT, '_mirror')
        call find_ldim_nptcls(vol2, ldim, nimgs)
        call vol%new(ldim, smpd)
        call vol%read(vol2)
        call vol%mirror('x')
        call vol%write(mirror)
        call vol%kill
        call docker%new(vol1, vol2, smpd, DOCK_HP_REFS, LP_ASSIGN, mskdiam)
        call docker%srch()
        call docker%get_dock_info(eul, shift, cc_direct)
        call docker%kill()
        call docker%new(vol1, mirror, smpd, DOCK_HP_REFS, LP_ASSIGN, mskdiam)
        call docker%srch()
        call docker%get_dock_info(eul, shift, cc_mirror)
        call docker%kill()
        call del_file(mirror)
        if( .not. ieee_is_finite(cc_direct) .or. .not. ieee_is_finite(cc_mirror) )then
            THROW_WARN('a docked correlation is not finite; the states count as identical')
            if( .not. ieee_is_finite(cc_direct) ) cc_direct = 1.
            if( .not. ieee_is_finite(cc_mirror) ) cc_mirror = 1.
        endif
        cc = max(cc_direct, cc_mirror)
    end function docked_pair_cc

    !> The docked correlation of the volumes @p vol1 and @p vol2 (docked_pair_cc) on a common grid:
    !! when their boxes differ, the larger is Fourier-cropped to the smaller's box (written beside it
    !! as <fbody>_common, removed after), which keeps its extent, so the pixel sizes must then agree.
    !! The mask diameter @p mskdiam (A) is fitted to the common box. -1, with a warning, when the
    !! volumes cover different extents.
    real function docked_cc_common_grid( vol1, vol2, mskdiam, cc_direct, cc_mirror ) result( cc )
        class(string), intent(in)  :: vol1, vol2
        real,          intent(in)  :: mskdiam
        real,          intent(out) :: cc_direct, cc_mirror
        type(string) :: v1, v2, tmp
        real    :: smpd1, smpd2
        integer :: ldim1(3), ldim2(3), n
        cc        = -1.
        cc_direct = -1.
        cc_mirror = -1.
        call find_ldim_nptcls(vol1, ldim1, n)
        call find_ldim_nptcls(vol2, ldim2, n)
        smpd1 = find_img_smpd(vol1)
        smpd2 = find_img_smpd(vol2)
        v1    = vol1
        v2    = vol2
        if( ldim1(1) > ldim2(1) )then
            tmp   = add2fbody(vol1, MRC_EXT, '_common')
            call fourier_crop(vol1, ldim1, smpd1, ldim2(1), tmp)
            smpd1 = smpd1 * real(ldim1(1)) / real(ldim2(1))
            v1    = tmp
        else if( ldim2(1) > ldim1(1) )then
            tmp   = add2fbody(vol2, MRC_EXT, '_common')
            call fourier_crop(vol2, ldim2, smpd2, ldim1(1), tmp)
            smpd2 = smpd2 * real(ldim2(1)) / real(ldim1(1))
            v2    = tmp
        endif
        if( abs(smpd1 - smpd2) <= 1.e-3 * smpd2 )then
            cc = docked_pair_cc(v1, v2, smpd2, fit_mskdiam(mskdiam, min(ldim1(1), ldim2(1)), smpd2), cc_direct, cc_mirror)
        else
            THROW_WARN('the volumes cover different extents; they are not docked')
        endif
        if( tmp%strlen() > 0 ) call del_file(tmp)

    contains

        ! the volume @p fname (@p ldim, @p smpd) Fourier-cropped to @p box, written to @p fname_out
        subroutine fourier_crop( fname, ldim, smpd, box, fname_out )
            class(string), intent(in) :: fname, fname_out
            integer,       intent(in) :: ldim(3), box
            real,          intent(in) :: smpd
            type(image) :: vol
            call vol%new(ldim, smpd)
            call vol%read(fname)
            call vol%fft()
            call vol%clip_inplace([box, box, box])
            call vol%ifft()
            call vol%set_smpd(smpd * real(ldim(1)) / real(box))
            call vol%write(fname_out)
            call vol%kill
        end subroutine fourier_crop

    end function docked_cc_common_grid

    !> The reference volumes recorded in @p dir (REFS_RECORD, check_state_refs): @p vols their
    !! absolute paths, @p box and @p smpd. .false., with none taken, when the record is missing or
    !! unreadable, holds more than @p nstates_max, or a volume is not on disk.
    logical function read_refs_record( dir, nstates_max, vols, box, smpd ) result( l_ok )
        class(string),             intent(in)    :: dir
        integer,                   intent(in)    :: nstates_max
        type(string), allocatable, intent(inout) :: vols(:)
        integer,                   intent(out)   :: box
        real,                      intent(out)   :: smpd
        type(string), allocatable :: names(:)
        character(len=1024) :: line
        integer :: funit, ios, n, s
        l_ok = .false.
        box  = 0
        smpd = 0.
        if( .not. file_exists(dir//'/'//REFS_RECORD) ) return
        open(newunit=funit, file=dir%to_char()//'/'//REFS_RECORD, status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        read(funit,*,iostat=ios) n, box, smpd
        if( ios /= 0 .or. n < 1 .or. n > nstates_max )then
            close(funit)
            return
        endif
        read(funit,'(A)',iostat=ios) line ! the even/odd values
        allocate(names(n))
        do s = 1,n
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            names(s) = dir//'/'//trim(adjustl(line))
            if( .not. file_exists(names(s)) ) ios = 1
            if( ios /= 0 ) exit
            names(s) = simple_abspath(names(s))
        enddo
        close(funit)
        if( ios /= 0 ) return
        if( allocated(vols) ) deallocate(vols)
        allocate(vols(n))
        do s = 1,n
            vols(s) = names(s)
        enddo
        l_ok = .true.
    end function read_refs_record

    !---------------- the assignment ----------------

    !> The matching in @p res (a copy of @p spproj's rows, with the matched classes' cls2D, cls3D
    !! and the particles' class labels; assign_cavgs_states' project) onto @p spproj, against
    !! @p nrefs references. Over the classes with a state and particles, those whose best
    !! correlation is an outlier are rejected (reject_by_fit: state 0) and the kept ones with a
    !! small margin flagged (flag_ambiguous: cls2D key 'ambiguous'); each kept class's pose and
    !! state are mapped onto its selected particles (map2ptcls), and each row selected before
    !! takes its particle's state, pose and correlation (state 0 when its class is rejected or it
    !! has none), in both particle segments: a particle follows its (sub-)class. @p desel marks
    !! the rows deselected here. The fit's statistics, each rejected class's correlation and
    !! population, and how the correlation follows the population are logged (@p what names the
    !! classes), with the counts beside the pool model's verdict @p model_rej (per row) when given;
    !! @p stats returns them. .false., with @p spproj and @p desel unchanged, when @p res holds
    !! other rows or no matching.
    logical function map_assignment( spproj, res, nrefs, what, stats, desel, model_rej ) result( l_ok )
        class(sp_project),    intent(inout) :: spproj, res
        integer,              intent(in)    :: nrefs
        character(len=*),     intent(in)    :: what
        type(assign_stats),   intent(out)   :: stats
        logical, allocatable, intent(inout) :: desel(:)
        logical, optional,    intent(in)    :: model_rej(:)
        integer, allocatable :: cls_states(:)
        real,    allocatable :: corrs(:), margins(:), pops(:)
        logical, allocatable :: l_cand(:), l_rej(:), l_amb(:), l_sel(:), l_mapped(:)
        integer :: n, ncls, icls, i, s
        l_ok = .false.
        n    = spproj%os_ptcl3D%get_noris()
        ncls = res%os_cls2D%get_noris()
        if( n == 0 .or. res%os_ptcl2D%get_noris() /= n ) return
        if( ncls == 0 .or. res%os_cls3D%get_noris() /= ncls ) return
        cls_states = res%os_cls3D%get_all_asint('state')
        corrs      = res%os_cls3D%get_all('corr')
        margins    = res%os_cls3D%get_all('margin')
        ! the candidates: classes with a state and particles
        allocate(l_cand(ncls))
        do icls = 1,ncls
            l_cand(icls) = cls_states(icls) > 0 .and. cls_states(icls) <= nrefs
            if( res%os_cls2D%isthere(icls, 'pop') ) l_cand(icls) = l_cand(icls) .and. res%os_cls2D%get_int(icls, 'pop') > 0
        enddo
        if( .not. any(l_cand) ) return
        call reject_by_fit(corrs, l_cand, l_rej, stats%med, stats%mad, stats%thres, stats%nmax)
        l_amb = flag_ambiguous(margins, l_cand .and. .not. l_rej)
        ! the populations, for the log: whether the rejection follows the classes' sizes
        allocate(pops(ncls), source=0.)
        do icls = 1,ncls
            if( res%os_cls2D%isthere(icls, 'pop') ) pops(icls) = real(res%os_cls2D%get_int(icls, 'pop'))
        enddo
        do icls = 1,ncls
            if( .not. l_cand(icls) .or. l_rej(icls) ) call res%os_cls3D%set_state(icls, 0)
            call res%os_cls2D%set(icls, 'ambiguous', merge(1, 0, l_amb(icls)))
        enddo
        ! the particles each kept class maps onto: selected, with a class of the matching
        allocate(l_sel(n), l_mapped(n))
        do i = 1,n
            l_sel(i)    = spproj%os_ptcl3D%get_state(i) > 0
            icls        = res%os_ptcl2D%get_class(i)
            l_mapped(i) = res%os_ptcl2D%get_state(i) > 0 .and. icls >= 1 .and. icls <= ncls
        enddo
        call res%map2ptcls()
        ! map2ptcls maps nothing onto particles without a 2D alignment (a virgin ptcl2D)
        if( res%is_virgin_field('ptcl3D') ) return
        stats%ncls  = ncls
        stats%ncand = count(l_cand)
        stats%nrej  = count(l_rej)
        stats%namb  = count(l_amb)
        allocate(stats%counts(nrefs), source=0)
        if( allocated(desel) ) deallocate(desel)
        allocate(desel(n), source=.false.)
        stats%l_model = present(model_rej)
        do i = 1,n
            if( .not. l_sel(i) ) cycle
            s = 0
            if( l_mapped(i) ) s = res%os_ptcl3D%get_state(i)
            if( s > nrefs ) s = 0
            if( s > 0 )then
                call spproj%os_ptcl3D%set_euler(i, res%os_ptcl3D%get_euler(i))
                call spproj%os_ptcl3D%set_shift(i, res%os_ptcl3D%get_2Dshift(i))
                call spproj%os_ptcl3D%set(i, 'corr', res%os_ptcl3D%get(i, 'corr'))
                stats%counts(s) = stats%counts(s) + 1
            else
                stats%ndesel = stats%ndesel + 1
                desel(i)     = .true.
            endif
            call spproj%os_ptcl3D%set_state(i, s)
            call spproj%os_ptcl2D%set_state(i, s)
            ! the fit's rejection beside the pool model's
            if( stats%l_model )then
                if( i <= size(model_rej) )then
                    if( s == 0 .and. model_rej(i) ) stats%nboth  = stats%nboth  + 1
                    if( s >  0 .and. model_rej(i) ) stats%nmodel = stats%nmodel + 1
                    if( s == 0 .and. .not. model_rej(i) ) stats%nfit = stats%nfit + 1
                endif
            endif
        enddo
        stats%nassigned = sum(stats%counts)
        write(logfhandle,'(A,I6,A)') '>>> ASSIGNMENT OF ', ncls, ' '//what
        write(logfhandle,'(A,I6,A,I6,A,I6,A)') '>>> ', stats%ncand, ' CANDIDATES: ', stats%nrej, ' REJECTED BY FIT, ',&
            &stats%namb, ' AMBIGUOUS (ASSIGNED)'
        call log_fit()
        do s = 1,nrefs
            write(logfhandle,'(A,I3,A,I8,A)') '>>>   STATE ', s, ': ', stats%counts(s), ' PARTICLES'
        enddo
        write(logfhandle,'(A,I8,A)') '>>> ', stats%ndesel, ' PARTICLES DESELECTED'
        if( stats%l_model )then
            write(logfhandle,'(A,I8,A,I8,A,I8)') '>>> PARTICLES REJECTED BY THE FIT AND THE POOL MODEL: ', stats%nboth,&
                &'; THE MODEL ALONE (NOT APPLIED): ', stats%nmodel, '; THE FIT ALONE: ', stats%nfit
        endif
        l_ok = .true.

    contains

        ! the fit's statistics; each rejected class's correlation, population and state; the
        ! median population of the kept and rejected candidates; and the correlation of the best
        ! correlation with the log of the population over the candidates (a strong positive one:
        ! the rejection follows the classes' sizes more than their fit)
        subroutine log_fit()
            real, allocatable :: x(:), y(:)
            real    :: r, sx, sy
            integer :: j
            write(logfhandle,'(A,F7.4,A,F7.4,A,F7.4,A,I6,A,I6,A)') '>>> FIT: MEDIAN BEST CORRELATION ', stats%med, ', MAD ',&
                &stats%mad, ', THRESHOLD ', stats%thres, ' (', count(l_cand .and. corrs < stats%thres),&
                &' CANDIDATES BELOW IT, AT MOST ', stats%nmax, ' REJECTED)'
            do j = 1,size(corrs)
                if( .not. l_rej(j) ) cycle
                write(logfhandle,'(A,I6,A,F7.4,A,I8,A,I3)') '>>>   REJECTED CLASS ', j, ': CORRELATION ', corrs(j),&
                    &', POPULATION ', nint(pops(j)), ', STATE ', cls_states(j)
            enddo
            if( .not. any(pops > 0.) ) return
            if( any(l_cand .and. .not. l_rej) .and. any(l_rej) )then
                write(logfhandle,'(A,F10.1,A,F10.1)') '>>> MEDIAN POPULATION OF THE KEPT / REJECTED CANDIDATES: ',&
                    &median_of(pack(pops, l_cand .and. .not. l_rej)), ' / ', median_of(pack(pops, l_rej))
            endif
            if( count(l_cand) < 3 ) return
            x  = pack(corrs, l_cand)
            y  = log(max(pack(pops, l_cand), 1.))
            x  = x - sum(x) / real(size(x))
            y  = y - sum(y) / real(size(y))
            sx = sqrt(sum(x**2))
            sy = sqrt(sum(y**2))
            r  = 0.
            if( sx > 0. .and. sy > 0. ) r = sum(x * y) / (sx * sy)
            write(logfhandle,'(A,F7.4)') '>>> CORRELATION OF THE BEST CORRELATION WITH LOG POPULATION OVER THE CANDIDATES: ', r
        end subroutine log_fit

    end function map_assignment

    !> The matched (sub-)class averages of @p res (map_assignment done) grouped into @p dir under
    !! @p label: a stack and a JPEG per state of the @p nrefs (the ambiguous ones included), one of
    !! those the fit rejected (a state given by the matching, @p states_matched, and none after the
    !! rejection, with particles), and one of the ambiguous ones (write_state_cavgs). None when the
    !! class averages are not on disk.
    subroutine write_matched_cavgs( res, states_matched, nrefs, dir, label )
        class(sp_project), intent(inout) :: res
        integer,           intent(in)    :: states_matched(:), nrefs
        class(string),     intent(in)    :: dir
        character(len=*),  intent(in)    :: label
        type(image), allocatable :: cavgs(:)
        integer,     allocatable :: states(:)
        logical,     allocatable :: l_rej(:), l_amb(:)
        type(string) :: stk
        real         :: smpd_stk
        integer      :: ncls, ncls_stk, icls
        ncls     = res%os_cls3D%get_noris()
        ncls_stk = 0
        call res%get_cavgs_stk(stk, ncls_stk, smpd_stk, fail=.false.)
        if( ncls < 1 .or. ncls_stk /= ncls .or. size(states_matched) /= ncls .or. .not. file_exists(stk) )then
            write(logfhandle,'(A)') '>>> NO CLASS AVERAGES OF THE '//label//' TO GROUP BY STATE'
            return
        endif
        states = res%os_cls3D%get_all_asint('state')
        allocate(l_rej(ncls), l_amb(ncls), source=.false.)
        do icls = 1,ncls
            l_rej(icls) = states_matched(icls) > 0 .and. states(icls) == 0
            if( res%os_cls2D%isthere(icls, 'pop') ) l_rej(icls) = l_rej(icls) .and. res%os_cls2D%get_int(icls, 'pop') > 0
            if( res%os_cls2D%isthere(icls, 'ambiguous') ) l_amb(icls) = res%os_cls2D%get_int(icls, 'ambiguous') == 1
        enddo
        cavgs = read_cavgs_into_imgarr(res)
        call write_state_cavgs(cavgs, states, nrefs, dir, label, l_rej=l_rej, l_amb=l_amb)
        call dealloc_imgarr(cavgs)
    end subroutine write_matched_cavgs

    !> The classes of @p proj (the assignment's project, map_assignment done) that state @p s holds
    !! (cls3D), selected for a solve3D_cavgs of that state alone: cls2D state 1 for them and 0 for
    !! every other class (the halves' own flags kept). The number selected; 0 when @p proj has no
    !! matching.
    integer function select_state_classes( proj, s ) result( nsel )
        class(sp_project), intent(inout) :: proj
        integer,           intent(in)    :: s
        integer :: icls, ncls
        nsel = 0
        ncls = proj%os_cls2D%get_noris()
        if( ncls == 0 .or. proj%os_cls3D%get_noris() /= ncls ) return
        do icls = 1,ncls
            if( proj%os_cls3D%get_state(icls) == s )then
                call proj%os_cls2D%set_state(icls, 1)
                nsel = nsel + 1
            else
                call proj%os_cls2D%set_state(icls, 0)
            endif
        enddo
    end function select_state_classes

    !> The poses of state @p s's particles from the solve3D_cavgs of that state's classes alone
    !! (@p res, its project), in place of the matching's: res's class poses mapped onto its
    !! particles (map2ptcls), and each row of @p spproj in state @p s whose class res holds in a
    !! state given that Euler angles and shift. The state's frame is that run's own: refine3D
    !! compares each particle with its own state's volume only. The number of rows given a pose; -1,
    !! with @p spproj unchanged, when @p res holds other rows or no class poses.
    integer function map_state_poses( spproj, res, s ) result( nposes )
        class(sp_project), intent(inout) :: spproj, res
        integer,           intent(in)    :: s
        integer :: i, n
        nposes = -1
        n      = spproj%os_ptcl3D%get_noris()
        if( n == 0 .or. res%os_ptcl2D%get_noris() /= n ) return
        if( res%os_cls3D%get_noris() == 0 ) return
        call res%map2ptcls()
        if( res%os_ptcl3D%get_noris() /= n ) return
        nposes = 0
        do i = 1,n
            if( spproj%os_ptcl3D%get_state(i) /= s ) cycle
            if( res%os_ptcl3D%get_state(i) <= 0 ) cycle
            call spproj%os_ptcl3D%set_euler(i, res%os_ptcl3D%get_euler(i))
            call spproj%os_ptcl3D%set_shift(i, res%os_ptcl3D%get_2Dshift(i))
            nposes = nposes + 1
        enddo
    end function map_state_poses

    !---------------- the local refinement ----------------

    !> The seeds of the local refinement on @p spproj's rows (refine=cont continues from a discrete
    !! pose: state, half-set, projection direction): the half-sets are the project's partition
    !! (eo_transaction: one half per micrograph, recorded on the stack rows and never changed, so
    !! every refinement has the same halves and a new stack is labelled once), and the projection
    !! direction closest to each mapped pose in the job's grid (@p nspace directions, c1).
    !! refine3D clears the update counts at its first iteration, so its starting reconstruction
    !! takes every row with a state.
    subroutine set_refine_seeds( spproj, nspace )
        class(sp_project), intent(inout) :: spproj
        integer,           intent(in)    :: nspace
        type(oris) :: eulspace
        type(sym)  :: c1
        logical    :: l_labelled
        call eo_transaction(spproj, l_labelled)
        call c1%new('c1')
        call eulspace%new(nspace, is_ptcl=.false.)
        call c1%build_refspiral(eulspace)
        call spproj%os_ptcl3D%set_projs(eulspace)
        call eulspace%kill
        call c1%kill
    end subroutine set_refine_seeds

    !> The local refinement's project: @p spproj written to @p projfile with its particles' states
    !! renumbered from 1 (compact_states), since refine3D reconstructs no volume for an empty state;
    !! @p spproj itself unchanged (map_refine_result maps the job's states back). @p nheld: the
    !! states the particles hold; nothing is written when none.
    subroutine write_refine_project( spproj, projfile, nheld )
        class(sp_project), intent(inout) :: spproj
        class(string),     intent(in)    :: projfile
        integer,           intent(out)   :: nheld
        integer, allocatable :: states(:), states_job(:), newind(:)
        states     = spproj%os_ptcl3D%get_all_asint('state')
        states_job = states
        call compact_states(states_job, newind, nheld)
        if( nheld < 1 ) return
        call spproj%os_ptcl3D%set_all('state', states_job)
        call spproj%write(projfile)
        call spproj%os_ptcl3D%set_all('state', states)
    end subroutine write_refine_project

    !> The local refinement's project @p res onto @p spproj, its states mapped back to
    !! @p spproj's (write_refine_project renumbered them): each state's volume and FSC it
    !! registered, in place of every earlier state's (of the @p nstates_max), and each selected
    !! row's refined pose (Euler angles and shift) where its state is the row's. The number of
    !! states registered; 0, with @p spproj unchanged, when @p res holds other rows or no volume.
    integer function map_refine_result( spproj, res, nstates_max ) result( nreg )
        class(sp_project), intent(inout) :: spproj, res
        integer,           intent(in)    :: nstates_max
        integer, allocatable :: states_job(:), newind(:)
        logical, allocatable :: l_vol(:)
        type(string) :: vol, fsc
        real         :: smpd
        integer      :: n, s, i, box, fsc_box, nposes, nheld
        nreg = 0
        n    = spproj%os_ptcl3D%get_noris()
        if( n == 0 .or. res%os_ptcl3D%get_noris() /= n ) return
        ! the job's state of each state of spproj (0: none)
        states_job = spproj%os_ptcl3D%get_all_asint('state')
        call compact_states(states_job, newind, nheld)
        allocate(l_vol(size(newind)), source=.false.)
        do s = 1,size(newind)
            if( newind(s) < 1 ) cycle
            if( .not. res%isthere_in_osout('vol', newind(s)) .or. .not. res%isthere_in_osout('fsc', newind(s)) ) cycle
            call res%get_vol('vol', newind(s), vol, smpd, box)
            call res%get_fsc(newind(s), fsc, fsc_box)
            l_vol(s) = file_exists(vol) .and. file_exists(fsc)
        enddo
        if( .not. any(l_vol) ) return
        do s = 1,nstates_max
            call spproj%remove_state_artifacts_from_osout(s)
        enddo
        do s = 1,size(l_vol)
            if( .not. l_vol(s) ) cycle
            call res%get_vol('vol', newind(s), vol, smpd, box)
            call res%get_fsc(newind(s), fsc, fsc_box)
            call spproj%add_vol2os_out(vol, smpd, s, 'vol')
            call spproj%add_fsc2os_out(fsc, s, fsc_box)
            nreg = nreg + 1
        enddo
        nposes = 0
        do i = 1,n
            s = spproj%os_ptcl3D%get_state(i)
            if( s <= 0 .or. s > size(newind) ) cycle
            if( res%os_ptcl3D%get_state(i) /= newind(s) ) cycle
            call spproj%os_ptcl3D%set_euler(i, res%os_ptcl3D%get_euler(i))
            call spproj%os_ptcl3D%set_shift(i, res%os_ptcl3D%get_2Dshift(i))
            nposes = nposes + 1
        enddo
        write(logfhandle,'(A,I3,A,I8,A)') '>>> LOCAL REFINEMENT: ', nreg, ' STATE VOLUMES, ', nposes, ' REFINED POSES'
    end function map_refine_result

    !---------------- class averages by state ----------------

    !> The class averages @p cavgs grouped by @p states (1..@p nstates) into @p dir: per state a
    !! stack and a tiled JPEG (cavgs_fbody), and, when given, one of the classes @p l_rej marks
    !! (cavgs_rejected) and one of those @p l_amb marks (cavgs_ambiguous); each JPEG logged
    !! (log_jpeg) under @p label. A group without class averages writes nothing and is logged.
    subroutine write_state_cavgs( cavgs, states, nstates, dir, label, l_rej, l_amb )
        class(image),      intent(inout) :: cavgs(:)
        integer,           intent(in)    :: states(:), nstates
        class(string),     intent(in)    :: dir
        character(len=*),  intent(in)    :: label
        logical, optional, intent(in)    :: l_rej(:), l_amb(:)
        integer :: s
        if( size(states) /= size(cavgs) ) THROW_HARD('# states /= # class averages; write_state_cavgs')
        do s = 1,nstates
            call write_cavgs_jpeg(states == s, dir//'/'//cavgs_fbody(s), label//' STATE '//int2str_pad(s, 2))
        enddo
        if( present(l_rej) ) call write_cavgs_jpeg(l_rej, dir//'/cavgs_rejected',  label//' REJECTED')
        if( present(l_amb) ) call write_cavgs_jpeg(l_amb, dir//'/cavgs_ambiguous', label//' AMBIGUOUS')

    contains

        subroutine write_cavgs_jpeg( mask, fbody, what )
            logical,          intent(in) :: mask(:)
            class(string),    intent(in) :: fbody
            character(len=*), intent(in) :: what
            type(string) :: stk, jpg
            integer      :: icls, k
            stk = fbody//MRC_EXT
            jpg = fbody//JPG_EXT
            if( file_exists(stk) ) call del_file(stk)
            if( file_exists(jpg) ) call del_file(jpg)
            if( size(mask) /= size(cavgs) .or. .not. any(mask) )then
                write(logfhandle,'(A)') '>>> '//what//' CLASS AVERAGES: NONE'
                return
            endif
            k = 0
            do icls = 1,size(cavgs)
                if( .not. mask(icls) ) cycle
                k = k + 1
                call cavgs(icls)%write(stk, k)
            enddo
            call mrc2jpeg_tiled(stk, jpg)
            write(logfhandle,'(A,I6)') '>>> '//what//' CLASS AVERAGES: ', k
            call log_jpeg(what//' CLASS AVERAGES', jpg)
        end subroutine write_cavgs_jpeg

    end subroutine write_state_cavgs

    !> The file body of state @p s's stack of class averages (cavgs_state<NN>).
    function cavgs_fbody( s ) result( fbody )
        integer, intent(in) :: s
        character(len=:), allocatable :: fbody
        fbody = 'cavgs_state'//int2str_pad(s, 2)
    end function cavgs_fbody

    !> The JPEG @p jpg logged under @p label with its absolute path, then once more as a line of its
    !! own ('>>> JPEG <path>'); nothing when it is not on disk.
    subroutine log_jpeg( label, jpg )
        character(len=*), intent(in) :: label
        class(string),    intent(in) :: jpg
        type(string) :: path
        if( .not. file_exists(jpg) ) return
        path = simple_abspath(jpg, check_exists=.false.)
        write(logfhandle,'(A,A)') '>>> '//label//' JPEG ', path%to_char()
        write(logfhandle,'(A,A)') '>>> JPEG ', path%to_char()
    end subroutine log_jpeg

    !---------------- the jobs' command lines ----------------
    ! Each sets the program and the keys that make the method on a fresh @p cline, in the job's
    ! folder (mkdir=no, @p projfile beside it); the caller adds the threads and parts. A mask
    ! diameter @p mskdiam not positive is left to the job's default.

    !> The route's solve2D: the first set classified into NCLS2D_FIRST classes, cache=yes.
    subroutine set_solve2D_route_cline( cline, projfile, mskdiam )
        class(cmdline),   intent(inout) :: cline
        character(len=*), intent(in)    :: projfile
        real,             intent(in)    :: mskdiam
        call cline%kill
        call cline%set('prg',      'solve2D')
        call cline%set('mkdir',    'no')
        call cline%set('projfile', projfile)
        call cline%set('ncls',     NCLS2D_FIRST)
        if( mskdiam > 0. ) call cline%set('mskdiam', mskdiam)
        call cline%set('cache',    'yes')
    end subroutine set_solve2D_route_cline

    !> cls_expansion: the classes of the selected particles split into @p ncls sub-classes in all
    !! (expansion_budget), their averages kernel pre-images of members (its defaults).
    subroutine set_clsexp_cline( cline, projfile, ncls, mskdiam )
        class(cmdline),   intent(inout) :: cline
        character(len=*), intent(in)    :: projfile
        integer,          intent(in)    :: ncls
        real,             intent(in)    :: mskdiam
        call cline%kill
        call cline%set('prg',      'cls_expansion')
        call cline%set('mkdir',    'no')
        call cline%set('projfile', projfile)
        call cline%set('ncls',     ncls)
        if( mskdiam > 0. ) call cline%set('mskdiam', mskdiam)
    end subroutine set_clsexp_cline

    !> solve3D_cavgs on the project's class averages and their halves: @p nstates volumes, c1, no
    !! pruning, the fixed low-pass ramp LPSTART_CAVGS3D to LPSTOP_CAVGS3D (force_lp_range, not one
    !! derived from the class FRCs).
    subroutine set_cavgs3D_cline( cline, projfile, nstates, mskdiam )
        class(cmdline),   intent(inout) :: cline
        character(len=*), intent(in)    :: projfile
        integer,          intent(in)    :: nstates
        real,             intent(in)    :: mskdiam
        call cline%kill
        call cline%set('prg',               'solve3D_cavgs')
        call cline%set('mkdir',             'no')
        call cline%set('projfile',          projfile)
        call cline%set('pgrp',              'c1')
        call cline%set('nstates',           nstates)
        call cline%set('prune',             'no')
        call cline%set('lpstart',           LPSTART_CAVGS3D)
        call cline%set('lpstop',            LPSTOP_CAVGS3D)
        call cline%set('force_lp_range',    'yes')
        if( mskdiam > 0. ) call cline%set('mskdiam', mskdiam)
    end subroutine set_cavgs3D_cline

    !> solve3D_cavgs of one state's class averages (a project selecting them, select_state_classes):
    !! one state, c1, no pruning, and the low-pass limits solve3D_cavgs derives from the class FRCs
    !! (its defaults, no fixed ramp), so each state has its own.
    subroutine set_state_cavgs3D_cline( cline, projfile, mskdiam )
        class(cmdline),   intent(inout) :: cline
        character(len=*), intent(in)    :: projfile
        real,             intent(in)    :: mskdiam
        call cline%kill
        call cline%set('prg',      'solve3D_cavgs')
        call cline%set('mkdir',    'no')
        call cline%set('projfile', projfile)
        call cline%set('pgrp',     'c1')
        call cline%set('nstates',  1)
        call cline%set('prune',    'no')
        if( mskdiam > 0. ) call cline%set('mskdiam', mskdiam)
    end subroutine set_state_cavgs3D_cline

    !> assign_cavgs_states: every class average of the project matched against the references
    !! @p vols (one per state) at LP_ASSIGN over NSPACE_ASSIGN directions per state.
    subroutine set_assign_cline( cline, projfile, vols, mskdiam )
        class(cmdline),   intent(inout) :: cline
        character(len=*), intent(in)    :: projfile
        class(string),    intent(in)    :: vols(:)
        real,             intent(in)    :: mskdiam
        integer :: s
        call cline%kill
        call cline%set('prg',      ASSIGN_LABEL)
        call cline%set('mkdir',    'no')
        call cline%set('projfile', projfile)
        call cline%set('nstates',  size(vols))
        do s = 1,size(vols)
            call cline%set('vol'//int2str(s), vols(s))
        enddo
        if( mskdiam > 0. ) call cline%set('mskdiam', mskdiam)
        call cline%set('lp',       LP_ASSIGN)
        call cline%set('nspace',   NSPACE_ASSIGN)
    end subroutine set_assign_cline

    !> The local refinement: one iteration (maxits=1) of refine3D refine=cont from the seeds
    !! (set_refine_seeds), @p nstates fixed states, @p nspace directions (the seeds' grid), c1, no
    !! sigma2 (objfun=cc). The caller runs it distributed (nparts, at 1 too): the distributed
    !! refine3D reconstructs its starting volumes, the shared-memory one needs them given.
    subroutine set_refine_cont_cline( cline, projfile, nstates, nspace, mskdiam )
        class(cmdline),   intent(inout) :: cline
        character(len=*), intent(in)    :: projfile
        integer,          intent(in)    :: nstates, nspace
        real,             intent(in)    :: mskdiam
        call cline%kill
        call cline%set('prg',      'refine3D')
        call cline%set('mkdir',    'no')
        call cline%set('projfile', projfile)
        call cline%set('refine',   'cont')
        call cline%set('maxits',   1)
        call cline%set('pgrp',     'c1')
        call cline%set('nstates',  nstates)
        call cline%set('nspace',   nspace)
        call cline%set('objfun',   'cc')
        if( mskdiam > 0. ) call cline%set('mskdiam', mskdiam)
    end subroutine set_refine_cont_cline

end module simple_state_assignment
