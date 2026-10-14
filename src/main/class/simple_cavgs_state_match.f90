!@descr: class averages matched against state volumes: each aligned with every reprojection of every state, its best state, 3D pose, correlation and margin
! Each state volume is reprojected over an even grid of directions on the whole sphere, so every
! mirrored view is a view of the grid and no mirror is searched. Every class average is aligned with
! every reprojection by a full in-plane search with shifts (match_imgs), at a common low-pass limit,
! after both sets are brought to one pixel size and box for that limit, and masked. A class's state is the one whose best
! reprojection correlates best; its margin is that correlation less the best of any other state;
! its 3D pose is that reprojection's direction with the in-plane angle and shift, in refine3D's
! convention (euler(3) = 360 - rotation), so map2ptcls can compose it with its particles' 2D poses.
module simple_cavgs_state_match
use simple_defs,             only: logfhandle, COSMSKHALFWIDTH
use simple_error,            only: simple_exception
use simple_math,             only: round2even
use simple_type_defs,        only: inpl_struct
use simple_image,            only: image
use simple_oris,             only: oris
use simple_cmdline,          only: cmdline
use simple_parameters,       only: parameters
use simple_simple_volinterp, only: reproject
use simple_strategy2D_utils, only: match_imgs
use simple_imgarr_utils,     only: dealloc_imgarr
implicit none

public :: match_cavgs_to_states, reduce_state_matches, matched_pose
private
#include "simple_local_flags.inc"

real,    parameter :: LP_SAMPLING   = 3.   ! the matched images' pixel size is lp/LP_SAMPLING (or the inputs' when coarser)
real,    parameter :: HP_MATCH      = 100. ! high-pass limit of the alignment (A)
integer, parameter :: BOX_MATCH_MIN = 32   ! the matched images' box at least (or the input's when smaller)

contains

    !> The class averages @p cavgs matched against the state volumes @p vols (one per state, one box
    !! and pixel size, not necessarily the class averages'): @p nspace reprojections per state over the whole sphere,
    !! aligned with every class average at @p lp (A) with shifts up to @p trs (pixels at the input
    !! sampling), within a soft mask of diameter @p mskdiam (A), on @p nthr threads. Per class: its
    !! @p states (the best), @p eulers (3 x ncls), @p shifts (2 x ncls, pixels at the input
    !! sampling), @p corrs (the best correlation) and @p margins (the best less the best of another
    !! state; 0 with one state).
    subroutine match_cavgs_to_states( cavgs, vols, nspace, lp, trs, mskdiam, nthr, states, eulers, shifts,&
            &corrs, margins )
        class(image),         intent(inout) :: cavgs(:), vols(:)
        integer,              intent(in)    :: nspace, nthr
        real,                 intent(in)    :: lp, trs, mskdiam
        integer, allocatable, intent(inout) :: states(:)
        real,    allocatable, intent(inout) :: eulers(:,:), shifts(:,:), corrs(:), margins(:)
        type(image),       allocatable :: cavgs_crop(:), vols_crop(:), refs(:), projs(:)
        type(inpl_struct), allocatable :: algninfo(:,:)
        real,              allocatable :: corrmat(:,:)
        logical,           allocatable :: l_msk(:,:,:)
        integer,           allocatable :: best_ref(:)
        type(parameters), allocatable  :: params_match
        type(cmdline) :: cline_match
        type(oris)    :: eulspace
        type(image)   :: img_msk
        integer     :: ncls, nstates, ldim(3), box_c, box_v, box_c2, box_v2, box_crop, s, icls, i, iref
        real        :: smpd_c, smpd_v, smpd_crop, mskrad, scale
        ncls    = size(cavgs)
        nstates = size(vols)
        if( ncls < 1 )    THROW_HARD('no class averages to match; match_cavgs_to_states')
        if( nstates < 1 ) THROW_HARD('no state volumes to match against; match_cavgs_to_states')
        if( nspace < 1 )  THROW_HARD('no projection direction; match_cavgs_to_states')
        ldim   = cavgs(1)%get_ldim()
        box_c  = ldim(1)
        smpd_c = cavgs(1)%get_smpd()
        ldim   = vols(1)%get_ldim()
        box_v  = ldim(1)
        smpd_v = vols(1)%get_smpd()
        do s = 2,nstates
            ldim = vols(s)%get_ldim()
            if( ldim(1) /= box_v .or. abs(vols(s)%get_smpd() - smpd_v) > 1.e-3 ) THROW_HARD('the state volumes differ in box or pixel size')
        enddo
        ! both sets brought to one pixel size, the coarsest of lp/LP_SAMPLING and the two sets' own
        ! (with boxes of BOX_MATCH_MIN at least), then to one box, the smaller, by clipping
        smpd_crop = max(lp / LP_SAMPLING, smpd_c, smpd_v)
        smpd_crop = max(min(smpd_crop, real(box_c) * smpd_c / real(BOX_MATCH_MIN), real(box_v) * smpd_v / real(BOX_MATCH_MIN)),&
            &smpd_c, smpd_v)
        box_c2    = min(box_c, round2even(real(box_c) * smpd_c / smpd_crop))
        box_v2    = min(box_v, round2even(real(box_v) * smpd_v / smpd_crop))
        box_crop  = min(box_c2, box_v2)
        scale     = real(box_c) / real(box_c2) ! the class averages' input pixels per matched pixel
        allocate(cavgs_crop(ncls), vols_crop(nstates))
        ! serial: large images transform on threads of their own
        do icls = 1,ncls
            call cavgs_crop(icls)%copy(cavgs(icls))
            call resample(cavgs_crop(icls), box_c, box_c2, [box_crop,box_crop,1])
        enddo
        do s = 1,nstates
            call vols_crop(s)%copy(vols(s))
            call resample(vols_crop(s), box_v, box_v2, [box_crop,box_crop,box_crop])
        enddo
        ! the alignment's parameters, for the matched images' box and pixel size (the polar transforms
        ! are sized by the box and the mask, as cluster_cavgs sizes its own); it also sets the threads
        ! the reprojection runs on
        call cline_match%set('box',     box_crop)
        call cline_match%set('smpd',    smpd_crop)
        call cline_match%set('mskdiam', mskdiam)
        call cline_match%set('nthr',    max(1, nthr))
        call cline_match%set('ctf',     'no')
        call cline_match%set('objfun',  'cc')
        call cline_match%set('oritype', 'cls2D')
        call cline_match%set('mkdir',   'no')
        call cline_match%set('lp',      lp)
        call cline_match%set('trs',     trs / scale)
        allocate(params_match)
        call params_match%new(cline_match, silent=.true.)
        ! the reprojections of every state over the whole sphere: references 1..nspace are state 1's
        call eulspace%new(nspace, is_ptcl=.false.)
        call eulspace%spiral()
        allocate(refs(nstates * nspace))
        do s = 1,nstates
            projs = reproject(vols_crop(s), eulspace)
            do i = 1,nspace
                call refs((s - 1) * nspace + i)%copy(projs(i))
            enddo
            call dealloc_imgarr(projs)
        enddo
        call dealloc_imgarr(vols_crop)
        ! both sets normalised within, and softly masked to, the particle
        mskrad = min(real(box_crop / 2) - COSMSKHALFWIDTH - 1., 0.5 * mskdiam / smpd_crop)
        call img_msk%new([box_crop,box_crop,1], smpd_crop)
        img_msk = 1.
        call img_msk%memoize_mask_coords
        call img_msk%mask2D_hard(mskrad)
        l_msk = img_msk%bin2logical()
        call img_msk%kill
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1,nstates * nspace
            call refs(i)%norm_within(l_msk)
            call refs(i)%mask2D_soft(mskrad, backgr=0.)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls = 1,ncls
            call cavgs_crop(icls)%norm_within(l_msk)
            call cavgs_crop(icls)%mask2D_soft(mskrad, backgr=0.)
        enddo
        !$omp end parallel do
        write(logfhandle,'(A,I6,A,I3,A,I6,A,I4,A,F6.2,A)') '>>> MATCHING ', ncls, ' CLASS AVERAGES AGAINST ', nstates,&
            &' STATES x ', nspace, ' DIRECTIONS (BOX ', box_crop, ', ', smpd_crop, ' A/PX)'
        ! every class average aligned with every reprojection, without mirrors (the grid holds every view)
        algninfo = match_imgs(params_match, HP_MATCH, lp, trs / scale, refs, cavgs_crop, l_nomirr=.true.)
        deallocate(params_match)
        call cline_match%kill
        call dealloc_imgarr(refs)
        call dealloc_imgarr(cavgs_crop)
        allocate(corrmat(nstates * nspace, ncls))
        corrmat = algninfo(:,:)%corr
        if( allocated(states)  ) deallocate(states)
        if( allocated(eulers)  ) deallocate(eulers)
        if( allocated(shifts)  ) deallocate(shifts)
        if( allocated(corrs)   ) deallocate(corrs)
        if( allocated(margins) ) deallocate(margins)
        allocate(states(ncls), best_ref(ncls), eulers(3,ncls), shifts(2,ncls), corrs(ncls), margins(ncls))
        call reduce_state_matches(corrmat, nstates, nspace, states, best_ref, corrs, margins)
        do icls = 1,ncls
            iref = best_ref(icls)
            eulers(:,icls) = matched_pose(eulspace%get_euler(iref - (states(icls) - 1) * nspace), algninfo(iref,icls)%e3)
            shifts(:,icls) = [algninfo(iref,icls)%x, algninfo(iref,icls)%y] * scale
        enddo
        call eulspace%kill
        deallocate(algninfo, corrmat, best_ref, l_msk)

    contains

        ! @p img (box @p box_in) Fourier-cropped to @p box_out, its pixel size set to the matched
        ! one, then clipped in real space to @p ldim_out
        subroutine resample( img, box_in, box_out, ldim_out )
            class(image), intent(inout) :: img
            integer,      intent(in)    :: box_in, box_out, ldim_out(3)
            integer :: ldim_img(3)
            ldim_img = img%get_ldim()
            if( box_out < box_in )then
                call img%fft()
                if( ldim_img(3) > 1 )then
                    call img%clip_inplace([box_out,box_out,box_out])
                else
                    call img%clip_inplace([box_out,box_out,1])
                endif
                call img%ifft()
            endif
            call img%set_smpd(smpd_crop)
            if( ldim_out(1) < box_out ) call img%clip_inplace(ldim_out)
        end subroutine resample

    end subroutine match_cavgs_to_states

    !> From the correlations @p corr of every reference (nstates blocks of @p nspace, state by state)
    !! with every class (columns), per class: its best state, best reference, best correlation, and
    !! the margin of that correlation over the best of any other state (0 with one state).
    pure subroutine reduce_state_matches( corr, nstates, nspace, best_state, best_ref, best_corr, margin )
        real,    intent(in)  :: corr(:,:)
        integer, intent(in)  :: nstates, nspace
        integer, intent(out) :: best_state(:), best_ref(:)
        real,    intent(out) :: best_corr(:), margin(:)
        real    :: state_best(nstates)
        logical :: l_other(nstates)
        integer :: icls, s, i0
        do icls = 1,size(corr, 2)
            do s = 1,nstates
                i0 = (s - 1) * nspace
                state_best(s) = maxval(corr(i0 + 1:i0 + nspace, icls))
            enddo
            best_state(icls) = maxloc(state_best, dim=1)
            i0 = (best_state(icls) - 1) * nspace
            best_ref(icls)   = i0 + maxloc(corr(i0 + 1:i0 + nspace, icls), dim=1)
            best_corr(icls)  = state_best(best_state(icls))
            margin(icls)     = 0.
            if( nstates > 1 )then
                l_other = .true.
                l_other(best_state(icls)) = .false.
                margin(icls) = best_corr(icls) - maxval(state_best, mask=l_other)
            endif
        enddo
    end subroutine reduce_state_matches

    !> The 3D pose of an image matched to the projection along @p eul_dir (its in-plane angle
    !! ignored) at in-plane rotation @p e3_match (the alignment's angle, degrees), in refine3D's
    !! convention: the direction's first two angles, and 360 less the rotation.
    pure function matched_pose( eul_dir, e3_match ) result( euls )
        real, intent(in) :: eul_dir(3), e3_match
        real :: euls(3)
        euls(1) = eul_dir(1)
        euls(2) = eul_dir(2)
        euls(3) = modulo(360. - e3_match, 360.)
    end function matched_pose

end module simple_cavgs_state_match
