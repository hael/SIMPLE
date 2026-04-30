!@descr: reference-section preparation helpers for matcher workflows
! POLAR_REFS*.bin lifecycle:
! - producers: centralized materializers in simple_matcher_refvol_utils, exec_polar_assembly
! - consumers: prep_pftc4align3D, prob_tab, prob_tab_neigh
module simple_matcher_pftc_prep
use simple_pftc_srch_api
use simple_builder,              only: builder
use simple_classaverager,        only: cavgs_merged, cavgs_even, cavgs_odd
use simple_matcher_ptcl_batch,   only: prep_sigmas_objfun
use simple_matcher_refvol_utils, only: pick_lp_est_state, estimate_lp_from_refs, polar_ref_sections_available
use simple_refine3D_fnames,      only: refine3D_polar_refs_fname
implicit none

public :: prep_pftc4align3D, prep_pftc4align2D, polar_ref_sections_available
private
#include "simple_local_flags.inc"

contains

    !>  \brief  prepares the polarft corrcalc object for search and imports the references
    subroutine prep_pftc4align2D( params, build, ptcl_match_imgs_pad, batchsz_max, which_iter, l_stream )
        use simple_matcher_2Dprep, only: prep2dref, calc_2Dref_offset
        class(parameters),          intent(inout) :: params
        class(builder),             intent(inout) :: build
        type(image),                intent(inout) :: ptcl_match_imgs_pad(:)
        integer,                    intent(in)    :: batchsz_max, which_iter
        logical,                    intent(in)    :: l_stream
        class(image), pointer :: cavgs_m(:), cavgs_e(:), cavgs_o(:)
        type(image), allocatable :: match_imgs(:)
        complex,     allocatable :: pft(:,:)
        real         :: xyz(3)
        integer      :: icls, pop, pop_even, pop_odd, centype, ithr
        logical      :: do_center, has_been_searched, input_center
        has_been_searched = .not.build%spproj%is_virgin_field(params%oritype)
        input_center      = trim(params%center) .eq. 'yes'
        ! create the polarft_calc object
        call build%pftc%new(params, params%ncls, [1,batchsz_max], params%kfromto)
        ! objective functions & sigma
        call prep_sigmas_objfun(params, build, l_stream)
        ! prepare the polarizer images
        call ptcl_match_imgs_pad(1)%memoize4polarize_oversamp(build%pftc%get_pdim_srch())
        allocate(match_imgs(params%ncls))
        cavgs_m => cavgs_merged
        cavgs_e => cavgs_even
        cavgs_o => cavgs_odd
        call cavgs_m(1)%construct_thread_safe_tmp_imgs(nthr_glob)
        ! mask memoization
        call cavgs_m(1)%memoize_mask_coords
        ! mode of cavg centering
        centype = get_centype(params%center_type)
        ! PREPARATION OF REFERENCES IN pftc
        ! read references and transform into polar coordinates
        !$omp parallel do default(shared) private(icls,ithr,pop,pop_even,pop_odd,do_center,xyz,pft)&
        !$omp schedule(static) proc_bind(close)
        do icls=1,params%ncls
            pop      = 1
            pop_even = 0
            pop_odd  = 0
            if( has_been_searched )then
                pop      = build%spproj_field%get_pop(icls, 'class'      )
                pop_even = build%spproj_field%get_pop(icls, 'class', eo=0)
                pop_odd  = build%spproj_field%get_pop(icls, 'class', eo=1)
            endif
            ithr = omp_get_thread_num() + 1
            if( pop > 0 )then
                call match_imgs(icls)%new([params%box_crop, params%box_crop, 1], params%smpd_crop, wthreads=.false.)
                ! Calculate center
                do_center = has_been_searched .and. (pop > MINCLSPOPLIM) .and. (which_iter > 2)
                do_center = input_center .and. do_center
                if( do_center )then
                    call match_imgs(icls)%copy_fast(cavgs_m(icls))
                    call calc_2Dref_offset(params, build, match_imgs(icls), icls, centype, xyz)
                else
                    xyz = 0.0
                endif
                ! Prepare the references
                ! allocte pft
                pft = build%pftc%allocate_pft()
                if( params%l_lpset )then
                    ! merged class average in both even and odd positions
                    call match_imgs(icls)%copy_fast(cavgs_m(icls))
                    call prep2Dref(params, build, match_imgs(icls), icls, xyz, ptcl_match_imgs_pad(ithr))
                    call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft)
                    call build%pftc%set_ref_pft(icls, pft, iseven=.true.)
                    call build%pftc%cp_even2odd_ref(icls)
                else
                    if( pop_even >= MINCLSPOPLIM .and. pop_odd >= MINCLSPOPLIM )then
                        ! even & odd
                        call match_imgs(icls)%copy_fast(cavgs_e(icls))
                        call prep2Dref(params, build, match_imgs(icls), icls, xyz, ptcl_match_imgs_pad(ithr))
                        call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft)
                        call build%pftc%set_ref_pft(icls, pft, iseven=.true.)
                        call match_imgs(icls)%copy_fast(cavgs_o(icls))
                        call prep2Dref(params, build, match_imgs(icls), icls, xyz, ptcl_match_imgs_pad(ithr))
                        call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft)
                        call build%pftc%set_ref_pft(icls, pft, iseven=.false.)
                    else
                        ! merged class average in both even and odd positions
                        call match_imgs(icls)%copy_fast(cavgs_m(icls))
                        call prep2Dref(params, build, match_imgs(icls), icls, xyz, ptcl_match_imgs_pad(ithr))
                        call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft)
                        call build%pftc%set_ref_pft(icls, pft, iseven=.true.)
                        call build%pftc%cp_even2odd_ref(icls)
                    endif
                endif
                deallocate(pft)
                call match_imgs(icls)%kill
            endif
        end do
        !$omp end parallel do
        call build%pftc%memoize_refs
        ! CLEANUP
        deallocate(match_imgs)
        call cavgs_m(1)%kill_thread_safe_tmp_imgs
        nullify(cavgs_m,cavgs_e,cavgs_o)
    end subroutine prep_pftc4align2D

    subroutine prep_pftc4align3D( params, build, cline, batchsz )
        class(parameters),        intent(inout) :: params
        class(builder),           intent(inout) :: build
        class(cmdline),           intent(in)    :: cline
        integer,                  intent(in)    :: batchsz
        real, allocatable :: gaufilter(:), default_frc(:)
        integer           :: iproj, nrefs, filtsz, state, istate
        logical           :: l_filtrefs
        ! Resolution limit estimation
        call pick_lp_est_state(params, build, state)
        if( trim(params%filt_mode).eq.'uniform' .and. &
            &cline%defined('lpstart') .and. cline%defined('lpstop') )then
            call estimate_lp_from_refs(params, build, cline, params%lpstart, params%lpstop, state)
        endif
        ! Calculator init
        ! Polar reference slots are state-major:
        ! (state - 1) * nspace + local_projection.
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1,batchsz], params%kfromto)
        ! objective functions & sigma
        call prep_sigmas_objfun(params, build, .false.)
        ! Read polar references
        call build%pftc%polar_cavger_new(.true.)
        call build%pftc%polar_cavger_read_all(refine3D_polar_refs_fname())
        if( file_exists(FRCS_FILE) )then
            call build%clsfrcs%read(string(FRCS_FILE))
        else
            call build%clsfrcs%new(params%nspace, params%box_crop, params%smpd_crop, params%nstates)
            ! Bootstrap before the first assembly has no empirical FSCs yet.
            ! Use a neutral FRC so reference preparation does not inherit a
            ! silent zero-resolution prior.
            filtsz = build%clsfrcs%get_filtsz()
            allocate(default_frc(filtsz), source=1.0)
            do istate = 1,params%nstates
                do iproj = 1,params%nspace
                    call build%clsfrcs%set_frc(iproj, default_frc, istate)
                enddo
            enddo
            deallocate(default_frc)
        endif
        ! prepare filter
        l_filtrefs = .false.
        if(trim(params%gauref).eq.'yes')then
            l_filtrefs = .true.
            filtsz = build%clsfrcs%get_filtsz()
            allocate(gaufilter(filtsz),source=0.)
            call gaussian_filter(params%gaufreq, params%smpd, params%box, gaufilter)
        endif
        ! PREPARATION OF REFERENCES IN pftc
        !$omp parallel do default(shared) private(iproj)&
        !$omp schedule(static) proc_bind(close)
        do iproj = 1,nrefs
            if( l_filtrefs ) call build%pftc%polar_filterrefs(iproj, gaufilter)
            ! transfer to pftc
            if( params%l_lpset )then
                ! put the merged class average in both even and odd positions
                call build%pftc%polar_cavger_set_ref_pft(iproj, 'merged')
                call build%pftc%cp_even2odd_ref(iproj)
            else
                ! transfer e/o refs to pftc
                call build%pftc%polar_cavger_set_ref_pft(iproj, 'even')
                call build%pftc%polar_cavger_set_ref_pft(iproj, 'odd')
            endif
        end do
        !$omp end parallel do
        if( allocated(gaufilter) ) deallocate(gaufilter)
    end subroutine prep_pftc4align3D

end module simple_matcher_pftc_prep
