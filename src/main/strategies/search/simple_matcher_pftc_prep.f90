!@descr: polar-reference preparation helpers for matcher workflows
module simple_matcher_pftc_prep
use simple_pftc_srch_api
use simple_builder,              only: builder
use simple_classaverager,        only: cavgs_merged, cavgs_even, cavgs_odd
use simple_matcher_ptcl_batch,   only: prep_sigmas_objfun
use simple_matcher_refvol_utils, only: report_resolution, estimate_lp_from_refs
implicit none

public :: prep_pftc4align3D_polar, prep_pftc4align2D, prep_pftc4align2D_polar
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
                do_center = (has_been_searched .and. (pop > MINCLSPOPLIM) .and. (which_iter > 2)&
                    &.and. .not.params%l_update_frac)
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
                    call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft, mask=build%l_resmsk)
                    call build%pftc%set_ref_pft(icls, pft, iseven=.true.)
                    call build%pftc%cp_even2odd_ref(icls)
                else
                    if( pop_even >= MINCLSPOPLIM .and. pop_odd >= MINCLSPOPLIM )then
                        ! even & odd
                        call match_imgs(icls)%copy_fast(cavgs_e(icls))
                        call prep2Dref(params, build, match_imgs(icls), icls, xyz, ptcl_match_imgs_pad(ithr))
                        call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft, mask=build%l_resmsk)
                        call build%pftc%set_ref_pft(icls, pft, iseven=.true.)
                        call match_imgs(icls)%copy_fast(cavgs_o(icls))
                        call prep2Dref(params, build, match_imgs(icls), icls, xyz, ptcl_match_imgs_pad(ithr))
                        call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft, mask=build%l_resmsk)
                        call build%pftc%set_ref_pft(icls, pft, iseven=.false.)
                    else
                        ! merged class average in both even and odd positions
                        call match_imgs(icls)%copy_fast(cavgs_m(icls))
                        call prep2Dref(params, build, match_imgs(icls), icls, xyz, ptcl_match_imgs_pad(ithr))
                        call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft, mask=build%l_resmsk)
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

    !>  \brief  prepares the polarft corrcalc object for search and imports polar references
    subroutine prep_pftc4align2D_polar( params, build, batchsz_max, which_iter, l_stream )
        use simple_matcher_2Dprep, only: calc_2Dref_offset
        class(parameters),          intent(inout) :: params
        class(builder),             intent(inout) :: build
        integer,                    intent(in)    :: batchsz_max, which_iter
        logical,                    intent(in)    :: l_stream
        type(image), allocatable :: tmp_imgs(:)
        real    :: xyz(3)
        integer :: icls, pop, pop_even, pop_odd, centype
        logical :: has_been_searched, do_center, l_center, l_gaufilt
        ! pftc instantiation
        call build%pftc%new(params, params%ncls, [1,batchsz_max], params%kfromto)
        ! Sigma2
        call prep_sigmas_objfun(params, build, l_stream)
        ! Read polar references
        call build%pftc%polar_cavger_new(.false.)
        call build%pftc%polar_cavger_read_all(string(POLAR_REFS_FBODY)//BIN_EXT)
        has_been_searched = .not.build%spproj%is_virgin_field(params%oritype)
        ! Centering-related objects
        do_center = (params%center .eq. 'yes') .and. has_been_searched&
             &.and. (which_iter > 2) .and. (.not.params%l_update_frac)
        if( do_center )then
            allocate(tmp_imgs(params%ncls))
            call build%pftc%polar_cavger_refs2cartesian(tmp_imgs, 'merged')
            call tmp_imgs(1)%construct_thread_safe_tmp_imgs(nthr_glob)
            ! mask memoization
            call tmp_imgs(1)%memoize_mask_coords
        endif
        ! Filtering
        l_gaufilt = trim(params%gauref)=='yes'
        ! Mode of cavg centering
        centype = get_centype(params%center_type)
        ! PREPARATION OF REFERENCES IN pftc
        !$omp parallel do default(shared) private(icls,pop,pop_even,pop_odd,xyz,l_center)&
        !$omp schedule(static) proc_bind(close)
        do icls=1,params%ncls
            ! populations
            pop      = 1
            pop_even = 0
            pop_odd  = 0
            if( has_been_searched )then
                pop      = build%spproj_field%get_pop(icls, 'class'      )
                pop_even = build%spproj_field%get_pop(icls, 'class', eo=0)
                pop_odd  = build%spproj_field%get_pop(icls, 'class', eo=1)
            endif
            if( pop > 0 )then
                ! centering
                l_center = do_center .and. (pop > MINCLSPOPLIM)
                xyz      = 0.
                if( l_center ) call calc_2Dref_offset(params, build, tmp_imgs(icls), icls, centype, xyz)
                ! Prep for alignment
                call build%pftc%polar_prep2Dref(build%clsfrcs, icls, l_gaufilt)
                ! transfer to pftc
                if( params%l_lpset )then
                    ! merged class average in both even and odd positions
                    call build%pftc%polar_cavger_set_ref_pft(icls, 'merged')
                    call build%pftc%cp_even2odd_ref(icls)
                else
                    if( pop_even >= MINCLSPOPLIM .and. pop_odd >= MINCLSPOPLIM )then
                        ! transfer e/o refs to pftc
                        call build%pftc%polar_cavger_set_ref_pft(icls, 'even')
                        call build%pftc%polar_cavger_set_ref_pft(icls, 'odd')
                    else
                        ! merged class average in both even and odd positions
                        call build%pftc%polar_cavger_set_ref_pft(icls, 'merged')
                        call build%pftc%cp_even2odd_ref(icls)
                    endif
                endif
                ! centering cavg & particles within the pftc
                if( l_center .and. (arg(xyz) > CENTHRESH) )then
                    call build%spproj_field%add_shift2class(icls, -xyz(1:2))
                    call build%pftc%shift_ref(icls, xyz(1:2))
                endif
            endif
            if( do_center ) call tmp_imgs(icls)%kill
        end do
        !$omp end parallel do
        call build%pftc%memoize_refs
        ! cleanup
        if( do_center )then
            call tmp_imgs(1)%kill_thread_safe_tmp_imgs
            deallocate(tmp_imgs)
        endif
    end subroutine prep_pftc4align2D_polar

    subroutine prep_pftc4align3D_polar( params, build, cline, batchsz )
        class(parameters),        intent(inout) :: params
        class(builder),           intent(inout) :: build
        class(cmdline),           intent(in)    :: cline
        integer,                  intent(in)    :: batchsz
        real, allocatable :: gaufilter(:)
        integer           :: iproj, nrefs, filtsz, state
        logical           :: l_filtrefs
        ! Resolution limit estimation
        call report_resolution(params, build, state)
        if( cline%defined('lpstart') .and. cline%defined('lpstop') )then
            call estimate_lp_from_refs(params, build, cline, params%lpstart, params%lpstop, state)
        endif
        ! Calculator init
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1,batchsz], params%kfromto)
        ! Read polar references
        call build%pftc%polar_cavger_new(.true.)
        call build%pftc%polar_cavger_read_all(string(POLAR_REFS_FBODY//BIN_EXT))
        call build%clsfrcs%read(string(FRCS_FILE))
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
        do iproj = 1,params%nspace
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
    end subroutine prep_pftc4align3D_polar

end module simple_matcher_pftc_prep