!@descr: reference preparation helpers for matcher workflows
module simple_matcher_pftc_prep
use simple_core_module_api
use simple_builder,              only: builder
use simple_classaverager,        only: cavgs_merged, cavgs_even, cavgs_odd
use simple_image,                only: image
use simple_matcher_ptcl_batch,   only: prep_sigmas_objfun
use simple_parameters,           only: parameters
implicit none

public :: prep_pftc4align2D
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

end module simple_matcher_pftc_prep
