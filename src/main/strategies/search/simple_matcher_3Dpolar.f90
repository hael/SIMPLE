!@descr: polar partial reference generation for reconstruct3D workers
module simple_matcher_3Dpolar
use simple_core_module_api
use simple_builder,            only: builder
use simple_matcher_3Drec,      only: init_rec, prep_imgs4rec, update_rec, write_partial_recs, finalize_rec_objs
use simple_matcher_ptcl_batch, only: alloc_ptcl_imgs, build_batch_particles3D, clean_batch_particles3D
use simple_parameters,         only: parameters
use simple_pftc_srch_api
implicit none

public :: calc_polar_partials
private
#include "simple_local_flags.inc"

contains

    subroutine calc_polar_partials( params, build, cline, nptcls, pinds )
        use simple_image, only: image
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: nptcls
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable   :: fpls(:)
        type(image),       allocatable   :: ptcl_match_imgs(:), ptcl_match_imgs_pad(:), ptcl_rec_imgs(:)
        integer :: batchlims(2), ibatch, batchsz, i, maxbatchsz, nrefs
        maxbatchsz = MAXIMGBATCHSZ
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1,maxbatchsz], params%kfromto)
        call build%pftc%polar_cavger_new(.true., nrefs=nrefs)
        call alloc_ptcl_imgs(params, build, ptcl_match_imgs, ptcl_match_imgs_pad, maxbatchsz)
        if( trim(params%polar) == 'obsfield' )then
            call init_rec(params, build, maxbatchsz, fpls)
            call alloc_imgarr(maxbatchsz, [params%box,params%box,1], params%smpd, ptcl_rec_imgs)
        endif
        do ibatch = 1,nptcls,maxbatchsz
            batchlims = [ibatch, min(nptcls, ibatch+maxbatchsz-1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            select case(trim(params%polar))
                case('obsfield')
                    call build_batch_particles3D(params, build, batchsz, pinds(batchlims(1):batchlims(2)), &
                        &ptcl_match_imgs, ptcl_match_imgs_pad, imgs4rec=ptcl_rec_imgs(:batchsz))
                case('yes')
                    call build_batch_particles3D(params, build, batchsz, pinds(batchlims(1):batchlims(2)), &
                        &ptcl_match_imgs, ptcl_match_imgs_pad)
                case default
                    THROW_HARD('unsupported POLAR mode for calc_polar_partials: '//trim(params%polar))
            end select
            select case(trim(params%polar))
                case('obsfield')
                    call prep_imgs4rec(params, build, batchsz, ptcl_rec_imgs(:batchsz), &
                        &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
                    call update_rec(params, build, batchsz, pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
                case('yes')
                    call build%pftc%polar_cavger_update_sums(batchsz, pinds(batchlims(1):batchlims(2)), &
                        &build%spproj, is3D=.true.)
                case default
                    THROW_HARD('unsupported POLAR mode for calc_polar_partials: '//trim(params%polar))
            end select
        end do
        select case(trim(params%polar))
            case('obsfield')
                call write_partial_recs(params, build, cline, fpls)
            case('yes')
                call build%pftc%polar_cavger_readwrite_partial_sums('write')
        end select
        if( allocated(fpls) )then
            do i = 1,size(fpls)
                if( allocated(fpls(i)%cmplx_plane) ) deallocate(fpls(i)%cmplx_plane)
                if( allocated(fpls(i)%ctfsq_plane) ) deallocate(fpls(i)%ctfsq_plane)
            end do
            deallocate(fpls)
        endif
        call clean_batch_particles3D(build, ptcl_match_imgs, ptcl_match_imgs_pad, &
            &imgs4rec=ptcl_rec_imgs)
        if( trim(params%polar) == 'obsfield' ) call finalize_rec_objs(params, build)
        call build%pftc%polar_cavger_kill
        call build%pftc%kill
    end subroutine calc_polar_partials

end module simple_matcher_3Dpolar
