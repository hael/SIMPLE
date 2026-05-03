!@descr: polar=yes partial reference generation for reconstruct3D workers
module simple_matcher_3Dpolar
use simple_core_module_api
use simple_builder,            only: builder
use simple_matcher_ptcl_batch, only: alloc_ptcl_imgs, build_batch_particles3D, clean_batch_particles3D
use simple_parameters,         only: parameters
use simple_pftc_srch_api
implicit none

public :: calc_polar_partials
private
#include "simple_local_flags.inc"

contains

    subroutine calc_polar_partials( params, build, nptcls, pinds )
        use simple_image, only: image
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        integer,           intent(in)    :: pinds(nptcls)
        type(image),       allocatable   :: ptcl_match_imgs(:), ptcl_match_imgs_pad(:)
        integer :: batchlims(2), ibatch, batchsz, maxbatchsz, nrefs
        if( trim(params%polar) /= 'yes' )then
            THROW_HARD('calc_polar_partials requires POLAR=yes; polar=obsfield uses calc_3Drec')
        endif
        maxbatchsz = MAXIMGBATCHSZ
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1,maxbatchsz], params%kfromto)
        call build%pftc%polar_cavger_new(.true., nrefs=nrefs)
        call alloc_ptcl_imgs(params, build, ptcl_match_imgs, ptcl_match_imgs_pad, maxbatchsz)
        do ibatch = 1,nptcls,maxbatchsz
            batchlims = [ibatch, min(nptcls, ibatch+maxbatchsz-1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call build_batch_particles3D(params, build, batchsz, pinds(batchlims(1):batchlims(2)), &
                &ptcl_match_imgs, ptcl_match_imgs_pad)
            call build%pftc%polar_cavger_update_sums(batchsz, pinds(batchlims(1):batchlims(2)), &
                &build%spproj, is3D=.true.)
        end do
        call build%pftc%polar_cavger_readwrite_partial_sums('write')
        call clean_batch_particles3D(build, ptcl_match_imgs, ptcl_match_imgs_pad)
        call build%pftc%polar_cavger_kill
        call build%pftc%kill
    end subroutine calc_polar_partials

end module simple_matcher_3Dpolar
