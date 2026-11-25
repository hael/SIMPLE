!@descr: phase-local estimates and bounded particle-batch capacity for the matchers
module simple_matcher_memory
use, intrinsic :: iso_fortran_env, only: int64
use simple_parameters, only: parameters
use simple_type_defs,  only: ptcl_ref
implicit none
private
public :: matcher_memory_terms, matcher_batch_capacity, matcher_block_kernel_bytes

contains

    !> Dominant alignment arrays, including interpolation shells, scratch and pending candidates.
    !! The remaining process memory is reserved for project state, library plans and reconstruction.
    subroutine matcher_memory_terms( params, nrefs, nptcls, nthreads, box_sums, fixed_bytes, particle_bytes )
        class(parameters), intent(in)  :: params
        integer,           intent(in)  :: nrefs, nptcls, nthreads, box_sums
        integer(int64),    intent(out) :: fixed_bytes, particle_bytes
        type(ptcl_ref) :: entry
        integer(int64) :: p, nk, ni, nr, box, sum_box, crop, pad, nrows
        p    = int(max(1,params%pftsz),int64)
        nk   = int(max(1,params%kfromto(2)-params%kfromto(1)+1),int64)
        ni   = int(max(1,params%box_crop/2-params%kfromto(1)+1),int64)
        nr   = int(max(1,nrefs),int64)
        box  = int(params%box,int64)
        crop = int(params%box_crop,int64)
        pad  = int(params%box_croppd,int64)
        nrows = int(max(0,params%top-params%fromp+1),int64)
        fixed_bytes = 16_int64 * nr * p * ni + 32_int64 * nr * (p+1_int64) * nk
        fixed_bytes = fixed_bytes + int(nptcls,int64) * &
            &int(storage_size(entry)/8 + 16 + 64*max(1,params%topk),int64)
        fixed_bytes = fixed_bytes + 4_int64 * nrows * (box/2_int64)
        fixed_bytes = fixed_bytes + int(nthreads,int64) * &
            &(128_int64*(p+1_int64)*max(ni,nk) + 8_int64*(crop*crop+pad*pad))
        particle_bytes = 12_int64*p*ni + 16_int64*(p+1_int64)*nk + 8_int64*nk + 32_int64
        particle_bytes = particle_bytes + 4_int64*(box+2_int64)*box
        if( box_sums > 0 )then
            sum_box = int(box_sums,int64)
            fixed_bytes = fixed_bytes + 24_int64*nr*(sum_box/2_int64+1_int64)*sum_box
            particle_bytes = particle_bytes + 4_int64*(sum_box+2_int64)*sum_box + 72_int64*nr
        else
            particle_bytes = particle_bytes + 32_int64*nr
        endif
    end subroutine matcher_memory_terms

    !> The block kernel's memory, for a table that scores with it: the bin-major reference memo of both halves
    !! (memoize_ref_bins) and each thread's particle terms, tile, products, transposed products and transforms
    pure integer(int64) function matcher_block_kernel_bytes( params, nrefs, nthreads ) result( nbytes )
        use simple_polarft_calc, only: PFTC_BLOCK_NB, PFTC_BLOCK_NR
        class(parameters), intent(in) :: params
        integer,           intent(in) :: nrefs, nthreads
        integer(int64) :: p1, n2k, nb, nr
        p1  = int(max(1,params%pftsz),int64) + 1_int64
        n2k = 2_int64 * int(max(1,params%kfromto(2)-params%kfromto(1)+1),int64)
        nb  = int(PFTC_BLOCK_NB,int64)
        nr  = int(PFTC_BLOCK_NR,int64)
        nbytes = 16_int64 * n2k * int(max(1,nrefs),int64) * p1
        nbytes = nbytes + int(max(1,nthreads),int64) * (8_int64 * (n2k*nb*p1 + n2k*nr + 2_int64*nb*nr*p1) +&
            &8_int64 * (p1 - 1_int64) * nr * nb)
    end function matcher_block_kernel_bytes

    !> Keep the preferred capacity when it fits; an empty part still needs one calculator slot.
    !! A missing budget leaves capacity unchanged. A tight budget never reduces it below one particle
    !! per thread (within the sample and the preferred batch): a batch's parallel loops give each
    !! thread a particle, and the fixed storage does not scale with the batch, so a smaller batch
    !! saves a few slots of memory and idles the rest of the team. Fixed storage over budget leaves
    !! that floor and no cache.
    pure integer function matcher_batch_capacity( nptcls, preferred, nthreads, budget, fixed_bytes, particle_bytes ) &
        &result(capacity)
        integer,        intent(in) :: nptcls, preferred, nthreads
        integer(int64), intent(in) :: budget, fixed_bytes, particle_bytes
        integer(int64) :: available
        integer        :: team_floor
        capacity = max(1,min(nptcls,preferred))
        if( budget <= 0_int64 .or. particle_bytes <= 0_int64 ) return
        team_floor = max(1,min(capacity,nthreads))
        available  = max(0_int64,budget-max(0_int64,fixed_bytes))
        capacity   = int(max(int(team_floor,int64),min(int(capacity,int64),available/particle_bytes)))
    end function matcher_batch_capacity

end module simple_matcher_memory
