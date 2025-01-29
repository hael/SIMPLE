! optimization(search)-based masking
module simple_opt_mask
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image, image_ptr
use simple_parameters, only: params_glob
implicit none
#include "simple_local_flags.inc"

public :: estimate_spher_mask
private

contains

    subroutine estimate_spher_mask( ref, targ, mskimg, mskfromto, best_msk )
        class(image), intent(inout) :: ref, targ, mskimg
        integer,      intent(in)    :: mskfromto(2)
        integer,      intent(out)   :: best_msk
        logical, allocatable :: lmask(:,:,:)
        type(image) :: targ_copy, targ_msk
        integer     :: box, ldim(3), imsk
        real        :: smpd, cur_cost, best_cost
        if( mskfromto(1) == mskfromto(2) )then  
            best_msk = mskfromto(1)
            return
        endif
        if( targ%is_ft() .or. ref%is_ft() ) THROW_HARD('Input ref & targ has to be in real-space representation')
        ldim  = targ%get_ldim()
        box   = ldim(1)
        smpd  = targ%get_smpd()
        lmask = mskimg%bin2logical()
        call targ_copy%copy(targ)
        call targ_msk%new(ldim, smpd)
        best_msk  = mskfromto(1)
        best_cost = huge(best_cost)
        if( ldim(3) > 1 )then ! 3D, so parallelize
            call targ_msk%set_wthreads(.true.)
            call ref%set_wthreads(.true.)
        endif
        do imsk = mskfromto(1), mskfromto(2)
            call targ_msk%copy_fast(targ_copy)
            call targ_msk%mask(real(imsk), 'hard')
            ! calculate squared Euclidean distance
            cur_cost = targ_msk%sqeuclid(ref, lmask)
            if( cur_cost <= best_cost )then
                best_cost = cur_cost
                best_msk  = imsk
            endif
        enddo
        call targ_copy%kill
        call targ_msk%kill
    end subroutine estimate_spher_mask

end module simple_opt_mask
