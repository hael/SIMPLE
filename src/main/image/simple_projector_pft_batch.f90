!@descr: fused projection helpers that operate on projector data without polarft class dependencies
module simple_projector_pft_batch
use simple_core_module_api
use simple_projector, only: projector
implicit none

public :: fproject_polar_batch
private
#include "simple_local_flags.inc"

contains

    !> Fused oversampled projection of one state's full reference bank.
    !> refs_state is indexed (irot, k_local, iproj) where k_local = k-kfromto(1)+1.
    subroutine fproject_polar_batch( vol_pad, eulspace, nspace, kfromto, mask, polar_x, polar_y, refs_state )
        class(projector), intent(in)    :: vol_pad
        class(oris),      intent(inout) :: eulspace
        integer,          intent(in)    :: nspace, kfromto(2)
        logical,          intent(in)    :: mask(:)
        real(sp),         intent(in)    :: polar_x(:,:), polar_y(:,:)
        complex(sp),      intent(inout) :: refs_state(:,:,:)
        type(ori)   :: o_tmp
        integer     :: iproj, irot, k, kloc, kfrom, kto
        real        :: e_rotmat(3,3), loc(3)
        kfrom = kfromto(1)
        kto   = kfromto(2)
        if( size(refs_state,1) /= size(polar_x,1) ) THROW_HARD('refs_state dim(1) mismatch in fproject_polar_batch')
        if( size(refs_state,2) /= size(polar_x,2) ) THROW_HARD('refs_state dim(2) mismatch in fproject_polar_batch')
        if( size(polar_x,2)    /= (kto-kfrom+1)   ) THROW_HARD('polar_x k-span mismatch in fproject_polar_batch')
        if( size(polar_y,1)    /= size(polar_x,1) .or. size(polar_y,2) /= size(polar_x,2) )then
            THROW_HARD('polar_x/polar_y shape mismatch in fproject_polar_batch')
        endif
        if( size(refs_state,3) /= nspace ) THROW_HARD('refs_state dim(3) mismatch in fproject_polar_batch')
        if( size(mask)         <  kto    ) THROW_HARD('mask too short in fproject_polar_batch')
        !$omp parallel do default(shared) private(iproj,o_tmp,e_rotmat,irot,k,kloc,loc) schedule(static) proc_bind(close)
        do iproj = 1, nspace
            call eulspace%get_ori(iproj, o_tmp)
            e_rotmat = o_tmp%get_mat()
            do irot = 1, size(polar_x,1)
                do k = kfrom, kto
                    kloc = k - kfrom + 1 ! this is the correct indexing since we pass the arrays without explicit shape 
                    if( mask(k) )then
                        loc = matmul([polar_x(irot,kloc), polar_y(irot,kloc), 0.0], e_rotmat)
                        refs_state(irot,kloc,iproj) = vol_pad%interp_fcomp_oversamp(loc)
                    else
                        refs_state(irot,kloc,iproj) = CMPLX_ZERO
                    endif
                enddo
            enddo
            call o_tmp%kill
        enddo
        !$omp end parallel do
    end subroutine fproject_polar_batch

end module simple_projector_pft_batch
