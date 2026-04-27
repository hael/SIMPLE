!@descr: fused projection helpers that operate on projector data without polarft class dependencies
module simple_projector_pft_batch
use simple_core_module_api
use simple_projector, only: projector
implicit none

public :: fproject_polar_batch, fproject_polar_batch_mirr
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
        integer :: iproj, irot, k, kloc, kfrom, kto, pftsz
        real    :: e_rotmat(3,3), loc(3), px, py
        kfrom = kfromto(1)
        kto   = kfromto(2)
        pftsz = size(polar_x,1)
        if( kfrom < 1 .or. kfrom > kto  )           THROW_HARD('invalid kfromto in fproject_polar_batch')
        if( size(refs_state,1) /= pftsz           ) THROW_HARD('refs_state dim(1) mismatch in fproject_polar_batch')
        if( size(refs_state,2) /= size(polar_x,2) ) THROW_HARD('refs_state dim(2) mismatch in fproject_polar_batch')
        if( size(polar_x,2)    /= (kto-kfrom+1)   ) THROW_HARD('polar_x k-span mismatch in fproject_polar_batch')
        if( size(polar_y,1)    /= pftsz .or. size(polar_y,2) /= size(polar_x,2) )then
            THROW_HARD('polar_x/polar_y shape mismatch in fproject_polar_batch')
        endif
        if( size(refs_state,3) /= nspace ) THROW_HARD('refs_state dim(3) mismatch in fproject_polar_batch')
        if( size(mask)         <  kto    ) THROW_HARD('mask too short in fproject_polar_batch')
        !$omp parallel do default(shared) private(iproj,e_rotmat,irot,k,kloc,loc,px,py) schedule(static) proc_bind(close)
        do iproj = 1, nspace
            e_rotmat = eulspace%get_mat(iproj)
            do k = kfrom, kto
                kloc = k - kfrom + 1
                if( mask(k) )then
                    do irot = 1, pftsz
                        ! 3D coordinates
                        ! simplified loc = matmul([polar_x(irot,kloc), polar_y(irot,kloc), 0.0], e_rotmat)
                        px     = real(polar_x(irot,kloc))
                        py     = real(polar_y(irot,kloc))
                        loc(1) = px*e_rotmat(1,1) + py*e_rotmat(2,1)
                        loc(2) = px*e_rotmat(1,2) + py*e_rotmat(2,2)
                        loc(3) = px*e_rotmat(1,3) + py*e_rotmat(2,3)
                        ! 3D interpolation
                        refs_state(irot,kloc,iproj) = vol_pad%interp_fcomp_oversamp(loc)
                    enddo
                else
                    refs_state(:,kloc,iproj) = CMPLX_ZERO
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine fproject_polar_batch

    !> Extract unique central sections and mirror them
    subroutine fproject_polar_batch_mirr( vol_pad, eulspace, nspace, kfromto, mask, polar_x, polar_y, refs_state )
        class(projector), intent(in)    :: vol_pad
        class(oris),      intent(inout) :: eulspace
        integer,          intent(in)    :: nspace, kfromto(2)
        logical,          intent(in)    :: mask(:)
        real(sp),         intent(in)    :: polar_x(:,:), polar_y(:,:)
        complex(sp),      intent(inout) :: refs_state(:,:,:)
        complex(sp) :: pftm(size(polar_x,1), kfromto(1):kfromto(2))
        integer     :: iproj, irot, k, kloc, kfrom, kto, m, l, pftsz
        real        :: e_rotmat(3,3), loc(3), px, py, psi
        kfrom = kfromto(1)
        kto   = kfromto(2)
        pftsz = size(polar_x,1)
        if( kfrom < 1 .or. kfrom > kto  )           THROW_HARD('invalid kfromto in fproject_polar_batch')
        if( size(refs_state,1) /= pftsz           ) THROW_HARD('refs_state dim(1) mismatch in fproject_polar_batch')
        if( size(refs_state,2) /= size(polar_x,2) ) THROW_HARD('refs_state dim(2) mismatch in fproject_polar_batch')
        if( size(polar_x,2)    /= (kto-kfrom+1)   ) THROW_HARD('polar_x k-span mismatch in fproject_polar_batch')
        if( size(polar_y,1)    /= pftsz .or. size(polar_y,2) /= size(polar_x,2) )then
            THROW_HARD('polar_x/polar_y shape mismatch in fproject_polar_batch')
        endif
        if( size(refs_state,3) /= nspace  ) THROW_HARD('refs_state dim(3) mismatch in fproject_polar_batch')
        if( size(mask)         <  kto     ) THROW_HARD('mask too short in fproject_polar_batch')
        if( .not.eulspace%isthere('mirr') ) THROW_HARD('Mirror index missing, use fproject_polar_batch instead')
        !$omp parallel do private(iproj,e_rotmat,irot,k,kloc,loc,m,px,py,l,pftm,psi)&
        !$omp& schedule(static) proc_bind(close) default(shared)
        do iproj = 1, nspace/2
            ! extract iproj section
            e_rotmat = eulspace%get_mat(iproj)
            do k = kfrom, kto
                kloc = k - kfrom + 1
                if( mask(k) )then
                    do irot = 1,pftsz
                        ! 3D coordinates
                        ! simplified loc = matmul([polar_x(irot,kloc), polar_y(irot,kloc), 0.0], e_rotmat)
                        px     = real(polar_x(irot,kloc))
                        py     = real(polar_y(irot,kloc))
                        loc(1) = px*e_rotmat(1,1) + py*e_rotmat(2,1)
                        loc(2) = px*e_rotmat(1,2) + py*e_rotmat(2,2)
                        loc(3) = px*e_rotmat(1,3) + py*e_rotmat(2,3)
                        ! 3D interpolation
                        refs_state(irot,kloc,iproj) = vol_pad%interp_fcomp_oversamp(loc)
                    enddo
                else
                    refs_state(:,kloc,iproj) = CMPLX_ZERO
                endif
            enddo
            ! mirror corresponding section
            m = eulspace%get_int(iproj, 'mirr')
            pftm(1,:) = conjg(refs_state(1,:,iproj))
            do k = 2,pftsz/2
                l = pftsz-k+2
                pftm(k,:) = refs_state(l,:,iproj)
                pftm(l,:) = refs_state(k,:,iproj)
            enddo
            k = pftsz/2 + 1
            pftm(k,:) = refs_state(k,:,iproj)
            ! update array
            psi = eulspace%get(m, 'psi')
            if( (psi > 0.1) .and. (psi < 359.9) )then
                refs_state(:,:,m) = conjg(pftm)
            else
                refs_state(:,:,m) = pftm
            endif
        enddo
        !$omp end parallel do
    end subroutine fproject_polar_batch_mirr

end module simple_projector_pft_batch
