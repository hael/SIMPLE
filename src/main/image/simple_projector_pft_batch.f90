!@descr: fused projection helpers that operate on projector data without polarft class dependencies
module simple_projector_pft_batch
use simple_core_module_api
use simple_projector, only: projector
implicit none

public :: fproject_polar_batch, fproject_polar_batch_mirr, fproject_polar_batch_opt
private
#include "simple_local_flags.inc"

contains

    !> Fused oversampled projection of one state's full reference bank.
    !> refs_state is indexed (irot, k_local, iproj) where k_local = k-kfromto(1)+1.
    subroutine fproject_polar_batch( vol_pad, eulspace, nspace, kfromto, polar_x, polar_y, refs_state )
        class(projector), intent(in)    :: vol_pad
        class(oris),      intent(inout) :: eulspace
        integer,          intent(in)    :: nspace, kfromto(2)
        real(sp),         intent(in)    :: polar_x(:,:), polar_y(:,:)
        complex(sp),      intent(inout) :: refs_state(:,:,:)
        integer :: iproj, irot, k, kloc, kfrom, kto, pftsz, noris
        real    :: e_rotmat(3,3), loc(3), px, py
        kfrom = kfromto(1)
        kto   = kfromto(2)
        pftsz = size(polar_x,1)
        noris = eulspace%get_noris()
        if( kfrom < 1 .or. kfrom > kto  )           THROW_HARD('invalid kfromto in fproject_polar_batch')
        if( noris < nspace )                        THROW_HARD('eulspace smaller than requested nspace in fproject_polar_batch')
        if( size(refs_state,1) /= pftsz           ) THROW_HARD('refs_state dim(1) mismatch in fproject_polar_batch')
        if( size(refs_state,2) /= size(polar_x,2) ) THROW_HARD('refs_state dim(2) mismatch in fproject_polar_batch')
        if( size(polar_x,2)    /= (kto-kfrom+1)   ) THROW_HARD('polar_x k-span mismatch in fproject_polar_batch')
        if( size(polar_y,1)    /= pftsz .or. size(polar_y,2) /= size(polar_x,2) )then
            THROW_HARD('polar_x/polar_y shape mismatch in fproject_polar_batch')
        endif
        if( size(refs_state,3) /= nspace ) THROW_HARD('refs_state dim(3) mismatch in fproject_polar_batch')
        !$omp parallel do default(shared) private(iproj,e_rotmat,irot,k,kloc,loc,px,py) schedule(static) proc_bind(close)
        do iproj = 1, nspace
            e_rotmat = eulspace%get_mat(iproj)
            do k = kfrom, kto
                kloc = k - kfrom + 1
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
            enddo
        enddo
        !$omp end parallel do
    end subroutine fproject_polar_batch

    !> Extract unique central sections and mirror them
    subroutine fproject_polar_batch_mirr( vol_pad, eulspace, nspace, kfromto, polar_x, polar_y, refs_state )
        class(projector), intent(in)    :: vol_pad
        class(oris),      intent(inout) :: eulspace
        integer,          intent(in)    :: nspace, kfromto(2)
        real(sp),         intent(in)    :: polar_x(:,:), polar_y(:,:)
        complex(sp),      intent(inout) :: refs_state(:,:,:)
        complex(sp) :: pftm(size(polar_x,1), kfromto(1):kfromto(2))
        integer     :: iproj, irot, k, kloc, kfrom, kto, m, l, pftsz, noris
        real        :: e_rotmat(3,3), loc(3), px, py, psi
        kfrom = kfromto(1)
        kto   = kfromto(2)
        pftsz = size(polar_x,1)
        noris = eulspace%get_noris()
        if( kfrom < 1 .or. kfrom > kto  )           THROW_HARD('invalid kfromto in fproject_polar_batch')
        if( noris < nspace )                        THROW_HARD('eulspace smaller than requested nspace in fproject_polar_batch_mirr')
        if( size(refs_state,1) /= pftsz           ) THROW_HARD('refs_state dim(1) mismatch in fproject_polar_batch')
        if( size(refs_state,2) /= size(polar_x,2) ) THROW_HARD('refs_state dim(2) mismatch in fproject_polar_batch')
        if( size(polar_x,2)    /= (kto-kfrom+1)   ) THROW_HARD('polar_x k-span mismatch in fproject_polar_batch')
        if( size(polar_y,1)    /= pftsz .or. size(polar_y,2) /= size(polar_x,2) )then
            THROW_HARD('polar_x/polar_y shape mismatch in fproject_polar_batch')
        endif
        if( size(refs_state,3) /= nspace  ) THROW_HARD('refs_state dim(3) mismatch in fproject_polar_batch')
        if( .not.eulspace%isthere('mirr') ) THROW_HARD('Mirror index missing, use fproject_polar_batch instead')
        !$omp parallel do private(iproj,e_rotmat,irot,k,kloc,loc,m,px,py,l,pftm,psi)&
        !$omp& schedule(static) proc_bind(close) default(shared)
        do iproj = 1, nspace/2
            ! extract iproj section
            e_rotmat = eulspace%get_mat(iproj)
            do k = kfrom, kto
                kloc = k - kfrom + 1
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

    !> Fully fused oversampled projection of one state's full reference bank.
    !> OpenMP CPU & GPU accelerated
    !> refs_state is indexed (irot, k_local, iproj) where k_local = k-kfromto(1)+1.
    subroutine fproject_polar_batch_opt( vol_pad, eulspace, nspace, pftc, refs_state )
        use simple_kbinterpol,   only: apod_fast_device
        use simple_polarft_calc, only: polarft_calc
        class(projector),    intent(in)    :: vol_pad
        class(oris),         intent(in)    :: eulspace
        integer,             intent(in)    :: nspace
        class(polarft_calc), intent(inout) :: pftc
        complex(sp),         intent(inout) :: refs_state(:,:,:)
        integer,       parameter :: PAD_FAC         = OSMPL_PAD_FAC
        real,          parameter :: PAD_FAC_SCALING = real(PAD_FAC**3)
        type(kbinterpol)         :: kb
        real,        allocatable :: rotmats(:,:,:),hkcoords(:,:,:), tmp(:,:)
        real    :: rotmat(3,3)
        integer :: clb(3), cub(3), cdim(3), kfromto(2)
        integer :: kfrom, kto, iproj, pftsz, noris, iwinsz, wdim, nk
        ! Dimensions and boundaries
        pftsz   = pftc%get_pftsz()
        kfromto = pftc%get_kfromto()
        kfrom   = kfromto(1)
        kto     = kfromto(2)
        nk      = kto - kfrom + 1
        noris   = eulspace%get_noris()
        if( kfrom < 1 .or. kfrom > kto  )  THROW_HARD('invalid kfromto in fproject_polar_batch_opt')
        if( noris < nspace )               THROW_HARD('eulspace smaller than requested nspace in fproject_polar_batch_opt')
        if( size(refs_state,1) /= pftsz  ) THROW_HARD('refs_state dim(1) mismatch in fproject_polar_batch_opt')
        if( size(refs_state,2) /= nk     ) THROW_HARD('refs_state dim(2) mismatch in fproject_polar_batch_opt')
        if( size(refs_state,3) /= nspace ) THROW_HARD('refs_state dim(3) mismatch in fproject_polar_batch_opt')
        kb     = vol_pad%get_kbwin()
        iwinsz = ceiling(kb%get_winsz() - 0.5)
        wdim   = kb%get_wdim()
        if( wdim > 3 ) THROW_HARD('increase fixed KB weight storage in fproject_polar_batch_opt')
        allocate(tmp(kfrom:kto,pftsz), rotmats(2,3,nspace), hkcoords(2,pftsz,nk))
        ! format polar coordinates
        call pftc%get_polar_coords(1, tmp)
        hkcoords(1,:,:) = transpose(tmp)
        call pftc%get_polar_coords(2, tmp)
        hkcoords(2,:,:) = transpose(tmp)
        deallocate(tmp)
        ! padded polar coordinates
        hkcoords = real(PAD_FAC) * hkcoords
        ! format rotation matrices
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iproj,rotmat)
        do iproj = 1, nspace
            rotmat = eulspace%get_mat(iproj)
            rotmats(:,:,iproj) = rotmat(1:2,:)
        enddo
        !$omp end parallel do
        ! 3D array boundaries
        clb  = lbound(vol_pad%cmat_exp)
        cub  = ubound(vol_pad%cmat_exp)
        cdim = cub - clb + 1
        ! Interpolation kernel
        call fproject_polar_batch_opt_kernel(vol_pad%cmat_exp, refs_state)
        ! cleanup
        deallocate(rotmats, hkcoords)
        contains

            subroutine fproject_polar_batch_opt_kernel( cmat_exp, refs )
                complex(sp), intent(in)    :: cmat_exp(cdim(1),cdim(2),cdim(3))
                complex(sp), intent(inout) :: refs(:,:,:)
                complex(sp) :: fc
                real        :: base(3), locpd(3), wx(3), wy(3), wz(3), px, py, sx, sy, sz
                integer     :: irot, i, j, kk, l, i0p(3), ix, iy, iz
#ifdef USE_OPENMP_OFFLOAD
                !$omp target teams distribute parallel do collapse(3) &
                !$omp& map(to:cmat_exp(1:cdim(1),1:cdim(2),1:cdim(3)), &
                !$omp& rotmats(1:2,1:3,1:nspace), hkcoords(1:2,1:pftsz,1:nk)) &
                !$omp& map(from:refs(1:pftsz,1:nk,1:nspace)) &
#else
                !$omp parallel do collapse(3) proc_bind(close) default(shared) schedule(static) &
#endif
                !$omp& firstprivate(kb,nspace,nk,pftsz,iwinsz,wdim,clb) &
                !$omp& private(iproj,kk,irot,px,py,locpd,i0p,base,sx,sy,sz,wx,wy,wz,fc,i,j,l,ix,iy,iz)
                do iproj = 1, nspace
                    do kk = 1, nk
                        do irot = 1, pftsz
                            ! padded coordinates
                            px = hkcoords(1,irot,kk)
                            py = hkcoords(2,irot,kk)
                            ! rotated padded coordinates
                            locpd(1) = px*rotmats(1,1,iproj) + py*rotmats(2,1,iproj)
                            locpd(2) = px*rotmats(1,2,iproj) + py*rotmats(2,2,iproj)
                            locpd(3) = px*rotmats(1,3,iproj) + py*rotmats(2,3,iproj)
                            ! 3D grid indices of interpolation window
                            i0p  = nint(locpd) - iwinsz
                            base = real(i0p, sp) - locpd
                            ! weights pre-computation and normalization
                            sx = 0.0
                            sy = 0.0
                            sz = 0.0
                            do i = 1, wdim
                                ! weights
                                wx(i) = apod_fast_device(kb, base(1) + real(i-1, sp))
                                wy(i) = apod_fast_device(kb, base(2) + real(i-1, sp))
                                wz(i) = apod_fast_device(kb, base(3) + real(i-1, sp))
                                ! accumulate sums
                                sx = sx + wx(i)
                                sy = sy + wy(i)
                                sz = sz + wz(i)
                            enddo
                            ! weights normalization
                            wx(:wdim) = wx(:wdim) / sx
                            wy(:wdim) = wy(:wdim) / sy
                            wz(:wdim) = wz(:wdim) / sz
                            ! weighted Fourier interpolation
                            i0p = i0p - clb ! adjust for base-1 reindexing
                            fc  = CMPLX_ZERO
                            do l = 1, wdim
                                do j = 1, wdim
                                    do i = 1, wdim
                                        ix = i0p(1) + i
                                        iy = i0p(2) + j
                                        iz = i0p(3) + l
                                        fc = fc + wx(i) * wy(j) * wz(l) * cmat_exp(ix, iy, iz)
                                    enddo
                                enddo
                            enddo
                            refs(irot,kk,iproj) = fc * PAD_FAC_SCALING
                        enddo
                    enddo
                enddo
#ifdef USE_OPENMP_OFFLOAD
                !$omp end target teams distribute parallel do
#else
                !$omp end parallel do
#endif
            end subroutine fproject_polar_batch_opt_kernel

    end subroutine fproject_polar_batch_opt

end module simple_projector_pft_batch
